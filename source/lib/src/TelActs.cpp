#include "TelActs.hh"


#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"


#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"




Acts::FreeToBoundMatrix
TelActs::freeToCurvilinearJacobian(const Acts::Vector3D &direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0); // == cos(phi) * sin(theta)
  const double y = direction(1); // == sin(phi) * sin(theta)
  const double z = direction(2); // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  Acts::FreeToBoundMatrix jacToCurv = Acts::FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < Acts::s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}


std::unique_ptr<TelActs::TelEvent> TelActs::createTelEvent(
  Acts::GeometryContext& gctx,
  const Acts::CombinatorialKalmanFilterResult<TelActs::TelSourceLink>& ckfResult,
  const std::map<Acts::GeometryIdentifier, size_t>&  mapSurId2DetId,
  size_t runN, size_t eventN, size_t detSetupN) {

  auto& fittedStates = ckfResult.fittedStates;
  auto& trackTips = ckfResult.trackTips;
  auto& fittedParameters = ckfResult.fittedParameters;

  std::unique_ptr<TelActs::TelEvent> telEvent(new TelActs::TelEvent{runN, eventN, detSetupN, {}, {}, {}});

  // std::cout<< "================result of fitting"<<std::endl;
  for(size_t indexTrack=0; indexTrack<trackTips.size(); indexTrack++){
    std::shared_ptr<TelActs::TelTrajectory> telTraj(new TelActs::TelTrajectory);
    telTraj->TN=indexTrack;

    // std::cout<< "================indexTrack"<< indexTrack <<std::endl;
    // std::cout<< "================tip"<< trackTips[indexTrack] <<std::endl;
    fittedStates.visitBackwards(
      trackTips[indexTrack], [&](const auto &state) {
                                 // only fill the track states with non-outlier measurement
                                 auto typeFlags = state.typeFlags();
                                 if (!state.hasSmoothed() && !state.hasPredicted()) {
                                   return true;
                                 }

                                 std::shared_ptr<TelActs::TelHit> telHit;
                                 auto telSurface = state.referenceSurface().getSharedPtr();
                                 // std::cout<< "telSurface->geometryId() "<<telSurface->geometryId()<<std::endl;
                                 // std::cout<< telSurface->center(gctx)<<std::endl;
                                 // for(auto [aSurId, aDetId]:mapSurId2DetId){
                                 //   std::cout<< aSurId <<" sur=det "<< aDetId<<std::endl;
                                 // }
                                 size_t id = mapSurId2DetId.at(telSurface->geometryId());

                                 Acts::Vector2D fit_pos_local;
                                 Acts::FreeVector freeParams;
                                 if(state.hasSmoothed()){
                                   fit_pos_local=Acts::Vector2D(state.smoothed()[Acts::eBoundLoc0],
                                                                state.smoothed()[Acts::eBoundLoc1]);
                                   freeParams = Acts::detail::transformBoundToFreeParameters(
                                     *telSurface, gctx, state.smoothed());
                                 }
                                 else{
                                   fit_pos_local=Acts::Vector2D(state.predicted()[Acts::eBoundLoc0],
                                                                state.predicted()[Acts::eBoundLoc1]);
                                   freeParams = Acts::detail::transformBoundToFreeParameters(
                                     *telSurface, gctx, state.predicted());
                                 }

                                 Acts::Vector3D fit_dir_global(freeParams[Acts::eFreeDir0],
                                                               freeParams[Acts::eFreeDir1],
                                                               freeParams[Acts::eFreeDir2]);

                                 Acts::Vector3D fit_pos_global(freeParams[Acts::eFreePos0],
                                                               freeParams[Acts::eFreePos1],
                                                               freeParams[Acts::eFreePos2]);

                                 Acts::Vector3D fit_dir_local = telSurface->transform(gctx).rotation().transpose() * fit_dir_global;

                                 std::shared_ptr<TelActs::TelHitMeasure> telHitMeas;
                                 if(state.hasUncalibrated()){
                                   Acts::Vector2D meas_pos_local;
                                   Acts::Vector2D residual_local;
                                   meas_pos_local = state.uncalibrated().value();
                                   residual_local = meas_pos_local -fit_pos_local;

                                   telHitMeas.reset(new TelActs::TelHitMeasure{id, {meas_pos_local(0), meas_pos_local(1)}, {}});
                                 }

                                 std::shared_ptr<TelActs::TelHitFit> telHitFit;
                                 telHitFit.reset(new TelActs::TelHitFit{
                                     id,
                                     {fit_pos_local(0), fit_pos_local(1)},
                                     {fit_dir_local(0), fit_dir_local(1), fit_dir_local(2)}, // local dir
                                     {fit_pos_global(0), fit_pos_global(1), fit_pos_global(2)},
                                     {fit_dir_global(0), fit_dir_global(1), fit_dir_global(2)},
                                     telHitMeas});

                                 telHit.reset(new TelActs::TelHit{id, telHitFit, nullptr});
                                 telTraj->Hs.push_back(telHit);
                                 return true;
                               });
    std::reverse(telTraj->Hs.begin(), telTraj->Hs.end());
    telEvent->Ts.push_back(telTraj);
  }
  return telEvent;
}



std::unique_ptr<TelActs::TelEvent> TelActs::createTelEvent(
  const JsonValue& js,
  std::vector<std::shared_ptr<Acts::PlaneLayer>>& planeLayers,
  const std::map<Acts::GeometryIdentifier, size_t>&  mapGeoId2DetId,
  size_t runN, size_t eventN, size_t detSetupN){

  std::map<size_t, std::shared_ptr<Acts::PlaneLayer>> mapDetId2PlaneLayer;
  for(auto &aPlaneLayer: planeLayers){
    Acts::GeometryIdentifier geoId= aPlaneLayer->geometryId();
    size_t detId = mapGeoId2DetId.at(geoId);
    mapDetId2PlaneLayer[detId] = aPlaneLayer;
  }

  std::unique_ptr<TelActs::TelEvent> telEvent(new TelActs::TelEvent{runN, eventN, detSetupN, {}, {}, {}});

  const auto &layers = js["layers"];
  for (const auto &layer : layers.GetArray()) {
    size_t detId = layer["ext"].GetUint();
    auto it = mapDetId2PlaneLayer.find(detId);
    if(it == mapDetId2PlaneLayer.end()){
      continue;
    }
    uint16_t tri = layer["tri"].GetUint();

    for (const auto &hit : layer["hit"].GetArray()) {
      double hitMeasU = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
      double hitMeasV = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;

      std::vector<std::shared_ptr<TelActs::TelRawMeasure>> rawMeasCol;
      for(const auto &pix :  hit["pix"].GetArray()){
        std::shared_ptr<TelActs::TelRawMeasure> rawMeas(new TelActs::TelRawMeasure);
        rawMeas->uvdcS[0] = int16_t(pix[0].GetInt());
        rawMeas->uvdcS[1] = int16_t(pix[0].GetInt());
        rawMeas->uvdcS[2] = int16_t(detId);
        rawMeas->uvdcS[3] = int16_t(tri);
        rawMeasCol.push_back(rawMeas);
      }
      telEvent->HMs.emplace_back(new TelActs::TelHitMeasure{detId, {hitMeasU, hitMeasV}, rawMeasCol});
    }
  }
  return telEvent;
}



using namespace Acts::UnitLiterals;

void TelActs::matchAddExtraHitMeas(
  Acts::GeometryContext& gctx,
  const std::map<Acts::GeometryIdentifier, size_t>&  mapSurId2DetId,
  std::shared_ptr<TelActs::TelEvent> telEvent,
  const std::vector<TelActs::TelSourceLink>& sourcelinksTargets
  ){
    for(auto &telTraj: telEvent->Ts){
      size_t fittedHitNum = telTraj->numberHitFitByMeas();
      if(fittedHitNum<5){
        continue;
      }

      for(const auto &sl : sourcelinksTargets){
        size_t id = mapSurId2DetId.at(sl.referenceSurface().geometryId());

        auto aTelHit = telTraj->hit(id);
        if(!aTelHit || !aTelHit->hasHitFit() || aTelHit->isFittedFromMeasure()){
          continue;
        }

        Acts::Vector2D xy_fit(aTelHit->HF->PLs[0], aTelHit->HF->PLs[1]);
        Acts::Vector2D xy_meas = sl.value();
        Acts::Vector2D xy_resid = xy_meas - xy_fit;
        if(xy_resid.norm()>0.1_mm){
          continue;
        }

        aTelHit->HM.reset(new TelActs::TelHitMeasure{id, {xy_meas(0), xy_meas(1)}, {}});
      }
    }
}

std::vector<TelActs::TelSourceLink> TelActs::createSourceLink(
  const std::map<Acts::GeometryIdentifier, size_t>&  mapGeoId2DetId,
  std::vector<std::shared_ptr<Acts::PlaneLayer>>& planeLayers,
  std::shared_ptr<TelActs::TelEvent> telEvent){

  std::map<size_t, std::shared_ptr<Acts::PlaneLayer>> mapDetId2PlaneLayer;
  for(auto &aPlaneLayer: planeLayers){
    Acts::GeometryIdentifier geoId= aPlaneLayer->geometryId();
    size_t detId = mapGeoId2DetId.at(geoId);
    mapDetId2PlaneLayer[detId] = aPlaneLayer;
  }

  std::vector<TelActs::TelSourceLink> sourceLinks;
  auto & hitsMeas= telEvent->HMs;
  for(auto &aHitMeas : hitsMeas){
    auto detN = aHitMeas->DN;
    auto it = mapDetId2PlaneLayer.find(detN);
    if(it==mapDetId2PlaneLayer.end()){
      continue;
    }

    auto u = aHitMeas->PLs[0];
    auto v = aHitMeas->PLs[1];

    Acts::Vector2D hitLoc;
    hitLoc << u, v;
    Acts::BoundMatrix hitCov = Acts::BoundMatrix::Zero();
    double hitResolU = 6_um;
    double hitResolV = 6_um;
    hitCov(0, 0) = hitResolU * hitResolU;
    hitCov(1, 1) = hitResolV * hitResolV;
    sourceLinks.emplace_back(*(it->second), hitLoc, hitCov);
  }
  return sourceLinks;
}
