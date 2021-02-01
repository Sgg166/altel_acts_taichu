#include "TelActs.hh"


#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"


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



void TelActs::fillTelTrajectories(Acts::GeometryContext& gctx,
                                  const Acts::CombinatorialKalmanFilterResult<TelActs::TelSourceLink>& ckfResult,
                                  std::shared_ptr<altel::TelEvent> telEvent,
                                  const std::map<Acts::GeometryIdentifier, size_t>&  mapSurId2DetId){

  auto& fittedStates = ckfResult.fittedStates;
  auto& trackTips = ckfResult.trackTips;
  auto& fittedParameters = ckfResult.fittedParameters;

  // std::cout<< "================result of fitting"<<std::endl;
  for(size_t indexTrack=0; indexTrack<trackTips.size(); indexTrack++){
    std::shared_ptr<altel::TelTrajectory> telTraj(new altel::TelTrajectory);
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

                                 std::shared_ptr<altel::TelTrajHit> telTrajHit;
                                 auto telSurface = state.referenceSurface().getSharedPtr();
                                 size_t id = mapSurId2DetId.at(telSurface->geometryId());


                                 Acts::Vector2D fit_pos_local;
                                 Acts::FreeVector freeParams;
                                 if(state.hasSmoothed()){
                                   fit_pos_local=Acts::Vector2D(state.smoothed()[Acts::eBoundLoc0],
                                                                state.smoothed()[Acts::eBoundLoc1]);
                                   freeParams = Acts::detail::transformBoundToFreeParameters(
                                     *telSurface, gctx, state.smoothed());
                                 }
                                 // if(state.hasFiltered()){
                                 //   fit_pos_local=Acts::Vector2D(state.filtered()[Acts::eBoundLoc0],
                                 //                                state.filtered()[Acts::eBoundLoc1]);
                                 //   freeParams = Acts::detail::transformBoundToFreeParameters(
                                 //     *telSurface, gctx, state.filtered());
                                 // }
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

                                 std::shared_ptr<altel::TelMeasHit> telMeasHit;
                                 if(state.hasUncalibrated() && typeFlags.test(Acts::TrackStateFlag::MeasurementFlag) ){
                                   // do not get outlier meas
                                   Acts::Vector2D meas_pos_local;
                                   Acts::Vector2D residual_local;
                                   meas_pos_local = state.uncalibrated().value();
                                   residual_local = meas_pos_local -fit_pos_local;
                                   telMeasHit = state.uncalibrated().m_hitMeas;
                                 }
                                 std::shared_ptr<altel::TelFitHit> telFitHit;
                                 telFitHit.reset(new altel::TelFitHit(
                                     uint16_t(id),
                                     fit_pos_local(0), fit_pos_local(1),
                                     fit_pos_global(0), fit_pos_global(1), fit_pos_global(2),
                                     fit_dir_global(0), fit_dir_global(1), fit_dir_global(2),
                                     telMeasHit));

                                 telTrajHit.reset(new altel::TelTrajHit{uint16_t(id), telFitHit, nullptr});
                                 telTraj->THs.push_back(telTrajHit);
                                 return true;
                               });
    std::reverse(telTraj->THs.begin(), telTraj->THs.end());
    telEvent->TJs.push_back(telTraj);
  }
}

std::unique_ptr<altel::TelEvent> TelActs::createTelEvent(
  const JsonValue& js,
  size_t runN, size_t eventN, size_t detSetupN,
  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer
){

  std::unique_ptr<altel::TelEvent> telEvent(new altel::TelEvent(runN, eventN, detSetupN, 0));

  const auto &layers = js["layers"];
  for (const auto &layer : layers.GetArray()) {
    size_t detId = layer["ext"].GetUint();
    auto it = mapDetId2PlaneLayer.find(detId);
    if(it == mapDetId2PlaneLayer.end()){
      continue;
    }
    uint16_t tri = layer["tri"].GetUint();
    telEvent->CK=tri;

    for (const auto &hit : layer["hit"].GetArray()) {
      double hitMeasU = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
      double hitMeasV = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
      // if(hitMeasV>5){
      //   continue;
      // }
      std::vector<altel::TelMeasRaw> rawMeasCol;
      for(const auto &pix :  hit["pix"].GetArray()){
        altel::TelMeasRaw measRaw(uint16_t(pix[0].GetInt()),
                                  uint16_t(pix[1].GetInt()),
                                  uint16_t(detId),
                                  uint16_t(tri));
        rawMeasCol.push_back(measRaw);
      }
      telEvent->MHs.emplace_back(new altel::TelMeasHit(detId, hitMeasU, hitMeasV, rawMeasCol));
    }
  }
  return telEvent;
}


std::unique_ptr<altel::TelEvent> TelActs::createTelEvent(
  const JsonValue& js,
  size_t runN, size_t eventN, size_t detSetupN
){

  std::unique_ptr<altel::TelEvent> telEvent(new altel::TelEvent(runN, eventN, detSetupN, 0));

  const auto &layers = js["layers"];
  for (const auto &layer : layers.GetArray()) {
    size_t detId = layer["ext"].GetUint();
    uint16_t tri = layer["tri"].GetUint();
    telEvent->CK=tri;

    for (const auto &hit : layer["hit"].GetArray()) {
      double hitMeasU = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
      double hitMeasV = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
      std::vector<altel::TelMeasRaw> rawMeasCol;
      for(const auto &pix :  hit["pix"].GetArray()){
        altel::TelMeasRaw measRaw(uint16_t(pix[0].GetInt()),
                                  uint16_t(pix[1].GetInt()),
                                  uint16_t(detId),
                                  uint16_t(tri));
        rawMeasCol.push_back(measRaw);
      }
      telEvent->MHs.emplace_back(new altel::TelMeasHit(detId, hitMeasU, hitMeasV, rawMeasCol));
    }
  }
  return telEvent;
}




using namespace Acts::UnitLiterals;


void TelActs::mergeAndMatchExtraTelEvent(std::shared_ptr<altel::TelEvent> aEvent,
                                         std::shared_ptr<altel::TelEvent> extraEvent,
                                         double maxMatchDist,
                                         double minFitHitsPerTraj){

  // TODO: test if existing, however it does not hurt ttree write
  aEvent->MRs.insert(aEvent->MRs.end(), extraEvent->MRs.begin(), extraEvent->MRs.end());
  aEvent->MHs.insert(aEvent->MHs.end(), extraEvent->MHs.begin(), extraEvent->MHs.end());

  auto &extraMeasHits = extraEvent->MHs;

  for(auto &aMeasHit : extraMeasHits){
    size_t id = aMeasHit->DN;

    double bestMatchedDist = HUGE_VAL;
    std::shared_ptr<altel::TelTrajHit> bestMatchedTrajHit;
    for(auto &aTraj: aEvent->TJs){
      size_t fittedHitNum = aTraj->numOriginMeasHit();
      if(fittedHitNum<minFitHitsPerTraj){
        continue;
      }
      auto aTrajHit = aTraj->trajHit(id);
      if(!aTrajHit || !aTrajHit->hasFitHit() || aTrajHit->hasOriginMeasHit()){
        continue;
      }
      Acts::Vector2D xy_fit(aTrajHit->FH->PLs[0], aTrajHit->FH->PLs[1]);
      Acts::Vector2D xy_meas (aMeasHit->PLs[0], aMeasHit->PLs[1]);
      Acts::Vector2D xy_resid = xy_meas - xy_fit;
      double dist = xy_resid.norm();
      if(dist>maxMatchDist || dist>bestMatchedDist){
        continue;
      }

      bestMatchedDist = dist;
      bestMatchedTrajHit = aTrajHit;
    }

    if(bestMatchedTrajHit){
      bestMatchedTrajHit->MM=aMeasHit;
    }

  }
}



void TelActs::mergeAndMatchExtraTelEventForTraj(std::shared_ptr<altel::TelEvent> aEvent,
                                                std::shared_ptr<altel::TelEvent> extraEvent,
                                                double maxMatchDist,
                                                double minFitHitsPerTraj){

  // TODO: test if existing, however it does not hurt ttree write
  aEvent->MRs.insert(aEvent->MRs.end(), extraEvent->MRs.begin(), extraEvent->MRs.end());
  aEvent->MHs.insert(aEvent->MHs.end(), extraEvent->MHs.begin(), extraEvent->MHs.end());

  auto &extraMeasHits = extraEvent->MHs;

  for(auto &aTraj: aEvent->TJs){
    size_t fittedHitNum = aTraj->numOriginMeasHit();
    if(fittedHitNum<minFitHitsPerTraj){
      continue;
    }

    for(auto &aTrajHit: aTraj->THs){
      if(!aTrajHit || !aTrajHit->hasFitHit() || aTrajHit->hasOriginMeasHit()){
        continue;
      }
      size_t aTrajHitId = aTrajHit->DN;
      double bestMatchedDist = HUGE_VAL;
      std::shared_ptr<altel::TelMeasHit> bestMatchedMeasHit;

      for(auto &aMeasHit : extraMeasHits){
        size_t aMeasHitId = aMeasHit->DN;
        if( aMeasHitId != aTrajHitId ){
          continue;
        }
        Acts::Vector2D xy_fit(aTrajHit->FH->PLs[0], aTrajHit->FH->PLs[1]);
        Acts::Vector2D xy_meas (aMeasHit->PLs[0], aMeasHit->PLs[1]);
        Acts::Vector2D xy_resid = xy_meas - xy_fit;
        double dist = xy_resid.norm();
        if(dist>maxMatchDist || dist>bestMatchedDist){
          continue;
        }
        bestMatchedDist = dist;
        bestMatchedMeasHit = aMeasHit;
      }
      if(bestMatchedMeasHit){
        aTrajHit->MM=bestMatchedMeasHit;
      }
    }
  }
}


std::vector<TelActs::TelSourceLink> TelActs::createSourceLinks(
  std::shared_ptr<altel::TelEvent> telEvent,
  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer){

  std::vector<TelActs::TelSourceLink> sourceLinks;
  auto & hitsMeas= telEvent->MHs;
  for(auto &aHitMeas : hitsMeas){
    sourceLinks.emplace_back(aHitMeas, mapDetId2PlaneLayer);
  }
  return sourceLinks;
}


std::pair<size_t, std::shared_ptr<Acts::PlaneLayer>> TelActs::createPlaneLayer(const JsonValue& js_det){
  size_t id = js_det["id"].GetUint();
  double cx = js_det["center"]["x"].GetDouble();
  double cy = js_det["center"]["y"].GetDouble();
  double cz = js_det["center"]["z"].GetDouble();
  double rx = js_det["rotation"]["x"].GetDouble();
  double ry = js_det["rotation"]["y"].GetDouble();
  double rz = js_det["rotation"]["z"].GetDouble();
  double sx = js_det["size"]["x"].GetDouble();
  double sy = js_det["size"]["y"].GetDouble();
  double sz = js_det["size"]["z"].GetDouble();

  double layerThickness = 80_um;

  std::shared_ptr<Acts::PlanarBounds> pBounds(new Acts::RectangleBounds(
                                                sx * Acts::UnitConstants::mm ,  // *0.5
                                                sy * Acts::UnitConstants::mm)); // *0.5
  //NOTE: workaround, enlarge sensor size x2 to prevent buggy decision of reaching end of tracker

  Acts::Vector3D translation(cx, cy, cz);
  Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
  Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
  Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
  Acts::Rotation3D rotation = rotZ * rotY * rotX;

  Acts::Transform3D xBeamRotation(Acts::Transform3D::Identity());
  xBeamRotation.linear()<< // y-z-x
    0, 0, 1,
    1, 0, 0,
    0, 1, 0;

  auto layerTransform = std::make_shared<Acts::Transform3D>(xBeamRotation * Acts::Translation3D(translation) * rotation);

  std::cout<< layerTransform->matrix()<<std::endl;

  // Acts::Transform3D layerTransform = xBeamRotation * Acts::Translation3D(translation) * rotation;
  std::shared_ptr<Acts::PlaneLayer> planeLayer = std::dynamic_pointer_cast<Acts::PlaneLayer>(
    Acts::PlaneLayer::create(*layerTransform, pBounds, nullptr, 0.5_mm));

  Acts::Material silicon = Acts::Material::fromMolarDensity(
    9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);
  auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
    Acts::MaterialSlab(silicon, layerThickness));

  planeLayer->surfaceRepresentation().assignSurfaceMaterial(surfaceMaterial);

  return {id, planeLayer};
}


std::shared_ptr<Acts::TrackingGeometry>
TelActs::createWorld(Acts::GeometryContext &gctx, double sizex, double sizey, double sizez,
                     const std::vector<std::shared_ptr<const Acts::PlaneLayer>>& planeLayers){

  auto tracker_cuboid = std::make_shared<Acts::CuboidVolumeBounds>(
      sizex/2., sizey/2.,  sizez/2.);

  std::vector<std::shared_ptr<const Acts::Layer>> layers;
  for(auto& pl : planeLayers){
    layers.push_back(pl);
  }

  Acts::LayerArrayCreator layerArrayCreator({});
  auto layer_array_binned = layerArrayCreator.layerArray(
    gctx, layers, sizex/-2, sizex/2,
    Acts::BinningType::arbitrary, Acts::BinningValue::binX);


  auto trackVolume = Acts::TrackingVolume::create(
      Acts::Transform3D::Identity(), tracker_cuboid, nullptr,
      std::move(layer_array_binned), nullptr, {}, "Tracker");

  std::shared_ptr<Acts::TrackingGeometry> geo_world;
  geo_world.reset(new Acts::TrackingGeometry(trackVolume));

  return geo_world;
}
