#include "TelMultiTrajectory.hpp"

using namespace Acts::UnitLiterals;


Acts::FreeToBoundMatrix
TelActs::TelMultiTrajectory::freeToCurvilinearJacobian(const Acts::Vector3D &direction) {
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


JsonValue TelActs::TelMultiTrajectory::createJsonValue(JsonAllocator& jsa, Acts::GeometryContext& gctx) const {
  JsonValue js_tracks(rapidjson::kArrayType);
  for (const size_t &trackTip : m_trackTips) {
    // size of trackTips should <= 1 for each seed
    JsonValue js_track(rapidjson::kObjectType);
    JsonValue js_states(rapidjson::kArrayType);
    JsonValue js_states_reverse(rapidjson::kArrayType); //tmp
    m_multiTrajectory.visitBackwards(trackTip, [&](const auto &state) {
                                                 // only fill the track states with non-outlier measurement
                                                 auto typeFlags = state.typeFlags();
                                                 if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
                                                   return true;
                                                 }
                                                 if (!state.hasSmoothed()) {
                                                   return true;
                                                 }

                                                 auto telSurface = state.referenceSurface().getSharedPtr();
                                                 auto telElement = dynamic_cast<const TelActs::TelElement*>(
                                                   telSurface->associatedDetectorElement());
                                                 if(!telElement){
                                                   return true;
                                                 }
                                                 size_t layerid = telElement->id();

                                                 // Get the source link info
                                                 auto meas = std::get<
                                                   Acts::Measurement<TelActs::PixelSourceLink, Acts::BoundIndices,
                                                                     Acts::eBoundLoc0, Acts::eBoundLoc1>>(
                                                                       *state.uncalibrated());
                                                 // Get local position
                                                 Acts::Vector2D pos_hit(meas.parameters()[Acts::eBoundLoc0],
                                                                        meas.parameters()[Acts::eBoundLoc1]);

                                                 // 1) Transform bound parameter to free parameter
                                                 // 1.1)Fist transform the smoothed bound parameters to free parameters
                                                 // to get the position and momentum
                                                 Acts::FreeVector freeParams =
                                                   Acts::detail::transformBoundToFreeParameters(
                                                     *telSurface, gctx, state.smoothed());
                                                 // 1.2)Get the global position, direction, q/p, t etc.
                                                 Acts::Vector3D pos(freeParams[Acts::eFreePos0],
                                                                    freeParams[Acts::eFreePos1],
                                                                    freeParams[Acts::eFreePos2]);
                                                 Acts::Vector3D dir(freeParams[Acts::eFreeDir0],
                                                                    freeParams[Acts::eFreeDir1],
                                                                    freeParams[Acts::eFreeDir2]);
                                                 double p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
                                                 int q = freeParams[Acts::eFreeQOverP]>0 ? 1:-1;
                                                 double t = freeParams[Acts::eFreeTime];

                                                 /// 1.3) Initialize the jacobian from local to the global frame
                                                 Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                                                 // Calculate the jacobian
                                                 telSurface->initJacobianToGlobal(gctx, jacToGlobal,
                                                                                  pos, dir, state.smoothed());
                                                 Acts::FreeSymMatrix freeCovariance =
                                                   jacToGlobal * state.smoothedCovariance() * jacToGlobal.transpose();

                                                 // 2) Transform free parameter to curvilinear parameter
                                                 Acts::FreeToBoundMatrix jacToCurv = freeToCurvilinearJacobian(dir);
                                                 Acts::BoundSymMatrix curvCovariance =
                                                   jacToCurv * freeCovariance * jacToCurv.transpose();

                                                 // Write out the curvilinear parameters and its covariance
                                                 JsonValue js_track_state(rapidjson::kObjectType); //
                                                 js_track_state.AddMember("id", layerid, jsa);

                                                 JsonValue js_hit(rapidjson::kObjectType);
                                                 js_hit.AddMember("x", pos_hit.x(), jsa);
                                                 js_hit.AddMember("y", pos_hit.y(), jsa);
                                                 js_track_state.AddMember("hit", std::move(js_hit), jsa);

                                                 JsonValue js_pos(rapidjson::kObjectType);
                                                 js_pos.AddMember("x", pos.x(), jsa);
                                                 js_pos.AddMember("y", pos.y(), jsa);
                                                 js_pos.AddMember("z", pos.z(), jsa);
                                                 js_track_state.AddMember("pos", std::move(js_pos), jsa);

                                                 JsonValue js_dir(rapidjson::kObjectType);
                                                 js_dir.AddMember("x", dir.x(), jsa);
                                                 js_dir.AddMember("y", dir.y(), jsa);
                                                 js_dir.AddMember("z", dir.z(), jsa);
                                                 js_track_state.AddMember("dir", std::move(js_dir), jsa);

                                                 js_track_state.AddMember("p", p, jsa);
                                                 js_track_state.AddMember("q", p, jsa);
                                                 js_track_state.AddMember("t", t, jsa);

                                                 const double *cov_data = curvCovariance.data();
                                                 JsonValue js_state_cov(rapidjson::kArrayType); //
                                                 js_state_cov.Reserve(Acts::eBoundSize * Acts::eBoundSize, jsa);
                                                 for (size_t n = 0; n < Acts::eBoundSize * Acts::eBoundSize; n++) {
                                                   js_state_cov.PushBack(JsonValue(*(cov_data + n)), jsa);
                                                 }
                                                 js_track_state.AddMember("cov", std::move(js_state_cov), jsa);
                                                 js_states_reverse.PushBack(std::move(js_track_state), jsa);
                                                 return true;
                                               });

    for (size_t i = js_states_reverse.Size(); i > 0; i--) {
      js_states.PushBack(std::move(js_states_reverse[i - 1]), jsa);
    }
    js_track.AddMember("states", std::move(js_states), jsa);
    js_tracks.PushBack(std::move(js_track), jsa);
  }
  return js_tracks;
}


std::unique_ptr<TelActs::TelEvent> TelActs::TelMultiTrajectory::createTelEvent(
  Acts::GeometryContext& gctx,
  const std::map<Acts::GeometryIdentifier, size_t>&  mapSurId2DetId,
  size_t runN, size_t eventN, size_t detSetupN) const{

  std::unique_ptr<TelActs::TelEvent> telEvent(new TelActs::TelEvent{runN, eventN, detSetupN, {}, {}, {}});

  // std::cout<< "================result of fitting"<<std::endl;
  for(size_t indexTrack=0; indexTrack<m_trackTips.size(); indexTrack++){
    std::shared_ptr<TelActs::TelTrajectory> telTraj(new TelActs::TelTrajectory);
    telTraj->TN=indexTrack;

    // std::cout<< "================indexTrack"<< indexTrack <<std::endl;
    // std::cout<< "================tip"<< m_trackTips[indexTrack] <<std::endl;
    m_multiTrajectory.visitBackwards(
      m_trackTips[indexTrack], [&](const auto &state) {
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

void TelActs::TelMultiTrajectory::fillSingleTrack(
  Acts::GeometryContext& gctx,
  const std::map<Acts::GeometryIdentifier, size_t>&  mapSurId2DetId,

  std::vector<size_t>& idMeas,
  std::vector<double>& xMeas,
  std::vector<double>& yMeas,
  std::vector<double>& xResidLocal,
  std::vector<double>& yResidLocal,

  std::vector<size_t>& idFit,
  std::vector<double>& xFitLocal,
  std::vector<double>& yFitLocal,
  std::vector<double>& xFitWorld,
  std::vector<double>& yFitWorld,
  std::vector<double>& zFitWorld,
  size_t indexTrack
  )const
{
  std::cout<< "fill track"<<std::endl;
  
  idMeas.clear();
  xMeas.clear();
  yMeas.clear();
  xResidLocal.clear();
  yResidLocal.clear();

  idFit.clear();
  xFitLocal.clear();
  yFitLocal.clear();
  xFitWorld.clear();
  yFitWorld.clear();
  zFitWorld.clear();

  if(indexTrack >= m_trackTips.size()){
    std::fprintf(stderr, "something wrong, request track #%d in %d tracks\n");
    throw;
  }

  m_multiTrajectory.visitBackwards(
    m_trackTips[indexTrack], [&](const auto &state) {
                      // only fill the track states with non-outlier measurement
                      auto typeFlags = state.typeFlags();
                      if (!state.hasSmoothed()) {
                        return true;
                      }

                      auto telSurface = state.referenceSurface().getSharedPtr();
                      // auto telElement = dynamic_cast<const TelActs::TelElement*>(
                      //   telSurface->associatedDetectorElement());
                      // if(!telElement){
                      //   return true;
                      // }
                      // size_t layerid = telElement->id();
                      size_t layerid = mapSurId2DetId.at(telSurface->geometryId());
                      // Get fit local
                      Acts::Vector2D fit_pos_local(state.smoothed()[Acts::eBoundLoc0],
                                                   state.smoothed()[Acts::eBoundLoc1]);
                      // Get fit global
                      Acts::FreeVector freeParams =
                        Acts::detail::transformBoundToFreeParameters(
                          *telSurface, gctx, state.smoothed());

                      Acts::Vector3D fit_dir_global(freeParams[Acts::eFreeDir0],
                                                    freeParams[Acts::eFreeDir1],
                                                    freeParams[Acts::eFreeDir2]);


                      Acts::Vector3D fit_pos_global(freeParams[Acts::eFreePos0],
                                                    freeParams[Acts::eFreePos1],
                                                    freeParams[Acts::eFreePos2]);
                      idFit.push_back(layerid);
                      xFitLocal.push_back(fit_pos_local(0));
                      yFitLocal.push_back(fit_pos_local(1));
                      xFitWorld.push_back(fit_pos_global(0));
                      yFitWorld.push_back(fit_pos_global(1));
                      zFitWorld.push_back(fit_pos_global(2));

                      if(state.hasUncalibrated()){
                        Acts::Vector2D meas_pos_local;
                        Acts::Vector2D residual_local;
                        meas_pos_local = state.uncalibrated().value();
                        residual_local = meas_pos_local -fit_pos_local;

                        idMeas.push_back(layerid);
                        xMeas.push_back(meas_pos_local(0));
                        yMeas.push_back(meas_pos_local(1));
                        xResidLocal.push_back(residual_local(0));
                        yResidLocal.push_back(residual_local(1));
                      }
                      return true;
                    });


  std::reverse(idMeas.begin(), idMeas.end());
  std::reverse(xMeas.begin(),  xMeas.end());
  std::reverse(yMeas.begin(),  yMeas.end());
  std::reverse(xResidLocal.begin(), xResidLocal.end());
  std::reverse(yResidLocal.begin(), yResidLocal.end());

  std::reverse(idFit.begin(), idFit.end());
  std::reverse(xFitLocal.begin(), xFitLocal.end());
  std::reverse(yFitLocal.begin(), yFitLocal.end());
  std::reverse(xFitWorld.begin(), xFitWorld.end());
  std::reverse(yFitWorld.begin(), yFitWorld.end());
}
