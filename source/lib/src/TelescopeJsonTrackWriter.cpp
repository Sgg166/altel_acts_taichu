// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeJsonTrackWriter.hpp"

#include "ActsExamples/Utilities/Paths.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
Acts::FreeToBoundMatrix freeToCurvilinearJacobian(const Acts::Vector3D& direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
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

namespace Telescope{
  using Measurement = Acts::Measurement<Telescope::PixelSourceLink,
        Acts::BoundIndices,                                 
	Acts::eBoundLoc0,
                                        Acts::eBoundLoc1>;
}

Telescope::TelescopeJsonTrackWriter::TelescopeJsonTrackWriter(
    const Telescope::TelescopeJsonTrackWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT(cfg.inputTrajectories, "TelescopeJsonTrackWriter", level),
      m_cfg(cfg) {

  for(auto &[id, surface]: m_cfg.trackSurfaces){
    m_surface_id_map[surface] = id;
  }

  std::cout<< "write jsfile_name "<<m_cfg.outputFileName<<std::endl;
  m_jsfp = std::fopen((m_cfg.outputDir+ "/"+ m_cfg.outputFileName).c_str(), "w");
  if(!m_jsfp) {
    std::fprintf(stderr, "File opening failed: %s \n", m_cfg.outputFileName.c_str());
    throw;
  }

  m_jsos.reset( new rapidjson::FileWriteStream(m_jsfp, m_jsbuffer, sizeof(m_jsbuffer)) );
  m_jsw.reset( new rapidjson::Writer< rapidjson::FileWriteStream>(*m_jsos) );
  // m_jsw->SetMaxDecimalPlaces(5);
  m_jsw->StartArray();
}

ActsExamples::ProcessCode Telescope::TelescopeJsonTrackWriter::endRun() {
  m_jsw->EndArray();
  m_jsw->Flush();
  std::fclose(m_jsfp);
  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode Telescope::TelescopeJsonTrackWriter::writeT
(const ActsExamples::AlgorithmContext& context,
 const std::vector<Telescope::PixelMultiTrajectory>& trajectories) {

  JsonAllocator jsa;
  Telescope::JsonValue js_output(rapidjson::kObjectType);
  js_output.AddMember("eventNum", JsonValue(context.eventNumber), jsa );
  Telescope::JsonValue js_tracks(rapidjson::kArrayType);

  for (const auto& traj : trajectories) {
    const auto& [trackTips, mj] = traj.trajectory();
    // Loop over all trajectories in a multiTrajectory
    // size of trackTips should <= 1
    for (const size_t& trackTip : trackTips) {
      // Collect the trajectory summary info
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
      if(trajState.nMeasurements != m_cfg.trackSurfaces.size()){
        continue;
      }
      Telescope::JsonValue js_track(rapidjson::kObjectType);
      Telescope::JsonValue js_states(rapidjson::kArrayType);
      Telescope::JsonValue js_states_reverse(rapidjson::kArrayType);
      mj.visitBackwards(trackTip,
                        [&](const auto& state) {
                          // only fill the track states with non-outlier measurement
                          auto typeFlags = state.typeFlags();
                          if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
                            return true;
                          }
                          if(!state.hasSmoothed()) {
                            return true;
                          }

                          auto state_surface = state.referenceSurface().getSharedPtr();
                          auto surface_it = m_surface_id_map.find(state_surface);
                          if(surface_it == m_surface_id_map.end() ){
                            return true;
                          }
                          size_t layerid = surface_it->second;


                          // Get the source link info
                          auto meas = std::get<Acts::Measurement
                                               <Telescope::PixelSourceLink, Acts::BoundIndices, Acts::eBoundLoc0,Acts::eBoundLoc1>
                                               >(*state.uncalibrated());

                          // Get local position
                          Acts::Vector2D pos_local(meas.parameters()[Acts::eBoundLoc0],
                                                   meas.parameters()[Acts::eBoundLoc1]);

                          // std::cout<<"local position of hit on layer " << layerid<<" : " << local.x()<<", "<<local.y()<<std::endl;

                          //The bound parameters info
                          //Acts::BoundTrackParameters boundpara(state_surface,
                          //                                state.smoothed(),
                          //                                state.smoothedCovariance()
	                  //						  );

			  // 1) Transform bound parameter to free parameter
                          // 1.1)Fist transform the smoothed bound parameters to free parameters to get the position and momentum 
			  Acts::FreeVector freeParams =
                          Acts::detail::transformBoundToFreeParameters(*state_surface, context.geoContext,
                                                         state.smoothed());
			  // 1.2)Get the global position, direction, q/p, t etc. 
			  Acts::Vector3D pos(freeParams[Acts::eFreePos0], freeParams[Acts::eFreePos1], freeParams[Acts::eFreePos2]);
                          Acts::Vector3D dir(freeParams[Acts::eFreeDir0], freeParams[Acts::eFreeDir1], freeParams[Acts::eFreeDir2]);
			  double p = std::abs(1 / freeParams[Acts::eFreeQOverP]);;
                          double q = p*freeParams[Acts::eFreeQOverP];
                          double t = freeParams[Acts::eFreeTime];

                          /// 1.3) Initialize the jacobian from local to the global frame
                          Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                          // Calculate the jacobian
                          state_surface->initJacobianToGlobal(context.geoContext, jacToGlobal, pos, dir,
                                                              state.smoothed());
                          Acts::FreeSymMatrix freeCovariance = jacToGlobal * state.smoothedCovariance()* jacToGlobal.transpose();

                          // 2) Transform free parameter to curvilinear parameter
                          Acts::FreeToBoundMatrix jacToCurv = freeToCurvilinearJacobian(dir);
                          Acts::BoundSymMatrix curvCovariance = jacToCurv * freeCovariance * jacToCurv.transpose();

                          // Write out the curvilinear parameters and its covariance
                          Telescope::JsonValue js_track_state(rapidjson::kObjectType); //
                          js_track_state.AddMember("id", JsonValue(layerid), jsa);
                          js_track_state.AddMember("x", JsonValue(pos.x()), jsa);
                          js_track_state.AddMember("y", JsonValue(pos.y()), jsa);
                          js_track_state.AddMember("z", JsonValue(pos.z()), jsa);
                          js_track_state.AddMember("lx", JsonValue(pos_local.x()), jsa);
                          js_track_state.AddMember("ly", JsonValue(pos_local.y()), jsa);
                          js_track_state.AddMember("dx", JsonValue(dir.x()), jsa);
                          js_track_state.AddMember("dy", JsonValue(dir.y()), jsa);
                          js_track_state.AddMember("dz", JsonValue(dir.z()), jsa);
                          js_track_state.AddMember("p", JsonValue(p), jsa);
                          js_track_state.AddMember("q", JsonValue(q), jsa);
                          js_track_state.AddMember("t", JsonValue(t), jsa);

                          const double *cov_data = curvCovariance.data();
                          Telescope::JsonValue js_state_cov(rapidjson::kArrayType); //
                          js_state_cov.Reserve(Acts::eBoundSize*Acts::eBoundSize, jsa);
                          for(size_t n=0; n<Acts::eBoundSize*Acts::eBoundSize; n++){
                            js_state_cov.PushBack(JsonValue(*(cov_data+n)), jsa);
                          }
                          js_track_state.AddMember("cov", std::move(js_state_cov), jsa);
                          js_states_reverse.PushBack(std::move(js_track_state), jsa);
                          return true;
                        });

      for(size_t i = js_states_reverse.Size(); i>0; i-- ){
        js_states.PushBack(std::move(js_states_reverse[i-1]), jsa);
      }
      js_track.AddMember("states", std::move(js_states), jsa);
      js_tracks.PushBack(std::move(js_track), jsa);
    }
  }
  js_output.AddMember("tracks", std::move(js_tracks), jsa);

  {
    const std::lock_guard<std::mutex> lock(m_mtx);
    js_output.Accept(*m_jsw);
    rapidjson::PutN(*m_jsos, '\n', 2);
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
