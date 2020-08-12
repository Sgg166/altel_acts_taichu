// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeJsonTrackWriter.hpp"

#include "ACTFW/Utilities/Paths.hpp"
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
                                        Acts::ParDef::eLOC_0,
                                        Acts::ParDef::eLOC_1>;
}

Telescope::TelescopeJsonTrackWriter::TelescopeJsonTrackWriter(
    const Telescope::TelescopeJsonTrackWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT(cfg.inputTrajectories, "TelescopeJsonTrackWriter", level),
      m_cfg(cfg) {

  std::string jsfile_name= m_cfg.outputDir+"/jswrite.json";
  std::cout<< "write jsfile_name "<<jsfile_name<<std::endl;
  m_jsfp = std::fopen((m_cfg.outputDir+"/jswrite.json").c_str(), "w");
  if(!m_jsfp) {
    std::fprintf(stderr, "File opening failed: %s \n", jsfile_name.c_str());
    throw;
  }

  m_jsos.reset( new rapidjson::FileWriteStream(m_jsfp, m_jsbuffer, sizeof(m_jsbuffer)) );
  m_jsw.reset( new rapidjson::Writer< rapidjson::FileWriteStream>(*m_jsos) );
  // m_jsw->SetMaxDecimalPlaces(5);
  m_jsw->StartArray();
}

FW::ProcessCode Telescope::TelescopeJsonTrackWriter::endRun() {
  m_jsw->EndArray();
  m_jsw->Flush();
  std::fclose(m_jsfp);
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode Telescope::TelescopeJsonTrackWriter::writeT
(const FW::AlgorithmContext& context,
 const std::vector<Telescope::PixelMultiTrajectory>& trajectories) {

  JsonAllocator jsa;
  Telescope::JsonValue js_output(rapidjson::kObjectType);
  js_output.AddMember("eventNum", JsonValue(context.eventNumber), jsa );
  Telescope::JsonValue js_tracks(rapidjson::kArrayType);

  for (const auto& traj : trajectories) {
    const auto& [trackTips, mj] = traj.trajectory();
    // Loop over all trajectories in a multiTrajectory
    for (const size_t& trackTip : trackTips) {
      // Collect the trajectory summary info
      auto trajState = Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

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
                          auto psurface = state.referenceSurface().getSharedPtr();
                          if(psurface != m_cfg.outputParaSurface){
                            return true;
                          }

			  // 1) The bound parameters info
                          Acts::BoundParameters boundpara(context.geoContext,
                                                          state.smoothedCovariance(),
                                                          state.smoothed(),
                                                          psurface);
                          Acts::Vector3D pos = boundpara.position();
                          Acts::Vector3D mom = boundpara.momentum();
                          Acts::Vector3D dir = mom.normalized();
                          double q = boundpara.charge();
                          double t = boundpara.time();
			  const auto& boundCovariance = *boundpara.covariance(); 

			  // 2) Transform bound parameter to free parameter
			  //Acts::FreeVector freepara;
                          //freepara << pos[0], pos[1], pos[2], t, dir[0],
                          //dir[1], dir[2], q / mom.norm();

			  /// Initialize the jacobian from local to the global frame
			  Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                          // Calculate the jacobian 
			  psurface->initJacobianToGlobal(context.geoContext, jacToGlobal, pos, dir,
                                     state.smoothed());
                          Acts::FreeSymMatrix freeCovariance = jacToGlobal * boundCovariance * jacToGlobal.transpose();

                          // 3) Transform free parameter to curvilinear parameter
                          Acts::FreeToBoundMatrix jacToCurv = freeToCurvilinearJacobian(dir);
                          Acts::BoundSymMatrix curvCovariance = jacToCurv * freeCovariance * jacToCurv.transpose();

			  // Write out the curvilinear parameters and its covariance
                          Telescope::JsonValue js_track_state(rapidjson::kObjectType); //
                          js_track_state.AddMember("x", JsonValue(pos.x()), jsa);
                          js_track_state.AddMember("y", JsonValue(pos.y()), jsa);
                          js_track_state.AddMember("z", JsonValue(pos.z()), jsa);
                          js_track_state.AddMember("px", JsonValue(mom.x()), jsa);
                          js_track_state.AddMember("py", JsonValue(mom.y()), jsa);
                          js_track_state.AddMember("pz", JsonValue(mom.z()), jsa);
                          js_track_state.AddMember("q", JsonValue(q), jsa);
                          js_track_state.AddMember("t", JsonValue(t), jsa);

                          const double *cov_data = curvCovariance.data();
                          Telescope::JsonValue js_state_cov(rapidjson::kArrayType); //
                          js_state_cov.Reserve(Acts::eBoundParametersSize*Acts::eBoundParametersSize, jsa);
                          for(size_t n=0; n<Acts::eBoundParametersSize*Acts::eBoundParametersSize; n++){
                            js_state_cov.PushBack(JsonValue(*(cov_data+n)), jsa);
                          }
                          js_track_state.AddMember("cov", std::move(js_state_cov), jsa);
                          js_tracks.PushBack(std::move(js_track_state), jsa);
                          return true;
                        });
    }
  }
  js_output.AddMember("tracks", std::move(js_tracks), jsa);

  {
    const std::lock_guard<std::mutex> lock(m_mtx);
    js_output.Accept(*m_jsw);
    rapidjson::PutN(*m_jsos, '\n', 1);
  }

  return FW::ProcessCode::SUCCESS;
}
