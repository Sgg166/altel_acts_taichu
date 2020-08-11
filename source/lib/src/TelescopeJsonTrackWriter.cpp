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

#include <ios>
#include <iostream>
#include <stdexcept>

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

			  // The bound parameters info
                          Acts::BoundParameters boundpara(context.geoContext,
                                                          state.smoothedCovariance(),
                                                          state.smoothed(),
                                                          psurface);
                          Acts::Vector3D gpos = boundpara.position();
                          Acts::Vector3D dir = boundpara.momentum().normalized();
                          const auto& boundCovariance = *boundpara.covariance(); 

			  // Transform to free parameter
			  Acts::FreeVector freepara;
                          freepara << gpos[0], gpos[1], gpos[2], boundpara.time(), dir[0],
                          dir[1], dir[2], boundpara.charge() / boundpara.momentum().norm();

			  /// Initialize the jacobian from local to the global frame
			  Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                          // Calculate the jacobian 
			  psurface->initJacobianToGlobal(context.geoContext, jacToGlobal, gpos, dir,
                                     state.smoothed());
                          Acts::FreeSymMatrix freeCovariance = jacToGlobal * boundCovariance * jacToGlobal.transpose();

			  // Write out the free parameters and its covariance
                          Telescope::JsonValue js_track_state(rapidjson::kObjectType); //
                          js_track_state.AddMember("x", JsonValue(freepara(Acts::eFreePos0)), jsa);
                          js_track_state.AddMember("y", JsonValue(freepara(Acts::eFreePos1)), jsa);
                          js_track_state.AddMember("z", JsonValue(freepara(Acts::eFreePos2)), jsa);
                          js_track_state.AddMember("t", JsonValue(freepara(Acts::eFreeTime)), jsa);
                          js_track_state.AddMember("dirx", JsonValue(freepara(Acts::eFreeDir0)), jsa);
                          js_track_state.AddMember("diry", JsonValue(freepara(Acts::eFreeDir1)), jsa);
                          js_track_state.AddMember("dirz", JsonValue(freepara(Acts::eFreeDir2)), jsa);
                          js_track_state.AddMember("q/p", JsonValue(freepara(Acts::eFreeQOverP)), jsa);

                          const double *cov_data = freeCovariance.data();
                          Telescope::JsonValue js_state_cov(rapidjson::kArrayType); //
                          js_state_cov.Reserve(Acts::eFreeParametersSize*Acts::eFreeParametersSize, jsa);
                          for(size_t n=0; n<Acts::eFreeParametersSize*Acts::eFreeParametersSize; n++){
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
