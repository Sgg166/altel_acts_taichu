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

  std::FILE* m_jsfp = std::fopen((m_cfg.outputDir+"/jswrite.json").c_str(), "w");
  m_jsos.reset( new rapidjson::FileWriteStream(m_jsfp, m_jsbuffer, sizeof(m_jsbuffer)) );
  m_jsw.reset( new rapidjson::Writer< rapidjson::FileWriteStream>(*m_jsos) );
  m_jsw->SetMaxDecimalPlaces(5);
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


  Telescope::JsonValue js_output(rapidjson::kObjectType);
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

                          Acts::BoundParameters boundpara(context.geoContext,
                                                          state.smoothedCovariance(),
                                                          state.smoothed(),
                                                          psurface);

                          Acts::Vector3D gpos = boundpara.position();

                          Telescope::JsonValue js_track_state(rapidjson::kObjectType); //
                          js_track_state.AddMember("x", JsonValue(gpos.x()), *m_jsa);
                          js_track_state.AddMember("y", JsonValue(gpos.y()), *m_jsa);
                          js_track_state.AddMember("z", JsonValue(gpos.z()), *m_jsa);

                          const double *cov_data = boundpara.covariance()->data();
                          Telescope::JsonValue js_state_cov(rapidjson::kArrayType); //
                          js_state_cov.Reserve(36, *m_jsa);
                          for(size_t n=0; n<36; n++){
                            js_state_cov.PushBack(JsonValue(*(cov_data+n)), *m_jsa);
                          }
                          js_track_state.AddMember("cov", std::move(js_state_cov), *m_jsa);
                          js_tracks.PushBack(std::move(js_track_state), *m_jsa);
                          return true;
                        });
    }
  }

  return FW::ProcessCode::SUCCESS;
}
