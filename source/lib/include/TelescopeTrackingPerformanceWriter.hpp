// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "ACTFW/Validation/TrackSummaryPlotTool.hpp"
#include "PixelMultiTrajectory.hpp"

class TFile;
class TTree;

namespace Telescope{

/// Write out the residual and pull of track parameters and efficiency.
///
/// Efficiency here is the fraction of smoothed tracks compared to all tracks.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class TelescopeTrackingPerformanceWriter final
  : public FW::WriterT<std::vector<PixelMultiTrajectory>> {
 public:
  struct Config {
    /// Input (fitted) trajectories collection.
    std::string inputTrajectories;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_telescope_tracking.root";
    /// Plot tool configurations.
    FW::ResPlotTool::Config resPlotToolConfig;
    FW::EffPlotTool::Config effPlotToolConfig;
    FW::TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
  };

  /// Construct from configuration and log level.
  TelescopeTrackingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TelescopeTrackingPerformanceWriter() override;

  /// Finalize plots.
  FW::ProcessCode endRun() final override;

 private:
  FW::ProcessCode writeT(
                     const FW::AlgorithmContext& ctx,
      const std::vector<PixelMultiTrajectory>& trajectories) final override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for residuals and pulls.
  FW::ResPlotTool m_resPlotTool;
  FW::ResPlotTool::ResPlotCache m_resPlotCache;
  /// Plot tool for efficiency
  FW::EffPlotTool m_effPlotTool;
  FW::EffPlotTool::EffPlotCache m_effPlotCache;
  /// Plot tool for track hit info
  FW::TrackSummaryPlotTool m_trackSummaryPlotTool;
  FW::TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache;
};

}  // namespace FW
