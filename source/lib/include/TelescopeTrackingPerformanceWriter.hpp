// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>

#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/Validation/EffPlotTool.hpp"
#include "ActsExamples/Validation/ResPlotTool.hpp"
#include "ActsExamples/Validation/TrackSummaryPlotTool.hpp"
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
  : public ActsExamples::WriterT<std::vector<PixelMultiTrajectory>> {
 public:
  struct Config {
    /// Input (fitted) trajectories collection.
    std::string inputTrajectories;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_telescope_tracking.root";
    /// Plot tool configurations.
    ActsExamples::ResPlotTool::Config resPlotToolConfig;
    ActsExamples::EffPlotTool::Config effPlotToolConfig;
    ActsExamples::TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
  };

  /// Construct from configuration and log level.
  TelescopeTrackingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TelescopeTrackingPerformanceWriter() override;

  /// Finalize plots.
  ActsExamples::ProcessCode endRun() final override;

 private:
  ActsExamples::ProcessCode writeT(
                     const ActsExamples::AlgorithmContext& ctx,
      const std::vector<PixelMultiTrajectory>& trajectories) final override;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for residuals and pulls.
  ActsExamples::ResPlotTool m_resPlotTool;
  ActsExamples::ResPlotTool::ResPlotCache m_resPlotCache;
  /// Plot tool for efficiency
  ActsExamples::EffPlotTool m_effPlotTool;
  ActsExamples::EffPlotTool::EffPlotCache m_effPlotCache;
  /// Plot tool for track hit info
  ActsExamples::TrackSummaryPlotTool m_trackSummaryPlotTool;
  ActsExamples::TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache;
};

}  // namespace FW
