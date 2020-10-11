 
// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "PixelMultiTrajectory.hpp"
#include "PixelSourceLink.hpp"

namespace Telescope {

  class TelescopeTrackFindingAlgorithm final : public ActsExamples::BareAlgorithm {
 public:
    using TrackFinderResult =
      Acts::Result<Acts::CombinatorialKalmanFilterResult<PixelSourceLink>>;
  /// Track finding function that takes input measurements, initial trackstate
  /// and track finder options and returns some track-finding-specific result.
  using CKFOptions =
      Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>;
  using TrackFinderFunction = std::function<TrackFinderResult(
      const std::vector<PixelSourceLink>&, const ActsExamples::TrackParameters&,
      const CKFOptions&)>;

  using SourceLinkTrack = std::vector<PixelSourceLink>;
  using SourceLinkTrackReader =
      std::function<std::vector<SourceLinkTrack>(const std::string&, size_t)>;

  /// Create the track finder function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static TrackFinderFunction makeTrackFinderFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      ActsExamples::Options::BFieldVariant magneticField);

  struct Config {
    /// Inout soucelinks
    std::string inputSourcelinks;
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased track finder function.
    TrackFinderFunction findTracks;
    /// CKF source link selector config
    Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg;
    uint64_t seedSurfaceGeoIDStart;
    uint64_t seedSurfaceGeoIDEnd;
    double seedResX{15 * Acts::UnitConstants::mm};
    double seedResY{15 * Acts::UnitConstants::mm};
    double seedResPhi{0.7};
    double seedResTheta{0.7};
    double seedEnergy{5 * Acts::UnitConstants::GeV};
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TelescopeTrackFindingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(const ActsExamples::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
