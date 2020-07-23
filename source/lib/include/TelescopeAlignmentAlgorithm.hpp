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

#include "Alignment.hpp"
#include "PixelSourceLink.hpp"

#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"

namespace Telescope{

  class TelescopeAlignmentAlgorithm final : public FW::BareAlgorithm {
 public:
    using AlignResult = Acts::Result<FW::AlignmentResult>;
  /// Fit function that takes input measurements, initial trackstate and fitter
  /// options and returns some fit-specific result.
  using AlignmentFunction = std::function<AlignResult(
      const std::vector<std::vector<PixelSourceLink>>&,
      const std::vector<Acts::CurvilinearParameters>&,
      const FW::AlignmentOptions<
          Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>&)>;
  using SourceLinkTrack = std::vector<PixelSourceLink>;
  using SourceLinkTrackReader =
      std::function<std::vector<SourceLinkTrack>(const std::string&, size_t)>;

  /// Create the fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static AlignmentFunction makeAlignmentFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      FW::Options::BFieldVariant magneticField, Acts::Logging::Level lvl);

  struct Config {
    /// Input data file name.
    std::string inputFileName;
    /// Input data reader
    SourceLinkTrackReader trackReader;
    /// Number of tracks
    size_t maxNumTracks = std::numeric_limits<size_t>::max();
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased fitter function.
    AlignmentFunction align;
    /// The alignd transform updater
    FW::AlignedTransformUpdater alignedTransformUpdater;
    /// The surfaces (or detector elements?) to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// The alignment mask at each iteration
    std::map<unsigned int, std::bitset<6>> iterationState;
    /// Maximum number of iterations
    size_t maxNumIterations = 100;
        /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<size_t, double> deltaChi2ONdfCutOff = {10, 0.00001};
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff = 0.10;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TelescopeAlignmentAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the alignment algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  FW::ProcessCode execute(const FW::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
