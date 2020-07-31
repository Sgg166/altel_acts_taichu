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

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Units.hpp"
#include "PixelMultiTrajectory.hpp"
#include "PixelSourceLink.hpp"

namespace Telescope {

  class TelescopeFittingAlgorithm final : public FW::BareAlgorithm {
 public:
  using FitterResult = Acts::Result<Acts::KalmanFitterResult<PixelSourceLink>>;
  /// Fit function that takes input measurements, initial trackstate and fitter
  /// options and returns some fit-specific result.
  using FitterFunction = std::function<FitterResult(
                                                    const std::vector<PixelSourceLink>&, const FW::TrackParameters&,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>&)>;
  using SourceLinkTrack = std::vector<PixelSourceLink>;
  using SourceLinkTrackReader =
      std::function<std::vector<SourceLinkTrack>(const std::string&, size_t)>;

  /// Create the fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static FitterFunction makeFitterFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      FW::Options::BFieldVariant magneticField, Acts::Logging::Level lvl);

  struct Config {
    /// Input data.
    std::string inputSourcelinks{"sourcelinks_to_fit"};
    /// Output fitted trajectories collection.
    std::string outputTrajectories{"trajectories_from_fit"};
    /// Type erased fitter function.
    FitterFunction fit;

    double seedResX{10 * Acts::UnitConstants::um};
    double seedResY{10 * Acts::UnitConstants::um};
    double seedResPhi{0.7};
    double seedResTheta{0.7};
    double beamEnergy{5 * Acts::UnitConstants::GeV};
  };

  /// Constructor of the fitting algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TelescopeFittingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the fitting algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  FW::ProcessCode execute(const FW::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
