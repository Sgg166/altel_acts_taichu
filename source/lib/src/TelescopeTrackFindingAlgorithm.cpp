// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeTrackFindingAlgorithm.hpp"

#include <stdexcept>

#include <boost/program_options.hpp>
#include <iostream>
#include <map>
#include <random>

#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"



#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"



Telescope::TelescopeTrackFindingAlgorithm::TelescopeTrackFindingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TelescopeTrackFindingAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputFileName.empty()) {
    throw std::invalid_argument("Missing input data file");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

FW::ProcessCode Telescope::TelescopeTrackFindingAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // Read input data
  // Note: one entry here is actually all the source links in one event
  const std::vector<SourceLinkTrack> sourcelinkTracks =
      m_cfg.trackReader(m_cfg.inputFileName, m_cfg.maxNumTracks);

  std::cout << "There are " << sourcelinkTracks.size() << " events read-in"
            << std::endl;

  // Prepare the output data with MultiTrajectory
  std::vector<PixelMultiTrajectory> trajectories;
  trajectories.reserve(sourcelinkTracks.size());

  // Construct a plane surface centered around (0., 0., 0) and has a normal
  // vector (1., 0., 0.) as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});

  // Loop over the tracks
  for (std::size_t ievent = 0; ievent < sourcelinkTracks.size(); ++ievent) {
    // The list of hits and the initial start parameters
    const auto& trackSourcelinks = sourcelinkTracks[ievent];

    // We can have empty tracks which must give empty fit results
    if (trackSourcelinks.empty()) {
      trajectories.push_back(PixelMultiTrajectory());
      ACTS_WARNING("Empty event " << ievent << " found.");
      continue;
    }

    // @Todo: create seeds and loop over CKF for different seeds
    // Set initial parameters for the particle track
    Acts::BoundSymMatrix cov;
    cov << std::pow(10_mm, 2), 0., 0., 0., 0., 0., 0., std::pow(10_mm, 2), 0.,
        0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 0., 0., 0.1, 0., 0., 0.,
        0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;

    Acts::Vector3D rPos(-120_mm, 0, 0);
    Acts::Vector3D rMom(4_GeV, 0, 0);
    Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(
        cov, rPos, rMom, 1., 0);

    // Set the CombinatorialKalmanFilter options
    TelescopeTrackFindingAlgorithm::CKFOptions ckfOptions(
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        m_cfg.sourcelinkSelectorCfg, &(*pSurface));
    
    ACTS_DEBUG("Invoke CKF");
    auto result = m_cfg.findTracks(trackSourcelinks, rStart, ckfOptions);
     if (result.ok()) {
      // Get the track finding output object
      const auto& trackFindingOutput = result.value();
      // Create a PixelMultiTrajectory
      trajectories.emplace_back(std::move(trackFindingOutput.fittedStates),
                                std::move(trackFindingOutput.trackTips),
                                std::move(trackFindingOutput.fittedParameters));
    } else {
      ACTS_WARNING("Track finding failed for truth seed "
                   << ievent << " with error" << result.error());
      // Track finding failed, but still create an empty PixelMultiTrajectory
      trajectories.push_back(PixelMultiTrajectory());
    }
 
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}


namespace {
  template <typename TrackFinder>
  struct TrackFinderFunctionImpl {
    TrackFinder trackFinder;

    TrackFinderFunctionImpl(TrackFinder&& f) : trackFinder(std::move(f)) {}

    Telescope::TelescopeTrackFindingAlgorithm::TrackFinderResult operator()
    (
     const std::vector<Telescope::PixelSourceLink>& sourceLinks,
     const FW::TrackParameters& initialParameters,
     const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>&
     options) const {
      return trackFinder.findTracks(sourceLinks, initialParameters, options);
    }
  };
}

Telescope::TelescopeTrackFindingAlgorithm::TrackFinderFunction
Telescope::TelescopeTrackFindingAlgorithm::makeTrackFinderFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    FW::Options::BFieldVariant magneticField, Acts::Logging::Level lvl) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding track
  // finder.
  return std::visit(
      [trackingGeometry, lvl](auto&& inputField) -> TrackFinderFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::StraightLineStepper;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using CKF =
            Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother,
                                            Acts::CKFSourceLinkSelector>;

        // construct all components for the track finder
        MagneticField field(std::move(inputField));
        Stepper stepper;
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        CKF trackFinder(
            std::move(propagator),
            Acts::getDefaultLogger("CombinatorialKalmanFilter", lvl));

        // build the track finder functions. owns the track finder object.
        return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));
      },
      std::move(magneticField));
}
