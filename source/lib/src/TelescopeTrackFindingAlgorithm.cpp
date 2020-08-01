// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeTrackFindingAlgorithm.hpp"

#include <stdexcept>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

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
