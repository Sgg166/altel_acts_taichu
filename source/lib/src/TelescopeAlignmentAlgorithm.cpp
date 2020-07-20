// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeAlignmentAlgorithm.hpp"
#include "PixelSourceLink.hpp"

#include <stdexcept>

#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"

FW::TelescopeAlignmentAlgorithm::TelescopeAlignmentAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TelescopeAlignmentAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputFileName.empty()) {
    throw std::invalid_argument("Missing input data file");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

FW::ProcessCode FW::TelescopeAlignmentAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;
  using Measurement =
      Acts::Measurement<FW::PixelSourceLink, Acts::ParDef::eLOC_0,
                        Acts::ParDef::eLOC_1>;

  // Read input data
  std::vector<SourceLinkTrack> sourcelinkTracks =
      m_cfg.trackReader(m_cfg.inputFileName, m_cfg.maxNumTracks);

  std::cout << "There are " << sourcelinkTracks.size() << " tracks read-in"
            << std::endl;

  // Prepare the initial track parameters collection
  std::vector<Acts::CurvilinearParameters> initialParameters;
  initialParameters.reserve(sourcelinkTracks.size());
  unsigned int iTrack = 0;
  while (iTrack < sourcelinkTracks.size()) {
    // Create initial parameters
    // The position is taken from the first measurement
    const auto& sourcelinks = sourcelinkTracks.at(iTrack);
    const Acts::Vector3D global0 =
        sourcelinks.at(0).globalPosition(ctx.geoContext);
    const Acts::Vector3D global1 =
        sourcelinks.at(1).globalPosition(ctx.geoContext);
    Acts::Vector3D distance = global1 - global0;

    const double phi = Acts::VectorHelpers::phi(distance);
    const double theta = Acts::VectorHelpers::theta(distance);

    // shift along the beam by 100_mm
    Acts::Vector3D rPos = global0 - distance / 2;
    Acts::Vector3D rMom(4_GeV * sin(theta) * cos(phi),
                        4_GeV * sin(theta) * sin(phi), 4_GeV * cos(theta));

    Acts::BoundSymMatrix cov;
    cov << std::pow(50_um, 2), 0., 0., 0., 0., 0., 0., std::pow(50_um, 2), 0.,
        0., 0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 0.0001, 0., 0., 0.,
        0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;

    Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(
        cov, rPos, rMom, 1., 0);
    initialParameters.push_back(rStart);
    iTrack++;
  }

  // Prepare the output data with MultiTrajectory
  TrajectoryContainer trajectories;
  trajectories.reserve(sourcelinkTracks.size());


  auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
      Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});

  
  // Set the KalmanFitter options
  Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
      ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
      Acts::VoidOutlierFinder(), pSurface.get());

  // Set the alignment options
  AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>
      alignOptions(kfOptions, m_cfg.alignedTransformUpdater,
                   m_cfg.alignedDetElements, m_cfg.chi2ONdfCutOff,
                   m_cfg.maxNumIterations, m_cfg.iterationState);

  ACTS_DEBUG("Invoke alignment");
  auto result = m_cfg.align(sourcelinkTracks, initialParameters, alignOptions);
  if (result.ok()) {
    ACTS_VERBOSE(
        "Alignment finished with deltaChi2 = " << result.value().deltaChi2);
  } else {
    ACTS_WARNING("Alignment failed with " << result.error());
  }

  // ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
