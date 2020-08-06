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

using namespace Acts::UnitLiterals;


Telescope::TelescopeTrackFindingAlgorithm::TelescopeTrackFindingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TelescopeTrackFindingAlgorithm", level),
      m_cfg(std::move(cfg)) {

}

FW::ProcessCode Telescope::TelescopeTrackFindingAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {

  const auto& sourcelinks
    = ctx.eventStore.get<std::vector<Telescope::PixelSourceLink>>(m_cfg.inputSourcelinks);

  std::vector<Telescope::PixelMultiTrajectory> trajectories;
  if (!sourcelinks.empty()) {
    // Start to find seeds for this event using the source links on the first two layers
    // @todo: add raw layer id in PixelSourceLink
    std::vector<Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy>> initialParameters;
    for(const auto& sl0 : sourcelinks){
      const auto& surface0 = sl0.referenceSurface();
      const auto& layer0 = surface0.geoID().layer();
      if(layer0 != 2) {
        continue;
      }
      const Acts::Vector3D global0 = sl0.globalPosition(ctx.geoContext);
      for(const auto& sl1 : sourcelinks){
        const auto& surface1 = sl1.referenceSurface();
        const auto& layer1 = surface1.geoID().layer();
        if(layer1 != 4) {
          continue;
        }
        const Acts::Vector3D global1 = sl1.globalPosition(ctx.geoContext);
        Acts::Vector3D distVec = global1 - global0;
        // compare their distance in x-y (r-phi) plane
        const double rDist = std::abs(Acts::VectorHelpers::perp(distVec));
        // @todo: add options for the seed cuts
        if(rDist<=3_mm) {
          Acts::BoundSymMatrix cov;
          cov <<
            m_cfg.seedResX * m_cfg.seedResX, 0., 0., 0., 0., 0.,
            0., m_cfg.seedResY * m_cfg.seedResY, 0., 0., 0., 0.,
            0., 0., m_cfg.seedResPhi*m_cfg.seedResPhi, 0., 0., 0.,
            0., 0., 0., m_cfg.seedResTheta*m_cfg.seedResTheta, 0., 0.,
            0., 0., 0., 0., 0.0001, 0.,
            0., 0., 0., 0.,     0., 1.;

          const double phi = Acts::VectorHelpers::phi(distVec);
          const double theta = Acts::VectorHelpers::theta(distVec);
          Acts::Vector3D rPos = global0 - distVec / 2;
          Acts::Vector3D rMom(m_cfg.beamEnergy * sin(theta) * cos(phi),
                              m_cfg.beamEnergy * sin(theta) * sin(phi),
                              m_cfg.beamEnergy * cos(theta));

          initialParameters.emplace_back(cov, rPos, rMom, 1., 0);
        }
      }
    }

    auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>
      (Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{0, 0., 1.});

    // Set the CombinatorialKalmanFilter options
    // @Todo: add options for CKF
    Telescope::TelescopeTrackFindingAlgorithm::CKFOptions ckfOptions
      (ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
       Acts::CKFSourceLinkSelector::Config{{Acts::GeometryID(),{500, 1}}}, refSurface.get());
    // 500 chi2cut,   1 max number of selected sourcelinks in a single surface;

    //Loop ever the seeds
    size_t iseed = 0;
    size_t nTracks = 0;
    for(const auto& rStart: initialParameters){
      auto result = m_cfg.findTracks(sourcelinks, rStart, ckfOptions);
      if (result.ok()) {
        // Get the track finding output object
        const auto& trackFindingOutput = result.value();
        // Create a PixelMultiTrajectory
        nTracks += trackFindingOutput.trackTips.size();
        trajectories.emplace_back(std::move(trackFindingOutput.fittedStates),
                                  std::move(trackFindingOutput.trackTips),
                                  std::move(trackFindingOutput.fittedParameters));
      } else {
        std::printf("Track finding failed in Event<%lu> seed<%lu>, with error \n",
                    ctx.eventNumber, iseed, result.error().message().c_str());
      }
      iseed++;
    }// end of the loop for all seeds

    std::printf("<< eventNumber: %lu    sourcelinks.size(): %lu   initialParameters.size(): %lu   nTracks %lu \n",
                ctx.eventNumber, sourcelinks.size(), initialParameters.size(), nTracks);

  }
  else{
    std::cout<<"Empty event <" << ctx.eventNumber << "> found."<<std::endl;
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
