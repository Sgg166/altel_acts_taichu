// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeFittingAlgorithm.hpp"

#include <stdexcept>

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"


Telescope::TelescopeFittingAlgorithm::TelescopeFittingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TelescopeFittingAlgorithm", level),
      m_cfg(std::move(cfg)) {

}

ActsExamples::ProcessCode Telescope::TelescopeFittingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // Read input data
  const auto& trackSourcelinks
    = ctx.eventStore.get<std::vector<Telescope::PixelSourceLink>>(m_cfg.inputSourcelinks);

  // Prepare the output data with MultiTrajectory
  std::vector<PixelMultiTrajectory> trajectories;

  // We can have empty tracks which must give empty fit results
  if (!trackSourcelinks.empty()) {
    // Set initial parameters for the particle track

    // @Below is what is used when the detector aligned along global x in the first beginning. 
    // But this is not working when it's aligned along global z. Needs investigation of the reason.
    //Acts::BoundSymMatrix cov;
    //cov << std::pow(10_mm, 2), 0., 0., 0., 0., 0., 0., std::pow(10_mm, 2), 0.,
    //    0., 0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 0.0001, 0., 0., 0.,
    //    0., 0., 0., 0.0001, 0., 0., 0., 0., 0., 0., 1.;
    //Acts::Vector3D rPos(-120_mm, 0, 0);
    //Acts::Vector3D rMom(4_GeV, 0, 0);
    double resX2 = m_cfg.seedResX * m_cfg.seedResX;
    double resY2 = m_cfg.seedResY * m_cfg.seedResY;
    double resPhi2 = m_cfg.seedResPhi * m_cfg.seedResPhi;
    double resTheta2 = m_cfg.seedResTheta * m_cfg.seedResTheta;

    Acts::BoundSymMatrix cov;
    cov <<
      resX2,0., 0., 0., 0., 0.,
      0., resY2, 0., 0., 0., 0.,
      0., 0., resPhi2, 0., 0., 0.,
      0., 0., 0., resTheta2, 0., 0.,
      0., 0., 0., 0., 0.0001, 0.,
      0., 0., 0., 0.,     0., 1.;

    const Acts::Vector3D global0 = trackSourcelinks.at(0).globalPosition(ctx.geoContext);
    const Acts::Vector3D global1 = trackSourcelinks.at(1).globalPosition(ctx.geoContext);
    Acts::Vector3D distance = global1 - global0;
    const double phi = Acts::VectorHelpers::phi(distance);
    const double theta = Acts::VectorHelpers::theta(distance);
    Acts::Vector3D rPos = global0 - distance / 2;
    Acts::Vector4D rPos4(rPos.x(), rPos.y(), rPos.z(), 0);
    Acts::Vector3D rDir(sin(theta) * cos(phi),
                        sin(theta) * sin(phi),
                        cos(theta));
    double q = 1;

    Acts::CurvilinearTrackParameters rStart(rPos4, rDir, q/m_cfg.beamEnergy, cov);

    // Set the KalmanFitter options
    // TODO, remove target surface, note the vector3d
    auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>
      (Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions
      (ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
       Acts::VoidOutlierFinder(), Acts::LoggerWrapper{logger()},
      Acts::PropagatorPlainOptions(), &(*pSurface));

    ACTS_DEBUG("Invoke fitter");
    auto result = m_cfg.fit(trackSourcelinks, rStart, kfOptions);
    if (result.ok()) {
      // Get the fit output object
      const auto& fitOutput = result.value();
      // The track entry indices container. One element here.
      std::vector<size_t> trackTips;
      trackTips.reserve(1);
      trackTips.emplace_back(fitOutput.trackTip);
      // The fitted parameters container. One element (at most) here.
      IndexedParams indexedParams;
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.value();
        ACTS_VERBOSE("Fitted paramemeters for track " << ctx.eventNumber);
        ACTS_VERBOSE("  " << params.parameters().transpose());
        // Push the fitted parameters to the container
        indexedParams.emplace(fitOutput.trackTip, std::move(params));
      } else {
        ACTS_DEBUG("No fitted paramemeters for track " << ctx.eventNumber);
      }
      // Create a PixelMultiTrajectory
      trajectories.emplace_back(std::move(fitOutput.fittedStates),
                                std::move(trackTips), std::move(indexedParams));
    } else {
      ACTS_WARNING("Fit failed for track " << ctx.eventNumber << " with error"
                   << result.error());
      trajectories.push_back(PixelMultiTrajectory());
    }
  }
  else{
    ACTS_WARNING("Empty event " << ctx.eventNumber << " found.");
    trajectories.push_back(PixelMultiTrajectory());
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}


namespace {
template <typename Fitter>
struct FitterFunctionImpl {
  Fitter fitter;

  FitterFunctionImpl(Fitter&& f) : fitter(std::move(f)) {}

  Telescope::TelescopeFittingAlgorithm::FitterResult operator()
  (
   const std::vector<Telescope::PixelSourceLink>& sourceLinks,
   const ActsExamples::TrackParameters& initialParameters,
   const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& options) const {
    return fitter.fit(sourceLinks, initialParameters, options);
  };
};
}  // namespace

Telescope::TelescopeFittingAlgorithm::FitterFunction
Telescope::TelescopeFittingAlgorithm::makeFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    ActsExamples::Options::BFieldVariant magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [trackingGeometry](auto&& inputField) -> FitterFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::StraightLineStepper;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

        // construct all components for the fitter
        MagneticField field(std::move(inputField));
        Stepper stepper;
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter fitter(std::move(propagator));

        // build the fitter functions. owns the fitter object.
        return FitterFunctionImpl<Fitter>(std::move(fitter));
      },
      std::move(magneticField));
}
