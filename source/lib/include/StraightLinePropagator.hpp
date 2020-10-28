
#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace TelActs {

struct StraightLinePropagator {
  using Stepper = Acts::StraightLineStepper;
  using State = Acts::StraightLineStepper::State;

  /// Nested propagator options struct
  struct PropagatorOptions {
    /// Propagation direction
    Acts::NavigationDirection direction = Acts::NavigationDirection::forward;

    /// The mass of electron
    double mass = 0.5109989461 * Acts::UnitConstants::MeV;

    /// Absolute maximum step size
    double maxStepSize = std::numeric_limits<double>::max();

    /// Absolute maximum path length
    double pathLimit = std::numeric_limits<double>::max();

    /// Required tolerance to reach target (surface, pathlength)
    double targetTolerance = 10 * Acts::UnitConstants::um;

    /// Maxium number of target trials
    size_t maxTargetTrials = 2;
  };

  ///
  ///@brief Constructor using an Acts StraightLine stepper
  ///
  StraightLinePropagator(Stepper &&stepper) : m_stepper(std::move(stepper)) {}

  ///
  ///@brief transport (bound) track parameters to target surface
  ///
  template <typename parameters_t>
  Acts::BoundTrackParameters
  transport(const Acts::GeometryContext &gctx,
            const Acts::MagneticFieldContext &mctx,
            const PropagatorOptions &options, const parameters_t &par,
            const Acts::Surface &targetSurface) const {
    // Construct a stepping state with bound/curvilinear track parameter
    State stepping(gctx, mctx, par, options.direction, options.maxStepSize,
                   options.targetTolerance);
    // Try to target the surface by a few trials
    bool targetReached = false;
    for (unsigned int i = 0; i < options.maxTargetTrials; i++) {
      if (target(gctx, stepping, options, targetSurface)) {
        targetReached = true;
        break;
      }
    }
    if (not targetReached) {
      throw std::runtime_error(
          "Target surface is not reached with provided starting parameters and "
          "navigation direction. Good luck.");
    }
    // get the bound parameters at target surface
    auto targetState = m_stepper.boundState(stepping, targetSurface);
    auto targetBoundParams = std::get<Acts::BoundTrackParameters>(targetState);
    return targetBoundParams;
  }

  ///
  ///@brief target at the provided target surface
  ///
  bool target(const Acts::GeometryContext &gctx, State &stepping,
              const PropagatorOptions &options,
              const Acts::Surface &targetSurface) const {
    const auto sIntersection = targetSurface.intersect(
        gctx, m_stepper.position(stepping),
        stepping.navDir * m_stepper.direction(stepping), true);
    // The target is reached or not
    bool targetReached = (sIntersection.intersection.status ==
                          Acts::Intersection3D::Status::onSurface);
    double distance = sIntersection.intersection.pathLength;

    if (not targetReached) {
      // Target is not reached, update the step size
      const double overstepLimit = m_stepper.overstepLimit(stepping);
      // Check the alternative solution
      if (distance < overstepLimit and sIntersection.alternative) {
        // Update the distance to the alternative solution
        distance = sIntersection.alternative.pathLength;
      }
      // Update the step size
      stepping.stepSize.update(stepping.navDir * distance,
                               Acts::ConstrainedStep::aborter);

      // Update the position/momentum, derivative and jacTransport
      step(stepping, options);
    }
    return targetReached;
  }

  ///
  ///@brief update the stepping state using the step size
  ///
  double step(State &stepping, const PropagatorOptions &options) const {
    // use the adjusted step size
    const auto h = stepping.stepSize;
    // time propagates along distance as 1/b = sqrt(1 + m²/p²)
    const auto dtds = std::hypot(1., options.mass / stepping.p);
    // Update the track parameters according to the equations of motion
    stepping.pos += h * stepping.dir;
    stepping.t += h * dtds;
    // Propagate the jacobian
    if (stepping.covTransport) {
      // The step transport matrix in global coordinates
      Acts::FreeMatrix D = Acts::FreeMatrix::Identity();
      D.block<3, 3>(0, 4) = Acts::ActsSymMatrixD<3>::Identity() * h;
      // Extend the calculation by the time propagation
      // Evaluate dt/dlambda
      D(3, 7) = h * options.mass * options.mass *
                (stepping.q == 0. ? 1. : stepping.q) / (stepping.p * dtds);
      // Set the derivative factor the time
      stepping.derivative(3) = dtds;
      // Update jacobian and derivative
      stepping.jacTransport = D * stepping.jacTransport;
      stepping.derivative.template head<3>() = stepping.dir;
    }
    // state the path length
    stepping.pathAccumulated += h;

    // return h
    return h;
  }

private:
  Stepper m_stepper;
};

} // namespace Telescope
