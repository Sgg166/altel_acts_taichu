// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TelescopeTrackFindingAlgorithm.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Telescope{
  template <typename TrackFinder>
  struct TrackFinderFunctionImpl {
    TrackFinder trackFinder;
  
    TrackFinderFunctionImpl(TrackFinder&& f) : trackFinder(std::move(f)) {}

    Telescope::TelescopeTrackFindingAlgorithm::TrackFinderResult operator()(
									    const std::vector<PixelSourceLink>& sourceLinks,
									    const FW::TrackParameters& initialParameters,
									    const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>& options) const {
      return trackFinder.findTracks(sourceLinks, initialParameters, options);
    }
  };
}  // namespace

Telescope::TelescopeTrackFindingAlgorithm::TrackFinderFunction
Telescope::TelescopeTrackFindingAlgorithm::makeTrackFinderFunction(
								   std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
								   FW::Options::BFieldVariant magneticField, Acts::Logging::Level lvl) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding track finder.
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
