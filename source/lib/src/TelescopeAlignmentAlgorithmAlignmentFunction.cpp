// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "TelescopeTrack.hpp"


#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace {
  template <typename Alignment>
  struct AlignmentFunctionImpl {
    Alignment align;

    AlignmentFunctionImpl(Alignment&& a) : align(std::move(a)) {}

    Telescope::AlignResult operator()
    (
     const std::vector<std::vector<Telescope::PixelSourceLink>>& sourceLinks,
     const std::vector<Acts::CurvilinearTrackParameters>& initialParameters,
     const ActsAlignment::AlignmentOptions< Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>& options) const
    {
      return align.align(sourceLinks, initialParameters, options);
    };
  };
}  // namespace




Telescope::AlignmentFunction Telescope::makeAlignmentFunction
(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
 ActsExamples::Options::BFieldVariant magneticField, Acts::Logging::Level lvl){
  
  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
                    [trackingGeometry, lvl](auto&& inputField) -> AlignmentFunction {
                      // each entry in the variant is already a shared_ptr
                      // need ::element_type to get the real magnetic field type
                      using InputMagneticField =
                        typename std::decay_t<decltype(inputField)>::element_type;
                      using MagneticField = Acts::SharedBField<InputMagneticField>;
                      using Stepper = Acts::EigenStepper<MagneticField>;
                      using Navigator = Acts::Navigator;
                      using Propagator = Acts::Propagator<Stepper, Navigator>;
                      using Updater = Acts::GainMatrixUpdater;
                      using Smoother = Acts::GainMatrixSmoother;
                      using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;
                      using Alignment = ActsAlignment::Alignment<Fitter>;

                      // construct all components for the fitter
                      MagneticField field(std::move(inputField));
                      Stepper stepper(std::move(field));
                      Navigator navigator(trackingGeometry);
                      navigator.resolvePassive = false;
                      navigator.resolveMaterial = true;
                      navigator.resolveSensitive = true;
                      Propagator propagator(std::move(stepper), std::move(navigator));
                      Fitter fitter(std::move(propagator));
                      Alignment alignment(std::move(fitter),
                                                      Acts::getDefaultLogger("Alignment", lvl));
                      
                      // build the alignment functions. owns the alignment object.
                      return AlignmentFunctionImpl<ActsAlignment::Alignment<Fitter>>(std::move(alignment));
                    },
                    
                    std::move(magneticField)
                    );
  
}
