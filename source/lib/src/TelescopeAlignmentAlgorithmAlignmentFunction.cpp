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

#include "TelescopeAlignmentAlgorithm.hpp"

#include "ACTFW/Alignment/Alignment.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryID.hpp"
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

  FW::TelescopeAlignmentAlgorithm::AlignResult operator()(
      std::vector<std::vector<FW::PixelSourceLink>>& sourceLinks,
      const std::vector<Acts::CurvilinearParameters>& initialParameters,
      const FW::AlignmentOptions<
          Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>& options) const {
    return align.align(sourceLinks, initialParameters, options);
  };
};
}  // namespace

FW::TelescopeAlignmentAlgorithm::AlignmentFunction
FW::TelescopeAlignmentAlgorithm::makeAlignmentFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    Options::BFieldVariant magneticField, Acts::Logging::Level lvl) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

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
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;
        using Alignment = FW::Alignment<Fitter>;

        // construct all components for the fitter
        MagneticField field(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter fitter(std::move(propagator),
                      Acts::getDefaultLogger("KalmanFitter", lvl));
        Alignment alignment(std::move(fitter),
                            Acts::getDefaultLogger("Alignment", lvl));

        // build the alignment functions. owns the alignment object.
        return AlignmentFunctionImpl<Alignment>(std::move(alignment));
      },
      std::move(magneticField));
}
