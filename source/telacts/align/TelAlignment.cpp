#include "TelActs.hh"


#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"


#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"



namespace {
template <typename Alignment> struct AlignmentFunctionImpl {
  Alignment align;
  AlignmentFunctionImpl(Alignment &&a) : align(std::move(a)) {}
  TelActs::AlignResult operator()
  (
   const std::vector<std::vector<TelActs::PixelSourceLink>> &sourceLinks,
   const std::vector<Acts::CurvilinearTrackParameters> &initialParameters,
   const ActsAlignment::AlignmentOptions<
   Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>> &options) const {
    return align.align(sourceLinks, initialParameters, options);
  };
};
} // namespace


TelActs::AlignmentFunction TelActs::makeAlignmentFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<Acts::ConstantBField> magneticField,
    Acts::Logging::Level lvl) {

  using InputMagneticField =
    typename std::decay_t<decltype(magneticField)>::element_type;
  using MagneticField = Acts::SharedBField<InputMagneticField>;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;
  using Alignment = ActsAlignment::Alignment<Fitter>;

  // construct all components for the fitter
  MagneticField field(std::move(magneticField));
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
}
