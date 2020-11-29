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
template <typename TrackFinder> struct TrackFinderFunctionImpl {
  TrackFinder trackFinder;

  TrackFinderFunctionImpl(TrackFinder &&f) : trackFinder(std::move(f)) {}

  TelActs::TrackFinderResult operator()(
      const std::vector<TelActs::PixelSourceLink> &sourceLinks,
      const Acts::BoundTrackParameters &initialParameters,
      const Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>
          &options) const {
    return trackFinder.findTracks(sourceLinks, initialParameters, options);
  }
};
} // namespace

TelActs::TrackFinderFunction
TelActs::makeTrackFinderFunction(
                                 std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                 std::shared_ptr<Acts::ConstantBField> magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // each entry in the variant is already a shared_ptr
  // need ::element_type to get the real magnetic field type
  using InputMagneticField =
    typename std::decay_t<decltype(magneticField)>::element_type;
  using MagneticField = Acts::SharedBField<InputMagneticField>;
  using Stepper = Acts::StraightLineStepper;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using CKF =
    Acts::CombinatorialKalmanFilter<Propagator, Updater, Smoother,
                                    Acts::CKFSourceLinkSelector>;

  // construct all components for the track finder
  MagneticField field(std::move(magneticField));
  Stepper stepper;
  Navigator navigator(trackingGeometry);
  navigator.resolvePassive = false;
  navigator.resolveMaterial = true;
  navigator.resolveSensitive = true;
  Propagator propagator(std::move(stepper), std::move(navigator));
  CKF trackFinder(std::move(propagator));

  // build the track finder functions. owns the track finder object.
  return TrackFinderFunctionImpl<CKF>(std::move(trackFinder));

}
