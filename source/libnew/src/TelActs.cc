#include "TelActs.hh"

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

std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
          std::vector<std::shared_ptr<TelActs::TelescopeDetectorElement>>>
TelActs::buildGeometry(Acts::GeometryContext &nominal_gctx,
                       const JsonValue& js) {

  std::shared_ptr<const Acts::TrackingGeometry> geo_world;
  std::vector<std::shared_ptr<TelActs::TelescopeDetectorElement>> element_col;

  std::vector<std::shared_ptr<const Acts::Surface>> surface_col;
  std::vector<Acts::LayerPtr> layer_col;

  if (!js.HasMember("geometry")) {
    throw;
  }
  const auto &js_geo = js["geometry"];
  const auto &js_dets = js_geo["detectors"];
  std::map<size_t, std::array<double, 6>> id_geo_map;
  for(const auto& js_det: js_dets.GetArray()){
    size_t id = js_det["id"].GetUint();
    double cx = js_det["center"]["x"].GetDouble();
    double cy = js_det["center"]["y"].GetDouble();
    double cz = js_det["center"]["z"].GetDouble();
    double rx = js_det["rotation"]["x"].GetDouble();
    double ry = js_det["rotation"]["y"].GetDouble();
    double rz = js_det["rotation"]["z"].GetDouble();
    double ptx = js_det["pitch"]["x"].GetDouble();
    double pty = js_det["pitch"]["y"].GetDouble();
    double ptz = js_det["pitch"]["z"].GetDouble();
    double px = js_det["pixel"]["x"].GetDouble();
    double py = js_det["pixel"]["y"].GetDouble();
    double pz = js_det["pixel"]["z"].GetDouble();
    double sx = js_det["size"]["x"].GetDouble();
    double sy = js_det["size"]["y"].GetDouble();
    double sz = js_det["size"]["z"].GetDouble();

    Acts::Vector3D translation(cx, cy, cz);
    // The rotation around global z axis
    Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
    // The rotation around global y axis
    Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
    // The rotation around global x axis
    Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
    Acts::Rotation3D rotation = rotZ * rotY * rotX;
    // Acts::Transform3D
    auto trafo = std::make_shared<Acts::Transform3D>(
        Acts::Translation3D(translation) * rotation);

    // Create the detector element
    // Boundaries of the surfaces (ALPIDE real SIZE: 29.941760
    // * 13.762560_mm*mm)
    using namespace Acts::UnitLiterals;

    auto detElement = std::make_shared<TelActs::TelescopeDetectorElement>(
        id, trafo, 50_mm, 25_mm, 80_um);

    surface_col.push_back(detElement->surface().getSharedPtr());
    layer_col.push_back(detElement->layer());
    element_col.push_back(detElement);
  }

  // Build tracking volume
  float tracker_halfX = 0.1_m;
  float tracker_halfY = 0.1_m;
  float tracker_halfZ = 1.0_m;
  auto tracker_cuboid = std::make_shared<Acts::CuboidVolumeBounds>(
      tracker_halfX, tracker_halfY, tracker_halfZ);

  auto layer_array_binned = Acts::LayerArrayCreator({}).layerArray(
      nominal_gctx, layer_col, -tracker_halfZ, tracker_halfZ,
      Acts::BinningType::arbitrary, Acts::BinningValue::binZ);

  auto trackVolume = Acts::TrackingVolume::create(
      Acts::Transform3D::Identity(), tracker_cuboid, nullptr,
      std::move(layer_array_binned), nullptr, {}, "Tracker");

  // Build world volume
  // assumming single tracker in the world
  auto tracker_array_binned =
      Acts::TrackingVolumeArrayCreator({}).trackingVolumeArray(
          nominal_gctx, {trackVolume}, Acts::BinningValue::binZ);
  auto mtvpWorld = Acts::TrackingVolume::create(Acts::Transform3D::Identity(),
                                                tracker_cuboid,
                                                tracker_array_binned, "World");
  geo_world.reset(new Acts::TrackingGeometry(mtvpWorld));

  return std::make_pair(geo_world, element_col);
}
