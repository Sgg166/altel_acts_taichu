#pragma once

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

#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
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

#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "TelescopeDetectorElement.hpp"
#include "PixelSourceLink.hpp"
#include "myrapidjson.h"

namespace TelActs{

  using AlignResult
  = Acts::Result<ActsAlignment::AlignmentResult>;
  using AlignmentFunction
  = std::function<AlignResult(const std::vector<std::vector<PixelSourceLink>> &,
                              const std::vector<Acts::CurvilinearTrackParameters> &,
                              const ActsAlignment::AlignmentOptions<
                              Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>> &)>;
  AlignmentFunction makeAlignmentFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                          std::shared_ptr<Acts::ConstantBField> magneticField,
                                          Acts::Logging::Level lvl);

  using CKFOptions
  =  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>;
  using TrackFinderResult
  = Acts::Result<Acts::CombinatorialKalmanFilterResult<PixelSourceLink>>;
  using TrackFinderFunction
  = std::function<TrackFinderResult(const std::vector<PixelSourceLink> &,
                                    const Acts::BoundTrackParameters &, const CKFOptions &)>;
  TrackFinderFunction makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                              std::shared_ptr<Acts::ConstantBField> magneticField);

  std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
            std::vector<std::shared_ptr<TelescopeDetectorElement>>>
  buildGeometry(Acts::GeometryContext &nominal_gctx, const JsonValue &js);

};
