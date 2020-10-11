#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "ActsAlignment/Kernel/Alignment.hpp"
#include "PixelSourceLink.hpp"

#include "Acts/TrackFitting/KalmanFitter.hpp"


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

#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

#include "TelescopeDetectorElement.hpp"

#include "myrapidjson.h"

namespace Telescope{
  using AlignResult = Acts::Result<ActsAlignment::AlignmentResult>;

  using AlignmentFunction = std::function<
    AlignResult(const std::vector<std::vector<PixelSourceLink>>&,
                const std::vector<Acts::CurvilinearTrackParameters>&,
                const ActsAlignment::AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>&)
    >;

  AlignmentFunction makeAlignmentFunction (std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                           ActsExamples::Options::BFieldVariant magneticField, Acts::Logging::Level lvl);




  void BuildGeometry(
                     Acts::GeometryContext& nominal_gctx,
                     std::shared_ptr<const Acts::TrackingGeometry>& geo_world,
                     std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>>& element_col,
                     const std::map<size_t, std::array<double, 6>>& opts,
                     double widthX, double heightY, double thickZ
                     );

}
