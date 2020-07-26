#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "Alignment.hpp"
#include "PixelSourceLink.hpp"

#include "Acts/Fitter/KalmanFitter.hpp"


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

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"

#include "myrapidjson.h"

namespace Telescope{
  using AlignResult = Acts::Result<FW::AlignmentResult>;

  using AlignmentFunction = std::function<
    AlignResult(const std::vector<std::vector<PixelSourceLink>>&,
                const std::vector<Acts::CurvilinearParameters>&,
                const FW::AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>&)
    >;

  AlignmentFunction makeAlignmentFunction (std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                           FW::Options::BFieldVariant magneticField, Acts::Logging::Level lvl);




  void BuildGeometry(
                     Acts::GeometryContext& nominal_gctx,
                     std::shared_ptr<const Acts::TrackingGeometry>& geo_world,
                     std::vector<std::shared_ptr<Acts::DetectorElementBase>>& element_col,
                     std::vector<std::shared_ptr<const Acts::Surface>>& surface_col,
                     std::vector<Acts::LayerPtr>& layer_col,
                     const rapidjson::GenericValue<rapidjson::UTF8<>, rapidjson::CrtAllocator> &js_opt
                     );

}
