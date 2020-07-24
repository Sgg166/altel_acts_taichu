#include <array>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

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

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include "TelescopeTrack.hpp"
#include "TelescopeDetectorElement.hpp"

using namespace Acts::UnitLiterals;

void Telescope::BuildGeometry(
                              Acts::GeometryContext& nominal_gctx,
                              std::shared_ptr<const Acts::TrackingGeometry>& geo_world,
                              std::vector<std::shared_ptr<Acts::DetectorElementBase>>& element_col,
                              std::vector<std::shared_ptr<const Acts::Surface>>& surface_col,
                              std::vector<Acts::LayerPtr>& layer_col,
                              std::vector<Acts::Vector3D> translations
                              ){
  
  // Construct the rotation
  Acts::RotationMatrix3D rotation = Acts::RotationMatrix3D::Identity();
  // Euler angle around y will be -90_degree
  double rotationAngle = 90_degree;
  rotation.col(0) = Acts::Vector3D(std::cos(rotationAngle), 0., std::sin(rotationAngle));
  rotation.col(1) = Acts::Vector3D(0., 1., 0.);
  rotation.col(2) = Acts::Vector3D(-std::sin(rotationAngle), 0., std::cos(rotationAngle));
  
  // Boundaries of the surfaces (ALPIDE SIZE: 27.52512 * 13.76256_mm*mm)
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>
    (Acts::RectangleBounds(1000_mm, 1000_mm));
  // (Acts::RectangleBounds(27.52512_mm, 13.76256_mm));

  // Material of the surfaces
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>
    (Acts::MaterialProperties(95.7, 465.2, 28.03, 14., 2.32e-3, 50_um));
  
  // Set translation vectors
  for (unsigned int i = 0; i < translations.size(); i++) {
    auto trafo = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity() * rotation);
    trafo->translation() = translations[i];
    // Create the detector element
    auto detElement = std::make_shared<Telescope::TelescopeDetectorElement>(trafo, rBounds, 1_um, surfaceMaterial);
    auto detSurface = detElement->surface().getSharedPtr();
    surface_col.push_back(detSurface);
    element_col.push_back(detElement);
    
    std::unique_ptr<Acts::SurfaceArray> surArray(new Acts::SurfaceArray(detSurface));
    auto detLayer = Acts::PlaneLayer::create(trafo, rBounds, std::move(surArray), 1.5_cm);
    layer_col.push_back(detLayer);
    
    auto mutableSurface = std::const_pointer_cast<Acts::Surface>(detSurface);
    mutableSurface->associateLayer(*detLayer);
  }

  // Build tracking volume

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig);  
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(nominal_gctx,
                                                                          layer_col,
                                                                          -0.5_m - 1._mm, 0.5_m + 1._mm,
                                                                          Acts::BinningType::arbitrary,
                                                                          Acts::BinningValue::binX));
  
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);
  auto trafoVol = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  trafoVol->translation() = Acts::Vector3D(0., 0., 0.);
  auto trackVolume = Acts::TrackingVolume::create
    (trafoVol, boundsVol, nullptr, std::move(layArr), nullptr, {}, "Tracker");
  
  // Build world volume
  std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3D>> tapVec;
  tapVec.push_back(std::make_pair(trackVolume, Acts::Vector3D(0., 0., 0.)));

  Acts::BinningData binData(Acts::BinningOption::open, Acts::BinningValue::binX,
                            std::vector<float>{-0.6_m, 0.6_m});
  std::unique_ptr<const Acts::BinUtility> bu(new Acts::BinUtility(binData));

  std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
                                                            new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));

  auto trafoWorld = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  trafoWorld->translation() = Acts::Vector3D(0., 0., 0.);
  auto worldVol = std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);
  Acts::MutableTrackingVolumePtr mtvpWorld(Acts::TrackingVolume::create
                                           (trafoWorld, worldVol, trVolArr, "World"));
  
  geo_world.reset(new Acts::TrackingGeometry(mtvpWorld));
}
