
#include "TelElement.hpp"


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
#include "Acts/Material/MaterialSlab.hpp"

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



#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"


#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"


using namespace Acts::UnitLiterals;

TelActs::TelElement::TelElement(const JsonValue &js_det){
    size_t id = js_det["id"].GetUint();
    double cx = js_det["center"]["x"].GetDouble();
    double cy = js_det["center"]["y"].GetDouble();
    double cz = js_det["center"]["z"].GetDouble();
    double rx = js_det["rotation"]["x"].GetDouble();
    double ry = js_det["rotation"]["y"].GetDouble();
    double rz = js_det["rotation"]["z"].GetDouble();
    double sx = js_det["size"]["x"].GetDouble();
    double sy = js_det["size"]["y"].GetDouble();
    double sz = js_det["size"]["z"].GetDouble();

    m_tel_det_id = id;

    m_elementThickness = 80_um;
    // auto rBounds = std::make_shared<Acts::RectangleBounds>(
    //   sx / 2.0 * Acts::UnitConstants::mm ,
    //   sy / 2.0 * Acts::UnitConstants::mm);

    std::shared_ptr<Acts::PlanarBounds> pBounds(new Acts::RectangleBounds(
                                                  sx * Acts::UnitConstants::mm ,
                                                  sy * Acts::UnitConstants::mm));
    //NOTE: workaround, enlarge sensor size x2 to prevent buggy decision of reaching end of tracker

    Acts::Vector3D translation(cx, cy, cz);
    Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
    Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
    Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
    Acts::Rotation3D rotation = rotZ * rotY * rotX;

    Acts::Transform3D xBeamRotation(Acts::Transform3D::Identity());
    xBeamRotation.linear()<< // y-z-x
      0, 0, 1,
      1, 0, 0,
      0, 1, 0;

    m_elementTransform = std::make_shared<Acts::Transform3D>(xBeamRotation * Acts::Translation3D(translation) * rotation);
    m_layer = Acts::PlaneLayer::create(*m_elementTransform, pBounds,
                                       nullptr, 5 * Acts::UnitConstants::mm);


    Acts::Material silicon = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, (2.329 / 28.0855) * 1_mol / 1_cm3);
    auto material = std::make_shared<Acts::HomogeneousSurfaceMaterial>(
      Acts::MaterialSlab(silicon, m_elementThickness));

    m_layer->surfaceRepresentation().assignSurfaceMaterial(material);

    m_surface = m_layer->surfaceRepresentation().getSharedPtr();
}


TelActs::TelElement::~TelElement(){

}



std::shared_ptr<Acts::TrackingGeometry>
TelActs::TelElement::buildWorld(Acts::GeometryContext &gctx, double sizex, double sizey, double sizez,
                                std::vector<std::shared_ptr<TelActs::TelElement>> element_col){

  std::vector<Acts::LayerPtr> layer_col;
  for(auto &ele: element_col){
    layer_col.push_back(ele->layer());
  }

  auto tracker_cuboid = std::make_shared<Acts::CuboidVolumeBounds>(
      sizex/2., sizey/2.,  sizez/2.);

  Acts::LayerArrayCreator layerArrayCreater({});
  auto layer_array_binned = layerArrayCreater.layerArray(
    gctx, layer_col, sizex/-2, sizex/2,
    Acts::BinningType::arbitrary, Acts::BinningValue::binX);


  auto trackVolume = Acts::TrackingVolume::create(
      Acts::Transform3D::Identity(), tracker_cuboid, nullptr,
      std::move(layer_array_binned), nullptr, {}, "Tracker");

  // // Build world volume
  // // assumming single tracker in the world
  // auto tracker_array_binned =
  //     Acts::TrackingVolumeArrayCreator({}).trackingVolumeArray(
  //       gctx, {trackVolume}, Acts::BinningValue::binX); // by x axis


  // std::cout<< trafo<<std::endl;

  // auto mtvpWorld = Acts::TrackingVolume::create(Acts::Transform3D::Identity(),
  //                                               tracker_cuboid,
  //                                               tracker_array_binned, "World");
  std::shared_ptr<Acts::TrackingGeometry> geo_world;
  geo_world.reset(new Acts::TrackingGeometry(trackVolume));

  return geo_world;
}


std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
          std::vector<std::shared_ptr<TelActs::TelescopeDetectorElement>>>
TelActs::TelElement::buildGeometry(Acts::GeometryContext &nominal_gctx,
                                   const JsonValue& js) {
  std::vector<std::shared_ptr<TelActs::TelescopeDetectorElement>> element_col;
  if (!js.HasMember("geometry")) {
    throw;
  }
  const auto &js_geo = js["geometry"];
  const auto &js_dets = js_geo["detectors"];
  for(const auto& js_det: js_dets.GetArray()){
    auto detElement = std::make_shared<TelActs::TelElement>(js_det);
    element_col.push_back(detElement);
  }


  std::shared_ptr<const Acts::TrackingGeometry> geo_world =
    TelActs::TelElement::buildWorld(nominal_gctx, 4.0_m, 0.1_m, 0.1_m,  element_col);

  return std::make_pair(geo_world, element_col);
}
