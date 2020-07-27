// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "TelescopeDetectorElement.hpp"

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

namespace Acts {
class TrackingGeometry;
}

namespace Telescope {

/// Global method to build the telescope tracking geometry
///
/// @tparam detector_element_t is the actual type of the detector
/// element, each derivative of a TelescopeDetectorElement can be used
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param matDecorator is an optional decorator for the material
template <typename detector_element_t>
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename detector_element_t::ContextType& gctx,
    std::vector<std::shared_ptr<detector_element_t>>& detectorStore,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator = nullptr) {
  using namespace Acts::UnitLiterals;

  // Construct the rotation
  Acts::RotationMatrix3D rotation = Acts::RotationMatrix3D::Identity();
  // Euler angle around y will be -90_degree
  double rotationAngle = 90_degree;
  Acts::Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Acts::Vector3D yPos(0., 1., 0.);
  Acts::Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Boundaries of the surfaces (ALPIDE SIZE: 27.52512 * 13.76256_mm*mm)
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>(
      // Acts::RectangleBounds(27.52512_mm, 13.76256_mm));
      Acts::RectangleBounds(1000_mm, 1000_mm));

  // Material of the surfaces
  Acts::MaterialProperties matProp(95.7, 465.2, 28.03, 14., 2.32e-3, 80_um);
  const auto surfaceMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Set translation vectors
  std::vector<Acts::Vector3D> translations;
  translations.reserve(6);
  translations.push_back({-95_mm, 0., 0.});
  translations.push_back({-57_mm, 0., 0.});
  translations.push_back({-19_mm, 0., 0.});
  translations.push_back({19_mm, 0., 0.});
  translations.push_back({57_mm, 0., 0.});
  translations.push_back({95_mm, 0., 0.});

  // Construct surfaces
  std::array<std::shared_ptr<const Acts::Surface>, 6> surfaces;
  unsigned int i;
  for (i = 0; i < translations.size(); i++) {
    Acts::Transform3D trafo(Acts::Transform3D::Identity() * rotation);
    trafo.translation() = translations[i];
    // Create the detector element
    auto detElement = std::make_shared<TelescopeDetectorElement>(
        std::make_shared<const Acts::Transform3D>(trafo), rBounds, 1_um,
        surfaceMaterial);
    // And remember the surface
    surfaces[i] = detElement->surface().getSharedPtr();
    // Add it to the event store
    detectorStore.push_back(std::move(detElement));
  }

  // Construct layers
  std::array<Acts::LayerPtr, 6> layers;
  for (i = 0; i < 6; i++) {
    Acts::Transform3D trafo(Acts::Transform3D::Identity() * rotation);
    trafo.translation() = translations[i];

    std::unique_ptr<Acts::SurfaceArray> surArray(
        new Acts::SurfaceArray(surfaces[i]));

    layers[i] = Acts::PlaneLayer::create(
        std::make_shared<const Acts::Transform3D>(trafo), rBounds,
        std::move(surArray), 1.5_cm);

    auto mutableSurface = const_cast<Acts::Surface*>(surfaces[i].get());
    mutableSurface->associateLayer(*layers[i]);
  }

  // Build volume for surfaces
  Acts::Transform3D trafoVol(Acts::Transform3D::Identity());
  trafoVol.translation() = Acts::Vector3D(0., 0., 0.);

  auto boundsVol =
      std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));

  Acts::LayerVector layVec;
  for (i = 0; i < 6; i++) {
    layVec.push_back(layers[i]);
  }
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      gctx, layVec, -0.5_m - 1._mm, 0.5_m + 1._mm, Acts::BinningType::arbitrary,
      Acts::BinningValue::binX));

  auto trackVolume = Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(trafoVol), boundsVol, nullptr,
      std::move(layArr), nullptr, {}, "Tracker");

  // Build world volume
  Acts::Transform3D trafoWorld(Acts::Transform3D::Identity());
  trafoWorld.translation() = Acts::Vector3D(0., 0., 0.);

  auto worldVol =
      std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);

  std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3D>> tapVec;
  tapVec.push_back(std::make_pair(trackVolume, Acts::Vector3D(0., 0., 0.)));

  std::vector<float> binBoundaries = {-0.6_m, 0.6_m};
  Acts::BinningData binData(Acts::BinningOption::open, Acts::BinningValue::binX,
                            binBoundaries);
  std::unique_ptr<const Acts::BinUtility> bu(new Acts::BinUtility(binData));

  std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
      new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));

  Acts::MutableTrackingVolumePtr mtvpWorld(Acts::TrackingVolume::create(
      std::make_shared<const Acts::Transform3D>(trafoWorld), worldVol, trVolArr,
      "World"));

  // Build and return tracking geometry
  return std::unique_ptr<Acts::TrackingGeometry>(
      new Acts::TrackingGeometry(mtvpWorld));
}

}  // end of namespace Telescope
