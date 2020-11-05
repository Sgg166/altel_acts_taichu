
#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Geometry/GeometryObjectSorter.hpp"

#include "myrapidjson.h"

namespace Acts {

class Surface;
class PlanarBounds;
class ISurfaceMaterial;
} // namespace Acts


namespace TelActs {

class TelElement;
using TelescopeDetectorElement = TelElement;

class TelElement : public Acts::DetectorElementBase {
public:
  // TelElement(size_t telDetID,
  //            std::shared_ptr<Acts::Transform3D> transform,
  //            double widthX, double heightY, double thickZ);

  TelElement(const JsonValue& js);

  ~TelElement() override;

  /// Return surface associated with this detector element
  const Acts::Surface &surface() const override{ return *m_surface;}

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override {  return m_elementThickness;}

  /// Return local to global transform associated with this identifier
  const Acts::Transform3D &
  transform(const Acts::GeometryContext &gctx) const override;

  /// Return the nominal local to global transform
  const Acts::Transform3D &
  nominalTransform(const Acts::GeometryContext &gctx) const;
  void addAlignedTransform(std::unique_ptr<Acts::Transform3D> alignedTransform);

  std::shared_ptr<const Acts::Layer> layer() const { return m_layer;}
  size_t id() const { return m_tel_det_id; }


  static std::shared_ptr<Acts::TrackingGeometry>
  buildWorld(Acts::GeometryContext &gctx, double sizex, double sizey, double sizez,
             std::vector<std::shared_ptr<TelActs::TelElement>> element_col);

  static std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
                   std::vector<std::shared_ptr<TelescopeDetectorElement>>>
  buildGeometry(Acts::GeometryContext &nominal_gctx, const JsonValue &js);

private:
  std::shared_ptr<Acts::Transform3D> m_elementTransform;
  std::shared_ptr<Acts::Surface> m_surface;
  double m_elementThickness;

  std::unique_ptr<Acts::Transform3D> m_alignedTransforms;
  std::shared_ptr<Acts::Layer> m_layer;
  size_t m_tel_det_id;

};


inline const Acts::Transform3D &
TelElement::transform(const Acts::GeometryContext &gctx) const {
  // Check if a different transform than the nominal exists
  if (m_alignedTransforms) {
    return (*m_alignedTransforms);
  }
  // Return the standard transform if not found
  return nominalTransform(gctx);
}

inline const Acts::Transform3D &TelElement::nominalTransform(
    const Acts::GeometryContext & /*gctx*/) const {
  return *m_elementTransform;
}

inline void TelElement::addAlignedTransform(
    std::unique_ptr<Acts::Transform3D> alignedTransform) {
  m_alignedTransforms = std::move(alignedTransform);
}

} // namespace Telescope
