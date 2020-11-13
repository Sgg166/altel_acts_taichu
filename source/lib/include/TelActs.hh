#pragma once

#include "Acts/MagneticField/ConstantBField.hpp"

#include "Acts/TrackFinding/CKFSourceLinkSelector.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"

#include "ActsAlignment/Kernel/Alignment.hpp"

#include "TelElement.hpp"
#include "TelEvent.hpp"
#include "TelSourceLink.hpp"
#include "myrapidjson.h"

namespace TelActs{

  using AlignResult
  = Acts::Result<ActsAlignment::AlignmentResult>;
  using AlignmentFunction
  = std::function<AlignResult(const std::vector<std::vector<TelSourceLink>> &,
                              const std::vector<Acts::CurvilinearTrackParameters> &,
                              const ActsAlignment::AlignmentOptions<
                              Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>> &)>;
  AlignmentFunction makeAlignmentFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                          std::shared_ptr<Acts::ConstantBField> magneticField,
                                          Acts::Logging::Level lvl);

  using CKFOptions
  =  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector>;
  using TrackFinderResult
  = Acts::Result<Acts::CombinatorialKalmanFilterResult<TelSourceLink>>;
  using TrackFinderFunction
  = std::function<TrackFinderResult(const std::vector<TelSourceLink> &,
                                    const Acts::BoundTrackParameters &, const CKFOptions &)>;
  TrackFinderFunction makeTrackFinderFunction(std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                                              std::shared_ptr<Acts::ConstantBField> magneticField);


  Acts::FreeToBoundMatrix
  freeToCurvilinearJacobian(const Acts::Vector3D &direction);


  std::unique_ptr<TelEvent>
  createTelEvent(Acts::GeometryContext& gctx,
                 const Acts::CombinatorialKalmanFilterResult<TelActs::TelSourceLink>& ckfResult,
                 size_t runN, size_t eventN, size_t detSetupN,
                 const std::map<Acts::GeometryIdentifier, size_t>&  mapGeoId2DetId);

  void matchAddExtraHitMeas(std::shared_ptr<TelEvent> telEvent,
                            const std::vector<TelSourceLink>& sourcelinksTargets,
                            const std::map<Acts::GeometryIdentifier, size_t>&  mapGeoId2DetId);

  std::unique_ptr<TelActs::TelEvent>
  createTelEvent(const JsonValue& js, size_t runN, size_t eventN, size_t detSetupN,
                 std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer);

  std::vector<TelActs::TelSourceLink>
  createSourceLinks(std::shared_ptr<TelActs::TelEvent> telEvent,
                   std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer);

  std::pair<size_t, std::shared_ptr<Acts::PlaneLayer>>
  createPlaneLayer(const JsonValue& js_det);

  std::shared_ptr<Acts::TrackingGeometry>
  createWorld(Acts::GeometryContext &gctx, double sizex, double sizey, double sizez,
              const std::vector<std::shared_ptr<const Acts::PlaneLayer>>& planeLayers);

};
