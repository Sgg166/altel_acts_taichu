// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <unordered_map>
#include <utility>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

#include "PixelSourceLink.hpp"

#include "TelescopeDetectorElement.hpp"
#include "myrapidjson.h"

namespace TelActs {
using IndexedParams = std::unordered_map<size_t, Acts::BoundTrackParameters>;

/// @brief Struct for truth track fitting/finding result with
/// Acts::KalmanFitter/Acts::CombinatorialKalmanFilter
///
/// It contains a MultiTrajectory with a vector of entry indices for individual
/// trajectories, and a map of fitted parameters indexed by the entry index.
/// In case of track fitting, there is at most one trajectory in the
/// MultiTrajectory; In case of track finding, there could be multiple
/// trajectories in the MultiTrajectory.
///
/// @note The MultiTrajectory is thought to be empty if there is no entry index
struct PixelMultiTrajectory {
public:
  /// @brief Default constructor
  ///
  PixelMultiTrajectory() = default;

  /// @brief Constructor from multiTrajectory and fitted track parameters
  ///
  /// @param multiTraj The multiTrajectory
  /// @param tTips The entry indices for trajectories in multiTrajectory
  /// @param parameters The fitted track parameters indexed by trajectory entry
  /// index
  PixelMultiTrajectory(const Acts::MultiTrajectory<PixelSourceLink> &multiTraj,
                       const std::vector<size_t> &tTips,
                       const IndexedParams &parameters)
      : m_multiTrajectory(multiTraj), m_trackTips(tTips),
        m_trackParameters(parameters) {}

  /// @brief Copy constructor
  ///
  /// @param rhs The source PixelMultiTrajectory
  PixelMultiTrajectory(const PixelMultiTrajectory &rhs)
      : m_multiTrajectory(rhs.m_multiTrajectory), m_trackTips(rhs.m_trackTips),
        m_trackParameters(rhs.m_trackParameters) {}

  /// @brief Copy move constructor
  ///
  /// @param rhs The source PixelMultiTrajectory
  PixelMultiTrajectory(PixelMultiTrajectory &&rhs)
      : m_multiTrajectory(std::move(rhs.m_multiTrajectory)),
        m_trackTips(std::move(rhs.m_trackTips)),
        m_trackParameters(std::move(rhs.m_trackParameters)) {}

  /// @brief Default destructor
  ///
  ~PixelMultiTrajectory() = default;

  /// @brief assignment operator
  ///
  /// @param rhs The source PixelMultiTrajectory
  PixelMultiTrajectory &operator=(const PixelMultiTrajectory &rhs) {
    m_multiTrajectory = rhs.m_multiTrajectory;
    m_trackTips = rhs.m_trackTips;
    m_trackParameters = rhs.m_trackParameters;
    return *this;
  }

  /// @brief assignment move operator
  ///
  /// @param rhs The source PixelMultiTrajectory
  PixelMultiTrajectory &operator=(PixelMultiTrajectory &&rhs) {
    m_multiTrajectory = std::move(rhs.m_multiTrajectory);
    m_trackTips = std::move(rhs.m_trackTips);
    m_trackParameters = std::move(rhs.m_trackParameters);
    return *this;
  }

  /// @brief Indicator if a trajectory exists
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return Whether there is trajectory with provided entry index
  bool hasTrajectory(const size_t &entryIndex) const {
    return std::count(m_trackTips.begin(), m_trackTips.end(), entryIndex) > 0;
  }

  /// @brief Indicator if there is fitted track parameters for one trajectory
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return Whether having fitted track parameters or not
  bool hasTrackParameters(const size_t &entryIndex) const {
    return m_trackParameters.count(entryIndex) > 0;
  }

  /// @brief Getter for multiTrajectory
  ///
  /// @return The multiTrajectory with trajectory entry indices
  ///
  /// @note It could return an empty multiTrajectory
  std::pair<std::vector<size_t>, Acts::MultiTrajectory<PixelSourceLink>>
  trajectory() const {
    return std::make_pair(m_trackTips, m_multiTrajectory);
  }

  /// @brief Getter of fitted track parameters for one trajectory
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return The fitted track parameters of the trajectory
  const Acts::BoundTrackParameters &
  trackParameters(const size_t &entryIndex) const {
    auto it = m_trackParameters.find(entryIndex);
    if (it != m_trackParameters.end()) {
      return it->second;
    } else {
      throw std::runtime_error(
          "No fitted track parameters for trajectory with entry index = " +
          std::to_string(entryIndex));
    }
  }

  JsonValue createJsonValue(JsonAllocator& jsa, Acts::GeometryContext& gctx) const {
    JsonValue js_tracks(rapidjson::kArrayType);
    for (const size_t &trackTip : m_trackTips) {
      // size of trackTips should <= 1 for each seed
      JsonValue js_track(rapidjson::kObjectType);
      JsonValue js_states(rapidjson::kArrayType);
      JsonValue js_states_reverse(rapidjson::kArrayType); //tmp
      m_multiTrajectory.visitBackwards(trackTip, [&](const auto &state) {
                                                   // only fill the track states with non-outlier measurement
                                                   auto typeFlags = state.typeFlags();
                                                   if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
                                                     return true;
                                                   }
                                                   if (!state.hasSmoothed()) {
                                                     return true;
                                                   }

                                                   auto telSurface = state.referenceSurface().getSharedPtr();
                                                   auto telElement = dynamic_cast<const TelActs::TelescopeDetectorElement*>(
                                                     telSurface->associatedDetectorElement());
                                                   if(!telElement){
                                                     return true;
                                                   }
                                                   size_t layerid = telElement->id();

                                                   // Get the source link info
                                                   auto meas = std::get<
                                                     Acts::Measurement<TelActs::PixelSourceLink, Acts::BoundIndices,
                                                                       Acts::eBoundLoc0, Acts::eBoundLoc1>>(
                                                                         *state.uncalibrated());
                                                   // Get local position
                                                   Acts::Vector2D pos_hit(meas.parameters()[Acts::eBoundLoc0],
                                                                            meas.parameters()[Acts::eBoundLoc1]);

                                                   // 1) Transform bound parameter to free parameter
                                                   // 1.1)Fist transform the smoothed bound parameters to free parameters
                                                   // to get the position and momentum
                                                   Acts::FreeVector freeParams =
                                                     Acts::detail::transformBoundToFreeParameters(
                                                       *telSurface, gctx, state.smoothed());
                                                   // 1.2)Get the global position, direction, q/p, t etc.
                                                   Acts::Vector3D pos(freeParams[Acts::eFreePos0],
                                                                      freeParams[Acts::eFreePos1],
                                                                      freeParams[Acts::eFreePos2]);
                                                   Acts::Vector3D dir(freeParams[Acts::eFreeDir0],
                                                                      freeParams[Acts::eFreeDir1],
                                                                      freeParams[Acts::eFreeDir2]);
                                                   double p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
                                                   int q = freeParams[Acts::eFreeQOverP]>0 ? 1:-1;
                                                   double t = freeParams[Acts::eFreeTime];

                                                   /// 1.3) Initialize the jacobian from local to the global frame
                                                   Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                                                   // Calculate the jacobian
                                                   telSurface->initJacobianToGlobal(gctx, jacToGlobal,
                                                                                    pos, dir, state.smoothed());
                                                   Acts::FreeSymMatrix freeCovariance =
                                                     jacToGlobal * state.smoothedCovariance() * jacToGlobal.transpose();

                                                   // 2) Transform free parameter to curvilinear parameter
                                                   Acts::FreeToBoundMatrix jacToCurv = freeToCurvilinearJacobian(dir);
                                                   Acts::BoundSymMatrix curvCovariance =
                                                     jacToCurv * freeCovariance * jacToCurv.transpose();

                                                   // Write out the curvilinear parameters and its covariance
                                                   JsonValue js_track_state(rapidjson::kObjectType); //
                                                   js_track_state.AddMember("id", layerid, jsa);

                                                   JsonValue js_hit(rapidjson::kObjectType);
                                                   js_hit.AddMember("x", pos_hit.x(), jsa);
                                                   js_hit.AddMember("y", pos_hit.y(), jsa);
                                                   js_track_state.AddMember("hit", std::move(js_hit), jsa);

                                                   JsonValue js_pos(rapidjson::kObjectType);
                                                   js_pos.AddMember("x", pos.x(), jsa);
                                                   js_pos.AddMember("y", pos.y(), jsa);
                                                   js_pos.AddMember("z", pos.z(), jsa);
                                                   js_track_state.AddMember("pos", std::move(js_pos), jsa);

                                                   JsonValue js_dir(rapidjson::kObjectType);
                                                   js_dir.AddMember("x", dir.x(), jsa);
                                                   js_dir.AddMember("y", dir.y(), jsa);
                                                   js_dir.AddMember("z", dir.z(), jsa);
                                                   js_track_state.AddMember("dir", std::move(js_dir), jsa);

                                                   js_track_state.AddMember("p", p, jsa);
                                                   js_track_state.AddMember("q", p, jsa);
                                                   js_track_state.AddMember("t", t, jsa);

                                                   const double *cov_data = curvCovariance.data();
                                                   JsonValue js_state_cov(rapidjson::kArrayType); //
                                                   js_state_cov.Reserve(Acts::eBoundSize * Acts::eBoundSize, jsa);
                                                   for (size_t n = 0; n < Acts::eBoundSize * Acts::eBoundSize; n++) {
                                                     js_state_cov.PushBack(JsonValue(*(cov_data + n)), jsa);
                                                   }
                                                   js_track_state.AddMember("cov", std::move(js_state_cov), jsa);
                                                   js_states_reverse.PushBack(std::move(js_track_state), jsa);
                                                   return true;
                                                 });

      for (size_t i = js_states_reverse.Size(); i > 0; i--) {
        js_states.PushBack(std::move(js_states_reverse[i - 1]), jsa);
      }
      js_track.AddMember("states", std::move(js_states), jsa);
      js_tracks.PushBack(std::move(js_track), jsa);
    }
    return js_tracks;
  }


/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
  static Acts::FreeToBoundMatrix
  freeToCurvilinearJacobian(const Acts::Vector3D &direction) {
    // Optimized trigonometry on the propagation direction
    const double x = direction(0); // == cos(phi) * sin(theta)
    const double y = direction(1); // == sin(phi) * sin(theta)
    const double z = direction(2); // == cos(theta)
    // can be turned into cosine/sine
    const double cosTheta = z;
    const double sinTheta = sqrt(x * x + y * y);
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = x * invSinTheta;
    const double sinPhi = y * invSinTheta;
    // prepare the jacobian to curvilinear
    Acts::FreeToBoundMatrix jacToCurv = Acts::FreeToBoundMatrix::Zero();
    if (std::abs(cosTheta) < Acts::s_curvilinearProjTolerance) {
      // We normally operate in curvilinear coordinates defined as follows
      jacToCurv(0, 0) = -sinPhi;
      jacToCurv(0, 1) = cosPhi;
      jacToCurv(1, 0) = -cosPhi * cosTheta;
      jacToCurv(1, 1) = -sinPhi * cosTheta;
      jacToCurv(1, 2) = sinTheta;
    } else {
      // Under grazing incidence to z, the above coordinate system definition
      // becomes numerically unstable, and we need to switch to another one
      const double c = sqrt(y * y + z * z);
      const double invC = 1. / c;
      jacToCurv(0, 1) = -z * invC;
      jacToCurv(0, 2) = y * invC;
      jacToCurv(1, 0) = c;
      jacToCurv(1, 1) = -x * y * invC;
      jacToCurv(1, 2) = -x * z * invC;
    }
    // Time parameter
    jacToCurv(5, 3) = 1.;
    // Directional and momentum parameters for curvilinear
    jacToCurv(2, 4) = -sinPhi * invSinTheta;
    jacToCurv(2, 5) = cosPhi * invSinTheta;
    jacToCurv(3, 4) = cosPhi * cosTheta;
    jacToCurv(3, 5) = sinPhi * cosTheta;
    jacToCurv(3, 6) = -sinTheta;
    jacToCurv(4, 7) = 1.;

    return jacToCurv;
  }


private:
  // The multiTrajectory
  Acts::MultiTrajectory<PixelSourceLink> m_multiTrajectory;

  // The entry indices of trajectories stored in multiTrajectory
  std::vector<size_t> m_trackTips = {};

  // The fitted parameters at the provided surface for individual trajectories
  IndexedParams m_trackParameters = {};
};

} // namespace Telescope
