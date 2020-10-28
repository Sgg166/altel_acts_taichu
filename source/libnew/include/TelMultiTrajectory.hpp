#pragma once

#include <unordered_map>
#include <utility>

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

#include "TelSourceLink.hpp"
#include "TelElement.hpp"
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
  class PixelMultiTrajectory;
  using TelMultiTrajectory = PixelMultiTrajectory;

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

  JsonValue createJsonValue(JsonAllocator& jsa, Acts::GeometryContext& gctx) const;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
  static Acts::FreeToBoundMatrix
  freeToCurvilinearJacobian(const Acts::Vector3D &direction);

private:
  // The multiTrajectory
  Acts::MultiTrajectory<PixelSourceLink> m_multiTrajectory;

  // The entry indices of trajectories stored in multiTrajectory
  std::vector<size_t> m_trackTips = {};

  // The fitted parameters at the provided surface for individual trajectories
  IndexedParams m_trackParameters = {};
};

} // namespace Telescope
