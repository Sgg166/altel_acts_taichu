#pragma once


#include <vector>
#include <array>
#include <string>
#include <algorithm>

#include "Acts/Utilities/Units.hpp"
#include "PixelSourceLink.hpp"

namespace Acts{
  class Surface;
};

using namespace Acts::UnitLiterals;


///
/// @brief Struct to read and create the source link tracks
///
struct TelescopeTrackReader {
  /// The number of pixels in local x direction
  size_t nPixX = 1024;

  /// The number of pixels in local y direction
  size_t nPixY = 512;

  /// The size of pixel pitch in local x direction
  double pitchX = 29.24_um;

  /// The size of pixel pitch in local y direction
  double pitchY = 26.88_um;

  /// The pixel detector resolution
  std::array<double, 2> resolution = {150_um, 150_um};

  /// The ordered detector surfaces
  std::vector<const Acts::Surface*> detectorSurfaces;

  /// Function to read and create source link tracks as input of fitter
  ///
  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  ///
  /// @return The created source link tracks
  std::vector<std::vector<FW::PixelSourceLink>> operator()
  (const std::string& fileName, size_t nTracks) const ;


};
