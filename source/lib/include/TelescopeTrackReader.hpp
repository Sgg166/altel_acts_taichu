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
/// Struct for 2D pixel hit
///
struct PixelHit {
  // The correspoinding hit index
  size_t surfaceIndex = 0;

  // The local x coordinate
  double locX = 0;

  // The local y cooridnate
  double locY = 0;
};

class ClusterPool {
 public:
  void addHit(uint64_t x, uint64_t y, uint64_t z) {
    uint64_t index = x + (y << 16) + (z << 32);
    m_hit_col.push_back(index);
  }

  void buildClusters() {
    std::vector<uint64_t> hit_col_remain = m_hit_col;

    while (!hit_col_remain.empty()) {
      std::vector<uint64_t> hit_col_this_cluster;
      std::vector<uint64_t> hit_col_this_cluster_edge;

      // get first edge seed hit
      // from un-identifed hit to edge hit
      hit_col_this_cluster_edge.push_back(hit_col_remain[0]);
      hit_col_remain.erase(hit_col_remain.begin());

      while (!hit_col_this_cluster_edge.empty()) {
        uint64_t e = hit_col_this_cluster_edge[0];
        uint64_t c = 0x00000001;  // LSB column  x
        uint64_t r = 0x00010000;  // LSB row     y

        //  8 sorround hits search,
        std::vector<uint64_t> sorround_col{e - c + r, e + r,    e + c + r,
                                           e - c,     e + c,    e - c - r,
                                           e - r,     e + c - r};

        for (auto& sr : sorround_col) {
          // only search on un-identifed hits
          auto sr_found_it =
              std::find(hit_col_remain.begin(), hit_col_remain.end(), sr);
          if (sr_found_it != hit_col_remain.end()) {
            // move the found sorround hit
            // from un-identifed hit to an edge hit
            hit_col_this_cluster_edge.push_back(sr);
            hit_col_remain.erase(sr_found_it);
          }
        }

        // after sorround search
        // move from edge hit to cluster hit
        hit_col_this_cluster.push_back(e);
        hit_col_this_cluster_edge.erase(hit_col_this_cluster_edge.begin());
      }

      double cx = 0;
      double cy = 0;
      uint64_t cz = 0;
      for (auto& hit : hit_col_this_cluster) {
        cx += (hit & 0xffff);
        cy += (hit & 0xffff0000) >> 16;
        cz = (hit & 0xffff00000000) >> 32;
      }
      cx /= hit_col_this_cluster.size();
      cy /= hit_col_this_cluster.size();

      m_ccenter_col.push_back(PixelHit{cz, cx, cy});
      m_cluster_col.push_back(std::move(hit_col_this_cluster));
    }
  }

  std::vector<uint64_t> m_hit_col;
  std::vector<std::vector<uint64_t>> m_cluster_col;
  std::vector<PixelHit> m_ccenter_col;
};


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
 private:
  /// Function to read and create raw tracks from a json file
  ///
  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  ///
  /// @return The created raw tracks
  std::vector<std::vector<PixelHit>>
  jsonTrackReader(const std::string& fileName, size_t nTracks) const;
};
