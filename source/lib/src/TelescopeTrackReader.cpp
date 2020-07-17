

#include "Acts/Surfaces/Surface.hpp"
#include "TelescopeTrackReader.hpp"
#include "myrapidjson.h"

std::vector<std::vector<FW::PixelSourceLink>>
TelescopeTrackReader::operator()(const std::string& fileName, size_t nTracks) const
{
  // Create a container for the output tracks
  std::vector<std::vector<FW::PixelSourceLink>> sourcelinkTracks;
  sourcelinkTracks.reserve(nTracks);

  // Read in the raw tracks
  std::vector<std::vector<PixelHit>> rawTracks =
    jsonTrackReader(fileName, nTracks);

  std::cout << "There are " << rawTracks.size()
            << " tracks from jsonTrackReader" << std::endl;
  // Loop over the raw track to create the source link track
  for (const auto& rtrack : rawTracks) {
    // The number of hits should be less or equal to number of provided
    // surfaces?
    assert(rtrack.size() <= detectorSurfaces.size());

    // Setup local covariance
    Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
    cov(0, 0) = resolution[0] * resolution[0];
    cov(1, 1) = resolution[1] * resolution[1];
    // Create the track sourcelinks
    std::vector<FW::PixelSourceLink> sourcelinks;
    sourcelinks.reserve(rtrack.size());
    for (const auto& hit : rtrack) {
      Acts::Vector2D loc;
      loc << hit.locX, hit.locY;
      // push a hit
      sourcelinks.emplace_back(*detectorSurfaces.at(hit.surfaceIndex), loc,
                               cov);
    }
    // push the sourcelinks into the trajectory container
    sourcelinkTracks.push_back(sourcelinks);
  }

  return sourcelinkTracks;
}


std::vector<std::vector<PixelHit>>
TelescopeTrackReader::jsonTrackReader(const std::string& fileName, size_t nTracks) const
{
  
  std::FILE* fp = std::fopen(fileName.c_str(), "r");
  if (!fp) {
    std::fprintf(stderr, "File %s opening failed\n", fileName.c_str());
    throw std::ios_base::failure("Could not open '" + fileName);
  }

  char readBuffer[65536];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document data;
  data.ParseStream(is);
  if (!data.IsArray()) {
    std::fprintf(stderr, "data file is not json array\n");
    exit(2);
  }

  rapidjson::Value::ConstValueIterator ev_it = data.Begin();
  rapidjson::Value::ConstValueIterator ev_it_end = data.End();

  std::vector<std::vector<PixelHit>> rawTracks;
  rawTracks.reserve(nTracks);
  size_t itrack = 0;
  while (ev_it != ev_it_end and itrack < nTracks) {
    ClusterPool cpool;

    for (auto& subev : ev_it->GetArray()) {
      for (auto& hit : subev["hit_xyz_array"].GetArray()) {
        uint16_t pixX = hit[0].GetInt();
        uint16_t pixY = hit[1].GetInt();
        uint16_t layerIndex = hit[2].GetInt();
        cpool.addHit(pixX, pixY, layerIndex);
        // track.emplace_back(PixelHit{layerIndex,
        //                             (pixX - (nPixX - 1) / 2.) * pitchX,
        //                             (pixY - (nPixY - 1) / 2.) * pitchY});
      }
    }
    cpool.buildClusters();
    ++ev_it;

    std::vector<PixelHit> track = cpool.m_ccenter_col;
    // std::cout << "track size " << track.size() << std::endl;

    std::vector<uint16_t> cluster_counter(6, 0);
    for (auto& c : track) {
      cluster_counter[c.surfaceIndex]++;
      c.locX = (c.locX - (nPixX - 1) / 2.) * pitchX;
      c.locY = (c.locY - (nPixY - 1) / 2.) * pitchY;
    }

    bool isgood_data = true;
    for (auto n : cluster_counter) {
      if (n != 1) {
        isgood_data = false;
        break;
      }
    }
    if (!isgood_data) {
      continue;
    }

    // std::fprintf(stdout, "\nadd track clusters:");
    // for (auto& h : track) {
    //  std::fprintf(stdout, " [%f, %f, %lu] ", h.locX, h.locY,
    //  h.surfaceIndex);
    //}
    // push the track to the track container
    rawTracks.push_back(track);
    itrack++;
  }
  return rawTracks;
}
