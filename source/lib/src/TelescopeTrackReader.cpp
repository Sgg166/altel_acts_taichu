

#include "Acts/Surfaces/Surface.hpp"
#include "TelescopeTrackReader.hpp"
#include "myrapidjson.h"

std::vector<std::vector<FW::PixelSourceLink>>
TelescopeTrackReader::operator()(const std::string& fileName, size_t nTracks) const
{
  // Create a container for the output tracks
  std::vector<std::vector<FW::PixelSourceLink>> sourcelinkTracks;
  sourcelinkTracks.reserve(nTracks);
  
  std::FILE* fp = std::fopen(fileName.c_str(), "r");
  if(!fp) {
    std::fprintf(stderr, "File opening failed: %s \n", fileName.c_str());
    throw;
  }

  char readBuffer[1000000];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document doc;
  doc.ParseStream(is);
  if(!doc.IsArray()){
    std::fprintf(stderr, "no, it is not data array\n");
    throw;
  }
  rapidjson::Value::ConstValueIterator ev_it = doc.Begin();
  rapidjson::Value::ConstValueIterator ev_it_end = doc.End();  
  if(ev_it == ev_it_end){
    std::fprintf(stderr, "empty array\n");
    throw;
  }
  
  uint64_t processed_datapack_count = 0;
  uint64_t good_datapack_count = 0;
  while(ev_it != ev_it_end && good_datapack_count < nTracks){
    const auto &evpack = *ev_it;
    ev_it++;
    processed_datapack_count ++;
    
    const auto &frames = evpack["layers"];
    
    bool is_good_datapack = true;
    for(const auto& layer : evpack["layers"].GetArray()){
      uint64_t l_hit_n = layer["hit"].GetArray().Size();
      if(l_hit_n != 1){
        is_good_datapack = false;
        continue;
      }
    }
    if(!is_good_datapack){
      continue;
    }    
    good_datapack_count ++;
    

    // The number of hits should be less or equal to number of provided
    // surfaces?
    // std::cout<< "frames.Size()" << frames.Size()<<std::endl;
    // assert(frames.Size() <= detectorSurfaces.size());
    
    // Setup local covariance
    Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
    cov(0, 0) = resolution[0] * resolution[0];
    cov(1, 1) = resolution[1] * resolution[1];
    // Create the track sourcelinks
    std::vector<FW::PixelSourceLink> sourcelinks; 
    
    std::printf("%u ", good_datapack_count);
    for(size_t i= 0; i< 6; i++){
      double x = frames[i]["hit"][0]["pos"][0].GetDouble();
      double y = frames[i]["hit"][0]["pos"][1].GetDouble();
      Acts::Vector2D loc;
      loc << x, y;
      std::printf("<%f, %f, %u> ", x, y, i);
      sourcelinks.emplace_back(*detectorSurfaces.at(i), loc, cov);
    }
    std::printf("\n");
    
    sourcelinkTracks.push_back(sourcelinks);
  }

  return sourcelinkTracks;
}
