#include "TelescopeJsonTrackReader.hpp"

#include <Acts/Utilities/Units.hpp>

#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Utilities/Paths.hpp"


#include "PixelSourceLink.hpp"

#include "myrapidjson.h"

using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

Telescope::TelescopeJsonTrackReader::TelescopeJsonTrackReader(
    const Telescope::TelescopeJsonTrackReader::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(0, 0),
      m_logger(Acts::getDefaultLogger("TelescopeJsonTrackReader", lvl))
{

  std::string datafile_name = m_cfg.inputDataFile;
  size_t n_datapack_select_opt = m_cfg.maxSelectDatapackNum;

  m_jsa.reset(new rapidjson::CrtAllocator);
  m_js_selected_datapack_col.reset(new JsonValue(rapidjson::kArrayType) );
  {
    std::FILE* fp = std::fopen(datafile_name.c_str(), "r");
    if(!fp) {
      std::fprintf(stderr, "File opening failed: %s \n", datafile_name.c_str());
      throw std::system_error(EIO, std::generic_category(), "File opening failed");;
    }

    char readBuffer[UINT16_MAX];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc(m_jsa.get());
    doc.ParseStream(is);
    std::fclose(fp);

    if(!doc.IsArray() || !doc.GetArray().Size()){
      std::fprintf(stderr, "no, it is not data array\n");
      throw std::system_error(EDOM, std::generic_category(), "File is not valid json array");;
    }

    uint64_t processed_datapack_count = 0;
    for(auto datapack_it = doc.Begin(); datapack_it != doc.End() && m_js_selected_datapack_col->GetArray().Size() < n_datapack_select_opt ; datapack_it++ ){
      auto &evpack = *datapack_it;
      processed_datapack_count ++;
      auto &frames = evpack["layers"];
      bool is_good_datapack = true;
      for(auto& layer : evpack["layers"].GetArray()){
        uint64_t l_hit_n = layer["hit"].GetArray().Size();
        if(l_hit_n != 1){
          is_good_datapack = false;
          continue;
        }
      }
      if(!is_good_datapack){
        continue;
      }
      m_js_selected_datapack_col->PushBack(evpack, *m_jsa); // moved
    }
  }

  if(cfg.oneEventMode){
    m_eventsRange = {0, 1};
  }
  else{
    m_eventsRange = {0, m_js_selected_datapack_col->GetArray().Size()};
  }
}

std::string Telescope::TelescopeJsonTrackReader::name() const {
  return "TelescopeJsonTrackReader";
}

std::pair<size_t, size_t> Telescope::TelescopeJsonTrackReader::availableEvents() const {
  return m_eventsRange;
}

FW::ProcessCode Telescope::TelescopeJsonTrackReader::read(const FW::AlgorithmContext& ctx) {

  Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
  cov_hit(0, 0) =m_cfg.resX*m_cfg.resX;
  cov_hit(1, 1) =m_cfg.resY*m_cfg.resY;

  if(m_cfg.oneEventMode){
    std::vector<std::vector<Telescope::PixelSourceLink>> sourcelinkTracks;
    sourcelinkTracks.reserve(m_js_selected_datapack_col->Size());
    for(auto datapack_it = m_js_selected_datapack_col->Begin(); datapack_it != m_js_selected_datapack_col->End() ; datapack_it++ ){
      std::vector<Telescope::PixelSourceLink> sourcelinks;
      const auto &datapack = *datapack_it;
      const auto &layers = datapack["layers"];
      {
        for(const auto& [i, s] : m_cfg.surfaces){
          double x_hit = layers[i]["hit"][0]["pos"][0].GetDouble() - 0.02924*1024/2.0;
          double y_hit = layers[i]["hit"][0]["pos"][1].GetDouble() - 0.02688*512/2.0;
          Acts::Vector2D loc_hit;
          loc_hit << x_hit, y_hit;
          sourcelinks.emplace_back(*s, loc_hit, cov_hit);
        }
      }
      sourcelinkTracks.push_back(std::move(sourcelinks));
    }
    ctx.eventStore.add(m_cfg.outputTracks, std::move(sourcelinkTracks));
  }
  else{
    std::vector<Telescope::PixelSourceLink> sourcelinks;
    const auto &datapack = (*m_js_selected_datapack_col)[ctx.eventNumber];
    const auto &layers = datapack["layers"];
    {
      for(const auto& [i, s] : m_cfg.surfaces){
        double x_hit = layers[i]["hit"][0]["pos"][0].GetDouble() - 0.02924*1024/2.0;
        double y_hit = layers[i]["hit"][0]["pos"][1].GetDouble() - 0.02688*512/2.0;
        Acts::Vector2D loc_hit;
        loc_hit << x_hit, y_hit;
        sourcelinks.emplace_back(*s, loc_hit, cov_hit);
      }
    }
    ctx.eventStore.add(m_cfg.outputTracks, std::move(sourcelinks));
  }
  return FW::ProcessCode::SUCCESS;
}
