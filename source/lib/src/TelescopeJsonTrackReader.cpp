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

#include "myrapidjson.h"

using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

Telescope::TelescopeJsonTrackReader::TelescopeJsonTrackReader(
    const Telescope::TelescopeJsonTrackReader::Config& cfg, Acts::Logging::Level lvl)
    : m_cfg(cfg),
      m_eventsRange(0, 0),
      m_logger(Acts::getDefaultLogger("TelescopeJsonTrackReader", lvl))
{

  std::string datafile_name = m_cfg.inputDataFile;
  size_t n_datapack_select_opt = m_cfg.maxSelectEventNum;
 
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
    for(auto ev_it = doc.Begin(); ev_it != doc.End() && m_js_selected_datapack_col->GetArray().Size() < n_datapack_select_opt ; ev_it++ ){
      auto &evpack = *ev_it;
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
      m_js_selected_datapack_col->PushBack(evpack, *m_jsa);
    }
  }
  
  m_eventsRange = {0, m_js_selected_datapack_col->GetArray().Size()};
  
}

std::string Telescope::TelescopeJsonTrackReader::name() const {
  return "TelescopeJsonTrackReader";
}

std::pair<size_t, size_t> Telescope::TelescopeJsonTrackReader::availableEvents() const {

  
  return m_eventsRange;
}

FW::ProcessCode Telescope::TelescopeJsonTrackReader::read(const FW::AlgorithmContext& ctx) {

  // SimParticleContainer particles;
  // particles.adopt_sequence(std::move(unordered));
  // ctx.eventStore.add(m_cfg.outputParticles, std::move(particles));

  return FW::ProcessCode::SUCCESS;
}
