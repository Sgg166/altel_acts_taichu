#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <memory>
#include <string>
#include <map>

#include "ACTFW/Framework/IReader.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"



namespace rapidjson{
  template<typename CharType>
  struct UTF8;

  template <typename Encoding, typename Allocator>
  class GenericValue;

  class CrtAllocator;
}

namespace Telescope{

class TelescopeJsonTrackReader final : public FW::IReader {
 public:
  struct Config {
    std::string inputDataFile;
    std::string outputTracks;
    size_t oneEventMode;
    size_t maxSelectDatapackNum;
    double resX;
    double resY;
    std::map<size_t, std::shared_ptr<const Acts::Surface>> surfaces;
  };

  /// Construct the particle reader.
  ///
  /// @params cfg is the configuration object
  /// @params lvl is the logging level
  TelescopeJsonTrackReader(const Config& cfg, Acts::Logging::Level lvl);

  std::string name() const final override;

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream.
  FW::ProcessCode read(const FW::AlgorithmContext& ctx) final override;

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<rapidjson::CrtAllocator> m_jsa;
  std::unique_ptr<rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>>  m_js_selected_datapack_col;

  const Acts::Logger& logger() const { return *m_logger; }
};

}
