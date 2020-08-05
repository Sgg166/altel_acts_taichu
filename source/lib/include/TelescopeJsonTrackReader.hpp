#pragma once

#include <Acts/Utilities/Logger.hpp>
#include <memory>
#include <string>
#include <map>
#include <mutex>

#include "ACTFW/Framework/IReader.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include "JsonGenerator.hpp"

#include "PixelSourceLink.hpp"


namespace Telescope{

class TelescopeJsonTrackReader final : public FW::IReader {
 public:
  struct Config {
    std::string inputDataFile;
    std::string outputTracks;
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

  static bool createSourcelinksFromJSON(const Telescope::JsonValue& js_evpack,
                                        const std::map<size_t, std::shared_ptr<const Acts::Surface>>& surfaces,
                                        const Acts::BoundMatrix& cov_hit,
                                        std::vector<Telescope::PixelSourceLink>& sourcelinks);

 private:
  Config m_cfg;
  std::pair<size_t, size_t> m_eventsRange;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Telescope::JsonAllocator> m_jsa;
  std::unique_ptr<Telescope::JsonGenerator> m_jsgen;
  std::unique_ptr<Telescope::JsonDocument> m_jsdoc;
  Acts::BoundMatrix m_cov_hit;
  std::mutex m_mtx_read;

  const Acts::Logger& logger() const { return *m_logger; }
};

}
