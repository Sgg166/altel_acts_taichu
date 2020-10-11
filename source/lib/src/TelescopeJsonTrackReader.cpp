#include "TelescopeJsonTrackReader.hpp"

#include <chrono>
#include <fstream>
#include <ios>
#include <ratio>
#include <stdexcept>
#include <string>
#include <vector>

#include <Acts/Utilities/Units.hpp>

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

Telescope::TelescopeJsonTrackReader::TelescopeJsonTrackReader(
    const Telescope::TelescopeJsonTrackReader::Config &cfg,
    Acts::Logging::Level lvl)
    : m_cfg(cfg), m_eventsRange(0, SIZE_MAX),
      m_logger(Acts::getDefaultLogger("TelescopeJsonTrackReader", lvl)) {
  m_jsa.reset(new Telescope::JsonAllocator);
  m_jsgen.reset(new Telescope::JsonGenerator(m_cfg.inputDataFile));
  m_jsdoc.reset(new Telescope::JsonDocument(m_jsa.get()));

  m_cov_hit = Acts::BoundMatrix::Zero();
  m_cov_hit(0, 0) = m_cfg.resX * m_cfg.resX;
  m_cov_hit(1, 1) = m_cfg.resY * m_cfg.resY;
}

std::string Telescope::TelescopeJsonTrackReader::name() const {
  return "TelescopeJsonTrackReader";
}

std::pair<size_t, size_t>
Telescope::TelescopeJsonTrackReader::availableEvents() const {
  return m_eventsRange;
}

ActsExamples::ProcessCode Telescope::TelescopeJsonTrackReader::read(
    const ActsExamples::AlgorithmContext &ctx) {
  Telescope::JsonValue evpack;
  {
    const std::lock_guard<std::mutex> lock(m_mtx_read);
    m_jsdoc->Populate(*m_jsgen);
    if (!m_jsgen->isvalid) {
      return ActsExamples::ProcessCode::ABORT;
    }
    m_jsdoc->Swap(evpack);
  }

  std::vector<Telescope::PixelSourceLink> sourcelinks;
  bool re =
      createSourcelinksFromJSON(evpack, m_cfg.surfaces, m_cov_hit, sourcelinks);
  if (!re) {
    return ActsExamples::ProcessCode::ABORT;
  }

  ctx.eventStore.add(m_cfg.outputSourcelinks, std::move(sourcelinks));
  return ActsExamples::ProcessCode::SUCCESS;
}

bool Telescope::TelescopeJsonTrackReader::createSourcelinksFromJSON(
    const Telescope::JsonValue &js_evpack,
    const std::map<size_t, std::shared_ptr<const Acts::Surface>> &surfaces,
    const Acts::BoundMatrix &cov_hit,
    std::vector<Telescope::PixelSourceLink> &sourcelinks) {
  const auto &layers = js_evpack["layers"];
  for (const auto &layer : layers.GetArray()) {
    size_t id_ext = layer["ext"].GetUint();
    auto surface_it = surfaces.find(id_ext);
    if (surface_it == surfaces.end()) {
      continue;
    }
    for (const auto &hit : layer["hit"].GetArray()) {
      double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
      double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
      Acts::Vector2D loc_hit;
      loc_hit << x_hit, y_hit;
      sourcelinks.emplace_back(*(surface_it->second), loc_hit, cov_hit);
    }
  }
  return true;
}
