#include "TelSourceLink.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

std::vector<TelActs::TelSourceLink> TelActs::TelSourceLink::CreateSourceLinks(
  const JsonValue &js,
  const std::vector<std::shared_ptr<TelActs::TelElement>> eles){

  std::vector<TelActs::TelSourceLink> sourcelinks;
  const auto &layers = js["layers"];
  JsonAllocator jsa;
  JsonValue js_hits(rapidjson::kArrayType);
  for (const auto &layer : layers.GetArray()) {
    size_t id_ext = layer["ext"].GetUint();
    for (const auto &hit : layer["hit"].GetArray()) {
      double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
      double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
      JsonValue js_hit(rapidjson::kObjectType);
      js_hit.AddMember("id", id_ext, jsa);
      js_hit.AddMember("x", x_hit, jsa);
      js_hit.AddMember("y", y_hit, jsa);
      js_hits.PushBack(std::move(js_hit), jsa);
    }
  }

  for(const auto &js_hit : js_hits.GetArray()) {
    size_t id = js_hit["id"].GetUint();
    std::shared_ptr<TelActs::TelElement> ele;
    for(const auto& anEle: eles){
      if(anEle->id() == id){
        ele = anEle;
        break;
      }
    }
    if(ele){
      Acts::Vector2D loc_hit;
      double x = js_hit["x"].GetDouble();
      double y = js_hit["y"].GetDouble();
      loc_hit << x, y;
      //////////// hit data
      Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
      double resX = 5_um;
      double resY = 5_um;
      cov_hit(0, 0) = resX * resX;
      cov_hit(1, 1) = resY * resY;
      sourcelinks.emplace_back(ele->surface(), loc_hit, cov_hit);
    }
  }
  return sourcelinks;
}
