#include "TelSourceLink.hpp"
#include "Acts/Utilities/Units.hpp"
#include <random>
using namespace Acts::UnitLiterals;

/*static std::random_device rd;
static std::mt19937 gen(rd());
static std::normal_distribution<double> dist_resX(6.0, 1.0);  // 均值6.0，标准差1.0
static std::normal_distribution<double> dist_resY(6.0, 1.0);  // 均值6.0，标准差1.0
*/
TelActs::TelSourceLink::TelSourceLink(const Acts::PlaneLayer &planeLayer, std::shared_ptr<altel::TelMeasHit> hitMeas)
  : m_hitMeas(hitMeas), m_cov(Acts::BoundMatrix::Zero()), m_surface(&planeLayer){
  if(!hitMeas){
    std::fprintf(stderr, "very wrong\n");
    throw;
  }
  m_values << m_hitMeas->PLs[0], m_hitMeas->PLs[1];

  double resX = 6_um;
  double resY = 6_um;

  m_cov(0, 0) = resX * resX;
  m_cov(1, 1) = resY * resY;


/*  double pitchU = 25_um;  // 像素 pitch = 25 µm 
  double pitchV = 25_um;

  // 统计 cluster 在 U 和 V 方向上的扩展
  std::set<uint16_t> unique_u, unique_v;
  for (const auto& raw : hitMeas->measRaws()) {
    unique_u.insert(raw.u());
    unique_v.insert(raw.v());
  }
  size_t cluSizeU = unique_u.size();
  size_t cluSizeV = unique_v.size();

  // 分辨率估算: 单像素 = pitch/√12, 多像素 = pitch/(√12 * √N)
  double resU = (cluSizeU > 0)
                  ? pitchU / (std::sqrt(12.0) * std::sqrt((double)cluSizeU))
                  : pitchU / std::sqrt(12.0);
  double resV = (cluSizeV > 0)
                  ? pitchV / (std::sqrt(12.0) * std::sqrt((double)cluSizeV))
                  : pitchV / std::sqrt(12.0);

  m_cov(0, 0) = resU * resU; 
  m_cov(1, 1) = resV * resV;


  double resX = dist_resX(gen) ;
  double resY = dist_resY(gen) ;
  m_cov(0, 0) = (resX*1_um) * (resX*1_um);
  m_cov(1, 1) = (resY*1_um) * (resY*1_um);
*/
}


TelActs::TelSourceLink::TelSourceLink(std::shared_ptr<altel::TelMeasHit> hitMeas,
                                      const std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer)
  :m_hitMeas(hitMeas), m_cov(Acts::BoundMatrix::Zero()){
  if(!hitMeas){
    std::fprintf(stderr, "very wrong\n");
    throw;
  }
  size_t detId = hitMeas->DN;
  auto it = mapDetId2PlaneLayer.find(detId);
  if(it==mapDetId2PlaneLayer.end()){
    std::fprintf(stderr, "very wrong\n");
    throw;
  }
  m_surface = it->second.get();

  m_values << hitMeas->PLs[0], hitMeas->PLs[1];

  double resX = 6_um;
  double resY = 6_um;
  m_cov(0, 0) = resX * resX;
  m_cov(1, 1) = resY * resY;



  /*double pitchU = 25_um;  // 像素 pitch = 25 µm 
  double pitchV = 25_um;

  // 统计 cluster 在 U 和 V 方向上的扩展
  std::set<uint16_t> unique_u;
  std::set<uint16_t> unique_v;
  for (const auto& raw : hitMeas->measRaws()) {
    unique_u.insert(raw.u());
    unique_v.insert(raw.v());
  }
  size_t cluSizeU = unique_u.size();
  size_t cluSizeV = unique_v.size();

  // 分辨率估算: 单像素 = pitch/√12, 多像素 = pitch/(√12 * √N)
  double resU = (cluSizeU > 0)
                  ? pitchU / (std::sqrt(12.0) * std::sqrt((double)cluSizeU))
                  : pitchU / std::sqrt(12.0);
  double resV = (cluSizeV > 0)
                  ? pitchV / (std::sqrt(12.0) * std::sqrt((double)cluSizeV))
                  : pitchV / std::sqrt(12.0);

  // 设置协方差矩阵 (只在局部平面位置上的 2x2 部分非零)
  m_cov(0, 0) = resU * resU;
  m_cov(1, 1) = resV * resV;

  double resX = dist_resX(gen) ;
  double resY = dist_resY(gen) ;
  m_cov(0, 0) = (resX*1_um) * (resX*1_um);
  m_cov(1, 1) = (resY*1_um) * (resY*1_um);*/
  m_hitMeas = hitMeas;

}


/*
std::vector<TelActs::TelSourceLink> TelActs::TelSourceLink::CreateSourceLinks(
  const JsonValue &js,
  const std::vector<std::shared_ptr<TelActs::TelElement>> eles)
{

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
      double resX = 6_um;
      double resY = 6_um;
      cov_hit(0, 0) = resX * resX;
      cov_hit(1, 1) = resY * resY;
      sourcelinks.emplace_back(ele->surface(), loc_hit, cov_hit);
    }
  }
  return sourcelinks;
}
*/

