#include "TelSourceLink.hpp"
#include "Acts/Utilities/Units.hpp"
#include "myrapidjson.h"
#include <random>
using namespace Acts::UnitLiterals;
 
namespace {
   struct ResolutionConfig {
    double u_size1 = 6_um, u_size2 = 6_um, u_sizeOther = 6_um;
    double v_size1 = 6_um, v_size2 = 6_um, v_sizeOther = 6_um;
    
    static ResolutionConfig& getInstance() {
      static ResolutionConfig instance;
      return instance;
    }
    
    void loadFromFile(const std::string& filename) {
      std::string content = JsonUtils::readFile(filename);
      JsonDocument doc = JsonUtils::createJsonDocument(content);
      
      if (doc.HasMember("resolutions")) {
        const auto& res = doc["resolutions"];
        if (res.HasMember("u_direction")) {
          const auto& u = res["u_direction"];
        //  if (u.HasMember("u")) u = std::stod(u["u"].GetString()) * 1_um;
          if (u.HasMember("size_1")) u_size1 = std::stod(u["size_1"].GetString()) * 1_um;
          if (u.HasMember("size_2")) u_size2 = std::stod(u["size_2"].GetString()) * 1_um;
          if (u.HasMember("size_other")) u_sizeOther = std::stod(u["size_other"].GetString()) * 1_um;
        }
        if (res.HasMember("v_direction")) {
          const auto& v = res["v_direction"];
        //  if (v.HasMember("v")) v = std::stod(v["v"].GetString()) * 1_um;
          if (v.HasMember("size_1")) v_size1 = std::stod(v["size_1"].GetString()) * 1_um;
          if (v.HasMember("size_2")) v_size2 = std::stod(v["size_2"].GetString()) * 1_um;
          if (v.HasMember("size_other")) v_sizeOther = std::stod(v["size_other"].GetString()) * 1_um;
        }
      }
    }
  };
  //Calculate the resolution corresponding to the cluster extension
  double calculateResolution(size_t clusterSize, bool isUdirection) {
    auto& config = ResolutionConfig::getInstance();
    if (isUdirection) {
   //   return config.u;
      if (clusterSize == 1) return config.u_size1;
      else if (clusterSize == 2) return config.u_size2;
      else return config.u_sizeOther;
    } else {
     // return config.v;
      if (clusterSize == 1) return config.v_size1;
      else if (clusterSize == 2) return config.v_size2;
      else return config.v_sizeOther;
    }
  }
 
  // Calculate the expansion of the cluster in the U and V directions and set the covariance matrix
  void setCovarianceFromClusterSize(Acts::BoundMatrix& cov, 
                                   std::shared_ptr<altel::TelMeasHit> hitMeas) {
    // Statistically analyze the expansion of the cluster in the U and V directions
    std::set<uint16_t> unique_u, unique_v;
    for (const auto& raw : hitMeas->measRaws()) {
      unique_u.insert(raw.u());
      unique_v.insert(raw.v());
    }
    
    size_t cluSizeU = unique_u.size();
    size_t cluSizeV = unique_v.size();
    
    double resU = calculateResolution(cluSizeU, true);
    double resV = calculateResolution(cluSizeV, false);
    
    cov(0, 0) = resU * resU;
    cov(1, 1) = resV * resV;
  }
 
  // Verify the effectiveness of hitMeas and set the detector ID
  void validateAndInitializeHitMeas(std::shared_ptr<altel::TelMeasHit> hitMeas, 
                                   size_t& detId) {
    if (!hitMeas) {
      std::fprintf(stderr, "very wrong\n");
      throw;
    }
    detId = hitMeas->DN;
  }
}
void TelActs::initializeResolutionConfig(const std::string& configFile) {
  if (!configFile.empty()) {
    ResolutionConfig::getInstance().loadFromFile(configFile);
  }
}
 
TelActs::TelSourceLink::TelSourceLink(const Acts::PlaneLayer &planeLayer, 
                                     std::shared_ptr<altel::TelMeasHit> hitMeas)
  : m_hitMeas(hitMeas), m_cov(Acts::BoundMatrix::Zero()), m_surface(&planeLayer) {
  validateAndInitializeHitMeas(hitMeas, m_detId);
  m_values << m_hitMeas->PLs[0], m_hitMeas->PLs[1];
  setCovarianceFromClusterSize(m_cov, hitMeas);
}
 
TelActs::TelSourceLink::TelSourceLink(std::shared_ptr<altel::TelMeasHit> hitMeas,
                                      const std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>>& mapDetId2PlaneLayer)
  : m_hitMeas(hitMeas), m_cov(Acts::BoundMatrix::Zero()) {
  validateAndInitializeHitMeas(hitMeas, m_detId);
  
  auto it = mapDetId2PlaneLayer.find(m_detId);
  if(it == mapDetId2PlaneLayer.end()){
    std::fprintf(stderr, "very wrong\n");
    throw;
  }
  m_surface = it->second.get();
  
  m_values << hitMeas->PLs[0], hitMeas->PLs[1];
  setCovarianceFromClusterSize(m_cov, hitMeas);
}
