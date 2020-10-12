#pragma once

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <map>
#include <ratio>

#include "myrapidjson.h"

namespace Telescope {

  using JsonAllocator = ::JsonAllocator;
  using JsonValue = ::JsonValue;
  using JsonDocument = ::JsonDocument;
  using JsonReader = ::JsonReader;
  using JsonGenerator = ::JsonFileDeserializer;


  static double ReadBeamEnergyFromDataFile(const std::string &datafile_name) {
    std::printf("\n------------------------------------\n");
    double beamEnergy = -1;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    JsonFileDeserializer jsf(datafile_name);
    while (jsf) {
      read_datapack_count ++;
      if (read_datapack_count > read_datapack_max) {
        throw;
      }
      JsonDocument jsd= jsf.getNextJsonDocument();
      const auto &js_evpack = jsd;
      if (js_evpack.HasMember("testbeam")) {
        const auto &js_testbeam = js_evpack["testbeam"];
        beamEnergy = js_testbeam["energy"].GetDouble();
        std::printf("beam energy  %f \n", beamEnergy);
        break;
      }
    }
    if(beamEnergy < -1){
      beamEnergy = 5;
      std::printf("unable to extract beam energy, set to default  %f \n", beamEnergy);
    }
    return beamEnergy;
  }

  static std::map<size_t, std::array<double, 6>>
  ReadGeoFromDataFile(const std::string &datafile_name) {
    std::printf("\n------------------------------------\n");
    std::map<size_t, std::array<double, 6>> id_geo_map;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;

    JsonFileDeserializer jsf(datafile_name);
    while (jsf) {
      read_datapack_count ++;
      if (read_datapack_count > read_datapack_max) {
        throw;
      }
      JsonDocument jsd= jsf.getNextJsonDocument();
      const auto &js_evpack = jsd;
      if (js_evpack.HasMember("telescope")) {
        const auto &js_telescope = js_evpack["telescope"];
        size_t ln = 0;
        for (const auto &l : js_telescope["locations"].GetObject()) {
          std::string name = l.name.GetString();
          double loc = l.value.GetDouble();
          std::printf("layer %s  %lu,  locationZ %f \n", name.c_str(), ln, loc);
          id_geo_map[ln] = {0, 0, loc, 0, 0, 0};
          ln++;
        }
        break;
      }
    }
    return id_geo_map;
  }


  static std::map<size_t, std::array<double, 6>>
  ReadGeoFromGeoFile(const std::string &geofile_name) {
    JsonFileDeserializer jsf(geofile_name);
    JsonDocument doc_geo = jsf.getNextJsonDocument();
    if (!doc_geo.HasMember("geometry")) {
      throw;
    }
    const auto &js_geo = doc_geo["geometry"];
    const auto &js_dets = js_geo["detectors"];
    std::map<size_t, std::array<double, 6>> id_geo_map;
    for(const auto& js_det: js_dets.GetArray()){
      size_t id = js_det["id"].GetUint();
      double cx = js_det["center"]["x"].GetDouble();
      double cy = js_det["center"]["y"].GetDouble();
      double cz = js_det["center"]["z"].GetDouble();
      double rx = js_det["rotation"]["x"].GetDouble();
      double ry = js_det["rotation"]["y"].GetDouble();
      double rz = js_det["rotation"]["z"].GetDouble();
      double ptx = js_det["pitch"]["x"].GetDouble();
      double pty = js_det["pitch"]["y"].GetDouble();
      double ptz = js_det["pitch"]["z"].GetDouble();
      double px = js_det["pixel"]["x"].GetDouble();
      double py = js_det["pixel"]["y"].GetDouble();
      double pz = js_det["pixel"]["z"].GetDouble();
      double sx = js_det["size"]["x"].GetDouble();
      double sy = js_det["size"]["y"].GetDouble();
      double sz = js_det["size"]["z"].GetDouble();
      id_geo_map[id] = {cx, cy, cz, rx, ry, rz};
    }
    return id_geo_map;
  }

  static void
  WriteGeoToGeoFile(const std::string &geofile_name,
                    std::map<size_t, std::array<double, 6>> id_geo_map) {
    JsonDocument jsdoc;
    jsdoc.SetObject();
    JsonAllocator& jsa= jsdoc.GetAllocator();
    using rapidjson::kObjectType;
    using rapidjson::kArrayType;
    JsonValue js_geo(kObjectType);
    JsonValue js_dets(kArrayType);
    for (const auto &id_geo : id_geo_map) {
      JsonValue js_det(kObjectType);
      js_det.AddMember("id", id_geo.first, jsa);

      JsonValue js_size(kObjectType);
      js_size.AddMember("x", 0.02924*1024., jsa);
      js_size.AddMember("y", 0.02688*512., jsa);
      js_size.AddMember("z", 1., jsa);
      js_det.AddMember("size", std::move(js_size), jsa);

      JsonValue js_pitch(kObjectType);
      js_pitch.AddMember("x", 0.02924, jsa);
      js_pitch.AddMember("y", 0.02688, jsa);
      js_pitch.AddMember("z", 1., jsa);
      js_det.AddMember("pitch", std::move(js_pitch), jsa);

      JsonValue js_pixel(kObjectType);
      js_pixel.AddMember("x", 1024, jsa);
      js_pixel.AddMember("y", 512, jsa);
      js_pixel.AddMember("z", 1, jsa);
      js_det.AddMember("pixel", std::move(js_pixel), jsa);

      // JsonValue js_color(kObjectType);
      // js_color.AddMember("r", colors[i%ncolor][0], jsa);
      // js_color.AddMember("g", colors[i%ncolor][1], jsa);
      // js_color.AddMember("b", colors[i%ncolor][2], jsa);
      // js_det.AddMember("color", std::move(js_color), jsa);

      JsonValue js_center(kObjectType);
      js_center.AddMember("x", (id_geo.second)[0], jsa);
      js_center.AddMember("y", (id_geo.second)[1], jsa);
      js_center.AddMember("z", (id_geo.second)[2], jsa);
      js_det.AddMember("center", std::move(js_center), jsa);

      JsonValue js_rotation(kObjectType);
      js_rotation.AddMember("x", (id_geo.second)[3], jsa);
      js_rotation.AddMember("y", (id_geo.second)[4], jsa);
      js_rotation.AddMember("z", (id_geo.second)[5], jsa);
      js_det.AddMember("rotation", std::move(js_rotation), jsa);

      js_dets.PushBack(js_det, jsa);
    }
    js_geo.AddMember("detectors", std::move(js_dets), jsa);
    jsdoc.AddMember("geometry", std::move(js_geo), jsa);

    std::string jsstr = JsonUtils::stringJsonValue(jsdoc, true);

    std::FILE *fp = std::fopen(geofile_name.c_str(), "w");
    std::fwrite(jsstr.data(), 1, jsstr.size(), fp);
    std::fclose(fp);
  }

} // namespace Telescope
