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

  using JsonGenerator = ::JsonGeneratorArrayUnwrap;

  static double ReadBeamEnergyFromDataFile(const std::string &datafile_name) {
    std::printf("\n------------------------------------\n");
    double beamEnergy;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_data(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    while (1) {
      if (read_datapack_count > read_datapack_max) {
        throw;
      }

      doc_data.Populate(gen);
      if (!gen.isvalid) {
        throw;
      }

      const auto &js_evpack = doc_data;
      if (js_evpack.HasMember("testbeam")) {
        const auto &js_testbeam = js_evpack["testbeam"];
        beamEnergy = js_testbeam["energy"].GetDouble();
        std::printf("beam energy  %f \n", beamEnergy);
        break;
      }
    }
    return beamEnergy;
  }

  static std::map<size_t, std::array<double, 6>>
  ReadGeoFromDataFile(const std::string &datafile_name) {
    std::printf("\n------------------------------------\n");
    std::map<size_t, std::array<double, 6>> id_geo_map;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_data(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    while (1) {
      if (read_datapack_count > read_datapack_max) {
        throw;
      }

      doc_data.Populate(gen);
      if (!gen.isvalid) {
        throw;
      }

      const auto &js_evpack = doc_data;
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
    std::map<size_t, std::array<double, 6>> id_geo_map;

    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_geo(&jsa);
    Telescope::JsonGenerator gen(geofile_name);
    doc_geo.Populate(gen);
    if (!gen.isvalid) {
      throw;
    }

    if (!doc_geo.HasMember("geo")) {
      throw;
    }

    const auto &js_geo = doc_geo["geo"];
    for (const auto &l : js_geo.GetArray()) {
      size_t id = l["id"].GetUint();
      double cx = l["centerX"].GetDouble();
      double cy = l["centerY"].GetDouble();
      double cz = l["centerZ"].GetDouble();
      double rx = l["rotX"].GetDouble();
      double ry = l["rotY"].GetDouble();
      double rz = l["rotZ"].GetDouble();
      id_geo_map[id] = {cx, cy, cz, rx, ry, rz};
    }
    return id_geo_map;
  }

  static void
  WriteGeoToGeoFile(const std::string &geofile_name,
                    std::map<size_t, std::array<double, 6>> id_geo_map) {
    Telescope::JsonAllocator jsa;
    Telescope::JsonValue js_output(rapidjson::kObjectType);
    Telescope::JsonValue js_geo(rapidjson::kArrayType);
    for (const auto &id_geo : id_geo_map) {
      Telescope::JsonValue js_ele(rapidjson::kObjectType);
      js_ele.AddMember("id", JsonValue(id_geo.first), jsa);
      js_ele.AddMember("centerX", JsonValue((id_geo.second)[0]), jsa);
      js_ele.AddMember("centerY", JsonValue((id_geo.second)[1]), jsa);
      js_ele.AddMember("centerZ", JsonValue((id_geo.second)[2]), jsa);
      js_ele.AddMember("rotX", JsonValue((id_geo.second)[3]), jsa);
      js_ele.AddMember("rotY", JsonValue((id_geo.second)[4]), jsa);
      js_ele.AddMember("rotZ", JsonValue((id_geo.second)[5]), jsa);
      js_geo.PushBack(std::move(js_ele), jsa);
    }

    js_output.AddMember("geo", std::move(js_geo), jsa);
    std::FILE *fp = std::fopen(geofile_name.c_str(), "w");
    char writeBuffer[UINT16_MAX];
    rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
    rapidjson::PrettyWriter<rapidjson::FileWriteStream> writer(os);
    js_output.Accept(writer);
    std::fclose(fp);
  }

} // namespace Telescope
