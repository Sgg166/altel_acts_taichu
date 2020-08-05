#pragma once

#include <ratio>
#include <chrono>
#include <filesystem>
#include <cstdio>
#include <map>

#include "myrapidjson.h"

namespace Telescope{

using JsonAllocator = rapidjson::CrtAllocator;
using JsonStackAllocator = rapidjson::CrtAllocator;
using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;
using JsonDocument = rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator, rapidjson::CrtAllocator>;
using JsonReader = rapidjson::GenericReader<rapidjson::UTF8<char>, rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

struct JsonGenerator{

  JsonGenerator(const std::filesystem::path &filepath)
  {

    std::fprintf(stdout, "input file  %s\n", filepath.c_str());
    std::filesystem::file_status st_file = std::filesystem::status(filepath);
    if(!std::filesystem::exists(st_file)){
      std::fprintf(stderr, "File < %s > does not exist.\n", filepath.c_str());
      throw;
    }
    if(!std::filesystem::is_regular_file(st_file)){
      std::fprintf(stderr, "File < %s > is not regular file.\n", filepath.c_str());
      throw;
    }
    filesize = std::filesystem::file_size(filepath);

    fp = std::fopen(filepath.c_str(), "r");
    if(!fp) {
      std::fprintf(stderr, "File opening failed: %s \n", filepath.c_str());
      throw;
    }
    is.reset(new rapidjson::FileReadStream(fp, readBuffer, sizeof(readBuffer)));
    rapidjson::SkipWhitespace(*is);
    if(is->Peek() == '['){
      is->Take();
      isarray = true;
    }
    else{
      isarray = false;
    }
    isvalid = true;
  }

  ~JsonGenerator(){
    if(fp){
      std::fclose(fp);
    }
  }

  bool operator()(JsonDocument& doc){
    reader.Parse<rapidjson::kParseStopWhenDoneFlag>(*is, doc);
    if(reader.HasParseError()) {
      if(is->Tell() + 10 < filesize){
        std::fprintf(stderr, "rapidjson error<%s> when parsing input data at offset %llu\n",
                     rapidjson::GetParseError_En(reader.GetParseErrorCode()), reader.GetErrorOffset());
        throw;
      }//otherwise, it reaches almost end of file. mute the errer message
      isvalid = false;
      return false;
    }
    if(isarray){
      if(is->Peek()==','){
        is->Take();
      }
      else{
        rapidjson::SkipWhitespace(*is);
        if(is->Peek()==','){
          is->Take();
        }
      }
    }
    isvalid = true;
    return true;
  }

  size_t filesize;
  std::FILE* fp;
  char readBuffer[UINT16_MAX+1];
  std::unique_ptr<rapidjson::FileReadStream> is;
  JsonReader reader;
  bool isvalid;
  bool isarray;

  static void example(const std::string& datafile_name){
    JsonAllocator s_jsa;
    JsonDocument doc(&s_jsa);
    JsonGenerator gen(datafile_name);
    int n = 0;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    while(1){// doc is cleared at beginning of each loop
      doc.Populate(gen);
      if(!gen.isvalid){
        break;
      }
      n++;
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    double time_sec_total = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
    double time_sec_per = time_sec_total / n;
    std::printf("datapack %llu, time_sec_total %f, time_sec_per %f,  data_req %f \n",  n, time_sec_total , time_sec_per, 1./time_sec_per);
  }


  static double ReadBeamEnergy(const std::string& datafile_name){
    std::printf("\n------------------------------------\n");
    double beamEnergy;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_data(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    while(1){
      if(read_datapack_count > read_datapack_max ){
        throw;
      }

      doc_data.Populate(gen);
      if(!gen.isvalid  ){
        throw;
      }

      const auto &js_evpack = doc_data;
      if(js_evpack.HasMember("testbeam")){
        const auto &js_testbeam = js_evpack["testbeam"];
        beamEnergy = js_testbeam["energy"].GetDouble();
        std::printf("beam energy  %f \n", beamEnergy );
        break;
      }
    }
    return beamEnergy;
  }

  static std::map<size_t, std::array<double, 6>>  ReadGeoFromDataFile(const std::string& datafile_name){
    std::printf("\n------------------------------------\n");
    std::map<size_t, std::array<double, 6>> id_geo_map;
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_data(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    while(1){
      if(read_datapack_count > read_datapack_max ){
        throw;
      }

      doc_data.Populate(gen);
      if(!gen.isvalid  ){
        throw;
      }

      const auto &js_evpack = doc_data;
      if(js_evpack.HasMember("telescope")){
        const auto &js_telescope = js_evpack["telescope"];
        size_t ln = 0;
        for(const auto& l: js_telescope["locations"].GetObject()){
          std::string name = l.name.GetString();
          double loc = l.value.GetDouble();
          std::printf("layer %s  %lu,  locationZ %f \n", name.c_str(), ln, loc);
          id_geo_map[ln]= {0, 0, loc, 0, 0, 0};
          ln ++;
        }
        break;
      }
    }
    return id_geo_map;
  }

  static std::map<size_t, std::array<double, 6>>  ReadGeoFromGeoFile(const std::string& geofile_name){
    std::map<size_t, std::array<double, 6>> id_geo_map;

    Telescope::JsonAllocator jsa;
    Telescope::JsonDocument doc_geo(&jsa);
    Telescope::JsonGenerator gen(geofile_name);
    doc_geo.Populate(gen);
    if(!gen.isvalid  ){
      throw;
    }

    if(!doc_geo.HasMember("geo")){
      throw;
    }

    const auto &js_geo = doc_geo["geo"];
    for(const auto& l: js_geo.GetArray()){
      size_t id = l["id"].GetUint();
      double cx = l["centerX"].GetDouble();
      double cy = l["centerY"].GetDouble();
      double cz = l["centerZ"].GetDouble();
      double rx = l["rotX"].GetDouble();
      double ry = l["rotY"].GetDouble();
      double rz = l["rotZ"].GetDouble();
      id_geo_map[id]= {cx, cy, cz, rx, ry, rz};
    }
    return id_geo_map;
  }

  static void  WriteGeoToGeoFile(const std::string& geofile_name, std::map<size_t, std::array<double, 6>> id_geo_map){
    Telescope::JsonAllocator jsa;
    Telescope::JsonValue js_output(rapidjson::kObjectType);
    Telescope::JsonValue js_geo(rapidjson::kArrayType);
    for (const auto& id_geo : id_geo_map) {
      Telescope::JsonValue js_ele(rapidjson::kObjectType);
      js_ele.AddMember("id",         JsonValue(id_geo.first), jsa);
      js_ele.AddMember("centerX",    JsonValue((id_geo.second)[0]), jsa);
      js_ele.AddMember("centerY",    JsonValue((id_geo.second)[1]), jsa);
      js_ele.AddMember("centerZ",    JsonValue((id_geo.second)[2]), jsa);
      js_ele.AddMember("rotX",       JsonValue((id_geo.second)[3]), jsa);
      js_ele.AddMember("rotY",       JsonValue((id_geo.second)[4]), jsa);
      js_ele.AddMember("rotZ",       JsonValue((id_geo.second)[5]), jsa);
      js_geo.PushBack(std::move(js_ele), jsa);
    }

    js_output.AddMember("geo", std::move(js_geo), jsa);
    std::FILE* fp = std::fopen(geofile_name.c_str(), "w");
    char writeBuffer[UINT16_MAX];
    rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));
    rapidjson::PrettyWriter< rapidjson::FileWriteStream> writer(os);
    js_output.Accept(writer);
    std::fclose(fp);
  }


  template<typename T>
  static void PrintJson(const T& o){
    rapidjson::StringBuffer sb;
    // rapidjson::PrettyWriter<rapidjson::StringBuffer> w(sb);
    rapidjson::Writer<rapidjson::StringBuffer> w(sb);
    o.Accept(w);
    std::fwrite(sb.GetString(), 1, sb.GetSize(), stdout);
    std::fputc('\n', stdout);
  }
};

}
