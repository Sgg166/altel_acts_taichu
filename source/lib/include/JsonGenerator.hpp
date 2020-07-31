#pragma once

#include <ratio>
#include <chrono>
#include <filesystem>
#include <cstdio>

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
