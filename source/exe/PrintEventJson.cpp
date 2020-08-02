#include "JsonGenerator.hpp"
#include "getopt.h"

#include <cstdio>
#include <iostream>

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -file       [jsonfile]   name of data json file
  -eventNumber [int] eventNumber to be printed

)";

static Telescope::JsonAllocator jsa;

int main(int argc, char* argv[]) {

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "file",      required_argument, NULL,           'f' },
     { "eventNumber",   required_argument, NULL,      'w' },
     { 0, 0, 0, 0 }};


  std::string datafile_name;
  size_t eventNumber = 0;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(0);
      break;
    case 'f':
      datafile_name = optarg;
      break;
    case 'w':
      eventNumber = std::stoul(optarg);
      break;

      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr,"case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr,"case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr,"case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if(datafile_name.empty()){
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "datafile:         %s\n", datafile_name.c_str());
  std::fprintf(stdout, "eventNumber:      %lu\n", eventNumber);
  std::fprintf(stdout, "\n");

  {
    uint64_t read_datapack_count = 0;
    uint64_t read_datapack_max = 20000;
    Telescope::JsonDocument doc_data(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    while(1){
      if(read_datapack_count > read_datapack_max ){
        break;
      }

      doc_data.Populate(gen);
      if(!gen.isvalid  ){
        break;
      }

      const auto &evpack = doc_data;
      if(read_datapack_count == eventNumber){
        std::printf("----------eventNuber: %lu------------ \n", read_datapack_count);
        Telescope::JsonGenerator::PrintJson(doc_data);

        size_t ln = 0;
        for(const auto&layer : evpack["layers"].GetArray()){
          std::printf("-layer %lu \n", ln);
          std::printf("--rawJson");
          Telescope::JsonGenerator::PrintJson(layer);
          std::printf("--hitSize:  %u,  hits:  ", layer["hit"].GetArray().Size());
          for(const auto&hit : layer["hit"].GetArray()){
            double hx = hit["pos"][0].GetDouble();
            double hy = hit["pos"][1].GetDouble();
            std::printf(" [%f, %f] ", hx, hy);
          }
          std::printf("\n");
          ln++;
        }
        return 0;
      }
      read_datapack_count ++;
    }
  }
  return 0;
}
