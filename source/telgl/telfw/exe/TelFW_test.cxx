#include "getopt.h"

#include <string>
#include <fstream>
#include <sstream>
#include <memory>
#include <chrono>
#include <thread>
#include <regex>
#include <atomic>
#include <future>

#include <cstdio>
#include <csignal>

#include <chrono>
#include <thread>
#include <iostream>

#include "TelGL.hh"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linenoise.h"
#include "myrapidjson.h"

#include "TelFW.hh"
#include "TelFWTest.hh"



static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -data     [PATH]   path to data file (input)
  -geometry [PATH]   path to geometry file (input)
)";

static sig_atomic_t g_done = 0;
int main(int argc, char **argv){
  signal(SIGINT, [](int){g_done+=1;});

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "data",       required_argument, NULL,           'd' },
     { "geometry",   required_argument, NULL,           'g' },
     { 0, 0, 0, 0 }};

  std::string data_path;
  std::string geometry_path;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'd':
      data_path = optarg;
      break;
    case 'g':
      geometry_path = optarg;
      break;
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      std::fprintf(stderr,"case 1\n");
      exit(1);
      break;
    case ':':
      std::fprintf(stderr,"case :\n");
      exit(1);
      break;
    case '?':
      std::fprintf(stderr,"case ?\n");
      exit(1);
      break;
    default:
      std::fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if(do_help){
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    exit(0);
  }

  if(geometry_path.empty() || data_path.empty()){
    std::fprintf(stdout, "Error: input option.\n");
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    std::exit(1);
  }


  //////////////////////////////////////////

  TelFW telfw(800, 400, "test");
  TelFWTest telfwtest(geometry_path, data_path);

  telfw.startAsync<TelFWTest>(&telfwtest, &TelFWTest::beginHook, &TelFWTest::clearHook, &TelFWTest::drawHook);
  using namespace std::chrono_literals;
  while(!g_done){
    std::this_thread::sleep_for(1s);
  }
  telfw.stopAsync();
  
  return 0;
}
