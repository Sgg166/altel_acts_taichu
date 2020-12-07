#include "TelEventTTreeReader.hpp"
#include "getopt.h"
#include "myrapidjson.h"

#include <numeric>
#include <chrono>

#include <TFile.h>
#include <TTree.h>

#include "TelGL.hh"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linenoise.h"
#include "myrapidjson.h"

#include "TelFW.hh"
#include "glfw_test.hh"


static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -wait                             wait for user keyboard input per event
  -eventSkip      <INT>             number of events to skip before start processing
  -eventMax       <INT>             max number of events to process
  -geometryFile   <PATH>            path to geometry input file (input)
  -rootFile       <PATH>            path to root file (input)

examples:
 ./bin/altelTelEventViewer -w -geo calice_geo_align4.json  -r  detresid.root
)";


int main(int argc, char *argv[]) {
  int64_t eventMaxNum = 0;
  int64_t eventSkipNum = 0;
  std::string geometryFilePath;
  std::string rootFilePath;

  int do_wait = 0;
  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"wait", no_argument, NULL, 'w'},
                                {"eventSkip", required_argument, NULL, 's'},
                                {"eventMax", required_argument, NULL, 'm'},
                                {"rootFile", required_argument, NULL, 'b'},
                                {"geometryFile", required_argument, NULL, 'g'},
                                {0, 0, 0, 0}};

    if(argc == 1){
      std::fprintf(stderr, "%s\n", help_usage.c_str());
      std::exit(1);
    }
    int c;
    int longindex;
    opterr = 1;
    while ((c = getopt_long_only(argc, argv, "-", longopts, &longindex)) != -1) {
      // // "-" prevents non-option argv
      // if(!optopt && c!=0 && c!=1 && c!=':' && c!='?'){ //for debug
      //   std::fprintf(stdout, "opt:%s,\targ:%s\n", longopts[longindex].name, optarg);;
      // }
      switch (c) {
      case 's':
        eventSkipNum = std::stoul(optarg);
        break;
      case 'm':
        eventMaxNum = std::stoul(optarg);
        break;
      case 'b':
        rootFilePath = optarg;
        break;
      case 'g':
        geometryFilePath = optarg;
        break;
      case 'w':
        do_wait=1;
        break;
        // help and verbose
      case 'v':
        do_verbose=1;
        //option is set to no_argument
        if(optind < argc && *argv[optind] != '-'){
          do_verbose = std::stoul(argv[optind]);
          optind++;
        }
        break;
      case 'h':
        std::fprintf(stdout, "%s\n", help_usage.c_str());
        std::exit(0);
        break;
        /////generic part below///////////
      case 0:
        // getopt returns 0 for not-NULL flag option, just keep going
        break;
      case 1:
        // If the first character of optstring is '-', then each nonoption
        // argv-element is handled as if it were the argument of an option
        // with character code 1.
        std::fprintf(stderr, "%s: unexpected non-option argument %s\n",
                     argv[0], optarg);
        std::exit(1);
        break;
      case ':':
        // If getopt() encounters an option with a missing argument, then
        // the return value depends on the first character in optstring:
        // if it is ':', then ':' is returned; otherwise '?' is returned.
        std::fprintf(stderr, "%s: missing argument for option %s\n",
                     argv[0], longopts[longindex].name);
        std::exit(1);
        break;
      case '?':
        // Internal error message is set to print when opterr is nonzero (default)
        std::exit(1);
        break;
      default:
        std::fprintf(stderr, "%s: missing getopt branch %c for option %s\n",
                     argv[0], c, longopts[longindex].name);
        std::exit(1);
        break;
      }
    }
  }/////////getopt end////////////////

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "geometryFile:  %s\n", geometryFilePath.c_str());
  std::fprintf(stdout, "rootFile:      %s\n", rootFilePath.c_str());
  std::fprintf(stdout, "\n");

  //////////// geometry
  std::string str_geo = JsonUtils::readFile(geometryFilePath);
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFilePath.c_str() );
    throw;
  }

  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }
  altel::TelEventTTreeReader ttreeReader;

  TFile tfile(rootFilePath.c_str(),"READ");
  if(!tfile.IsOpen()){
    std::fprintf(stderr, "tfile is not open\n");
    throw;
  }

  TTree *pTree = 0;
  tfile.GetObject("eventTree",pTree);
  if(!pTree){
    std::fprintf(stderr, "pTree is invalid\n");
    throw;
  }
  ttreeReader.setTTree(pTree);
  TelFW telfw(800, 400, "test");
  glfw_test telfwtest(geometryFilePath);
  telfw.startAsync<glfw_test>(&telfwtest, &glfw_test::beginHook, &glfw_test::clearHook, &glfw_test::drawHook);

  auto tp_start = std::chrono::system_clock::now();

  size_t totalNumEvents  = ttreeReader.numEvents();
  for(size_t eventNum = eventSkipNum; eventNum<totalNumEvents; eventNum++){
    if(eventNum> eventMaxNum && eventMaxNum>0){
      break;
    }

    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    std::fprintf(stdout, "FileEvent #%d, event #%d, clock/trigger #%d\n", eventNum, telEvent->eveN(), telEvent->clkN());
    telfwtest.pushBufferEvent(telEvent);

    if(do_wait){
      std::fprintf(stdout, "waiting, press any key to next event\n");
      std::getc(stdin);
    }
  }

  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;

  tfile.Close();
  if(do_wait){
    std::fprintf(stdout, "waiting, press any key to exit\n");
    std::getc(stdin);
  }

  telfw.stopAsync();
  return 0;
}
