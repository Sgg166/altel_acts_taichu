#include "getopt.h"

#include <numeric>
#include <chrono>

#include "TelGL.hh"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linenoise.h"
#include "myrapidjson.h"

#include "TelFW.hh"
#include "glfw_test.hh"

#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

static const std::string default_geometry = R"(
{"geometry": {"detectors": [
    {"id": 0, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 0.0}, "rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 1, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 40.0}, "rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 2, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1}, "center": {"x": 0.0,"y": 0.0,"z":80.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 3, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 120.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 4, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 160.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 5, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 200.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}}
]}}
)";


static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -wait                             wait for user keyboard input per event
  -eventSkip      <INT>             number of events to skip before start processing
  -eventMax       <INT>             max number of events to process
  -geometryFile   <PATH>            path to geometry input file (input)
  -eudaqRawFile   <PATH>            path to eudaq raw file (input)

examples:
 ./bin/altelEudaqRawViewer -w -geo calice_geo_align4.json  -eudaq run00001.raw
)";

int main(int argc, char *argv[]) {
  int64_t eventMaxNum = 0;
  int64_t eventSkipNum = 0;
  std::string geometryFilePath;
  std::string eudaqRawFilePath;

  int do_wait = 0;
  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"wait", no_argument, NULL, 'w'},
                                {"eventSkip", required_argument, NULL, 's'},
                                {"eventMax", required_argument, NULL, 'm'},
                                {"eudaqRawFile", required_argument, NULL, 'e'},
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
      case 'e':
        eudaqRawFilePath = optarg;
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
  std::fprintf(stdout, "geometryFile:  <%s>\n", geometryFilePath.c_str());
  std::fprintf(stdout, "eudaqRawFile:  <%s>\n", eudaqRawFilePath.c_str());
  std::fprintf(stdout, "\n");

  //////////// geometry

  std::string str_geo;
  if(geometryFilePath.empty()){
    str_geo = default_geometry;
  }
  else{
    str_geo = JsonUtils::readFile(geometryFilePath);
  }
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFilePath.c_str() );
    throw;
  }

  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }

  auto reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), eudaqRawFilePath);

  if(!reader){
    std::fprintf(stderr, "raw file is not open\n");
    throw;
  }

  TelFW telfw(800, 400, "test");
  glfw_test telfwtest(str_geo, true);
  telfw.startAsync<glfw_test>(&telfwtest, &glfw_test::beginHook, &glfw_test::clearHook, &glfw_test::drawHook);

  for(size_t eventNum = 0;;eventNum++){
    auto eudaqEvent = reader->GetNextEvent();
    if(!eudaqEvent){
      std::fprintf(stdout, "reach end of data file\n");
      break;
    }
    if(eventNum> eventMaxNum && eventMaxNum>0){
      std::fprintf(stdout, "reach to eventMaxNum set by option\n");
      break;
    }
    if(eventNum < eventSkipNum){
      continue;
    }
    std::shared_ptr<altel::TelEvent> telEvent = altel::createTelEvent(eudaqEvent);
    std::fprintf(stdout, "FileEvent #%d, event #%d, clock/trigger #%d\n", eventNum, telEvent->eveN(), telEvent->clkN());
    telfwtest.pushBufferEvent(telEvent);

    if(do_wait){
      std::fprintf(stdout, "waiting, press any key to next event\n");
      std::getc(stdin);
    }
  }

  reader.reset();
  if(do_wait){
    std::fprintf(stdout, "finished, press any key to exit\n");
    std::getc(stdin);
  }

  telfw.stopAsync();
  return 0;
}
