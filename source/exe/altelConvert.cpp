#include "TelEventTTreeWriter.hpp"
#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

#include <numeric>
#include <chrono>
#include <regex>

#include <TFile.h>
#include <TTree.h>
#include <Math/SpecFunc.h>
#include <Math/DistFunc.h>

#include "linenoise.h"
#include "myrapidjson.h"

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -eventSkip      <INT>             number of events to skip before start processing (default 0, disabled)
  -eventMax       <INT>             max number of events to process  (default -1, disabled)
  -daqFiles  <<PATH0> [PATH1]...>   paths to input daq data files (input). old option -eudaqFiles
  -rootFile       <PATH>            path to out root file of reconstructed trajactories (output)

examples:
./altelConvert  -daqFiles eudaqRaw/altel_Run069017_200824002945.raw  -rootFile detresid.root -eventMax 10000
)";

int main(int argc, char *argv[]) {
  int64_t eventMaxNum = 0;
  int64_t eventSkipNum = 0;
  std::set<int64_t> includeDetId;
  std::set<int64_t> excludeDetId;
  std::set<int64_t> targetDetId;

  std::vector<std::string> rawFilePathCol;
  std::string geometryFilePath;
  std::string rootFilePath;


  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"eventSkip", required_argument, NULL, 's'},
                                {"eventMax", required_argument, NULL, 'm'},
                                {"daqFiles", required_argument, NULL, 'f'},
                                {"rootFile", required_argument, NULL, 'b'},
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
      case 'f':{
        optind--;
        for( ;optind < argc && *argv[optind] != '-'; optind++){
          const char* fileStr = argv[optind];
          rawFilePathCol.push_back(std::string(fileStr));
        }
        break;
      }
      case 'b':
        rootFilePath = optarg;
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

  std::fprintf(stdout, "%d daqFiles:\n", rawFilePathCol.size());
  for(auto &rawfilepath: rawFilePathCol){
    std::fprintf(stdout, "                %s\n", rawfilepath.c_str());
  }

  if (rawFilePathCol.empty() ||
      rootFilePath.empty() ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(1);
  }
  /////////////////////////////////////

  altel::TelEventTTreeWriter ttreeWriter;
  TTree *pTree = new TTree("eventTree", "eventTree");
  ttreeWriter.setTTree(pTree);


  uint32_t rawFileNum=0;
  size_t readEventNum = 0;
  size_t emptyEventNum = 0;
  size_t eventNum = 0;
  auto tp_start = std::chrono::system_clock::now();

  eudaq::FileReaderUP reader;
  std::unique_ptr<JsonFileDeserializer> jsreader;

  bool is_eudaq_raw = true;
  if(std::regex_match(rawFilePathCol.front(), std::regex("\\S+.json")) ){
    is_eudaq_raw= false;
  }
  while(1){
    if(eventNum> eventMaxNum && eventMaxNum>0){
      break;
    }
    std::shared_ptr<altel::TelEvent> fullEvent;
    if(is_eudaq_raw){
      if(!reader){
        if(rawFileNum<rawFilePathCol.size()){
          std::fprintf(stdout, "processing raw file: %s\n", rawFilePathCol[rawFileNum].c_str());
          reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), rawFilePathCol[rawFileNum]);
          rawFileNum++;
        }
        else{
          std::fprintf(stdout, "processed %d raw files, quit\n", rawFileNum);
          break;
        }
      }
      auto eudaqEvent = reader->GetNextEvent();
      if(!eudaqEvent){
        reader.reset();
        continue; // goto for next raw file
      }
      readEventNum++;
      if(readEventNum<=eventSkipNum){
        continue;
      }
      eventNum++;
      fullEvent  = altel::createTelEvent(eudaqEvent);
    }
    else{
      if(!jsreader){
        if(rawFileNum<rawFilePathCol.size()){
          std::fprintf(stdout, "processing js file: %s\n", rawFilePathCol[rawFileNum].c_str());
          jsreader.reset(new JsonFileDeserializer(rawFilePathCol[rawFileNum]));
          rawFileNum++;
        }
        else{
          std::fprintf(stdout, "processed %d raw files, quit\n", rawFileNum);
          break;
        }
      }
      auto evpack = jsreader->getNextJsonDocument();
      if(evpack.IsNull()){
        jsreader.reset();
        continue;
      }
      eventNum++;
      fullEvent  = TelActs::createTelEvent(evpack, 0, eventNum, 0);
    }


    ttreeWriter.fillTelEvent(fullEvent);
    eventNum ++;
  }

  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;
  double time_s = dur_diff.count();

  TFile tfile(rootFilePath.c_str(),"recreate");
  pTree->Write();
  tfile.Close();
  return 0;
}
