#include "TelEventTTreeReader.hpp"
#include "getopt.h"
#include "myrapidjson.h"
 
#include <numeric>
#include <chrono>
#include <iostream>
#include <stdexcept>
 
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
 ./altelTelEventViewer -w -geo ../../testbeam_data_2507/RUN/geo_setup2_align3_0p04.json  -r  detresid.root
)";
 
int main(int argc, char *argv[]) {
    int64_t eventMaxNum = 0;
    int64_t eventSkipNum = 0;
    std::string geometryFilePath;
    std::string rootFilePath;
    int totalGoodTrajectories = 0;
    int totalAllTrajectories = 0;
    int totalMatchedHits = 0;
    int goodTrajHits = 0;
    int totalFitHits = 0;
    int totaloriginMeasHits=0;
    int alltotalTrajHits=0;
    int do_wait = 0;
    int do_verbose = 0;
 
    //////////////getopt begin//////////////////
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", no_argument, NULL, 'v'},
        {"wait", no_argument, NULL, 'w'},
        {"eventSkip", required_argument, NULL, 's'},
        {"eventMax", required_argument, NULL, 'm'},
        {"rootFile", required_argument, NULL, 'b'},
        {"geometryFile", required_argument, NULL, 'g'},
        {0, 0, 0, 0}
    };
 
    if(argc == 1){
        std::fprintf(stderr, "%s\n", help_usage.c_str());
        return 1;
    }
 
    int c;
    int longindex;
    opterr = 1;
    while ((c = getopt_long_only(argc, argv, "-", longopts, &longindex)) != -1) {
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
            do_wait = 1;
            break;
        case 'v':
            do_verbose = 1;
            if(optind < argc && *argv[optind] != '-'){
                do_verbose = std::stoul(argv[optind]);
                optind++;
            }
            break;
        case 'h':
            std::fprintf(stdout, "%s\n", help_usage.c_str());
            return 0;
        case ':':
            std::fprintf(stderr, "%s: missing argument for option %s\n",
                         argv[0], longopts[longindex].name);
            return 1;
        case '?':
            return 1;
        default:
            std::fprintf(stderr, "%s: missing getopt branch %c for option %s\n",
                         argv[0], c, longopts[longindex].name);
            return 1;
        }
    }
    ///////////getopt end////////////////
 
    std::fprintf(stdout, "\n");
    std::fprintf(stdout, "geometryFile:  %s\n", geometryFilePath.c_str());
    std::fprintf(stdout, "rootFile:      %s\n", rootFilePath.c_str());
    std::fprintf(stdout, "\n");
 
    // Check if required files are provided
    if(geometryFilePath.empty() || rootFilePath.empty()){
        std::fprintf(stderr, "Error: Both geometryFile and rootFile must be provided\n");
        return 1;
    }
 
    //////////// geometry loading
    std::string str_geo = JsonUtils::readFile(geometryFilePath);
    if(str_geo.empty()){
        std::fprintf(stderr, "Error: Cannot read geometry file <%s>\n", geometryFilePath.c_str());
        return 1;
    }
 
    JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
    if(jsd_geo.IsNull()){
        std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFilePath.c_str());
        return 1;
    }
 
    if (!jsd_geo.HasMember("geometry")) {
        std::fprintf(stderr, "Error: Geometry file missing 'geometry' member\n");
        return 1;
    }
 
    altel::TelEventTTreeReader ttreeReader;
 
    // ROOT file handling
    TFile tfile(rootFilePath.c_str(),"READ");
    if(!tfile.IsOpen()){
        std::fprintf(stderr, "Error: Cannot open ROOT file <%s>\n", rootFilePath.c_str());
        return 1;
    }
 
    TTree *pTree = nullptr;
    tfile.GetObject("eventTree",pTree);
    if(!pTree){
        std::fprintf(stderr, "Error: 'eventTree' not found in ROOT file\n");
        tfile.Close();
        return 1;
    }
 
    if(pTree->GetEntries() == 0){
        std::fprintf(stderr, "Error: Tree has no entries\n");
        tfile.Close();
        return 1;
    }
 
    ttreeReader.setTTree(pTree);
 
    // Initialize visualization
    TelFW telfw(800, 400, "test");
    glfw_test telfwtest(geometryFilePath);
    telfw.startAsync<glfw_test>(&telfwtest, &glfw_test::beginHook, &glfw_test::clearHook, &glfw_test::drawHook);
 
    auto tp_start = std::chrono::system_clock::now();
 
    size_t totalNumEvents = ttreeReader.numEvents();
    std::fprintf(stdout, "Total events in file: %zu\n", totalNumEvents);
 
    for(size_t eventNum = eventSkipNum; eventNum < totalNumEvents; eventNum++){
        if(eventMaxNum > 0 && eventNum >= eventSkipNum + eventMaxNum){
            break;
        }
 
        // 创建事件时进行错误检查
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        if(!telEvent){
            std::fprintf(stderr, "Warning: Failed to create TelEvent for event %zu\n", eventNum);
            continue;
        }
 
        std::fprintf(stdout, "FileEvent #%zu, event #%d, clock/trigger #%d\n",
                    eventNum, telEvent->eveN(), telEvent->clkN());
 
        // Process measurement hits
        std::map<uint32_t, std::vector<std::shared_ptr<altel::TelMeasHit>>> map_layer_measHits;
        for(auto& mh: telEvent->measHits()){
            if(!mh) continue;
            uint32_t detN = mh->detN();
            map_layer_measHits[detN].push_back(mh);
        }
 
        for(auto& [detN, mhs]: map_layer_measHits){
            std::fprintf(stdout, "detN : %u, cluster_size: %zu ", detN, mhs.size());
            for(auto &mh : mhs){
                std::fprintf(stdout, " cluster(u,v) = (%.3f,%.3f)", mh->u(), mh->v());
                std::fprintf(stdout, "[Raw_size:%zu :", mh->measRaws().size());
                for(auto &mr : mh->measRaws()){
                    std::fprintf(stdout, "Raw(u,v)(%u,%u)", mr.u(), mr.v());
                }
                std::fprintf(stdout, "]");
            }
            std::fprintf(stdout, "\n");
        }
 
        // Process trajectories
  //     std::fprintf(stdout, "=============================");
        int goodTrajectoriesCount = 0;
        if (!telEvent->trajs().empty()) {
            std::fprintf(stdout, "trajsNum: %zu\n", telEvent->trajs().size());
            size_t trajIndex = 0;
            for (auto aTraj: telEvent->trajs()) {
              trajIndex++;
              totalAllTrajectories++;
              alltotalTrajHits+=aTraj->numTrajHit();
              bool isGoodTrajectory = (aTraj->numOriginMeasHit() == 5);
              if(isGoodTrajectory){
                totalGoodTrajectories++;
                goodTrajectoriesCount++;
                std::fprintf(stdout, "  good traja[%d]: Including points=%zu\n", goodTrajectoriesCount, aTraj->numOriginMeasHit());
                std::fprintf(stdout, "  === Cluster Information for Good Trajectory %d ===\n", goodTrajectoriesCount);
                for (auto &aTrajHit: aTraj->trajHits()) {
                  goodTrajHits++;
                 // std::fprintf(stdout, "    TrajHit[%zu]:\n", j);
                  auto aFitHit = aTrajHit->fitHit();
                  auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
                  if(aFitHit){
                    totalFitHits++;
                    std::fprintf(stdout, "      FitHit: detN=%d,local pos=(%.3f, %.3f),global pos=(%.2f,%.2f,%.2f)\n",
                                        aFitHit->detN(), aFitHit->u(), aFitHit->v(),aFitHit->x(),aFitHit->y(),aFitHit->z());
                    auto aOriginMeasHit = aFitHit->originMeasHit();
                    if (aOriginMeasHit) {
                      totaloriginMeasHits++;
                      std::fprintf(stdout, "      OriginMeasHit: detN=%d, pos=(%.3f, %.3f), cluster_size=%zu\n",
                                               aOriginMeasHit->detN(), aOriginMeasHit->u(), aOriginMeasHit->v(),aOriginMeasHit->measRaws().size());
                      std::fprintf(stdout, "        Raw pixels: ");
                      for (const auto& raw : aOriginMeasHit->measRaws()) {
                        std::fprintf(stdout, "(%u,%u) ", raw.u(), raw.v());
                      }
                      std::fprintf(stdout, "\n");
                    }
                  }
                  if (aMatchedMeasHit) {
                    //totalMatchedHits++;
                    totalMatchedHits+=aTraj->numMatchedMeasHit();
                    std::fprintf(stdout, "      Status: MATCHED with measurement\n");
                    std::fprintf(stdout, "      MatchedMeasHit: detN=%d, pos=(%.3f, %.3f), cluster_size=%zu\n",
                                         aMatchedMeasHit->detN(), aMatchedMeasHit->u(), aMatchedMeasHit->v(),aMatchedMeasHit->measRaws().size());
                    std::fprintf(stdout, "        Matched raw pixels: ");
                    for (const auto& raw : aMatchedMeasHit->measRaws()) {
                      std::fprintf(stdout, "(%u,%u) ", raw.u(), raw.v());
                    }
                    std::fprintf(stdout, "\n");
                  }else {
                    std::fprintf(stdout, "      Status: NO MATCH (prediction only)\n");
                  }
                  std::fprintf(stdout, "    ---\n");
                }
              }else{
                std::fprintf(stdout, "  Bad trajectory[%zu]: only %zu origin hits\n",
                         trajIndex, aTraj->numOriginMeasHit());
              }
            }
            std::fprintf(stdout, "  === End Cluster Information ===\n");
        }
        std::fprintf(stdout, "The total number of good trajectories for this event: %d/%zu\n",
                        goodTrajectoriesCount, telEvent->trajs().size());
        telfwtest.pushBufferEvent(telEvent);
        if(do_wait){
          std::fprintf(stdout, "waiting, press any key to next event\n");
          std::getc(stdin);
        }
    } 
    auto tp_end = std::chrono::system_clock::now();
    std::chrono::duration<double> dur_diff = tp_end - tp_start;
    std::fprintf(stdout, "Processing time: %.2f seconds\n", dur_diff.count());
    std::fprintf(stdout, "\n=== MATCHING STATISTICS ===\n");
    std::fprintf(stdout, "Total trajectories: %d\n", totalAllTrajectories);
    std::fprintf(stdout, "Total good trajectories: %d\n",totalGoodTrajectories);
    std::fprintf(stdout, "allTotal trajectory hits: %d\n", alltotalTrajHits);
    std::fprintf(stdout, "good trajectory hits: %d\n", goodTrajHits);
    std::fprintf(stdout, "Total fit hits: %d\n", totalFitHits);
    std::fprintf(stdout, "Total origin hits: %d\n",totaloriginMeasHits );
    std::fprintf(stdout, "Total matched hits: %d\n", totalMatchedHits);
    tfile.Close();
 
    if(do_wait){
        std::fprintf(stdout, "waiting, press any key to exit\n");
        std::getc(stdin);
    }
 
    telfw.stopAsync();
    return 0;
  }
