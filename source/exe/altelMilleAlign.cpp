#include "TelMille.hh"
 
#include "getopt.h"
#include "Acts/Utilities/Units.hpp" 
#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"
#include "TelEventTTreeReader.hpp"
#include <iostream>
#include <algorithm>
#include <set>
 
#include <TFile.h>
#include <TTree.h> 
using namespace Acts::UnitLiterals;  // 解决_um问题

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -eudaqFiles  <PATH0 [PATH1]...>   paths to input eudaq raw data files (input)
  -rootfile          <PATH>          path to input ROOT file (input)
  -maxEventNumber    <int>           max number of events to be processed
  -maxTrackNumber    <int>           max number of tracks to be processed
  -pedeSteeringFile  <PATH>          path to pede steering file (output)
  -milleBinaryFile   <PATH>          path to mille binary file (output)
  -inputGeometryFile <PATH>          geometry input file
  -resolDefault   <  <float_UV>|<float_U float_V> >
                                    default U/V resolution for all detectors
  -resolDetector  <<int_ID>  <<float_UV>|<float_U float_V> >>
                                    U/V resolution(s) for a specific detector by int_ID
  
  -ResolDefault     [PATH] path to resolution config JSON file
  -rootMode raw     The ROOT file is processed according to the original measurement mode.
  -rootMode traj    The ROOT files are processed according to the trajectory filtering mode.
  -maxTrackChi2Ndf  <float_X float_Y>
                                    reject tracks with prefit chi2/ndf above the cut
  -maxTrackResidual <float_X float_Y>
                                    reject tracks with max abs prefit residual above the cut [mm]


example:
./altelMilleBin -pede pede.txt -mille mille.bin  -eudaqFiles  eudaqRaw/altel_Run069017_200824002945.raw eudaqRaw/altel_Run069018_200824003322.raw -input ../init_geo.json -maxE 1000000 -resolDefault 0.04 0.04
 
example with ROOT file:
./altelMilleBin -pede pede.txt -mille mille.bin -rootfile run333.root -input geo_run333_init.json -maxE 100000000 -resolDefault 0.1
)";

// Resolution configuration structure
struct ResolutionConfig {
  double u_size1 = 6_um, u_size2 = 6_um, u_sizeOther = 6_um;
  double v_size1 = 6_um, v_size2 = 6_um, v_sizeOther = 6_um;
  bool loaded = false;
 
  void loadFromFile(const std::string& filename) {
    std::string content = JsonUtils::readFile(filename);
    JsonDocument doc = JsonUtils::createJsonDocument(content);
 
    if (doc.HasMember("resolutions")) {
      const auto& res = doc["resolutions"];
      if (res.HasMember("u_direction")) {
        const auto& u = res["u_direction"];
        if (u.HasMember("size_1")) u_size1 = std::stod(u["size_1"].GetString()) * 1_um ;
        if (u.HasMember("size_2")) u_size2 = std::stod(u["size_2"].GetString()) * 1_um ;
        if (u.HasMember("size_other")) u_sizeOther = std::stod(u["size_other"].GetString()) * 1_um ;
      }
      if (res.HasMember("v_direction")) {
        const auto& v = res["v_direction"];
        if (v.HasMember("size_1")) v_size1 = std::stod(v["size_1"].GetString()) * 1_um ;
        if (v.HasMember("size_2")) v_size2 = std::stod(v["size_2"].GetString()) * 1_um ;
        if (v.HasMember("size_other")) v_sizeOther = std::stod(v["size_other"].GetString()) * 1_um ;
      }
      loaded = true;
      std::printf("Resolution config loaded from: %s\n", filename.c_str());
    }
  }
};

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"eudaqFiles", required_argument, NULL, 'f'},
                              {"inputGeometryFile", required_argument, NULL, 'g'},
                              {"pedeSteeringFile", required_argument, NULL, 'u'},
                              {"milleBinaryFile", required_argument, NULL, 'q'},
                              {"resolDefault", required_argument, NULL, 'r'},
                              {"resolDetector", required_argument, NULL, 's'},
                              {"maxEventNumber", required_argument, NULL, 'm'},
                              {"maxTrackNumber", required_argument, NULL, 'n'},
                              {"rootfile", required_argument, NULL, 'R'},
                              {"ResolDefault", required_argument,NULL, 'D'},
                              {"rootMode", required_argument, NULL, 'M'},
                              {"maxTrackChi2Ndf", required_argument, NULL, 'C'},
                              {"maxTrackResidual", required_argument, NULL, 'A'},
                             {0, 0, 0, 0}};
 
  std::vector<std::string> rawFilePathCol;
  std::string inputGeometryFile_path;
  std::string pedeSteeringFile_path;
  std::string milleBinaryFile_path;
  std::string rootfile_path;

  std::string rootMode = "raw";

  std::string resolDefault_path;
  ResolutionConfig resolConfig;
  
  size_t maxTrackNumber = -1;
  size_t maxEventNumber = -1;
 
  std::map<uint16_t, std::pair<double, double>> mapResolDet;
  double resolDefaultU =0.03;
  double resolDefaultV =0.03;
  double maxTrackChi2NdfX = -1.0;
  double maxTrackChi2NdfY = -1.0;
  double maxTrackResidualX = -1.0;
  double maxTrackResidualY = -1.0;
 
  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'f':{
      optind--;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        const char* fileStr = argv[optind];
        rawFilePathCol.push_back(std::string(fileStr));
      }
      break;
    }
    case 'g':
      inputGeometryFile_path = optarg;
      break;
    case 'u':
      pedeSteeringFile_path = optarg;
      break;
    case 'q':
      milleBinaryFile_path = optarg;
      break;
    case 'D':
      resolDefault_path = optarg;
      resolConfig.loadFromFile(resolDefault_path);
      break;
    case 'r':{
      optind--;
      std::vector<double> resolVec;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        resolVec.push_back(std::stod(argv[optind]));
      }
      if(resolVec.size()==1){
        resolDefaultU=resolVec[0];
        resolDefaultV=resolVec[0];
      }
      else if(resolVec.size()==2){
        resolDefaultU=resolVec[0];
        resolDefaultV=resolVec[1];
      }
      break;
    }
    case 's':{
      optind--;
      uint16_t detN = -1;
      std::vector<double> resolVec;
      bool isFirstArg = true;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        if(isFirstArg){
          detN = std::stoi(argv[optind]);
          isFirstArg = false;
        }
        else{
          resolVec.push_back(std::stod(argv[optind]));
        }
      }
      if(resolVec.size()==1){
        mapResolDet[detN]=std::pair<double, double>(resolVec[0], resolVec[0]);
      }
      else if(resolVec.size()==2){
        mapResolDet[detN]=std::pair<double, double>(resolVec[0], resolVec[1]);
      }
      break;
    }
    case 'm':
      maxEventNumber = std::stoull(optarg);
      break;
    case 'n':
      maxTrackNumber = std::stoull(optarg);
      break;
    case 'R':
      rootfile_path = std::string(optarg);
      break;
    case 'M':
      rootMode = std::string(optarg);
      break;
    case 'C':{
      optind--;
      std::vector<double> cutVec;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        cutVec.push_back(std::stod(argv[optind]));
      }
      if(cutVec.size()==1){
        maxTrackChi2NdfX = cutVec[0];
        maxTrackChi2NdfY = cutVec[0];
      }
      else if(cutVec.size()==2){
        maxTrackChi2NdfX = cutVec[0];
        maxTrackChi2NdfY = cutVec[1];
      }
      break;
    }
    case 'A':{
      optind--;
      std::vector<double> cutVec;
      for( ;optind < argc && *argv[optind] != '-'; optind++){
        cutVec.push_back(std::stod(argv[optind]));
      }
      if(cutVec.size()==1){
        maxTrackResidualX = cutVec[0];
        maxTrackResidualY = cutVec[0];
      }
      else if(cutVec.size()==2){
        maxTrackResidualX = cutVec[0];
        maxTrackResidualY = cutVec[1];
      }
      break;
    }
      /////generic part below///////////
    case 0: /* getopt_long() set a variable, just keep going */
      break;
    case 1:
      fprintf(stderr, "case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr, "case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr, "case ?\n");
      exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }
 
  if ((rawFilePathCol.empty() && rootfile_path.empty()) ||
      inputGeometryFile_path.empty() ||
      milleBinaryFile_path.empty() ||
      pedeSteeringFile_path.empty()) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(0);
  }
 
  std::fprintf(stdout, "\n");
  if(!rawFilePathCol.empty()){
    std::fprintf(stdout, "%zu eudaqFiles:\n", rawFilePathCol.size());
    for(auto &rawfilepath: rawFilePathCol){
      std::fprintf(stdout, "  %s\n", rawfilepath.c_str());
    }
  }
  if(!rootfile_path.empty()){
    std::fprintf(stdout, "ROOT file:  %s\n", rootfile_path.c_str());
  }
  if (rootMode != "raw" && rootMode != "traj") {
  std::fprintf(stderr, "Invalid -rootMode <%s>, must be raw or traj\n", rootMode.c_str());
  std::exit(1);
  }
  if (do_help) {
  std::fprintf(stdout, "%s\n", help_usage.c_str());
  return 0;
  }
  std::fprintf(stdout, "inputGeometryFile:  %s\n", inputGeometryFile_path.c_str());
  std::fprintf(stdout, "milleBinaryFile:    %s\n", milleBinaryFile_path.c_str());
  std::fprintf(stdout, "pedeSteeringFile:   %s\n", pedeSteeringFile_path.c_str());
  //++++++++++++++++++++
  if (resolConfig.loaded) {
    std::fprintf(stdout, "Resolution config: loaded from %s\n", resolDefault_path.c_str());
    std::fprintf(stdout, "  U: size_1=%.6f, size_2=%.6f, size_other=%.6f\n", 
                 resolConfig.u_size1, resolConfig.u_size2, resolConfig.u_sizeOther);
    std::fprintf(stdout, "  V: size_1=%.6f, size_2=%.6f, size_other=%.6f\n", 
                 resolConfig.v_size1, resolConfig.v_size2, resolConfig.v_sizeOther);
  } else {
    std::fprintf(stdout, "resolDefault:       [%f   %f]\n", resolDefaultU, resolDefaultV);//原数据
  }
  std::fprintf(stdout, "resolDetector:\n");//原数据
  //+++++++++++
  std::fprintf(stdout, "maxTrackChi2Ndf:    [%f   %f]\n", maxTrackChi2NdfX, maxTrackChi2NdfY);
  std::fprintf(stdout, "maxTrackResidual:   [%f   %f]\n", maxTrackResidualX, maxTrackResidualY);
 /* for(auto &[detN, resolUV]: mapResolDet){
    std::fprintf(stdout, "  det #%d:  [%f   %f]\n", detN, resolUV.first, resolUV.second);
  }*/
 
 
  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(inputGeometryFile_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", inputGeometryFile_path.c_str() );
    throw;
  }
  // JsonUtils::printJsonValue(jsd_geo, true);
 
  std::set<uint16_t> geoDetNs;
  for(const auto& geoJS: jsd_geo["geometry"]["detectors"].GetArray()){
    geoDetNs.insert(geoJS["id"].GetUint());
  }
  if(geoDetNs.size() != jsd_geo["geometry"]["detectors"].Size()){
    std::fprintf(stderr, "geometory id, something wrong\n");
    throw;
  }
  for(auto &detN: geoDetNs){
    if(mapResolDet.find(detN)==mapResolDet.end()){
      mapResolDet[detN]=std::pair<double, double>(resolDefaultU, resolDefaultV);
    }
  }
 
 
  altel::TelMille telmille;
  telmille.setGeometry(jsd_geo);
  telmille.setTrackQualityCuts(maxTrackChi2NdfX, maxTrackChi2NdfY,
                               maxTrackResidualX, maxTrackResidualY);
  for(auto& [detN, resolUV ]: mapResolDet){
   //++++++++++
    if (resolConfig.loaded) {
      telmille.setResolutionConfig(resolConfig.u_size1, resolConfig.u_size2, resolConfig.u_sizeOther,
                                   resolConfig.v_size1, resolConfig.v_size2, resolConfig.v_sizeOther);
    } else {
      telmille.setResolution(detN, resolUV.first, resolUV.second);//原数据
    }
  }//++++++++++
  telmille.startMilleBinary(milleBinaryFile_path);
 
  JsonAllocator jsa;
 
  size_t nTracks = 0;
  size_t nEvents = 0;
 
  // Check if using ROOT file or EUDAQ files
  bool useRootFile = !rootfile_path.empty();
  
  if(useRootFile){
    // Read from ROOT file
    TFile *tfile = new TFile(rootfile_path.c_str(), "READ");
    if(!tfile || !tfile->IsOpen()){
      std::fprintf(stderr, "ROOT file <%s> cannot be opened.\n", rootfile_path.c_str());
      std::exit(1);
    }
 
    TTree *pTree = 0;
    tfile->GetObject("eventTree", pTree);
    if(!pTree){
      std::fprintf(stderr, "eventTree not found in ROOT file.\n");
      std::exit(1);
    }
 
    altel::TelEventTTreeReader ttreeReader;
    ttreeReader.setTTree(pTree);
    size_t totalNumEvents = ttreeReader.numEvents();
 
    std::fprintf(stdout, "ROOT mode: %s\n", rootMode.c_str());

    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
      if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
        break;
      }
 
      std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
      if(!telEvent){
        continue;
      }
 
      nEvents++;

//============================
      if(rootMode == "raw"){
      JsonValue js_track_filtered(rapidjson::kArrayType);
      std::set<uint16_t> measDetNs;
      bool isMoreThanOneHitPerPlane = false;

      std::fprintf(stdout, "\n[RAW] Event #%zu ", nEvents);

      for(const auto &aMeasHit: telEvent->measHits()){
        if(!aMeasHit){
          continue;
        }

        uint16_t detN = aMeasHit->detN();
        double measU = aMeasHit->u();
        double measV = aMeasHit->v();

        std::fprintf(stdout, "[%f, %f, %d] ", measU, measV, detN);

        if(geoDetNs.find(detN) == geoDetNs.end()){
          continue;
        }

        auto [it_insert, isSuccess] = measDetNs.insert(detN);
        if(!isSuccess){
          isMoreThanOneHitPerPlane = true;
          break;
        }

        std::set<uint16_t> unique_u, unique_v;
        for(const auto& raw : aMeasHit->measRaws()){
          unique_u.insert(raw.u());
          unique_v.insert(raw.v());
        }

        JsonValue js_hit(rapidjson::kObjectType);
        js_hit.AddMember("id", detN, jsa);
        js_hit.AddMember("x", measU, jsa);
        js_hit.AddMember("y", measV, jsa);
        js_hit.AddMember("cluster_size_u", static_cast<uint32_t>(unique_u.size()), jsa);
        js_hit.AddMember("cluster_size_v", static_cast<uint32_t>(unique_v.size()), jsa);

        js_track_filtered.PushBack(std::move(js_hit), jsa);
      }

      if(isMoreThanOneHitPerPlane || measDetNs.size() != geoDetNs.size()){
        std::cout << "skipping event " << nEvents << std::endl;
        continue;
      }

      if(telmille.fillTrackXYRz(js_track_filtered)){
        nTracks++;
      }
    }

    // =========================================================
    // Mode 2: traj
    // 用 traj 挑好轨迹，但真正用于 alignment 的是原始 measHit
    // =========================================================
    
      else if(rootMode == "traj"){
        std::fprintf(stdout, "\n[TRAJ] Event #%zu ", nEvents);
        for(auto &aTraj : telEvent->trajs()){
          if(!aTraj){
            continue;
          }
          std::set<uint16_t> usedDetNs;
          JsonValue js_track_filtered(rapidjson::kArrayType);
          bool isValidTraj = true;
          for(const auto& trajHit : aTraj->trajHits()){
            if(!trajHit){
              isValidTraj = false;
              break;
            }
          // 这条轨迹在该 plane 上需要有拟合结果
            if(!trajHit->hasFitHit()){
              isValidTraj = false;
              break;
            }
            uint16_t detN = trajHit->detN();
            if(geoDetNs.find(detN) == geoDetNs.end()){
              continue;
            }

          // 每个 plane 只能有一个 hit

            auto [it_insert, isSuccess] = usedDetNs.insert(detN);
            if(!isSuccess){
              isValidTraj = false;
              break;
            }

          // 真正用于 alignment 的是原始测量 hit
          // 优先使用 matchedMeasHit，其次使用 fitHit 挂的 originMeasHit        
            std::shared_ptr<altel::TelMeasHit> measHit = nullptr;
            if(trajHit->hasMatchedMeasHit()){            
              measHit = trajHit->matchedMeasHit();
            }
            else if(trajHit->hasOriginMeasHit()){
              measHit = trajHit->fitHit()->originMeasHit();
            }
            if(!measHit){
              isValidTraj = false;
              break;
            }
            if(measHit->detN() != detN){
              std::fprintf(stderr,"detector mismatch: trajHit det=%u, measHit det=%u\n",detN, measHit->detN());
              isValidTraj = false;
              break;
            }
            double measU = measHit->u();
            double measV = measHit->v();
            std::fprintf(stdout, "[%f, %f, %d] ", measU, measV, detN);          
            JsonValue js_hit(rapidjson::kObjectType);
            js_hit.AddMember("id", detN, jsa);
            js_hit.AddMember("x", measU, jsa);
            js_hit.AddMember("y", measV, jsa);

          // 从原始 cluster/raw pixel 恢复 cluster size
          
            std::set<uint16_t> unique_u, unique_v;          
            for(const auto& raw : measHit->measRaws()){
              unique_u.insert(raw.u());
              unique_v.insert(raw.v());
            }
         
            js_hit.AddMember("cluster_size_u",static_cast<uint32_t>(unique_u.size()), jsa);
            js_hit.AddMember("cluster_size_v",static_cast<uint32_t>(unique_v.size()), jsa);
            js_track_filtered.PushBack(std::move(js_hit), jsa);
          }

        // 轨迹不完整，丢弃
          if(!isValidTraj){
            continue;
          }

        // 必须覆盖所有几何 plane    
          if(usedDetNs.size() == geoDetNs.size()){
            if(telmille.fillTrackXYRz(js_track_filtered)){
              nTracks++;
            }
          }
        }
      }
    }

    tfile->Close();
    delete tfile;
  }
  else{
    // Read from EUDAQ files
    uint32_t rawFileN=0;
    eudaq::FileReaderUP reader;
 
    while(1){
      if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
        break;
      }
 
      if(!reader){
        if(rawFileN<rawFilePathCol.size()){
          std::fprintf(stdout, "processing raw file: %s\n", rawFilePathCol[rawFileN].c_str());
          reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), rawFilePathCol[rawFileN]);
          rawFileN++;
        }
        else{
          std::fprintf(stdout, "processed %zu raw files, quit\n", rawFileN);
          break;
        }
      }
      auto eudaqEvent = reader->GetNextEvent();
      if(!eudaqEvent){
        reader.reset();
        continue;
      }
      nEvents++;
 
      std::shared_ptr<altel::TelEvent> telEvent = altel::createTelEvent(eudaqEvent);
 
      // TODO: TelEvent to json
      JsonValue js_track_filtered(rapidjson::kArrayType);
      std::set<uint16_t> measDetNs;
      bool isMoreThanOneHitPerPlane = false;
      std::fprintf(stdout, "\nEvent #%zu ", nEvents);
      for(const auto &aMeasHit: telEvent->measHits()){
        uint16_t detN = aMeasHit->detN();
        double measU = aMeasHit->u();
        double measV = aMeasHit->v();
 
        std::fprintf(stdout, "[%f, %f, %d] ", measU, measV, detN);
        if(geoDetNs.find(detN) == geoDetNs.end()){
          continue;
        }
        auto [it_insert, isSucess] = measDetNs.insert(detN);
        if(!isSucess){
          isMoreThanOneHitPerPlane = true;
          break;
        }
//++++++
        std::set<uint16_t> unique_u, unique_v;
        for(const auto& raw : aMeasHit->measRaws()){
          unique_u.insert(raw.u());
          unique_v.insert(raw.v());
        }
        size_t clusterSizeU = unique_u.size();
        size_t clusterSizeV = unique_v.size();
//++++++++
        JsonValue js_hit(rapidjson::kObjectType);
        js_hit.AddMember("id", detN, jsa);
        js_hit.AddMember("x", measU, jsa);
        js_hit.AddMember("y", measV, jsa);
//++++++++
        js_hit.AddMember("cluster_size_u", clusterSizeU, jsa);
        js_hit.AddMember("cluster_size_v", clusterSizeV, jsa);
        //++++++++
        js_track_filtered.PushBack(std::move(js_hit), jsa);
      }
 
      if(isMoreThanOneHitPerPlane || measDetNs.size() != geoDetNs.size()){
        std::cout<< "skipping event "<<nEvents <<std::endl;
        continue;
      }
 
      if(telmille.fillTrackXYRz(js_track_filtered)){
        nTracks++;
      }
    }
  }
 
  std::fprintf(stdout, "\n%zu tracks are picked from %zu events\n", nTracks, nEvents);
 
  telmille.endMilleBinary();
  telmille.createPedeStreeringModeXYRz(pedeSteeringFile_path);
 
  return 0;
}
