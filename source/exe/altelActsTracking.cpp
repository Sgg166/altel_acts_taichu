#include "TelEventTTreeWriter.hpp"
#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

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

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -wait                             wait for user keyboard input per event
  -eventSkip      <INT>             number of events to skip before start processing
  -eventMax       <INT>             max number of events to process
  -geometryFile   <PATH>            path to geometry input file (input)
  -eudaqFiles  <<PATH0> [PATH1]...> paths to input eudaq raw data files (input)
  -rootFile       <PATH>            path to root file (output)
  -particleEnergy <FLOAT>           energy of beam particle, electron, (Gev)
  -targetIds    <<INT0> [INT1]...>  IDs of target detector which are complectely excluded from track fitting. Residual are caculated.

examples:
./bin/TelDetectorResidual -eudaqFiles eudaqRaw/altel_Run069017_200824002945.raw -geometryFile calice_geo_align4.json -rootFile detresid.root -targetIds 5 -eventMax 10000
)";

int main(int argc, char *argv[]) {
  int64_t eventMaxNum = 0;
  int64_t eventSkipNum = 0;
  std::vector<int64_t> targetDetId;

  std::vector<std::string> rawFilePathCol;
  std::string geometryFilePath;
  std::string rootFilePath;

  double particleQ = 1;
  double particleMass = 0.511 * Acts::UnitConstants::MeV;
  double particleEnergy = 5.0 * Acts::UnitConstants::GeV;

  double distCollimator = 5_m;
  double widthCollimator = 4_cm;
  int do_wait = 0;

  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"wait", no_argument, NULL, 'w'},
                                {"eventSkip", required_argument, NULL, 's'},
                                {"eventMax", required_argument, NULL, 'm'},
                                {"eudaqFiles", required_argument, NULL, 'f'},
                                {"rootFile", required_argument, NULL, 'b'},
                                {"geometryFile", required_argument, NULL, 'g'},
                                {"particleEnergy", required_argument, NULL, 'e'},
                                {"targetIds", required_argument, NULL, 'd'},
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
      case 'g':
        geometryFilePath = optarg;
        break;
      case 'e':
        particleEnergy = std::stod(optarg) * Acts::UnitConstants::GeV;
        break;
      case 'd':{
        //optind is increased by 2 when option is set to required_argument
        for(int i = optind-1; i < argc && *argv[i] != '-'; i++){
          targetDetId.push_back(std::stol(argv[i]));
          optind = i+1;
        }
        break;
      }
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
  std::fprintf(stdout, "%d eudaqFiles:\n", rawFilePathCol.size());
  for(auto &rawfilepath: rawFilePathCol){
    std::fprintf(stdout, "  %s\n", rawfilepath.c_str());
  }
  std::fprintf(stdout, "geometryFile:  %s\n", geometryFilePath.c_str());
  std::fprintf(stdout, "rootFile:      %s\n", rootFilePath.c_str());
  std::fprintf(stdout, "\n");

  if (rawFilePathCol.empty() ||
      rootFilePath.empty() ||
      geometryFilePath.empty()) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(1);
  }

  /////////////////////////////////////
  Acts::Logging::Level logLevel =
    do_verbose ? (Acts::Logging::VERBOSE) : (Acts::Logging::INFO);

  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  /////////////  track seed
  double seedResX = 0.5*widthCollimator;
  double seedResY = 0.5*widthCollimator;
  double seedResPhi = 0.5*widthCollimator/distCollimator;
  double seedResTheta = 0.5*widthCollimator/distCollimator;

  double seedResX2 = seedResX * seedResX;
  double seedResY2 = seedResY * seedResY;
  double seedResPhi2 = seedResPhi * seedResPhi;
  double seedResTheta2 = seedResTheta * seedResTheta;

  Acts::BoundSymMatrix seedCov;
  seedCov <<
    seedResX2,0.,       0.,         0.,           0.,     0.,
    0.,       seedResY2,0.,         0.,           0.,     0.,
    0.,       0.,       seedResPhi2,0.,           0.,     0.,
    0.,       0.,       0.,         seedResTheta2,0.,     0.,
    0.,       0.,       0.,         0.,           0.0001, 0.,
    0.,       0.,       0.,         0.,           0.,     1.;

  double seedPhi = 0;
  double seedTheta = 0.5*M_PI;
  Acts::Vector4D seedPos4(-distCollimator, 0, 0, 0);
  Acts::CurvilinearTrackParameters seedParameters(seedPos4, seedPhi, seedTheta,
                                                  particleEnergy, particleQ, seedCov);

  //////////// geometry
  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geometryFilePath);
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFilePath.c_str() );
    throw;
  }
  // JsonUtils::printJsonValue(jsd_geo, true);

  std::printf("--------create acts geo object-----\n");
  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }
  const auto &js_geo = jsd_geo["geometry"];
  const auto &js_dets = js_geo["detectors"];

  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>> mapDetId2PlaneLayer;
  std::map<std::shared_ptr<const Acts::PlaneLayer>, size_t> mapPlaneLayer2DetId;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> allPlaneLayers;
  for(const auto& js_det: js_dets.GetArray()){
    auto [detId, planeLayer] = TelActs::createPlaneLayer(js_det);
    mapDetId2PlaneLayer[detId] = planeLayer;
    mapPlaneLayer2DetId[planeLayer] = detId;
    allPlaneLayers.push_back(planeLayer);
  }
  if(allPlaneLayers.size()!=mapDetId2PlaneLayer.size()){
    std::fprintf(stderr, "Error: there must be duplicated detID\n");
    throw;
  }

  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>> mapDetId2PlaneLayer_dets;
  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>> mapDetId2PlaneLayer_targets;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> layerDets;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> layerTargets;


  for(auto &[detId, planeLayer] :mapDetId2PlaneLayer){
    auto it = find (targetDetId.begin(), targetDetId.end(), detId);
    if(it == targetDetId.end()){
      layerDets.push_back(planeLayer);
      mapDetId2PlaneLayer_dets[detId] = planeLayer;
    }
    else{
      layerTargets.push_back(planeLayer);
      mapDetId2PlaneLayer_targets[detId] = planeLayer;
    }
  }

  std::vector<uint16_t> detId_dets;
  for(const auto& [detId, layer]: mapDetId2PlaneLayer_dets){
    detId_dets.push_back(detId);
  }
  std::vector<uint16_t> detId_targets;
  for(const auto& [detId, layer]: mapDetId2PlaneLayer_targets){
    detId_targets.push_back(detId);
  }


  if(layerDets.size()<3){
    std::fprintf(stdout, "error: number of detector elements is only %d.", layerDets.size());
    throw;
  }

  Acts::GeometryObjectSorterT<std::shared_ptr<const Acts::PlaneLayer>> layerSorter(gctx, Acts::BinningValue::binX);
  std::sort(layerDets.begin(), layerDets.end(), layerSorter);

  std::fprintf(stdout, "elementN = %d, detN = %d, targetN = %d\n",
               mapDetId2PlaneLayer.size(), mapDetId2PlaneLayer_dets.size(), mapDetId2PlaneLayer_targets.size());

  std::shared_ptr<const Acts::TrackingGeometry> worldGeo =
    TelActs::createWorld(gctx, 11.0_m, 0.1_m, 0.1_m,  allPlaneLayers); //include collimator at -5_m
// geometry closed, geometry ID only valid after geometry closing

  std::map<Acts::GeometryIdentifier, size_t> mapGeoId2DetId;
  for(auto& [detId, aPlaneLayer]: mapDetId2PlaneLayer){
    Acts::GeometryIdentifier geoId = aPlaneLayer->geometryId();
    mapGeoId2DetId[geoId] = detId;
  }

  ///////////// trackfind conf
  auto trackFindFun = TelActs::makeTrackFinderFunction(worldGeo, magneticField);

  std::vector<Acts::CKFSourceLinkSelector::Config::InputElement> ckfConfigEle_vec;
  for(auto& [detId, aPlaneLayer]: mapDetId2PlaneLayer){
    if(aPlaneLayer==layerDets[0]){// x sorted layers
      ckfConfigEle_vec.push_back({aPlaneLayer->geometryId(), {10,10}}); //first layer can have multi-branches
    }
    else{
      ckfConfigEle_vec.push_back({aPlaneLayer->geometryId(), {10,1}}); //chi2, max branches
    }
  }
  Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg(ckfConfigEle_vec);

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Vector3D{-4_m, 0., 0.}, Acts::Vector3D{1., 0., 0.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;
  pOptions.mass = particleMass;

  auto kfLogger = Acts::getDefaultLogger("CKF", logLevel);
  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector> ckfOptions(
    gctx, mctx, cctx, sourcelinkSelectorCfg, Acts::LoggerWrapper{*kfLogger}, pOptions,
    refSurface.get());

  altel::TelEventTTreeWriter ttreeWriter;
  TTree *pTree = new TTree("eventTree", "eventTree");
  ttreeWriter.setTTree(pTree);

  TelFW telfw(800, 400, "test");
  glfw_test telfwtest(geometryFilePath);
  telfw.startAsync<glfw_test>(&telfwtest, &glfw_test::beginHook, &glfw_test::clearHook, &glfw_test::drawHook);

  uint32_t rawFileNum=0;
  size_t readEventNum = 0;
  size_t emptyEventNum = 0;
  size_t eventNum = 0;
  size_t trackNum = 0;
  size_t droppedTrackNum = 0;
  auto tp_start = std::chrono::system_clock::now();

  eudaq::FileReaderUP reader;

  while(1){
    if(eventNum> eventMaxNum && eventMaxNum>0){
      break;
    }
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
    // std::cout<<"reading"<<std::endl;

    readEventNum++;
    if(readEventNum<=eventSkipNum){
      continue;
    }
    eventNum++;

    std::shared_ptr<altel::TelEvent> fullEvent  = altel::createTelEvent(eudaqEvent);
    // TODO test nullptr
    std::shared_ptr<altel::TelEvent> detEvent(new altel::TelEvent(fullEvent->runN(),
                                                                  fullEvent->eveN(),
                                                                  fullEvent->detN(),
                                                                  fullEvent->clkN()));
    detEvent->measHits()=fullEvent->measHits(detId_dets);
    // std::shared_ptr<altel::TelEvent> detEvent  = TelActs::createTelEvent(evpack, runN, eventNum, setupN, mapDetId2PlaneLayer_dets);
    std::vector<TelActs::TelSourceLink> sourcelinks  = TelActs::createSourceLinks(detEvent, mapDetId2PlaneLayer_dets);

    if(sourcelinks.empty()) {
      emptyEventNum ++;
      eventNum ++;
      // std::cout<<"empty"<<std::endl;
      continue;
    }

    ////////////////////////////////
    auto result = trackFindFun(sourcelinks, seedParameters, ckfOptions);
    if (!result.ok()){
      std::fprintf(stderr, "Track finding failed in Event<%lu> , with error \n",
                   eventNum, result.error().message().c_str());
      throw;
    }

    TelActs::fillTelTrajectories(gctx, result.value(), detEvent, mapGeoId2DetId);

    std::shared_ptr<altel::TelEvent> targetEvent(new altel::TelEvent(fullEvent->runN(),
                                                                     fullEvent->eveN(),
                                                                     fullEvent->detN(),
                                                                     fullEvent->clkN()));
    targetEvent->measHits()=fullEvent->measHits(detId_targets);

    TelActs::mergeAndMatchExtraTelEvent(detEvent, targetEvent, 400_um, 3);

    for(auto &aTraj: detEvent->TJs){
      size_t fittedHitNum = aTraj->numOriginMeasHit();
      if(fittedHitNum<3){
        if(fittedHitNum != 1)
          droppedTrackNum++;
        continue;
      }
      trackNum ++;
    }

    ttreeWriter.fillTelEvent(detEvent);
    eventNum ++;

    telfwtest.pushBufferEvent(detEvent);

    if(do_wait){
      std::cout<<"waiting, press any key to conitnue"<<std::endl;
      std::getc(stdin);
    }
  }

  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;
  double time_s = dur_diff.count();
  std::fprintf(stdout, "total time: %.6fs, \nprocessed total %d events,  include %d empty events,\nfound %d good tracks, dropped %d tracks\n",
               time_s, eventNum, emptyEventNum, trackNum, droppedTrackNum);
  std::fprintf(stdout, "event rate: %.0fhz, non-empty event rate: %.0fhz, empty event rate: %.0fhz, track rate: %.0fhz\n",
               eventNum/time_s, (eventNum-emptyEventNum)/time_s, emptyEventNum/time_s, trackNum/time_s);

  TFile tfile(rootFilePath.c_str(),"recreate");
  pTree->Write();
  tfile.Close();

  if(do_wait){
    std::cout<<"waiting, press any key to conitnue"<<std::endl;
    std::getc(stdin);
  }

  telfw.stopAsync();
  return 0;
}
