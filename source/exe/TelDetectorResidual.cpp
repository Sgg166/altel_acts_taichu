#include "TelEventTTreeWriter.hpp"
#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

#include <numeric>
#include <chrono>

#include <TFile.h>
#include <TTree.h>

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                        help message
  -verbose                     verbose flag
  -eventSkip      <int>        number of events to skip before start processing
  -eventMax       <int>        max number of events to process
  -geometryFile   <path>       path to geometry input file (input)
  -hitFile        <path>       path data input file (input)
  -rootFile       <path>       path to root file (output)
  -particleEnergy <float>      energy of beam particle, electron, (Gev)
  -targetId       <int>...     IDs of target detector which are complectely excluded from track fitting. Residual are caculated.

examples:
./bin/TelDetectorResidual -hitFile /work/data/TB2006/alpide_200629033515.json  -geo /work/testbeam/altel_align/runspace/test313/align_313_geo.json -r detresid.root -eventM 10
)";


int main(int argc, char *argv[]) {
  int64_t eventMaxNum = -1;
  int64_t eventSkipNum = 0;
  std::vector<int64_t> dutIDs;
  std::string hitFilePath;
  std::string geometryFilePath;
  std::string rootFilePath;

  double particleQ = 1;
  double particleMass = 0.511 * Acts::UnitConstants::MeV;
  double particleEnergy = 5.0 * Acts::UnitConstants::GeV;

  double distCollimator = 5_m;
  double widthCollimator = 4_cm;

  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"eventSkip", required_argument, NULL, 's'},
                                {"eventMax", required_argument, NULL, 'm'},
                                {"hitFile", required_argument, NULL, 'f'},
                                {"rootFile", required_argument, NULL, 'b'},
                                {"geometryFile", required_argument, NULL, 'g'},
                                {"particleEnergy", required_argument, NULL, 'e'},
                                {"targetId", required_argument, NULL, 'd'},
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
      case 'f':
        hitFilePath = optarg;
        break;
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
          dutIDs.push_back(std::stol(argv[i]));
          optind = i+1;
        }
        break;
      }
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
  std::fprintf(stdout, "hitFileFile:   %s\n", hitFilePath.c_str());
  std::fprintf(stdout, "geometryFile:  %s\n", geometryFilePath.c_str());
  std::fprintf(stdout, "rootFile:      %s\n", rootFilePath.c_str());
  std::fprintf(stdout, "\n");

  if (hitFilePath.empty() || rootFilePath.empty() ||
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
  JsonUtils::printJsonValue(jsd_geo, true);

  std::printf("--------create acts geo object-----\n");
  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }
  const auto &js_geo = jsd_geo["geometry"];
  const auto &js_dets = js_geo["detectors"];

  std::map<size_t, std::shared_ptr<const Acts::PlaneLayer>> mapDetId2PlaneLayer;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> layerDets;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> layerTargets;
  std::vector<std::shared_ptr<const Acts::PlaneLayer>> allPlaneLayers;

  for(const auto& js_det: js_dets.GetArray()){
    auto [detId, planeLayer] = TelActs::createPlaneLayer(js_det);
    layerDets.push_back(planeLayer);
    mapDetId2PlaneLayer[detId] = planeLayer;
  }
  if(layerDets.size()<3){
    std::fprintf(stdout, "error: number of detector elements is only %d.", layerDets.size());
    throw;
  }

  Acts::GeometryObjectSorterT<std::shared_ptr<const Acts::PlaneLayer>> layerSorter(gctx, Acts::BinningValue::binX);
  std::sort(layerDets.begin(), layerDets.end(), layerSorter);

  const auto &js_targets = js_geo["targets"];
  for(const auto& js_target: js_targets.GetArray()){
    auto [detId, planeLayer] = TelActs::createPlaneLayer(js_target);
    layerTargets.push_back(planeLayer);
    mapDetId2PlaneLayer[detId] = planeLayer;
  }

  std::fprintf(stdout, "elementN = %d, detN = %d, targetN = %d\n",
               mapDetId2PlaneLayer.size(), layerDets.size(), layerTargets.size());

  for(auto& [detId, aPlaneLayer]: mapDetId2PlaneLayer){
    allPlaneLayers.push_back(aPlaneLayer);
  }

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
      ckfConfigEle_vec.push_back({aPlaneLayer->geometryId(), {20,10}}); //first layer can have multi-branches
    }
    else{
      ckfConfigEle_vec.push_back({aPlaneLayer->geometryId(), {20,1}}); //chi2, max branches
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

  TelActs::TelEventTTreeWriter ttreeWriter;
  ttreeWriter.pTree.reset(new TTree("eventTree", "eventTree"));
  ttreeWriter.createBranch();

  JsonFileDeserializer jsfd(hitFilePath);

  for(size_t i=0; i< eventSkipNum; i++){
    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object after skip %d event, possible end of file\n", i);
    }
  }

  size_t emptyEventNum = 0;
  size_t eventNum = 0;
  size_t trackNum = 0;
  size_t droppedTrackNum = 0;
  auto tp_start = std::chrono::system_clock::now();
  while(jsfd && (eventNum< eventMaxNum || eventMaxNum<0)){

    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, possible end of file\n");
      break;
    }

    size_t runN = 0;
    size_t setupN = 0;
    std::shared_ptr<TelActs::TelEvent> detEvent  = TelActs::createTelEvent(evpack, layerDets,  mapGeoId2DetId, runN, eventNum, setupN);
    std::vector<TelActs::TelSourceLink> sourcelinks  = TelActs::createSourceLink(mapGeoId2DetId, layerDets, detEvent);

    if(sourcelinks.empty()) {
      emptyEventNum ++;
      eventNum ++;
      continue;
    }

    if(do_verbose){
      std::fprintf(stdout, "\n\n\n");
      for(const auto& l: evpack["layers"].GetArray()){
        JsonUtils::printJsonValue(l, false);
      }
    }

    ////////////////////////////////
    auto result = trackFindFun(sourcelinks, seedParameters, ckfOptions);
    if (!result.ok()){
      std::fprintf(stderr, "Track finding failed in Event<%lu> , with error \n",
                   eventNum, result.error().message().c_str());
      throw;
    }

    std::shared_ptr<TelActs::TelEvent> targetEvent  = TelActs::createTelEvent(evpack, layerTargets,  mapGeoId2DetId,
                                                                              runN, eventNum, setupN);
    std::vector<TelActs::TelSourceLink> sourcelinksTargets  = TelActs::createSourceLink(mapGeoId2DetId, layerTargets, targetEvent);

    std::shared_ptr<TelActs::TelEvent> telEvent = TelActs::createTelEvent(gctx, result.value(), mapGeoId2DetId,
                                                                          runN, eventNum, setupN);

    for(auto &telTraj: telEvent->Ts){
      size_t fittedHitNum = telTraj->numberHitFitByMeas();
      if(fittedHitNum<5){
        if(fittedHitNum != 1)
          droppedTrackNum++;
        continue;
      }
      trackNum ++;
    }
    TelActs::matchAddExtraHitMeas(gctx, mapGeoId2DetId, telEvent, sourcelinksTargets);
    // todo: match targetEvent directly
    ttreeWriter.fill(telEvent);
    eventNum ++;
  }

  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;
  double time_s = dur_diff.count();
  std::fprintf(stdout, "total time: %.6fs, \nprocessed total %d events,  include %d empty events,\nfound %d good tracks, dropped %d tracks\n",
               time_s, eventNum, emptyEventNum, trackNum, droppedTrackNum);
  std::fprintf(stdout, "event rate: %.0fhz, non-empty event rate: %.0fhz, empty event rate: %.0fhz, track rate: %.0fhz\n",
               eventNum/time_s, (eventNum-emptyEventNum)/time_s, emptyEventNum/time_s, trackNum/time_s);

  TFile tfile(rootFilePath.c_str(),"recreate");
  ttreeWriter.pTree->Write();
  tfile.Close();

  return 0;
}
