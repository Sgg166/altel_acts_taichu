#include "TelMultiTrajectory.hpp"
#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

#include <numeric>
#include <chrono>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2.h>
#include <TImage.h>
#include <TCanvas.h>


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
  -dutId          <int>...     DUT IDs which are complectely excluded from track fitting. Residual are caculated.

examples:
./bin/TelTrackResidual -hitFile /work/data/TB2006/alpide_200629033515.json  -geo /work/testbeam/altel_align/runspace/test313/align_313_geo.json -r actsfit.root -even 10
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
                                {"dutId", required_argument, NULL, 'd'},
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
  std::vector<std::shared_ptr<TelActs::TelElement>> eleDets;
  for(const auto& js_det: js_dets.GetArray()){
    auto ele = std::make_shared<TelActs::TelElement>(js_det);
    eleDets.push_back(ele);
  }
  if(eleDets.size()<3){
    std::fprintf(stdout, "error: number of detector elements is only %d.", eleDets.size());
    throw;
  }

  std::shared_ptr<TelActs::TelElement> firstEleDet;
  std::vector<Acts::LayerPtr> layerDets;
  for(auto &ele: eleDets){
    layerDets.push_back(ele->layer());
  }
  Acts::GeometryObjectSorterT<Acts::LayerPtr> layerSorter(gctx, Acts::BinningValue::binX);
  std::sort(layerDets.begin(), layerDets.end(), layerSorter);
  for(auto &ele: eleDets){
    if(ele->layer() == layerDets[0] ){
      firstEleDet = ele;
    }
  }
  if(!firstEleDet){
    std::fprintf(stderr, "error: assign first element");
  }

  const auto &js_targets = js_geo["targets"];
  std::vector<std::shared_ptr<TelActs::TelElement>> eleTargets;
  for(const auto& js_target: js_targets.GetArray()){
    auto ele = std::make_shared<TelActs::TelElement>(js_target);
    eleTargets.push_back(ele);
  }

  std::vector<std::shared_ptr<TelActs::TelElement>> eleDetsAndTargets;
  for(auto &e: eleDets){
    eleDetsAndTargets.push_back(e);
  }
  for(auto &e: eleTargets){
    eleDetsAndTargets.push_back(e);
  }

  std::fprintf(stdout, "elementN = %d, detN = %d, targetN = %d\n",
               eleDetsAndTargets.size(), eleDets.size(), eleTargets.size());

  std::shared_ptr<const Acts::TrackingGeometry> worldGeo =
    TelActs::TelElement::buildWorld(gctx, 11.0_m, 0.1_m, 0.1_m,  eleDetsAndTargets); //include collimator at -5_m
  std::map<Acts::GeometryIdentifier, size_t> mapSurId2DetId;
  std::map<size_t, Acts::GeometryIdentifier> mapDetId2SurId;
  for(auto &ele: eleDetsAndTargets){
    size_t detId = ele->id();
    Acts::GeometryIdentifier surId = ele->surface().geometryId();
    mapDetId2SurId[detId] = surId;
    mapSurId2DetId[surId] = detId;
  }


  std::shared_ptr<const Acts::Surface> firstSurfaceDet = firstEleDet->surface().getSharedPtr();

  ///////////// trackfind conf
  auto trackFindFun = TelActs::makeTrackFinderFunction(worldGeo, magneticField);

  std::vector<Acts::CKFSourceLinkSelector::Config::InputElement> ckfConfigEle_vec;
  for(auto& e: eleDetsAndTargets){
    if(e==firstEleDet){
      ckfConfigEle_vec.push_back({e->surface().geometryId(), {20,10}}); //first layer can have multi-branches
    }
    else{
      ckfConfigEle_vec.push_back({e->surface().geometryId(), {20,1}}); //chi2, max branches
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

//////////////////////////////////////////////////
  uint32_t rRunN;
  uint32_t rEventN;
  int16_t rConfigN;
  int16_t rNumTraj_PerEvent;

  ///rawMeas
  std::vector<int16_t> rRawMeasVec_DetN, *pRawMeasVec_DetN = &rRawMeasVec_DetN;
  std::vector<int16_t> rRawMeasVec_U, *pRawMeasVec_U = &rRawMeasVec_U;
  std::vector<int16_t> rRawMeasVec_V, *pRawMeasVec_V = &rRawMeasVec_V;
  std::vector<int16_t> rRawMeasVec_Clk, *pRawMeasVec_Clk = &rRawMeasVec_Clk;

  //hitMeas
  std::vector<int16_t> rHitMeasVec_DetN, *pHitMeasVec_DetN = &rHitMeasVec_DetN;
  std::vector<double>  rHitMeasVec_U, *pHitMeasVec_U = &rHitMeasVec_U;
  std::vector<double>  rHitMeasVec_V, *pHitMeasVec_V = &rHitMeasVec_V;
  std::vector<int16_t> rHitMeasVec_NumRawMeas_PerHitMeas,
    *pHitMeasVec_NumRawMeas_PerHitMeas = &rHitMeasVec_NumRawMeas_PerHitMeas;
  std::vector<int16_t> rHitMeasVec_Index_To_RawMeas,
    *pHitMeasVec_Index_To_RawMeas = &rHitMeasVec_Index_To_RawMeas;

  //hitFit
  std::vector<int16_t> rHitFitVec_TrajN, *pHitFitVec_TrajN = &rHitFitVec_TrajN;
  std::vector<int16_t> rHitFitVec_DetN, *pHitFitVec_DetN = &rHitFitVec_DetN;
  std::vector<double> rHitFitVec_U, *pHitFitVec_U = &rHitFitVec_U;
  std::vector<double> rHitFitVec_V, *pHitFitVec_V = &rHitFitVec_V;
  std::vector<double> rHitFitVec_X, *pHitFitVec_X = &rHitFitVec_X;
  std::vector<double> rHitFitVec_Y, *pHitFitVec_Y = &rHitFitVec_Y;
  std::vector<double> rHitFitVec_Z, *pHitFitVec_Z = &rHitFitVec_Z;
  std::vector<double> rHitFitVec_DirX, *pHitFitVec_DirX = &rHitFitVec_DirX;
  std::vector<double> rHitFitVec_DirY, *pHitFitVec_DirY = &rHitFitVec_DirY;
  std::vector<double> rHitFitVec_DirZ, *pHitFitVec_DirZ = &rHitFitVec_DirZ;

  std::vector<int16_t> rHitFitVec_Index_To_Origin_HitMeas,
    *pHitFitVec_Index_To_Origin_HitMeas = &rHitFitVec_Index_To_Origin_HitMeas;

  std::vector<int16_t> rHitFitVec_Index_To_Matched_HitMeas,
    *pHitFitVec_Index_To_Matched_HitMeas = &rHitFitVec_Index_To_Matched_HitMeas;

  //traj
  std::vector<int16_t> rTrajVec_NumHitFit_PerTraj,
    *pTrajVec_NumHitFit_PerTraj = &rTrajVec_NumHitFit_PerTraj;
  std::vector<int16_t> rTrajVec_NumHitMeas_Origin_PerTraj,
    *pTrajVec_NumHitMeas_Origin_PerTraj = &rTrajVec_NumHitMeas_Origin_PerTraj;
  std::vector<int16_t> rTrajVec_NumHitMeas_Matched_PerTraj,
    *pTrajVec_NumHitMeas_Matched_PerTraj = &rTrajVec_NumHitMeas_Matched_PerTraj;
  std::vector<int16_t> rTrajVec_Index_To_HitFit,
    *pTrajVec_Index_To_HitFit = &rTrajVec_Index_To_HitFit;


  //ana
  std::vector<int16_t> rAnaVec_Matched_DetN, *pAnaVec_Matched_DetN = &rAnaVec_Matched_DetN;
  std::vector<double> rAnaVec_Matched_ResdU, *pAnaVec_Matched_ResdU = &rAnaVec_Matched_ResdU;
  std::vector<double> rAnaVec_Matched_ResdV, *pAnaVec_Matched_ResdV = &rAnaVec_Matched_ResdV;


  //////////////////////////////////////////////////////////
  TTree tree("eventTree", "eventTree");

  auto bRunN = tree.Branch("RunN", &rRunN);
  auto bEventN = tree.Branch("EventN", &rEventN);
  auto bConfigN = tree.Branch("ConfigN", &rConfigN);
  auto bNumTraj_PerEvent = tree.Branch("NumTraj_PerEvent", &rNumTraj_PerEvent);

  auto bRawMeasVec_DetN = tree.Branch("RawMeasVec_DetN", &pRawMeasVec_DetN);
  auto bRawMeasVec_U = tree.Branch("RawMeasVec_U", &pRawMeasVec_U);
  auto bRawMeasVec_V = tree.Branch("RawMeasVec_V", &pRawMeasVec_V);
  auto bRawMeasVec_Clk = tree.Branch("RawMeasVec_Clk", &pRawMeasVec_Clk);

  auto bHitMeasVec_DetN = tree.Branch("HitMeasVec_DetN", &pHitMeasVec_DetN);
  auto bHitMeasVec_U = tree.Branch("HitMeasVec_U", &pHitMeasVec_U);
  auto bHitMeasVec_V = tree.Branch("HitMeasVec_V", &pHitMeasVec_V);
  auto bHitMeasVec_NumRawMeas_PerHitMeas =
    tree.Branch("HitMeasVec_NumRawMeas_PerHitMeas", &pHitMeasVec_NumRawMeas_PerHitMeas);
  auto bHitMeasVec_Index_To_RawMeas =
    tree.Branch("HitMeasVec_Index_To_RawMeas", &pHitMeasVec_Index_To_RawMeas);

  auto bHitFitVec_DetN = tree.Branch("HitFitVec_DetN", &pHitFitVec_DetN);
  auto bHitFitVec_U = tree.Branch("HitFitVec_U", &pHitFitVec_U);
  auto bHitFitVec_V = tree.Branch("HitFitVec_V", &pHitFitVec_V);
  auto bHitFitVec_X = tree.Branch("HitFitVec_X", &pHitFitVec_X);
  auto bHitFitVec_Y = tree.Branch("HitFitVec_Y", &pHitFitVec_Y);
  auto bHitFitVec_Z = tree.Branch("HitFitVec_Z", &pHitFitVec_Z);
  auto bHitFitVec_DirX = tree.Branch("HitFitVec_DirX", &pHitFitVec_DirX);
  auto bHitFitVec_DirY = tree.Branch("HitFitVec_DirY", &pHitFitVec_DirY);
  auto bHitFitVec_DirZ = tree.Branch("HitFitVec_DirZ", &pHitFitVec_DirZ);
  auto bHitFitVec_Index_To_Origin_HitMeas =
    tree.Branch("HitFitVec_Index_To_Origin_HitMeas", &pHitFitVec_Index_To_Origin_HitMeas);
  auto bHitFitVec_Index_To_Matched_HitMeas =
    tree.Branch("HitFitVec_Index_To_Matched_HitMeas", &pHitFitVec_Index_To_Matched_HitMeas);

  auto bTrajVec_NumHitFit_PerTraj = tree.Branch("TrajVec_NumHitFit_PerTraj", &pTrajVec_NumHitFit_PerTraj);
  auto bTrajVec_NumHitMeas_Origin_PerTraj = tree.Branch("TrajVec_NumHitMeas_Origin_PerTraj", &pTrajVec_NumHitMeas_Origin_PerTraj);
  auto bTrajVec_NumHitMeas_Matched_PerTraj = tree.Branch("TrajVec_NumHitMeas_Matched_PerTraj", &pTrajVec_NumHitMeas_Matched_PerTraj);

  auto bTrajVec_Index_To_HitFit =
    tree.Branch("TrajVec_Index_To_HitFit", &pTrajVec_Index_To_HitFit);


  // ana
  auto bAnaVec_Matched_DetN = tree.Branch("AnaVec_Matched_DetN", &pAnaVec_Matched_DetN);
  auto bAnaVec_Matched_ResdU = tree.Branch("AnaVec_Matched_ResdU", &pAnaVec_Matched_ResdU);
  auto bAnaVec_Matched_ResdV = tree.Branch("AnaVec_Matched_ResdV", &pAnaVec_Matched_ResdV);

/////////////////////////////////////////////////

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

    auto sourcelinks = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleDets);
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
    std::unique_ptr<TelActs::TelMultiTrajectory> mj;
    auto result = trackFindFun(sourcelinks, seedParameters, ckfOptions);
    if (!result.ok()){
      std::fprintf(stderr, "Track finding failed in Event<%lu> , with error \n",
                   eventNum, result.error().message().c_str());
      throw;
    }

    const auto &trackFindingOutput = result.value();
    mj.reset(new TelActs::TelMultiTrajectory(
        std::move(trackFindingOutput.fittedStates),
        std::move(trackFindingOutput.trackTips),
        std::move(trackFindingOutput.fittedParameters)));

    auto sourcelinksTargets = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleTargets);


    auto telEvent = mj->createTelEvent(gctx, mapSurId2DetId, 0, eventNum, 0);

///////////////branch vector clear////////////////////////
    ///rawMeas
    rRawMeasVec_DetN.clear();
    rRawMeasVec_U.clear();
    rRawMeasVec_V.clear();
    rRawMeasVec_Clk.clear();
    //hitMeas
    rHitMeasVec_DetN.clear();
    rHitMeasVec_U.clear();
    rHitMeasVec_V.clear();
    rHitMeasVec_NumRawMeas_PerHitMeas.clear();
    rHitMeasVec_Index_To_RawMeas.clear();

    //hitFit
    rHitFitVec_TrajN.clear();
    rHitFitVec_DetN.clear();
    rHitFitVec_U.clear();
    rHitFitVec_V.clear();
    rHitFitVec_X.clear();
    rHitFitVec_Y.clear();
    rHitFitVec_Z.clear();
    rHitFitVec_DirX.clear();
    rHitFitVec_DirY.clear();
    rHitFitVec_DirZ.clear();
    rHitFitVec_Index_To_Origin_HitMeas.clear();
    rHitFitVec_Index_To_Matched_HitMeas.clear();

    //traj
    rTrajVec_NumHitFit_PerTraj.clear();
    rTrajVec_NumHitMeas_Origin_PerTraj.clear();
    rTrajVec_NumHitMeas_Matched_PerTraj.clear();
    rTrajVec_Index_To_HitFit.clear();


    //matched analysis
    rAnaVec_Matched_DetN.clear();
    rAnaVec_Matched_ResdU.clear();
    rAnaVec_Matched_ResdV.clear();

///////////////////////////////////////////////



    for(auto &telTraj: telEvent->Ts){
      size_t fittedHitNum = telTraj->numberHitFitByMeas();
      if(fittedHitNum<5){
        // std::cout<<"skip telTraj with meas "<< telTraj->numberHitFitByMeas()<< std::endl;
        if(fittedHitNum != 1)
          droppedTrackNum++;
        continue;
      }

      for(const auto &sl : sourcelinksTargets){
        size_t id = mapSurId2DetId.at(sl.referenceSurface().geometryId());

        auto aTelHit = telTraj->hit(id);
        if(!aTelHit || !aTelHit->hasHitFit() || aTelHit->isFittedFromMeasure()){
          continue;
        }

        Acts::Vector2D xy_fit(aTelHit->HF->PLs[0], aTelHit->HF->PLs[1]);
        Acts::Vector2D xy_meas = sl.value();
        Acts::Vector2D xy_resid = xy_meas - xy_fit;
        if(xy_resid.norm()>0.1_mm){
          continue;
        }

        aTelHit->HM.reset(new TelActs::TelHitMeasure{id, {xy_meas(0), xy_meas(1)}, {}});
      }
      trackNum ++;
    }

////////////////////////////////////////fill tree////////////////////
    std::map<std::shared_ptr<TelActs::TelRawMeasure>, int16_t> poolMapRawMeas;
    std::map<std::shared_ptr<TelActs::TelHitMeasure>, int16_t> poolMapHitMeas;
    std::map<std::shared_ptr<TelActs::TelHitFit>, int16_t> poolMapHitFit;

    for(auto &aRawMeas: telEvent->Ms){
      auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()) );
      if(inserted){
        // inserted and return true, when no exist
        rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
        rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
        rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
        rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
      }
    }

    for(auto &aHitMeas: telEvent->HMs){
      auto [it, inserted] = poolMapHitMeas.emplace( aHitMeas, int16_t(rHitMeasVec_DetN.size()) );
      if(inserted){
        // inserted and return true, when no exist
        rHitMeasVec_DetN.push_back(aHitMeas->DN);
        rHitMeasVec_U.push_back(aHitMeas->PLs[0]);
        rHitMeasVec_V.push_back(aHitMeas->PLs[1]);

        rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
        for(auto &aRawMeas: aHitMeas->Ms){
          auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
          if(inserted){
            rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
            rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
            rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
            rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
          }
          rHitMeasVec_Index_To_RawMeas.push_back(it->second);
        }
        rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeas->Ms.size()));
      }
    }

    rNumTraj_PerEvent = 0;
    for(auto &aTraj: telEvent->Ts){
      size_t fittedHitNum = aTraj->numberHitFitByMeas();
      if(fittedHitNum<5){
        continue;
      }

      rNumTraj_PerEvent ++;
      rTrajVec_Index_To_HitFit.push_back(-1);//traj: -1 {M-N} -1 {M-N}

      int16_t aNumHitFit_PerTraj = 0;
      int16_t aNumHitMeas_Origin_PerTraj = 0;
      int16_t aNumHitMeas_Matched_PerTraj = 0;

      for(auto &aHit: aTraj->Hs){
        aNumHitFit_PerTraj++;
        auto aHitFit = aHit->HF;
        auto [it, inserted] = poolMapHitFit.emplace( aHitFit, int16_t(rHitFitVec_DetN.size()) );
        if(!inserted){
          //TODO: handle this case for ambiguility case
          std::fprintf(stderr, "a shared hitfit by multi-traj, WRONG\n");
          throw;
        }

        rTrajVec_Index_To_HitFit.push_back(it->second);
        rHitFitVec_TrajN.push_back(int16_t(rTrajVec_NumHitFit_PerTraj.size()));
        rHitFitVec_DetN.push_back(aHitFit->DN);
        rHitFitVec_U.push_back(aHitFit->PLs[0]);
        rHitFitVec_V.push_back(aHitFit->PLs[1]);
        rHitFitVec_X.push_back(aHitFit->PGs[0]);
        rHitFitVec_Y.push_back(aHitFit->PGs[1]);
        rHitFitVec_Z.push_back(aHitFit->PGs[2]);
        rHitFitVec_DirX.push_back(aHitFit->DGs[0]);
        rHitFitVec_DirY.push_back(aHitFit->DGs[1]);
        rHitFitVec_DirZ.push_back(aHitFit->DGs[2]);

        int16_t indexHitMeasOri = -1;
        auto aHitMeasOri = aHitFit->HM;
        if(aHitMeasOri){
          aNumHitMeas_Origin_PerTraj ++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasOri, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasOri->DN);
            rHitMeasVec_U.push_back(aHitMeasOri->PLs[0]);
            rHitMeasVec_V.push_back(aHitMeasOri->PLs[1]);

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasOri->Ms){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
                rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
                rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
                rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasOri->Ms.size()));
          }
          indexHitMeasOri = it->second;
        }
        rHitFitVec_Index_To_Origin_HitMeas.push_back(indexHitMeasOri);

        int16_t  indexHitMeasMatched = -1;
        auto aHitMeasMatched = aHit->HM;
        if(aHitMeasMatched){
          aNumHitMeas_Matched_PerTraj++;
          auto [it, inserted] = poolMapHitMeas.emplace( aHitMeasMatched, int16_t(rHitMeasVec_DetN.size()));
          if(inserted){
            rHitMeasVec_DetN.push_back(aHitMeasMatched->DN);
            rHitMeasVec_U.push_back(aHitMeasMatched->PLs[0]);
            rHitMeasVec_V.push_back(aHitMeasMatched->PLs[1]);

            rHitMeasVec_Index_To_RawMeas.push_back(-1);//: -1 {M-N} -1 {M-N}
            for(auto &aRawMeas: aHitMeasMatched->Ms){
              auto [it, inserted] = poolMapRawMeas.emplace(aRawMeas, int16_t(rRawMeasVec_DetN.size()));
              if(inserted){
                rRawMeasVec_U.push_back(aRawMeas->uvdcS[0]);
                rRawMeasVec_V.push_back(aRawMeas->uvdcS[1]);
                rRawMeasVec_DetN.push_back(aRawMeas->uvdcS[2]);
                rRawMeasVec_Clk.push_back(aRawMeas->uvdcS[3]);
              }
              rHitMeasVec_Index_To_RawMeas.push_back(it->second);
            }
            rHitMeasVec_NumRawMeas_PerHitMeas.push_back(int16_t(aHitMeasMatched->Ms.size()));
            //

          }
          indexHitMeasMatched = it->second;
        }
        rHitFitVec_Index_To_Matched_HitMeas.push_back(indexHitMeasMatched);

        //ana matched
        if(aHitMeasMatched){
          rAnaVec_Matched_DetN.push_back(aHitMeasMatched->DN);
          rAnaVec_Matched_ResdU.push_back(aHitMeasMatched->PLs[0] - aHitFit->PLs[0]);
          rAnaVec_Matched_ResdV.push_back(aHitMeasMatched->PLs[1] - aHitFit->PLs[1]);
        }
      }
      rTrajVec_NumHitFit_PerTraj.push_back(aNumHitFit_PerTraj);
      rTrajVec_NumHitMeas_Origin_PerTraj.push_back(aNumHitMeas_Origin_PerTraj);
      rTrajVec_NumHitMeas_Matched_PerTraj.push_back(aNumHitMeas_Matched_PerTraj);
    }

    rRunN = telEvent->RN;
    rEventN = telEvent->EN;
    rConfigN = telEvent->DN;

    tree.Fill();
////////////////////////////////////////fill tree////////////////////
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
  tree.Write();
  tfile.Close();

  return 0;
}
