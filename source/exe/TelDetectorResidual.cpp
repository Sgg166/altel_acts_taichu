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

///////////////////////////////////////////////////

  TTree tree("actsfit", "actsfit");

  size_t nEvent;
  size_t nTracks;

  std::vector<size_t> idMeas;
  std::vector<double> xMeas;
  std::vector<double> yMeas;
  std::vector<double> xResidLocal;
  std::vector<double> yResidLocal;

  std::vector<size_t> idFit;
  std::vector<double> xFitLocal;
  std::vector<double> yFitLocal;
  std::vector<double> xFitWorld;
  std::vector<double> yFitWorld;
  std::vector<double> zFitWorld;

  std::vector<size_t>* p_idMeas = &idMeas;
  std::vector<double>* p_xMeas = &xMeas;
  std::vector<double>* p_yMeas = &yMeas;
  std::vector<double>* p_xResidLocal = &xResidLocal;
  std::vector<double>* p_yResidLocal = &yResidLocal;

  std::vector<size_t>* p_idFit = &idFit;
  std::vector<double>* p_xFitLocal = &xFitLocal;
  std::vector<double>* p_yFitLocal = &yFitLocal;
  std::vector<double>* p_xFitWorld = &xFitWorld;
  std::vector<double>* p_yFitWorld = &yFitWorld;
  std::vector<double>* p_zFitWorld = &zFitWorld;

  auto br_idMeas = tree.Branch("idMeas", &p_idMeas);
  auto br_xMeas = tree.Branch("xMeas", &p_xMeas);
  auto br_YMeas = tree.Branch("yMeas", &p_yMeas);
  auto br_XResidLocal = tree.Branch("xResidLocal", &p_xResidLocal);
  auto br_YResidLocal = tree.Branch("yResidLocal", &p_yResidLocal);

  auto br_idFit = tree.Branch("idFit", &p_idFit);
  auto br_xFitLocal = tree.Branch("xFitLocal", &p_xFitLocal);
  auto br_yFitLocal = tree.Branch("yFitLocal", &p_yFitLocal);
  auto br_xFitWorld = tree.Branch("xFitWorld", &p_xFitWorld);
  auto br_yFitWorld = tree.Branch("yFitWorld", &p_yFitWorld);
  auto br_zFitWorld = tree.Branch("zFitWorld", &p_zFitWorld);

///////////////////////////////////////////////////

  JsonFileDeserializer jsfd(hitFilePath);
  size_t emptyEventNum = 0;
  size_t eventNum = 0;
  size_t trackNum = 0;
  auto tp_start = std::chrono::system_clock::now();
  while(jsfd && (eventNum< eventMaxNum || eventMaxNum<0)){
    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }

    auto sourcelinks = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleDets);
    if(sourcelinks.empty()) {
      emptyEventNum ++;
      eventNum ++;
      continue;
    }

    if(do_verbose){
      for(const auto& l: evpack["layers"].GetArray()){
        JsonUtils::printJsonValue(l, false);
      }
    }

    ////////////////////////////////
    std::unique_ptr<TelActs::TelMultiTrajectory> mj;
    // Loop ever the seeds
    size_t nTracks = 0;
    auto result = trackFindFun(sourcelinks, seedParameters, ckfOptions);
    if (result.ok()) {
      // Get the track finding output object
      const auto &trackFindingOutput = result.value();
      // Create a PixelMultiTrajectory
      nTracks += trackFindingOutput.trackTips.size();

      mj.reset(new TelActs::TelMultiTrajectory(
        std::move(trackFindingOutput.fittedStates),
        std::move(trackFindingOutput.trackTips),
        std::move(trackFindingOutput.fittedParameters)));
    }
    else {
      std::printf("Track finding failed in Event<%lu> , with error \n",
                  eventNum, result.error().message().c_str());
      // for(const auto& l: evpack["layers"].GetArray()){
      //   JsonUtils::printJsonValue(l, false);
      // }
      throw;
    }

    auto sourcelinksTargets = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleTargets);

    size_t nTrackPerEvent= mj?mj->trackNumber():0;
    for(size_t indexTrack=0; indexTrack<nTrackPerEvent; indexTrack++){
      mj->fillSingleTrack(gctx,
                          idMeas, xMeas, yMeas,
                          xResidLocal, yResidLocal,
                          idFit, xFitLocal, yFitLocal,
                          zFitWorld, xFitWorld, yFitWorld,
                          indexTrack);

      if(idMeas.size()<4){
        continue;
      }
      for(const auto &sl : sourcelinksTargets){
        auto baseEle = sl.referenceSurface().associatedDetectorElement();
        const TelActs::TelElement* ele = dynamic_cast<const TelActs::TelElement*>(baseEle);
        if(!ele){
          std::cout<< "too wrong"<<std::endl;
          continue;
        }
        size_t id = ele->id();
        auto it = std::find (idFit.begin(), idFit.end(), id);
        if(it==idFit.end()){
          // std::fprintf(stdout, "not able to find fitted position for id %d\n", id);
          // for(auto &e: idFit){
          //   std::cout<< e<<std::endl;
          // }
          continue;
        }
        size_t n= it-idFit.begin();
        Acts::Vector2D xy_fit(xFitLocal[n], yFitLocal[n]);
        Acts::Vector2D xy_meas = sl.value();
        Acts::Vector2D xy_resid = xy_meas - xy_fit;
        if(xy_resid.norm()>0.1_mm){
          // std::fprintf(stdout, "large distance %f\n,  xy_fit = %f, %f   xy_meas = %f, %f\n",
          //              xy_resid.norm(), xy_fit(0), xy_fit(1), xy_meas(0), xy_meas(1));
          // JsonUtils::printJsonValue(evpack, false);
          continue;
        }
        idMeas.push_back(id);
        xMeas.push_back(xy_meas(0));
        yMeas.push_back(xy_meas(1));
        xResidLocal.push_back(xy_resid(0));
        yResidLocal.push_back(xy_resid(1));
      }
      trackNum ++;
      tree.Fill();
    }

    eventNum ++;
    //TODO gl
  }
  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;
  double time_s = dur_diff.count();
  std::fprintf(stdout, "processed total %d events, include %d empty events, found %d good tracks\n", eventNum, emptyEventNum, trackNum);
  std::fprintf(stdout, "total time: %fs, event rate: %fhz,  non-empty event rate: %fhz,  empty event rate: %fhz,  track rate: %fhz",
               time_s, eventNum/time_s, (eventNum-emptyEventNum)/time_s, emptyEventNum/time_s, trackNum/time_s);

  TFile tfile(rootFilePath.c_str(),"recreate");
  tree.Write();
  tfile.Close();

  return 0;
}
