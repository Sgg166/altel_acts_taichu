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
  -eventMax       [int]        max number of events to process
  -geometryFile   [PATH]       path to geometry input file (input)
  -hitFile        [PATH]       path data input file (input)
  -rootFile       [PATH]       path to root file (output)
  -energy         [float]      beam energy, GeV
  -dutID          [int]        ID of DUT which is excluded from track fit

examples:
./bin/TelTrackResidual -hitFile /work/data/TB2006/alpide_200629033515.json  -geo /work/testbeam/altel_align/runspace/test313/align_313_geo.json -r actsfit.root -even 10
)";


int main(int argc, char *argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"eventMax", required_argument, NULL, 'm'},
                              {"hitFile", required_argument, NULL, 'f'},
                              {"rootFile", required_argument, NULL, 'b'},
                              {"geomertyFile", required_argument, NULL, 'g'},
                              {"energy", required_argument, NULL, 'e'},
                              {"dutID", required_argument, NULL, 'd'},
                              {0, 0, 0, 0}};

  int64_t eventMaxNum = -1;
  int64_t dutID = -1;
  std::string dataFile_path;
  std::string rootFile_path;
  std::string geoFile_path;
  double beamEnergy = -1;

  double seedResX = 5_um;
  double seedResY = 5_um;
  double seedResPhi = 0.03_rad;
  double seedResTheta = 0.03_rad;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'm':
      eventMaxNum = std::stoul(optarg);
      break;
    case 'f':
      dataFile_path = optarg;
      break;
    case 'b':
      rootFile_path = optarg;
      break;
    case 'g':
      geoFile_path = optarg;
      break;
    case 'e':
      beamEnergy = std::stod(optarg) * Acts::UnitConstants::GeV;
      break;
    case 'd':
      dutID = std::stol(optarg);
      break;

      /////generic part below///////////
    case 0:
      //getopt_long() set a variable, just keep going
      break;
    case '?':
      //internal unrecognized/ambiguous
      std::exit(1);
      break;
    case 1:
      fprintf(stderr, "case 1\n");
      std::exit(1);
      break;
    case ':':
      fprintf(stderr, "case :\n");
      std::exit(1);
      break;
    default:
      fprintf(stderr, "case default, missing branch in getopt\n");
      std::exit(1);
      break;
    }
  }

  if(do_help){
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    std::exit(0);
  }

  if (beamEnergy < 0) {
    beamEnergy = 5.0 * Acts::UnitConstants::GeV;
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "dataFile:         %s\n", dataFile_path.c_str());
  std::fprintf(stdout, "geoFile:          %s\n", geoFile_path.c_str());
  std::fprintf(stdout, "rootFile:         %s\n", rootFile_path.c_str());
  std::fprintf(stdout, "seedResPhi:       %f\n", seedResPhi);
  std::fprintf(stdout, "seedResTheta:     %f\n", seedResTheta);
  std::fprintf(stdout, "\n");

  if (dataFile_path.empty() || rootFile_path.empty() ||
      geoFile_path.empty()) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  /////////////////////////////////////
  Acts::Logging::Level logLevel =
    do_verbose ? (Acts::Logging::VERBOSE) : (Acts::Logging::INFO);

  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  /////////////  track seed conf
  Acts::BoundSymMatrix cov_seed;
  double seedResX2 = seedResX * seedResX;
  double seedResY2 = seedResY * seedResY;
  double seedResPhi2 = seedResPhi * seedResPhi;
  double seedResTheta2 = seedResTheta * seedResTheta;
  cov_seed <<
    seedResX2,0.,       0.,         0.,           0.,     0.,
    0.,       seedResY2,0.,         0.,           0.,     0.,
    0.,       0.,       seedResPhi2,0.,           0.,     0.,
    0.,       0.,       0.,         seedResTheta2,0.,     0.,
    0.,       0.,       0.,         0.,           0.0001, 0.,
    0.,       0.,       0.,         0.,           0.,     1.;

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geoFile_path);
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geoFile_path.c_str() );
    throw;
  }
  JsonUtils::printJsonValue(jsd_geo, true);

  std::printf("--------create acts geo object-----\n");
  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }
  const auto &js_geo = jsd_geo["geometry"];


  std::vector<std::shared_ptr<TelActs::TelElement>> eleDets;
  const auto &js_dets = js_geo["detectors"];
  for(const auto& js_det: js_dets.GetArray()){
    auto ele = std::make_shared<TelActs::TelElement>(js_det);
    eleDets.push_back(ele);
  }

  std::vector<std::shared_ptr<TelActs::TelElement>> eleTargets;
  const auto &js_targets = js_geo["targets"];
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
    TelActs::TelElement::buildWorld(gctx, 4.0_m, 0.1_m, 0.1_m,  eleDetsAndTargets);

  std::map<Acts::GeometryIdentifier, size_t> mapSurId2DetId;
  std::map<size_t, Acts::GeometryIdentifier> mapDetId2SurId;
  for(auto &ele: eleDetsAndTargets){
    size_t detId = ele->id();
    Acts::GeometryIdentifier surId = ele->surface().geometryId();
    mapDetId2SurId[detId] = surId;
    mapSurId2DetId[surId] = detId;
  }

  //40_mm, 20_mm, 80_um

  // Set up surfaces
  size_t id_seed = 0;
  std::shared_ptr<const Acts::Surface> surface_seed;
  for(auto& e: eleDets){
    if(e->id() == id_seed){
      surface_seed = e->surface().getSharedPtr();
    }
  }

  ///////////// trackfind conf
  auto trackFindFun = TelActs::makeTrackFinderFunction(worldGeo, magneticField);

  std::vector<Acts::CKFSourceLinkSelector::Config::InputElement> ckfConfigEle_vec;
  for(auto& e: eleDetsAndTargets){
    ckfConfigEle_vec.push_back({e->surface().geometryId(), {20,1}}); //chi2, max branches
  }

  Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg(ckfConfigEle_vec);

  //   = Acts::CKFSourceLinkSelector::Config{
  //   // {Acts::GeometryIdentifier(), {20, 1}}
  //   {surface_seed0->geometryIdentifier(), {20, 1}},
  //   {surface_seed0->geometryIdentifier(), {20, 1}},
  //   {surface_seed0->geometryIdentifier(), {20, 1}}
  // }; //max chi2, max number of selected hits

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{0, 0., 1.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;
  pOptions.mass = 0.511 * Acts::UnitConstants::MeV;


  auto kfLogger = Acts::getDefaultLogger("KalmanFilter", logLevel);
  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector> ckfOptions(
    gctx, mctx, cctx, sourcelinkSelectorCfg, Acts::LoggerWrapper{*kfLogger}, pOptions,
    refSurface.get());

///////////////////////////////////////////////////

  TTree tree("actsfit", "actsfit");

  // size_t nEvent;
  // size_t nTracks;

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

  JsonFileDeserializer jsfd(dataFile_path);

  size_t eventNumber = 0;
  size_t trackNumber = 0;
  size_t seedNumber = 0;
  auto tp_start = std::chrono::system_clock::now();
  while(jsfd && (eventNumber< eventMaxNum || eventMaxNum<0)){
    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }

    auto sourcelinks = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleDets);
    if(sourcelinks.empty()) {
      std::fprintf(stdout, "Empty event <%d>.\n", eventNumber);
    }

    std::vector<Acts::CurvilinearTrackParameters> initialParameters;
    for (auto &sl : sourcelinks) {
      if( surface_seed != sl.referenceSurface().getSharedPtr() ){
        continue;
      }
      Acts::Vector3D seedOriginWorld = sl.globalPosition(gctx);
      const double phi = 0;
      const double theta = 0.5*M_PI;
      Acts::Vector4D rPos4(seedOriginWorld.x(), seedOriginWorld.y(), seedOriginWorld.z(), 0);
      double q = 1;
      initialParameters.emplace_back(rPos4, phi, theta, beamEnergy, q, cov_seed);
    }

    ////////////////////////////////
    std::vector<TelActs::TelMultiTrajectory> trajectories;
    // Loop ever the seeds
    size_t iseed = 0;
    size_t nTracks = 0;
    for (const auto &rStart : initialParameters) {
      auto result = trackFindFun(sourcelinks, rStart, ckfOptions);
      if (result.ok()) {
        // Get the track finding output object
        const auto &trackFindingOutput = result.value();
        // Create a PixelMultiTrajectory
        nTracks += trackFindingOutput.trackTips.size();
        trajectories.emplace_back(
          std::move(trackFindingOutput.fittedStates),
          std::move(trackFindingOutput.trackTips),
          std::move(trackFindingOutput.fittedParameters));
      }
      else {
        std::printf("Track finding failed in Event<%lu> seed<%lu>, with error \n",
                    eventNumber, iseed, result.error().message().c_str());
        // for(const auto& l: evpack["layers"].GetArray()){
        //   JsonUtils::printJsonValue(l, false);
        // }
      }
      iseed++;
    }
    seedNumber+=iseed;

    auto sourcelinksTargets = TelActs::TelSourceLink::CreateSourceLinks(evpack, eleTargets);
    for(const auto &mj: trajectories){
      if(mj.trackNumber()==0){
        continue;
      }
      mj.fillSingleTrack(gctx, mapSurId2DetId,
                         idMeas, xMeas, yMeas,
                         xResidLocal, yResidLocal,
                         idFit, xFitLocal, yFitLocal,
                         zFitWorld, xFitWorld, yFitWorld,
                         0);

      if(idMeas.size()<4){
        std::cout<<idMeas.size()<< " measus is less than"<< 4<<std::endl;
        continue;
      }
      for(const auto &sl : sourcelinksTargets){
        auto baseEle = sl.referenceSurface().associatedDetectorElement();
        // const TelActs::TelElement* ele = dynamic_cast<const TelActs::TelElement*>(baseEle);
        // if(!ele){
        //   std::cout<< "too wrong"<<std::endl;
        //   continue;
        // }
        // size_t id = ele->id();
        size_t id = mapSurId2DetId.at(sl.referenceSurface().geometryId());


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
      trackNumber ++;
      tree.Fill();
    }

    eventNumber ++;
    //TODO gl
  }
  auto tp_end = std::chrono::system_clock::now();
  std::chrono::duration<double> dur_diff = tp_end-tp_start;
  double time_s = dur_diff.count();
  std::fprintf(stdout, "processed %d events,  %d seeds, %d good tracks\n", eventNumber, seedNumber, trackNumber);
  std::fprintf(stdout, "total time: %fs, event rate: %fhz,  seed rate: %fhz,  track rate: %fhz",
               time_s, eventNumber/time_s, seedNumber/time_s, trackNumber/time_s);

  TFile tfile(rootFile_path.c_str(),"recreate");
  tree.Write();
  tfile.Close();

  return 0;
}
