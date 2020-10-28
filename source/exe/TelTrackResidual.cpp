#include "TelMultiTrajectory.hpp"
#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"


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
  std::string datafile_name;
  std::string rootfile_name;
  std::string geofile_name;
  double beamEnergy = -1;

  double seedResX = 5_um;
  double seedResY = 5_um;
  double seedResPhi = 0.1_rad;
  double seedResTheta = 0.1_rad;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'm':
      eventMaxNum = std::stoul(optarg);
      break;
    case 'f':
      datafile_name = optarg;
      break;
    case 'b':
      rootfile_name = optarg;
      break;
    case 'g':
      geofile_name = optarg;
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
  std::fprintf(stdout, "datafile:         %s\n", datafile_name.c_str());
  std::fprintf(stdout, "geofile:          %s\n", geofile_name.c_str());
  std::fprintf(stdout, "rootFile:         %s\n", rootfile_name.c_str());
  std::fprintf(stdout, "seedResPhi:       %f\n", seedResPhi);
  std::fprintf(stdout, "seedResTheta:     %f\n", seedResTheta);
  std::fprintf(stdout, "\n");

  if (datafile_name.empty() || rootfile_name.empty() ||
      geofile_name.empty()) {
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
  cov_seed << seedResX * seedResX, 0., 0., 0., 0., 0., 0.,
    seedResY * seedResY, 0., 0., 0., 0., 0., 0.,
    seedResPhi * seedResPhi, 0., 0., 0., 0., 0., 0.,
    seedResTheta * seedResTheta, 0., 0., 0., 0., 0., 0.,
    0.0001, 0., 0., 0., 0., 0., 0., 1.;

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geofile_name);
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geofile_name.c_str() );
    throw;
  }
  JsonUtils::printJsonValue(jsd_geo, true);

  std::printf("--------create acts geo object-----\n");
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  std::vector<std::shared_ptr<TelActs::TelElement>> element_col;
  std::tie(trackingGeometry, element_col)  = TelActs::TelElement::buildGeometry(gctx, jsd_geo);
  //40_mm, 20_mm, 80_um

  // Set up surfaces
  size_t id_seed = 0;
  std::shared_ptr<const Acts::Surface> surface_seed;
  for(auto& e: element_col){
    if(e->id() == id_seed){
      surface_seed = e->surface().getSharedPtr();
    }
  }

  ///////////// trackfind conf
  auto trackFindFun = TelActs::makeTrackFinderFunction(trackingGeometry, magneticField);

  Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg = Acts::CKFSourceLinkSelector::Config{
    {Acts::GeometryIdentifier(), {10, 1}}}; //max chi2, max number of selected hits

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
  JsonFileDeserializer jsfd(datafile_name);

  size_t eventNumber = 0;
  while(jsfd && (eventNumber< eventMaxNum || eventMaxNum<0)){
    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }

    auto sourcelinks = TelActs::TelSourceLink::CreateSourceLinks(evpack, element_col);// all element, TODO: select
    if(sourcelinks.empty()) {
      std::fprintf(stdout, "Empty event <%d>.\n", eventNumber);
    }

    std::vector<Acts::CurvilinearTrackParameters> initialParameters;
    for (auto &sl : sourcelinks) {
      if( surface_seed != sl.referenceSurface().getSharedPtr() ){
        continue;
      }
      Acts::Vector3D global = sl.globalPosition(gctx);
      const double phi = 0.0000001;
      const double theta = 0.0000001;
      Acts::Vector4D rPos4(global.x(), global.y(), global.z(), 0);
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
      }
      iseed++;
    }

    for(const auto &mj: trajectories){
      mj.fillSingleTrack(gctx,
                         idMeas, xMeas, yMeas,
                         xResidLocal, yResidLocal,
                         idFit, xFitLocal, yFitLocal,
                         xFitWorld, yFitWorld, zFitWorld);


    }

    eventNumber ++;
    //TODO gl
  }
  return 0;
}
