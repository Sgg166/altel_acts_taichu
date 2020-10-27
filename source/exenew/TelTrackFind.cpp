// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include "TelMultiTrajectory.hpp"

#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                        help message
  -verbose                     verbose flag
  -eventMax       [int]        max number of events to process
  -geometryFile   [PATH]       path to geometry input file (input)
  -hitFile        [PATH]       path data input file (input)
  -trackFile      [PATH]       path to fitted track data file (output)
  -energy         [float]      beam energy, GeV
  -dutID          [int]        ID of DUT which is excluded from track fit
  -seedResPhi     [float]      preset seed track resolution Phi
  -seedResTheta   [float]      preset seed track resolution Theta
)";


  // -hitResX        [float]      preset detector hit resolution X
  // -hitResY        [float]      preset detector hit resolution Y
  // -seedResX       [float]      preset seed track resolution X
  // -seedResY       [float]      preset seed track resolution Y


int main(int argc, char *argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"eventMax", required_argument, NULL, 'm'},
                              {"hitFile", required_argument, NULL, 'f'},
                              {"trackFile", required_argument, NULL, 'b'},
                              {"geomertyFile", required_argument, NULL, 'g'},
                              {"energy", required_argument, NULL, 'e'},
                              {"dutID", required_argument, NULL, 'd'},
                              // {"hitResX", required_argument, NULL, 'r'},
                              // {"hitResY", required_argument, NULL, 's'},
                              // {"seedResX", required_argument, NULL, 'j'},
                              // {"seedResY", required_argument, NULL, 'k'},
                              {"seedResPhi", required_argument, NULL, 't'},
                              {"seedResTheta", required_argument, NULL, 'w'},
                              {0, 0, 0, 0}};

  int64_t eventMaxNum = -1;
  int64_t dutID = -1;
  std::string datafile_name;
  std::string fitted_datafile_name;
  std::string geofile_name;
  double beamEnergy = -1;
  // double resX = 5_um;
  // double resY = 5_um;
  // Use large starting parameter covariance

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
      fitted_datafile_name = optarg;
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
    // case 'r':
    //   resX = std::stod(optarg);
    //   break;
    // case 's':
    //   resY = std::stod(optarg);
    //   break;
    case 'j':
      seedResX = std::stod(optarg);
      break;
    case 'k':
      seedResY = std::stod(optarg);
      break;
    case 't':
      seedResPhi = std::stod(optarg);
      break;
    case 'w':
      seedResTheta = std::stod(optarg);
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
  std::fprintf(stdout, "fittedFile:       %s\n", fitted_datafile_name.c_str());
  // std::fprintf(stdout, "resX:             %f\n", resX);
  // std::fprintf(stdout, "resY:             %f\n", resY);
  // std::fprintf(stdout, "seedResX:         %f\n", seedResX);
  // std::fprintf(stdout, "seedResY:         %f\n", seedResY);
  std::fprintf(stdout, "seedResPhi:       %f\n", seedResPhi);
  std::fprintf(stdout, "seedResTheta:     %f\n", seedResTheta);
  std::fprintf(stdout, "\n");

  if (datafile_name.empty() || fitted_datafile_name.empty() ||
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
  std::tie(trackingGeometry, element_col)  = TelActs::buildGeometry(gctx, jsd_geo);
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
  // 500 chi2cut,   1 max number of selected sourcelinks in a single surface;

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{0, 0., 1.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;
  pOptions.mass = 0.511 * Acts::UnitConstants::MeV;

  // Set the CombinatorialKalmanFilter options
  // @Todo: add options for CKF  Acts::LoggerWrapper{kfLogger}
  auto kfLogger = Acts::getDefaultLogger("KalmanFilter", logLevel);
  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector> ckfOptions(
    gctx, mctx, cctx, sourcelinkSelectorCfg, Acts::LoggerWrapper{*kfLogger}, pOptions,
    refSurface.get());

///////////////////////////////////////////////////
  JsonFileDeserializer jsfd(datafile_name);
  JsonFileSerializer jsfs(fitted_datafile_name);

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

    //////////////////////////////////////////////////////////////
    // JsonValue js_output(rapidjson::kObjectType);
    // js_output.AddMember("eventNum", JsonValue(eventNumber), jsa);
    // JsonValue js_tracks(rapidjson::kArrayType);
    // for (const auto &traj : trajectories) {
    //   auto js_track_mj = traj.createJsonValue(jsa, gctx);
    //   for(auto &js_track : js_track_mj.GetArray()){
    //     js_tracks.PushBack(std::move(js_track), jsa);
    //   }
    // }
    // js_output.AddMember("hits", std::move(js_hits), jsa);
    // js_output.AddMember("tracks", std::move(js_tracks), jsa);
    // jsfs.putNextJsonValue(js_output);

    eventNumber ++;
    //TODO gl
  }
  return 0;
}
