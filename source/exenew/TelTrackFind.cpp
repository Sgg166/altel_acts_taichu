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

#include "PixelMultiTrajectory.hpp"

#include "TelActs.hh"
#include "getopt.h"
#include "myrapidjson.h"

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                        help message
  -verbose                     verbose flag
  -thread         [int]        multiple threads
  -eventMax       [int]        max number of events
  -data           [jsonfile]   data input file
  -fittedFile     [jsonfile]   path to fitted data file (output)
  -geo            [jsonfile]   geometry input file
  -energy         [float]      beam energy, GeV
  -outputdir      [path]       output dir path
  -resX           [float]      preset detector hit resolution X
  -resY           [float]      preset detector hit resolution Y
  -seedResX       [float]      preset seed track resolution X
  -seedResY       [float]      preset seed track resolution Y
  -seedResPhi     [float]      preset seed track resolution Phi
  -seedResTheta   [float]      preset seed track resolution Theta
)";

int main(int argc, char *argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"thread", required_argument, NULL, 'd'},
                              {"eventMax", required_argument, NULL, 'm'},
                              {"data", required_argument, NULL, 'f'},
                              {"fittedFile", required_argument, NULL, 'b'},
                              {"outputdir", required_argument, NULL, 'o'},
                              {"geomerty", required_argument, NULL, 'g'},
                              {"energy", required_argument, NULL, 'e'},
                              {"resX", required_argument, NULL, 'r'},
                              {"resY", required_argument, NULL, 's'},
                              {"seedResX", required_argument, NULL, 'j'},
                              {"seedResY", required_argument, NULL, 'k'},
                              {"seedResPhi", required_argument, NULL, 't'},
                              {"seedResTheta", required_argument, NULL, 'w'},
                              {0, 0, 0, 0}};

  size_t threadNum = 1;
  size_t eventMaxNum = -1;
  std::string datafile_name;
  std::string fitted_datafile_name;
  std::string geofile_name;
  std::string outputDir;
  double beamEnergy = -1;
  double resX = 5_um;
  double resY = 5_um;
  // Use large starting parameter covariance

  double seedResX = 15_mm;
  double seedResY = 15_mm;
  double seedResPhi = 0.7_rad;
  double seedResTheta = 0.7_rad;


  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(0);
      break;
    case 'd':
      threadNum = std::stoul(optarg);
      break;
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
    case 'o':
      outputDir = optarg;
      break;
    case 'r':
      resX = std::stod(optarg);
      break;
    case 's':
      resY = std::stod(optarg);
      break;
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

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "datafile:         %s\n", datafile_name.c_str());
  std::fprintf(stdout, "geofile:          %s\n", geofile_name.c_str());
  std::fprintf(stdout, "outputDir:        %s\n", outputDir.c_str());
  std::fprintf(stdout, "fittedFile:       %s\n", fitted_datafile_name.c_str());
  std::fprintf(stdout, "resX:             %f\n", resX);
  std::fprintf(stdout, "resY:             %f\n", resY);
  std::fprintf(stdout, "seedResX:         %f\n", seedResX);
  std::fprintf(stdout, "seedResY:         %f\n", seedResY);
  std::fprintf(stdout, "seedResPhi:       %f\n", seedResPhi);
  std::fprintf(stdout, "seedResTheta:     %f\n", seedResTheta);
  std::fprintf(stdout, "\n");

  if (datafile_name.empty() || fitted_datafile_name.empty() ||
      outputDir.empty() || geofile_name.empty()) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::filesystem::path path_dir_output = std::filesystem::absolute(outputDir);
  std::filesystem::file_status st_dir_output =
      std::filesystem::status(path_dir_output);
  if (!std::filesystem::exists(st_dir_output)) {
    std::fprintf(stderr, "Output folder does not exist: %s\n",
                 path_dir_output.c_str());
    std::filesystem::file_status st_parent =
        std::filesystem::status(path_dir_output.parent_path());
    if (std::filesystem::exists(st_parent) &&
        std::filesystem::is_directory(st_parent)) {
      if (std::filesystem::create_directory(path_dir_output)) {
        std::printf("Create output folder: %s\n", path_dir_output.c_str());
      } else {
        throw;
      }
    } else {
      throw;
    }
  }


  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(geofile_name);
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geofile_name.c_str() );
    throw;
  }
  JsonUtils::printJsonValue(jsd_geo, true);

  if (beamEnergy < 0) {
    beamEnergy = 5.0 * Acts::UnitConstants::GeV;
  }
  std::fprintf(stdout, "beamEnergy:       %f\n", beamEnergy);

  std::printf("--------create acts geo object-----\n");
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  std::vector<std::shared_ptr<TelActs::TelescopeDetectorElement>> element_col;
  std::tie(trackingGeometry, element_col)  = TelActs::buildGeometry(gctx, jsd_geo);
  //40_mm, 20_mm, 80_um

  // Set up surfaces
  std::map<size_t, std::shared_ptr<const Acts::Surface>> surfaces_selected;
  for (const auto &e : element_col) {
    auto id = e->telDetectorID();
    surfaces_selected[id] = e->surface().getSharedPtr();
  }

  /////////////////////////////////////
  Acts::Logging::Level logLevel =
      do_verbose ? (Acts::Logging::VERBOSE) : (Acts::Logging::INFO);

  /////////////  track seed conf
  size_t seedSurfaceGeoIDStart = surfaces_selected[0]->geometryId().value();
  size_t seedSurfaceGeoIDEnd = surfaces_selected[1]->geometryId().value();

  Acts::BoundSymMatrix cov;
  cov << seedResX * seedResX, 0., 0., 0., 0., 0., 0.,
    seedResY * seedResY, 0., 0., 0., 0., 0., 0.,
    seedResPhi * seedResPhi, 0., 0., 0., 0., 0., 0.,
    seedResTheta * seedResTheta, 0., 0., 0., 0., 0., 0.,
    0.0001, 0., 0., 0., 0., 0., 0., 1.;


  ///////////// trackfind conf
  auto trackFindFun = TelActs::makeTrackFinderFunction(trackingGeometry, magneticField);

  Acts::CKFSourceLinkSelector::Config sourcelinkSelectorCfg = Acts::CKFSourceLinkSelector::Config{
    {Acts::GeometryIdentifier(), {500, 1}}}; //max chi2, max number of selected hits

  auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
    Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{0, 0., 1.});

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  // Set the CombinatorialKalmanFilter options
  // @Todo: add options for CKF  Acts::LoggerWrapper{kfLogger}
  auto kfLogger = Acts::getDefaultLogger("KalmanFilter", logLevel);
  Acts::CombinatorialKalmanFilterOptions<Acts::CKFSourceLinkSelector> ckfOptions(
    gctx, mctx, cctx, sourcelinkSelectorCfg, Acts::LoggerWrapper{*kfLogger}, pOptions,
    refSurface.get());
  // 500 chi2cut,   1 max number of selected sourcelinks in a single surface;


///////////////////////////////////////////////////
  JsonFileDeserializer jsf(datafile_name);
  size_t eventNumber = 0;
  while(jsf && (eventNumber< eventMaxNum || eventMaxNum<0)){
    auto evpack = jsf.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }

    std::vector<TelActs::PixelSourceLink> sourcelinks;
    // TODO, read source link


    if (sourcelinks.empty()) {
      std::fprintf(stdout, "Empty event <%d> found.\n", eventNumber);
      eventNumber ++;
      continue;
    }

    // Start to find seeds for this event using the source links on the first
    // two layers
    // @todo: add raw layer id in PixelSourceLink
    std::vector<Acts::CurvilinearTrackParameters> initialParameters;
    for (const auto &sl0 : sourcelinks) {
      const auto &surface0 = sl0.referenceSurface();
      if (surface0.geometryId().value() != seedSurfaceGeoIDStart) {
        continue;
      }
      const Acts::Vector3D global0 = sl0.globalPosition(gctx);
      for (const auto &sl1 : sourcelinks) {
        const auto &surface1 = sl1.referenceSurface();
        if (surface1.geometryId().value() != seedSurfaceGeoIDEnd) {
          continue;
        }
        const Acts::Vector3D global1 = sl1.globalPosition(gctx);
        Acts::Vector3D distVec = global1 - global0;
        // compare their distance in x-y (r-phi) plane
        const double rDist = std::abs(Acts::VectorHelpers::perp(distVec));
        // @todo: add options for the seed cuts
        if (rDist <= 3_mm) {
          const double phi = Acts::VectorHelpers::phi(distVec);
          const double theta = Acts::VectorHelpers::theta(distVec);
          Acts::Vector3D rPos = global0 - distVec / 2;
          Acts::Vector4D rPos4(rPos.x(), rPos.y(), rPos.z(), 0);
          double q = 1;
          initialParameters.emplace_back(rPos4, phi, theta, beamEnergy, q,cov);
        }
      }
    }

    std::vector<TelActs::PixelMultiTrajectory> trajectories;
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
    eventNumber ++;
    //TODO: write  fitted_datafile_name
  }

  return 0;
}
