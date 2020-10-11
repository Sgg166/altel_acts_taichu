// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/ConstantBField.hpp"
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

#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Framework/IWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"

#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrack.hpp"
#include "TelescopeTrackingPerformanceWriter.hpp"
#include "TelescopeTrackFindingAlgorithm.hpp"
#include "ObjTelescopeTrackWriter.hpp"
#include "RootTelescopeTrackWriter.hpp"
#include "TelescopeJsonTrackReader.hpp"
#include "TelescopeJsonTrackWriter.hpp"
#include "JsonGenerator.hpp"


#include <filesystem>
#include <memory>
#include "getopt.h"

using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;
using JsonDocument = rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

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

int main(int argc, char* argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",           no_argument,       &do_help,      1  },
     { "verbose",        no_argument,       &do_verbose,   1  },
     { "thread",         required_argument, NULL,           'd' },
     { "eventMax",       required_argument, NULL,           'm' },
     { "data",           required_argument, NULL,           'f' },
     { "fittedFile",     required_argument, NULL,           'b' },
     { "outputdir",      required_argument, NULL,      'o' },
     { "geomerty",       required_argument, NULL,      'g' },
     { "energy",       required_argument, NULL,        'e' },
     { "resX",           required_argument, NULL,          'r' },
     { "resY",           required_argument, NULL,          's' },
     { "seedResX",       required_argument, NULL,        'j' },
     { "seedResY",       required_argument, NULL,      'k' },
     { "seedResPhi",     required_argument, NULL,        't' },
     { "seedResTheta",   required_argument, NULL,      'w' },
     { 0, 0, 0, 0 }};


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
  double seedResPhi  = 0.7_rad;
  double seedResTheta = 0.7_rad;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
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
      fprintf(stderr,"case 1\n");
      exit(1);
      break;
    case ':':
      fprintf(stderr,"case :\n");
      exit(1);
      break;
    case '?':
      fprintf(stderr,"case ?\n");
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


  if(datafile_name.empty() || fitted_datafile_name.empty() || outputDir.empty() || geofile_name.empty() ){
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::filesystem::path path_dir_output = std::filesystem::absolute(outputDir);
  std::filesystem::file_status st_dir_output = std::filesystem::status(path_dir_output);
  if(!std::filesystem::exists(st_dir_output)){
    std::fprintf(stderr, "Output folder does not exist: %s\n", path_dir_output.c_str());
    std::filesystem::file_status st_parent = std::filesystem::status(path_dir_output.parent_path());
    if(std::filesystem::exists(st_parent) && std::filesystem::is_directory(st_parent)){
      if(std::filesystem::create_directory(path_dir_output)){
        std::printf("Create output folder: %s\n", path_dir_output.c_str());
      }
      else{
        throw;
      }
    }
    else{
      throw;
    }
  }

  std::map<size_t, std::array<double, 6>> geoconf;
  geoconf = Telescope::JsonGenerator::ReadGeoFromGeoFile(geofile_name);
  for(const auto& [id, lgeo] : geoconf){
    std::printf("layer: %lu   centerX: %f   centerY: %f   centerZ: %f  rotationX: %f   rotationY: %f   rotationZ: %f\n",
                id, lgeo[0], lgeo[1], lgeo[2], lgeo[3], lgeo[4], lgeo[5]);
  }

  if(beamEnergy<0){
    beamEnergy = Telescope::JsonGenerator::ReadBeamEnergyFromDataFile(datafile_name) * Acts::UnitConstants::GeV;
  }

  std::fprintf(stdout, "beamEnergy:       %f\n", beamEnergy);

  Acts::GeometryContext gctx;
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  // Setup detector geometry
  std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>> element_col;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

  Telescope::BuildGeometry(gctx, trackingGeometry, element_col, geoconf, 40_mm, 20_mm, 80_um);
  std::map<size_t, std::shared_ptr<const Acts::Surface>> surfaces_selected;
  for(const auto& e: element_col){
    auto id = e->telDetectorID();
    surfaces_selected[id] =  e->surface().getSharedPtr();
  }

  /////////////////////////////////////
  Acts::Logging::Level logLevel = do_verbose? (Acts::Logging::VERBOSE):(Acts::Logging::INFO);

  // logLevel = Acts::Logging::WARNING;

  Telescope::TelescopeJsonTrackReader::Config conf_reader;
  conf_reader.inputDataFile = datafile_name;
  conf_reader.outputSourcelinks = "sourcelinks";
  conf_reader.resX = resX;
  conf_reader.resY = resY;
  conf_reader.surfaces = surfaces_selected;

  Telescope::TelescopeTrackFindingAlgorithm::Config conf_trackfinding;
  conf_trackfinding.inputSourcelinks="sourcelinks";
  conf_trackfinding.outputTrajectories="trajectories";
  conf_trackfinding.findTracks=Telescope::TelescopeTrackFindingAlgorithm::makeTrackFinderFunction
    (trackingGeometry, magneticField);
  conf_trackfinding.sourcelinkSelectorCfg = Acts::CKFSourceLinkSelector::Config{{Acts::GeometryIdentifier(),{500, 1}}};
  conf_trackfinding.seedSurfaceGeoIDStart = surfaces_selected[0]->geometryId().value();
  conf_trackfinding.seedSurfaceGeoIDEnd = surfaces_selected[1]->geometryId().value();
  conf_trackfinding.seedResX = seedResX;
  conf_trackfinding.seedResY = seedResY;
  conf_trackfinding.seedResPhi = seedResPhi;
  conf_trackfinding.seedResTheta = seedResTheta;
  conf_trackfinding.seedEnergy = beamEnergy;


  // write the tracks (measurements only for the moment) as Obj
  Telescope::ObjTelescopeTrackWriter::Config conf_trackObjWriter;
  conf_trackObjWriter.inputTrajectories = "trajectories";
  conf_trackObjWriter.outputDir = path_dir_output.c_str();
  // The number of tracks you want to show (in default, all of tracks will be
  // shown)
  conf_trackObjWriter.maxNumTracks = 100;

  ////////////////seq////////////////////////
  ActsExamples::Sequencer::Config conf_seq;
  conf_seq.logLevel  = logLevel;
  conf_seq.outputDir = path_dir_output.c_str();
  conf_seq.numThreads = threadNum;
  conf_seq.events = eventMaxNum;
  ActsExamples::Sequencer seq(conf_seq);

  seq.addReader(std::make_shared<Telescope::TelescopeJsonTrackReader>(conf_reader, logLevel));
  seq.addAlgorithm(std::make_shared<Telescope::TelescopeTrackFindingAlgorithm>(conf_trackfinding, logLevel) );


  // write tracks as root tree
  // @Todo: adapt the writer to write trajectories produced by CKF
  Telescope::RootTelescopeTrackWriter::Config conf_trackRootWriter;
  conf_trackRootWriter.inputTrajectories = "trajectories";
  conf_trackRootWriter.outputDir = path_dir_output.c_str();
  conf_trackRootWriter.outputFilename = "telescope_tracks.root";
  conf_trackRootWriter.outputTreename = "tracks";
  seq.addWriter(std::make_shared<Telescope::RootTelescopeTrackWriter>(conf_trackRootWriter,  logLevel));


  Telescope::TelescopeJsonTrackWriter::Config conf_trackJsonWriter;
  conf_trackJsonWriter.inputTrajectories = "trajectories";
  conf_trackJsonWriter.outputDir = path_dir_output.c_str();
  conf_trackJsonWriter.outputFileName = fitted_datafile_name;
  conf_trackJsonWriter.trackSurfaces = surfaces_selected;
  seq.addWriter(std::make_shared<Telescope::TelescopeJsonTrackWriter>(conf_trackJsonWriter,  logLevel));

  seq.run();

  return 0;
}
