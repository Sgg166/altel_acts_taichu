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

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Framework/IWriter.hpp"

#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrackReader.hpp"
#include "TelescopeTrack.hpp"
#include "TelescopeTrackingPerformanceWriter.hpp"
#include "TelescopeTrackFindingAlgorithm.hpp"
#include "ObjTelescopeTrackWriter.hpp"
#include "RootTelescopeTrackWriter.hpp"
#include "TelescopeTrackReader.hpp"

#include <memory>
#include "getopt.h"
#include "myrapidjson.h"

using JsonValue = rapidjson::GenericValue<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;
using JsonDocument = rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>;

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -file       [jsonfile]   name of data json file
  -energy     [float]      beam energy
  -geo        [jsonfile]   geometry input file
  -resX       [float]      preset detector hit resolution X
  -resY       [float]      preset detector hit resolution Y
  -resPhi     [float]      preset seed track resolution Phi
  -resTheta   [float]      preset seed track resolution Theta

)";

int main(int argc, char* argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "file",      required_argument, NULL,           'f' },
     { "outputdir",      required_argument, NULL,      'o' },
     { "energy",       required_argument, NULL,        'e' },
     { "geomerty",       required_argument, NULL,      'g' },
     { "resX",       required_argument, NULL,          'r' },
     { "resY",       required_argument, NULL,          's' },
     { "resPhi",       required_argument, NULL,        't' },
     { "resTheta",       required_argument, NULL,      'w' },
     { 0, 0, 0, 0 }};


  std::string datafile_name;
  std::string geofile_name;
  std::string outputDir = "./";
  double beamEnergy = 4 * Acts::UnitConstants::GeV;
  double resX = 5_um;
  double resY = 5_um;
  // Use large starting parameter covariance
  double resLoc1 = 15_mm;
  double resLoc2 = 15_mm;
  double resPhi = 0.7;
  double resTheta = 0.7;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(0);
      break;
    case 'f':
      datafile_name = optarg;
      break;
    case 'e':
      beamEnergy = std::stod(optarg) * Acts::UnitConstants::GeV;
      break;
    case 'g':
      geofile_name = optarg;
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
    case 't':
      resPhi = std::stod(optarg);
      break;
    case 'w':
      resTheta = std::stod(optarg);
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

  if(datafile_name.empty() || outputDir.empty()){
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "datafile:         %s\n", datafile_name.c_str());
  std::fprintf(stdout, "geofile:          %s\n", geofile_name.c_str());
  std::fprintf(stdout, "outputDir:        %s\n", outputDir.c_str());
  std::fprintf(stdout, "beamEnergy:       %f\n", beamEnergy);
  std::fprintf(stdout, "resX:             %f\n", resX);
  std::fprintf(stdout, "resY:             %f\n", resY);
  std::fprintf(stdout, "resPhi:           %f\n", resPhi);
  std::fprintf(stdout, "resTheta:         %f\n", resTheta);
  std::fprintf(stdout, "\n");


  ///////select datapacks/////////////////////
  JsonValue js_selected_datapack_col(rapidjson::kArrayType);
  size_t n_datapack_select_opt = 20000;
  {
    std::FILE* fp = std::fopen(datafile_name.c_str(), "r");
    if(!fp) {
      std::fprintf(stderr, "File opening failed: %s \n", datafile_name.c_str());
      throw std::system_error(EIO, std::generic_category(), "File opening failed");;
    }

    char readBuffer[UINT16_MAX];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc(&jsa);
    doc.ParseStream(is);
    std::fclose(fp);

    if(!doc.IsArray() || !doc.GetArray().Size()){
      std::fprintf(stderr, "no, it is not data array\n");
      throw std::system_error(EDOM, std::generic_category(), "File is not valid json array");;
    }

    uint64_t processed_datapack_count = 0;
    for(auto ev_it = doc.Begin(); ev_it != doc.End() && js_selected_datapack_col.Size() < n_datapack_select_opt ; ev_it++ ){
      auto &evpack = *ev_it;
      processed_datapack_count ++;
      auto &frames = evpack["layers"];
      bool is_good_datapack = true;
      for(auto& layer : evpack["layers"].GetArray()){
        uint64_t l_hit_n = layer["hit"].GetArray().Size();
        if(l_hit_n != 1){
          is_good_datapack = false;
          continue;
        }
      }
      if(!is_good_datapack){
        continue;
      }
      js_selected_datapack_col.PushBack(evpack, jsa);
    }

    std::fprintf(stdout,
                 "Select %lu datapacks out of %lu processed datapackes. The data file contains %lu datapacks.\n",
                 js_selected_datapack_col.Size(), processed_datapack_count, doc.GetArray().Size() );

  }
  ///////end of select datapacks/////////////////////

  JsonValue js_geometry(rapidjson::kArrayType);
  if(!geofile_name.empty())
    {
      std::FILE* fp = std::fopen(geofile_name.c_str(), "r");
      if(!fp) {
	std::fprintf(stderr, "File opening failed: %s \n", geofile_name.c_str());
	throw std::system_error(EIO, std::generic_category(), "File opening failed");;
      }
      char readBuffer[UINT16_MAX];
      rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
      rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc(&jsa);
      doc.ParseStream(is);
      std::fclose(fp);
      js_geometry.CopyFrom<rapidjson::CrtAllocator>(doc["alignment_result"],jsa);
    }
  else{
    std::vector<std::vector<double>> positions{{0., 0., -95_mm},
					       {0., 0., -57_mm},
					       {0., 0., -19_mm},
					       {0., 0., 19_mm},
					       {0., 0., 57_mm},
					       {0., 0., 95_mm}};
	  
    for(auto &p: positions){
      JsonValue js_ele(rapidjson::kObjectType);
      js_ele.AddMember("centerX",    JsonValue(p[0]), jsa);
      js_ele.AddMember("centerY",    JsonValue(p[1]), jsa);
      js_ele.AddMember("centerZ",    JsonValue(p[2]), jsa);
      js_ele.AddMember("rotX", JsonValue(0), jsa);
      js_ele.AddMember("rotY", JsonValue(0), jsa);
      js_ele.AddMember("rotZ", JsonValue(0), jsa);
      js_geometry.PushBack(std::move(js_ele), jsa);
    }
  }

  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  // Setup detector geometry
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> element_col;
  std::vector<std::shared_ptr<const Acts::Surface>> surface_col;
  std::vector<Acts::LayerPtr> layer_col;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

  Telescope::BuildGeometry(gctx, trackingGeometry, element_col, surface_col, layer_col, js_geometry);
  // Set up the detector elements to be aligned (fix the first one)
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> dets;
  unsigned int idet = 0;
  for (const auto& det : element_col) {
    idet++;
    // Skip the first detector element
    if (idet == 1) {
      continue;
    }
    dets.push_back(det);
  }
  //

  /////////////////////////////////////
  Acts::Logging::Level logLevel = do_verbose? (Acts::Logging::VERBOSE):(Acts::Logging::INFO); 
  auto trackFinderFun = Telescope::TelescopeTrackFindingAlgorithm::makeTrackFinderFunction
    (trackingGeometry, magneticField,  logLevel);

  // write tracks as root tree
  // @Todo: adapt the writer to write trajectories produced by CKF
  Telescope::RootTelescopeTrackWriter::Config conf_trackRootWriter;
  conf_trackRootWriter.inputTrajectories = "trajectories";
  conf_trackRootWriter.outputDir = outputDir;
  conf_trackRootWriter.outputFilename = "telescope_tracks.root";
  conf_trackRootWriter.outputTreename = "tracks";

  // write the tracks (measurements only for the moment) as Csv
  Telescope::ObjTelescopeTrackWriter::Config conf_trackObjWriter;
  conf_trackObjWriter.inputTrajectories = "trajectories";
  conf_trackObjWriter.outputDir = outputDir;
  // The number of tracks you want to show (in default, all of tracks will be
  // shown)
  conf_trackObjWriter.maxNumTracks = 100;
  
  // write reconstruction performance data
  // @Todo: adapt the writer to write trajectories produced by CKF
  Telescope::TelescopeTrackingPerformanceWriter::Config conf_perfFitter;
  conf_perfFitter.inputTrajectories = "trajectories";
  conf_perfFitter.outputDir = outputDir;

  ////////////////seq////////////////////////
  FW::Sequencer::Config conf_seq;
  conf_seq.logLevel  = Acts::Logging::INFO;
  conf_seq.outputDir =  outputDir;
  conf_seq.numThreads = 1;
  conf_seq.events = 1;


  std::vector<std::shared_ptr<FW::IWriter>> writer_col;
  writer_col.push_back(std::make_shared<Telescope::RootTelescopeTrackWriter>(conf_trackRootWriter,  Acts::Logging::INFO));
  writer_col.push_back(std::make_shared<Telescope::ObjTelescopeTrackWriter>(conf_trackObjWriter,  Acts::Logging::INFO));
  writer_col.push_back(std::make_shared<Telescope::TelescopeTrackingPerformanceWriter>(conf_perfFitter,  Acts::Logging::INFO));

  //run
  {
    FW::WhiteBoard eventStore;
    FW::AlgorithmContext context(0, 0, eventStore);

    // Setup local covariance
    Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
    cov_hit(0, 0) =resX*resX;
    cov_hit(1, 1) =resY*resY;

    // One entry in the container below is the sorted source links (the layer as the key) in one event
    std::vector<std::vector<Telescope::PixelSourceLink>> sourcelinkTracks;
    sourcelinkTracks.reserve(js_selected_datapack_col.Size());
    for(auto ev_it = js_selected_datapack_col.Begin(); ev_it != js_selected_datapack_col.End() ; ev_it++ ){
      const auto &evpack = *ev_it;
      std::vector<Telescope::PixelSourceLink> sourcelinks;
      auto &frames = evpack["layers"];
      for(size_t i= 0; i< 6; i++){
        double x = frames[i]["hit"][0]["pos"][0].GetDouble() - 0.02924*1024/2.0;
        double y = frames[i]["hit"][0]["pos"][1].GetDouble() - 0.02688*512/2.0;
        Acts::Vector2D loc;
        loc << x, y;
        sourcelinks.emplace_back(*surface_col.at(i), loc, cov_hit);
      }
      sourcelinkTracks.push_back(sourcelinks);
    }

    std::vector<Telescope::PixelMultiTrajectory> trajectories;
    trajectories.reserve(js_selected_datapack_col.Size());
    size_t itrack = 0;
    for(const auto& trackSourcelinks : sourcelinkTracks) {
      itrack++;
      if (trackSourcelinks.empty()) {
        trajectories.push_back(Telescope::PixelMultiTrajectory());
        std::cout<<"Empty track " << itrack << " found."<<std::endl;
        continue;
      }
      std::cout<<"There are " << trackSourcelinks.size() << " source links in this event "<<std::endl;
     
      // Start to find seeds for this event using the source links on the first two layers
      // @todo: add raw layer id in PixelSourceLink
      std::vector<Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy>> initialParameters; 
      for(const auto& sl0 : trackSourcelinks){
	const auto& surface0 = sl0.referenceSurface();
	const auto& layer0 = surface0.geoID().layer();
	if(layer0 != 2) {
	  continue;
	}
	const Acts::Vector3D global0 = sl0.globalPosition(gctx); 
        for(const auto& sl1 : trackSourcelinks){
	  const auto& surface1 = sl1.referenceSurface();
	  const auto& layer1 = surface1.geoID().layer();
	  if(layer1 != 4) {
	    continue;
	  }
          const Acts::Vector3D global1 = sl1.globalPosition(gctx);
	  Acts::Vector3D distVec = global1 - global0; 
          // compare their distance in x-y (r-phi) plane
          const double rDist = std::abs(Acts::VectorHelpers::perp(distVec));
	  // @todo: add options for the seed cuts
	  if(rDist<=3_mm) {
	    Acts::BoundSymMatrix cov;
	    cov <<
	      resLoc1 * resLoc1, 0., 0., 0., 0., 0.,
	      0., resLoc2 * resLoc2, 0., 0., 0., 0.,
	      0., 0., resPhi*resPhi, 0., 0., 0.,
	      0., 0., 0., resTheta*resTheta, 0., 0.,
	      0., 0., 0., 0., 0.0001, 0.,
	      0., 0., 0., 0.,     0., 1.;

	    const double phi = Acts::VectorHelpers::phi(distVec);
	    const double theta = Acts::VectorHelpers::theta(distVec);
	    Acts::Vector3D rPos = global0 - distVec / 2;
	    Acts::Vector3D rMom(beamEnergy * sin(theta) * cos(phi),
				beamEnergy * sin(theta) * sin(phi),
				beamEnergy * cos(theta));

            initialParameters.emplace_back(cov, rPos, rMom, 1., 0); 
	  }
	}
      }
      std::cout<<"There are " << initialParameters.size() << " created seeds"<<std::endl; 

      auto refSurface = Acts::Surface::makeShared<Acts::PlaneSurface>
        (Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{0, 0., 1.});
       
      // Set the CombinatorialKalmanFilter options
      // @Todo: add options for CKF
      Telescope::TelescopeTrackFindingAlgorithm::CKFOptions ckfOptions(
								       gctx, mctx, cctx,
								       Acts::CKFSourceLinkSelector::Config{{Acts::GeometryID(),{30, 5}}}, refSurface.get());
 
      //Loop ever the seeds
      size_t iseed = 0;
      for(const auto& rStart: initialParameters){ 
	auto result = trackFinderFun(trackSourcelinks, rStart, ckfOptions);
	if (result.ok()) {
	  // Get the track finding output object
	  const auto& trackFindingOutput = result.value();
	  // Create a PixelMultiTrajectory
	  trajectories.emplace_back(std::move(trackFindingOutput.fittedStates),
				    std::move(trackFindingOutput.trackTips),
				    std::move(trackFindingOutput.fittedParameters));
	} else {
	  std::cout<<"Track finding failed for track " << itrack <<" seed "
		   << iseed << " with error" << result.error()<<std::endl;
	  // Track finding failed, but still create an empty PixelMultiTrajectory
	  trajectories.push_back(Telescope::PixelMultiTrajectory());
	}
	iseed++;
      }// end of the loop for all seeds
    } // end of the loop for all read-in events

    eventStore.add("trajectories", std::move(trajectories));

    for (auto& wrt : writer_col) {
      if (wrt->write(++context) != FW::ProcessCode::SUCCESS) {
        throw std::runtime_error("Failed to write output data");
      }
    }

    for (auto& wrt : writer_col) {
      if (wrt->endRun() != FW::ProcessCode::SUCCESS) {
        throw std::runtime_error("Failed to write output data");
      }
    }

  }
  return 0;
}
