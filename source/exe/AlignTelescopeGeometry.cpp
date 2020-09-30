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

#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrack.hpp"
#include "TelescopeJsonTrackReader.hpp"
#include "JsonGenerator.hpp"

#include <memory>
#include "getopt.h"

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -file       [jsonfile]   name of data json file
  -energy     [float]      beam energy, GeV
  -out        [jsonfile]   alignment result
  -geo        [jsonfile]   geometry input file
  -resX       [float]      preset detector hit resolution X
  -resY       [float]      preset detector hit resolution Y
  -resPhi     [float]      preset seed track resolution Phi
  -resTheta   [float]      preset seed track resolution Theta
  -maxItera   [integer]    max time of alignement iterations
  -deltaItera [integer]    converge condition: delta iteration
  -deltaChi2  [float]      converge condition: delta Chi2ONdf

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
     { "energy",       required_argument, NULL,        'e' },
     { "out",       required_argument, NULL,           'o' },
     { "geomerty",       required_argument, NULL,      'g' },
     { "resX",       required_argument, NULL,          'r' },
     { "resY",       required_argument, NULL,          's' },
     { "resPhi",       required_argument, NULL,        't' },
     { "resTheta",       required_argument, NULL,      'w' },
     { "maxItera",       required_argument, NULL,      'x' },
     { "deltaItera",       required_argument, NULL,    'y' },
     { "deltaChi2",       required_argument, NULL,     'z' },
     { 0, 0, 0, 0 }};


  std::string datafile_name;
  std::string outputfile_name;
  std::string geofile_name;
  double resX = 150_um;
  double resY = 150_um;
  // Use large starting parameter covariance
  double resLoc1 = 50_um;
  double resLoc2 = 50_um;
  double resPhi = 0.7;
  double resTheta = 0.7;
  // Iterations converge criteria
  size_t maxNumIterations = 400;
  size_t nIterations = 10;
  double deltaChi2ONdf = 1e-5;
  double beamEnergy = -1;

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
    case 'o':
      outputfile_name = optarg;
      break;
    case 'e':
      beamEnergy = std::stod(optarg) * Acts::UnitConstants::GeV;
      break;
    case 'g':
      geofile_name = optarg;
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
    case 'x':
      maxNumIterations = std::stoul(optarg);
      break;
    case 'y':
      nIterations = std::stoul(optarg);
      break;
    case 'z':
      deltaChi2ONdf = std::stod(optarg);
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

  if(datafile_name.empty() || outputfile_name.empty()){
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "datafile:         %s\n", datafile_name.c_str());
  std::fprintf(stdout, "outfile:          %s\n", outputfile_name.c_str());
  std::fprintf(stdout, "geofile:          %s\n", geofile_name.c_str());
  std::fprintf(stdout, "resX:             %f\n", resX);
  std::fprintf(stdout, "resY:             %f\n", resY);
  std::fprintf(stdout, "resPhi:           %f\n", resPhi);
  std::fprintf(stdout, "resTheta:         %f\n", resTheta);
  std::fprintf(stdout, "maxNumIterations: %lu\n", maxNumIterations);
  std::fprintf(stdout, "nIterations:      %lu\n", nIterations);
  std::fprintf(stdout, "deltaChi2ONdf:    %f\n", deltaChi2ONdf);
  std::fprintf(stdout, "\n");

  Acts::GeometryContext gctx;
  Acts::MagneticFieldContext mctx;
  Acts::CalibrationContext cctx;

  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);

  // Setup detector geometry
  std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>> element_col;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;


  std::printf("--------read geo and beam energy\n");
  std::map<size_t, std::array<double, 6>> geoconf;
  if(!geofile_name.empty())
    geoconf = Telescope::JsonGenerator::ReadGeoFromGeoFile(geofile_name);
  else{
    geoconf = Telescope::JsonGenerator::ReadGeoFromDataFile(datafile_name);
  }
  for(const auto& [id, lgeo] : geoconf){
  std::printf("layer: %lu   centerX: %f   centerY: %f   centerZ: %f  rotationX: %f   rotationY: %f   rotationZ: %f\n",
              id, lgeo[0], lgeo[1], lgeo[2], lgeo[3], lgeo[4], lgeo[5]);
  }
  if(beamEnergy<0){
    beamEnergy = Telescope::JsonGenerator::ReadBeamEnergyFromDataFile(datafile_name) * Acts::UnitConstants::GeV;
  }

  std::fprintf(stdout, "beamEnergy:       %f\n", beamEnergy);

  Telescope::BuildGeometry(gctx, trackingGeometry, element_col, geoconf, 50_mm, 25_mm, 80_um);
  std::map<size_t, std::shared_ptr<const Acts::Surface>> surfaces_selected;
  for(const auto& e: element_col){
    auto id = e->telDetectorID();
    surfaces_selected[id] =  e->surface().getSharedPtr();
  }

  // Set up the detector elements to be aligned (fix the first one)
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> element_col_align;
  for (auto it = element_col.begin() ; it != element_col.end(); ++it){
    if(it != element_col.begin())
      element_col_align.push_back(*it);
  }
  //

  /////////////////////////////////////
  // The criteria to determine if the iteration has converged. @Todo: to use
  // delta chi2 instead
  double chi2ONdfCutOff = 0.01;
  std::pair<size_t, double> deltaChi2ONdfCutOff = {nIterations,deltaChi2ONdf};
  // set up the alignment dnf for each iteration
  std::map<unsigned int, std::bitset<6>> iterationState;
  for (unsigned int iIter = 0; iIter < maxNumIterations; iIter++) {
    std::bitset<6> mask(std::string("111111"));
    if(iIter%2 ==0){
      // only align offset along x/y
      mask = std::bitset<6>(std::string("000011"));
    }else{
      // only align rotation around beam
      mask = std::bitset<6>(std::string("100000"));
    }
    iterationState.emplace(iIter, mask);
  }

  FW::AlignedTransformUpdater alignedTransformUpdaterFun =
    [](Acts::DetectorElementBase& detElement,
       const Acts::GeometryContext& gctx,
       const Acts::Transform3D& aTransform) {
      Telescope::TelescopeDetectorElement* telescopeDetElement =
        dynamic_cast<Telescope::TelescopeDetectorElement*>(&detElement);
      if (telescopeDetElement) {
        telescopeDetElement->addAlignedTransform
          (std::make_unique<Acts::Transform3D>(aTransform));
        return true;
      }
      return false;
    };

  auto alignFun = Telescope::makeAlignmentFunction
    (trackingGeometry, magneticField, do_verbose? (Acts::Logging::VERBOSE):(Acts::Logging::INFO));

  {
    // Setup local covariance
    Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
    cov_hit(0, 0) =resX*resX;
    cov_hit(1, 1) =resY*resY;

    //
    size_t n_datapack_select_opt = 20000;
    Telescope::JsonDocument doc(&jsa);
    Telescope::JsonGenerator gen(datafile_name);
    std::vector<std::vector<Telescope::PixelSourceLink>> sourcelinkTracks;
    while(1){
      if(sourcelinkTracks.size()> n_datapack_select_opt){
        break;
      }
      doc.Populate(gen);
      if(!gen.isvalid  ){
        break;
      }

      Telescope::JsonValue evpack;
      doc.Swap(evpack);
      std::vector<Telescope::PixelSourceLink> sourcelinks;
      Telescope::TelescopeJsonTrackReader::createSourcelinksFromJSON(evpack, surfaces_selected, cov_hit, sourcelinks);

      //drop multiple hits events
      if(sourcelinks.size() != surfaces_selected.size()){
        continue;
      }
      std::vector<size_t> geo_ids;
      bool found_same_geo_id = false ;
      for(const auto &sl: sourcelinks){
        size_t geo_id = sl.referenceSurface().geoID().value();
        if(std::find( geo_ids.begin(), geo_ids.end(), geo_id ) != geo_ids.end()){
          found_same_geo_id = true;
          break;
        }
        geo_ids.push_back(geo_id);
      }
      if(found_same_geo_id){
        continue;
      }
      sourcelinkTracks.push_back(sourcelinks);
    }

    //seeding
    std::vector<Acts::CurvilinearParameters> initialParameters;
    initialParameters.reserve(sourcelinkTracks.size());
    Acts::BoundSymMatrix cov_seed;
    cov_seed <<
      resLoc1 * resLoc1,  0.,  0.,   0.,   0.,  0.,
      0., resLoc2 * resLoc2, 0.,   0.,   0.,  0.,
      0.,                 0., resPhi*resPhi, 0.,   0.,  0.,
      0.,                 0., 0.,   resTheta*resTheta, 0.,  0.,
      0.,                 0., 0.,   0.,   0.0001,  0.,
      0.,                 0., 0.,   0.,   0.,  1.;

    for( const auto& sourcelinks : sourcelinkTracks) {
      // Create initial parameters
      // The position is taken from the first measurement
      const Acts::Vector3D global0 = sourcelinks.at(0).globalPosition(gctx);
      const Acts::Vector3D global1 = sourcelinks.at(1).globalPosition(gctx);
      Acts::Vector3D distance = global1 - global0;

      const double phi = Acts::VectorHelpers::phi(distance);
      const double theta = Acts::VectorHelpers::theta(distance);
      Acts::Vector3D rPos = global0 - distance / 2;
      Acts::Vector3D rMom(beamEnergy * sin(theta) * cos(phi),
                          beamEnergy * sin(theta) * sin(phi),
                          beamEnergy * cos(theta));

      Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(cov_seed, rPos, rMom, 1., 0);
      initialParameters.push_back(rStart);
    }

    // Set the KalmanFitter options
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions
      (gctx, mctx, cctx, Acts::VoidOutlierFinder()); //pSurface default nullptr

    // Set the alignment options
    FW::AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>> alignOptions
      (kfOptions,
       alignedTransformUpdaterFun,
       element_col_align,
       chi2ONdfCutOff,
       deltaChi2ONdfCutOff,
       maxNumIterations,
       iterationState
       );

    std::printf("Invoke alignment");
    auto result = alignFun(sourcelinkTracks, initialParameters, alignOptions);

    if (!result.ok()){
      std::printf("Alignment failed with %s \n", result.error().message().c_str());
    }

    std::printf("Alignment finished with deltaChi2 = %f \n", result.value().deltaChi2);
    std::map<size_t, std::array<double, 6>> geo_result;
    for(const auto& det : element_col) {
      const auto& surface = &det->surface();
      const auto& transform = det->transform(gctx);
      const auto& translation = transform.translation();
      const auto& rotation = transform.rotation();
      const Acts::Vector3D rotAngles = rotation.eulerAngles(2, 1, 0);

      size_t id = det->telDetectorID();
      double cx = translation.x();
      double cy = translation.y();
      double cz = translation.z();
      double rx = rotAngles(2);
      double ry = rotAngles(1);
      double rz = rotAngles(0);
      geo_result[id]={cx, cy, cz, rx, ry, rz};

      std::printf("layer: %lu   centerX: %f   centerY: %f   centerZ: %f  rotationX: %f   rotationY: %f   rotationZ: %f\n",
                  id, cx, cy, cz, rx, ry, rz);
    }
    std::printf("select %lu datapacks (single track events) \n", sourcelinkTracks.size());
    Telescope::JsonGenerator::WriteGeoToGeoFile(outputfile_name, geo_result);
  }
  return 0;
}
