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

#include "TelescopeAlignmentAlgorithm.hpp"
#include "TelescopeDetectorElement.hpp"
#include "TelescopeTrackReader.hpp"


#include <memory>
#include "getopt.h"
#include "myrapidjson.h"



using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -file [jsonfile]   name of data json file
)";

int main(int argc, char* argv[]) {  
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] =
    {
     { "help",       no_argument,       &do_help,      1  },
     { "verbose",    no_argument,       &do_verbose,   1  },
     { "file",     required_argument, NULL,           'f' },
     { 0, 0, 0, 0 }};
  
  std::string datafile_name;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      break;
    case 'f':
      datafile_name = optarg;
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

  if(do_help){
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    exit(0);
  }
  
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  std::string inputDir = "./";
  std::string outputDir = "./";
  
  // Setup the magnetic field
  auto magneticField = std::make_shared<Acts::ConstantBField>(0_T, 0_T, 0_T);  

  
  // Setup detector geometry
  /////////////////////////////////////
  Telescope::TelescopeDetectorElement::ContextType gctx;
  std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>> element_col;

  // Construct the rotation
  Acts::RotationMatrix3D rotation = Acts::RotationMatrix3D::Identity();
  // Euler angle around y will be -90_degree
  double rotationAngle = 90_degree;
  rotation.col(0) = Acts::Vector3D(std::cos(rotationAngle), 0., std::sin(rotationAngle));
  rotation.col(1) = Acts::Vector3D(0., 1., 0.);
  rotation.col(2) = Acts::Vector3D(-std::sin(rotationAngle), 0., std::cos(rotationAngle));
  
  // Boundaries of the surfaces (ALPIDE SIZE: 27.52512 * 13.76256_mm*mm)
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>
    (Acts::RectangleBounds(1000_mm, 1000_mm));
  // (Acts::RectangleBounds(27.52512_mm, 13.76256_mm));

  // Material of the surfaces
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>
    (Acts::MaterialProperties(95.7, 465.2, 28.03, 14., 2.32e-3, 50_um));
  
  // Set translation vectors
  std::vector<Acts::Vector3D> translations{
                                           {-95_mm, 0., 0.},
                                           {-57_mm, 0., 0.},
                                           {-19_mm, 0., 0.},
                                           {19_mm, 0., 0.},
                                           {57_mm, 0., 0.},
                                           {95_mm, 0., 0.}};

  std::vector<std::shared_ptr<const Acts::Surface>> surface_col;
  std::vector<Acts::LayerPtr> layer_col;
  
  for (unsigned int i = 0; i < translations.size(); i++) {
    auto trafo = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity() * rotation);
    trafo->translation() = translations[i];
    // Create the detector element
    auto detElement = std::make_shared<Telescope::TelescopeDetectorElement>(trafo, rBounds, 1_um, surfaceMaterial);
    auto detSurface = detElement->surface().getSharedPtr();
    surface_col.push_back(detSurface);
    element_col.push_back(detElement);
    
    std::unique_ptr<Acts::SurfaceArray> surArray(new Acts::SurfaceArray(detSurface));
    auto detLayer = Acts::PlaneLayer::create(trafo, rBounds, std::move(surArray), 1.5_cm);
    layer_col.push_back(detLayer);
    
    auto mutableSurface = std::const_pointer_cast<Acts::Surface>(detSurface);
    mutableSurface->associateLayer(*detLayer);
  }

  // Build tracking volume

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig);  
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      gctx, layer_col, -0.5_m - 1._mm, 0.5_m + 1._mm, Acts::BinningType::arbitrary,
      Acts::BinningValue::binX));
  
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);
  auto trafoVol = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  trafoVol->translation() = Acts::Vector3D(0., 0., 0.);
  auto trackVolume = Acts::TrackingVolume::create
    (trafoVol, boundsVol, nullptr, std::move(layArr), nullptr, {}, "Tracker");

  // Build world volume
  std::vector<std::pair<Acts::TrackingVolumePtr, Acts::Vector3D>> tapVec;
  tapVec.push_back(std::make_pair(trackVolume, Acts::Vector3D(0., 0., 0.)));

  Acts::BinningData binData(Acts::BinningOption::open, Acts::BinningValue::binX,
                            std::vector<float>{-0.6_m, 0.6_m});
  std::unique_ptr<const Acts::BinUtility> bu(new Acts::BinUtility(binData));

  std::shared_ptr<const Acts::TrackingVolumeArray> trVolArr(
      new Acts::BinnedArrayXD<Acts::TrackingVolumePtr>(tapVec, std::move(bu)));

  auto trafoWorld = std::make_shared<Acts::Transform3D>(Acts::Transform3D::Identity());
  trafoWorld->translation() = Acts::Vector3D(0., 0., 0.);
  auto worldVol = std::make_shared<const Acts::CuboidVolumeBounds>(1.2_m, 0.1_m, 0.1_m);
  Acts::MutableTrackingVolumePtr mtvpWorld(Acts::TrackingVolume::create
                                           (trafoWorld, worldVol, trVolArr, "World"));

  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry(new Acts::TrackingGeometry(mtvpWorld));

  // setup the alignment algorithm
  Telescope::TelescopeAlignmentAlgorithm::Config conf_alignment;
  //@Todo: add run number information in the file name
  // The number of tracks you want to process (in default, all of tracks will be
  // read and fitted)
  conf_alignment.alignedTransformUpdater
    = [](Acts::DetectorElementBase* detElement,
         const Acts::GeometryContext& gctx,
         const Acts::Transform3D& aTransform) {
        Telescope::TelescopeDetectorElement* telescopeDetElement =
          dynamic_cast<Telescope::TelescopeDetectorElement*>(detElement);
        if (telescopeDetElement) {
          telescopeDetElement->addAlignedTransform
            (std::make_unique<Acts::Transform3D>(aTransform));
          return true;
        }
        return false;
      };
  
  // Set up the detector elements to be aligned (fix the first one)
  std::vector<Acts::DetectorElementBase*> dets;
  unsigned int idet = 0;
  for (const auto& det : element_col) {
    idet++;
    // Skip the first detector element
    if (idet == 1) {
      continue;
    }
    dets.push_back(det.get());
  }
  conf_alignment.alignedDetElements = std::move(dets);
  
   // The criteria to determine if the iteration has converged. @Todo: to use
  // delta chi2 instead
  conf_alignment.chi2ONdfCutOff = 0.05;
  conf_alignment.deltaChi2ONdfCutOff = {10,0.00001};
  // The maximum number of iterations
  conf_alignment.maxNumIterations = 400;
  // set up the alignment dnf for each iteration
  std::map<unsigned int, std::bitset<6>> iterationState;
  for (unsigned int iIter = 0; iIter < conf_alignment.maxNumIterations; iIter++) {
    std::bitset<6> mask(std::string("111111"));
    if (iIter % 2 == 0 ) {
      // fix the x offset (i.e. offset along the beam) and rotation around y
      mask = std::bitset<6>(std::string("010110"));
    }else {
      // fix the x offset and rotation aroundi x, z
      mask = std::bitset<6>(std::string("101001"));
    }
    iterationState.emplace(iIter, mask);
  }
  conf_alignment.iterationState = std::move(iterationState);
  conf_alignment.align = Telescope::TelescopeAlignmentAlgorithm::makeAlignmentFunction
                        (trackingGeometry, magneticField, Acts::Logging::INFO);

  
  auto algAlign = std::make_shared<Telescope::TelescopeAlignmentAlgorithm>(conf_alignment, Acts::Logging::INFO);


  std::FILE* fp = std::fopen(datafile_name.c_str(), "r");
  if(!fp) {
    std::fprintf(stderr, "File opening failed: %s \n", datafile_name.c_str());
    throw;
  }


  char readBuffer[1000000];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  rapidjson::GenericDocument<rapidjson::UTF8<char>, rapidjson::CrtAllocator>  doc;
  // rapidjson::Document doc;
  doc.ParseStream(is);
  if(!doc.IsArray() || !doc.GetArray().Size()){
    std::fprintf(stderr, "no, it is not data array\n");
    throw;
  }

  while(1){
    FW::WhiteBoard wb;
    FW::AlgorithmContext ctx(0, 0, wb);

    // new sourcelink
    // Read input data
    std::vector<std::vector<Telescope::PixelSourceLink>> sourcelinkTracks;

    auto ev_it = doc.Begin();
    auto ev_it_end = doc.End();  
    uint64_t processed_datapack_count = 0;
    uint64_t good_datapack_count = 0;
    while(ev_it != ev_it_end && good_datapack_count < 20000){
      const auto &evpack = *ev_it;
      ev_it++;
      processed_datapack_count ++;
    
      const auto &frames = evpack["layers"];
    
      bool is_good_datapack = true;
      for(const auto& layer : evpack["layers"].GetArray()){
        uint64_t l_hit_n = layer["hit"].GetArray().Size();
        if(l_hit_n != 1){
          is_good_datapack = false;
          continue;
        }
      }
      if(!is_good_datapack){
        continue;
      }    
      good_datapack_count ++;

      // Setup local covariance
      Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
      cov(0, 0) = 50_um * 50_um;
      cov(1, 1) = 50_um * 50_um;
      // Create the track sourcelinks
      std::vector<Telescope::PixelSourceLink> sourcelinks;
      
      for(size_t i= 0; i< 6; i++){
        double x = frames[i]["hit"][0]["pos"][0].GetDouble() - 0.02924*1024/2.0;
        double y = frames[i]["hit"][0]["pos"][1].GetDouble() - 0.02688*512/2.0;
        Acts::Vector2D loc;
        loc << x, y;
        sourcelinks.emplace_back(*surface_col.at(i), loc, cov);
      }
      sourcelinkTracks.push_back(sourcelinks);
    }
    
    std::cout << "There are " << sourcelinkTracks.size() << " tracks read-in"<< std::endl;    
    // ctx.eventStore.add("TeleSourcelinkTracks", std::move(sourcelinkTracks));



    // Prepare the initial track parameters collection
    std::vector<Acts::CurvilinearParameters> initialParameters;
    initialParameters.reserve(sourcelinkTracks.size());
    unsigned int iTrack = 0;
    while (iTrack < sourcelinkTracks.size()) {
      // Create initial parameters
      // The position is taken from the first measurement
      const auto& sourcelinks = sourcelinkTracks.at(iTrack);
      const Acts::Vector3D global0 =
        sourcelinks.at(0).globalPosition(ctx.geoContext);
      const Acts::Vector3D global1 =
        sourcelinks.at(1).globalPosition(ctx.geoContext);
      Acts::Vector3D distance = global1 - global0;

      const double phi = Acts::VectorHelpers::phi(distance);
      const double theta = Acts::VectorHelpers::theta(distance);

      // shift along the beam by 100_mm
      Acts::Vector3D rPos = global0 - distance / 2;
      Acts::Vector3D rMom(6_GeV * sin(theta) * cos(phi),
                          6_GeV * sin(theta) * sin(phi), 6_GeV * cos(theta));

      Acts::BoundSymMatrix cov;
      cov <<
        std::pow(50_um, 2), 0., 0.,   0.,   0.,     0.,
        0., std::pow(50_um, 2), 0.,   0.,   0.,     0.,
        0.,                 0., 0.01, 0.,   0.,     0.,
        0.,                 0., 0.,   0.01, 0.,     0.,
        0.,                 0., 0.,   0.,   0.0001, 0.,
        0.,                 0., 0.,   0.,   0.,     1.;

    
      Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(
                                                                         cov, rPos, rMom, 1., 0);
      initialParameters.push_back(rStart);
      iTrack++;
    }

    auto pSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(
                                                                  Acts::Vector3D{0., 0., 0.}, Acts::Vector3D{1., 0., 0.});

  
    // Set the KalmanFitter options
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
                                                                 ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
                                                                 Acts::VoidOutlierFinder(), pSurface.get());

    // Set the alignment options
    FW::AlignmentOptions<Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>
      alignOptions(kfOptions, conf_alignment.alignedTransformUpdater,
                   conf_alignment.alignedDetElements, conf_alignment.chi2ONdfCutOff, conf_alignment.deltaChi2ONdfCutOff,
                   conf_alignment.maxNumIterations, conf_alignment.iterationState);
  
    std::printf("Invoke alignment");
    auto result = conf_alignment.align(sourcelinkTracks, initialParameters, alignOptions);

    if (result.ok()) {
      std::cout<<"Alignment finished with deltaChi2 = " << result.value().deltaChi2;
      // ctx.eventStore.add("TeleAlignResult", std::move( result.value() ) );
    
    } else {
      std::cout<<"Alignment failed with " << result.error();
    }    
    
    if (algAlign->execute(++ctx) != FW::ProcessCode::SUCCESS) {
      throw std::runtime_error("Failed to process event data");
    }

    // auto &result = wb.get<FW::AlignmentResult>("TeleAlignResult");
    // std::cout<< "main Alignment finished with deltaChi2 = " << result.deltaChi2<< std::endl;    

    /////////////////////////////
    
    break;
  }  
  return 0;
}
