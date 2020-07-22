// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include "Acts/Utilities/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"

#include "ACTFW/EventData/Track.hpp"

#include "ACTFW/Geometry/CommonGeometry.hpp"

#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"


#include "ObjTelescopeTrackWriter.hpp"
#include "RootTelescopeTrackWriter.hpp"
#include "TelescopeAlignmentAlgorithm.hpp"
#include "BuildTelescopeDetector.hpp"
#include "TelescopeDetectorElement.hpp"
#include "TelescopeFittingAlgorithm.hpp"
#include "TelescopeTrackingPerformanceWriter.hpp"

#include "TelescopeTrackReader.hpp"

#include "getopt.h"

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
  Telescope::TelescopeDetectorElement::ContextType nominalContext;
  std::vector<std::shared_ptr<Telescope::TelescopeDetectorElement>> detectorStore;
  std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry
    = Telescope::buildDetector<Telescope::TelescopeDetectorElement>(nominalContext, detectorStore, matDeco);

  
  
  // Get the surfaces;
  std::vector<const Acts::Surface*> surfaces;
  trackingGeometry->visitSurfaces
    ([&](const Acts::Surface* s) {
       if (s and s->associatedDetectorElement()) {
         surfaces.push_back(s);
       }
     });

  // The source link tracks reader
  Telescope::TelescopeTrackReader trackReader;
  trackReader.detectorSurfaces = surfaces;  

  // setup the alignment algorithm
  Telescope::TelescopeAlignmentAlgorithm::Config conf_alignment;
  //@Todo: add run number information in the file name
  conf_alignment.inputFileName = inputDir + datafile_name;
  conf_alignment.outputTrajectories = "trajectories";
  conf_alignment.trackReader = trackReader;
  // The number of tracks you want to process (in default, all of tracks will be
  // read and fitted)
  conf_alignment.maxNumTracks = 20000;
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
  for (const auto& det : detectorStore) {
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
  conf_alignment.chi2ONdfCutOff = 0.005;
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

  FW::Sequencer::Config conf_seq;
  conf_seq.logLevel  = logLevel;
  conf_seq.outputDir =  outputDir;  
  conf_seq.numThreads = 1;
  conf_seq.events = 1;

  FW::Sequencer seq(conf_seq);
  
  seq.addAlgorithm(std::make_shared<Telescope::TelescopeAlignmentAlgorithm>(conf_alignment, Acts::Logging::INFO));
  // seq.addAlgorithm(std::make_shared<Telescope::TelescopeFittingAlgorithm>(conf_fitter, logLevel));

  seq.run();
  
  
  return 0;
}


  // // setup the fitting algorithm
  // TelescopeFittingAlgorithm::Config conf_fitter;
  // //@Todo: add run number information in the file name
  // conf_fitter.inputFileName = inputDir + "/alpide-data.json";
  // conf_fitter.outputTrajectories = "trajectories";
  // conf_fitter.trackReader = trackReader;
  // // The number of tracks you want to process (in default, all of tracks will be
  // // read and fitted)
  // conf_fitter.maxNumTracks = 20000;
  // conf_fitter.fit = TelescopeFittingAlgorithm::makeFitterFunction
  //   (trackingGeometry, magneticField, logLevel);


  // // write tracks as root tree
  // RootTelescopeTrackWriter::Config conf_trackRootWriter;
  // conf_trackRootWriter.inputTrajectories = conf_fitter.outputTrajectories;
  // conf_trackRootWriter.outputDir = outputDir;
  // conf_trackRootWriter.outputFilename = "telescope_tracks.root";
  // conf_trackRootWriter.outputTreename = "tracks";

  // // write the tracks (measurements only for the moment) as Csv
  // Obj::ObjTelescopeTrackWriter::Config conf_trackObjWriter;
  // conf_trackObjWriter.inputTrajectories = conf_fitter.outputTrajectories;
  // conf_trackObjWriter.outputDir = outputDir;
  // // The number of tracks you want to show (in default, all of tracks will be
  // // shown)
  // conf_trackObjWriter.maxNumTracks = 100;
  
  // // write reconstruction performance data
  // TelescopeTrackingPerformanceWriter::Config conf_perfFitter;
  // conf_perfFitter.inputTrajectories = conf_fitter.outputTrajectories;
  // conf_perfFitter.outputDir = outputDir;



  // seq.addWriter(std::make_shared<RootTelescopeTrackWriter>(conf_trackRootWriter, logLevel));
  // seq.addWriter(std::make_shared<Obj::ObjTelescopeTrackWriter>(conf_trackObjWriter, logLevel));
  // seq.addWriter(std::make_shared<TelescopeTrackingPerformanceWriter>(conf_perfFitter, logLevel));
