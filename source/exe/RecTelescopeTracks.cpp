// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Utilities/Units.hpp"

#include "ObjTelescopeTrackWriter.hpp"
#include "RootTelescopeTrackWriter.hpp"
#include "TelescopeAlignmentAlgorithm.hpp"
#include "TelescopeDetector.hpp"
#include "TelescopeDetectorElement.hpp"
#include "TelescopeFittingAlgorithm.hpp"
#include "TelescopeTrackingPerformanceWriter.hpp"

#include "TelescopeTrackReader.hpp"

using namespace Acts::UnitLiterals;
using namespace FW;

int main(int argc, char* argv[]) {
  TelescopeDetector detector;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  //  detector.addOptions(desc);
  Options::addBFieldOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd =
      std::make_shared<FW::RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Setup detector geometry
  auto geometry = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }

  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);

  // Get the surfaces;
  std::vector<const Acts::Surface*> surfaces;
  surfaces.reserve(6);
  trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // The source link tracks reader
  TelescopeTrackReader trackReader;
  trackReader.detectorSurfaces = surfaces;

  // setup the alignment algorithm
  TelescopeAlignmentAlgorithm::Config alignment;
  //@Todo: add run number information in the file name
  alignment.inputFileName = inputDir + "/alpide-data.json";
  alignment.outputTrajectories = "trajectories";
  alignment.trackReader = trackReader;
  // The number of tracks you want to process (in default, all of tracks will be
  // read and fitted)
  alignment.maxNumTracks = 20000;
  alignment.alignedTransformUpdater = [](Acts::DetectorElementBase* detElement,
                                         const Acts::GeometryContext& gctx,
                                         const Acts::Transform3D& aTransform) {
    Telescope::TelescopeDetectorElement* telescopeDetElement =
        dynamic_cast<Telescope::TelescopeDetectorElement*>(detElement);
    auto alignContext =
        std::any_cast<Telescope::TelescopeDetectorElement::ContextType>(gctx);
    if (telescopeDetElement) {
      telescopeDetElement->addAlignedTransform(
          std::make_unique<Acts::Transform3D>(aTransform), alignContext.iov);
      return true;
    }
    return false;
  };
  // Set up the detector elements to be aligned (fix the first one)
  std::vector<Acts::DetectorElementBase*> dets;
  dets.reserve(detector.detectorStore.size());
  unsigned int idet = 0;
  for (const auto& det : detector.detectorStore) {
    idet++;
    // Skip the first detector element
    if (idet == 1) {
      continue;
    }
    dets.push_back(det.get());
  }
  std::cout << "There are " << dets.size() << " detector elements to be aligned"
            << std::endl;
  alignment.alignedDetElements = std::move(dets);
   // The criteria to determine if the iteration has converged. @Todo: to use
  // delta chi2 instead
  alignment.chi2ONdfCutOff = 0.1;
  // The maximum number of iterations
  alignment.maxNumIterations = 160;
  // set up the alignment dnf for each iteration
  std::map<unsigned int, std::bitset<6>> iterationState;
  for (unsigned int iIter = 0; iIter < alignment.maxNumIterations; iIter++) {
    std::bitset<6> mask(std::string("111111"));
    if (iIter % 4 == 0 or iIter % 4 == 1) {
      // fix the x offset (i.e. offset along the beam) and rotation around y
      mask = std::bitset<6>(std::string("101110"));
    } else if (iIter % 4 == 2) {
      // align only the x offset (the x offset and y, z offset could not be
      // aligned together)
      mask = std::bitset<6>(std::string("000001"));
    } else {
      // fix the x offset and rotation around x, z
      mask = std::bitset<6>(std::string("010110"));
    }
    iterationState.emplace(iIter, mask);
  }
  alignment.iterationState = std::move(iterationState);
  alignment.align = TelescopeAlignmentAlgorithm::makeAlignmentFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<TelescopeAlignmentAlgorithm>(alignment, logLevel));

  // setup the fitting algorithm
  TelescopeFittingAlgorithm::Config fitter;
  //@Todo: add run number information in the file name
  fitter.inputFileName = inputDir + "/alpide-data.json";
  fitter.outputTrajectories = "trajectories";
  fitter.trackReader = trackReader;
  // The number of tracks you want to process (in default, all of tracks will be
  // read and fitted)
  fitter.maxNumTracks = 20000;
  fitter.fit = TelescopeFittingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<TelescopeFittingAlgorithm>(fitter, logLevel));

  // write tracks as root tree
  RootTelescopeTrackWriter::Config trackRootWriter;
  trackRootWriter.inputTrajectories = fitter.outputTrajectories;
  trackRootWriter.outputDir = outputDir;
  trackRootWriter.outputFilename = "telescope_tracks.root";
  trackRootWriter.outputTreename = "tracks";
  sequencer.addWriter(
      std::make_shared<RootTelescopeTrackWriter>(trackRootWriter, logLevel));

  if (vm["output-obj"].template as<bool>()) {
    // write the tracks (measurements only for the moment) as Csv
    Obj::ObjTelescopeTrackWriter::Config trackObjWriter;
    trackObjWriter.inputTrajectories = fitter.outputTrajectories;
    trackObjWriter.outputDir = outputDir;
    // The number of tracks you want to show (in default, all of tracks will be
    // shown)
    trackObjWriter.maxNumTracks = 100;
    sequencer.addWriter(std::make_shared<Obj::ObjTelescopeTrackWriter>(
        trackObjWriter, logLevel));
  }

  // write reconstruction performance data
  TelescopeTrackingPerformanceWriter::Config perfFitter;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.outputDir = outputDir;
  sequencer.addWriter(std::make_shared<TelescopeTrackingPerformanceWriter>(
      perfFitter, logLevel));

  return sequencer.run();
}
