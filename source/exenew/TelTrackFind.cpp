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




/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] direction Normalised direction vector
///
/// @return Projection Jacobian
Acts::FreeToBoundMatrix
freeToCurvilinearJacobian(const Acts::Vector3D &direction) {
  // Optimized trigonometry on the propagation direction
  const double x = direction(0); // == cos(phi) * sin(theta)
  const double y = direction(1); // == sin(phi) * sin(theta)
  const double z = direction(2); // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  Acts::FreeToBoundMatrix jacToCurv = Acts::FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < Acts::s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 4) = cosPhi * cosTheta;
  jacToCurv(3, 5) = sinPhi * cosTheta;
  jacToCurv(3, 6) = -sinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}



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
  std::map<std::shared_ptr<const Acts::Surface>, size_t> surface_id_map;
  for (const auto &e : element_col) {
    auto id = e->telDetectorID();
    auto surface = e->surface().getSharedPtr();
    surfaces_selected[id] = surface;
    surface_id_map[surface] = id;
  }

  /////////////////////////////////////
  Acts::Logging::Level logLevel =
      do_verbose ? (Acts::Logging::VERBOSE) : (Acts::Logging::INFO);

  //////////// hit data
  Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
  cov_hit(0, 0) = resX * resX;
  cov_hit(1, 1) = resY * resY;

  /////////////  track seed conf
  size_t seedSurfaceGeoIDStart = surfaces_selected[0]->geometryId().value();
  size_t seedSurfaceGeoIDEnd = surfaces_selected[1]->geometryId().value();

  Acts::BoundSymMatrix cov_seed;
  cov_seed << seedResX * seedResX, 0., 0., 0., 0., 0., 0.,
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
  JsonFileDeserializer jsfd(datafile_name);
  JsonFileSerializer jsfs(fitted_datafile_name);

  size_t eventNumber = 0;
  while(jsfd && (eventNumber< eventMaxNum || eventMaxNum<0)){
    auto evpack = jsfd.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }

    std::vector<TelActs::PixelSourceLink> sourcelinks;
    const auto &layers = evpack["layers"];
    for (const auto &layer : layers.GetArray()) {
      size_t id_ext = layer["ext"].GetUint();
      auto surface_it = surfaces_selected.find(id_ext);
      if (surface_it == surfaces_selected.end()) {
        continue;
      }
      for (const auto &hit : layer["hit"].GetArray()) {
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        Acts::Vector2D loc_hit;
        loc_hit << x_hit, y_hit;
        sourcelinks.emplace_back(*(surface_it->second), loc_hit, cov_hit);
      }
    }

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
          initialParameters.emplace_back(rPos4, phi, theta, beamEnergy, q,cov_seed);
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


    //////////////////////////////////////////////////////////////

    JsonAllocator jsa;
    JsonValue js_output(rapidjson::kObjectType);
    js_output.AddMember("eventNum", JsonValue(eventNumber), jsa);
    JsonValue js_tracks(rapidjson::kArrayType);

    for (const auto &traj : trajectories) {
      const auto &[trackTips, mj] = traj.trajectory();
      // Loop over all trajectories in a multiTrajectory
      // size of trackTips should <= 1
      for (const size_t &trackTip : trackTips) {
        // Collect the trajectory summary info
        auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
        if (trajState.nMeasurements != surfaces_selected.size()) {
          continue;
        }
        JsonValue js_track(rapidjson::kObjectType);
        JsonValue js_states(rapidjson::kArrayType);
        JsonValue js_states_reverse(rapidjson::kArrayType);
        mj.visitBackwards(trackTip, [&](const auto &state) {
                                      // only fill the track states with non-outlier measurement
                                      auto typeFlags = state.typeFlags();
                                      if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
                                        return true;
                                      }
                                      if (!state.hasSmoothed()) {
                                        return true;
                                      }

                                      auto state_surface = state.referenceSurface().getSharedPtr();
                                      auto surface_it = surface_id_map.find(state_surface);
                                      if (surface_it == surface_id_map.end()) {
                                        return true;
                                      }
                                      size_t layerid = surface_it->second;

                                      // Get the source link info
                                      auto meas = std::get<
                                        Acts::Measurement<TelActs::PixelSourceLink, Acts::BoundIndices,
                                                          Acts::eBoundLoc0, Acts::eBoundLoc1>>(
                                                            *state.uncalibrated());

                                      // Get local position
                                      Acts::Vector2D pos_local(meas.parameters()[Acts::eBoundLoc0],
                                                               meas.parameters()[Acts::eBoundLoc1]);

                                      // The bound parameters info
                                      // Acts::BoundTrackParameters boundpara(state_surface,
                                      //                                state.smoothed(),
                                      //                                state.smoothedCovariance()
                                      //						  );

                                      // 1) Transform bound parameter to free parameter
                                      // 1.1)Fist transform the smoothed bound parameters to free parameters
                                      // to get the position and momentum
                                      Acts::FreeVector freeParams =
                                        Acts::detail::transformBoundToFreeParameters(
                                          *state_surface, gctx, state.smoothed());
                                      // 1.2)Get the global position, direction, q/p, t etc.
                                      Acts::Vector3D pos(freeParams[Acts::eFreePos0],
                                                         freeParams[Acts::eFreePos1],
                                                         freeParams[Acts::eFreePos2]);
                                      Acts::Vector3D dir(freeParams[Acts::eFreeDir0],
                                                         freeParams[Acts::eFreeDir1],
                                                         freeParams[Acts::eFreeDir2]);
                                      double p = std::abs(1 / freeParams[Acts::eFreeQOverP]);
                                      ;
                                      double q = p * freeParams[Acts::eFreeQOverP];
                                      double t = freeParams[Acts::eFreeTime];

                                      /// 1.3) Initialize the jacobian from local to the global frame
                                      Acts::BoundToFreeMatrix jacToGlobal = Acts::BoundToFreeMatrix::Zero();
                                      // Calculate the jacobian
                                      state_surface->initJacobianToGlobal(gctx, jacToGlobal,
                                                                          pos, dir, state.smoothed());
                                      Acts::FreeSymMatrix freeCovariance =
                                        jacToGlobal * state.smoothedCovariance() * jacToGlobal.transpose();

                                      // 2) Transform free parameter to curvilinear parameter
                                      Acts::FreeToBoundMatrix jacToCurv = freeToCurvilinearJacobian(dir);
                                      Acts::BoundSymMatrix curvCovariance =
                                        jacToCurv * freeCovariance * jacToCurv.transpose();

                                      // Write out the curvilinear parameters and its covariance
                                      JsonValue js_track_state(rapidjson::kObjectType); //
                                      js_track_state.AddMember("id", JsonValue(layerid), jsa);
                                      js_track_state.AddMember("x", JsonValue(pos.x()), jsa);
                                      js_track_state.AddMember("y", JsonValue(pos.y()), jsa);
                                      js_track_state.AddMember("z", JsonValue(pos.z()), jsa);
                                      js_track_state.AddMember("lx", JsonValue(pos_local.x()), jsa);
                                      js_track_state.AddMember("ly", JsonValue(pos_local.y()), jsa);
                                      js_track_state.AddMember("dx", JsonValue(dir.x()), jsa);
                                      js_track_state.AddMember("dy", JsonValue(dir.y()), jsa);
                                      js_track_state.AddMember("dz", JsonValue(dir.z()), jsa);
                                      js_track_state.AddMember("p", JsonValue(p), jsa);
                                      js_track_state.AddMember("q", JsonValue(q), jsa);
                                      js_track_state.AddMember("t", JsonValue(t), jsa);

                                      const double *cov_data = curvCovariance.data();
                                      JsonValue js_state_cov(rapidjson::kArrayType); //
                                      js_state_cov.Reserve(Acts::eBoundSize * Acts::eBoundSize, jsa);
                                      for (size_t n = 0; n < Acts::eBoundSize * Acts::eBoundSize; n++) {
                                        js_state_cov.PushBack(JsonValue(*(cov_data + n)), jsa);
                                      }
                                      js_track_state.AddMember("cov", std::move(js_state_cov), jsa);
                                      js_states_reverse.PushBack(std::move(js_track_state), jsa);
                                      return true;
                                    });

        for (size_t i = js_states_reverse.Size(); i > 0; i--) {
          js_states.PushBack(std::move(js_states_reverse[i - 1]), jsa);
        }
        js_track.AddMember("states", std::move(js_states), jsa);
        js_tracks.PushBack(std::move(js_track), jsa);
      }
    }
    js_output.AddMember("tracks", std::move(js_tracks), jsa);

    jsfs.putNextJsonValue(js_output);

  }

  return 0;
}
