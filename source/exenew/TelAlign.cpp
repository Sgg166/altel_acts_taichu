#include "TelActs.hh"
#include "getopt.h"

using namespace Acts::UnitLiterals;

static const std::string help_usage = R"(
Usage:
  -help                    help message
  -verbose                 verbose flag
  -hitFile        [PATH]   name of data json file
  -energy         [float]  beam energy, GeV
  -outputGeometry [PATH]   alignment result file
  -inputGeometry  [PATH]   geometry input file
  -hitResX        [float]  preset detector hit resolution X
  -hitResY        [float]  preset detector hit resolution Y
  -seedResX       [float]  preset detector seed resolution X
  -seedResY       [float]  preset detector seed resolution Y
  -seedResPhi     [float]  preset seed track resolution Phi
  -seedResTheta   [float]  preset seed track resolution Theta
  -maxItera       [int]    max time of alignement iterations
  -deltaItera     [int]    converge condition: delta iteration
  -deltaChi2      [float]  converge condition: delta Chi2ONdf
  -chi2ONdfCutOff [float]
)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"hitFile", required_argument, NULL, 'f'},
                              {"energy", required_argument, NULL, 'e'},
                              {"outputGeometry", required_argument, NULL, 'o'},
                              {"inputGeomerty", required_argument, NULL, 'g'},
                              {"hitResX", required_argument, NULL, 'r'},
                              {"hitResY", required_argument, NULL, 's'},
                              {"seedResX", required_argument, NULL, 'u'},
                              {"seedResY", required_argument, NULL, 'q'},
                              {"seedResPhi", required_argument, NULL, 't'},
                              {"seedResTheta", required_argument, NULL, 'w'},
                              {"maxItera", required_argument, NULL, 'x'},
                              {"deltaItera", required_argument, NULL, 'y'},
                              {"deltaChi2", required_argument, NULL, 'z'},
                              {"chi2ONdfCutOff", required_argument, NULL, 'p'},
                             {0, 0, 0, 0}};

  std::string datafile_name;
  std::string outputfile_name;
  std::string geofile_name;
  double resX = 150_um;
  double resY = 150_um;

 // Use large starting parameter covariance
  double resLoc1 = 20_um;
  double resLoc2 = 20_um;
  double resPhi = 0.02;
  double resTheta = 0.02;

  // Iterations converge criteria
  size_t maxNumIterations = 400;
  size_t nIterations = 10;
  double deltaChi2ONdf = 1e-5;
  double chi2ONdfCutOff = 0.0001;

  double beamEnergy = -1;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
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
    case 'u':
      resLoc1 = std::stod(optarg);
      break;
    case 'q':
      resLoc2 = std::stod(optarg);
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
    case 'p':
      chi2ONdfCutOff = std::stod(optarg);
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

  if (datafile_name.empty() || outputfile_name.empty()) {
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

  // Set up surfaces
  std::map<size_t, std::shared_ptr<const Acts::Surface>> surfaces_selected;
  for (const auto &e : element_col) {
    auto id = e->telDetectorID();
    surfaces_selected[id] = e->surface().getSharedPtr();
  }

  // Set up the detector elements to be aligned (fix the first one)
  std::vector<Acts::DetectorElementBase *> element_col_align;
  for (auto it = element_col.begin(); it != element_col.end(); ++it) {
    if (it != element_col.begin())
      element_col_align.push_back((*it).get());
  }
  //

  /////////////////////////////////////
  // The criteria to determine if the iteration has converged. @Todo: to use
  // delta chi2 instead
  std::pair<size_t, double> deltaChi2ONdfCutOff = {nIterations, deltaChi2ONdf};
  // set up the alignment dnf for each iteration
  std::map<unsigned int, std::bitset<6>> iterationState;
  for (unsigned int iIter = 0; iIter < maxNumIterations; iIter++) {
    std::bitset<6> mask(std::string("111111"));

    // mask = std::bitset<6>(std::string("000011"));
    if (iIter % 2 == 0) {
      // only align offset along x/y
      mask = std::bitset<6>(std::string("000011"));
    } else {
      // only align rotation around beam
      mask = std::bitset<6>(std::string("100000"));
    }
    iterationState.emplace(iIter, mask);
  }

  ActsAlignment::AlignedTransformUpdater alignedTransformUpdaterFun =
      [](Acts::DetectorElementBase *detElement,
         const Acts::GeometryContext &gctx,
         const Acts::Transform3D &aTransform) {
        TelActs::TelescopeDetectorElement *telescopeDetElement =
            dynamic_cast<TelActs::TelescopeDetectorElement *>(detElement);
        if (telescopeDetElement) {
          telescopeDetElement->addAlignedTransform(
              std::make_unique<Acts::Transform3D>(aTransform));
          return true;
        }
        return false;
      };

  auto alignFun = TelActs::makeAlignmentFunction(
      trackingGeometry, magneticField,
      do_verbose ? (Acts::Logging::VERBOSE) : (Acts::Logging::INFO));

  {
    // Setup local covariance
    Acts::BoundMatrix cov_hit = Acts::BoundMatrix::Zero();
    cov_hit(0, 0) = resX * resX;
    cov_hit(1, 1) = resY * resY;

    //
    size_t n_datapack_select_opt = 20000;
    JsonFileDeserializer jsf(datafile_name);
    std::vector<std::vector<TelActs::PixelSourceLink>> sourcelinkTracks;
    while (jsf) {
      if (sourcelinkTracks.size() > n_datapack_select_opt) {
        break;
      }
      auto evpack = jsf.getNextJsonDocument();
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

      // drop multiple hits events
      if (sourcelinks.size() != surfaces_selected.size()) {
        continue;
      }
      std::vector<size_t> geo_ids;
      bool found_same_geo_id = false;
      for (const auto &sl : sourcelinks) {
        size_t geo_id = sl.referenceSurface().geometryId().value();
        if (std::find(geo_ids.begin(), geo_ids.end(), geo_id) !=
            geo_ids.end()) {
          found_same_geo_id = true;
          break;
        }
        geo_ids.push_back(geo_id);
      }
      if (found_same_geo_id) {
        continue;
      }
      sourcelinkTracks.push_back(sourcelinks);
    }

    // seeding
    std::vector<Acts::CurvilinearTrackParameters> initialParameters;
    initialParameters.reserve(sourcelinkTracks.size());

    for (const auto &sourcelinks : sourcelinkTracks) {
      // Create initial parameters
      // The position is taken from the first measurement
      const Acts::Vector3D global0 = sourcelinks.at(0).globalPosition(gctx);
      const Acts::Vector3D global1 = sourcelinks.at(1).globalPosition(gctx);
      Acts::Vector3D distance = global1 - global0;

      const double phi = Acts::VectorHelpers::phi(distance);
      const double theta = Acts::VectorHelpers::theta(distance);
      Acts::Vector3D rPos = global0;
      Acts::Vector4D rPos4(rPos.x(), rPos.y(), rPos.z(), 0);

      Acts::BoundSymMatrix cov_seed;
      cov_seed << resLoc1 * resLoc1, 0., 0., 0., 0., 0., 0., resLoc2 * resLoc2,
        0., 0., 0., 0., 0., 0., resPhi * resPhi, 0., 0., 0., 0., 0., 0.,
        resTheta * resTheta, 0., 0., 0., 0., 0., 0., 0.0001, 0., 0., 0., 0., 0.,
        0., 1.;

      Acts::CurvilinearTrackParameters rStart(rPos4, phi, theta, 1., beamEnergy,
                                              cov_seed);
      initialParameters.push_back(rStart);
    }

    auto kfLogger = Acts::getDefaultLogger("KalmanFilter",
                                           do_verbose ? (Acts::Logging::VERBOSE)
                                                      : (Acts::Logging::INFO));

    // Set the KalmanFitter options
    Acts::PropagatorPlainOptions pOptions;
    pOptions.mass = 0.511 * Acts::UnitConstants::MeV;
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
        gctx, mctx, cctx, Acts::VoidOutlierFinder(),
        Acts::LoggerWrapper{*kfLogger},pOptions
        ); // pSurface default nullptr

    // Set the alignment options
    ActsAlignment::AlignmentOptions<
        Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>>
        alignOptions(kfOptions, alignedTransformUpdaterFun, element_col_align,
                     chi2ONdfCutOff, deltaChi2ONdfCutOff, maxNumIterations,
                     iterationState);

    std::printf("Invoke alignment");
    auto result = alignFun(sourcelinkTracks, initialParameters, alignOptions);

    if (!result.ok()) {
      std::printf("Alignment failed with %s \n",
                  result.error().message().c_str());
    }

    std::map<size_t, std::array<double, 6>> geo_result;
    auto &js_dets = jsd_geo["geometry"]["detectors"];
    for (const auto &det : element_col) {
      const auto &surface = &det->surface();
      const auto &transform = det->transform(gctx);
      const auto &translation = transform.translation();
      const auto &rotation = transform.rotation();
      const Acts::Vector3D rotAngles = rotation.eulerAngles(2, 1, 0);

      size_t id = det->telDetectorID();
      double cx = translation.x();
      double cy = translation.y();
      double cz = translation.z();
      double rx = rotAngles(2);
      double ry = rotAngles(1);
      double rz = rotAngles(0);
      geo_result[id] = {cx, cy, cz, rx, ry, rz};

      std::printf("layer: %lu   centerX: %f   centerY: %f   centerZ: %f  "
                  "rotationX: %f   rotationY: %f   rotationZ: %f\n",
                  id, cx, cy, cz, rx, ry, rz);

      //update geo js
      for(auto &js_det : js_dets.GetArray()){
        if(js_det["id"].GetUint() == id){
          js_det["center"]["x"]=cx;
          js_det["center"]["y"]=cy;
          js_det["center"]["z"]=cz;
          js_det["rotation"]["x"]=rx;
          js_det["rotation"]["y"]=ry;
          js_det["rotation"]["z"]=rz;
          break;
        }
      }
    }
    std::printf("select %lu datapacks (single track events) \n",
                sourcelinkTracks.size());
  }

  // JsonUtils::printJsonValue(jsd_geo, true);
  std::string jsstr = JsonUtils::stringJsonValue(jsd_geo, true);
  std::FILE *fp = std::fopen(outputfile_name.c_str(), "w");
  std::fwrite(jsstr.data(), 1, jsstr.size(), fp);
  std::fclose(fp);
  return 0;
}
