#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "StraightLinePropagator.hpp"

#include "JsonGenerator.hpp"

#include <filesystem>
#include "getopt.h"


using namespace Acts::UnitLiterals;


static const std::string help_usage = R"(
Usage:
  -help                        help message
  -verbose                     verbose flag
  -eventMax       [int]        max number of events
  -targetLayerID  [int]        target layer id
  -fitData        [jsonfile]   input fitted data file
  -oriData        [jsonfile]   input origin data file
  -geomerty       [jsonfile]   geometry input file
  -outputdir      [path]       output dir path
)";



int main(int argc, char* argv[]){

  int do_help = false;
  int do_verbose = false;

  struct option longopts[] =
    {
     { "help",           no_argument,       &do_help,      1  },
     { "verbose",        no_argument,       &do_verbose,   1  },
     { "eventMax",       required_argument, NULL,           'm' },
     { "targetLayerID",  required_argument, NULL,           't' },
     { "fitData",        required_argument, NULL,           'f' },
     { "oriData",        required_argument, NULL,           'r' },
     { "geomerty",       required_argument, NULL,      'g' },
     { "outputdir",      required_argument, NULL,      'o' },
     { 0, 0, 0, 0 }};


  size_t target_layer_id = 6;
  size_t eventMaxNum = 10;
  std::string oridatafile_name;
  std::string fitdatafile_name;
  std::string geofile_name;
  std::string outputDir;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL))!= -1) {
    switch (c) {
    case 'm':
      eventMaxNum = std::stoul(optarg);
      break;
    case 't':
      target_layer_id = std::stoul(optarg);
      break;
    case 'f':
      fitdatafile_name = optarg;
      break;
    case 'r':
      oridatafile_name = optarg;
      break;
    case 'g':
      geofile_name = optarg;
      break;
    case 'o':
      outputDir = optarg;
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
      fprintf(stderr,"unrecognized option\n  please check the options usage\n");
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(1);
      break;
    default:
      fprintf(stderr,"case default, missing branch in switch-case\n");
      exit(1);
      break;
    }
  }

  if(do_help){
    std::fprintf(stdout, "help message\n");
    std::fprintf(stdout, "%s\n", help_usage.c_str());
    exit(0);
  }

  Telescope::JsonGenerator gen_ori(oridatafile_name);
  Telescope::JsonGenerator gen_fit(fitdatafile_name);

  Acts::StraightLineStepper stepper;
  Telescope::StraightLinePropagator propagator(std::move(stepper));
  Telescope::StraightLinePropagator::PropagatorOptions options;
  Acts::GeometryContext geoContext;
  Acts::MagneticFieldContext magContext;

  // Construct a target plane surface
  std::map<size_t, std::array<double, 6>> geoconf;
  geoconf = Telescope::JsonGenerator::ReadGeoFromGeoFile(geofile_name);
  std::array<double, 6> target_layer_geoconf = geoconf.at(target_layer_id);
  double cx = target_layer_geoconf[0];
  double cy = target_layer_geoconf[1];
  double cz = target_layer_geoconf[2];
  Acts::Vector3D translation(cx,cy,cz);

  double rx = target_layer_geoconf[3];
  double ry = target_layer_geoconf[4];
  double rz = target_layer_geoconf[5];
  // The rotation around global z axis
  Acts::AngleAxis3D rotZ(rz, Acts::Vector3D::UnitZ());
  // The rotation around global y axis
  Acts::AngleAxis3D rotY(ry, Acts::Vector3D::UnitY());
  // The rotation around global x axis
  Acts::AngleAxis3D rotX(rx, Acts::Vector3D::UnitX());
  Acts::Rotation3D rotation = rotZ * rotY * rotX;

  auto trafo = std::make_shared<Acts::Transform3D>(Acts::Translation3D(translation) * rotation);
  auto target_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo,
                                                                      std::make_shared<const Acts::RectangleBounds>(30_mm, 15_mm));

  for(size_t n = 0; n<eventMaxNum; n++){
    std::cout<<"\n\n"<<std::endl;
    std::vector<Acts::CurvilinearParameters> track_curPara_v;
    {
      Telescope::JsonAllocator jsa;
      Telescope::JsonDocument jsdoc(&jsa);
      // while(1){// doc is cleared at beginning of each loop
      Telescope::JsonValue js_pack;
      jsdoc.Populate(gen_fit);
      if(!gen_fit.isvalid){
        return 0;
      }
      jsdoc.Swap(js_pack);

      Acts::BoundSymMatrix cov;
      Acts::Vector3D pos;
      Acts::Vector3D mom;
      double charge;
      double time;
      const auto& js_tracks = js_pack["tracks"];
      for(const auto& js_track: js_tracks.GetArray()){
        const auto& js_states = js_track["states"];
        for(const auto& js_state: js_states.GetArray()){
          size_t id = js_state["id"].GetUint();
          if(id != target_layer_id -1 ){
            continue;
          }

          double x = js_state["x"].GetDouble();
          double y = js_state["y"].GetDouble();
          double z = js_state["z"].GetDouble();
          pos = Acts::Vector3D(x, y, z);

          double px = js_state["px"].GetDouble();
          double py = js_state["py"].GetDouble();
          double pz = js_state["pz"].GetDouble();
          mom = Acts::Vector3D(px, py, pz);

          time = js_state["t"].GetDouble();
          charge = js_state["q"].GetDouble();

          const auto& js_cov = js_state["cov"];
          std::vector<double> cov_data;
          for(const auto& js_e : js_cov.GetArray()){
            double e = js_e.GetDouble();
            cov_data.push_back(e);
          }
          cov = Acts::BoundSymMatrix(cov_data.data());
          track_curPara_v.emplace_back(cov, pos, mom, charge, time);
        }
      }
    }

    // get the bound parameters at the target surface
    std::vector<Acts::BoundParameters>  target_boundPara_v;
    for(auto & curPara : track_curPara_v){
      auto targetParams = propagator.transport(geoContext, magContext, options, curPara, *target_surface);
      target_boundPara_v.push_back( std::move(targetParams) );
    }

    std::vector<Acts::Vector2D> hit_local_v;
    std::vector<Acts::Vector3D> hit_global_v;
    {

      Telescope::JsonAllocator jsa;
      Telescope::JsonDocument jsdoc(&jsa);
      // while(1){// doc is cleared at beginning of each loop
      Telescope::JsonValue js_pack;
      jsdoc.Populate(gen_ori);
      if(!gen_ori.isvalid){
        return 0;
      }
      jsdoc.Swap(js_pack);

      const auto &layers = js_pack["layers"];
      for(const auto&layer : layers.GetArray()){
        size_t id_ext = layer["ext"].GetUint();
        if(id_ext != target_layer_id){
          continue;
        }
        for(const auto&hit : layer["hit"].GetArray()){
          double x_hit = hit["pos"][0].GetDouble() - 0.02924*1024/2.0;
          double y_hit = hit["pos"][1].GetDouble() - 0.02688*512/2.0;
          Acts::Vector2D lpos;
          Acts::Vector3D mom;
          Acts::Vector3D gpos;
          lpos<< x_hit, y_hit;
          target_surface->localToGlobal(geoContext, lpos, mom, gpos);
          hit_local_v.push_back(lpos);
          hit_global_v.push_back(gpos);
        }
      }
    }

    // std::cout<<"=======starting parameter info:====="<<std::endl;
    // std::cout<<curPara<<std::endl;
    // std::cout<<"position(): \n"<< curPara.position()<<std::endl;
    // std::cout<<"momentum(): \n"<< curPara.momentum()<<std::endl;
    // std::cout<<"covariance: \n"<< *curPara.covariance()<<std::endl;
    std::cout<<"======target parameter info:====="<<std::endl;
    // std::cout<<"position(): \n"<< targetParams.position()<<std::endl;
    // std::cout<<"momentum(): \n"<< targetParams.momentum()<<std::endl;
    // std::cout<<"covariance: \n"<< *targetParams.covariance()<<std::endl;
    for(auto &boundPara : target_boundPara_v){
      std::cout<< boundPara<<std::endl;
    }

    std::cout<<"======ori hit info:====="<<std::endl;
    for(auto &gpos : hit_global_v){
      std::cout<<"global hit postion\n "<< gpos<<std::endl;
    }
  }
  return 0;
}
