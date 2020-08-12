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

using namespace Acts::UnitLiterals;


int main(int argc, char* argv[]){
  Acts::GeometryContext geoContext;

  Acts::MagneticFieldContext magContext;

  Acts::StraightLineStepper stepper;
  Telescope::StraightLinePropagator propagator( std::move(stepper));
  Telescope::StraightLinePropagator::PropagatorOptions options;

  Telescope::JsonAllocator jsa;
  Telescope::JsonDocument jsdoc(&jsa);
  Telescope::JsonGenerator gen("/home/yiliu/testbeam/altel_acts/INSTALL/bin/out7/jswrite.json");
  // while(1){// doc is cleared at beginning of each loop
  Telescope::JsonValue js_pack;
  jsdoc.Populate(gen);
  if(!gen.isvalid){
    return 0;
  }
  jsdoc.Swap(js_pack);

  Acts::BoundSymMatrix cov;
  Acts::Vector3D pos;
  Acts::Vector3D mom;
  double charge;
  double time;
  const auto& js_tracks = js_pack["tracks"];
  for(const auto& js_state: js_tracks.GetArray()){
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
  }

  Acts::CurvilinearParameters curPara(cov, pos, mom, charge, time);

  // Construct a target plane surface
  auto tBounds = std::make_shared<const Acts::RectangleBounds>(27.52512_mm, 13.76256_mm);
  Acts::Translation3D tTranslation{0, 0, 490_mm};
  auto tTransform = std::make_shared<const Acts::Transform3D>(tTranslation);
  auto tSurface = Acts::Surface::makeShared<Acts::PlaneSurface>(tTransform, tBounds);

  // get the bound parameters at the target surface
  auto targetParams = propagator.transport(geoContext, magContext, options, curPara, *tSurface);

  std::cout<<"=======starting parameter info:====="<<std::endl;
  std::cout<<"position(): \n"<< curPara.position()<<std::endl;
  std::cout<<"momentum(): \n"<< curPara.momentum()<<std::endl;
  std::cout<<"covariance: \n"<< *curPara.covariance()<<std::endl;

  std::cout<<"======target parameter info:====="<<std::endl;
  std::cout<<"position(): \n"<< targetParams.position()<<std::endl;
  std::cout<<"momentum(): \n"<< targetParams.momentum()<<std::endl;
  std::cout<<"covariance: \n"<< *targetParams.covariance()<<std::endl;

  return 0;
}
