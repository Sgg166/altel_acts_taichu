#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Units.hpp"
#include "PixelMultiTrajectory.hpp"
#include "PixelSourceLink.hpp"

namespace Telescope {

class TelescopeSeeddingAlgorithm final : public FW::BareAlgorithm {
public:

  struct Config {
    /// Input data.
    std::string inputSourcelinks{"sourcelinks_to_seed"};
    /// Output fitted trajectories collection.
    std::string outputSeed{"seeds"};

    double seedResX{10 * Acts::UnitConstants::um};
    double seedResY{10 * Acts::UnitConstants::um};
    double seedResPhi{0.7};
    double seedResTheta{0.7};
    double seedBeamEnergy{5 * Acts::UnitConstants::GeV};
  };

  /// Constructor of the seeding algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  TelescopeSeedingAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the seeding algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  FW::ProcessCode execute(const FW::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace FW
