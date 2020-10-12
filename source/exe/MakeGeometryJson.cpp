#include "JsonGenerator.hpp"
#include "getopt.h"

#include <cstdio>
#include <iostream>

static const std::string help_usage = R"(
Usage:
  -help              help message
  -verbose           verbose flag
  -dataHit  [PATH]   path to detector hit data (input json file)
  -geometry [PATH]   path to geometry file (output json file)
)";

int main(int argc, char *argv[]) {
  rapidjson::CrtAllocator jsa;

  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"data", required_argument, NULL, 'f'},
                              {"geometry", required_argument, NULL, 'g'},
                              {0, 0, 0, 0}};

  std::string dataHit_path;
  std::string geometry_path;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'h':
      do_help = 1;
      std::cout << "here" << std::endl;
      std::fprintf(stdout, "%s\n", help_usage.c_str());
      exit(0);
      break;
    case 'f':
      dataHit_path = optarg;
      break;
    case 'g':
      geometry_path = optarg;
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

  if (dataHit_path.empty()) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    exit(0);
  }

  std::map<size_t, std::array<double, 6>> geoconf;
  geoconf = Telescope::ReadGeoFromDataFile(dataHit_path);

  for (const auto &[id, lgeo] : geoconf) {
    std::printf("layer: %lu   centerX: %f   centerY: %f   centerZ: %f  "
                "rotationX: %f   rotationY: %f   rotationZ: %f\n",
                id, lgeo[0], lgeo[1], lgeo[2], lgeo[3], lgeo[4], lgeo[5]);
  }

  double beamEnergy =
      Telescope::ReadBeamEnergyFromDataFile(datafile_name) *
      Acts::UnitConstants::GeV;

  std::fprintf(stdout, "beamEnergy:       %f\n", beamEnergy);

  return 0;
}
