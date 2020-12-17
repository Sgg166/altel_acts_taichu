#include "TelMille.hh"

#include "getopt.h"

#include <iostream>
#include <algorithm>


static const std::string help_usage = R"(
Usage:
  -help                    help message
  -verbose                 verbose flag
  -hitFile          [PATH]   name of data json file
  -pedeSteeringFile [PATH]   path to pede steering file (output)
  -milleBinaryFile  [PATH]   path to mille binary file (output)
  -inputGeometry    [PATH]   geometry input file
  -hitResX          [float]  preset detector hit resolution X
  -hitResY          [float]  preset detector hit resolution Y
  -maxEventNumber   [int]    max number of events to be processed
  -maxTrackNumber   [int]    max number of tracks to be processed

example:

./altelMilleBin  -hitFile /work/data/TB2008_CALICE/jsondata/altel_Run069000.json  -pede pede.txt -mille mille.bin -input ../init_geo.json -maxE 100000 -hitResX 0.1 -hitResY 0.1
./altelMilleBin  -hitFile /work/data/TB2006/alpide_200629033515.json -pede pede.txt -mille mille.bin -input ../init_313_geo.json -maxE 1000000 -hitResX 0.014 -hitResY 0.014
./altelMilleBin  -hitFile /work/data/TB2008_CALICE/test/altel_Run069001.json -pede steerfile -mille mille.bin -input out2_91.json -maxE 100000 -hitResX 0.027 -hitResY 0.027

)";

int main(int argc, char *argv[]) {
  int do_help = false;
  int do_verbose = false;
  struct option longopts[] = {{"help", no_argument, &do_help, 1},
                              {"verbose", no_argument, &do_verbose, 1},
                              {"hitFile", required_argument, NULL, 'f'},
                              {"inputGeomerty", required_argument, NULL, 'g'},
                              {"pedeSteeringFile", required_argument, NULL, 'u'},
                              {"milleBinaryFile", required_argument, NULL, 'q'},
                              {"hitResX", required_argument, NULL, 'r'},
                              {"hitResY", required_argument, NULL, 's'},
                              {"maxEventNumber", required_argument, NULL, 'm'},
                              {"maxTrackNumber", required_argument, NULL, 'n'},
                             {0, 0, 0, 0}};

  std::string hitFile_path;
  std::string outputGeometry_path;
  std::string inputGeomerty_path;
  std::string pedeSteeringFile_path;
  std::string milleBinaryFile_path;
  size_t maxTrackNumber = -1;
  size_t maxEventNumber = -1;

  double hitResX =0.02;
  double hitResY =0.02;

  int c;
  opterr = 1;
  while ((c = getopt_long_only(argc, argv, "", longopts, NULL)) != -1) {
    switch (c) {
    case 'f':
      hitFile_path = optarg;
      break;
    case 'o':
      outputGeometry_path = optarg;
      break;
    case 'g':
      inputGeomerty_path = optarg;
      break;
    case 'u':
      pedeSteeringFile_path = optarg;
      break;
    case 'q':
      milleBinaryFile_path = optarg;
      break;
    case 'r':
      hitResX = std::stod(optarg);
      break;
    case 's':
      hitResY = std::stod(optarg);
      break;
    case 'm':
      maxEventNumber = std::stoull(optarg);
      break;
    case 'n':
      maxTrackNumber = std::stoull(optarg);
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

  if (hitFile_path.empty() ||
      inputGeomerty_path.empty()||
      milleBinaryFile_path.empty()||
      pedeSteeringFile_path.empty()
    ) {
    std::fprintf(stderr, "%s\n", help_usage.c_str());
    std::exit(0);
  }

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "hitFile:            %s\n", hitFile_path.c_str());
  std::fprintf(stdout, "inputGeomerty:      %s\n", inputGeomerty_path.c_str());
  std::fprintf(stdout, "milleBinaryFile:    %s\n", milleBinaryFile_path.c_str());
  std::fprintf(stdout, "pedeSteeringFile:   %s\n", pedeSteeringFile_path.c_str());
  std::fprintf(stdout, "hitResX:            %f\n", hitResX);
  std::fprintf(stdout, "hitResY:            %f\n", hitResY);
  std::fprintf(stdout, "\n");

  std::printf("--------read geo-----\n");
  std::string str_geo = JsonUtils::readFile(inputGeomerty_path.c_str());
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", inputGeomerty_path.c_str() );
    throw;
  }
  // JsonUtils::printJsonValue(jsd_geo, true);

  altel::TelMille telmille;
  telmille.setGeometry(jsd_geo);
  telmille.setResolution(hitResX, hitResY);
  telmille.startMilleBinary(milleBinaryFile_path);

  size_t nGeoLayers = jsd_geo["geometry"]["detectors"].Size();

  size_t n_datapack_select_opt = 20000;
  JsonFileDeserializer jsf(hitFile_path);
  JsonAllocator jsa;

  size_t nTracks = 0;
  size_t nEvents = 0;
  while (jsf) {
    if (nTracks >= maxTrackNumber || nEvents >= maxEventNumber ) {
      break;
    }
    auto evpack = jsf.getNextJsonDocument();
    if(evpack.IsNull()){
      std::fprintf(stdout, "reach null object, end of file\n");
      break;
    }
    nEvents++;

    // std::cout<< ">>>>"<<std::endl;

    JsonValue js_track_filtered(rapidjson::kArrayType);
    const auto &layers = evpack["layers"];
    for (const auto &layer : layers.GetArray()) {
      size_t id_ext = layer["ext"].GetUint();
      for (const auto &hit : layer["hit"].GetArray()) {
        double x_hit = hit["pos"][0].GetDouble() - 0.02924 * 1024 / 2.0;
        double y_hit = hit["pos"][1].GetDouble() - 0.02688 * 512 / 2.0;
        JsonValue js_hit(rapidjson::kObjectType);
        js_hit.AddMember("id", id_ext, jsa);
        js_hit.AddMember("x", x_hit, jsa);
        js_hit.AddMember("y", y_hit, jsa);
        js_track_filtered.PushBack(std::move(js_hit), jsa);
      }
    }
    // JsonUtils::printJsonValue(js_track_filtered, false);

    // drop multiple hits events
    if (js_track_filtered.Size() != nGeoLayers) { // TODO
      // std::cout<< "skipping event "<<nEvents <<std::endl;
      continue;
    }
    std::vector<size_t> geo_ids;
    bool found_same_geo_id = false;
    for (const auto &js_hit : js_track_filtered.GetArray()) {
      size_t geo_id = js_hit["id"].GetUint();
      if (std::find(geo_ids.begin(), geo_ids.end(), geo_id) !=
          geo_ids.end()) {
        found_same_geo_id = true;
        break;
      }
      geo_ids.push_back(geo_id);
    }
    if (found_same_geo_id) {
      // std::cout<< "skipping muilt-tracks event "<<nEvents <<std::endl;
      // JsonUtils::printJsonValue(js_track_filtered, false);
      continue;
    }

    nTracks++;
    telmille.fillTrackXYRz(js_track_filtered);
  }
  std::fprintf(stdout, "%i tracks are picked from %i events\n", nTracks, nEvents);

  telmille.endMilleBinary();
  telmille.createPedeStreeringModeXYRz(pedeSteeringFile_path);

  // JsonUtils::printJsonValue(jsd_geo, true);
  // std::string jsstr = JsonUtils::stringJsonValue(jsd_geo, true);
  // std::FILE *fp = std::fopen(outputfile_name.c_str(), "w");
  // std::fwrite(jsstr.data(), 1, jsstr.size(), fp);
  // std::fclose(fp);
  return 0;
}
