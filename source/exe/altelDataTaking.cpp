#include "eudaq/Producer.hh"
#include "Telescope.hh"

#include <list>
#include <iostream>
#include <chrono>
#include <thread>
#include <csignal>

#include "getopt.h"

#include <numeric>
#include <chrono>

#include "TelGL.hh"
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linenoise.h"
#include "myrapidjson.h"

#include "TelFW.hh"
#include "glfw_test.hh"

#include "eudaq/FileReader.hh"
#include "CvtEudaqAltelRaw.hh"

static const std::string default_geometry = R"(
{"geometry": {"detectors": [
    {"id": 0, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 0.0}, "rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 1, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 40.0}, "rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 2, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1}, "center": {"x": 0.0,"y": 0.0,"z":80.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 3, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 120.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 4, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 160.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}},
    {"id": 5, "size": {"x": 29.94176,"y": 13.76256,"z": 1.0}, "pitch": {"x": 0.02924,"y": 0.02688,"z": 1.0},
     "pixel": {"x": 1024,"y": 512,"z": 1},"center": {"x": 0.0,"y": 0.0,"z": 200.0},"rotation": {"x": 0.0,"y": 0.0,"z": 0.0}}
]}}
)";


static const std::string default_rbcpconf =
#include "altel_run_conf.hh"
  ;

static const std::string help_usage = R"(
Usage:
  -help                             help message
  -verbose                          verbose flag
  -geometryFile   <PATH>            path to geometry input file (input)
  -rbcpConfFile   <PATH>            path to eudaq rbcp configure file (input)

examples:
 ./bin/altelDataTaking -geo calice_geo_align4.json
)";


static sig_atomic_t g_done = 0;
int main(int argc, char *argv[]) {
  signal(SIGINT, [](int){g_done+=1;});
  std::string geometryFilePath;
  std::string rbcpConfFilePath;

  int do_wait = 0;
  int do_verbose = 0;
  {////////////getopt begin//////////////////
    struct option longopts[] = {{"help", no_argument, NULL, 'h'},//option -W is reserved by getopt
                                {"verbose", no_argument, NULL, 'v'},//val
                                {"rbcpConfFile", required_argument, NULL, 'e'},
                                {"geometryFile", required_argument, NULL, 'g'},
                                {0, 0, 0, 0}};

    // if(argc == 1){
    //   std::fprintf(stderr, "%s\n", help_usage.c_str());
    //   std::exit(1);
    // }
    int c;
    int longindex;
    opterr = 1;
    while ((c = getopt_long_only(argc, argv, "-", longopts, &longindex)) != -1) {
      // // "-" prevents non-option argv
      // if(!optopt && c!=0 && c!=1 && c!=':' && c!='?'){ //for debug
      //   std::fprintf(stdout, "opt:%s,\targ:%s\n", longopts[longindex].name, optarg);;
      // }
      switch (c) {
      case 'e':
        rbcpConfFilePath = optarg;
        break;
      case 'g':
        geometryFilePath = optarg;
        break;
      case 'w':
        do_wait=1;
        break;
        // help and verbose
      case 'v':
        do_verbose=1;
        //option is set to no_argument
        if(optind < argc && *argv[optind] != '-'){
          do_verbose = std::stoul(argv[optind]);
          optind++;
        }
        break;
      case 'h':
        std::fprintf(stdout, "%s\n", help_usage.c_str());
        std::exit(0);
        break;
        /////generic part below///////////
      case 0:
        // getopt returns 0 for not-NULL flag option, just keep going
        break;
      case 1:
        // If the first character of optstring is '-', then each nonoption
        // argv-element is handled as if it were the argument of an option
        // with character code 1.
        std::fprintf(stderr, "%s: unexpected non-option argument %s\n",
                     argv[0], optarg);
        std::exit(1);
        break;
      case ':':
        // If getopt() encounters an option with a missing argument, then
        // the return value depends on the first character in optstring:
        // if it is ':', then ':' is returned; otherwise '?' is returned.
        std::fprintf(stderr, "%s: missing argument for option %s\n",
                     argv[0], longopts[longindex].name);
        std::exit(1);
        break;
      case '?':
        // Internal error message is set to print when opterr is nonzero (default)
        std::exit(1);
        break;
      default:
        std::fprintf(stderr, "%s: missing getopt branch %c for option %s\n",
                     argv[0], c, longopts[longindex].name);
        std::exit(1);
        break;
      }
    }
  }/////////getopt end////////////////

  std::fprintf(stdout, "\n");
  std::fprintf(stdout, "geometryFile:  <%s>\n", geometryFilePath.c_str());
  std::fprintf(stdout, "rbcpConfFileFile:  <%s>\n", rbcpConfFilePath.c_str());
  std::fprintf(stdout, "\n");

  //////////// geometry

  std::string str_geo;
  if(geometryFilePath.empty()){
    str_geo = default_geometry;
  }
  else{
    str_geo = JsonUtils::readFile(geometryFilePath);
  }
  JsonDocument jsd_geo = JsonUtils::createJsonDocument(str_geo);
  if(jsd_geo.IsNull()){
    std::fprintf(stderr, "Geometry file <%s> does not contain any json objects.\n", geometryFilePath.c_str() );
    throw;
  }

  if (!jsd_geo.HasMember("geometry")) {
    throw;
  }
  TelFW telfw(800, 400, "test");
  glfw_test telfwtest(str_geo, true);
  telfw.startAsync<glfw_test>(&telfwtest, &glfw_test::beginHook, &glfw_test::clearHook, &glfw_test::drawHook);

  std::string str_rbcpconf;
  if(rbcpConfFilePath.empty()){
    str_rbcpconf = default_rbcpconf;
  }
  else{
    str_rbcpconf = JsonUtils::readFile(rbcpConfFilePath);
  }

  std::unique_ptr<altel::Telescope> m_tel;
  m_tel.reset(new altel::Telescope(str_rbcpconf)); // todo
  m_tel->Init();

  m_tel->Start_no_tel_reading();

  while(!g_done){
    auto ev_tel = m_tel->ReadEvent();
    if(ev_tel.empty()){
      std::this_thread::sleep_for(std::chrono::microseconds(100));
      continue;
    }
    std::shared_ptr<eudaq::Event> eudaqEvent = eudaq::Event::MakeUnique("AltelRaw");
    uint64_t trigger_n = ev_tel.front()->GetTrigger();
    eudaqEvent->SetTriggerN(trigger_n);

    std::map<uint32_t, uint32_t> map_layer_clusterN;

    for(auto& e: ev_tel){
      uint32_t word32_count  = 2; // layerID_uint32, cluster_n_uint32
      for(auto &ch : e->m_clusters){
        word32_count += 3; // x_float, y_float , pixel_n_uint32
        word32_count += ch.pixelHits.size(); // pixel_xy_uint32,
      }
      std::vector<uint32_t> layer_block(word32_count);
      uint32_t* p_block = layer_block.data();
      uint32_t layerID = e->GetExtension();
      *p_block =  layerID;
      p_block++;
      uint32_t clusters_size = e->m_clusters.size();
      *p_block = e->m_clusters.size();
      p_block ++;
      for(auto &ch : e->m_clusters){
        *(reinterpret_cast<float*>(p_block)) = ch.x();
        p_block ++;
        *(reinterpret_cast<float*>(p_block)) = ch.y();
        p_block ++;
        *p_block = ch.pixelHits.size();
        p_block ++;
        for(auto &ph : ch.pixelHits){
          // Y<< 16 + X
          *p_block =  uint32_t(ph.x()) + (uint32_t(ph.y())<<16);
          p_block ++;
        }
      }
      if(p_block - layer_block.data() != layer_block.size()){
        std::cerr<<"error altel data block"<<std::endl;
        throw;
      }
      map_layer_clusterN[layerID]= clusters_size;
      eudaqEvent->AddBlock(layerID, layer_block);
    }
    std::shared_ptr<altel::TelEvent> telEvent = altel::createTelEvent(eudaqEvent);
    if(!telEvent){
      continue;
    }
    if(telEvent->measRaws().empty() && telEvent->measHits().empty() && telEvent->trajs().empty()){
      continue;
    }
    if(telEvent->measHits().size()<4){
      continue;
    }
    
    telfwtest.pushBufferEvent(telEvent);
  }
  m_tel->Stop();
  m_tel.reset();
  telfw.stopAsync();
  return 0;
}
