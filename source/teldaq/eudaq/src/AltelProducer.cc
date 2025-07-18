#include "eudaq/Producer.hh"

#include <list>
#include <iostream>
#include <chrono>
#include <thread>
#include <regex>

#include "Telescope.hh"

template<typename ... Args>
static std::string FormatString( const std::string& format, Args ... args ){
  std::size_t size = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1;
  std::unique_ptr<char[]> buf( new char[ size ] );
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 );
}


namespace altel{
  class AltelProducer : public eudaq::Producer {
  public:
    using eudaq::Producer::Producer;
    ~AltelProducer() override {};
    void DoInitialise() override;
    void DoConfigure() override;
    void DoStartRun() override;
    void DoStopRun() override;
    void DoTerminate() override;
    void DoReset() override;
    void DoStatus() override;

    static const uint32_t m_id_factory = eudaq::cstr2hash("AltelProducer");

    void RunLoop() override;
  private:
    bool m_exit_of_run;
    std::unique_ptr<altel::Telescope> m_tel;

    std::atomic<uint64_t> m_tg_n_begin;
    std::atomic<uint64_t> m_tg_n_last;

    std::chrono::system_clock::time_point m_tp_run_begin;

    uint64_t m_st_n_tg_old;
    std::chrono::system_clock::time_point m_st_tp_old;
  };
}

namespace{
  auto _reg_ = eudaq::Factory<eudaq::Producer>::
    Register<altel::AltelProducer, const std::string&, const std::string&>(altel::AltelProducer::m_id_factory);
}


namespace{
  std::string LoadFileToString(const std::string& path){
    std::ifstream ifs(path);
    if(!ifs.good()){
        std::cerr<<"LoadFileToString:: ERROR, unable to load file<"<<path<<">\n";
        throw;
    }

    std::string str;
    str.assign((std::istreambuf_iterator<char>(ifs) ),
               (std::istreambuf_iterator<char>()));
    return str;
  }
}

void altel::AltelProducer::DoInitialise(){
  m_tel.reset();
  const eudaq::Configuration &param = *GetInitConfiguration();
  param.Print();

  std::vector<std::string> vecLayerName;
  std::string tel_json_str;

  if(param.Has("GEOMETRY_SETUP")){
    std::map<std::string, double> mapLayerPos;
    std::string str_GEOMETRY_SETUP;
    str_GEOMETRY_SETUP = param.Get("GEOMETRY_SETUP", "");
    std::regex block_regex("([a-zA-Z0-9]+)\\:([0-9]+)"); // sm[1]  name, sm[2]  pos
    auto blocks_begin = std::sregex_iterator(str_GEOMETRY_SETUP.begin(), str_GEOMETRY_SETUP.end(), block_regex);
    auto blocks_end = std::sregex_iterator();
    std::cout << "Ini file: found " << std::distance(blocks_begin, blocks_end) << " telescope layers"<<std::endl;
    for (std::sregex_iterator ism = blocks_begin; ism != blocks_end; ++ism){
      // std::smatch &sm= *ism;
      std::string sm_str = (*ism).str();
      std::string layer_name = (*ism)[1].str();
      vecLayerName.push_back(layer_name);
      double layer_pos = std::stod((*ism)[2].str());
      mapLayerPos[layer_name] = layer_pos;
    }

    rapidjson::Document jsdoc;
    rapidjson::Document::AllocatorType& a = jsdoc.GetAllocator();
    jsdoc.SetObject();
    jsdoc.AddMember("telescope", rapidjson::Value(rapidjson::kObjectType), a);
    jsdoc["telescope"].AddMember("locations", rapidjson::Value(rapidjson::kObjectType), a);
    jsdoc["telescope"].AddMember("config", rapidjson::Value(rapidjson::kObjectType), a);
    for(const auto pairLayerPos : mapLayerPos){
      rapidjson::Value name_js(pairLayerPos.first.c_str(), a);
      rapidjson::Value pos_js(pairLayerPos.second);
      jsdoc["telescope"]["locations"].AddMember(name_js, pos_js, a);
    }
    rapidjson::StringBuffer sb;
    // rapidjson::Writer<rapidjson::StringBuffer> writer(sb);
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
    jsdoc.Accept(writer);
    tel_json_str = sb.GetString();
  }

  std::map<std::string,  std::set<std::pair<uint16_t, uint16_t>>> mask_col;
  for(const auto & lname : vecLayerName){
    std::string pmask_para_key("PIXEL_MASK_OVERRIDE_");
    pmask_para_key+=lname;
    if(param.Has(pmask_para_key)){
      std::set<std::pair<uint16_t, uint16_t>> maskXYs;
      std::string str_PIXEL_MASK_OVERRIDE_x;
      str_PIXEL_MASK_OVERRIDE_x = param.Get(pmask_para_key, "");
      {
        std::regex block_regex("([0-9]+)\\:([0-9]+)"); // X:Y
        auto blocks_begin = std::sregex_iterator(str_PIXEL_MASK_OVERRIDE_x.begin(), str_PIXEL_MASK_OVERRIDE_x.end(), block_regex);
        auto blocks_end = std::sregex_iterator();
        std::cout << "found layer <"<<lname<<"> mask size: "  << std::distance(blocks_begin, blocks_end) << " "<<std::endl;
        for (std::sregex_iterator ism = blocks_begin; ism != blocks_end; ++ism){
          const std::smatch &sm= *ism;
          uint16_t maskx = (uint16_t)std::stoul(sm[1].str());
          uint16_t masky = (uint16_t)std::stoul(sm[2].str());
          maskXYs.emplace(maskx, masky);
        }
      }
      {
        std::regex block_regex("([0-9]+)\\:([NnYyXx])"); // X:Y
        auto blocks_begin = std::sregex_iterator(str_PIXEL_MASK_OVERRIDE_x.begin(), str_PIXEL_MASK_OVERRIDE_x.end(), block_regex);
        auto blocks_end = std::sregex_iterator();
        std::cout << "found layer <"<<lname<<"> mask size: "  << std::distance(blocks_begin, blocks_end) << " "<<std::endl;
        for (std::sregex_iterator ism = blocks_begin; ism != blocks_end; ++ism){
          const std::smatch &sm= *ism;
          uint16_t maskx = (uint16_t)std::stoul(sm[1].str());
          for(uint16_t masky = 0; masky<512; masky++ ){
            maskXYs.emplace(maskx, masky);
          }
        }
      }

      {
        std::regex block_regex("([NnYyXx])\\:([0-9]+)"); // X:Y
        auto blocks_begin = std::sregex_iterator(str_PIXEL_MASK_OVERRIDE_x.begin(), str_PIXEL_MASK_OVERRIDE_x.end(), block_regex);
        auto blocks_end = std::sregex_iterator();
        std::cout << "found layer <"<<lname<<"> mask size: "  << std::distance(blocks_begin, blocks_end) << " "<<std::endl;
        for (std::sregex_iterator ism = blocks_begin; ism != blocks_end; ++ism){
          const std::smatch &sm= *ism;
          uint16_t masky = (uint16_t)std::stoul(sm[2].str());
          for(uint16_t maskx = 0; maskx<1024; maskx++ ){
            maskXYs.emplace(maskx, masky);
          }
        }
      }
      mask_col.emplace(lname, std::move(maskXYs));
    }
  }

  if(!tel_json_str.empty()){
    m_tel.reset(new Telescope(tel_json_str,""));
  }else{
    std::cout<<"not able to create tele_json from eudaq init file"<<std::endl;
    m_tel.reset(new Telescope("", ""));
  }
  if(m_tel)  m_tel->Init();
  if(m_tel && !mask_col.empty() )  m_tel->FlushPixelMask(mask_col);

}

void altel::AltelProducer::DoConfigure(){
  //do nothing here


}

void altel::AltelProducer::DoStartRun(){
  m_tel->Start_no_tel_reading();
}

void altel::AltelProducer::DoStopRun(){
  m_tel->Stop();
  m_exit_of_run = true;
}

void altel::AltelProducer::DoReset(){
  m_tel.reset();
}

void altel::AltelProducer::DoTerminate(){
  std::terminate();
}

void altel::AltelProducer::RunLoop(){

  m_tp_run_begin = std::chrono::system_clock::now();
  auto tp_start_run = std::chrono::steady_clock::now();
  m_exit_of_run = false;
  bool is_first_event = true;
  while(!m_exit_of_run){
    auto telev = m_tel->ReadEvent();
    if(!telev){
      std::this_thread::sleep_for(std::chrono::microseconds(100));
      continue;
    }

    uint64_t trigger_n = telev->clkN();

    auto ev_eudaq = eudaq::Event::MakeUnique("AltelRaw");
    ev_eudaq->SetTriggerN(trigger_n);

    std::map<uint32_t,  std::vector<std::shared_ptr<altel::TelMeasHit>>> map_layer_measHits;
    std::vector<uint32_t> detNs={9,7,5,3,2,32};
    for(auto& detN: detNs){
      map_layer_measHits[detN];
    }
    for(auto& mh: telev->measHits()){
      if(!mh){
        continue;
      }
      uint32_t detN = mh->detN();
      map_layer_measHits[detN].push_back(mh);
    }

    for(auto& [detN, mhs]: map_layer_measHits){
      uint32_t word32_count  = 2; // layerID_uint32, cluster_n_uint32
      for(auto& mh : mhs){
        word32_count += 3; // x_float, y_float , pixel_n_uint32
        word32_count += mh->measRaws().size(); // pixel_xy_uint32
      }

      std::vector<uint32_t>  layer_block(word32_count);
      uint32_t* p_block = layer_block.data();
      *p_block =  detN;

      p_block++;
      *p_block = mhs.size(); // cluster_n

      for(auto &mh : mhs){
        p_block ++;
        *(reinterpret_cast<float*>(p_block)) = mh->u();

        p_block ++;
        *(reinterpret_cast<float*>(p_block)) = mh->v();

        p_block ++;
        *p_block = mh->measRaws().size();

        for(auto &mr : mh->measRaws()){
          // Y<< 16 + X
          p_block ++;
          *p_block =  uint32_t(mr.u()) + (uint32_t(mr.v())<<16);
        }
      }
      ev_eudaq->AddBlock(detN, layer_block);
    }

    SendEvent(std::move(ev_eudaq));
    if(is_first_event){
      is_first_event = false;
      m_tg_n_begin = trigger_n;
    }
    m_tg_n_last = trigger_n;
  }
}

void altel::AltelProducer::DoStatus() {
  if(m_exit_of_run){
    auto tp_now = std::chrono::system_clock::now();
    m_st_tp_old = tp_now;
    m_st_n_tg_old = 0;
  }

  if (!m_exit_of_run) {
    auto tp_now = std::chrono::system_clock::now();
    std::chrono::duration<double> dur_period_sec = tp_now - m_st_tp_old;
    std::chrono::duration<double> dur_accu_sec = tp_now - m_tp_run_begin;
    uint64_t st_n_tg_now = m_tg_n_last;
    uint64_t st_n_tg_begin = m_tg_n_begin;

    double sec_accu = dur_accu_sec.count();
    double sec_period = dur_period_sec.count();

    double st_hz_tg_accu = (st_n_tg_now - st_n_tg_begin) / sec_accu ;
    double st_hz_tg_period = (st_n_tg_now - m_st_n_tg_old) / sec_period ;

    SetStatusTag("TriggerID(latest:first)", FormatString("%u:%u", st_n_tg_now, st_n_tg_begin));
    SetStatusTag("TriggerHz(per:avg)", FormatString("%.1f:%.1f", st_hz_tg_period, st_hz_tg_accu));
    // SetStatusTag("EventHz(per,avg)", std::to_string());
    // SetStatusTag("Cluster([layer:avg:per])", std::to_string());

    m_st_tp_old = tp_now;
    m_st_n_tg_old = st_n_tg_now;
  }
}
