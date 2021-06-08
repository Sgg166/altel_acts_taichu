#include "eudaq/Producer.hh"

#include <list>
#include <iostream>
#include <chrono>
#include <thread>

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
  std::cout<< "it is AltelProducer "<<std::endl;
  auto ini = GetInitConfiguration();
  m_tel.reset(new altel::Telescope("builtin", "builtin"));
  m_tel->Init();
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
