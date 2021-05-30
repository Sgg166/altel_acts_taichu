
#include <regex>
#include "Layer.hh"
#include "FirmwarePortal.hh"

//using namespace std::chrono_literals;
using namespace altel;

Layer::~Layer(){
  m_conn.reset();

  m_is_async_watching = false;
  if(m_fut_async_watch.valid())
    m_fut_async_watch.get();
}

void Layer::start(){
  m_vec_ring_ev.clear();
  m_vec_ring_ev.resize(m_size_ring);
  m_count_ring_write = 0;
  m_count_ring_read = 0;
  m_hot_p_read = m_size_ring -1; // tail

  m_tg_expected = 0;
  m_flag_wait_first_event = true;

  m_st_n_tg_ev_now =0;
  m_st_n_ev_input_now =0;
  m_st_n_ev_output_now =0;
  m_st_n_ev_bad_now =0;
  m_st_n_ev_overflow_now =0;
  m_st_n_tg_ev_begin = 0;

  std::stringstream ssbuf;
  std::string strbuf;
  NetMsg daqMsg{NetMsg::daqstart, 0, 0, 0, {}};
  ssbuf.str(std::string());
  msgpack::pack(ssbuf, daqMsg);
  strbuf = ssbuf.str();
  m_isDataAccept= true;
  m_conn->sendRaw(strbuf.data(), strbuf.size());

  if(!m_is_async_watching){
    m_fut_async_watch = std::async(std::launch::async, &Layer::AsyncWatchDog, this);
  }
}

void Layer::stop(){
  m_isDataAccept= false;

  std::stringstream ssbuf;
  std::string strbuf;
  NetMsg daqMsg{NetMsg::daqstop, 0, 0, 0, {}};
  ssbuf.str(std::string());
  msgpack::pack(ssbuf, daqMsg);
  strbuf = ssbuf.str();
  m_conn->sendRaw(strbuf.data(), strbuf.size());

  m_is_async_watching = false;
  if(m_fut_async_watch.valid())
    m_fut_async_watch.get();
}

void Layer::init(){
  m_isDataAccept= false;
  m_conn =  TcpConnection::connectToServer(m_host,  m_port, reinterpret_cast<FunProcessMessage>(&Layer::perConnProcessRecvMesg), nullptr, this);

  std::stringstream ssbuf;
  std::string strbuf;
  NetMsg daqMsg{NetMsg::daqinit, 0, 0, 0, {}};
  ssbuf.str(std::string());
  msgpack::pack(ssbuf, daqMsg);
  strbuf = ssbuf.str();
  m_conn->sendRaw(strbuf.data(), strbuf.size());

}

int Layer::perConnProcessRecvMesg(void* pconn, msgpack::object_handle& oh){ // IMPROVE IT AS A RING
  if(!m_isDataAccept){
    std::cout<< "msg is dropped"<<std::endl;
    return 0;
  }

  msgpack::object msg = oh.get();
  unique_zone& life = oh.zone();

  NetMsg netmsg;
  try{
    netmsg = std::move(msg.as<NetMsg>());
  }
  catch(...){
    std::cout<< "msg parsing error"<<std::endl;
    return 0;
  }

  std::cout<< "bin "<<TcpConnection::binToHexString(netmsg.bin.data(),netmsg.bin.size())<<std::endl;
  if(netmsg.type!=NetMsg::Type::data){
    std::cout<< "unknown msg type"<<std::endl;
  }

  std::string datastr(netmsg.bin.data(),netmsg.bin.size());
  auto df = std::make_shared<DataFrame>(std::move(datastr));

  m_st_n_ev_input_now ++;
  uint64_t next_p_ring_write = m_count_ring_write % m_size_ring;
  if(next_p_ring_write == m_hot_p_read){
    // buffer full, permanent data lose
    m_st_n_ev_overflow_now ++;
    return 0;
  }
  uint16_t tg_l16 = 0xffff & df->GetCounter();
  //std::cout<< "id "<< tg_l16 <<"  ";
  if(m_flag_wait_first_event){
    m_flag_wait_first_event = false;
    m_extension = df->GetExtension() ;
    m_tg_expected = tg_l16;
    m_st_n_tg_ev_begin = m_tg_expected;
  }
  if(tg_l16 != (m_tg_expected & 0xffff)){
    // std::cout<<(tg_expected & 0x7fff)<< " " << tg_l16<<"\n";
    uint32_t tg_guess_0 = (m_tg_expected & 0xffff0000) + tg_l16;
    uint32_t tg_guess_1 = (m_tg_expected & 0xffff0000) + 0x10000 + tg_l16;
    if(tg_guess_0 > m_tg_expected && tg_guess_0 - m_tg_expected < 200){
      // std::cout<< "missing trigger, expecting : provided "<< (tg_expected & 0xffff) << " : "<< tg_l16<<" ("<< m_extension <<") \n";
      m_tg_expected =tg_guess_0;
    }
    else if (tg_guess_1 > m_tg_expected && tg_guess_1 - m_tg_expected < 200){
      // std::cout<< "missing trigger, expecting : provided "<< (tg_expected & 0xffff) << " : "<< tg_l16<<" ("<< m_extension <<") \n";
      m_tg_expected =tg_guess_1;
    }
    else{
      // std::cout<< "broken trigger ID, expecting : provided "<< (tg_expected & 0xffff) << " : "<< tg_l16<<" ("<<df->GetExtension() <<") \n";
      m_tg_expected ++;
      m_st_n_ev_bad_now ++;
      // permanent data lose
      return 0;
    }
  }
  //TODO: fix tlu firmware, mismatch between modes AIDA start at 1, EUDET start at 0
  df->SetTrigger(m_tg_expected);
  m_st_n_tg_ev_now = m_tg_expected;

  m_vec_ring_ev[next_p_ring_write] = df;
  m_count_ring_write ++;
  m_tg_expected ++;

  return 1;
}

std::string  Layer::GetStatusString(){
  std::unique_lock<std::mutex> lk(m_mtx_st);
  return m_st_string;
}

DataFrameSP& Layer::Front(){
  if(m_count_ring_write > m_count_ring_read) {
    uint64_t next_p_ring_read = m_count_ring_read % m_size_ring;
    m_hot_p_read = next_p_ring_read;
    // keep hot read to prevent write-overlapping
    return m_vec_ring_ev[next_p_ring_read];
  }
  else{
    return m_ring_end;
  }
}

void Layer::PopFront(){
  if(m_count_ring_write > m_count_ring_read) {
    uint64_t next_p_ring_read = m_count_ring_read % m_size_ring;
    m_hot_p_read = next_p_ring_read;
    // keep hot read to prevent write-overlapping
    m_vec_ring_ev[next_p_ring_read].reset();
    m_count_ring_read ++;
  }
}

uint64_t Layer::Size(){
  return  m_count_ring_write - m_count_ring_read;
}

void Layer::ClearBuffer(){
  m_count_ring_write = m_count_ring_read;
  m_vec_ring_ev.clear();
}

uint64_t Layer::AsyncWatchDog(){
  std::chrono::system_clock::time_point m_tp_old;
  std::chrono::system_clock::time_point m_tp_run_begin;

  m_tp_run_begin = std::chrono::system_clock::now();
  m_tp_old = m_tp_run_begin;
  m_is_async_watching = true;

  m_st_n_tg_ev_old =0;
  m_st_n_ev_input_old = 0;
  m_st_n_ev_bad_old =0;
  m_st_n_ev_overflow_old = 0;

  while(m_is_async_watching){
    std::this_thread::sleep_for(std::chrono::seconds(1));
    uint64_t st_n_tg_ev_begin = m_st_n_tg_ev_begin;
    uint64_t st_n_tg_ev_now = m_st_n_tg_ev_now;
    uint64_t st_n_ev_input_now = m_st_n_ev_input_now;
    uint64_t st_n_ev_bad_now = m_st_n_ev_bad_now;
    uint64_t st_n_ev_overflow_now = m_st_n_ev_overflow_now;

    // time
    auto tp_now = std::chrono::system_clock::now();
    std::chrono::duration<double> dur_period_sec = tp_now - m_tp_old;
    std::chrono::duration<double> dur_accu_sec = tp_now - m_tp_run_begin;
    double sec_period = dur_period_sec.count();
    double sec_accu = dur_accu_sec.count();

    // period
    uint64_t st_n_tg_ev_period = st_n_tg_ev_now - m_st_n_tg_ev_old;
    uint64_t st_n_ev_input_period = st_n_ev_input_now - m_st_n_ev_input_old;
    uint64_t st_n_ev_bad_period = st_n_ev_bad_now - m_st_n_ev_bad_old;
    uint64_t st_n_ev_overflow_period = st_n_ev_overflow_now - m_st_n_ev_overflow_old;

    // ratio
    //double st_output_vs_input_accu = st_n_ev_input_now? st_ev_output_now / st_ev_input_now : 1;
    double st_bad_vs_input_accu = st_n_ev_input_now? 1.0 * st_n_ev_bad_now / st_n_ev_input_now : 0;
    double st_overflow_vs_input_accu = st_n_ev_input_now? 1.0 *  st_n_ev_overflow_now / st_n_ev_input_now : 0;
    double st_input_vs_trigger_accu = st_n_ev_input_now? 1.0 * st_n_ev_input_now / (st_n_tg_ev_now - st_n_tg_ev_begin + 1) : 1;
    //double st_output_vs_input_period = st_ev_input_period? st_ev_output_period / st_ev_input_period : 1;
    double st_bad_vs_input_period = st_n_ev_input_period? 1.0 * st_n_ev_bad_period / st_n_ev_input_period : 0;
    double st_overflow_vs_input_period = st_n_ev_input_period? 1.0 *  st_n_ev_overflow_period / st_n_ev_input_period : 0;
    double st_input_vs_trigger_period = st_n_tg_ev_period? 1.0 *  st_n_ev_input_period / st_n_tg_ev_period : 1;

    // hz
    double st_hz_tg_accu = (st_n_tg_ev_now - st_n_tg_ev_begin + 1) / sec_accu ;
    double st_hz_input_accu = st_n_ev_input_now / sec_accu ;

    double st_hz_tg_period = st_n_tg_ev_period / sec_period ;
    double st_hz_input_period = st_n_ev_input_period / sec_period ;

    std::string st_string_new =
      FirmwarePortal::FormatString("L<%u> event(%d)/trigger(%d - %d)=Ev/Tr(%.4f) dEv/dTr(%.4f) tr_accu(%.2f hz) ev_accu(%.2f hz) tr_period(%.2f hz) ev_period(%.2f hz)",
                                   m_extension, st_n_ev_input_now, st_n_tg_ev_now, st_n_tg_ev_begin, st_input_vs_trigger_accu, st_input_vs_trigger_period,
                                   st_hz_tg_accu, st_hz_input_accu, st_hz_tg_period, st_hz_input_period
                                   );

    {
      std::unique_lock<std::mutex> lk(m_mtx_st);
      m_st_string = std::move(st_string_new);
    }

    //write to old
    m_st_n_tg_ev_old = st_n_tg_ev_now;
    m_st_n_ev_input_old = st_n_ev_input_now;
    m_st_n_ev_bad_old = st_n_ev_bad_now;
    m_st_n_ev_overflow_old = st_n_ev_overflow_now;
    m_tp_old = tp_now;
  }
  return 0;
}
