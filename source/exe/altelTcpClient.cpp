#include "Telescope.hh"

#include <list>
#include <iostream>
#include <chrono>
#include <thread>
#include <csignal>

#include <TFile.h>
#include <TTree.h>

#include "getopt.h"

#include <numeric>
#include <chrono>

#include <TelEvent.hpp>
#include <TelEventTTreeWriter.hpp>

static std::vector<std::unique_ptr<altel::Layer>> g_layers;
static sig_atomic_t g_done = 0;

int watchDog(){
  while(!g_done){
    std::this_thread::sleep_for(std::chrono::seconds(1));
    for(auto &l: g_layers){
      std::string l_status = l->GetStatusString();
      std::fprintf(stdout, "%s\n", l_status.c_str());
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {
  signal(SIGINT, [](int){g_done+=1;});

  g_layers.push_back(std::make_unique<altel::Layer>("test0", "131.169.133.170", 9000));
  g_layers.push_back(std::make_unique<altel::Layer>("test1", "131.169.133.171", 9000));

  for(auto &l: g_layers){
    l->init();
  }

  for(auto &l: g_layers){
    l->start();
  }

  altel::TelEventTTreeWriter ttreeWriter;
  TTree *pTree = new TTree("eventTree", "eventTree");
  ttreeWriter.setTTree(pTree);

  std::future<int> fut_watchdog = std::async(std::launch::async, &watchDog);
  std::cout<<"main loop"<<std::endl;

  uint32_t eventN= 0;
  while(!g_done){
    uint32_t trigger_n = -1;
    for(auto &l: g_layers){
      if( l->Size() == 0){
        // TODO check cached size of all layers
        continue;
      }
      else{
        uint32_t trigger_n_ev = l->Front()->GetTrigger();
        if(trigger_n_ev< trigger_n)
          trigger_n = trigger_n_ev;
      }
    }

    std::vector<altel::DataFrameSP> ev_sync;
    for(auto &l: g_layers){
      auto &ev_front = l->Front();
      // if(ev_front) l->PopFront();
      if(ev_front && ev_front->GetTrigger() == trigger_n){
        ev_sync.push_back(ev_front);
        l->PopFront();
      }
    }

    auto telev_sync = std::make_shared<altel::TelEvent>(0, eventN, 0, trigger_n);
    for(auto &ev: ev_sync){
      if(ev->m_raw.size()!=16){
        std::cout<< "bin "<<TcpConnection::binToHexString(ev->m_raw.data(),ev->m_raw.size())<<std::endl;
        ev->Print(std::cout, 0);
        std::cout<<std::endl;
      }
      auto telev = altel::Layer::createTelEvent(ev->m_raw);
      telev_sync->MRs.insert(telev_sync->MRs.end(), telev->MRs.begin(),telev->MRs.end());
      telev_sync->MHs.insert(telev_sync->MHs.end(), telev->MHs.begin(),telev->MHs.end());
    }

    ttreeWriter.fillTelEvent(telev_sync);
    eventN ++;
  }

  if(fut_watchdog.valid())
    fut_watchdog.get();

  for(auto &l: g_layers){
    l->stop();
  }

  TFile tfile("testTcpClient.root","recreate");
  pTree->Write();
  tfile.Close();

  return 0;
}

