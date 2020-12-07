// #include "TCanvas.h"
// #include "TROOT.h"

  // gROOT->ProcessLine(".include /home/yiliu/testbeam/altel_acts/source/teldata/event/include/");
  // gROOT->ProcessLine(".include /home/yiliu/testbeam/altel_acts/source/teldata/root/include/");
  // gROOT->ProcessLine(".include /home/yiliu/testbeam/altel_acts/source/teldata/root/src/");
  // gROOT->ProcessLine(".L /home/yiliu/testbeam/altel_acts/source/teldata/root/src/TelEventTTreeReader.cpp");

// #include "TelEvent.hpp"
// #include "TelEventTTreeReader.hpp"

#include "TelEventTTreeReader.cpp"

void altelTTreeExampleReader(){

  std::string rootFilePath("/run/media/yiliu/WIN/runspace/calice1908/detresid.root");
  std::cout<< rootFilePath<<std::endl;
  TFile *tfile = new TFile(rootFilePath.c_str(),"READ");
  if(!tfile || !tfile->IsOpen()){
    std::fprintf(stderr, "tfile is not open\n");
    throw;
  }
  TTree *pTree = 0;
  tfile->GetObject("eventTree",pTree);
  if(!pTree){
    std::fprintf(stderr, "pTree is invalid\n");
    throw;
  }
  altel::TelEventTTreeReader ttreeReader;
  ttreeReader.setTTree(pTree);

  size_t totalNumEvents  = ttreeReader.numEvents();
  for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    std::fprintf(stdout, "FileEvent #%zu, event #%u, clock/trigger #%lu\n", eventNum, telEvent->eveN(), telEvent->clkN());
    for(auto traj: telEvent.trajs()){
      traj.
    }

    std::fprintf(stdout, "waiting, press any key to next event\n");
    std::getc(stdin);
  }

  tfile->Close();
  delete tfile;
  return;
}
