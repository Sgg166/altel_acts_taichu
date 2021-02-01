
// 1. create root init file:   rootlogon.C   (absolute path is preferred)
// {
//   gROOT->ProcessLine(".include /work/testbeam/altel_acts/source/teldata/event/include/");
//   gROOT->ProcessLine(".include /work/testbeam/altel_acts/source/teldata/root/include/");
//   gROOT->ProcessLine(".L /work/testbeam/altel_acts/source/teldata/root/src/TelEventTTreeReader.cpp");
// }
// 2. run root macro in the same folder of rootlogon.C
//       >> root -l 'altelTTreeExampleReader.C("../geo7/geo7_gev4p4_200629005112_t3.root")'

void altelTTreeExampleReader(const std::string& rootFilePath){
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
    std::fprintf(stdout, "Event:  FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
                 eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
    for(auto aTraj: telEvent->trajs()){
      std::fprintf(stdout, "  Traj: numOriginMeasHit %zu, numFitHit %zu\n", aTraj->numOriginMeasHit(), aTraj->numFitHit());
      if(aTraj->numOriginMeasHit()<6){
        // not so good, skip
        continue;
      }
      for(auto &aTrajHit: aTraj->trajHits()){
        std::fprintf(stdout, "    Hit: ");
        auto aFitHit = aTrajHit->fitHit();
        if(aFitHit){
          std::fprintf(stdout, " fit-xyzuv(%f, %f, %f,     %f, %f)",
                       aFitHit->x(), aFitHit->y(), aFitHit->z(), aFitHit->u(), aFitHit->v());
          // uv are local, xyz are global
          auto aOriginMeasHit = aFitHit->originMeasHit();
          if(aOriginMeasHit){
            std::fprintf(stdout, " origin-uv(%f, %f)", aOriginMeasHit->u(), aOriginMeasHit->v());
          }
        }
        std::fprintf(stdout, "\n");
        // please see TelEvent.hh to get  more detailed traj values
      }
    }
    std::fprintf(stdout, "waiting, press any key to next event\n");
    std::getc(stdin);
  }

  tfile->Close();
  delete tfile;
  return;
}
