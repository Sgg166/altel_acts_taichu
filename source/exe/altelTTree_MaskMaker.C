
// 1. create root init file:   rootlogon.C   (absolute path is preferred)
// {
//   gROOT->ProcessLine(".include /work/testbeam/altel_acts/source/teldata/event/include/");
//   gROOT->ProcessLine(".include /work/testbeam/altel_acts/source/teldata/root/include/");
//   gROOT->ProcessLine(".L /work/testbeam/altel_acts/source/teldata/root/src/TelEventTTreeReader.cpp");
// }
// 2. run root macro in the same folder of rootlogon.C
//       >> root -l 'altelTTreeExampleReader.C("../geo7/geo7_gev4p4_200629005112_t3.root")'

void altelTTree_MaskMaker(const std::string& rootFilePath){
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

  std::map<std::array<uint16_t, 3>, uint64_t> fire_count;

  size_t totalNumEvents  = ttreeReader.numEvents();
  for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    std::fprintf(stdout, "Event:  FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu, numMeasRaw %zu \n",
                 eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size(), telEvent->measRaws().size());

    size_t nmr = 0;
    for(auto &mh: telEvent->measHits()){
        nmr += mh->measRaws().size();
    }
    if(nmr>2){
      //skip sparking event
      continue;
    }
    for(auto &mh: telEvent->measHits()){
      for(auto &mr : mh->measRaws()){
        uint16_t x = mr.u();
        uint16_t y = mr.v();
        uint16_t id = mr.detN();
        if(fire_count.count({id, x, y})==0){
          fire_count[{id, x, y}]= 0;
        }
        fire_count[{id, x, y}]+=1;
      }
    }
    std::fprintf(stdout, "numMeasRaw(in MeasHits) %zu \n", nmr);

    // std::fprintf(stdout, "waiting, press any key to next event\n");
    // std::getc(stdin);
  }

  std::cout<<"================="<<fire_count.size()<<std::endl;
  for(auto &[idxy, count]: fire_count){
    //totalNumEvents/20000;
    if(count>totalNumEvents/50000.){
      printf("%hu [%hu, %hu]  ,  %f \n", idxy[0], idxy[1], idxy[2], count*1.0/totalNumEvents);
    }
  }

  tfile->Close();
  delete tfile;
  return;
}
