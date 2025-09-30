void Efficiency(const std::string& rootFilePath){
  std::cout << rootFilePath << std::endl;
  TFile *tfile = new TFile(rootFilePath.c_str(), "READ");
  if(!tfile || !tfile->IsOpen()){
    std::fprintf(stderr, "tfile is not open\n");
    throw;
  }
  TTree *pTree = 0;
  tfile->GetObject("eventTree", pTree);
  if(!pTree){
      std::fprintf(stderr, "pTree is invalid\n");
      throw;
  }
  altel::TelEventTTreeReader ttreeReader;
  ttreeReader.setTTree(pTree);
 
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits = 0;
  size_t totalFitHits = 0;
  size_t Invalid_tracks= 0;
  size_t hit_dut_goodTracks=0;
  size_t eventsWithGoodTracks = 0;
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
    bool eventHasGoodTrack = false;
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);
      if(!isGoodTrack){
        continue;
      }
     // eventHasGoodTrack = true;
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        auto aMatchedMeasHit = aTrajHit->matchedMeasHit();
        if(!aFitHit){
          continue;
        }
        if (aFitHit->detN() != 32){
          continue;
        }
        goodTracks++;
        totalFitHits++;
        double hit_u = aFitHit->u();
        double hit_v = aFitHit->v();
        if(hit_u < -12.8 || hit_u > 12.8 || hit_v < -6.4 || hit_v > 6.4 ){
          Invalid_tracks++;
          continue;
         //  std::fprintf(stdout,"(hit_u,hit_v) :(%f,%f)\n", hit_u,hit_v);
        }
        hit_dut_goodTracks++;
        eventHasGoodTrack = true;
        if(aMatchedMeasHit){
          totalMatchedHits++; 
        }
      }
    }
    if(eventHasGoodTrack){
      eventsWithGoodTracks++;
    }
  }
  double efficiency = (goodTracks > 0) ? (double)totalMatchedHits / (double)hit_dut_goodTracks : 0.0;
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Events with good tracks: %zu\n", eventsWithGoodTracks);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good tracks (==5 hits ): %zu\n", goodTracks);
  std::fprintf(stdout, "Good tracks (==5 hits) and hit dut: %zu\n", hit_dut_goodTracks);
  std::fprintf(stdout, "Good tracks but no hit dut: %zu\n", Invalid_tracks);
  std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  std::fprintf(stdout, "Efficiency: %.3f\n",efficiency);
 
  std::fprintf(stdout, "Analysis complete. Results saved to residuals_det32.root\n");
  return;
}
