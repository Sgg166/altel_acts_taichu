void  uu_vv(const std::string& rootFilePath){
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
 
  
  TH1F *hGoodTrajectoriesPerEvent = new TH1F("hGoodTrajectoriesPerEvent", 
                                            "Number of Good Trajectories per Event;Good Trajectories per Event;Number of Events", 
                                            20, -1, 9);
 
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t goodEventNum = 0;
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits=0;
  size_t totalFitHits =0;
 
  size_t eventsWith_One_GoodTrack = 0;        
  size_t eventsWithGoodAndBadTracks = 0;      
  size_t eventsWithMultipleGoodTracks = 0;    
  size_t eventsWithMultipleGoodAndBad = 0;    
  size_t eventsWithNoGoodTracks = 0;          
 
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    bool hasGoodTrack = false;
    size_t goodTracksInEvent = 0;      
    size_t badTracksInEvent = 0;       
 
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==5 );
 
      if(isGoodTrack){
        hasGoodTrack = true;
        goodTracks++;
        goodTracksInEvent++;
 
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
          if(aFitHit && aMatchedMeasHit){
            totalFitHits++;
            totalMatchedHits++;
          }
        }
      } else {
        badTracksInEvent++;  
      }
    }
 
    hGoodTrajectoriesPerEvent->Fill(goodTracksInEvent);
 
    if(hasGoodTrack) {
      goodEventNum++;
      
      if(goodTracksInEvent == 1 && badTracksInEvent == 0) {
        eventsWith_One_GoodTrack++;  //
      } else if(goodTracksInEvent == 1 && badTracksInEvent > 0) {
        eventsWithGoodAndBadTracks++;  
      } else if(goodTracksInEvent > 1 && badTracksInEvent == 0) {
        eventsWithMultipleGoodTracks++;  
      } else if(goodTracksInEvent > 1 && badTracksInEvent > 0) {
        eventsWithMultipleGoodAndBad++;  
      }
    } else {
      eventsWithNoGoodTracks++;  
    }
  }
 
  double efficiency = (goodTracks > 0) ? (double)totalMatchedHits / (double)goodTracks : 0.0;
  std::fprintf(stdout, "========================================\n");
  std::fprintf(stdout, "Event statistics results:\n");
  std::fprintf(stdout, "========================================\n");
  std::fprintf(stdout, "total event number: %zu\n", totalNumEvents);
  std::fprintf(stdout, "The number of events containing good trajectories: %zu\n", goodEventNum);
  std::fprintf(stdout, "The number of events without a good trajectory %zu\n", eventsWithNoGoodTracks);
  std::fprintf(stdout, "total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "good tracks (==5 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "bad tracks: %zu\n", totalTracks - goodTracks);
  std::fprintf(stdout, "----------------------------------------\n");
  std::fprintf(stdout, "Contains a detailed distribution of good trajectory events:\n");
  std::fprintf(stdout, "An event with one good track : %zu\n", eventsWith_One_GoodTrack);
  std::fprintf(stdout, "An event with one good trajectory and multiple bad trajectories: %zu\n", eventsWithGoodAndBadTracks);
  std::fprintf(stdout, "Events with multiple good trajectories: %zu\n", eventsWithMultipleGoodTracks);
  std::fprintf(stdout, "Events with multiple good trajectories and multiple bad trajectories: %zu\n", eventsWithMultipleGoodAndBad);
  std::fprintf(stdout, "----------------------------------------\n");
  std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  std::fprintf(stdout, "Efficiency: %.3f\n",efficiency);
  std::fprintf(stdout, "========================================\n");
 
  // 打印好轨迹数量分布统计
  std::fprintf(stdout, "\nGood track number distribution statistics:\n");
  std::fprintf(stdout, "========================================\n");
  for(int i = 0; i <= 20; i++) {
    int count = hGoodTrajectoriesPerEvent->GetBinContent(i+1); // ROOT bin从1开始
    if(count > 0) {
      std::fprintf(stdout, "The number of events containing %d good trajectories: %d\n", i/2, count);
    }
  }
  std::fprintf(stdout, "========================================\n");
 
  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 2400, 1600);
 
  hGoodTrajectoriesPerEvent->SetFillColor(kBlue);
  hGoodTrajectoriesPerEvent->SetLineColor(kBlack);
  hGoodTrajectoriesPerEvent->SetLineWidth(2);
  hGoodTrajectoriesPerEvent->SetStats(1);  
  hGoodTrajectoriesPerEvent->Draw("HIST TEXT");
 
 /* TLatex statText;
  statText.SetNDC();
  statText.SetTextSize(0.05);
  statText.SetTextColor(kBlack);
  double meanGoodTraj = hGoodTrajectoriesPerEvent->GetMean();
  double rmsGoodTraj = hGoodTrajectoriesPerEvent->GetRMS();
  statText.DrawLatex(0.15, 0.85, Form("Mean = %.2f", meanGoodTraj));
  statText.DrawLatex(0.15, 0.78, Form("RMS = %.2f", rmsGoodTraj));
 */
  c1->SaveAs("residuals_det32.svg");
 
  TFile *outFile = new TFile("residuals_det32.root", "RECREATE");
  hGoodTrajectoriesPerEvent->Write();  
  c1->Write();
  std::fprintf(stdout, "Analysis complete. Results saved to residuals_det32.root\n");
 // std::getc(stdin);
 
  return;
}
