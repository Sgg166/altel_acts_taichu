void analyzeOriginMeasHits(const std::string& rootFilePath){
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
 
  // 创建直方图
  TH1D *hNumOriginHits = new TH1D("hNumOriginHits", "Number of Original Measured Hits per Fitted Hit", 20, 0, 20);
  hNumOriginHits->GetXaxis()->SetTitle("Number of Original Measured Hits");
  hNumOriginHits->GetYaxis()->SetTitle("Entries");
 
  TH1D *hNumOriginHitsGoodTracks = new TH1D("hNumOriginHitsGoodTracks", "Number of Original Measured Hits per Fitted Hit (Good Tracks Only)", 20, 0, 20);
  hNumOriginHitsGoodTracks->GetXaxis()->SetTitle("Number of Original Measured Hits");
  hNumOriginHitsGoodTracks->GetYaxis()->SetTitle("Entries");
 
  size_t totalNumEvents  = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
 
  for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() >= 3); // 好轨迹标准
      
      if(isGoodTrack){
        goodTracks++;
      }
 
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        if(aFitHit){
          auto aOriginMeasHit = aFitHit->originMeasHit();
          if(aOriginMeasHit){
            size_t numOrigin = aOriginMeasHit->measRaws().size();
            hNumOriginHits->Fill(numOrigin); // 所有轨迹的统计
            
            if(isGoodTrack){
              hNumOriginHitsGoodTracks->Fill(numOrigin); // 只有好轨迹的统计
            }
          }
        }
      }
    }
  }
 
  // 输出统计信息
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good tracks (>=6 hits): %zu\n", goodTracks);
 
  // 创建画布并绘制直方图
  TCanvas *c1 = new TCanvas("c1", "Origin Measured Hits Analysis", 1200, 600);
  c1->Divide(2,1);
 
  c1->cd(1);
  hNumOriginHits->Draw();
  hNumOriginHits->SetFillColor(kBlue);
  hNumOriginHits->SetFillStyle(3001);
 
  c1->cd(2);
  hNumOriginHitsGoodTracks->Draw();
  hNumOriginHitsGoodTracks->SetFillColor(kRed);
  hNumOriginHitsGoodTracks->SetFillStyle(3001);
 
  c1->Update();
 
  // 保存结果
  TFile *outFile = new TFile("origin_hits_analysis.root", "RECREATE");
  hNumOriginHits->Write();
  hNumOriginHitsGoodTracks->Write();
  outFile->Close();
 
  tfile->Close();
  delete tfile;
 
  std::fprintf(stdout, "Analysis complete. Results saved to origin_hits_analysis.root\n");
  std::fprintf(stdout, "Press any key to exit...\n");
  std::getc(stdin);
 
  return;
}
