void after_Hit_U_V_TH2F(const std::string& rootFilePath){
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
 
  TH2F* huv_cluster_det2  = new TH2F("huv_cluster_det2" , "Origin_Cluster UV Distribution det2"  , 100, -12.8, 12.8, 100, -6.4, 6.4);
  TH2F* huv_cluster_det3  = new TH2F("huv_cluster_det3" , "Origin_Cluster UV Distribution det3"  , 100, -12.79, 12.81, 100, -6.61, 6.19);
  TH2F* huv_cluster_det32 = new TH2F("huv_cluster_det32", "Matched_Cluster UV Distribution det32", 100, -11.6, 14, 100, -5.5, 7.3);
  TH2F* huv_cluster_det5  = new TH2F("huv_cluster_det5" , "Origin_Cluster UV Distribution det5"  , 100, -12.82, 12.78, 100, -7.5, 5.3);
  TH2F* huv_cluster_det7  = new TH2F("huv_cluster_det7" , "Origin_Cluster UV Distribution det7"  , 100, -12.97, 12.68, 100, -6, 6.8);
  TH2F* huv_cluster_det9  = new TH2F("huv_cluster_det9" , "Origin_Cluster UV Distribution det9"  , 100, -12.8, 12.8, 100, -6.4, 6.4);
 
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalFitHits =0;
  size_t totalMatchedHits=0;
  // std::fprintf(stdout,"11111111111111totalNumEvents: %zu",totalNumEvents);
  for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);

      if(isGoodTrack){
        goodTracks++;
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
          int detID = aTrajHit->detN();
          if(aFitHit){
            totalFitHits++;
            auto aOriginMeasHit = aFitHit->originMeasHit();
            if(aOriginMeasHit){
              double u=aOriginMeasHit->u();
              double v=aOriginMeasHit->v();
               if(detID == 2){
                 huv_cluster_det2->Fill(u, v);
               }
               else if(detID == 3) {
                 huv_cluster_det3->Fill(u, v);
               }
              /* else if(detID == 32){
                 huv_cluster_det32->Fill(u, v);
               }*/
               else if(detID == 5) {
                 huv_cluster_det5->Fill(u, v);
               }
               else if(detID == 7) {
                 huv_cluster_det7->Fill(u, v);
               }
               else if(detID == 9) {
                 huv_cluster_det9->Fill(u, v);
               }
            }           
          }
          if(aMatchedMeasHit){
            totalMatchedHits++;
            double u = aMatchedMeasHit->u();
            double v = aMatchedMeasHit->v();
            if(detID == 32){
              huv_cluster_det32->Fill(u, v);
            }
          }
        }
      }
    }
  }
    
  std::fprintf(stdout, "Good  tracks (==5 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
 
 
  TCanvas *c1 = new TCanvas("c1", "Cluster UV Distribution by Detector", 1800, 1200);
  c1->Divide(3, 3);
 
  c1->cd(1);
  huv_cluster_det2->SetMarkerColor(kRed);
  huv_cluster_det2->SetStats(1);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.9);
  huv_cluster_det2->GetXaxis()->SetTitle("U");
  huv_cluster_det2->GetYaxis()->SetTitle("V");

  huv_cluster_det2->Draw("COLZ");
 
  c1->cd(2);
  huv_cluster_det3->SetMarkerColor(kBlue);
  huv_cluster_det3->GetXaxis()->SetTitle("U");
  huv_cluster_det3->GetYaxis()->SetTitle("V");
  huv_cluster_det3->Draw("COLZ");
 
  c1->cd(3);
  huv_cluster_det32->SetMarkerColor(kViolet);
  huv_cluster_det32->GetXaxis()->SetTitle("U");
  huv_cluster_det32->GetYaxis()->SetTitle("V");
  huv_cluster_det32->Draw("COLZ");
 
  c1->cd(4);
  huv_cluster_det5->SetMarkerColor(kGreen);
  huv_cluster_det5->GetXaxis()->SetTitle("U");
  huv_cluster_det5->GetYaxis()->SetTitle("V");
  huv_cluster_det5->Draw("COLZ");
 
  c1->cd(5);
  huv_cluster_det7->SetMarkerColor(kMagenta);
  huv_cluster_det7->GetXaxis()->SetTitle("U");
  huv_cluster_det7->GetYaxis()->SetTitle("V");
  huv_cluster_det7->Draw("COLZ");
 
  c1->cd(6);
  huv_cluster_det9->SetMarkerColor(kOrange);
  huv_cluster_det9->GetXaxis()->SetTitle("U");
  huv_cluster_det9->GetYaxis()->SetTitle("V");
  huv_cluster_det9->Draw("COLZ");
 
  c1->SaveAs("Cluster_UV_Distribution.svg");
 
  TFile *outFile = new TFile("cluster_uv_distribution.root", "RECREATE");
  huv_cluster_det2->Write();
  huv_cluster_det3->Write();
  huv_cluster_det32->Write();
  huv_cluster_det5->Write();
  huv_cluster_det7->Write();
  huv_cluster_det9->Write();
  c1->Write();
//  outFile->Close();

  //tfile->Close();
 // delete tfile;
 // std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to origin_hits_by_detector.root\n");
  std::fprintf(stdout, "Press any key to exit...\n");
//  std::getc(stdin);
  return;
}
