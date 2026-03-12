void Match(const std::string& rootFilePath){
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
 

  TH1F *hNumMatchedHits_det32  = new TH1F("hNumMatchedHits_det32", "good_Matched_Cluster_size_det 32", 50, 0,25);

  TH1F *hResidualU_det32_num1  = new TH1F("hResidualU_det32_num1", "Residual U (det32, numMatched=1);u_{meas}-u_{fit} [um];Entries [100bin]", 100, -100, 100);
  TH1F *hResidualV_det32_num1  = new TH1F("hResidualV_det32_num1", "Residual V (det32, numMatched=1);v_{meas}-v_{fit} [um];Entries [100bin]", 100, -100, 100);
  TH1F *hResidualU_det32_num2  = new TH1F("hResidualU_det32_num2", "Residual U (det32, numMatched=2);u_{meas}-u_{fit} [um];Entries [100bin]", 100, -100, 100);
  TH1F *hResidualV_det32_num2  = new TH1F("hResidualV_det32_num2", "Residual V (det32, numMatched=2);v_{meas}-v_{fit} [um];Entries [100bin]", 100, -100, 100);
  TH1F *hResidualU_det32_other = new TH1F("hResidualU_det32_other", "Residual U (det32, numMatched>2);u_{meas}-u_{fit} [um];Entries [100bin]", 100, -100, 100);
  TH1F *hResidualV_det32_other = new TH1F("hResidualV_det32_other", "Residual V (det32, numMatched>2);v_{meas}-v_{fit} [um];Entries [100bin]", 100, -100, 100);

  TH1F *hClusterWidthU_det32 = new TH1F("hClusterWidthU_det32", "Cluster Width  in U direction (det32); Cluster Width_U ;Entries", 30, 0, 10);
  TH1F *hClusterWidthV_det32 = new TH1F("hClusterWidthV_det32", "Cluster Height in V direction (det32); Cluster Width_V ;Entries", 30, 0, 10);

  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits=0;

  size_t centerPixelResponses = 0;
  size_t neighborPixelResponses = 0;
  size_t totalResponses = 0;
  
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
  
   /* int goodTrajectoryCount = 0;
    for(auto aTraj: telEvent->trajs()){
      if(aTraj->numOriginMeasHit() == 5){
        goodTrajectoryCount++;
      }
    }
    if(goodTrajectoryCount != 1){
      continue;
    }
  */
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() >=3); 
 
      if(!isGoodTrack){
        continue;
      }
   //   goodTracks++;                
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
        if(!aFitHit ){
          continue;
        }
        if (aFitHit->detN() != 2){
          continue;
        }
        double fit_u=aFitHit->u();
        double fit_v=aFitHit->v();
        if(fit_u < 0 || fit_u > 4.85 || fit_v < 0 || fit_v > 8.19 ){
          //Invalid_tracks++;
          continue;
        //  std::fprintf(stdout,"(hit_u,hit_v) :(%f,%f)\n", hit_u,hit_v);
        }
        goodTracks++;
        if(aMatchedMeasHit){          
          totalMatchedHits++;
          size_t numMatched = aMatchedMeasHit->measRaws().size();
          double cluster_u = aMatchedMeasHit->u();
          double cluster_v = aMatchedMeasHit->v();
          double Residual_u=(cluster_u-fit_u)*1000.0;
          double Residual_v=(cluster_v-fit_v)*1000.0;
          hNumMatchedHits_det32->Fill(numMatched);
        //  hResidualU_det32_num1->Fill(Residual_u);
        //  hResidualV_det32_num1->Fill(Residual_v);
          //----------------------------------------------------------
          std::set<uint16_t> unique_u_coords;
          std::set<uint16_t> unique_v_coords;
          for(const auto& raw : aMatchedMeasHit->measRaws()){
            unique_u_coords.insert(raw.u());
            unique_v_coords.insert(raw.v());
          }
          size_t num_unique_u = unique_u_coords.size();
          size_t num_unique_v = unique_v_coords.size();

          hClusterWidthU_det32->Fill(num_unique_u);
          hClusterWidthV_det32->Fill(num_unique_v);
         //-----------------------------------------------------------
          if(numMatched == 1){
            hResidualU_det32_num1->Fill(Residual_u);
            hResidualV_det32_num1->Fill(Residual_v);
          } else if(numMatched == 2){
            std::set<uint16_t> unique_u_coords1;
            std::set<uint16_t> unique_v_coords1;
            for(const auto& raw : aMatchedMeasHit->measRaws()){
              unique_u_coords1.insert(raw.u());
              unique_v_coords1.insert(raw.v());
      //        std::printf(" Original pixel coordinates: u=%d, v=%d\n", raw.u(), raw.v());
            }
            size_t num_unique_u1 = unique_u_coords1.size();
            size_t num_unique_v1 = unique_v_coords1.size();
            if(num_unique_u1 == 1 && num_unique_v1 == 2){
              hResidualU_det32_num1->Fill(Residual_u);
              hResidualV_det32_num2->Fill(Residual_v);
            }
            else if(num_unique_u1 == 2 && num_unique_v1 == 1){
              hResidualU_det32_num2->Fill(Residual_u);
              hResidualV_det32_num1->Fill(Residual_v);
            }
            else{
              hResidualU_det32_num2->Fill(Residual_u);
              hResidualV_det32_num2->Fill(Residual_v);
            }
          }else if(numMatched>2) {
            std::set<uint16_t> unique_u_coords2;
            std::set<uint16_t> unique_v_coords2;
            for(const auto& raw : aMatchedMeasHit->measRaws()){
              unique_u_coords2.insert(raw.u());
              unique_v_coords2.insert(raw.v());
            }
            size_t num_unique_u2 = unique_u_coords2.size();
            size_t num_unique_v2 = unique_v_coords2.size();
            if(num_unique_u2 == 1 && num_unique_v2>2){
              hResidualU_det32_num1->Fill(Residual_u);
              hResidualV_det32_other->Fill(Residual_v);
            }
            else if(num_unique_u2>2 && num_unique_v2==1){
              hResidualU_det32_other->Fill(Residual_u);
              hResidualV_det32_num1->Fill(Residual_v);
            }
            else if(num_unique_u2==2 && num_unique_v2>2){
              hResidualU_det32_num2->Fill(Residual_u);
              hResidualV_det32_other->Fill(Residual_v);
            }
            else if(num_unique_u2>2 && num_unique_v2==2){
              hResidualU_det32_other->Fill(Residual_u);
              hResidualV_det32_num2->Fill(Residual_v);
            }
            else{
            hResidualU_det32_other->Fill(Residual_u);
            hResidualV_det32_other->Fill(Residual_v);
            }                                                
          }
        }
      }
    }
  }
 
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good  tracks (=5 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  

  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 1800, 1200);
  c1->Divide(2, 4);
 
  c1->cd(1);
  gStyle->SetPaintTextFormat(".2f");
  hNumMatchedHits_det32->SetFillColor(kGreen);
  if(hNumMatchedHits_det32->GetEntries() > 0){
    hNumMatchedHits_det32->Scale(1.0 / hNumMatchedHits_det32->GetEntries());
  }
  hNumMatchedHits_det32->SetStats(1);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.9);
  hNumMatchedHits_det32->GetXaxis()->SetTitle("good_Matched_Cluster_size");
  hNumMatchedHits_det32->GetYaxis()->SetTitle("Entries(Normalization)");
  hNumMatchedHits_det32->SetMinimum(0);
//  hNumMatchedHits_det32->SetMaximum(1.0);
  hNumMatchedHits_det32->SetFillStyle(3001);
  hNumMatchedHits_det32->SetMarkerSize(2.5);
  hNumMatchedHits_det32->Draw("HIST TEXT0");

  c1->cd(3);
  hResidualU_det32_num1->SetLineColor(kBlue+2);
  hResidualU_det32_num1->SetLineWidth(2);
  hResidualU_det32_num1->Draw();
  TF1 *fitU1 = new TF1("fitU1", "gaus", -50, 50);
  fitU1->SetLineColor(kRed);
  hResidualU_det32_num1->Fit(fitU1, "R");

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.07);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitU1->GetParameter(1)));
  latex1.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitU1->GetParameter(2)));
  
  c1->cd(4);
  hResidualV_det32_num1->SetLineColor(kBlue+2);
  hResidualV_det32_num1->SetLineWidth(2);
  hResidualV_det32_num1->Draw();
  TF1 *fitV1 = new TF1("fitV1", "gaus", -50, 50);
  fitV1->SetLineColor(kRed);
  hResidualV_det32_num1->Fit(fitV1, "R");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.07);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitV1->GetParameter(1)));
  latex2.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitV1->GetParameter(2)));  
  
  c1->cd(5);
  hResidualU_det32_num2->SetLineColor(kBlue+2);
  hResidualU_det32_num2->SetLineWidth(2);
  hResidualU_det32_num2->Draw();
  TF1 *fitU2 = new TF1("fitU2", "gaus", -50, 50);
  fitU2->SetLineColor(kRed);
  hResidualU_det32_num2->Fit(fitU2, "R");

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.07);
  latex3.SetTextColor(kRed);
  latex3.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitU2->GetParameter(1)));
  latex3.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitU2->GetParameter(2)));

  c1->cd(6);
  hResidualV_det32_num2->SetLineColor(kBlue+2);
  hResidualV_det32_num2->SetLineWidth(2);
  hResidualV_det32_num2->Draw();
  TF1 *fitV2 = new TF1("fitV2", "gaus", -50, 50);
  fitV2->SetLineColor(kRed);
  hResidualV_det32_num2->Fit(fitV2, "R");

  TLatex latex4;
  latex4.SetNDC();
  latex4.SetTextSize(0.07);
  latex4.SetTextColor(kRed);
  latex4.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitV2->GetParameter(1)));
  latex4.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitV2->GetParameter(2)));

  c1->cd(7);
  hResidualU_det32_other->SetLineColor(kBlue+2);
  hResidualU_det32_other->SetLineWidth(2);
  hResidualU_det32_other->Draw();
  TF1 *fitU3 = new TF1("fitU3", "gaus", -50, 50);
  fitU3->SetLineColor(kRed);
  hResidualU_det32_other->Fit(fitU3, "R");

  TLatex latex5;
  latex5.SetNDC();
  latex5.SetTextSize(0.07);
  latex5.SetTextColor(kRed);
  latex5.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitU3->GetParameter(1)));
  latex5.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitU3->GetParameter(2)));

  c1->cd(8);
  hResidualV_det32_other->SetLineColor(kBlue+2);
  hResidualV_det32_other->SetLineWidth(2);
  hResidualV_det32_other->Draw();
  TF1 *fitV3 = new TF1("fitV3", "gaus", -50, 50);
  fitV3->SetLineColor(kRed);
  hResidualV_det32_other->Fit(fitV3, "R");

  TLatex latex6;
  latex6.SetNDC();
  latex6.SetTextSize(0.07);
  latex6.SetTextColor(kRed);
  latex6.DrawLatex(0.55, 0.80, Form("Mean = %.1fum", fitV3->GetParameter(1)));
  latex6.DrawLatex(0.55, 0.70, Form("Sigma = %.1fum", fitV3->GetParameter(2)));

  c1->Update();

  double sigmaU1 = fitU1->GetParameter(2); // for numMatched==1
  double sigmaV1 = fitV1->GetParameter(2);

  double sigmaU2 = fitU2->GetParameter(2); // for numMatched==2
  double sigmaV2 = fitV2->GetParameter(2);

  double sigmaU3 = fitU3->GetParameter(2); // for numMatched>2
  double sigmaV3 = fitV3->GetParameter(2);
  
  std::ofstream sigmaFile("/home/pub/B/altel_acts/Draw/root_file/match_sigma_output.txt");
  if (sigmaFile.is_open()) {
    std::cout << std::fixed << std::setprecision(1);
    sigmaFile << "SigmaU1 " << sigmaU1 << std::endl;
    sigmaFile << "SigmaV1 " << sigmaV1 << std::endl;
    sigmaFile << "SigmaU2 " << sigmaU2 << std::endl;
    sigmaFile << "SigmaV2 " << sigmaV2 << std::endl;
    sigmaFile << "SigmaU3 " << sigmaU3 << std::endl;
    sigmaFile << "SigmaV3 " << sigmaV3 << std::endl;
    sigmaFile.close();
  }
  std::fprintf(stdout, "SigmaU1 = %.1f um, SigmaV1 = %.1f um\n", sigmaU1, sigmaV1);
  std::fprintf(stdout, "SigmaU2 = %.1f um, SigmaV2 = %.1f um\n", sigmaU2, sigmaV2);
  std::fprintf(stdout, "SigmaU3 = %.1f um, SigmaV3 = %.1f um\n", sigmaU3, sigmaV3);
  
  //c1->SaveAs("/home/pub/B/altel_acts/Draw/save_svg/match_hits_by_detector.svg");

  TCanvas *c2 = new TCanvas("c2", "Origin Hits by Detector Analysis", 1800, 1200);
  c2->Divide(1, 2);
  c2->cd(1);
  hClusterWidthU_det32->SetFillColor(kOrange);
 // hClusterWidthU_det32->SetLineColor(kBlack);
  hClusterWidthU_det32->SetMarkerSize(2.5);
  hClusterWidthU_det32->Draw("HIST TEXT0");

  c2->cd(2);
  hClusterWidthV_det32->SetFillColor(kAzure+1);
//  hClusterWidthV_det32->SetLineColor(kBlack);
  hClusterWidthV_det32->SetMarkerSize(2.5);
  //gStyle->SetTextAngle(0);
  hClusterWidthV_det32->Draw("HIST TEXT0");
  c2->Update();

  TFile *outFile = new TFile("/home/pub/B/altel_acts/Draw/root_file/match_hits_by_detector.root", "RECREATE");
  hNumMatchedHits_det32->Write();
  hResidualU_det32_num1->Write();
  hResidualV_det32_num1->Write();
  hResidualU_det32_num2->Write();
  hResidualV_det32_num2->Write();
  hResidualU_det32_other->Write();
  hResidualV_det32_other->Write();
  hClusterWidthU_det32->Write();
  hClusterWidthV_det32->Write();
  c1->Write();
  c2->Write();
  
  std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to /home/pub/B/altel_acts/Draw/save_svg/match_hits_by_detector.root\n");
 
  return;
}
