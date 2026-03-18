void match(const std::string& rootFilePath){
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
 

  TH1F *hNumMatchedHits_det2  = new TH1F("hNumMatchedHits_det2", "good_Matched_Cluster_size_det 2", 15, 0,15);

  TH1F *hResidualU_det2_num1  = new TH1F("hResidualU_det2_num1", "Residual U (det2, numMatched=1);u_{meas}-u_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualV_det2_num1  = new TH1F("hResidualV_det2_num1", "Residual V (det2, numMatched=1);v_{meas}-v_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualU_det2_num2  = new TH1F("hResidualU_det2_num2", "Residual U (det2, numMatched=2);u_{meas}-u_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualV_det2_num2  = new TH1F("hResidualV_det2_num2", "Residual V (det2, numMatched=2);v_{meas}-v_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualU_det2_num3  = new TH1F("hResidualU_det2_num3", "Residual U (det2, numMatched=3);u_{meas}-u_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualV_det2_num3  = new TH1F("hResidualV_det2_num3", "Residual V (det2, numMatched=3);v_{meas}-v_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualU_det2_num4  = new TH1F("hResidualU_det2_num4", "Residual U (det2, numMatched=4);u_{meas}-u_{fit} [um];Entries [200bin]", 200, -50, 50);
  TH1F *hResidualV_det2_num4  = new TH1F("hResidualV_det2_num4", "Residual V (det2, numMatched=4);v_{meas}-v_{fit} [um];Entries [200bin]", 200, -50, 50);



  TH1F *hClusterWidthU_det2 = new TH1F("hClusterWidthU_det2", "Cluster Width  in U/V direction (det2); Cluster Width_U/V ;Entries", 10, 0, 10);
  TH1F *hClusterWidthV_det2 = new TH1F("hClusterWidthV_det2", "Cluster Height in V direction (det2); Cluster Width_V ;Entries", 10, 0, 10);

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
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==4); 
 
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
          hNumMatchedHits_det2->Fill(numMatched);
        //  hResidualU_det2_num1->Fill(Residual_u);
        //  hResidualV_det2_num1->Fill(Residual_v);
          //----------------------------------------------------------
          std::set<uint16_t> unique_u_coords;
          std::set<uint16_t> unique_v_coords;
          for(const auto& raw : aMatchedMeasHit->measRaws()){
            unique_u_coords.insert(raw.u());
            unique_v_coords.insert(raw.v());
          }
          size_t num_unique_u = unique_u_coords.size();
          size_t num_unique_v = unique_v_coords.size();

          hClusterWidthU_det2->Fill(num_unique_u);
          hClusterWidthV_det2->Fill(num_unique_v);
        
          if(num_unique_u == 1){
            hResidualU_det2_num1->Fill(Residual_u);
          }
          else if(num_unique_u == 2){
            hResidualU_det2_num2->Fill(Residual_u);
          }
          else if(num_unique_u == 3){
            hResidualU_det2_num3->Fill(Residual_u);
          }
          else if(num_unique_u == 4){
            hResidualU_det2_num4->Fill(Residual_u);
          }

          if(num_unique_v == 1){
            hResidualV_det2_num1->Fill(Residual_v);
          }
          else if(num_unique_v == 2){
            hResidualV_det2_num2->Fill(Residual_v);
          }
          else if(num_unique_v == 3){
            hResidualV_det2_num3->Fill(Residual_v);
          }
          else if(num_unique_v == 4) {
            hResidualV_det2_num4->Fill(Residual_v);
          }
        }
      }
    }
  }
 
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good  tracks (=4 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  

  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 1800, 1200);
 
  gStyle->SetPaintTextFormat(".3f");
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.9);
  hNumMatchedHits_det2->SetStats(1);
  hNumMatchedHits_det2->SetFillColor(kGreen);
  if(hNumMatchedHits_det2->GetEntries() > 0){
    hNumMatchedHits_det2->Scale(1.0 / hNumMatchedHits_det2->GetEntries());
  }
  hNumMatchedHits_det2->GetXaxis()->SetTitle("good_Matched_Cluster_size");
  hNumMatchedHits_det2->GetYaxis()->SetTitle("Entries(Normalization)");
  hNumMatchedHits_det2->SetMinimum(0);
//  hNumMatchedHits_det2->SetMaximum(1.0);
  hNumMatchedHits_det2->SetFillStyle(3001);
  hNumMatchedHits_det2->SetMarkerSize(1);
  hNumMatchedHits_det2->Draw("HIST TEXT0");
  
  TCanvas *c2 = new TCanvas("c2", "Origin Hits by Detector Analysis", 1000, 1200);
  c2->Divide(2, 4);

  c2->cd(1);
  hResidualU_det2_num1->SetLineColor(kBlue+2);
  hResidualU_det2_num1->SetLineWidth(2);
  hResidualU_det2_num1->SetMarkerStyle(20);
  hResidualU_det2_num1->SetMarkerSize(1);
  hResidualU_det2_num1->Draw("P");
  TF1 *fitU1 = new TF1("fitU1", "gaus", -50, 50);
  fitU1->SetLineColor(kRed);
  hResidualU_det2_num1->Fit(fitU1, "R");

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.05);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitU1->GetParameter(1)));
  latex1.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitU1->GetParameter(2)));
  
  c2->cd(2);
  hResidualV_det2_num1->SetLineColor(kBlue+2);
  hResidualV_det2_num1->SetLineWidth(2);
  hResidualV_det2_num1->SetMarkerStyle(20);
  hResidualV_det2_num1->SetMarkerSize(1);
  hResidualV_det2_num1->Draw("P");
  TF1 *fitV1 = new TF1("fitV1", "gaus", -30, 30);
  fitV1->SetLineColor(kRed);
  hResidualV_det2_num1->Fit(fitV1, "R");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.05);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitV1->GetParameter(1)));
  latex2.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitV1->GetParameter(2)));  
  
  c2->cd(3);
  hResidualU_det2_num2->SetLineColor(kBlue+2);
  hResidualU_det2_num2->SetLineWidth(2);
  hResidualU_det2_num2->SetMarkerStyle(20);
  hResidualU_det2_num2->SetMarkerSize(1);
  hResidualU_det2_num2->Draw("P");
  TF1 *fitU2 = new TF1("fitU2", "gaus", -30, 30);
  fitU2->SetLineColor(kRed);
  hResidualU_det2_num2->Fit(fitU2, "R");

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.05);
  latex3.SetTextColor(kRed);
  latex3.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitU2->GetParameter(1)));
  latex3.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitU2->GetParameter(2)));

  c2->cd(4);
  hResidualV_det2_num2->SetLineColor(kBlue+2);
  hResidualV_det2_num2->SetLineWidth(2);
  hResidualV_det2_num2->SetMarkerStyle(20);
  hResidualV_det2_num2->SetMarkerSize(1);
  hResidualV_det2_num2->Draw("P");
  TF1 *fitV2 = new TF1("fitV2", "gaus", -30, 30);
  fitV2->SetLineColor(kRed);
  hResidualV_det2_num2->Fit(fitV2, "R");

  TLatex latex4;
  latex4.SetNDC();
  latex4.SetTextSize(0.05);
  latex4.SetTextColor(kRed);
  latex4.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitV2->GetParameter(1)));
  latex4.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitV2->GetParameter(2)));

  c2->cd(5);
  hResidualU_det2_num3->SetLineColor(kBlue+2);
  hResidualU_det2_num3->SetLineWidth(2);
  hResidualU_det2_num3->SetMarkerStyle(20);
  hResidualU_det2_num3->SetMarkerSize(1);
  hResidualU_det2_num3->Draw("P");
  TF1 *fitU3 = new TF1("fitU3", "gaus", -30, 30);
  fitU3->SetLineColor(kRed);
  hResidualU_det2_num3->Fit(fitU3, "R");

  TLatex latex5;
  latex5.SetNDC();
  latex5.SetTextSize(0.05);
  latex5.SetTextColor(kRed);
  latex5.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitU3->GetParameter(1)));
  latex5.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitU3->GetParameter(2)));

  c2->cd(6);
  hResidualV_det2_num3->SetLineColor(kBlue+2);
  hResidualV_det2_num3->SetLineWidth(2);
  hResidualV_det2_num3->SetMarkerStyle(20);
  hResidualV_det2_num3->SetMarkerSize(1);
  hResidualV_det2_num3->Draw("P");
  TF1 *fitV3 = new TF1("fitV3", "gaus", -30, 30);
  fitV3->SetLineColor(kRed);
  hResidualV_det2_num3->Fit(fitV3, "R");

  TLatex latex6;
  latex6.SetNDC();
  latex6.SetTextSize(0.05);
  latex6.SetTextColor(kRed);
  latex6.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitV3->GetParameter(1)));
  latex6.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitV3->GetParameter(2)));


  c2->cd(7);
  hResidualU_det2_num4->SetLineColor(kBlue+2);
  hResidualU_det2_num4->SetLineWidth(2);
  hResidualU_det2_num4->SetMarkerStyle(20);
  hResidualU_det2_num4->SetMarkerSize(1);
  hResidualU_det2_num4->Draw("P");
  TF1 *fitU4 = new TF1("fitU4", "gaus", -30, 30);
  fitU4->SetLineColor(kRed);
  hResidualU_det2_num4->Fit(fitU4, "R");

  TLatex latex7;
  latex7.SetNDC();
  latex7.SetTextSize(0.05);
  latex7.SetTextColor(kRed);
  latex7.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitU4->GetParameter(1)));
  latex7.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitU4->GetParameter(2)));

  c2->cd(8);
  hResidualV_det2_num4->SetLineColor(kBlue+2);
  hResidualV_det2_num4->SetLineWidth(2);
  hResidualV_det2_num4->SetMarkerStyle(20);
  hResidualV_det2_num4->SetMarkerSize(1);
  hResidualV_det2_num4->Draw("P");
  TF1 *fitV4 = new TF1("fitV4", "gaus", -30, 30);
  fitV4->SetLineColor(kRed);
  hResidualV_det2_num4->Fit(fitV4, "R");

  TLatex latex8;
  latex8.SetNDC();
  latex8.SetTextSize(0.05);
  latex8.SetTextColor(kRed);
  latex8.DrawLatex(0.65, 0.60, Form("Mean = %.2fum", fitV4->GetParameter(1)));
  latex8.DrawLatex(0.65, 0.50, Form("Sigma = %.2fum", fitV4->GetParameter(2)));


  c2->Update();
  
  //c1->SaveAs("/home/pub/B/altel_acts/Draw/save_svg/match_hits_by_detector.svg");

  TCanvas *c3 = new TCanvas("c3", "Origin Hits by Detector Analysis", 1800, 1200);
  gStyle->SetOptStat(0); 
  hClusterWidthU_det2->SetFillColor(kOrange-9);
  hClusterWidthU_det2->SetLineColor(kOrange+2);
  hClusterWidthU_det2->SetLineWidth(2);
  hClusterWidthU_det2->SetFillStyle(3354);
  if(hClusterWidthU_det2->GetEntries() > 0){
    hClusterWidthU_det2->Scale(1.0 / hClusterWidthU_det2->GetEntries());
  }
  hClusterWidthU_det2->Draw("HIST ");

  hClusterWidthV_det2->SetFillColor(kAzure-9);
  hClusterWidthV_det2->SetLineColor(kAzure+2);
  hClusterWidthV_det2->SetLineWidth(2);
  hClusterWidthV_det2->SetFillStyle(3345);
  if(hClusterWidthV_det2->GetEntries() > 0){
    hClusterWidthV_det2->Scale(1.0 / hClusterWidthV_det2->GetEntries());
  }
  hClusterWidthV_det2->Draw("HIST  SAME");

  //再画一遍边框，让轮廓更清楚
  hClusterWidthU_det2->Draw("HIST SAME");
  hClusterWidthV_det2->Draw("HIST SAME");

  TLegend *leg = new TLegend(0.60, 0.72, 0.90, 0.88);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(hClusterWidthU_det2, "Cluster Width U", "f");
  leg->AddEntry(hClusterWidthV_det2, "Cluster Height V", "f");
  leg->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);
  latex.DrawLatex(0.50, 0.55, Form("Mean cluster width: %.1f", hClusterWidthU_det2->GetMean()));
  latex.DrawLatex(0.50, 0.47, Form("Mean cluster height: %.1f", hClusterWidthV_det2->GetMean()));
  c3->Update();

  TFile *outFile = new TFile("/home/sungaige/A/altel_acts/Draw/root_file/match_hits_by_detector.root", "RECREATE");
  hNumMatchedHits_det2->Write();
  hResidualU_det2_num1->Write();
  hResidualV_det2_num1->Write();
  hResidualU_det2_num2->Write();
  hResidualV_det2_num2->Write();
  hResidualU_det2_num3->Write();
  hResidualV_det2_num3->Write();
  hResidualU_det2_num4->Write();
  hResidualV_det2_num4->Write();
  hClusterWidthU_det2->Write();
  hClusterWidthV_det2->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to /home/sungaige/A/altel_acts/Draw/save_svg/match_hits_by_detector.root\n");
 
  return;
}
