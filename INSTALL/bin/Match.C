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
 

  TH1F *hNumMatchedHits_det32  = new TH1F("hNumMatchedHits_det32", "good_Matched_Cluster_size_det 32", 30, 0,15);
 
  TH1F *hResidualU_det32_num1  = new TH1F("hResidualU_det32_num1", "Residual U (det32, numMatched=1);u_{meas}-u_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  TH1F *hResidualV_det32_num1  = new TH1F("hResidualV_det32_num1", "Residual V (det32, numMatched=1);v_{meas}-v_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  TH1F *hResidualU_det32_num2  = new TH1F("hResidualU_det32_num2", "Residual U (det32, numMatched=2);u_{meas}-u_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  TH1F *hResidualV_det32_num2  = new TH1F("hResidualV_det32_num2", "Residual V (det32, numMatched=2);v_{meas}-v_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  TH1F *hResidualU_det32_other = new TH1F("hResidualU_det32_other", "Residual U (det32, numMatched>2);u_{meas}-u_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  TH1F *hResidualV_det32_other = new TH1F("hResidualV_det32_other", "Residual V (det32, numMatched>2);v_{meas}-v_{fit} [mm];Entries [100bin]", 100, -0.1, 0.1);
  
  TH1F* hu_err =new TH1F("hu_err","u error;u_err[um];Entries",150, 0,10);
  TH1F* hv_err =new TH1F("hv_err","v error;v_err[um];Entries",150, 0,10);


  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits=0;

  size_t centerPixelResponses = 0;
  size_t neighborPixelResponses = 0;
  size_t totalResponses = 0;
  std::map<uint32_t, size_t> pixelHitCount;  //The number of hits per pixel
  std::vector<uint32_t> allPixelIds;         //Record the numbers of all response pixels
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==5); 
 
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
        double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);
     //   std::fprintf(stdout, "hit_u_err : %.8f\n", hit_u_err);
    //    std::fprintf(stdout, "hit_v_err : %.8f\n", hit_v_err);
        if (aFitHit->detN() != 32){
          continue;
        }
        double fit_u=aFitHit->u();
        double fit_v=aFitHit->v();
      /*  double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);*/
        if(fit_u < -12.8 || fit_u > 12.8 || fit_v < -6.4 || fit_v > 6.4 ){
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
          double Residual_u=cluster_u-fit_u;
          double Residual_v=cluster_v-fit_v;
          hNumMatchedHits_det32->Fill(numMatched);
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
  
  std::fprintf(stdout, "\n=== Response Position Analysis ===\n");
  std::fprintf(stdout, "Total responses analyzed: %zu\n", totalResponses);
  std::fprintf(stdout, "Center pixel responses: %zu (%.2f%%)\n",
                centerPixelResponses, totalResponses > 0 ? (double)centerPixelResponses/totalResponses*100 : 0);
  std::fprintf(stdout, "Neighbor pixel responses: %zu (%.2f%%)\n",
                neighborPixelResponses, totalResponses > 0 ? (double)neighborPixelResponses/totalResponses*100 : 0);
  std::fprintf(stdout, "Far responses: %zu (%.2f%%)\n",
                totalResponses - centerPixelResponses - neighborPixelResponses,
                totalResponses > 0 ? (double)(totalResponses - centerPixelResponses - neighborPixelResponses)/totalResponses*100 : 0);


  std::fprintf(stdout, "\n=== Pixel ID Analysis ===\n");
  std::fprintf(stdout, "Total pixel responses: %zu\n", allPixelIds.size());
  std::fprintf(stdout, "Unique pixels fired: %zu\n", pixelHitCount.size());
  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 1800, 1200);
  c1->Divide(2, 4);
 
  c1->cd(1);
  hNumMatchedHits_det32->SetFillColor(kGreen);
 /* if(hNumMatchedHits_det32->GetEntries() > 0){
    hNumMatchedHits_det32->Scale(1.0 / hNumMatchedHits_det32->GetEntries());
  }*/
  hNumMatchedHits_det32->SetStats(1);
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.9);
  hNumMatchedHits_det32->GetXaxis()->SetTitle("good_Matched_Cluster_size");
  hNumMatchedHits_det32->GetYaxis()->SetTitle("Percentage");
  hNumMatchedHits_det32->SetMinimum(0);
//  hNumMatchedHits_det32->SetMaximum(1.0);
  hNumMatchedHits_det32->SetFillStyle(3001);
  hNumMatchedHits_det32->Draw("HIST TEXT");

  c1->cd(3);
  hResidualU_det32_num1->SetLineColor(kBlue+2);
  hResidualU_det32_num1->SetLineWidth(2);
  hResidualU_det32_num1->Draw();
  TF1 *fitU1 = new TF1("fitU1", "gaus", -0.05, 0.05);
  fitU1->SetLineColor(kRed);
  hResidualU_det32_num1->Fit(fitU1, "R");

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.07);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU1->GetParameter(1)));
  latex1.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU1->GetParameter(2)));
  
  c1->cd(4);
  hResidualV_det32_num1->SetLineColor(kBlue+2);
  hResidualV_det32_num1->SetLineWidth(2);
  hResidualV_det32_num1->Draw();
  TF1 *fitV1 = new TF1("fitV1", "gaus", -0.05, 0.05);
  fitV1->SetLineColor(kRed);
  hResidualV_det32_num1->Fit(fitV1, "R");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.07);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV1->GetParameter(1)));
  latex2.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV1->GetParameter(2)));  
  
  c1->cd(5);
  hResidualU_det32_num2->SetLineColor(kBlue+2);
  hResidualU_det32_num2->SetLineWidth(2);
  hResidualU_det32_num2->Draw();
  TF1 *fitU2 = new TF1("fitU2", "gaus", -0.05, 0.05);
  fitU2->SetLineColor(kRed);
  hResidualU_det32_num2->Fit(fitU2, "R");

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.07);
  latex3.SetTextColor(kRed);
  latex3.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU2->GetParameter(1)));
  latex3.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU2->GetParameter(2)));

  c1->cd(6);
  hResidualV_det32_num2->SetLineColor(kBlue+2);
  hResidualV_det32_num2->SetLineWidth(2);
  hResidualV_det32_num2->Draw();
  TF1 *fitV2 = new TF1("fitV2", "gaus", -0.05, 0.05);
  fitV2->SetLineColor(kRed);
  hResidualV_det32_num2->Fit(fitV2, "R");

  TLatex latex4;
  latex4.SetNDC();
  latex4.SetTextSize(0.07);
  latex4.SetTextColor(kRed);
  latex4.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV2->GetParameter(1)));
  latex4.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV2->GetParameter(2)));

  c1->cd(7);
  hResidualU_det32_other->SetLineColor(kBlue+2);
  hResidualU_det32_other->SetLineWidth(2);
  hResidualU_det32_other->Draw();
  TF1 *fitU3 = new TF1("fitU3", "gaus", -0.05, 0.05);
  fitU3->SetLineColor(kRed);
  hResidualU_det32_other->Fit(fitU3, "R");

  TLatex latex5;
  latex5.SetNDC();
  latex5.SetTextSize(0.07);
  latex5.SetTextColor(kRed);
  latex5.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU3->GetParameter(1)));
  latex5.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU3->GetParameter(2)));

  c1->cd(8);
  hResidualV_det32_other->SetLineColor(kBlue+2);
  hResidualV_det32_other->SetLineWidth(2);
  hResidualV_det32_other->Draw();
  TF1 *fitV3 = new TF1("fitV3", "gaus", -0.05, 0.05);
  fitV3->SetLineColor(kRed);
  hResidualV_det32_other->Fit(fitV3, "R");

  TLatex latex6;
  latex6.SetNDC();
  latex6.SetTextSize(0.07);
  latex6.SetTextColor(kRed);
  latex6.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV3->GetParameter(1)));
  latex6.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV3->GetParameter(2)));

  c1->Update();
 
  c1->SaveAs("match_hits_by_detector.svg");
  
  TFile *outFile = new TFile("match_hits_by_detector.root", "RECREATE");
  hNumMatchedHits_det32->Write();
  hResidualU_det32_num1->Write();
  hResidualV_det32_num1->Write();
  hResidualU_det32_num2->Write();
  hResidualV_det32_num2->Write();
  hResidualU_det32_other->Write();
  hResidualV_det32_other->Write();
  c1->Write();
  
  TCanvas *c2 = new TCanvas("c2", "Origin Hits by Detector Analysis", 1800, 1200);

  c2->Divide(2, 1);
  c2->cd(1);
  hu_err->SetLineColor(kBlue+2);
  hu_err->SetLineWidth(2);
  hu_err->Draw();
 /* TF1 *fitU4 = new TF1("fitU4", "gaus", 2, 8);
  fitU4->SetLineColor(kRed);
  hu_err->Fit(fitU4, "R");

  TLatex latex7;
  latex7.SetNDC();
  latex7.SetTextSize(0.04);
  latex7.SetTextColor(kRed);
  latex7.DrawLatex(0.13, 0.80, Form("Mean = %.4f", fitU4->GetParameter(1)));
  latex7.DrawLatex(0.13, 0.70, Form("Sigma = %.4f", fitU4->GetParameter(2)));
*/
  c2->cd(2);
  hv_err->SetLineColor(kBlue+2);
  hv_err->SetLineWidth(2);
  hv_err->Draw();
 /* TF1 *fitV4 = new TF1("fitV4", "gaus", 0, 10);
  fitV4->SetLineColor(kRed);
  hv_err->Fit(fitV4, "R");  

  TLatex latex8;
  latex8.SetNDC();
  latex8.SetTextSize(0.04);
  latex8.SetTextColor(kRed);
  latex8.DrawLatex(0.13, 0.80, Form("Mean = %.4f", fitV4->GetParameter(1)));
  latex8.DrawLatex(0.13, 0.70, Form("Sigma = %.4f", fitV4->GetParameter(2)));
*/
  std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to origin_hits_by_detector.root\n");
 
  return;
}
