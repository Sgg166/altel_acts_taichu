void u_v_kinkAngle_err(const std::string& rootFilePath){
 // gSystem->Load("/home/pub/B/altel_acts/INSTALL/lib64/libActsCore.so");
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
 

  TH1F* hu_err_0 =new TH1F("hu_err_0","u error ;u_err[um];Entries",200,3,5.5);
  TH1F* hv_err_0 =new TH1F("hv_err_0","v error ;v_err[um];Entries",200,3,5.5);
  
  TH1F* hu_err =new TH1F("hu_err","u error Angle of deflection[0,0.01415];u_err[um];Entries",200, 3,6);
  TH1F* hv_err =new TH1F("hv_err","v error Angle of deflection[0,0.01415];v_err[um];Entries",200, 3,6);
 
  TH1F* hu_err2 =new TH1F("hu_err2","u error Angle of deflection(0.01415,0.02831];u_err[um];Entries",160, 3,6);
  TH1F* hv_err2 =new TH1F("hv_err2","v error Angle of deflection(0.01415,0.02831];v_err[um];Entries",160, 3,6); 

  TH1F *hkink_angle =  new TH1F("hkink_angle", "Angle of deflection distribution;Angle of deflection[Rad];Entries [100bin]", 50, -0.0001, 0.0006);


  size_t totalNumEvents = ttreeReader.numEvents(); 
  size_t totalFitHits = 0;
  double sum_u_err=0;
  double sum_u_err2=0;
  double sum_v_err=0;
  double sum_v_err2=0;
  double sum_tot_err=0;
  double sum_tot_err2=0;
  size_t nErrPoints=0;
  double min_kinkAngle = 0;
  size_t min_eventNum = 0;
  double max_kinkAngle = 0;
  size_t max_eventNum = 0;
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);
      if(!isGoodTrack){
         continue;
      }
      uint16_t id_before = 9;
      uint16_t id_after = 2;
      auto trajHit_before = aTraj->trajHit(id_before);
      auto trajHit_after = aTraj->trajHit(id_after);
      if(!trajHit_before || !trajHit_after){
        continue;
      }
      auto fitHit_before = trajHit_before->FH;
      auto fitHit_after = trajHit_after->FH;
      if(!fitHit_before || !fitHit_after){
        continue;
      }
      if(!fitHit_before->OM || !fitHit_after->OM){
        continue;
      }
     // Acts::Vector3D dir_before(fitHit_before->DGs[0], fitHit_before->DGs[1], fitHit_before->DGs[2]);
     // Acts::Vector3D dir_after(fitHit_after->DGs[0], fitHit_after->DGs[1], fitHit_after->DGs[2]);
      double dir_before_x = fitHit_before->DGs[0];
      double dir_before_y = fitHit_before->DGs[1]; 
      double dir_before_z = fitHit_before->DGs[2];
      double dir_after_x = fitHit_after->DGs[0];
      double dir_after_y = fitHit_after->DGs[1];
      double dir_after_z = fitHit_after->DGs[2];
 
      double dot_product = dir_before_x*dir_after_x + dir_before_y*dir_after_y + dir_before_z*dir_after_z;
      double norm_before = sqrt(dir_before_x*dir_before_x + dir_before_y*dir_before_y + dir_before_z*dir_before_z);
      double norm_after = sqrt(dir_after_x*dir_after_x + dir_after_y*dir_after_y + dir_after_z*dir_after_z);
      double kinkAngle = acos(dot_product / (norm_before * norm_after));
     // double kinkAngle = norm_after-norm_before;
    //  hkink_angle->Fill(kinkAngle);
      if (kinkAngle > max_kinkAngle) {
        max_kinkAngle = kinkAngle;
        max_eventNum = eventNum;
      }
      if (kinkAngle< min_kinkAngle){
        min_kinkAngle= kinkAngle;
        min_eventNum = eventNum;
      }

     // double kinkAngle = (dir_after-dir_before).norm();
    //  double kinkAngle = acos(dir_before.dot(dir_after) / (dir_before.norm() * dir_after.norm()));
     // std::fprintf(stdout, "kinkAngle =  %.6f",kinkAngle);
      
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
        if(!aFitHit ){
          continue;
        }
        if (aFitHit->detN() != 32){
          continue;
        }
        double fit_u=aFitHit->u();
        double fit_v=aFitHit->v();

        if(fit_u < -12.8 || fit_u > 12.8 || fit_v < -6.4 || fit_v > 6.4 ){
          continue;
        }
        if(!aMatchedMeasHit){
          continue;
        }
        double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hkink_angle->Fill(kinkAngle);
        hu_err_0->Fill(hit_u_err);
        hv_err_0->Fill(hit_v_err);
        if(kinkAngle >= 0 && kinkAngle <= 0.000246995){
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);        
        }
        else if(kinkAngle > 0.000246995 && kinkAngle <= 0.00049399){
          hu_err2->Fill(hit_u_err);
          hv_err2->Fill(hit_v_err);
        }
        else{
          continue;
        }
        double tot_err = std::sqrt(hit_u_err*hit_u_err + hit_v_err*hit_v_err);
        sum_u_err+=hit_u_err; sum_u_err2+=hit_u_err*hit_u_err;
        sum_v_err+=hit_v_err; sum_v_err2+=hit_v_err*hit_v_err;
        sum_tot_err+=tot_err; sum_tot_err2+=tot_err*tot_err;
        nErrPoints++;

      //----------------------------------------------------------------------------------------
        if(eventNum < 20 && totalFitHits < 30) {  // 只打印前几个事件的前几个命中点
          std::printf("Event %zu, Hit %zu: u_err:%.8f,v_err: %.8f\n",eventNum, totalFitHits, hit_u_err, hit_v_err);
        }
      //------------------------------------------------------------------------------------------

      }
    }
  }

  std::fprintf(stdout, "\n=== Kink Angle Range ===\n");
  std::fprintf(stdout, "Maximum kink angle: %.8f radians (event %zu)\n", max_kinkAngle, max_eventNum);
  std::fprintf(stdout, "Minimum kink angle: %.8f radians (event %zu)\n", min_kinkAngle, min_eventNum);

  
 // if(nErrPoints>0){
    double mean_u_err = sum_u_err/nErrPoints;
    double rms_u_err  = std::sqrt((sum_u_err2 - sum_u_err*sum_u_err/nErrPoints)/nErrPoints);
    double mean_v_err = sum_v_err/nErrPoints;
    double rms_v_err  = std::sqrt((sum_v_err2 - sum_v_err*sum_v_err/nErrPoints)/nErrPoints);
    double mean_tot_err = sum_tot_err/nErrPoints;
    double rms_tot_err  = std::sqrt((sum_tot_err2 - sum_tot_err*sum_tot_err/nErrPoints)/nErrPoints);
    std::fprintf(stdout, "\n=== Track fitting precision at DUT ===\n");
    std::fprintf(stdout, "Number of points used for error analysis: %zu\n", nErrPoints);
    std::fprintf(stdout, "U direction error: mean = %.6f um, RMS = %.6f um\n", mean_u_err, rms_u_err);
    std::fprintf(stdout, "V direction error: mean = %.6f um, RMS = %.6f um\n", mean_v_err, rms_v_err);
    std::fprintf(stdout, "Total position error: mean = %.6f um, RMS = %.6f um\n", mean_tot_err, rms_tot_err);
  //}

  TCanvas *c1 = new TCanvas("c1", "U VS V Error Analysis ", 1200, 1000);
  c1->Divide(1,2);
  c1->cd(1);
  hu_err_0->SetStats(1);
  gStyle->SetStatX(0.9);
  hu_err_0->SetLineColor(kBlue+2);
  hu_err_0->SetLineWidth(2);
  hu_err_0->Draw();
  TF1 *fitU0 = new TF1("fitU0", "gaus", 4.863, 4.864);
  fitU0->SetLineColor(kRed);
  hu_err_0->Fit(fitU0, "R");

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.06);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.13, 0.80, Form("Mean = %.6f", mean_u_err));
 // latex1.DrawLatex(0.13, 0.80, Form("Mean = %.6f", fitU0->GetParameter(1)));
 // latex1.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitU0->GetParameter(2)));

  c1->cd(2);
  hv_err_0->SetStats(1);
  hv_err_0->SetLineColor(kBlue+2);
  hv_err_0->SetLineWidth(2);
  hv_err_0->Draw();
  TF1 *fitV0 = new TF1("fitV1", "gaus", 4.863, 4.864);
  fitV0->SetLineColor(kRed);
  hv_err_0->Fit(fitV0, "R");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.06);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.13, 0.80, Form("Mean = %.6f", mean_v_err));
//  latex2.DrawLatex(0.13, 0.80, Form("Mean = %.6f", fitV0->GetParameter(1)));
 // latex2.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitV0->GetParameter(2)));


  TCanvas *c2 = new TCanvas("c2", "Origin Hits by Detector Analysis", 1800, 1200);

  c2->Divide(2, 2);
  c2->cd(1);
  hu_err->SetStats(1);
  hu_err->SetLineColor(kBlue+2);
  hu_err->SetLineWidth(2);
  hu_err->Draw();
  /*TF1 *fitU1 = new TF1("fitU4", "gaus", 4.8637,4.8638);
  fitU1->SetLineColor(kRed);
  hu_err->Fit(fitU1, "R");

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.06);
  latex3.SetTextColor(kRed);
  latex3.DrawLatex(0.13, 0.80, Form("Mean = %.6f", fitU1->GetParameter(1)));
  //latex3.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitU1->GetParameter(2)));
*/
  c2->cd(3);
  hv_err->SetStats(1);
  hv_err->SetLineColor(kBlue+2);
  hv_err->SetLineWidth(2);
  hv_err->Draw();
  /*TF1 *fitV1 = new TF1("fitV1", "gaus", 4.86368,4.86375);
  fitV1->SetLineColor(kRed);
  hv_err->Fit(fitV1, "R");  

  TLatex latex4;
  latex4.SetNDC();
  latex4.SetTextSize(0.06);
  latex4.SetTextColor(kRed);
  latex4.DrawLatex(0.4, 0.80, Form("Mean = %.6f", fitV1->GetParameter(1)));
 // latex4.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitV1->GetParameter(2)));
*/
   c2->cd(2);
   hu_err2->SetStats(1);
   hu_err2->SetLineColor(kBlue+2);
   hu_err2->SetLineWidth(2);
   hu_err2->Draw();
  /* TF1 *fitU2 = new TF1("fitU2", "gaus", 4.86364,4.8638);
   fitU2->SetLineColor(kRed);
   hu_err2->Fit(fitU2, "R");  

   TLatex latex5;
   latex5.SetNDC();
   latex5.SetTextSize(0.06);
   latex5.SetTextColor(kRed);
   latex5.DrawLatex(0.13, 0.80, Form("Mean = %.6f", fitU2->GetParameter(1)));
  // latex5.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitU2->GetParameter(2)));
*/
   c2->cd(4);
   hv_err2->SetStats(1);
   hv_err2->SetLineColor(kBlue+2);
   hv_err2->SetLineWidth(2);
   hv_err2->Draw();
  /* TF1 *fitV2 = new TF1("fitV2", "gaus", 4.86364,4.8638);
   fitV2->SetLineColor(kRed);
   hv_err2->Fit(fitV2, "R");  

   TLatex latex6;
   latex6.SetNDC();
   latex6.SetTextSize(0.06);
   latex6.SetTextColor(kRed);
   latex6.DrawLatex(0.13, 0.80, Form("Mean = %.6f", fitV2->GetParameter(1)));
  // latex6.DrawLatex(0.13, 0.70, Form("Sigma = %.6f", fitV2->GetParameter(2)));
*/

   TCanvas *c3 = new TCanvas("c3", "Kink Angle Distribution", 800, 600);   
  // hkink_angle->Scale(1.0 / hkink_angle->Integral());
   hkink_angle->Draw("hist");



  std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to origin_hits_by_detector.root\n");
 
  return;
}
