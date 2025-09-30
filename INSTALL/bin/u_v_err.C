void u_v_err(const std::string& rootFilePath){
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
 

  
  TH1F* hu_err =new TH1F("hu_err","u error;u_err[um];Entries",100, 3.5,5.5);
  TH1F* hv_err =new TH1F("hv_err","v error;v_err[um];Entries",100, 3.5,5.5);

  size_t totalNumEvents = ttreeReader.numEvents(); 
  size_t totalFitHits = 0;
  double sum_u_err=0;
  double sum_u_err2=0;
  double sum_v_err=0;
  double sum_v_err2=0;
  double sum_tot_err=0;
  double sum_tot_err2=0;
  size_t nErrPoints=0;

  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);
      if(!isGoodTrack){
         continue;
      }
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
        if(!aFitHit ){
          continue;
        }
    /*    double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);

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
      */  
        if (aFitHit->detN() != 32){
          continue;
        }
        double fit_u=aFitHit->u();
        double fit_v=aFitHit->v();

       /* double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);        
        double tot_err = std::sqrt(hit_u_err*hit_u_err + hit_v_err*hit_v_err);
        sum_u_err+=hit_u_err; sum_u_err2+=hit_u_err*hit_u_err;
        sum_v_err+=hit_v_err; sum_v_err2+=hit_v_err*hit_v_err;
        sum_tot_err+=tot_err; sum_tot_err2+=tot_err*tot_err;
        nErrPoints++;*/
        if(fit_u < -12.8 || fit_u > 12.8 || fit_v < -6.4 || fit_v > 6.4 ){
          continue;
        }

        double hit_u_err =aFitHit->u_err()*1000.0;
        double hit_v_err =aFitHit->v_err()*1000.0;
        hu_err->Fill(hit_u_err);
        hv_err->Fill(hit_v_err);        
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
  if(nErrPoints>0){
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
  }

  TCanvas *c2 = new TCanvas("c2", "Origin Hits by Detector Analysis", 1800, 1200);

  c2->Divide(2, 1);
  c2->cd(1);
  hu_err->SetLineColor(kBlue+2);
  hu_err->SetLineWidth(2);
  hu_err->Draw();
  /*TF1 *fitU4 = new TF1("fitU4", "gaus", 2, 8);
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
  /*TF1 *fitV4 = new TF1("fitV4", "gaus", 0, 10);
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
