void residual(const std::string& rootFilePath){
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
 
  std::map<std::string, int> detCombinationCount; 
  
  TH1F *hResidual_u = new TH1F("hResidual_x", "Residual x;x_{meas}-x_{fit} [um];Entries", 200, -100, 100);
  TH1F *hResidual_v = new TH1F("hResidual_y", "Residual y;y_{meas}-y_{fit} [um];Entries", 200, -100, 100);
 
  TH2F *cluster_uv_Residual = new TH2F("cluster_uv_Residual","cluster_uv_Residual",100,-40, 40,100,-40,40);
 
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t outEventNum =0;
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==5);
      if(!isGoodTrack){
        continue;
      }
 
      std::vector<int> usedDetIds;
      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        if(aFitHit){
          usedDetIds.push_back(aFitHit->detN());
        }
      }
 
      std::sort(usedDetIds.begin(), usedDetIds.end());
      std::string detCombination;
      for(size_t i = 0; i < usedDetIds.size(); i++){
        detCombination += std::to_string(usedDetIds[i]);
        if(i < usedDetIds.size() - 1){
          detCombination += "-";
        }
      }
 
      if(usedDetIds.size() >= 3){
        detCombinationCount[detCombination]++;
        std::cout << "Event " << eventNum << ": Detectors [" << detCombination << "]" << std::endl;
      }
 
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
        outEventNum++;
 
        if(aMatchedMeasHit){
          double cluster_u = aMatchedMeasHit->u();
          double cluster_v = aMatchedMeasHit->v();
          double Residual_u=(cluster_u - fit_u)*1000.0;
          double Residual_v=(cluster_v - fit_v)*1000.0;
          hResidual_u->Fill(Residual_u);
          hResidual_v->Fill(Residual_v);
          cluster_uv_Residual->Fill(Residual_u,Residual_v);
        }
      }
    }
  }
 
  std::cout << "\n========================================" << std::endl;
  std::cout << "Detector Combination Statistics:" << std::endl;
  std::cout << "========================================" << std::endl;
  for(auto &pair: detCombinationCount){
    std::cout << "Detectors [" << pair.first << "]: " << pair.second << " tracks" << std::endl;
  }
  std::cout << "========================================" << std::endl;
 
  std::fprintf(stdout,"found %zu outEventNum\n",outEventNum);

  TCanvas *c1 = new TCanvas("c1", "Residuals", 2000, 1000);
  c1->Divide(2,1);

  c1->cd(1);
  hResidual_u->SetLineColor(kBlue+2);
  hResidual_u->SetLineWidth(2);
  hResidual_u->SetMarkerStyle(20);
  hResidual_u->Draw();

  TF1 *fitU = new TF1("fitU", "gaus", -30, 30);
  fitU->SetLineColor(kRed);
  hResidual_u->Fit(fitU, "R");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  latex.SetTextColor(kRed);
  latex.DrawLatex(0.7, 0.65, Form("Mean = %.2f", fitU->GetParameter(1)));
  latex.DrawLatex(0.7, 0.60, Form("Sigma = %.2f", fitU->GetParameter(2)));

  c1->cd(2);
  hResidual_v->SetLineColor(kGreen+2);
  hResidual_v->SetLineWidth(2);
  hResidual_v->Draw();

  TF1 *fitV = new TF1("fitV", "gaus", -30, 30);
  fitV->SetLineColor(kRed);
  hResidual_v->Fit(fitV, "R");

  latex.SetTextColor(kRed);
  latex.DrawLatex(0.7, 0.65, Form("Mean = %.2f", fitV->GetParameter(1)));
  latex.DrawLatex(0.7, 0.60, Form("Sigma = %.2f", fitV->GetParameter(2)));

  TCanvas *c2 = new TCanvas("c2", "Residuals", 1000, 1000);
  cluster_uv_Residual->SetStats(0);
  cluster_uv_Residual->GetXaxis()->SetTitle("Residual_U [um]");
  cluster_uv_Residual->GetYaxis()->SetTitle("Residual_V [um]");
  cluster_uv_Residual->Draw("COLZ");

//  c1->SaveAs("residuals_1D_fit.png");
 // c1->SaveAs("residuals_1D_fit.pdf");

  TFile *outFile = new TFile("/home/sungaige/A/altel_acts/Draw/root_file/residuals_u_v.root", "RECREATE");
  hResidual_u->Write();
  hResidual_v->Write();
  cluster_uv_Residual->Write();
  c1->Write();
  c2->Write();
//  outFile->Close();

  std::fprintf(stdout, "Analysis complete. Results saved to /home/sungaige/A/altel_acts/Draw/root_file/residuals_u_v.root\n");
}
