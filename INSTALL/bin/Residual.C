void Residual(const std::string& rootFilePath){
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

  TH1F *hResidual_u = new TH1F("hResidual_u", "Residual u;u_{fit}-u_{meas} [mm];Entries", 200, -0.1, 0.1);
  TH1F *hResidual_v = new TH1F("hResidual_v", "Residual v;v_{fit}-v_{meas} [mm];Entries", 200, -0.1, 0.1);

  size_t totalNumEvents = ttreeReader.numEvents();

  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);

    for(auto aTraj: telEvent->trajs()){
      bool isGoodTrack = (aTraj->numOriginMeasHit() >= 3);
      if(!isGoodTrack) continue;

      for(auto &aTrajHit: aTraj->trajHits()){
        auto aFitHit = aTrajHit->fitHit();
        auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
        if(aFitHit && aMatchedMeasHit){
          int detID = aTrajHit->detN();
          if(detID == 32){
            double cluster_u = aMatchedMeasHit->u();
            double cluster_v = aMatchedMeasHit->v();
            double fit_u=aFitHit->u();
            double fit_v=aFitHit->v();
            double Residual_u=fit_u - cluster_u; // consistent with x_fit - x_meas
            double Residual_v=fit_v - cluster_v;

            hResidual_u->Fill(Residual_u);
            hResidual_v->Fill(Residual_v);
          }
        }
      }
    }
  }

  TCanvas *c1 = new TCanvas("c1", "Residuals", 1600, 600);
  c1->Divide(2,1);

  c1->cd(1);
  hResidual_u->SetLineColor(kBlue+2);
  hResidual_u->SetLineWidth(2);
  hResidual_u->Draw();

  TF1 *fitU = new TF1("fitU", "gaus", -0.05, 0.05);
  fitU->SetLineColor(kRed);
  hResidual_u->Fit(fitU, "R");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.SetTextColor(kRed);
  latex.DrawLatex(0.55, 0.85, Form("Mean = %.4f", fitU->GetParameter(1)));
  latex.DrawLatex(0.55, 0.80, Form("Sigma = %.4f", fitU->GetParameter(2)));

  c1->cd(2);
  hResidual_v->SetLineColor(kGreen+2);
  hResidual_v->SetLineWidth(2);
  hResidual_v->Draw();

  TF1 *fitV = new TF1("fitV", "gaus", -0.05, 0.05);
  fitV->SetLineColor(kRed);
  hResidual_v->Fit(fitV, "R");

  latex.SetTextColor(kRed);
  latex.DrawLatex(0.55, 0.85, Form("Mean = %.4f", fitV->GetParameter(1)));
  latex.DrawLatex(0.55, 0.80, Form("Sigma = %.4f", fitV->GetParameter(2)));

  c1->SaveAs("residuals_1D_fit.png");
  c1->SaveAs("residuals_1D_fit.pdf");

  TFile *outFile = new TFile("residuals_1D.root", "RECREATE");
  hResidual_u->Write();
  hResidual_v->Write();
  c1->Write();
//  outFile->Close();

  std::fprintf(stdout, "Analysis complete. Results saved to residuals_1D.root\n");
}

