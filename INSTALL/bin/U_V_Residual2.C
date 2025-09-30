void U_V_Residual(const std::string& rootFilePath){
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
 
  TH2F *cluster_u_Residual_det32 = new TH2F("cluster_u_Residual_det32", "cluster_u_Residual_det 32", 200, -12.8,12.8, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det32 = new TH2F("cluster_v_Residual_det32", "cluster_v_Residual_det 32", 200, -6.4,6.4, 200, -0.5, 0.5 );
  TH1F *hResidual_u = new TH1F("hResidual_u", "Residual u;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v = new TH1F("hResidual_v", "Residual v;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  
  TProfile *pResidualU_vs_clusterU = new TProfile("pResidualU_vs_clusterU", "Residual_u_mean_cluster_u", 50, -12.8, 12.8);
  TProfile *pResidualV_vs_clusterV = new TProfile("pResidualV_vs_clusterV", "Residual_v_mean_cluster_v", 50, -6.4, 6.4);

  TH2F *cluster_uv_Residual_sqrt_32 =new TH2F("cluster_uv_Residual_sqrt_32","cluster_uv_Residual_sqrt_32",200,-12.8,12.8,200,-6.4,6.4);
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits=0;
  size_t totalFitHits =0;
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==5 );

 
      if(isGoodTrack){
        goodTracks++;
 
 
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
          if(aFitHit && aMatchedMeasHit){
            totalFitHits++;
            totalMatchedHits++;
           // int detID = aTrajHit->detN();
           // if(detID == 32){
            double cluster_u = aMatchedMeasHit->u();
            double cluster_v = aMatchedMeasHit->v();
            double fit_u=aFitHit->u();
            double fit_v=aFitHit->v();
            double Residual_u=cluster_u-fit_u;
            double Residual_v=cluster_v-fit_v;
           // double Residual_uv_sqrt=std::sqrt((Residual_u)*(Residual_u)+(Residual_v)*(Residual_v));
       
            int detID = aTrajHit->detN();
            if(detID == 32){
              cluster_u_Residual_det32->Fill(cluster_u,Residual_u);
              cluster_v_Residual_det32->Fill(cluster_v,Residual_v);
              hResidual_u->Fill(Residual_u);
              hResidual_v->Fill(Residual_v);
              pResidualU_vs_clusterU->Fill(cluster_u, Residual_u);
              pResidualV_vs_clusterV->Fill(cluster_v, Residual_v);
              cluster_uv_Residual_sqrt_32->Fill(cluster_u,cluster_v,Residual_u);
            }
          }
        }
      }
    }
  }
  double efficiency = (goodTracks > 0) ? (double)totalMatchedHits / (double)goodTracks : 0.0;
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good  tracks (==5 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  std::fprintf(stdout, "Efficiency: %.3f\n",efficiency); 
  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 2400, 1600);
  c1->Divide(2, 4);  
 
  c1->cd(1);
  cluster_u_Residual_det32->SetFillColor(kGreen);
  cluster_u_Residual_det32->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det32->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det32->SetMarkerColor(9);
  cluster_u_Residual_det32->SetMarkerStyle(20);
  cluster_u_Residual_det32->SetMarkerSize(0.8);
  cluster_u_Residual_det32->Draw("P");
 
  c1->cd(2);
  cluster_v_Residual_det32->SetFillColor(kGreen);
  cluster_v_Residual_det32->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det32->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det32->SetMarkerColor(4);
  cluster_v_Residual_det32->SetMarkerStyle(20);
  cluster_v_Residual_det32->SetMarkerSize(0.8);
  cluster_v_Residual_det32->Draw("P");
 
 
  c1->cd(3);
  hResidual_u->SetLineColor(kBlue+2);
  hResidual_u->SetLineWidth(2);
  hResidual_u->Draw();
 
  TF1 *fitU = new TF1("fitU", "gaus", -0.05, 0.05);
  fitU->SetLineColor(kRed);
  hResidual_u->Fit(fitU, "R");
  //gStyle->SetOptFit(1111);
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.SetTextColor(kRed);
  latex.DrawLatex(0.65, 0.85, Form("Mean = %.4f", mean));
  latex.DrawLatex(0.65, 0.80, Form("Sigma = %.4f", sigma));

  c1->cd(4);
  hResidual_v->SetLineColor(kGreen+2);
  hResidual_v->SetLineWidth(2);
  hResidual_v->Draw();
 
  TF1 *fitV = new TF1("fitV", "gaus", -0.05, 0.05);
  fitV->SetLineColor(kRed);
  hResidual_v->Fit(fitV, "R");
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);
  latex.SetTextColor(kRed);
  latex.DrawLatex(0.65, 0.85, Form("Mean = %.4f", mean));
  latex.DrawLatex(0.65, 0.80, Form("Sigma = %.4f", sigma));

  c1->cd(5);
  pResidualU_vs_clusterU->SetMarkerColor(kRed);
  pResidualU_vs_clusterU->GetXaxis()->SetTitle("cluster_u [mm]");
  pResidualU_vs_clusterU->GetYaxis()->SetTitle("<Residual_u> [mm]");
  pResidualU_vs_clusterU->SetMarkerStyle(20);
  pResidualU_vs_clusterU->SetMarkerSize(0.8);
  pResidualU_vs_clusterU->SetLineColor(kRed);
  pResidualU_vs_clusterU->SetLineWidth(2);
  pResidualU_vs_clusterU->Draw("PE");
  
//   添加一条y=0的参考线
  TLine *line1 = new TLine(-12.8, 0, 12.8, 0);
  line1->SetLineColor(kBlack);
  line1->SetLineStyle(2);
  line1->Draw("same");
   
  c1->cd(6);
  pResidualV_vs_clusterV->SetMarkerColor(kRed);
  pResidualV_vs_clusterV->GetXaxis()->SetTitle("cluster_v [mm]");
  pResidualV_vs_clusterV->GetYaxis()->SetTitle("<Residual_v> [mm]");
  pResidualV_vs_clusterV->SetMarkerStyle(20);
  pResidualV_vs_clusterV->SetMarkerSize(0.8);
  pResidualV_vs_clusterV->SetLineColor(8);
  pResidualV_vs_clusterV->SetLineWidth(2);
  pResidualV_vs_clusterV->Draw("PE");
  TLine *line2 = new TLine(-6.4, 0, 6.4, 0);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(2);
  line2->Draw("same");


  c1->cd(7);
  cluster_uv_Residual_sqrt_32->SetStats(0);
  cluster_uv_Residual_sqrt_32->GetXaxis()->SetTitle("cluster_u [mm]");
  cluster_uv_Residual_sqrt_32->GetYaxis()->SetTitle("cluster_v [mm]");
  cluster_uv_Residual_sqrt_32->GetZaxis()->SetTitle("sqrt{Residual_u^2 + Residual_v^2} [mm]");
  cluster_uv_Residual_sqrt_32->SetMinimum(-0.05); 
  cluster_uv_Residual_sqrt_32->SetMaximum(0.05);
  //cluster_uv_Residual_sqrt_32->Draw("LEGO");
  cluster_uv_Residual_sqrt_32->Draw("COLZ");

  c1->SaveAs("residuals_det32.svg");
 
  TFile *outFile = new TFile("residuals_det32.root", "RECREATE");
  cluster_u_Residual_det32->Write();
  cluster_v_Residual_det32->Write();
  hResidual_u->Write();
  hResidual_v->Write();
  pResidualU_vs_clusterU->Write();  
  pResidualV_vs_clusterV->Write();
  cluster_uv_Residual_sqrt_32->Write();
  c1->Write();
  std::fprintf(stdout, "Analysis complete. Results saved to residuals_det32.root\n");
 // std::getc(stdin);
 
  return;
}
