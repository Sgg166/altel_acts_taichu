void v_v_r(const std::string& rootFilePath){
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
 
  double region_u_min = -11.6;
  double region_u_max = 14;
  double region_v_min = -5.5;
  double region_v_max = 7.3;
 
  TH2F *cluster_u_Residual_det32 = new TH2F("cluster_u_Residual_det32", "cluster_u_Residual_det 32", 200, region_u_min, region_u_max, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det32 = new TH2F("cluster_v_Residual_det32", "cluster_v_Residual_det 32", 200, region_v_min, region_v_max, 200, -0.5, 0.5 );
 
  TH1F *hResidual_u =  new TH1F("hResidual_u", " Residual u;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v = new TH1F("hResidual_v", " Residual v;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
 
  TProfile *pResidualU_vs_clusterU = new TProfile("pResidualU_vs_clusterU", "Residual_u_mean_cluster_u", 50, region_u_min, region_u_max);
  TProfile *pResidualV_vs_clusterV = new TProfile("pResidualV_vs_clusterV", "Residual_v_mean_cluster_v", 50, region_v_min, region_v_max);
 
  TH2F *cluster_uv_Residual_sqrt_32 =new TH2F("cluster_uv_Residual_sqrt_32","cluster_uv_Residual_sqrt_32",200,region_u_min, region_u_max,200,region_v_min, region_v_max);
  size_t totalNumEvents = ttreeReader.numEvents();
  size_t goodEventNum = 0;
  size_t totalTracks = 0;
  size_t goodTracks = 0;
  size_t totalMatchedHits=0;
  size_t totalFitHits =0;
 
  // 新增统计变量
  size_t eventsWithOnlyGoodTrack = 0;        // 只有1条好轨迹的事件
  size_t eventsWithGoodAndBadTracks = 0;      // 1条好轨迹 + 多条坏轨迹的事件
  size_t eventsWithMultipleGoodTracks = 0;    // 多条好轨迹的事件
  size_t eventsWithMultipleGoodAndBad = 0;    // 多条好轨迹 + 多条坏轨迹的事件
  size_t eventsWithNoGoodTracks = 0;          // 没有好轨迹的事件
 
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    bool hasGoodTrack = false;
    size_t goodTracksInEvent = 0;      // 当前事件中的好轨迹数量
    size_t badTracksInEvent = 0;       // 当前事件中的坏轨迹数量
 
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() ==5 );
 
      if(isGoodTrack){
        goodTracks++;
        goodTracksInEvent++;
        hasGoodTrack = true;
 
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          auto aMatchedMeasHit=aTrajHit->matchedMeasHit();
          if(aFitHit && aMatchedMeasHit){
            totalFitHits++;
            totalMatchedHits++;
            int detID = aTrajHit->detN();
            if(detID == 32){
              double cluster_u = aMatchedMeasHit->u();
              double cluster_v = aMatchedMeasHit->v();
              double fit_u=aFitHit->u();
              double fit_v=aFitHit->v();
              double Residual_u=cluster_u-fit_u;
              double Residual_v=cluster_v-fit_v;
              double Residual_uv_sqrt=std::sqrt((Residual_u)*(Residual_u)+(Residual_v)*(Residual_v));
              bool isInRegion = (fit_u >= region_u_min && fit_u <= region_u_max &&
                                           fit_v >= region_v_min && fit_v <= region_v_max);
              if(isInRegion){
                cluster_u_Residual_det32->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det32->Fill(cluster_v,Residual_v);
                hResidual_u->Fill(Residual_u);
                hResidual_v->Fill(Residual_v);
                pResidualU_vs_clusterU->Fill(cluster_u, Residual_u);
                pResidualV_vs_clusterV->Fill(cluster_v, Residual_v);
                cluster_uv_Residual_sqrt_32->Fill(cluster_u,cluster_v,Residual_uv_sqrt);
              }
            }
          }
        }
      } else {
        badTracksInEvent++;  // 坏轨迹计数
      }
    }
 
    // 事件分类统计
    if(hasGoodTrack) {
      goodEventNum++;
      
      if(goodTracksInEvent == 1 && badTracksInEvent == 0) {
        eventsWithOnlyGoodTrack++;  // 只有1条好轨迹
      } else if(goodTracksInEvent == 1 && badTracksInEvent > 0) {
        eventsWithGoodAndBadTracks++;  // 1条好轨迹 + 多条坏轨迹
      } else if(goodTracksInEvent > 1 && badTracksInEvent == 0) {
        eventsWithMultipleGoodTracks++;  // 多条好轨迹
      } else if(goodTracksInEvent > 1 && badTracksInEvent > 0) {
        eventsWithMultipleGoodAndBad++;  // 多条好轨迹 + 多条坏轨迹
      }
    } else {
      eventsWithNoGoodTracks++;  // 没有好轨迹
    }
  }
 
  double efficiency = (goodTracks > 0) ? (double)totalMatchedHits / (double)goodTracks : 0.0;
  std::fprintf(stdout, "========================================\n");
  std::fprintf(stdout, "事件统计结果:\n");
  std::fprintf(stdout, "========================================\n");
  std::fprintf(stdout, "总事件数: %zu\n", totalNumEvents);
  std::fprintf(stdout, "包含好轨迹的事件数: %zu\n", goodEventNum);
  std::fprintf(stdout, "没有好轨迹的事件数: %zu\n", eventsWithNoGoodTracks);
  std::fprintf(stdout, "总轨迹数: %zu\n", totalTracks);
  std::fprintf(stdout, "好轨迹数 (==5 hits): %zu\n", goodTracks);
  std::fprintf(stdout, "坏轨迹数: %zu\n", totalTracks - goodTracks);
  std::fprintf(stdout, "----------------------------------------\n");
  std::fprintf(stdout, "包含好轨迹事件的详细分布:\n");
  std::fprintf(stdout, "只有1条好轨迹的事件: %zu\n", eventsWithOnlyGoodTrack);
  std::fprintf(stdout, "1条好轨迹 + 多条坏轨迹的事件: %zu\n", eventsWithGoodAndBadTracks);
  std::fprintf(stdout, "多条好轨迹的事件: %zu\n", eventsWithMultipleGoodTracks);
  std::fprintf(stdout, "多条好轨迹 + 多条坏轨迹的事件: %zu\n", eventsWithMultipleGoodAndBad);
  std::fprintf(stdout, "----------------------------------------\n");
  std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
  std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
  std::fprintf(stdout, "Efficiency: %.3f\n",efficiency);
  std::fprintf(stdout, "========================================\n");
 
  //double mean  = fitFunc->GetParameter(1);
  //double sigma = fitFunc->GetParameter(2);
  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 2400, 1600);
  c1->Divide(2, 4);
 
  c1->cd(1);
  cluster_u_Residual_det32->SetFillColor(kGreen);
  gStyle->SetStatX(0.9);
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
//  gStyle->SetOptFit(1111);
 
  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.07);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU->GetParameter(1)));
  latex1.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU->GetParameter(2)));
 
 
  c1->cd(4);
  hResidual_v->SetLineColor(kGreen+2);
  hResidual_v->SetLineWidth(2);
  hResidual_v->Draw();
 
  TF1 *fitV = new TF1("fitV", "gaus", -0.05, 0.05);
  fitV->SetLineColor(kRed);
  hResidual_v->Fit(fitV, "R");
 
  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.07);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV->GetParameter(1)));
  latex2.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV->GetParameter(2)));
 
 
  c1->cd(5);
  pResidualU_vs_clusterU->SetMarkerColor(kRed);
  pResidualU_vs_clusterU->GetXaxis()->SetTitle("cluster_u [mm]");
  pResidualU_vs_clusterU->GetYaxis()->SetTitle("<Residual_u> [mm]");
  pResidualU_vs_clusterU->SetMarkerStyle(20);
  pResidualU_vs_clusterU->SetMarkerSize(0.8);
  pResidualU_vs_clusterU->SetLineColor(kRed);
  pResidualU_vs_clusterU->SetLineWidth(2);
  pResidualU_vs_clusterU->GetYaxis()->SetRangeUser(-0.01,0.01);
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
  pResidualV_vs_clusterV->GetYaxis()->SetRangeUser(-0.01,0.01);
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
