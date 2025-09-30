void after_Hit_Reader(const std::string& rootFilePath){
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
 

  TH1F *hNumOriginHits_det2 = new TH1F("hNumOriginHits_det2", "good_fit_origin_Claster_size_det 2", 30, 0, 15);
  TH1F *hNumOriginHits_det3 = new TH1F("hNumOriginHits_det3", "good_fit_origin_Claster_size_det 3", 30, 0, 15);
  TH1F *hNumOriginHits_det32 = new TH1F("hNumOriginHits_det32", "good_fit_origin_Claster_size_det 32", 30, 0, 15);
  TH1F *hNumOriginHits_det5 = new TH1F("hNumOriginHits_det5", "good_fit_origin_Claster_size_det 5", 30, 0, 15);
  TH1F *hNumOriginHits_det7 = new TH1F("hNumOriginHits_det7", "good_fit_origin_Claster_size_det 7", 30, 0, 15);
  TH1F *hNumOriginHits_det9 = new TH1F("hNumOriginHits_det9", "good_fit_origin_Claster_size_det 9", 30, 0, 15);
 

  TH2F *cluster_u_Residual_det2 = new TH2F("cluster_u_Residual_det2", "cluster_u_Residual_det 2", 200, -12.8,12.8, 200, -0.5, 0.5 );
  TH2F *cluster_u_Residual_det3 = new TH2F("cluster_u_Residual_det3", "cluster_u_Residual_det 3", 200, -12.8,12.8, 200, -0.5, 0.5 );
  TH2F *cluster_u_Residual_det5 = new TH2F("cluster_u_Residual_det5", "cluster_u_Residual_det 5", 200, -12.8,12.8, 200, -0.5, 0.5 );
  TH2F *cluster_u_Residual_det7 = new TH2F("cluster_u_Residual_det7", "cluster_u_Residual_det 7", 200, -12.8,12.8, 200, -0.5, 0.5 );
  TH2F *cluster_u_Residual_det9 = new TH2F("cluster_u_Residual_det9", "cluster_u_Residual_det 9", 200, -12.8,12.8, 200, -0.5, 0.5 );

  TH2F *cluster_v_Residual_det2 = new TH2F("cluster_v_Residual_det2", "cluster_v_Residual_det 2", 200, -6.4,6.4, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det3 = new TH2F("cluster_v_Residual_det3", "cluster_v_Residual_det 3", 200, -6.4,6.4, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det5 = new TH2F("cluster_v_Residual_det5", "cluster_v_Residual_det 5", 200, -6.4,6.4, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det7 = new TH2F("cluster_v_Residual_det7", "cluster_v_Residual_det 7", 200, -6.4,6.4, 200, -0.5, 0.5 );
  TH2F *cluster_v_Residual_det9 = new TH2F("cluster_v_Residual_det9", "cluster_v_Residual_det 9", 200, -6.4,6.4, 200, -0.5, 0.5 );

  TH1F *hResidual_u_det2 = new TH1F("hResidual_u_det2", "Residual u det2;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_u_det3 = new TH1F("hResidual_u_det3", "Residual u det3;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_u_det5 = new TH1F("hResidual_u_det5", "Residual u det5;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_u_det7 = new TH1F("hResidual_u_det7", "Residual u det7;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_u_det9 = new TH1F("hResidual_u_det9", "Residual u det9;u_{meas}-u_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);

  TH1F *hResidual_v_det2 = new TH1F("hResidual_v_det2", "Residual v det2;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v_det3 = new TH1F("hResidual_v_det3", "Residual v det3;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v_det5 = new TH1F("hResidual_v_det5", "Residual v det5;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v_det7 = new TH1F("hResidual_v_det7", "Residual v det7;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);
  TH1F *hResidual_v_det9 = new TH1F("hResidual_v_det9", "Residual v det9;v_{meas}-v_{fit} [mm];Entries [200bin]", 200, -0.1, 0.1);




  size_t totalNumEvents = ttreeReader.numEvents();
  size_t totalTracks = 0;
  size_t goodTracks = 0;
 
  for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
    for(auto aTraj: telEvent->trajs()){
      totalTracks++;
      bool isGoodTrack = (aTraj->numOriginMeasHit() >= 3); 
 
      if(isGoodTrack){
        goodTracks++;
        
        
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          if(aFitHit){
            auto aOriginMeasHit = aFitHit->originMeasHit();
            if(aOriginMeasHit){
              size_t numOrigin = aOriginMeasHit->measRaws().size();
              double cluster_u = aOriginMeasHit->u();
              double cluster_v = aOriginMeasHit->v();
              double fit_u=aFitHit->u();
              double fit_v=aFitHit->v();
              double Residual_u=cluster_u-fit_u;
              double Residual_v=cluster_v-fit_v;
              int detID = aTrajHit->detN();
              
              if(detID == 2){
                hNumOriginHits_det2->Fill(numOrigin);
                cluster_u_Residual_det2->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det2->Fill(cluster_v,Residual_v);
                hResidual_u_det2->Fill(Residual_u);
                hResidual_v_det2->Fill(Residual_v);
              }
              else if(detID == 3){
                hNumOriginHits_det3->Fill(numOrigin);
                cluster_u_Residual_det3->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det3->Fill(cluster_v,Residual_v);
                hResidual_u_det3->Fill(Residual_u);
                hResidual_v_det3->Fill(Residual_v);
              }
    /*          else if(detID == 32){
                hNumOriginHits_det32->Fill(numOrigin);
              }*/
              else if(detID == 5){
                hNumOriginHits_det5->Fill(numOrigin);
                cluster_u_Residual_det5->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det5->Fill(cluster_v,Residual_v);
                hResidual_u_det5->Fill(Residual_u);
                hResidual_v_det5->Fill(Residual_v);
              }
              else if(detID == 7){
                hNumOriginHits_det7->Fill(numOrigin);
                cluster_u_Residual_det7->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det7->Fill(cluster_v,Residual_v);
                hResidual_u_det7->Fill(Residual_u);
                hResidual_v_det7->Fill(Residual_v);
              }
              else if(detID == 9){
                hNumOriginHits_det9->Fill(numOrigin);
                cluster_u_Residual_det9->Fill(cluster_u,Residual_u);
                cluster_v_Residual_det9->Fill(cluster_v,Residual_v);
                hResidual_u_det9->Fill(Residual_u);
                hResidual_v_det9->Fill(Residual_v);
              }
            }
          }
        }
      }
    }
  }
 
  std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
  std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
  std::fprintf(stdout, "Good tracks (>3 hits): %zu\n", goodTracks);
 
  TCanvas *c1 = new TCanvas("c1", "Origin Hits by Detector Analysis", 1800, 1200);
  c1->Divide(5, 5);
 
  c1->cd(1);
  hNumOriginHits_det2->SetFillColor(kRed);
  hNumOriginHits_det2->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  if(hNumOriginHits_det2->GetEntries() > 0){
    hNumOriginHits_det2->Scale(1.0 / hNumOriginHits_det2->GetEntries());
  }
  hNumOriginHits_det2->GetYaxis()->SetTitle("Percentage");
  hNumOriginHits_det2->SetFillStyle(3001);
  hNumOriginHits_det2->SetMinimum(0);
  hNumOriginHits_det2->SetMaximum(1.0);
  hNumOriginHits_det2->Draw("HIST");


  c1->cd(2);
  cluster_u_Residual_det2->SetFillColor(kGreen);
  cluster_u_Residual_det2->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det2->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det2->SetMarkerColor(9);
  cluster_u_Residual_det2->SetMarkerStyle(20);
  cluster_u_Residual_det2->SetMarkerSize(0.8);
  cluster_u_Residual_det2->Draw("P");

  c1->cd(3);
  cluster_v_Residual_det2->SetFillColor(kGreen);
  cluster_v_Residual_det2->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det2->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det2->SetMarkerColor(4);
  cluster_v_Residual_det2->SetMarkerStyle(20);
  cluster_v_Residual_det2->SetMarkerSize(0.8);
  cluster_v_Residual_det2->Draw("P");


  c1->cd(4);
  hResidual_u_det2->SetLineColor(kBlue+2);
  hResidual_u_det2->SetLineWidth(2);
  hResidual_u_det2->Draw();

  TF1 *fitU2 = new TF1("fitU2", "gaus", -0.05, 0.05);
  fitU2->SetLineColor(kRed);
  hResidual_u_det2->Fit(fitU2, "R");
//  gStyle->SetOptFit(1111);

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextSize(0.07);
  latex1.SetTextColor(kRed);
  latex1.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU2->GetParameter(1)));
  latex1.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU2->GetParameter(2)));


  c1->cd(5);
  hResidual_v_det2->SetLineColor(kGreen+2);
  hResidual_v_det2->SetLineWidth(2);
  hResidual_v_det2->Draw();

  TF1 *fitV2 = new TF1("fitV2", "gaus", -0.05, 0.05);
  fitV2->SetLineColor(kRed);
  hResidual_v_det2->Fit(fitV2, "R");

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextSize(0.07);
  latex2.SetTextColor(kRed);
  latex2.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV2->GetParameter(1)));
  latex2.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV2->GetParameter(2)));
 
  c1->cd(6);
  hNumOriginHits_det3->SetFillColor(kBlue);
  hNumOriginHits_det3->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  hNumOriginHits_det3->GetYaxis()->SetTitle("Percentage");
  if(hNumOriginHits_det3->GetEntries() > 0){
    hNumOriginHits_det3->Scale(1.0 / hNumOriginHits_det3->GetEntries());
  }
  hNumOriginHits_det3->SetMinimum(0);
  hNumOriginHits_det3->SetMaximum(1.0);
  hNumOriginHits_det3->SetFillStyle(3001);
  hNumOriginHits_det3->Draw("HIST");
 
  c1->cd(7);
  cluster_u_Residual_det3->SetFillColor(kGreen);
  cluster_u_Residual_det3->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det3->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det3->SetMarkerColor(9);
  cluster_u_Residual_det3->SetMarkerStyle(20);
  cluster_u_Residual_det3->SetMarkerSize(0.8);
  cluster_u_Residual_det3->Draw("P");

  c1->cd(8);
  cluster_v_Residual_det3->SetFillColor(kGreen);
  cluster_v_Residual_det3->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det3->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det3->SetMarkerColor(4);
  cluster_v_Residual_det3->SetMarkerStyle(20);
  cluster_v_Residual_det3->SetMarkerSize(0.8);
  cluster_v_Residual_det3->Draw("P");
  

  c1->cd(9);
  hResidual_u_det3->SetLineColor(kBlue+2);
  hResidual_u_det3->SetLineWidth(2);
  hResidual_u_det3->Draw();

  TF1 *fitU3 = new TF1("fitU3", "gaus", -0.05, 0.05);
  fitU3->SetLineColor(kRed);
  hResidual_u_det3->Fit(fitU3, "R");
//  gStyle->SetOptFit(1111);

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextSize(0.07);
  latex3.SetTextColor(kRed);
  latex3.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU3->GetParameter(1)));
  latex3.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU3->GetParameter(2)));


  c1->cd(10);
  hResidual_v_det3->SetLineColor(kGreen+2);
  hResidual_v_det3->SetLineWidth(2);
  hResidual_v_det3->Draw();

  TF1 *fitV3 = new TF1("fitV3", "gaus", -0.05, 0.05);
  fitV3->SetLineColor(kRed);
  hResidual_v_det3->Fit(fitV3, "R");

  TLatex latex4;
  latex4.SetNDC();
  latex4.SetTextSize(0.07);
  latex4.SetTextColor(kRed);
  latex4.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV3->GetParameter(1)));
  latex4.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV3->GetParameter(2)));

/*  c1->cd(7);
  hNumOriginHits_det32->SetFillColor(kGreen);
  if(hNumOriginHits_det32->GetEntries() > 0){
    hNumOriginHits_det32->Scale(1.0 / hNumOriginHits_det32->GetEntries());
  }
  hNumOriginHits_det32->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  hNumOriginHits_det32->GetYaxis()->SetTitle("Percentage");
  hNumOriginHits_det32->SetMinimum(0);
  hNumOriginHits_det32->SetMaximum(1.0);
  hNumOriginHits_det32->SetFillStyle(3001);
  hNumOriginHits_det32->Draw("HIST");
*/

  c1->cd(11);
  hNumOriginHits_det5->SetFillColor(kMagenta);
  if(hNumOriginHits_det5->GetEntries() > 0){
    hNumOriginHits_det5->Scale(1.0 / hNumOriginHits_det5->GetEntries());
  }
  hNumOriginHits_det5->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  hNumOriginHits_det5->GetYaxis()->SetTitle("Percentage");
  hNumOriginHits_det5->SetMinimum(0);
  hNumOriginHits_det5->SetMaximum(1.0);
  hNumOriginHits_det5->SetFillStyle(3001);
  hNumOriginHits_det5->Draw("HIST");


  c1->cd(12);
  cluster_u_Residual_det5->SetFillColor(kGreen);
  cluster_u_Residual_det5->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det5->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det5->SetMarkerColor(9);
  cluster_u_Residual_det5->SetMarkerStyle(20);
  cluster_u_Residual_det5->SetMarkerSize(0.8);
  cluster_u_Residual_det5->Draw("P");

  c1->cd(13);
  cluster_v_Residual_det5->SetFillColor(kGreen);
  cluster_v_Residual_det5->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det5->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det5->SetMarkerColor(4);
  cluster_v_Residual_det5->SetMarkerStyle(20);
  cluster_v_Residual_det5->SetMarkerSize(0.8);
  cluster_v_Residual_det5->Draw("P");
 
  
  c1->cd(14);
  hResidual_u_det5->SetLineColor(kBlue+2);
  hResidual_u_det5->SetLineWidth(2);
  hResidual_u_det5->Draw();

  TF1 *fitU5 = new TF1("fitU5", "gaus", -0.05, 0.05);
  fitU5->SetLineColor(kRed);
  hResidual_u_det5->Fit(fitU5, "R");
//  gStyle->SetOptFit(1111);

  TLatex latex5;
  latex5.SetNDC();
  latex5.SetTextSize(0.07);
  latex5.SetTextColor(kRed);
  latex5.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU5->GetParameter(1)));
  latex5.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU5->GetParameter(2)));


  c1->cd(15);
  hResidual_v_det5->SetLineColor(kGreen+2);
  hResidual_v_det5->SetLineWidth(2);
  hResidual_v_det5->Draw();

  TF1 *fitV5 = new TF1("fitV5", "gaus", -0.05, 0.05);
  fitV5->SetLineColor(kRed);
  hResidual_v_det5->Fit(fitV5, "R");

  TLatex latex6;
  latex6.SetNDC();
  latex6.SetTextSize(0.07);
  latex6.SetTextColor(kRed);
  latex6.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV5->GetParameter(1)));
  latex6.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV5->GetParameter(2)));

  

  c1->cd(16);
  hNumOriginHits_det7->SetFillColor(kOrange);
  if(hNumOriginHits_det7->GetEntries() > 0){
    hNumOriginHits_det7->Scale(1.0 / hNumOriginHits_det7->GetEntries());
  }
  hNumOriginHits_det7->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  hNumOriginHits_det7->GetYaxis()->SetTitle("Percentage");
  hNumOriginHits_det7->SetMinimum(0);
  hNumOriginHits_det7->SetMaximum(1.0);
  hNumOriginHits_det7->SetFillStyle(3001);
  hNumOriginHits_det7->Draw("HIST");
 
  c1->cd(17);
  cluster_u_Residual_det7->SetFillColor(kGreen);
  cluster_u_Residual_det7->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det7->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det7->SetMarkerColor(9);
  cluster_u_Residual_det7->SetMarkerStyle(20);
  cluster_u_Residual_det7->SetMarkerSize(0.8);
  cluster_u_Residual_det7->Draw("P");

  c1->cd(18);
  cluster_v_Residual_det7->SetFillColor(kGreen);
  cluster_v_Residual_det7->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det7->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det7->SetMarkerColor(4);
  cluster_v_Residual_det7->SetMarkerStyle(20);
  cluster_v_Residual_det7->SetMarkerSize(0.8);
  cluster_v_Residual_det7->Draw("P");

  
  c1->cd(19);
  hResidual_u_det7->SetLineColor(kBlue+2);
  hResidual_u_det7->SetLineWidth(2);
  hResidual_u_det7->Draw();

  TF1 *fitU7 = new TF1("fitU7", "gaus", -0.05, 0.05);
  fitU7->SetLineColor(kRed);
  hResidual_u_det7->Fit(fitU7, "R");
//  gStyle->SetOptFit(1111);

  TLatex latex7;
  latex7.SetNDC();
  latex7.SetTextSize(0.07);
  latex7.SetTextColor(kRed);
  latex7.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU7->GetParameter(1)));
  latex7.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU7->GetParameter(2)));


  c1->cd(20);
  hResidual_v_det7->SetLineColor(kGreen+2);
  hResidual_v_det7->SetLineWidth(2);
  hResidual_v_det7->Draw();

  TF1 *fitV7 = new TF1("fitV7", "gaus", -0.05, 0.05);
  fitV7->SetLineColor(kRed);
  hResidual_v_det7->Fit(fitV7, "R");

  TLatex latex8;
  latex8.SetNDC();
  latex8.SetTextSize(0.07);
  latex8.SetTextColor(kRed);
  latex8.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV7->GetParameter(1)));
  latex8.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV7->GetParameter(2)));

  c1->cd(21);
  hNumOriginHits_det9->SetFillColor(kViolet);
  if(hNumOriginHits_det9->GetEntries() > 0){
    hNumOriginHits_det9->Scale(1.0 / hNumOriginHits_det9->GetEntries());
  }
  hNumOriginHits_det9->GetXaxis()->SetTitle("good_fit_origin_Claster_size");
  hNumOriginHits_det9->GetYaxis()->SetTitle("Percentage");
  hNumOriginHits_det9->SetMinimum(0);
  hNumOriginHits_det9->SetMaximum(1.0);
  hNumOriginHits_det9->SetFillStyle(3001);
  hNumOriginHits_det9->Draw("HIST");

  c1->cd(22);
  cluster_u_Residual_det9->SetFillColor(kGreen);
  cluster_u_Residual_det9->GetXaxis()->SetTitle("good_Matched_Cluster_u [mm]");
  cluster_u_Residual_det9->GetYaxis()->SetTitle("Residual_u [mm]");
  cluster_u_Residual_det9->SetMarkerColor(9);
  cluster_u_Residual_det9->SetMarkerStyle(20);
  cluster_u_Residual_det9->SetMarkerSize(0.8);
  cluster_u_Residual_det9->Draw("P");

  c1->cd(23);
  cluster_v_Residual_det9->SetFillColor(kGreen);
  cluster_v_Residual_det9->GetXaxis()->SetTitle("good_Matched_Cluster_v [mm]");
  cluster_v_Residual_det9->GetYaxis()->SetTitle("Residual_v [mm]");
  cluster_v_Residual_det9->SetMarkerColor(4);
  cluster_v_Residual_det9->SetMarkerStyle(20);
  cluster_v_Residual_det9->SetMarkerSize(0.8);
  cluster_v_Residual_det9->Draw("P");

  c1->cd(24);
  hResidual_u_det9->SetLineColor(kBlue+2);
  hResidual_u_det9->SetLineWidth(2);
  hResidual_u_det9->Draw();

  TF1 *fitU9 = new TF1("fitU9", "gaus", -0.05, 0.05);
  fitU9->SetLineColor(kRed);
  hResidual_u_det9->Fit(fitU9, "R");
//  gStyle->SetOptFit(1111);

  TLatex latex9;
  latex9.SetNDC();
  latex9.SetTextSize(0.07);
  latex9.SetTextColor(kRed);
  latex9.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitU9->GetParameter(1)));
  latex9.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitU9->GetParameter(2)));


  c1->cd(25);
  hResidual_v_det9->SetLineColor(kGreen+2);
  hResidual_v_det9->SetLineWidth(2);
  hResidual_v_det9->Draw();

  TF1 *fitV9 = new TF1("fitV9", "gaus", -0.05, 0.05);
  fitV9->SetLineColor(kRed);
  hResidual_v_det9->Fit(fitV9, "R");

  TLatex latex10;
  latex10.SetNDC();
  latex10.SetTextSize(0.07);
  latex10.SetTextColor(kRed);
  latex10.DrawLatex(0.55, 0.80, Form("Mean = %.4f", fitV9->GetParameter(1)));
  latex10.DrawLatex(0.55, 0.70, Form("Sigma = %.4f", fitV9->GetParameter(2)));
 
  c1->Update();
 
  c1->SaveAs("origin_hits_by_detector.svg");
  
  TFile *outFile = new TFile("origin_hits_by_detector.root", "RECREATE");
  hNumOriginHits_det2->Write();
  hNumOriginHits_det3->Write();
  //hNumOriginHits_det32->Write();
  hNumOriginHits_det5->Write();
  hNumOriginHits_det7->Write();
  hNumOriginHits_det9->Write();
  
  cluster_u_Residual_det2->Write();
  cluster_v_Residual_det2->Write();
  cluster_u_Residual_det3->Write();
  cluster_v_Residual_det3->Write();
  cluster_u_Residual_det5->Write();
  cluster_v_Residual_det5->Write();
  cluster_u_Residual_det7->Write();
  cluster_v_Residual_det7->Write();
  cluster_u_Residual_det9->Write();
  cluster_v_Residual_det9->Write();
 
  hResidual_u_det2->Write();
  hResidual_v_det2->Write();
  hResidual_u_det3->Write();
  hResidual_v_det3->Write();
  hResidual_u_det5->Write();
  hResidual_v_det5->Write();
  hResidual_u_det7->Write();
  hResidual_v_det7->Write();
  hResidual_u_det9->Write();
  hResidual_v_det9->Write();

  c1->Write();
//  outFile->Close();
 
 // tfile->Close();
 // delete tfile;
 // std::fprintf(stdout,"OVER\n");
  std::fprintf(stdout, "Analysis complete. Results saved to origin_hits_by_detector.root\n");
  std::fprintf(stdout, "Press any key to exit...\n");
 // std::getc(stdin);
 
  return;
}
