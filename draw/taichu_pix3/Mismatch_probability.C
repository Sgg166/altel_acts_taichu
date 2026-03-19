void Mismatch_probability(const std::string& rootFilePath){
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

    TEfficiency* pEffU = new TEfficiency("pEffU", "Mismatched Probability vs U;u [mm];Mismatch_Probability", 256, -11.6, 14);
    TEfficiency* pEffV = new TEfficiency("pEffV", "Mismatched Probability vs V;v [mm];Mismatch_Probability", 128, -5.5, 7.3);
    TEfficiency* pEff2D = new TEfficiency("eff2D","Mismatched Probability vs Position;u [mm 32bin];v [mm 16bin];Mismatch_Probability", 32 , -11.6, 14, 16, -5.5, 7.3);

    TH2F* hTotalTracks = new TH2F("hTotalTracks", "Total Tracks; fit_u; fit_v", 32, -11.6, 14, 16, -5.5, 7.3);
    TH2F* huv_fit_det32 = new TH2F("huv_fit_det32", "huv_fit_det32; fit_u; fit_v", 32, -11.6, 14, 16, -5.5, 7.3);
    TH2F* huv_no_matched_fit_det32 = new TH2F("huv_no_matched_fit_det32", "huv_no_matched_fit_det32; no_matched_fit_u; no_matched_fit_v", 32, -11.6, 14, 16, -5.5, 7.3);
   
    size_t totalNumEvents = ttreeReader.numEvents();
    size_t outEventNum = 0;
    size_t noclusterNum = 0;

    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
      std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
      int goodTrajectoryCount = 0;
      for(auto aTraj: telEvent->trajs()){
        if(aTraj->numOriginMeasHit() == 5){
          goodTrajectoryCount++;
        }
      }
      if(goodTrajectoryCount != 1){
        continue;
      }
      
      for(auto aTraj: telEvent->trajs()){
        if(aTraj->numOriginMeasHit() != 5) continue;
      
        for(auto &aTrajHit: aTraj->trajHits()){
          auto aFitHit = aTrajHit->fitHit();
          auto aMatchedMeasHit = aTrajHit->matchedMeasHit();
          if(!aFitHit) continue;
          if(aFitHit->detN() != 32) continue;
          double fit_u = aFitHit->u();
          double fit_v = aFitHit->v();
          if(fit_u < -12.8 || fit_u > 12.8 || fit_v < -6.4 || fit_v > 6.4){
            continue;
          }
          outEventNum++;
          hTotalTracks->Fill(fit_u, fit_v);
          bool unmatched = (aMatchedMeasHit == nullptr);
          pEffU->Fill(unmatched, fit_u);
          pEffV->Fill(unmatched, fit_v);
          pEff2D->Fill(unmatched, fit_u, fit_v );

          if(!unmatched){
            huv_fit_det32->Fill(fit_u, fit_v);
          }else{
            noclusterNum++;
            huv_no_matched_fit_det32->Fill(fit_u, fit_v);
          }
        }
      }
    }

    double total_MismatchProbability = (outEventNum > 0)? double(noclusterNum)/double(outEventNum) : 0.0;

    std::fprintf(stdout,"found %.3f Mismatch Probability\n", total_MismatchProbability);
    std::fprintf(stdout,"found %zu outEventNum\n", outEventNum);
    std::fprintf(stdout,"found %zu noclusterNum\n", noclusterNum);

    TCanvas* c1 = new TCanvas("c1", "Mismatch Probability Projection in U or V", 1600, 1600);
    c1->Divide(1,2);

    c1->cd(1);
    pEffU->SetMarkerStyle(20);
    pEffU->Draw("AP");
    gPad->Update();
    pEffU->SetMarkerSize(0.6);
    auto grU = pEffU->GetPaintedGraph();
    if (grU) {
      grU->GetYaxis()->SetRangeUser(0.0, 1);
    }

    c1->cd(2);
    pEffV->SetMarkerStyle(20);
    pEffV->Draw("AP");
    gPad->Update();
    pEffV->SetMarkerSize(0.6);
    auto grV = pEffV->GetPaintedGraph();
    if (grV) {
      grV->GetYaxis()->SetRangeUser(0.0, 1);
    }

    TCanvas* c2 = new TCanvas("c2", "Mismatch Probability Projection", 1600, 1600);
    c2->SetFillStyle(1001);

//    gStyle->SetPaintTextFormat(".2f");
    pEff2D->Draw("TEXT colz");
    gPad->Update();
    //pEff2D->SetStats(1);
    pEff2D->SetMarkerStyle(20);
    pEff2D->SetMarkerSize(1);
    pEff2D->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.0, 1);
    TLatex  latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextColor(kRed);
    latex.DrawLatex(0.65, 0.90, Form("total MismatchProbability: %.3f", total_MismatchProbability));

    TCanvas *c3 = new TCanvas("c3", "Residuals", 1000, 1000);
    c3->Divide(2,2);

    c3->cd(1);
    hTotalTracks->Draw("colz text");

    c3->cd(2);
    huv_fit_det32->Draw("colz text");

    c3->cd(3);
    huv_no_matched_fit_det32->Draw("colz text");

    // Save
    TFile *outFile = new TFile("/home/pub/B/altel_acts/Draw/root_file/Mismatch_Probability.root", "RECREATE");
    pEff2D->Write();
    pEffU->Write();
    pEffV->Write();
    hTotalTracks->Write();
    huv_no_matched_fit_det32->Write();
    huv_fit_det32->Write();
    c1->Write();
    c2->Write();

    std::fprintf(stdout, "Analysis complete. Result saved.\n");
}

