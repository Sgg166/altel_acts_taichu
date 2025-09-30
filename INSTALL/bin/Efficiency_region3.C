void Efficiency_region3(const std::string& rootFilePath){
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

    float region_u_min = -6.4;
    float region_u_max = 6.4;
    float region_v_min = -3.2;
    float region_v_max = 3.2;

    int nBinsU = 27;
    int nBinsV = 27;

    size_t totalNumEvents = ttreeReader.numEvents();
    size_t totalTracks = 0;
    size_t goodTracks = 0;
    size_t hit_dut_goodTracks = 0;
    size_t totalMatchedHits = 0;
    size_t totalFitHits = 0;
    size_t totalMatchedTracks=0;
    size_t totalMatchedHits_inRegion = 0;
    size_t totalFitHits_inRegion = 0;
    size_t goodTracks_inRegion = 0;
    size_t Invalid_tracks =0;
    TEfficiency* pEff2D = new TEfficiency("eff2D","Efficiency vs Position;u [mm 10bin];v [mm 27bin];Efficiency",nBinsU, region_u_min, region_u_max, nBinsV, region_v_min, region_v_max);
    TEfficiency* pEffU = new TEfficiency("effU","Efficiency vs u;u [mm 10bin];Efficiency", nBinsU, region_u_min, region_u_max);
    TEfficiency* pEffV = new TEfficiency("effV","Efficiency vs v;v [mm 27bin];Efficiency", nBinsV, region_v_min, region_v_max);
 
    TH2I* hTrackCount = new TH2I("hTrackCount", "Track count per bin;u [mm 10bin];v [mm 27bin];Count",nBinsU, region_u_min, region_u_max,nBinsV, region_v_min, region_v_max);

    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        for(auto aTraj: telEvent->trajs()){
            totalTracks++;

            bool isGoodTrack = (aTraj->numOriginMeasHit() ==5 );
            //goodTracks++;            

            if(!isGoodTrack){
              continue;
            }
            bool hasMatchedInDUT = false;
            bool hasFitInRegion = false;
            for(auto &aTrajHit: aTraj->trajHits()){
              auto aFitHit = aTrajHit->fitHit();
              auto aMatchedMeasHit = aTrajHit->matchedMeasHit();

              if(!aFitHit){
                continue;
              }
              if (aFitHit->detN() != 32){
                continue;
              }
              goodTracks++;
              double hit_u = aFitHit->u();
              double hit_v = aFitHit->v();
              if(hit_u<=-12.8 || hit_u>=12.8 || hit_v<=-6.4 || hit_v>=6.4 ){
                Invalid_tracks++;
                continue;
               // std::fprintf(stdout,"(hit_u,hit_v) :(%f,%f)\n", hit_u,hit_v);
              }
              hit_dut_goodTracks++;
              if (aMatchedMeasHit) {
                totalMatchedHits++;
                hasMatchedInDUT = true;//测试totalMatchedHits是否=totalMatchedTracks
              }
              bool isInRegion = (hit_u >= region_u_min && hit_u <= region_u_max &&
                                       hit_v >= region_v_min && hit_v <= region_v_max);
              if(isInRegion){
                goodTracks_inRegion++;
                totalFitHits_inRegion++;
                hasFitInRegion = true;       
                hTrackCount->Fill(hit_u, hit_v);
                bool hasMatch = (aMatchedMeasHit !=nullptr);
                pEff2D->Fill(hasMatch, hit_u, hit_v );
                pEffU->Fill (hasMatch, hit_u);
                pEffV->Fill (hasMatch, hit_v);                                                               
                if(aMatchedMeasHit){
                  totalMatchedHits_inRegion++;
                }
              }
              if (hasMatchedInDUT) {
                totalMatchedTracks++;
              }
            }
        }
    }

    double efficiency = (hit_dut_goodTracks > 0) ?
                       (double)totalMatchedTracks / (double)hit_dut_goodTracks : 0.0;
    double Error=(hit_dut_goodTracks > 0) ?
                      std::sqrt((( (double)hit_dut_goodTracks*(double)totalMatchedTracks)-((double)totalMatchedTracks*(double)totalMatchedTracks))/((double)hit_dut_goodTracks*(double)hit_dut_goodTracks*(double)hit_dut_goodTracks)) : 0.0;
    double local_efficiency = (goodTracks_inRegion > 0) ?
                             (double)totalMatchedHits_inRegion / (double)goodTracks_inRegion : 0.0;
    double local_Error=(goodTracks_inRegion > 0) ?
                        std::sqrt((( (double)goodTracks_inRegion*(double)totalMatchedHits_inRegion)-((double)totalMatchedHits_inRegion*(double)totalMatchedHits_inRegion))/((double)goodTracks_inRegion*(double)goodTracks_inRegion*(double)goodTracks_inRegion)) : 0.0;

    std::fprintf(stdout, "=== Overall efficiency result ===\n");
    std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
    std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
    std::fprintf(stdout, "Good tracks (==5 hits ): %zu\n", goodTracks);
    std::fprintf(stdout, "Good tracks (==5 hits) and hit dut: %zu\n", hit_dut_goodTracks);
    std::fprintf(stdout, "Good tracks but no hit dut: %zu\n", Invalid_tracks);
    std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
    std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
    std::fprintf(stdout, "Total Matched Tracks: %zu\n", totalMatchedTracks);
    std::fprintf(stdout, "Overall Efficiency: %.3f\n", efficiency);
    std::fprintf(stdout, "Overall Efficiency Error: %.5f\n", Error);

    std::fprintf(stdout, "\n=== Local efficiency result ===\n");
    std::fprintf(stdout, "Region bounds: u[%.1f, %.1f], v[%.1f, %.1f]\n",
                 region_u_min, region_u_max, region_v_min, region_v_max);
    std::fprintf(stdout, "Good tracks in region: %zu\n", goodTracks_inRegion);
    std::fprintf(stdout, "Total fit hits in region: %zu\n", totalFitHits_inRegion);
    std::fprintf(stdout, "Total MatchedHits in region: %zu\n", totalMatchedHits_inRegion);
    std::fprintf(stdout, "Local Efficiency: %.3f\n", local_efficiency);
    std::fprintf(stdout, "Local Efficiency Error: %.5f\n", local_Error);

    TCanvas* c1 = new TCanvas("c1", "Efficiency vs Position", 1000, 1000);
    c1->Divide(2, 2);
    c1->cd(1);
    c1->SetFillStyle(1001);
    c1->SetFillColor(kWhite);
   
    pEff2D->Draw("TEXT colz");
    gStyle->SetPaintTextFormat(".2f");
    //gStyle->SetPaintTextSize(0.002);
    gPad->Update();
    pEff2D->SetMarkerStyle(20);
    pEff2D->SetMarkerSize(1);
    pEff2D->SetTitle("Efficiency vs Position");
   // pEff2D->GetXaxis()->SetTitle("u [mm]");
   // pEff2D->GetPaintedHistogram()->GetYaxis()->SetTitle("v [mm]");
   // pEff2D->GetPaintedHistogram()->GetZaxis()->SetTitle("Efficiency");
    pEff2D->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.65, 1);

   // c1->cd(2);
    TLatex  latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextColor(kRed);
   // latex.DrawLatex(0.1, 0.70, Form("Total Events: %zu", totalNumEvents));
   // latex.DrawLatex(0.1, 0.60, Form("Good Tracks: %zu", goodTracks));
   // latex.DrawLatex(0.65, 0.90, Form("Overall Efficiency: %.3f", efficiency));
    latex.DrawLatex(0.65, 0.90, Form("Local Efficiency: %.3f", local_efficiency));
   // latex.DrawLatex(0.1, 0.30, Form("Region: u[%.1f, %.1f], v[%.1f, %.1f]",
     //                                 region_u_min, region_u_max, region_v_min, region_v_max));

    c1->cd(2);
    hTrackCount->Draw("TEXT colz");
    hTrackCount->SetStats(0);
    hTrackCount->SetMarkerSize(1);
    hTrackCount->SetTitle("Track count per bin");
    //hTrackCount->GetZaxis()->SetRangeUser(30, 600);
    
    c1->cd(3);
    pEffU->Draw("AP");
    gPad->Update();
    pEffU->SetTitle("Efficiency vs u");
    pEffU->SetMarkerStyle(20);
    pEffU->SetMarkerSize(0.8);
    auto grU = pEffU->GetPaintedGraph();
    if (grU) {
      grU->GetYaxis()->SetRangeUser(0.0, 1);
    }

    c1->cd(4);
//    pEffV->SetStats(1);

    pEffV->Draw("  AP");
    gPad->Update();
    pEffV->SetTitle("Efficiency vs v");
    pEffV->SetMarkerStyle(20);
    pEffV->SetMarkerSize(0.8);
    auto grV = pEffV->GetPaintedGraph();
    if (grV) {
      grV->GetYaxis()->SetRangeUser(0.0, 1);
    }
   
    TFile *outFile = new TFile("efficiency_projections.root", "RECREATE");
    c1->Write();
    pEff2D->Write();
    hTrackCount->Write();
    pEffU->Write();
    pEffV->Write();
    outFile->Close();
    c1->SaveAs("efficiency_projections.svg");
    std::fprintf(stdout, "Efficiency projections saved as 'efficiency_projections.svg'\n");

    std::fprintf(stdout, "Analysis complete.\n");
    return;
}
  
