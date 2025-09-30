void Efficiency_u_region(const std::string& rootFilePath){
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

    float region_u_min = -11.6;
    float region_u_max = 14;
    float region_v_min = -4.5;
    float region_v_max = 6.3;

    int nBinsU = 256;
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
    
    TEfficiency* pEffU = new TEfficiency("effU","Efficiency vs u;u [mm 256bin];Efficiency", nBinsU, region_u_min, region_u_max);

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
                bool hasMatch = (aMatchedMeasHit !=nullptr);
                pEffU->Fill (hasMatch, hit_u);
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
    //c1->Divide(2, 2);
    pEffU->Draw("AP");
    gPad->Update();
    pEffU->SetTitle("Efficiency vs u");
    pEffU->SetMarkerStyle(20);
    pEffU->SetMarkerSize(0.8);
    auto grU = pEffU->GetPaintedGraph();
    if (grU) {
      grU->GetYaxis()->SetRangeUser(0.0, 1);
    }
    TFile *outFile = new TFile("efficiency_projections.root", "RECREATE");
    c1->Write();
    pEffU->Write();
    outFile->Close();
    c1->SaveAs("efficiency_projections.svg");
    std::fprintf(stdout, "Efficiency projections saved as 'efficiency_projections.svg'\n");

    std::fprintf(stdout, "Analysis complete.\n");
    return;
}
  
