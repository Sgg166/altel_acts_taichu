void ep(const std::string& rootFilePath){
    std::cout << rootFilePath << std::endl;
    TFile* tfile = new TFile(rootFilePath.c_str(), "READ");
    if(!tfile || !tfile->IsOpen()){
        std::fprintf(stderr, "tfile is not open\n");
        throw;
    }
    TTree* pTree = 0;
    tfile->GetObject("eventTree", pTree);
    if(!pTree){
        std::fprintf(stderr, "pTree is invalid\n");
        throw;
    }
    altel::TelEventTTreeReader ttreeReader;
    ttreeReader.setTTree(pTree);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    size_t totalTracks = 0;
    size_t goodTracks = 0;
    size_t hit_dut_goodTracks = 0;
    size_t totalMatchedHits = 0;
    size_t totalFitHits = 0;
    size_t totalMatchedTracks = 0;
    size_t totalMatchedHits_inPixel = 0;
    size_t totalFitHits_inPixel = 0;
    size_t goodTracks_inPixel = 0;
    size_t Invalid_tracks =0;
 
    //Create pixel-level efficiency histograms
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",5, -12.5, 12.5, 5, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 5, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 5, -12.5, 12.5);
 
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",5, -12.5, 12.5, 5, -12.5, 12.5);
 
    TH2I* hNumMatchedDist = new TH2I("hNumMatchedDist", "NumMatched Distribution per pixel bin;u [um];v [um];NumMatched", 5, -12.5, 12.5, 5, -12.5, 12.5);
    TH1I* hNumMatched1D = new TH1I("hNumMatched1D","NumMatched Distribution;NumMatched;Count", 10, 0, 10);
 
    // Create 25 histograms for cluster size distribution in each bin
    std::vector<TH1I*> hClusterSizeBins;
    for (int ix = 0; ix < 5; ix++) {
        for (int iy = 0; iy < 5; iy++) {
            TString name = Form("hClusterSize_bin_%d_%d", ix, iy);
            TString title = Form("Cluster Size Distribution (bin u=%d, v=%d);Cluster Size;Count", ix, iy);
            hClusterSizeBins.push_back(new TH1I(name, title, 30, 0, 15));
        }
    }
 
    // Coordinates map function: Maps any coordinate to the range [-12.5um, +12.5um]
    auto mapToPixel = [](double coord_mm) {
      double coord_um = coord_mm * 1000.0;
      return coord_um - 25.0 * std::round(coord_um / 25.0);
    };
 
    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
        for(auto aTraj: telEvent->trajs()){
            totalTracks++;
 
            bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);
 
            if(!isGoodTrack){
              continue;
            }
 
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
              totalFitHits++;
 
              double hit_u = aFitHit->u();
              double hit_v = aFitHit->v();
              if(hit_u < -12.8 || hit_u > 12.8 || hit_v < -6.4 || hit_v > 6.4 ){
                Invalid_tracks++;
                continue;
              }
 
              // Map to the pixel cycle
              double pixel_u = mapToPixel(hit_u);
              double pixel_v = mapToPixel(hit_v);
 
              // Check if it is within the pixel range
              bool isInPixel = (pixel_u >= -12.5 && pixel_u <= 12.5 && pixel_v >= -12.5 && pixel_v <= 12.5);
              if(isInPixel){
                hit_dut_goodTracks++;
                hTrackCount_pixel->Fill(pixel_u, pixel_v);
                bool hasMatch = (aMatchedMeasHit != nullptr);
                pEff2D_pixel->Fill(hasMatch, pixel_u, pixel_v);
                pEffU_pixel->Fill(hasMatch, pixel_u);
                pEffV_pixel->Fill(hasMatch, pixel_v);
                if(aMatchedMeasHit){
                  totalMatchedHits_inPixel++;
                  size_t numMatched = aMatchedMeasHit->measRaws().size();
                  hNumMatchedDist->Fill(pixel_u, pixel_v, numMatched);
                  hNumMatched1D->Fill(numMatched);
                  
                  // Fill cluster size for the corresponding bin
                  // 使用TH2I的FindBin函数确保bin计算一致性
                  int root_binx = hNumMatchedDist->GetXaxis()->FindBin(pixel_u);
                  int root_biny = hNumMatchedDist->GetYaxis()->FindBin(pixel_v);
                  
                  // ROOT bin从1开始，转换为0-based索引
                  int binx = root_binx - 1;
                  int biny = root_biny - 1;
                  
                  if (binx >= 0 && binx < 5 && biny >= 0 && biny < 5) {
                    int idx = biny * 5 + binx; // Linear index
                    hClusterSizeBins[idx]->Fill(numMatched);
                  }
                }
              }
            }
        }
    }
 
    // Verification: Compare hClusterSizeBins and hNumMatchedDist
    std::fprintf(stdout, "\n=== Verification: Cluster Size Statistics ===\n");
    bool all_match = true;
    for(int i=0; i<25; i++){
        // 直接计算hClusterSizeBins中所有值的总和，避免浮点数精度问题
        double sum_from_cluster = 0;
        for(int bin=1; bin<=hClusterSizeBins[i]->GetNbinsX(); bin++){
            double bin_center = hClusterSizeBins[i]->GetBinCenter(bin);
            double bin_content = hClusterSizeBins[i]->GetBinContent(bin);
            sum_from_cluster += bin_center * bin_content;
        }
        
        // Get corresponding hNumMatchedDist bin value
        int binx = i % 5;
        int biny = i / 5;
        double sum_from_dist = hNumMatchedDist->GetBinContent(binx+1, biny+1);
        
        bool match = (fabs(sum_from_cluster - sum_from_dist) < 0.001);
        all_match = all_match && match;
        
        std::fprintf(stdout, "Bin %2d (u=%d,v=%d): ClusterSum=%8.1f, DistSum=%8.1f, Match=%s\n", 
               i, binx, biny, sum_from_cluster, sum_from_dist, match ? "YES" : "NO");
    }
    std::fprintf(stdout, "All bins match: %s\n", all_match ? "YES" : "NO");
 
    double pixel_efficiency = (hit_dut_goodTracks > 0) ?
        (double)totalMatchedHits_inPixel / (double)hit_dut_goodTracks : 0.0;
 
    std::fprintf(stdout, "\n=== Pixel level efficiency result ===\n");
    std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
    std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
    std::fprintf(stdout, "Good tracks (==5 hits ): %zu\n", goodTracks);
    std::fprintf(stdout, "Good tracks (==5 hits) and hit dut: %zu\n", hit_dut_goodTracks);
    std::fprintf(stdout, "Good tracks but no hit dut: %zu\n", Invalid_tracks);
    std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
    std::fprintf(stdout, "Total MatchedHits in pixel: %zu\n", totalMatchedHits_inPixel);
    std::fprintf(stdout, "Pixel Level Efficiency: %.3f\n", pixel_efficiency);
 
    TCanvas* c1 = new TCanvas("c1", "Pixel Level Efficiency", 1000, 1000);
    c1->Divide(2, 2);
    c1->cd(1);
    c1->SetFillStyle(1001);
    c1->SetFillColor(kWhite);
 
    pEff2D_pixel->Draw("TEXT colz");
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetFixedAspectRatio();
    gPad->Update();
    pEff2D_pixel->SetMarkerStyle(20);
    pEff2D_pixel->SetMarkerSize(2);
    pEff2D_pixel->SetTitle("Pixel Efficiency vs Position");
    pEff2D_pixel->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.88, 1);
 
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextColor(kRed);
    latex.DrawLatex(0.5, 0.90, Form("Pixel Efficiency: %.3f", pixel_efficiency));
 
    c1->cd(2);
    hTrackCount_pixel->Draw("TEXT colz");
    hTrackCount_pixel->SetStats(0);
    hTrackCount_pixel->SetMarkerSize(2);
    hTrackCount_pixel->SetTitle("Track count per pixel bin");
 
    c1->cd(3);
    pEffU_pixel->Draw("AP");
    gPad->Update();
    pEffU_pixel->SetTitle("Pixel Efficiency vs u");
    pEffU_pixel->SetMarkerStyle(20);
    pEffU_pixel->SetMarkerSize(0.8);
    auto grU_pixel = pEffU_pixel->GetPaintedGraph();
    if (grU_pixel) {
        grU_pixel->GetYaxis()->SetRangeUser(0.88, 1);
    }
 
    c1->cd(4);
    pEffV_pixel->Draw("AP");
    gPad->Update();
    pEffV_pixel->SetTitle("Pixel Efficiency vs v");
    pEffV_pixel->SetMarkerStyle(20);
    pEffV_pixel->SetMarkerSize(0.8);
    auto grV_pixel = pEffV_pixel->GetPaintedGraph();
    if (grV_pixel) {
        grV_pixel->GetYaxis()->SetRangeUser(0.88, 1);
    }
 
    c1->SaveAs("pixel_efficiency.svg");
    std::fprintf(stdout, "Pixel efficiency plots saved as 'pixel_efficiency.svg'\n");
 
    TCanvas* c3 = new TCanvas("c3", "NumMatched Distribution", 1200, 500);
    c3->Divide(2, 1);
 
    c3->cd(1);
    hNumMatchedDist->Draw("TEXT colz");
    hNumMatchedDist->SetStats(0);
    hNumMatchedDist->SetMarkerSize(2);
    hNumMatchedDist->SetTitle("NumMatched Distribution per pixel bin");
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetFixedAspectRatio();
 
    c3->cd(2);
    hNumMatched1D->Draw();
    hNumMatched1D->SetTitle("NumMatched Distribution");
    hNumMatched1D->SetStats(1);
 
    c3->SaveAs("nummatched_distribution.svg");
 
    TCanvas* c5 = new TCanvas("c5", "Cluster Size Distribution per Pixel Bin", 2000, 2000);
    c5->Divide(5,5);
 
    for (int i = 0; i < 25; i++) {
      c5->cd(i+1);
      if (hClusterSizeBins[i]->GetEntries() > 0) {
        hClusterSizeBins[i]->Scale(1.0 / hClusterSizeBins[i]->GetEntries()); // Normalize to fraction
      }
 
      hClusterSizeBins[i]->Draw();
      gPad->Update();
      hClusterSizeBins[i]->SetStats(1); // Show statistics
      hClusterSizeBins[i]->SetFillColor(kBlue);
      hClusterSizeBins[i]->SetFillStyle(3001);
    }
 
    c5->SaveAs("cluster_size_per_bin.svg");
 
    std::fprintf(stdout, "Pixel level analysis complete.\n");
    return;
}
