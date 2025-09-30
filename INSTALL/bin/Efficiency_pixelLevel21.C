void Efficiency_pixelLevel21(const std::string& rootFilePath){
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
    size_t cluster_size=0;
 
    size_t centerPixelResponses = 0;
    size_t neighborPixelResponses = 0;
    size_t totalResponses = 0;
 
    
    std::map<uint32_t, size_t> pixelHitCount;  //The number of hits per pixel
    std::vector<uint32_t> allPixelIds;         //Record the numbers of all response pixels
    size_t uniquePixels = 0;                    // Hit the number of only one pixel
 
    //Create pixel-level efficiency histograms
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",5, -12.5, 12.5, 5, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 5, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 5, -12.5,12.5);
 
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",5, -12.5, 12.5, 5, -12.5,12.5);
 
    TH1I* hPixelIdDist = new TH1I("hPixelIdDist", "Pixel ID Distribution;Pixel ID;Count",100, 0, 524288); //Cover every pixel
    TH2I* hPixelUVMap = new TH2I("hPixelUVMap", "Pixel Hit Map (UV);u;v;Count",1024, 0, 1024, 512, 0, 512);
 
    //Response position distribution histogram
    TH2I* hResponseOffset = new TH2I("hResponseOffset", "Response Offset Distribution;#Delta u [um];#Delta v [um];Count",50, -25, 25, 50, -25, 25);
    TH1I* hResponseDistance = new TH1I("hResponseDistance", "Response Distance from Center;Distance [um];Count",50, 0, 35);
    TH1I* hResponseType = new TH1I("hResponseType", "Response Type;Type;Count", 3, 0, 3);
    // Set the response type label
    hResponseType->GetXaxis()->SetBinLabel(1, "Center");
    hResponseType->GetXaxis()->SetBinLabel(2, "Neighbor");
    hResponseType->GetXaxis()->SetBinLabel(3, "Far");
 
    // Coordinates map function: Maps any coordinate to the range [-12.5um, +12.5um]
    auto mapToPixel = [](double coord_mm) {
      double coord_um = coord_mm * 1000.0;
      double mapped = coord_um-25*std::floor((coord_um+12.5)/25);
      if (mapped < 0) mapped += 25.0;
      return mapped - 12.5;
    };
 
    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
        for(auto aTraj: telEvent->trajs()){
            totalTracks++;
 
            bool isGoodTrack = (aTraj->numOriginMeasHit() ==5 );
 
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
                  cluster_size += numMatched;
 
                  // Obtain the numbers of all response pixels
                  for(auto& raw : aMatchedMeasHit->measRaws()){
                    uint32_t pixelId = raw.pixelId();  // useing pixelId()
                    allPixelIds.push_back(pixelId);
                    pixelHitCount[pixelId]++;  // Count the number of hits of this pixel
                    
                    hPixelIdDist->Fill(pixelId);
                    hPixelUVMap->Fill(raw.u(), raw.v());
                  }
 
                  // Analyze the distribution of response positions
                  double fit_u = aFitHit->u() * 1000.0; // 转换为微米
                  double fit_v = aFitHit->v() * 1000.0;
 
                  double meas_u = aMatchedMeasHit->u() * 1000.0;
                  double meas_v = aMatchedMeasHit->v() * 1000.0;
 
                  double delta_u = meas_u - fit_u;
                  double delta_v = meas_v - fit_v;
                  double distance = std::sqrt(delta_u * delta_u + delta_v * delta_v);
                  
                  hResponseOffset->Fill(delta_u, delta_v);
                  hResponseDistance->Fill(distance);
 
                  //Determine the response type
                  totalResponses++;
                  if (distance <= 12.5) {
                      centerPixelResponses++;
                      hResponseType->Fill(0); // Center
                  } else if (distance <= 25.0) {
                      neighborPixelResponses++;
                      hResponseType->Fill(1); // Neighbor
                  } else {
                      hResponseType->Fill(2); // Far
                  }
                }
              }
            }
        }
    }
 
    double pixel_efficiency = (hit_dut_goodTracks > 0) ?
        (double)totalMatchedHits_inPixel / (double)hit_dut_goodTracks : 0.0;
 
    std::fprintf(stdout, "=== Pixel level efficiency result ===\n");
    std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
    std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
    std::fprintf(stdout, "Good tracks (==5 hits ): %zu\n", goodTracks);
    std::fprintf(stdout, "Good tracks (==5 hits) and hit dut: %zu\n", hit_dut_goodTracks);
    std::fprintf(stdout, "Good tracks but no hit dut: %zu\n", Invalid_tracks);
    std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
    std::fprintf(stdout, "Total MatchedHits in pixel: %zu\n", totalMatchedHits_inPixel);
    std::fprintf(stdout, "cluster_size: %zu\n", cluster_size);
    std::fprintf(stdout, "Pixel Level Efficiency: %.3f\n", pixel_efficiency);
 
   
    std::fprintf(stdout, "\n=== Response Position Analysis ===\n");
    std::fprintf(stdout, "Total responses analyzed: %zu\n", totalResponses);
    std::fprintf(stdout, "Center pixel responses: %zu (%.2f%%)\n",
                centerPixelResponses, totalResponses > 0 ? (double)centerPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Neighbor pixel responses: %zu (%.2f%%)\n",
                neighborPixelResponses, totalResponses > 0 ? (double)neighborPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Far responses: %zu (%.2f%%)\n",
                totalResponses - centerPixelResponses - neighborPixelResponses,
                totalResponses > 0 ? (double)(totalResponses - centerPixelResponses - neighborPixelResponses)/totalResponses*100 : 0);
 
   
    std::fprintf(stdout, "\n=== Pixel ID Analysis ===\n");
    std::fprintf(stdout, "Total pixel responses: %zu\n", allPixelIds.size());
    std::fprintf(stdout, "Unique pixels fired: %zu\n", pixelHitCount.size());
 
    //Find the most active pixel
    if(!pixelHitCount.empty()){
      auto maxPixel = std::max_element(pixelHitCount.begin(), pixelHitCount.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; });
      
      auto [u_max, v_max] = altel::TelMeasRaw::pixelIdToUV(maxPixel->first);
      std::fprintf(stdout, "Most active pixel: ID=%u, (u=%d, v=%d), hits=%zu\n",
                  maxPixel->first, u_max, v_max, maxPixel->second);
    }
 
    // Calculate the average number of responses per pixel
    if(pixelHitCount.size() > 0){
      double avgHitsPerPixel = (double)allPixelIds.size() / pixelHitCount.size();
      std::fprintf(stdout, "Average hits per pixel: %.2f\n", avgHitsPerPixel);
    }
 
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
 
    TCanvas* c4 = new TCanvas("c4", "Response Position Analysis", 1200, 800);
    c4->Divide(2, 2);
 
    c4->cd(1);
    hResponseOffset->Draw("TEXT colz");
    hResponseOffset->SetStats(0);
    hResponseOffset->SetTitle("Response Offset Distribution");
    hResponseOffset->GetXaxis()->SetTitle("#Delta u [um]");
    hResponseOffset->GetYaxis()->SetTitle("#Delta v [um]");
    gPad->SetFixedAspectRatio();
 
    c4->cd(2);
    hResponseDistance->Draw();
    hResponseDistance->SetTitle("Response Distance from Center");
    hResponseDistance->GetXaxis()->SetTitle("Distance [um]");
    hResponseDistance->GetYaxis()->SetTitle("Count");
    hResponseDistance->SetStats(1);
 
    c4->cd(3);
    hResponseType->Draw();
    hResponseType->SetTitle("Response Type Distribution");
    hResponseType->SetStats(0);
    hResponseType->SetFillStyle(3001);
    hResponseType->SetFillColor(kBlue);
 
    c4->SaveAs("response_position_analysis.svg");
    std::fprintf(stdout, "Response position analysis plots saved as 'response_position_analysis.svg'\n");
 
   
    TCanvas* c5 = new TCanvas("c5", "Pixel ID Analysis", 1200, 600);
    c5->Divide(2, 1);
 
    c5->cd(1);
    hPixelIdDist->Draw();
    hPixelIdDist->SetMarkerStyle(20);
    hPixelIdDist->SetMarkerColor(2);
    hPixelIdDist->SetMarkerSize(1.0);
    hPixelIdDist->SetTitle("Pixel ID Distribution");
    hPixelIdDist->GetXaxis()->SetTitle("Pixel ID");
    hPixelIdDist->GetYaxis()->SetTitle("Count");
 
    c5->cd(2);
    hPixelUVMap->Draw("colz");
    hPixelUVMap->SetStats(0);
    hPixelUVMap->SetTitle("Pixel Hit Map (UV Coordinates)");
    hPixelUVMap->GetXaxis()->SetTitle("u coordinate");
    hPixelUVMap->GetYaxis()->SetTitle("v coordinate");
    //hPixelUVMap->GetZaxis()->SetRangeUser(2, 6);
 
    c5->SaveAs("pixel_id_analysis.svg");
    std::fprintf(stdout, "Pixel ID analysis plots saved as 'pixel_id_analysis.svg'\n");
 
    std::fprintf(stdout, "Pixel level analysis complete.\n");
    return;
}
