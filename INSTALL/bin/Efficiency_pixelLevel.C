void Efficiency_pixelLevel(const std::string& rootFilePath){
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
    //Create pixel-level efficiency histograms
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",10, -12.5, 12.5, 10, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 10, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 10, -12.5,12.5);
 
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",10, -12.5, 12.5, 10, -12.5,12.5);
    
    TH2I* hNumMatchedDist = new TH2I("hNumMatchedDist", "NumMatched Distribution per pixel bin;u [um];v [um];NumMatched", 10, -12.5, 12.5, 10, -12.5, 12.5);
    TH1I* hNumMatched1D = new TH1I("hNumMatched1D","NumMatched Distribution;NumMatched;Count", 30, 0, 15); 
  //  TH1F *hClusterSize  = new TH1F("local_pixel_Cluster_size_det", "local_pixel_Cluster_size_det ", 30, 0,15);
    
    std::vector<TH1I*> hClusterSizeBins;
    for (int ix = 0; ix < 5; ix++) {
        for (int iy = 0; iy < 5; iy++) {
            TString name = Form("hClusterSize_bin_%d_%d", ix, iy);
            TString title = Form("Cluster Size Distribution (bin u=%d, v=%d);Cluster Size;Count", ix, iy);
            hClusterSizeBins.push_back(new TH1I(name, title, 30, 0, 15));
        }
    }
  
 
    // Coordinates map function: Maps any coordinate to the range [-12.5um, +12.5um]
 /*   auto mapToPixel = [](double coord_mm) {
      double coord_um = coord_mm * 1000.0;
 //      return coord_um;
      return coord_um - 25.0 * std::round(coord_um / 25.0);
    };*/
    auto mapToPixel = [](double coord_mm) {
      double coord_um = coord_mm * 1000.0;
      double mapped = coord_um-25*std::floor((coord_um+12.5)/25);//Replace the coordinate range with[0,25]
      if (mapped < 0) mapped += 25.0;//
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

           // bool hasMatchedInDUT = false;
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
               // std::fprintf(stdout,"(hit_u,hit_v) :(%f,%f)\n", hit_u,hit_v);
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
                  hNumMatchedDist->Fill(pixel_u, pixel_v, numMatched);
                  hNumMatched1D->Fill(numMatched);
                /*  bool is_local_pixel =(pixel_u >= -7.5 && pixel_u <= -2.5 && pixel_v >= -7.5 && pixel_v <= -2.5);
                  if(is_local_pixel)
                  {
                    hClusterSize->Fill(numMatched);
                  }*/
               
                  int binx = int((pixel_u + 12.5) / 5.0); // u方向 [0..4]
                  int biny = int((pixel_v + 12.5) / 5.0); // v方向 [0..4]
                  if (binx >= 0 && binx < 5 && biny >= 0 && biny < 5) {
                    int idx = biny * 5 + binx; 
                    hClusterSizeBins[idx]->Fill(numMatched);
                  }
                }
              }
            }
        }
    }
 
  
    double pixel_efficiency = (hit_dut_goodTracks > 0) ?  
                              (double)totalMatchedHits_inPixel / (double)hit_dut_goodTracks : 0.0;
    double Error=(hit_dut_goodTracks > 0) ?
                  std::sqrt((( (double)hit_dut_goodTracks*(double)totalMatchedHits_inPixel)-((double)totalMatchedHits_inPixel*(double)totalMatchedHits_inPixel))/((double)hit_dut_goodTracks*(double)hit_dut_goodTracks*(double)hit_dut_goodTracks)) : 0.0;

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
    std::fprintf(stdout, "Overall Efficiency Error: %.5f\n", Error);


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
    pEff2D_pixel->SetMarkerSize(1.6);
    pEff2D_pixel->SetTitle("Pixel Efficiency vs Position");
    pEff2D_pixel->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.86, 1);
 
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

 
    TCanvas* c2 = new TCanvas("c2", "NumMatched Distribution", 1200, 500);
    c2->Divide(2, 1);

    c2->cd(1);
    hNumMatchedDist->Draw("TEXT colz");
    hNumMatchedDist->SetStats(0);
    hNumMatchedDist->SetMarkerSize(2);
    hNumMatchedDist->SetTitle("NumMatched Distribution per pixel bin");
    gStyle->SetPaintTextFormat(".2f");
    gPad->SetFixedAspectRatio();

    c2->cd(2);
    hNumMatched1D->Draw();
    hNumMatched1D->SetTitle("NumMatched Distribution");
    hNumMatched1D->SetStats(1);

    c2->SaveAs("nummatched_distribution.svg");

    TCanvas* c3 = new TCanvas("c3", "Cluster Size Distribution per Pixel Bin", 2500, 2000);
    c3->Divide(5,5);

    for (int i = 0; i < 25; i++) {
      c3->cd(i+1);
      /*if (hClusterSizeBins[i]->GetEntries() > 0) {
        hClusterSizeBins[i]->Scale(1.0 / hClusterSizeBins[i]->GetEntries()); 
      }*/

      hClusterSizeBins[i]->Draw("HIST TEXT");
      gPad->Update();
      hClusterSizeBins[i]->SetStats(1);
      gStyle->SetStatX(0.9);
    //  hClusterSizeBins[i]->SetStatW(0.15); 
    //  hClusterSizeBins[i]->SetStatH(0.2);  
      hClusterSizeBins[i]->SetFillColor(kBlue);
      hClusterSizeBins[i]->SetFillStyle(3001);
    //  hClusterSizeBins[i]->SetMaximum(1300.0);
    //  hClusterSizeBins[i]->Draw("HIST TEXT");

    }
   // hClusterSize->SetStats(1);
   // hClusterSize->SetFillStyle(3001);
   // hClusterSize->Draw("HIST TEXT");
    c3->SaveAs("cluster_size_per_bin.svg");
 
 
    std::fprintf(stdout, "Pixel level analysis complete.\n");
    return;
}
