void Efficiency_pixel_level(const std::string& rootFilePath){
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
    float region_v_min = -5.5;
    float region_v_max = 7.3;
 
    int nBinsU = 26;
    int nBinsV = 13;
 
    size_t totalNumEvents = ttreeReader.numEvents();
    size_t totalTracks = 0;
    size_t goodTracks = 0;
    size_t totalMatchedHits = 0;
    size_t totalFitHits = 0;
    size_t totalMatchedTracks=0;
    size_t totalMatchedHits_inRegion = 0;
    size_t totalFitHits_inRegion = 0;
    size_t goodTracks_inRegion = 0;
    size_t totalMatchedHits_inPixel = 0;
    size_t totalFitHits_inPixel = 0;
 
    // 原始区域级效率直方图
    TEfficiency* pEff2D = new TEfficiency("eff2D","Efficiency vs Position;u [mm 26bin];v [mm 13bin];Efficiency",nBinsU, region_u_min, region_u_max, nBinsV, region_v_min, region_v_max);
    TEfficiency* pEffU = new TEfficiency("effU","Efficiency vs u;u [mm 26bin];Efficiency", nBinsU, region_u_min, region_u_max);
    TEfficiency* pEffV = new TEfficiency("effV","Efficiency vs v;v [mm 13bin];Efficiency", nBinsV, region_v_min, region_v_max);
 
    // 像素级效率直方图 (25um周期，分成5个bin)
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",
                                               5, -12.5, 12.5, 5, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 
                                              5, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 
                                              5, -12.5, 12.5);
 
    TH2I* hTrackCount = new TH2I("hTrackCount", "Track count per bin;u [mm 26bin];v [mm 13bin];Count",nBinsU, region_u_min, region_u_max,nBinsV, region_v_min, region_v_max);
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",
                                      5, -12.5, 12.5, 5, -12.5, 12.5);
 
    int minTracksPerBin = 0; // cut 阈值
 
    // 像素映射函数：将任意坐标映射到[-12.5, +12.5]um范围内
    auto mapToPixel = [](double coord_mm) {
        double coord_um = coord_mm * 1000.0; // 转换为um
        return coord_um - 25.0 * std::round(coord_um / 25.0); // 映射到[-12.5, +12.5]um
    };
 
    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
        for(auto aTraj: telEvent->trajs()){
            totalTracks++;
 
            bool isGoodTrack = (aTraj->numOriginMeasHit() == 5);
 
            if(isGoodTrack){
                goodTracks++;
                bool hasMatchedInDUT = false;
                bool hasFitInRegion = false;
 
                for(auto &aTrajHit: aTraj->trajHits()){
                    auto aFitHit = aTrajHit->fitHit();
                    auto aMatchedMeasHit = aTrajHit->matchedMeasHit();
 
                    if(aFitHit){
                        totalFitHits++;
                        if (aFitHit->detN() == 32){
                            if (aMatchedMeasHit) {
                              totalMatchedHits++;
                              hasMatchedInDUT = true;
                            }
                            double hit_u = aFitHit->u();
                            double hit_v = aFitHit->v();
 
                             bool isInRegion = (hit_u >= region_u_min && hit_u <= region_u_max &&
                                           hit_v >= region_v_min && hit_v <= region_v_max);
                             if(isInRegion){
                               goodTracks_inRegion++;
                               totalFitHits_inRegion++;
                               hasFitInRegion = true;
 
                               hTrackCount->Fill(hit_u, hit_v);
                               int binx = hTrackCount->GetXaxis()->FindBin(hit_u);
                               int biny = hTrackCount->GetYaxis()->FindBin(hit_v);
                               int nTracksInBin = hTrackCount->GetBinContent(binx, biny);
                               bool hasMatch = (aMatchedMeasHit != nullptr);
 
                               // 填充区域级直方图
                               if(nTracksInBin >= minTracksPerBin){
                                 pEff2D->Fill(hasMatch, hit_u, hit_v);
                                 pEffU->Fill (hasMatch, hit_u);
                                 pEffV->Fill (hasMatch, hit_v);
                               }
 
                               // 像素级分析
                               double u_mapped = mapToPixel(hit_u);
                               double v_mapped = mapToPixel(hit_v);
 
                               // 统计像素级track数
                               hTrackCount_pixel->Fill(u_mapped, v_mapped);
                               int binx_pixel = hTrackCount_pixel->GetXaxis()->FindBin(u_mapped);
                               int biny_pixel = hTrackCount_pixel->GetYaxis()->FindBin(v_mapped);
                               int nTracksInBin_pixel = hTrackCount_pixel->GetBinContent(binx_pixel, biny_pixel);
 
                               // 填充像素级直方图
                               if(nTracksInBin_pixel >= minTracksPerBin){
                                 pEff2D_pixel->Fill(hasMatch, u_mapped, v_mapped);
                                 pEffU_pixel->Fill(hasMatch, u_mapped);
                                 pEffV_pixel->Fill(hasMatch, v_mapped);
                               }
 
                               // 统计像素级数据
                               totalFitHits_inPixel++;
                               if(aMatchedMeasHit){
                                 totalMatchedHits_inPixel++;
                                 totalMatchedHits_inRegion++;
                               }
                             }
                             if (hasMatchedInDUT) {
                               totalMatchedTracks++;
                             }
                        }
                    }
                }
            }
        }
    }
 
    double efficiency = (goodTracks > 0) ?
                       (double)totalMatchedTracks / (double)goodTracks : 0.0;
    double local_efficiency = (goodTracks_inRegion > 0) ?
                             (double)totalMatchedHits_inRegion / (double)goodTracks_inRegion : 0.0;
    double pixel_efficiency = (totalFitHits_inPixel > 0) ?
                             (double)totalMatchedHits_inPixel / (double)totalFitHits_inPixel : 0.0;
 
    std::fprintf(stdout, "=== Overall efficiency result ===\n");
    std::fprintf(stdout, "Total events: %zu\n", totalNumEvents);
    std::fprintf(stdout, "Total tracks: %zu\n", totalTracks);
    std::fprintf(stdout, "Good tracks (==5 hits ): %zu\n", goodTracks);
    std::fprintf(stdout, "Total fit hits: %zu\n", totalFitHits);
    std::fprintf(stdout, "Total MatchedHits: %zu\n", totalMatchedHits);
    std::fprintf(stdout, "Total Matched Tracks: %zu\n", totalMatchedTracks);
    std::fprintf(stdout, "Overall Efficiency: %.3f\n", efficiency);
 
    std::fprintf(stdout, "\n=== Local efficiency result ===\n");
    std::fprintf(stdout, "Region bounds: u[%.1f, %.1f], v[%.1f, %.1f]\n",
                 region_u_min, region_u_max, region_v_min, region_v_max);
    std::fprintf(stdout, "Good tracks in region: %zu\n", goodTracks_inRegion);
    std::fprintf(stdout, "Total fit hits in region: %zu\n", totalFitHits_inRegion);
    std::fprintf(stdout, "Total MatchedHits in region: %zu\n", totalMatchedHits_inRegion);
    std::fprintf(stdout, "Local Efficiency: %.3f\n", local_efficiency);
 
    std::fprintf(stdout, "\n=== Pixel level efficiency result ===\n");
    std::fprintf(stdout, "Pixel size: 25 um\n");
    std::fprintf(stdout, "Pixel range: [-12.5, +12.5] um\n");
    std::fprintf(stdout, "Total fit hits in pixel: %zu\n", totalFitHits_inPixel);
    std::fprintf(stdout, "Total MatchedHits in pixel: %zu\n", totalMatchedHits_inPixel);
    std::fprintf(stdout, "Pixel Level Efficiency: %.3f\n", pixel_efficiency);
 
    // 创建画布：区域级效率
    TCanvas* c1 = new TCanvas("c1", "Regional Efficiency", 1200, 800);
    c1->Divide(2, 3);
    c1->cd(1);
    c1->SetFillStyle(1001);
    c1->SetFillColor(kWhite);
 
    pEff2D->Draw("TEXT colz");
    gStyle->SetPaintTextFormat(".2f");
    gPad->Update();
    pEff2D->SetMarkerStyle(20);
    pEff2D->SetMarkerSize(1);
    pEff2D->SetTitle("Regional Efficiency vs Position");
    pEff2D->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.8, 1);
 
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextColor(kRed);
    latex.DrawLatex(0.7, 0.85, Form("Regional Efficiency: %.3f", local_efficiency));
 
    c1->cd(3);
    pEffU->Draw("AP");
    gPad->Update();
    pEffU->SetTitle("Regional Efficiency vs u");
    pEffU->SetMarkerStyle(20);
    pEffU->SetMarkerSize(0.8);
    auto grU = pEffU->GetPaintedGraph();
    if (grU) {
      grU->GetYaxis()->SetRangeUser(0.8, 1);
    }
 
    c1->cd(4);
    pEffV->Draw("AP");
    gPad->Update();
    pEffV->SetTitle("Regional Efficiency vs v");
    pEffV->SetMarkerStyle(20);
    pEffV->SetMarkerSize(0.8);
    auto grV = pEffV->GetPaintedGraph();
    if (grV) {
      grV->GetYaxis()->SetRangeUser(0.8, 1);
    }
 
    c1->cd(2);
    hTrackCount->Draw("TEXT colz");
    hTrackCount->SetStats(0);
    hTrackCount->SetMarkerSize(1);
    hTrackCount->SetTitle("Regional Track count per bin");
 
    // 创建画布：像素级效率
    TCanvas* c2 = new TCanvas("c2", "Pixel Level Efficiency", 1200, 800);
    c2->Divide(2, 3);
    c2->cd(1);
    c2->SetFillStyle(1001);
    c2->SetFillColor(kWhite);
 
    pEff2D_pixel->Draw("TEXT colz");
    gStyle->SetPaintTextFormat(".2f");
    gPad->Update();
    pEff2D_pixel->SetMarkerStyle(20);
    pEff2D_pixel->SetMarkerSize(1);
    pEff2D_pixel->SetTitle("Pixel Efficiency vs Position");
    pEff2D_pixel->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.8, 1);
 
    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.05);
    latex2.SetTextColor(kRed);
    latex2.DrawLatex(0.7, 0.85, Form("Pixel Efficiency: %.3f", pixel_efficiency));
 
    c2->cd(3);
    pEffU_pixel->Draw("AP");
    gPad->Update();
    pEffU_pixel->SetTitle("Pixel Efficiency vs u");
    pEffU_pixel->SetMarkerStyle(20);
    pEffU_pixel->SetMarkerSize(0.8);
    auto grU_pixel = pEffU_pixel->GetPaintedGraph();
    if (grU_pixel) {
      grU_pixel->GetYaxis()->SetRangeUser(0.8, 1);
    }
 
    c2->cd(4);
    pEffV_pixel->Draw("AP");
    gPad->Update();
    pEffV_pixel->SetTitle("Pixel Efficiency vs v");
    pEffV_pixel->SetMarkerStyle(20);
    pEffV_pixel->SetMarkerSize(0.8);
    auto grV_pixel = pEffV_pixel->GetPaintedGraph();
    if (grV_pixel) {
      grV_pixel->GetYaxis()->SetRangeUser(0.8, 1);
    }
 
    c2->cd(2);
    hTrackCount_pixel->Draw("TEXT colz");
    hTrackCount_pixel->SetStats(0);
    hTrackCount_pixel->SetMarkerSize(1);
    hTrackCount_pixel->SetTitle("Pixel Track count per bin");
 
    // 保存结果
    TFile *outFile = new TFile("efficiency_analysis.root", "RECREATE");
    c1->SaveAs("regional_efficiency.svg");
    c2->SaveAs("pixel_efficiency.svg");
    std::fprintf(stdout, "Regional efficiency saved as 'regional_efficiency.svg'\n");
    std::fprintf(stdout, "Pixel efficiency saved as 'pixel_efficiency.svg'\n");
 
    std::fprintf(stdout, "Analysis complete.\n");
    return;
}
