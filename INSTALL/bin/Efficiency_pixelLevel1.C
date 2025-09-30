void Efficiency_pixelLevel1(const std::string& rootFilePath){
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
    
    // 创建像素级效率直方图
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",5, -12.5, 12.5, 5, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 5, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 5, -12.5,12.5);
 
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",5, -12.5, 12.5, 5, -12.5,12.5);
 
    TH2I* hNumMatchedDist = new TH2I("hNumMatchedDist", "NumMatched Distribution per pixel bin;u [um];v [um];NumMatched", 5, -12.5, 12.5, 5, -12.5, 12.5);
    TH1I* hNumMatched1D = new TH1I("hNumMatched1D","NumMatched Distribution;NumMatched;Count", 30, 0, 15);
 
    std::vector<TH1I*> hClusterSizeBins;
    for (int ix = 0; ix < 5; ix++) {
        for (int iy = 0; iy < 5; iy++) {
            TString name = Form("hClusterSize_bin_%d_%d", ix, iy);
            TString title = Form("Cluster Size Distribution (bin u=%d, v=%d);Cluster Size;Count", ix, iy);
            hClusterSizeBins.push_back(new TH1I(name, title, 30, 0, 15));
        }
    }
 
    // 新增：分析响应位置分布的直方图
    TH2I* hResponseOffset = new TH2I("hResponseOffset", "Response Offset;u offset [um];v offset [um];Count", 21, -10.5, 10.5, 21, -10.5, 10.5);
    TH1I* hResponseDistance = new TH1I("hResponseDistance", "Distance from Expected Position;Distance [um];Count", 50, 0, 25);
    TH2I* hResponseMap = new TH2I("hResponseMap", "Response Position Map;Expected u [um];Response u [um];Count", 5, -12.5, 12.5, 5, -12.5, 12.5);
    TH2I* hResponseMapV = new TH2I("hResponseMapV", "Response Position Map V;Expected v [um];Response v [um];Count", 5, -12.5, 12.5, 5, -12.5, 12.5);
    
    // 统计不同响应类型的计数器
    size_t samePixelResponse = 0;      // 响应在同一像素
    size_t adjacentPixelResponse = 0;   // 响应在相邻像素
    size_t distantPixelResponse = 0;    // 响应在远处像素
 
    // 坐标映射函数：将任意坐标映射到[-12.5um, +12.5um]范围
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
 
              // 映射到像素周期
              double pixel_u = mapToPixel(hit_u);
              double pixel_v = mapToPixel(hit_v);
 
              // 检查是否在像素范围内
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
 
                  // 新增：分析响应位置分布
                  if(numMatched > 0){
                    // 计算簇的平均位置（实际响应位置）
                    // 由于TelMeasRaw没有电荷信息，我们使用简单的几何平均
                    double cluster_u = 0.0;
                    double cluster_v = 0.0;
                    
                    for(auto& rawHit : aMatchedMeasHit->measRaws()){
                      // 获取原始像素坐标
                      uint16_t raw_u = rawHit.u();
                      uint16_t raw_v = rawHit.v();
                      
                      // 转换为毫米单位（假设像素间距为25um）
                      double raw_u_mm = raw_u * 0.025 - 0.025 * (1024/2. - 0.5);
                      double raw_v_mm = raw_v * 0.025 - 0.025 * (512/2. - 0.5);
                      
                      cluster_u += raw_u_mm;
                      cluster_v += raw_v_mm;
                    }
                    
                    // 计算几何平均位置
                    cluster_u /= numMatched;
                    cluster_v /= numMatched;
                    
                    // 转换为微米单位
                    cluster_u *= 1000.0;
                    cluster_v *= 1000.0;
                    
                    // 映射到像素坐标系
                    double cluster_pixel_u = mapToPixel(cluster_u / 1000.0);
                    double cluster_pixel_v = mapToPixel(cluster_v / 1000.0);
                    
                    // 计算偏移量
                    double offset_u = cluster_pixel_u - pixel_u;
                    double offset_v = cluster_pixel_v - pixel_v;
                    double distance = std::sqrt(offset_u * offset_u + offset_v * offset_v);
                    
                    // 填充响应分布直方图
                    hResponseOffset->Fill(offset_u, offset_v);
                    hResponseDistance->Fill(distance);
                    hResponseMap->Fill(pixel_u, cluster_pixel_u);
                    hResponseMapV->Fill(pixel_v, cluster_pixel_v);
                    
                    // 统计响应类型
                    if(distance < 2.5){ // 同一像素内
                      samePixelResponse++;
                    } else if(distance < 7.5){ // 相邻像素
                      adjacentPixelResponse++;
                    } else { // 远处像素
                      distantPixelResponse++;
                    }
                  }
 
                  int binx = int((pixel_u + 12.5) / 5.0);
                  int biny = int((pixel_v + 12.5) / 5.0);
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
    
    // 新增：输出响应分布统计
    std::fprintf(stdout, "\n=== Response Position Analysis ===\n");
    std::fprintf(stdout, "Same pixel response (< 2.5 um): %zu (%.1f%%)\n", 
                samePixelResponse, 100.0 * samePixelResponse / totalMatchedHits_inPixel);
    std::fprintf(stdout, "Adjacent pixel response (2.5-7.5 um): %zu (%.1f%%)\n", 
                adjacentPixelResponse, 100.0 * adjacentPixelResponse / totalMatchedHits_inPixel);
    std::fprintf(stdout, "Distant pixel response (> 7.5 um): %zu (%.1f%%)\n", 
                distantPixelResponse, 100.0 * distantPixelResponse / totalMatchedHits_inPixel);
 
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
      hClusterSizeBins[i]->Draw("HIST TEXT");
      gPad->Update();
      hClusterSizeBins[i]->SetStats(1);
      gStyle->SetStatX(0.9);
      hClusterSizeBins[i]->SetFillColor(kBlue);
      hClusterSizeBins[i]->SetFillStyle(3001);
      hClusterSizeBins[i]->SetMaximum(1300.0);
    }
    c3->SaveAs("cluster_size_per_bin.svg");
 
    // 新增：创建响应分布分析画布
    TCanvas* c4 = new TCanvas("c4", "Response Position Analysis", 1500, 1000);
    c4->Divide(3, 2);
 
    c4->cd(1);
    hResponseOffset->Draw("TEXT colz");
    hResponseOffset->SetStats(0);
    hResponseOffset->SetTitle("Response Offset Distribution");
    gStyle->SetPaintTextFormat(".0f");
    gPad->SetFixedAspectRatio();
 
    c4->cd(2);
    hResponseDistance->Draw();
    hResponseDistance->SetTitle("Distance from Expected Position");
    hResponseDistance->SetStats(1);
    hResponseDistance->SetFillColor(kBlue);
    hResponseDistance->SetFillStyle(3001);
 
    c4->cd(3);
    hResponseMap->Draw("TEXT colz");
    hResponseMap->SetStats(0);
    hResponseMap->SetTitle("Response Map (u direction)");
    gStyle->SetPaintTextFormat(".0f");
    gPad->SetFixedAspectRatio();
 
    c4->cd(4);
    hResponseMapV->Draw("TEXT colz");
    hResponseMapV->SetStats(0);
    hResponseMapV->SetTitle("Response Map (v direction)");
    gStyle->SetPaintTextFormat(".0f");
    gPad->SetFixedAspectRatio();
 
    c4->cd(5);
    // 创建饼图显示响应类型分布
    TH1F* hResponseType = new TH1F("hResponseType", "Response Type Distribution", 3, 0, 3);
    hResponseType->SetBinContent(1, samePixelResponse);
    hResponseType->SetBinContent(2, adjacentPixelResponse);
    hResponseType->SetBinContent(3, distantPixelResponse);
    hResponseType->GetXaxis()->SetBinLabel(1, "Same Pixel");
    hResponseType->GetXaxis()->SetBinLabel(2, "Adjacent");
    hResponseType->GetXaxis()->SetBinLabel(3, "Distant");
    hResponseType->SetBarWidth(0.8);
    hResponseType->SetFillColor(kBlue);
    hResponseType->SetStats(0);
    hResponseType->Draw("BAR");
 
    c4->SaveAs("response_position_analysis.svg");
    std::fprintf(stdout, "Response position analysis plots saved as 'response_position_analysis.svg'\n");
 
    std::fprintf(stdout, "Pixel level analysis complete.\n");
    return;
}
