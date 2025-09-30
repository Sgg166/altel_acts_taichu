void Efficiency_pixelLevel4(const std::string& rootFilePath){
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
 
    // 统计响应位置分布
    size_t centerPixelResponses = 0;      // 中心像素响应
    size_t neighborPixelResponses = 0;     // 周边像素响应
    size_t totalResponses = 0;             // 总响应数
 
    // 新增：存储像素响应具体坐标的容器
    std::vector<std::tuple<double, double, double, double, double, int>> responseCoordinates;
    // 格式：fit_u, fit_v, meas_u, meas_v, distance, response_type
 
    //Create pixel-level efficiency histograms
    TEfficiency* pEff2D_pixel = new TEfficiency("eff2D_pixel","Pixel Efficiency;u [um];v [um];Efficiency",5, -12.5, 12.5, 5, -12.5, 12.5);
    TEfficiency* pEffU_pixel = new TEfficiency("effU_pixel","Pixel Efficiency vs u;u [um];Efficiency", 5, -12.5, 12.5);
    TEfficiency* pEffV_pixel = new TEfficiency("effV_pixel","Pixel Efficiency vs v;v [um];Efficiency", 5, -12.5,12.5);
 
    TH2I* hTrackCount_pixel = new TH2I("hTrackCount_pixel", "Track count per pixel bin;u [um];v [um];Count",5, -12.5, 12.5, 5, -12.5,12.5);
 
    TH2I* hNumMatchedDist = new TH2I("hNumMatchedDist", "NumMatched Distribution per pixel bin;u [um];v [um];NumMatched", 5, -12.5, 12.5, 5, -12.5, 12.5);
    TH1I* hNumMatched1D = new TH1I("hNumMatched1D","NumMatched Distribution;NumMatched;Count", 30, 0, 15);
 
    // 新增：响应位置分布直方图
    TH2I* hResponseOffset = new TH2I("hResponseOffset", "Response Offset Distribution;#Delta u [um];#Delta v [um];Count",
                                    50, -25, 25, 50, -25, 25);
    TH1I* hResponseDistance = new TH1I("hResponseDistance", "Response Distance from Center;Distance [um];Count",
                                      50, 0, 35);
    TH1I* hResponseType = new TH1I("hResponseType", "Response Type;Type;Count", 3, 0, 3);
    // 设置响应类型标签
    hResponseType->GetXaxis()->SetBinLabel(1, "Center");
    hResponseType->GetXaxis()->SetBinLabel(2, "Neighbor");
    hResponseType->GetXaxis()->SetBinLabel(3, "Far");
 
    // 新增：坐标分布直方图
    TH2I* hFitCoordinates = new TH2I("hFitCoordinates", "Fit Hit Coordinates;fit_u [um];fit_v [um];Count",
                                   100, -15, 15, 100, -15, 15);
    TH2I* hMeasCoordinates = new TH2I("hMeasCoordinates", "Measured Hit Coordinates;meas_u [um];meas_v [um];Count",
                                    100, -15, 15, 100, -15, 15);
    TH2I* hCoordinateComparison = new TH2I("hCoordinateComparison", "Measured vs Fit Coordinates;fit_u [um];meas_u [um];Count",
                                          100, -15, 15, 100, -15, 15);
 
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
                  hNumMatchedDist->Fill(pixel_u, pixel_v, numMatched);
                  hNumMatched1D->Fill(numMatched);
 
                  // 分析响应位置分布
                  // 获取拟合击中位置（预期位置）
                  double fit_u = aFitHit->u() * 1000.0; // 转换为微米
                  double fit_v = aFitHit->v() * 1000.0;
 
                  // 获取测量击中位置（实际响应位置）
                  double meas_u = aMatchedMeasHit->u() * 1000.0;
                  double meas_v = aMatchedMeasHit->v() * 1000.0;
 
                  // 计算偏移量
                  double delta_u = meas_u - fit_u;
                  double delta_v = meas_v - fit_v;
                  double distance = std::sqrt(delta_u * delta_u + delta_v * delta_v);
 
                  // 填充坐标分布直方图
                  hFitCoordinates->Fill(fit_u, fit_v);
                  hMeasCoordinates->Fill(meas_u, meas_v);
                  hCoordinateComparison->Fill(fit_u, meas_u);
 
                  // 填充偏移分布直方图
                  hResponseOffset->Fill(delta_u, delta_v);
                  hResponseDistance->Fill(distance);
 
                  // 判断响应类型
                  totalResponses++;
                  int responseType = -1;
                  if (distance <= 12.5) {
                      // 中心像素响应（距离小于等于像素半宽）
                      centerPixelResponses++;
                      hResponseType->Fill(0); // Center
                      responseType = 0;
                  } else if (distance <= 25.0) {
                      // 周边像素响应（距离在1-2个像素宽度之间）
                      neighborPixelResponses++;
                      hResponseType->Fill(1); // Neighbor
                      responseType = 1;
                  } else {
                      // 远距离响应
                      hResponseType->Fill(2); // Far
                      responseType = 2;
                  }
 
                  // 存储坐标信息
                  responseCoordinates.emplace_back(fit_u, fit_v, meas_u, meas_v, distance, responseType);
 
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
 
    // 输出响应位置统计结果
    std::fprintf(stdout, "\n=== Response Position Analysis ===\n");
    std::fprintf(stdout, "Total responses analyzed: %zu\n", totalResponses);
    std::fprintf(stdout, "Center pixel responses: %zu (%.2f%%)\n",
                centerPixelResponses, totalResponses > 0 ? (double)centerPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Neighbor pixel responses: %zu (%.2f%%)\n",
                neighborPixelResponses, totalResponses > 0 ? (double)neighborPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Far responses: %zu (%.2f%%)\n",
                totalResponses - centerPixelResponses - neighborPixelResponses,
                totalResponses > 0 ? (double)(totalResponses - centerPixelResponses - neighborPixelResponses)/totalResponses*100 : 0);
 
    // 新增：输出坐标信息到文件
    std::ofstream coordFile("pixel_response_coordinates.txt");
    if (coordFile.is_open()) {
        coordFile << "# 像素响应坐标数据\n";
        coordFile << "# 格式: fit_u(um) fit_v(um) meas_u(um) meas_v(um) distance(um) response_type\n";
        coordFile << "# response_type: 0=Center, 1=Neighbor, 2=Far\n";
        coordFile << "# 总响应数: " << responseCoordinates.size() << "\n\n";
        
        for (const auto& coord : responseCoordinates) {
            coordFile << std::get<0>(coord) << " " << std::get<1>(coord) << " "
                     << std::get<2>(coord) << " " << std::get<3>(coord) << " "
                     << std::get<4>(coord) << " " << std::get<5>(coord) << "\n";
        }
        coordFile.close();
        std::fprintf(stdout, "坐标数据已保存到 'pixel_response_coordinates.txt'\n");
    } else {
        std::fprintf(stderr, "无法创建坐标输出文件\n");
    }
 
    // 新增：在控制台输出前10个响应的坐标信息
    std::fprintf(stdout, "\n=== 前10个像素响应坐标示例 ===\n");
    std::fprintf(stdout, "序号\t拟合位置(um)\t\t\t测量位置(um)\t\t\t距离(um)\t类型\n");
    std::fprintf(stdout, "    \t(u, v)\t\t\t\t(u, v)\t\t\t\t    \t\n");
    for (size_t i = 0; i < std::min(size_t(10), responseCoordinates.size()); ++i) {
        const auto& coord = responseCoordinates[i];
        std::string typeStr = (std::get<5>(coord) == 0) ? "Center" : 
                             (std::get<5>(coord) == 1) ? "Neighbor" : "Far";
        std::fprintf(stdout, "%zu\t(%.2f, %.2f)\t\t(%.2f, %.2f)\t\t%.2f\t%s\n", 
                    i+1, std::get<0>(coord), std::get<1>(coord), 
                    std::get<2>(coord), std::get<3>(coord), 
                    std::get<4>(coord), typeStr.c_str());
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
 
    // 创建响应位置分析画布
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
 
    c4->cd(4);
    // 显示坐标对比
    hCoordinateComparison->Draw("TEXT colz");
    hCoordinateComparison->SetStats(0);
    hCoordinateComparison->SetTitle("Measured vs Fit Coordinates");
    hCoordinateComparison->GetXaxis()->SetTitle("fit_u [um]");
    hCoordinateComparison->GetYaxis()->SetTitle("meas_u [um]");
    gPad->SetFixedAspectRatio();
 
    c4->SaveAs("response_position_analysis.svg");
    std::fprintf(stdout, "Response position analysis plots saved as 'response_position_analysis.svg'\n");
 
    // 新增：创建坐标分布画布
    TCanvas* c5 = new TCanvas("c5", "Coordinate Distribution Analysis", 1200, 500);
    c5->Divide(2, 1);
 
    c5->cd(1);
    hFitCoordinates->Draw("TEXT colz");
    hFitCoordinates->SetStats(0);
    hFitCoordinates->SetTitle("Fit Hit Coordinates Distribution");
    hFitCoordinates->GetXaxis()->SetTitle("fit_u [um]");
    hFitCoordinates->GetYaxis()->SetTitle("fit_v [um]");
    gPad->SetFixedAspectRatio();
 
    c5->cd(2);
    hMeasCoordinates->Draw("TEXT colz");
    hMeasCoordinates->SetStats(0);
    hMeasCoordinates->SetTitle("Measured Hit Coordinates Distribution");
    hMeasCoordinates->GetXaxis()->SetTitle("meas_u [um]");
    hMeasCoordinates->GetYaxis()->SetTitle("meas_v [um]");
    gPad->SetFixedAspectRatio();
 
    c5->SaveAs("coordinate_distribution.svg");
    std::fprintf(stdout, "Coordinate distribution plots saved as 'coordinate_distribution.svg'\n");
 
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
 
    std::fprintf(stdout, "Pixel level analysis complete.\n");
    return;
}
