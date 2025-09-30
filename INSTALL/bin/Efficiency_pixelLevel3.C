void Efficiency_pixelLevel3(const std::string& rootFilePath){
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
    
    // 新增：统计响应位置分布的变量
    size_t centerPixelResponses = 0;      // 中心像素响应
    size_t neighborPixelResponses = 0;     // 周边像素响应
    size_t totalResponses = 0;             // 总响应数
 
    // 新增：存储响应坐标的容器
    std::vector<std::pair<double, double>> fitPositions;    // 拟合位置 (u, v) in um
    std::vector<std::pair<double, double>> measPositions;   // 测量位置 (u, v) in um
    std::vector<double> responseDistances;                  // 响应距离
    std::vector<int> responseTypes;                         // 响应类型: 0=center, 1=neighbor, 2=far
 
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
 
    // 新增：坐标分布散点图
    TH2D* hFitPositions = new TH2D("hFitPositions", "Fit Hit Positions;u [mm];v [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
    TH2D* hMeasPositions = new TH2D("hMeasPositions", "Measured Hit Positions;u [mm];v [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
    TH2D* hFitVsMeas = new TH2D("hFitVsMeas", "Fit vs Measured Positions;Fit u [mm];Meas u [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
 
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
 
                  // 新增：分析响应位置分布
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
                  
                  // 填充偏移分布直方图
                  hResponseOffset->Fill(delta_u, delta_v);
                  hResponseDistance->Fill(distance);
                  
                  // 填充坐标散点图
                  hFitPositions->Fill(fit_u/1000.0, fit_v/1000.0);
                  hMeasPositions->Fill(meas_u/1000.0, meas_v/1000.0);
                  hFitVsMeas->Fill(fit_u/1000.0, meas_u/1000.0);
                  
                  // 判断响应类型
                  totalResponses++;
                  int responseType = 2; // 默认为远距离
                  if (distance <= 12.5) {
                      // 中心像素响应（距离小于等于像素半宽）
                      centerPixelResponses++;
                      responseType = 0; // Center
                  } else if (distance <= 25.0) {
                      // 周边像素响应（距离在1-2个像素宽度之间）
                      neighborPixelResponses++;
                      responseType = 1; // Neighbor
                  }
                  
               
                  fitPositions.push_back(std::make_pair(fit_u, fit_v));
                  measPositions.push_back(std::make_pair(meas_u, meas_v));
                  responseDistances.push_back(distance);
                  responseTypes.push_back(responseType);
                  
                  hResponseType->Fill(responseType);
 
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
 
    // 新增：输出响应位置统计结果
    std::fprintf(stdout, "\n=== Response Position Analysis ===\n");
    std::fprintf(stdout, "Total responses analyzed: %zu\n", totalResponses);
    std::fprintf(stdout, "Center pixel responses: %zu (%.2f%%)\n", 
                centerPixelResponses, totalResponses > 0 ? (double)centerPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Neighbor pixel responses: %zu (%.2f%%)\n", 
                neighborPixelResponses, totalResponses > 0 ? (double)neighborPixelResponses/totalResponses*100 : 0);
    std::fprintf(stdout, "Far responses: %zu (%.2f%%)\n", 
                totalResponses - centerPixelResponses - neighborPixelResponses, 
                totalResponses > 0 ? (double)(totalResponses - centerPixelResponses - neighborPixelResponses)/totalResponses*100 : 0);
 
    // 新增：输出具体坐标信息（前20个作为示例）
    std::fprintf(stdout, "\n=== Sample Response Coordinates (first 20) ===\n");
    std::fprintf(stdout, "Index\tFit(u,v)\t\tMeas(u,v)\t\tDistance\tType\n");
    std::fprintf(stdout, "----------------------------------------------------------------------------------------\n");
    for (size_t i = 0; i < std::min(size_t(20), fitPositions.size()); i++) {
        std::string typeStr = (responseTypes[i] == 0) ? "Center" : 
                             (responseTypes[i] == 1) ? "Neighbor" : "Far";
        std::fprintf(stdout, "%4zu\t(%.1f,%.1f)\t(%.1f,%.1f)\t%.2f\t%s\n", 
                    i, 
                    fitPositions[i].first, fitPositions[i].second,
                    measPositions[i].first, measPositions[i].second,
                    responseDistances[i], typeStr.c_str());
    }
 
    // 新增：计算并输出坐标统计信息
    if (!fitPositions.empty()) {
        std::fprintf(stdout, "\n=== Coordinate Statistics ===\n");
        
        // 计算拟合位置统计
        double fit_u_mean = 0, fit_v_mean = 0;
        double meas_u_mean = 0, meas_v_mean = 0;
        double dist_mean = 0;
        
        for (size_t i = 0; i < fitPositions.size(); i++) {
            fit_u_mean += fitPositions[i].first;
            fit_v_mean += fitPositions[i].second;
            meas_u_mean += measPositions[i].first;
            meas_v_mean += measPositions[i].second;
            dist_mean += responseDistances[i];
        }
        
        fit_u_mean /= fitPositions.size();
        fit_v_mean /= fitPositions.size();
        meas_u_mean /= measPositions.size();
        meas_v_mean /= measPositions.size();
        dist_mean /= responseDistances.size();
        
        std::fprintf(stdout, "Fit positions mean: (%.2f, %.2f) um\n", fit_u_mean, fit_v_mean);
        std::fprintf(stdout, "Meas positions mean: (%.2f, %.2f) um\n", meas_u_mean, meas_v_mean);
        std::fprintf(stdout, "Mean response distance: %.2f um\n", dist_mean);
        
        // 计算标准差
        double fit_u_std = 0, fit_v_std = 0;
        double meas_u_std = 0, meas_v_std = 0;
        double dist_std = 0;
         
        for (size_t i = 0; i < fitPositions.size(); i++) {
            fit_u_std += std::pow(fitPositions[i].first - fit_u_mean, 2);
            fit_v_std += std::pow(fitPositions[i].second - fit_v_mean, 2);
            meas_u_std += std::pow(measPositions[i].first - meas_u_mean, 2);
            meas_v_std += std::pow(measPositions[i].second - meas_v_mean, 2);
            dist_std += std::pow(responseDistances[i] - dist_mean, 2);
        }
        
        fit_u_std = std::sqrt(fit_u_std / fitPositions.size());
        fit_v_std = std::sqrt(fit_v_std / fitPositions.size());
        meas_u_std = std::sqrt(meas_u_std / measPositions.size());
        meas_v_std = std::sqrt(meas_v_std / measPositions.size());
        dist_std = std::sqrt(dist_std / responseDistances.size());
        
        std::fprintf(stdout, "Fit positions std: (%.2f, %.2f) um\n", fit_u_std, fit_v_std);
        std::fprintf(stdout, "Meas positions std: (%.2f, %.2f) um\n", meas_u_std, meas_v_std);
        std::fprintf(stdout, "Response distance std: %.2f um\n", dist_std);
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
    pEff2D_pixel->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.9, 1);
 
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
    double center_frac = totalResponses > 0 ? (double)centerPixelResponses / totalResponses : 0;
    double neighbor_frac = totalResponses > 0 ? (double)neighborPixelResponses / totalResponses : 0;
    double far_frac = totalResponses > 0 ? (double)(totalResponses - centerPixelResponses - neighborPixelResponses) / totalResponses : 0;
    
    double values[3] = {center_frac * 100, neighbor_frac * 100, far_frac * 100};
    int colors[3] = {kBlue, kGreen, kRed};
    TPie* pie = new TPie("responsePie", "Response Type Distribution", 3, values, colors);
    
    pie->SetEntryLabel(0, "Center");
    pie->SetEntryLabel(1, "Neighbor");
    pie->SetEntryLabel(2, "Far");
    
    pie->SetRadius(0.3);
    pie->Draw("rsc");
    
    c4->SaveAs("response_position_analysis.svg");
    std::fprintf(stdout, "Response position analysis plots saved as 'response_position_analysis.svg'\n");
 
    // 创建坐标分布画布
    TCanvas* c5 = new TCanvas("c5", "Coordinate Distribution Analysis", 1500, 1000);
    c5->Divide(3, 2);
 
    c5->cd(1);
    hFitPositions->Draw("P");
    hFitPositions->SetTitle("Fit Hit Positions");
    hFitPositions->GetXaxis()->SetTitle("u [um]");
    hFitPositions->GetYaxis()->SetTitle("v [um]");
    hFitPositions->SetMarkerStyle(20);
    hFitPositions->SetMarkerSize(0.5);
 
    c5->cd(2);
    hMeasPositions->Draw("P");
    hMeasPositions->SetTitle("Measured Hit Positions");
    hMeasPositions->GetXaxis()->SetTitle("u [um]");
    hMeasPositions->GetYaxis()->SetTitle("v [um]");
    hMeasPositions->SetMarkerStyle(20);
    hMeasPositions->SetMarkerSize(0.5);
 
    c5->cd(3);
    hFitVsMeas->Draw("P");
    hFitVsMeas->SetTitle("Fit vs Measured u-coordinates");
    hFitVsMeas->GetXaxis()->SetTitle("Fit u [um]");
    hFitVsMeas->GetYaxis()->SetTitle("Meas u [um]");
    hFitVsMeas->SetMarkerStyle(20);
    hFitVsMeas->SetMarkerSize(0.5);
    
    // 添加对角线参考线
    TLine* diagonal = new TLine(-15, -15, 15, 15);
    diagonal->SetLineColor(kRed);
    diagonal->SetLineStyle(2);
    diagonal->Draw("same");
 
    c5->cd(4);
    // 创建按响应类型分类的散点图
    TH2D* hFitCenter = new TH2D("hFitCenter", "Fit Positions - Center Responses;u [mm];v [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
    TH2D* hFitNeighbor = new TH2D("hFitNeighbor", "Fit Positions - Neighbor Responses;u [mm];v [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
    TH2D* hFitFar = new TH2D("hFitFar", "Fit Positions - Far Responses;u [mm];v [mm]", 26, -11.6, 14, 13, -5.5, 7.3);
    
    for (size_t i = 0; i < fitPositions.size(); i++) {
        if (responseTypes[i] == 0) {
            hFitCenter->Fill(fitPositions[i].first, fitPositions[i].second);
        } else if (responseTypes[i] == 1) {
            hFitNeighbor->Fill(fitPositions[i].first, fitPositions[i].second);
        } else {
            hFitFar->Fill(fitPositions[i].first, fitPositions[i].second);
        }
    }
    
    hFitCenter->Draw("P");
    hFitCenter->SetMarkerColor(kBlue);
    hFitCenter->SetMarkerStyle(20);
    hFitCenter->SetMarkerSize(0.8);
    hFitCenter->SetTitle("Fit Positions by Response Type");
 
    c5->cd(5);
    hFitNeighbor->Draw("P");
    hFitNeighbor->SetMarkerColor(kGreen);
    hFitNeighbor->SetMarkerStyle(20);
    hFitNeighbor->SetMarkerSize(0.8);
 
    c5->cd(6);
    hFitFar->Draw("P");
    hFitFar->SetMarkerColor(kRed);
    hFitFar->SetMarkerStyle(20);
    hFitFar->SetMarkerSize(0.8);
 
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
