void trajHit_per_detector(const std::string& rootFilePath){
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
 
    // 创建每个探测器的TrajHit数量直方图
    TH1F *h_trajHit_det2 = new TH1F("h_trajHit_det2", "TrajHit Distribution det2", 20, 0, 20);
    TH1F *h_trajHit_det3 = new TH1F("h_trajHit_det3", "TrajHit Distribution det3", 20, 0, 20);
    TH1F *h_trajHit_det32 = new TH1F("h_trajHit_det32", "TrajHit Distribution det32", 20, 0, 20);
    TH1F *h_trajHit_det5 = new TH1F("h_trajHit_det5", "TrajHit Distribution det5", 20, 0, 20);
    TH1F *h_trajHit_det7 = new TH1F("h_trajHit_det7", "TrajHit Distribution det7", 20, 0, 20);
    TH1F *h_trajHit_det9 = new TH1F("h_trajHit_det9", "TrajHit Distribution det9", 20, 0, 20);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    size_t totalEvents = 0;
    size_t totalTrajHits = 0;
 
    // 创建每个探测器的事件计数器
    std::map<int, size_t> detectorEvents;
    std::map<int, size_t> detectorTrajHits;
 
    for(size_t eventNum = 0; eventNum < totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        std::fprintf(stdout, "Event: FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
                     eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
        
        totalEvents++;
        
        // 为每个事件创建探测器标记
        std::map<int, int> eventDetectorTrajHits;
        
        if (!telEvent->trajs().empty()) {
            for (size_t i = 0; i < telEvent->trajs().size(); i++) {
                const auto& traj = telEvent->trajs()[i];
                if (!traj) {
                    std::fprintf(stderr, "Warning: Trajectory %zu is null\n", i);
                    continue;
                }
                
                for (size_t j = 0; j < traj->numTrajHit(); j++) {
                    auto trajHit = traj->trajHit(j);
                    if (!trajHit){
                        continue;
                    }
                    
                    totalTrajHits++;
                    int detID = trajHit->detN();
                    
                    // 统计每个探测器的TrajHit数量
                    eventDetectorTrajHits[detID]++;
                    detectorTrajHits[detID]++;
                    
                    // 根据探测器ID填充相应的直方图
                    switch(detID) {
                        case 2:
                            h_trajHit_det2->Fill(eventDetectorTrajHits[detID]);
                            break;
                        case 3:
                            h_trajHit_det3->Fill(eventDetectorTrajHits[detID]);
                            break;
                        case 32:
                            h_trajHit_det32->Fill(eventDetectorTrajHits[detID]);
                            break;
                        case 5:
                            h_trajHit_det5->Fill(eventDetectorTrajHits[detID]);
                            break;
                        case 7:
                            h_trajHit_det7->Fill(eventDetectorTrajHits[detID]);
                            break;
                        case 9:
                            h_trajHit_det9->Fill(eventDetectorTrajHits[detID]);
                            break;
                    }
                }
            }
            
            // 统计每个探测器是否有事件
            for(const auto& pair : eventDetectorTrajHits) {
                detectorEvents[pair.first]++;
            }
        }
    }
 
    // 创建画布
    TCanvas *c1 = new TCanvas("c1", "TrajHit per Detector Distribution", 1200, 900);
    c1->Divide(3, 3);
 
    // 绘制每个探测器的直方图
    c1->cd(1);
    h_trajHit_det2->SetLineColor(kRed);
    h_trajHit_det2->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det2->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det2->GetEntries() > 0){
        h_trajHit_det2->Scale(1.0 / h_trajHit_det2->GetEntries());
        h_trajHit_det2->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det2->SetMinimum(0);
    h_trajHit_det2->Draw("HIST");
 
    c1->cd(2);
    h_trajHit_det3->SetLineColor(kBlue);
    h_trajHit_det3->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det3->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det3->GetEntries() > 0){
        h_trajHit_det3->Scale(1.0 / h_trajHit_det3->GetEntries());
        h_trajHit_det3->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det3->SetMinimum(0);
    h_trajHit_det3->Draw("HIST");
 
    c1->cd(3);
    h_trajHit_det32->SetLineColor(kViolet);
    h_trajHit_det32->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det32->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det32->GetEntries() > 0){
        h_trajHit_det32->Scale(1.0 / h_trajHit_det32->GetEntries());
        h_trajHit_det32->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det32->SetMinimum(0);
    h_trajHit_det32->Draw("HIST");
 
    c1->cd(4);
    h_trajHit_det5->SetLineColor(kGreen);
    h_trajHit_det5->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det5->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det5->GetEntries() > 0){
        h_trajHit_det5->Scale(1.0 / h_trajHit_det5->GetEntries());
        h_trajHit_det5->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det5->SetMinimum(0);
    h_trajHit_det5->Draw("HIST");
 
    c1->cd(5);
    h_trajHit_det7->SetLineColor(kMagenta);
    h_trajHit_det7->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det7->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det7->GetEntries() > 0){
        h_trajHit_det7->Scale(1.0 / h_trajHit_det7->GetEntries());
        h_trajHit_det7->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det7->SetMinimum(0);
    h_trajHit_det7->Draw("HIST");
 
    c1->cd(6);
    h_trajHit_det9->SetLineColor(kOrange);
    h_trajHit_det9->GetXaxis()->SetTitle("Number of TrajHits per Event");
    h_trajHit_det9->GetYaxis()->SetTitle("Events");
    if(h_trajHit_det9->GetEntries() > 0){
        h_trajHit_det9->Scale(1.0 / h_trajHit_det9->GetEntries());
        h_trajHit_det9->GetYaxis()->SetTitle("Fraction");
    }
    h_trajHit_det9->SetMinimum(0);
    h_trajHit_det9->Draw("HIST");
 
    c1->SaveAs("trajHit_per_detector.svg");
 
    // 保存结果到ROOT文件
    TFile *outFile = new TFile("trajHit_analysis.root", "RECREATE");
    h_trajHit_det2->Write();
    h_trajHit_det3->Write();
    h_trajHit_det32->Write();
    h_trajHit_det5->Write();
    h_trajHit_det7->Write();
    h_trajHit_det9->Write();
    c1->Write();
 
    // 输出统计信息
    std::fprintf(stdout, "\n=== TRAJHIT STATISTICS ===\n");
    std::fprintf(stdout, "Total events: %zu\n", totalEvents);
    std::fprintf(stdout, "Total trajectory hits: %zu\n", totalTrajHits);
    std::fprintf(stdout, "\n=== PER DETECTOR STATISTICS ===\n");
    
    for(const auto& pair : detectorEvents) {
        int detID = pair.first;
        size_t events = pair.second;
        size_t trajHits = detectorTrajHits[detID];
        double avgTrajHitsPerEvent = events > 0 ? static_cast<double>(trajHits) / events : 0.0;
        
        std::fprintf(stdout, "Detector %d: %zu events, %zu total TrajHits, %.2f avg TrajHits per event\n", 
                     detID, events, trajHits, avgTrajHitsPerEvent);
    }
    
    return;
}
