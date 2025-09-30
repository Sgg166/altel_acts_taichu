void Hit_RawNum_test2(const std::string& rootFilePath){
    std::cout<< rootFilePath<<std::endl;
    TFile *tfile = new TFile(rootFilePath.c_str(),"READ");
    if(!tfile || !tfile->IsOpen()){
        std::fprintf(stderr, "tfile is not open\n");
        throw;
    }
    TTree *pTree = 0;
    tfile->GetObject("eventTree",pTree);
    if(!pTree){
        std::fprintf(stderr, "pTree is invalid\n");
        throw;
    }
    altel::TelEventTTreeReader ttreeReader;
    ttreeReader.setTTree(pTree);
 
    TH1F *hClaster_RawNum_det2 = new TH1F("hClaster_RawNum_det2", "The number of raw pixels contained in each cluster in detID 2", 30, 0, 15);
    TH1F *hClaster_RawNum_det3 = new TH1F("hClaster_RawNum_det3", "The number of raw pixels contained in each cluster in detID 3", 30, 0, 15);
    TH1F *hClaster_RawNum_det32 = new TH1F("hClaster_RawNum_det32", "The number of raw pixels contained in each cluster in detID 32",30,0, 15);
    TH1F *hClaster_RawNum_det5 = new TH1F("hClaster_RawNum_det5", "The number of raw pixels contained in each cluster in detID 5", 30, 0, 15);
    TH1F *hClaster_RawNum_det7 = new TH1F("hClaster_RawNum_det7", "The number of raw pixels contained in each cluster in detID 7", 30, 0, 15);
    TH1F *hClaster_RawNum_det9 = new TH1F("hClaster_RawNum_det9", "The number of raw pixels contained in each cluster in detID 9", 30, 0, 15);
 
    TH1F *hHitCount_det32 = new TH1F("hHitCount_det32", "Number of hits per event in detID 32", 30, 0, 15);
    
    // 新增：统计detID32中不同measHits.size()的事件数量
    // 假设最多20个hits，可以根据实际情况调整
    TH1F *hDet32HitCountDistribution = new TH1F("hDet32HitCountDistribution", "Distribution of events with different hit counts in detID 32", 21, -0.5, 20.5);
    hDet32HitCountDistribution->GetXaxis()->SetTitle("Number of Hits in detID 32");
    hDet32HitCountDistribution->GetYaxis()->SetTitle("Number of Events");
 
    size_t totalNumEvents = ttreeReader.numEvents();
    size_t total_det32_hits = 0;
    for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        std::fprintf(stdout, "Event: FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
        eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
 
        const auto& measHits = telEvent->measHits();
        size_t det32_event_hit_count = 0;
        if(!measHits.empty()){
            for(size_t i = 0; i < measHits.size(); i++){
                const auto& hit = measHits[i];
                unsigned int detID = hit->detN();
                size_t numRaws = hit->measRaws().size();
 
                if(detID == 2)
                {
                    hClaster_RawNum_det2->Fill(numRaws);
                    // std::cout<<" numRaws : "<<numRaws<<std::endl;
                }
                else if(detID == 3 ) hClaster_RawNum_det3 ->Fill(numRaws);
                else if(detID == 32)
                {
                    hClaster_RawNum_det32->Fill(numRaws);
                    // hHitCount_det32->Fill(det32_event_hit_count);
                    det32_event_hit_count++;
                    total_det32_hits++;
                }
                else if(detID == 5 ) hClaster_RawNum_det5 ->Fill(numRaws);
                else if(detID == 7 ) hClaster_RawNum_det7 ->Fill(numRaws);
                else if(detID == 9 ) hClaster_RawNum_det9 ->Fill(numRaws);
            }
        }
        if(det32_event_hit_count>0){
            hHitCount_det32->Fill(det32_event_hit_count);
            // 新增：统计不同hit count的事件数量
            hDet32HitCountDistribution->Fill(det32_event_hit_count);
        } else {
            // 统计0个hit的事件
            hDet32HitCountDistribution->Fill(0);
        }
    }
 
    TCanvas *c1 = new TCanvas("c1", "Raw Pixels per Hit Distribution 30.raw", 1600, 1200);
    c1->Divide(3, 4);
 
    c1->cd(1);
    hClaster_RawNum_det2->SetLineColor(kRed);
    hClaster_RawNum_det2->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det2 ->Scale(1.0 / hClaster_RawNum_det2 ->GetEntries());
    hClaster_RawNum_det2->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det2->SetMinimum(0);
    hClaster_RawNum_det2->SetMaximum(1.0);
    hClaster_RawNum_det2->Draw("HIST");
 
    c1->cd(2);
    hClaster_RawNum_det3->SetLineColor(kBlue);
    hClaster_RawNum_det3->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det3->Scale(1.0 / hClaster_RawNum_det3 ->GetEntries());
    hClaster_RawNum_det3->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det3->SetMinimum(0);
    hClaster_RawNum_det3->SetMaximum(1.0);
    hClaster_RawNum_det3->Draw("HIST");
 
    c1->cd(3);
    hClaster_RawNum_det32->SetLineColor(kViolet);
    hClaster_RawNum_det32->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det32->Scale(1.0 / hClaster_RawNum_det32->GetEntries());
    hClaster_RawNum_det32->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det32->SetMinimum(0);
    hClaster_RawNum_det32->SetMaximum(1.0);
    hClaster_RawNum_det32->Draw("HIST");
 
    c1->cd(7);
    hHitCount_det32->SetLineColor(kViolet);
    hHitCount_det32->GetXaxis()->SetTitle("Number of Hit count");
    hHitCount_det32->Scale(1.0 /hHitCount_det32 ->GetEntries());
    hHitCount_det32->GetYaxis()->SetTitle("Percentage");
    hHitCount_det32->SetMinimum(0);
    hHitCount_det32->SetMaximum(1.0);
    hHitCount_det32->Draw("HIST");
 
    c1->cd(4);
    hClaster_RawNum_det5->SetLineColor(kGreen);
    hClaster_RawNum_det5->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det5->Scale(1.0 / hClaster_RawNum_det5 ->GetEntries());
    hClaster_RawNum_det5->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det5->SetMinimum(0);
    hClaster_RawNum_det5->SetMaximum(1.0);
    hClaster_RawNum_det5->Draw("HIST");
 
    c1->cd(5);
    hClaster_RawNum_det7->SetLineColor(kMagenta);
    hClaster_RawNum_det7->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det7->Scale(1.0 / hClaster_RawNum_det7 ->GetEntries());
    hClaster_RawNum_det7->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det7->SetMinimum(0);
    hClaster_RawNum_det7->SetMaximum(1.0);
    hClaster_RawNum_det7->Draw("HIST");
 
    c1->cd(6);
    hClaster_RawNum_det9->SetLineColor(kOrange);
    hClaster_RawNum_det9->GetXaxis()->SetTitle("Number of Raw Pixels per Hit");
    hClaster_RawNum_det9->Scale(1.0 / hClaster_RawNum_det9->GetEntries());
    hClaster_RawNum_det9->GetYaxis()->SetTitle("Percentage");
    hClaster_RawNum_det9->SetMinimum(0);
    hClaster_RawNum_det9->SetMaximum(1.0);
    hClaster_RawNum_det9->Draw("HIST");
 
    // 新增：在画布的第8个位置绘制detID32 hit count分布的折线图
    c1->cd(8);
    hDet32HitCountDistribution->SetLineColor(kBlack);
    hDet32HitCountDistribution->SetLineWidth(2);
    hDet32HitCountDistribution->SetMarkerStyle(22);
    hDet32HitCountDistribution->SetMarkerSize(0.8);
    hDet32HitCountDistribution->Draw("LP"); // LP表示Line和Point，即折线图加数据点
 
    c1->SaveAs("Hit_Raw44.svg");
 
    TFile *outFile = new TFile("hits_analysis.root", "RECREATE");
    hClaster_RawNum_det2 ->Write();
    hClaster_RawNum_det3 ->Write();
    hClaster_RawNum_det32->Write();
    hClaster_RawNum_det5 ->Write();
    hClaster_RawNum_det7 ->Write();
    hClaster_RawNum_det9 ->Write();
    hHitCount_det32 ->Write();
    hDet32HitCountDistribution->Write(); // 保存新的histogram
    c1->Write();
    //outFile->Close();
 
    //tfile->Close();
    //delete tfile;
    std::cout << "Total hits for detID 32: " << total_det32_hits << std::endl;
    std::cout << "Average hits per event for detID 32: " << (double)total_det32_hits/totalNumEvents << std::endl;
    
    // 新增：输出detID32中不同hit count的事件数量统计
    std::cout << "\n=== detID32 Hit Count Distribution ===" << std::endl;
    for(int i = 0; i <= 20; i++){
        int count = hDet32HitCountDistribution->GetBinContent(i+1);
        if(count > 0){
            std::cout << "Events with " << i << " hits in detID32: " << count << std::endl;
        }
    }
    
    std::cout << "OVER,The data has been saved to hits_analysis.root and Hit_Raw44.svg" << std::endl;
    return;
}
