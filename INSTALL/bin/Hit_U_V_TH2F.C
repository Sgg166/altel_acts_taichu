void Hit_U_V_TH2F(const std::string& rootFilePath){
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
 
    TH2F *h2_det2 = new TH2F("h2_det2", "Hit Position (u,v) in detID 2", 100, -14, 14, 100, -8, 8);
    TH2F *h2_det3 = new TH2F("h2_det3", "Hit Position (u,v) in detID 3", 100, -14, 14, 100, -8, 8);
    TH2F *h2_det32 = new TH2F("h2_det32", "Hit Position (u,v) in detID 32",100, -14, 14, 100, -8, 8);
    TH2F *h2_det5 = new TH2F("h2_det5", "Hit Position (u,v) in detID 5", 100, -14, 14, 100, -8, 8);
    TH2F *h2_det7 = new TH2F("h2_det7", "Hit Position (u,v) in detID 7", 100, -14, 14, 100, -8, 8);
    TH2F *h2_det9 = new TH2F("h2_det9", "Hit Position (u,v) in detID 9", 100, -14, 14, 100, -8, 8);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    std::cout << "totalNumEvents: " << totalNumEvents << std::endl;
 
    for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        std::fprintf(stdout, "Event: FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
                    eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
 
        const auto& measHits = telEvent->measHits();
        if(!measHits.empty()){
            for(size_t i = 0; i < measHits.size(); i++){
                const auto& hit = measHits[i];
                unsigned int detID = hit->detN();
                float u = hit->u();
                float v = hit->v();
 
                if(detID == 2) {
                    h2_det2->Fill(u, v);
                }
                else if(detID == 3) {
                    h2_det3->Fill(u, v);
                }
                else if(detID == 32) {
                    h2_det32->Fill(u, v);
                }
                else if(detID == 5) {
                    h2_det5->Fill(u, v);
                }
                else if(detID == 7) {
                    h2_det7->Fill(u, v);
                }
                else if(detID == 9) {
                    h2_det9->Fill(u, v);
                }
            }
        }
    }
 
    /*
    std::cout << "detID 2: Collected points " << h2_det2->GetEntries() << std::endl;
    std::cout << "detID 3: Collected points " << h2_det3->GetEntries() << std::endl;
    std::cout << "detID 32: Collected points " << h2_det32->GetEntries() << std::endl;
    std::cout << "detID 5: Collected points " << h2_det5->GetEntries() << std::endl;
    std::cout << "detID 7: Collected points " << h2_det7->GetEntries() << std::endl;
    std::cout << "detID 9: Collected points " << h2_det9->GetEntries() << std::endl;
 */
  
   /* double globalMax = h2_det2->GetMaximum();
    globalMax = std::max(globalMax, h2_det3->GetMaximum());
    globalMax = std::max(globalMax, h2_det32->GetMaximum());
    globalMax = std::max(globalMax, h2_det5->GetMaximum());
    globalMax = std::max(globalMax, h2_det7->GetMaximum());
    globalMax = std::max(globalMax, h2_det9->GetMaximum());
    
    std::cout << "Global maximum value: " << globalMax << std::endl;
    */

    TCanvas *c1 = new TCanvas("c1", "Hit Position (u,v) 2D Histograms 43.raw", 1200, 900);
    c1->Divide(3, 3);
 
    // detID = 2
    c1->cd(1);
    h2_det2->GetXaxis()->SetTitle("u");
    h2_det2->GetYaxis()->SetTitle("v");
    h2_det2->SetStats(1); // 显示统计信息
    gStyle->SetOptStat(1111);
    gStyle->SetStatX(0.9);
   // h2_det2->Scale(1.0 / globalMax);  // 使用全局最大值归一化
   // h2_det2->Scale(1.0 / h2_det2->GetMaximum());按每个探测器上最高柱子归一化
   // h2_det2->Scale(100.0 / h2_det2->GetEntries());
    h2_det2->GetZaxis()->SetTitle("Hit count");
   // h2_det2->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det2->SetMaximum(1.0);
    h2_det2->Draw("LEGO1"); //LEGO1绘制二维直方图作为乐高图，(LEGO2Z 乐高图 + 色条)Z轴表示落在对应 (u, v) 区域的 hit 或数据点的数量(4行画乐高图)
    // h2_det2->Draw("COLZ"); // COLZ绘制热力图，用颜色表示密度，并显示颜色条
   
    // detID = 3
    c1->cd(2);
    h2_det3->GetXaxis()->SetTitle("u");
    h2_det3->GetYaxis()->SetTitle("v");
    h2_det3->SetStats(1);
   // h2_det3->Scale(1.0 / h2_det3->GetEntries());
   // h2_det3->Scale(1.0 / h2_det3->GetMaximum());
    h2_det3->GetZaxis()->SetTitle("Hit count");
   // h2_det3->Scale(100.0 / globalMax);  
   // h2_det3->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det3->SetMaximum(1.0);
    h2_det3->Draw("LEGO1");
    // h2_det3->Draw("COLZ");
    
 
    // detID = 32
    c1->cd(3);
    h2_det32->GetXaxis()->SetTitle("u");
    h2_det32->GetYaxis()->SetTitle("v");
    h2_det32->SetStats(1);
  //  h2_det32->Scale(1.0 / h2_det32->GetMaximum());
   // h2_det32->Scale(1.0 / globalMax);  
   // h2_det32->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det32->Scale(100.0 / h2_det32->GetEntries());
    h2_det32->GetZaxis()->SetTitle("Hit count");
   // h2_det32->SetMaximum(1.0);
    gStyle->SetPalette(kRainBow);
    h2_det32->Draw("LEGO1");
   // h2_det32->Draw("COLZ");
 
    // detID = 5
    c1->cd(4);
    h2_det5->GetXaxis()->SetTitle("u");
    h2_det5->GetYaxis()->SetTitle("v");
    h2_det5->SetStats(1);
   // h2_det5->Scale(1.0 / h2_det5->GetMaximum());
   // h2_det5->Scale(1.0 / globalMax);  
   // h2_det5->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det5->Scale(100.0 / h2_det5->GetEntries());
    h2_det5->GetZaxis()->SetTitle("Hit count ");
   // h2_det5->SetMaximum(1.0);    
    h2_det5->Draw("LEGO1");
   // h2_det5->Draw("COLZ");
 
    // detID = 7
    c1->cd(5);
    h2_det7->GetXaxis()->SetTitle("u");
    h2_det7->GetYaxis()->SetTitle("v");
    h2_det7->SetStats(1);
   // h2_det7->Scale(1.0 / h2_det7->GetMaximum());
    //h2_det7->Scale(1.0 / globalMax);  
   // h2_det7->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det7->Scale(100.0 / h2_det7->GetEntries());
    h2_det7->GetZaxis()->SetTitle("Hit count ");
    //h2_det7->SetMaximum(1.0);    
    h2_det7->Draw("LEGO1");
   // h2_det7->Draw("COLZ");
 
    // detID = 9
    c1->cd(6);
    h2_det9->GetXaxis()->SetTitle("u");
    h2_det9->GetYaxis()->SetTitle("v");
    h2_det9->SetStats(1);
   // h2_det9->Scale(1.0 / h2_det9->GetMaximum());
   // h2_det9->Scale(1.0 / globalMax);  
   // h2_det9->GetZaxis()->SetTitle("Normalized Hits");
   // h2_det9->Scale(100.0 / h2_det9->GetEntries());
    h2_det9->GetZaxis()->SetTitle("Hit count ");
   // h2_det9->SetMaximum(1.0);   
    h2_det9->Draw("LEGO1");
   // h2_det9->Draw("COLZ");
 
    //c1->SaveAs("hit_u_v44.svg");
 
    TFile *outFile = new TFile("hit_u_v.root", "RECREATE");
    h2_det2->Write();
    h2_det3->Write();
    h2_det32->Write();
    h2_det5->Write();
    h2_det7->Write();
    h2_det9->Write();
    c1->Write("canvas_histograms");
    std::cout <<"OVER,The histograms have been saved to the ROOT file" << std::endl;
    
   /* delete h2_det2;
    delete h2_det3;
    delete h2_det32;
    delete h2_det5;
    delete h2_det7;
    delete h2_det9;
    delete c1;
    delete outFile;
    tfile->Close();
    delete tfile;
    */
    return;
}
