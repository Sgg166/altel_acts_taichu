void Hit_U_V_TH3F(const std::string& rootFilePath){
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
 
    // 创建3D直方图，第三维度为命中点数量 (假设范围0-200)
    TH3F *h3_det2 = new TH3F("h3_det2", "Hit Position (u,v,numHits) in detID 2", 
                            100, -15, 15, 100, -15, 15, 50, 0, 200);
    TH3F *h3_det3 = new TH3F("h3_det3", "Hit Position (u,v,numHits) in detID 3", 
                            100, -15, 15, 100, -15, 15, 50, 0, 200);
    TH3F *h3_det32 = new TH3F("h3_det32", "Hit Position (u,v,numHits) in detID 32", 
                             100, -15, 15, 100, -15, 15, 50, 0, 200);
    TH3F *h3_det5 = new TH3F("h3_det5", "Hit Position (u,v,numHits) in detID 5", 
                            100, -15, 15, 100, -15, 15, 50, 0, 200);
    TH3F *h3_det7 = new TH3F("h3_det7", "Hit Position (u,v,numHits) in detID 7", 
                            100, -15, 15, 100, -15, 15, 50, 0, 200);
    TH3F *h3_det9 = new TH3F("h3_det9", "Hit Position (u,v,numHits) in detID 9", 
                            100, -15, 15, 100, -15, 15, 50, 0, 200);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    std::cout << "totalNumEvents: " << totalNumEvents << std::endl;
 
    for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
        std::fprintf(stdout, "Event: FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
                    eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
 
        const auto& measHits = telEvent->measHits();
        float numHits = telEvent->measHits().size(); // 使用命中点数量作为第三维度
        
        if(!measHits.empty()){
            for(size_t i = 0; i < measHits.size(); i++){
                const auto& hit = measHits[i];
                unsigned int detID = hit->detN();
                float u = hit->u();
                float v = hit->v();
 
                if(detID == 2) {
                    h3_det2->Fill(u, v, numHits);
                }
                else if(detID == 3) {
                    h3_det3->Fill(u, v, numHits);
                }
                else if(detID == 32) {
                    h3_det32->Fill(u, v, numHits);
                }
                else if(detID == 5) {
                    h3_det5->Fill(u, v, numHits);
                }
                else if(detID == 7) {
                    h3_det7->Fill(u, v, numHits);
                }
                else if(detID == 9) {
                    h3_det9->Fill(u, v, numHits);
                }
            }
        }
    }
 
    // 创建3D直方图的2D投影进行可视化
    TCanvas *c1 = new TCanvas("c1", "Hit Position 3D Histograms (numHits as 3rd dim)", 1200, 900);
    c1->Divide(3, 2);
 
    // detID = 2 - XY投影 (u vs v)
    c1->cd(1);
    h3_det2->Project3D("yx")->GetXaxis()->SetTitle("u");
    h3_det2->Project3D("yx")->GetYaxis()->SetTitle("v");
    h3_det2->Project3D("yx")->SetTitle("detID 2: Hit Position (u,v)");
    h3_det2->Project3D("yx")->SetStats(1);
    h3_det2->Project3D("yx")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 2 - XZ投影 (u vs numHits)
    c1->cd(2);
    h3_det2->Project3D("zx")->GetXaxis()->SetTitle("Number of Hits");
    h3_det2->Project3D("zx")->GetYaxis()->SetTitle("u");
    h3_det2->Project3D("zx")->SetTitle("detID 2: u vs NumHits");
    h3_det2->Project3D("zx")->SetStats(1);
    h3_det2->Project3D("zx")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 2 - YZ投影 (v vs numHits)
    c1->cd(3);
    h3_det2->Project3D("zy")->GetXaxis()->SetTitle("Number of Hits");
    h3_det2->Project3D("zy")->GetYaxis()->SetTitle("v");
    h3_det2->Project3D("zy")->SetTitle("detID 2: v vs NumHits");
    h3_det2->Project3D("zy")->SetStats(1);
    h3_det2->Project3D("zy")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 3 - XY投影
    c1->cd(4);
    h3_det3->Project3D("yx")->GetXaxis()->SetTitle("u");
    h3_det3->Project3D("yx")->GetYaxis()->SetTitle("v");
    h3_det3->Project3D("yx")->SetTitle("detID 3: Hit Position (u,v)");
    h3_det3->Project3D("yx")->SetStats(1);
    h3_det3->Project3D("yx")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 32 - XY投影
    c1->cd(5);
    h3_det32->Project3D("yx")->GetXaxis()->SetTitle("u");
    h3_det32->Project3D("yx")->GetYaxis()->SetTitle("v");
    h3_det32->Project3D("yx")->SetTitle("detID 32: Hit Position (u,v)");
    h3_det32->Project3D("yx")->SetStats(1);
    h3_det32->Project3D("yx")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 5 - XY投影
    c1->cd(6);
    h3_det5->Project3D("yx")->GetXaxis()->SetTitle("u");
    h3_det5->Project3D("yx")->GetYaxis()->SetTitle("v");
    h3_det5->Project3D("yx")->SetTitle("detID 5: Hit Position (u,v)");
    h3_det5->Project3D("yx")->SetStats(1);
    h3_det5->Project3D("yx")->Draw("COLZ");
    gPad->SetGrid();
 
    c1->SaveAs("hit_u_v_3d_numHits.svg");
 
    // 另一个画布显示其他探测器的numHits投影
    TCanvas *c2 = new TCanvas("c2", "NumHits Projections for other detectors", 1200, 600);
    c2->Divide(3, 1);
 
    // detID = 7 - XZ投影 (u vs numHits)
    c2->cd(1);
    h3_det7->Project3D("zx")->GetXaxis()->SetTitle("Number of Hits");
    h3_det7->Project3D("zx")->GetYaxis()->SetTitle("u");
    h3_det7->Project3D("zx")->SetTitle("detID 7: u vs NumHits");
    h3_det7->Project3D("zx")->SetStats(1);
    h3_det7->Project3D("zx")->Draw("COLZ");
    gPad->SetGrid();
 
    // detID = 9 - XZ投影 (u vs numHits)
    c2->cd(2);
    h3_det9->Project3D("zx")->GetXaxis()->SetTitle("Number of Hits");
    h3_det9->Project3D("zx")->GetYaxis()->SetTitle("u");
    h3_det9->Project3D("zx")->SetTitle("detID 9: u vs NumHits");
    h3_det9->Project3D("zx")->SetStats(1);
    h3_det9->Project3D("zx")->Draw("COLZ");
    gPad->SetGrid();
 
    // 所有探测器的numHits分布对比
    c2->cd(3);
    h3_det2->Project3D("z")->GetXaxis()->SetTitle("Number of Hits");
    h3_det2->Project3D("z")->SetTitle("detID 2");
    h3_det2->Project3D("z")->SetLineColor(kBlue);
    h3_det2->Project3D("z")->SetStats(0);
    h3_det2->Project3D("z")->DrawNormalized();
    
    h3_det3->Project3D("z")->SetLineColor(kRed);
    h3_det3->Project3D("z")->DrawNormalized("same");
    
    h3_det32->Project3D("z")->SetLineColor(kGreen);
    h3_det32->Project3D("z")->DrawNormalized("same");
    
    h3_det5->Project3D("z")->SetLineColor(kOrange);
    h3_det5->Project3D("z")->DrawNormalized("same");
    
    h3_det7->Project3D("z")->SetLineColor(kMagenta);
    h3_det7->Project3D("z")->DrawNormalized("same");
    
    h3_det9->Project3D("z")->SetLineColor(kCyan);
    h3_det9->Project3D("z")->DrawNormalized("same");
    
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry("z_h3_det2", "detID 2", "l");
    legend->AddEntry("z_h3_det3", "detID 3", "l");
    legend->AddEntry("z_h3_det32", "detID 32", "l");
    legend->AddEntry("z_h3_det5", "detID 5", "l");
    legend->AddEntry("z_h3_det7", "detID 7", "l");
    legend->AddEntry("z_h3_det9", "detID 9", "l");
    legend->Draw();
    
    gPad->SetGrid();
 
    c2->SaveAs("hit_u_v_3d_numHits_distributions.svg");
 
    // 保存3D直方图到ROOT文件
    TFile *outFile = new TFile("hit_u_v_3d_numHits.root", "RECREATE");
    h3_det2->Write();
    h3_det3->Write();
    h3_det32->Write();
    h3_det5->Write();
    h3_det7->Write();
    h3_det9->Write();
    c1->Write("canvas_3d_projections");
    c2->Write("canvas_numHits_distributions");
    std::cout <<"OVER,The 3D histograms with numHits as 3rd dimension have been saved to the ROOT file" << std::endl;
 
    return;
}
