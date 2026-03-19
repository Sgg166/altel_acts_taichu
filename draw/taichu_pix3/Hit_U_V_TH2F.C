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
 
    TH2F *h2_det2 = new TH2F("h2_det2",   "Hit Position (u,v) in detID 2 ; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);
    TH2F *h2_det3 = new TH2F("h2_det3",   "Hit Position (u,v) in detID 3 ; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);
    TH2F *h2_det32 = new TH2F("h2_det32", "Hit Position (u,v) in detID 32; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);
    TH2F *h2_det5 = new TH2F("h2_det5",   "Hit Position (u,v) in detID 5 ; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);
    TH2F *h2_det7 = new TH2F("h2_det7",   "Hit Position (u,v) in detID 7 ; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);
    TH2F *h2_det9 = new TH2F("h2_det9",   "Hit Position (u,v) in detID 9 ; u[64bin]; v[32bin] ", 64, -12.8, 12.8, 32, -6.4, 6.4);

    TH1F *h1_det32_u = new TH1F("h1_det32_u", "Hit Position_u in detID 32; u ", 50, -12.8, 12.8);
    TH1F *h1_det32_v = new TH1F("h1_det32_v", "Hit Position_v in detID 32; v ", 50, -6.4, 6.4);


    size_t totalNumEvents = ttreeReader.numEvents();
    std::cout << "totalNumEvents: " << totalNumEvents << std::endl;
 
    for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
       // std::fprintf(stdout, "Event: FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
         //           eventNum, telEvent->eveN(), telEvent->clkN(), telEvent->trajs().size(), telEvent->measHits().size());
 
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
                    h1_det32_u->Fill(u );
                    h1_det32_v->Fill(v );
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
 
    TCanvas *c1 = new TCanvas("c1", "Hit Position (u,v) 2D Histograms 43.raw", 1200, 900);
    c1->Divide(3, 3);
 
    // detID = 2
    c1->cd(1);
    h2_det2->SetStats(1); 
    gStyle->SetOptStat(1111);
    gStyle->SetStatX(0.9);
    h2_det2->GetZaxis()->SetTitle("Hit count");
    //h2_det3->SetMinimum(0);     
    //h2_det3->SetMaximum(500);  
    //h2_det2->Draw("LEGO1"); 
     h2_det2->Draw("COLZ"); 
   
    // detID = 3
    c1->cd(2);
    h2_det3->GetZaxis()->SetTitle("Hit count");
    h2_det3->Draw("COLZ");
    
 
    // detID = 32
    c1->cd(3);
    //h2_det32->SetStats(0);
    h2_det32->GetZaxis()->SetTitle("Hit count");
   // h2_det32->SetMinimum(10);     
   // h2_det32->SetMaximum(20);  
    h2_det32->Draw("COLZ");
   

    c1->cd(7);
    h1_det32_u->GetYaxis()->SetTitle("U_Hit count");
    h1_det32_u->Draw();
 

    c1->cd(8);
    h1_det32_v->GetYaxis()->SetTitle("V_Hit count");
    h1_det32_v->Draw();

    // detID = 5
    c1->cd(4);
    h2_det5->GetZaxis()->SetTitle("Hit count ");
    h2_det5->Draw("COLZ");

    // detID = 7
    c1->cd(5);
    h2_det7->GetZaxis()->SetTitle("Hit count ");
    h2_det7->Draw("COLZ");
 
    // detID = 9
    c1->cd(6);
    h2_det9->GetZaxis()->SetTitle("Hit count ");
    h2_det9->Draw("COLZ");
 
    c1->SaveAs("hit_u_v.svg");
 
    TFile *outFile = new TFile("hit_u_v.root", "RECREATE");
    h2_det2->Write();
    h2_det3->Write();
    h2_det32->Write();
    h1_det32_u->Write();
    h1_det32_v->Write();
    h2_det5->Write();
    h2_det7->Write();
    h2_det9->Write();
    c1->Write("canvas_histograms");
    std::cout <<"OVER,The histograms have been saved to the ROOT file" << std::endl;
    
    return;
}
