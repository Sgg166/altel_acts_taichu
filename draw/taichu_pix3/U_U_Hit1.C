void U_U_Hit1(const std::string& rootFilePath){
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
 
    TH2F *h2_det9_vs_det7_u = new TH2F("h2_det9_vs_det2_x", "Detector 9 X vs Detector 2 X;x_det9;x_det2;Hit count", 1024, -12.8, 12.8, 1024, -12.8, 12.8);
    TH2F *h2_det9_vs_det7_v = new TH2F("h2_det9_vs_det2_y", "Detector 9 Y vs Detector 2 Y;y_det9;y_det2;Hit count", 512, -6.4 , 6.4, 512, -6.4 , 6.4);
 
    size_t totalNumEvents = ttreeReader.numEvents();
    std::cout << "totalNumEvents: " << totalNumEvents << std::endl;
 
    std::vector<unsigned int> requiredDetectors = {9, 7,5,32,3,2};
    size_t validEventCount = 0;
    size_t skippedEventCount = 0;
 
    for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
        std::shared_ptr<altel::TelEvent> telEvent = ttreeReader.createTelEvent(eventNum);
 
        const auto& measHits = telEvent->measHits();
        if(measHits.empty()){
            continue;
        }
        
        
        std::map<unsigned int, bool> detHasHits;
        for(auto detID : requiredDetectors){
            detHasHits[detID] = false;
        }
        
        //Count whether each detector has been hit
        for(size_t i = 0; i < measHits.size(); i++){
            const auto& hit = measHits[i];
            unsigned int detID = hit->detN();
            
            if(detHasHits.count(detID)){
                detHasHits[detID] = true;
            }
        }
        
        //Check whether all the necessary detectors have been hit
        bool allDetectorsHaveHits = true;
        for(auto& [detID, hasHit] : detHasHits){
            if(!hasHit){
                allDetectorsHaveHits = false;
                break;
            }
        }
        
        if(!allDetectorsHaveHits){
            skippedEventCount++;
            continue;
        }
        
        validEventCount++;
        
        std::vector<float> u_det0;
        std::vector<float> v_det0;
        std::vector<float> u_det1;
        std::vector<float> v_det1;
        
        for(size_t i = 0; i < measHits.size(); i++){
            const auto& hit = measHits[i];
            unsigned int detID = hit->detN();
            float u = hit->u();
            float v = hit->v();
            
            if(detID == 9) {
                u_det0.push_back(u);
                v_det0.push_back(v);
            }
            else if(detID == 2) {
                u_det1.push_back(u);
                v_det1.push_back(v);
            }
        }
        
        if(!u_det0.empty() && !u_det1.empty()){
            for(float u0 : u_det0){
                for(float u1 : u_det1){
                    h2_det9_vs_det7_u->Fill(u0, u1);
                }
            }
        }
        
        if(!v_det0.empty() && !v_det1.empty()){
            for(float v0 : v_det0){
                for(float v1 : v_det1){
                    h2_det9_vs_det7_v->Fill(v0, v1);
                }
            }
        }
    }
    
 //   std::cout << "Total events processed: " << std::min((size_t)100000, totalNumEvents) << std::endl;
    std::cout << "Total events processed: " <<  totalNumEvents << std::endl;
    std::cout << "Valid events (all detectors have hits): " << validEventCount << std::endl;
    std::cout << "Skipped events (missing detector hits): " << skippedEventCount << std::endl;
    std::cout << "Efficiency: " << (double)validEventCount / (double) totalNumEvents   << std::endl;
 
    TCanvas *c1 = new TCanvas("c1", "Hit Position (u,v) 2D Histograms 43.raw", 1000, 1000);
    c1->Divide(1, 2);
 
    c1->cd(1);
    h2_det9_vs_det7_u->SetStats(1);
    gStyle->SetOptStat(1111);
    gStyle->SetStatX(0.9);
    h2_det9_vs_det7_u->Draw("COLZ");
 
    c1->cd(2);
    h2_det9_vs_det7_v->Draw("COLZ");
 
    TFile *outFile = new TFile("/home/pub/B/altel_acts/Draw/root_file/hit_u_u.root", "RECREATE");
    h2_det9_vs_det7_u->Write();
    h2_det9_vs_det7_v->Write();
    c1->Write();
    std::cout <<"OVER,The histograms have been saved to the ROOT file" << std::endl;
 
    return;
}
