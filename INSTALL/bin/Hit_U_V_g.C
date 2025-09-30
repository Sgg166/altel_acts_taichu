void Hit_U_V_g(const std::string& rootFilePath){
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
  
  // 创建存储不同detID的u和v值的向量
  std::vector<float> u_det2;
  std::vector<float> v_det2;
  std::vector<float> u_det3;
  std::vector<float> v_det3;
  std::vector<float> u_det32;
  std::vector<float> v_det32;
  std::vector<float> u_det5;
  std::vector<float> v_det5;
  std::vector<float> u_det7;
  std::vector<float> v_det7;
  std::vector<float> u_det9;
  std::vector<float> v_det9;
 
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
        
        // 根据detID存储坐标
        if(detID == 2) {
          u_det2.push_back(u);
          v_det2.push_back(v);
        }
        else if(detID == 3) {
          u_det3.push_back(u);
          v_det3.push_back(v);
        }
        else if(detID == 32) {
          u_det32.push_back(u);
          v_det32.push_back(v);
        }
        else if(detID == 5) {
          u_det5.push_back(u);
          v_det5.push_back(v);
        }
        else if(detID == 7) {
          u_det7.push_back(u);
          v_det7.push_back(v);
        }
        else if(detID == 9) {
          u_det9.push_back(u);
          v_det9.push_back(v);
        }
      }
    }
  }
  /*
  std::cout << "detID 2: Collected points " << u_det2.size()  << std::endl;
  std::cout << "detID 2: Collected points " << v_det2.size()  << std::endl;
  std::cout << "detID 3: Collected points " << u_det3.size()  << std::endl;
  std::cout << "detID 3: Collected points " << v_det3.size()  << std::endl;
  std::cout << "detID 32: Collected points " << u_det32.size()  << std::endl;
  std::cout << "detID 32: Collected points " << v_det32.size()  << std::endl;
  std::cout << "detID 5: Collected points " << u_det5.size()  << std::endl;
  std::cout << "detID 5: Collected points " << v_det5.size()  << std::endl;
  std::cout << "detID 7: Collected points " << u_det7.size()  << std::endl;
  std::cout << "detID 7: Collected points " << v_det7.size()  << std::endl;
  std::cout << "detID 9: Collected points " << u_det9.size()  << std::endl;
  std::cout << "detID 9: Collected points " << v_det9.size()  << std::endl;
 */
  TCanvas *c1 = new TCanvas("c1", "Hit Position (u,v) Scatter Plot", 1200, 900);
  c1->Divide(3, 3);
  
  // detID = 2
  c1->cd(1);
  TGraph *gr2 = new TGraph(u_det2.size(), &u_det2[0], &v_det2[0]);
  gr2->SetTitle("Hit Position (u,v) in detID 2");
  gr2->GetXaxis()->SetTitle("u");
  gr2->GetYaxis()->SetTitle("v");
  gr2->SetMarkerStyle(20);  // 实心圆点
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerSize(0.5);  // 稍小的点
  gr2->Draw("AP");  // A表示坐标轴，P表示绘制点
  
  // detID = 3
  c1->cd(2);
  TGraph *gr3 = new TGraph(u_det3.size(), &u_det3[0], &v_det3[0]);
  gr3->SetTitle("Hit Position (u,v) in detID 3");
  gr3->GetXaxis()->SetTitle("u");
  gr3->GetYaxis()->SetTitle("v");
  gr3->SetMarkerStyle(20);
  gr3->SetMarkerColor(kBlue);
  gr3->SetMarkerSize(0.5);
  gr3->Draw("AP");
  
  //detID = 32
  c1->cd(3);
  TGraph *gr32 = new TGraph(u_det32.size(), &u_det32[0], &v_det32[0]);
  gr32->SetTitle("Hit Position (u,v) in detID 32");
  gr32->GetXaxis()->SetTitle("u");
  gr32->GetYaxis()->SetTitle("v");
  gr32->SetMarkerStyle(20);
  gr32->SetMarkerColor(kBlack);
  gr32->SetMarkerSize(0.5);
  gr32->Draw("AP");


  // detID = 5
  c1->cd(4);
  TGraph *gr5 = new TGraph(u_det5.size(), &u_det5[0], &v_det5[0]);
  gr5->SetTitle("Hit Position (u,v) in detID 5");
  gr5->GetXaxis()->SetTitle("u");
  gr5->GetYaxis()->SetTitle("v");
  gr5->SetMarkerStyle(20);
  gr5->SetMarkerColor(kGreen+2);
  gr5->SetMarkerSize(0.5);
  gr5->Draw("AP");
  
  // detID = 7
  c1->cd(5);
  TGraph *gr7 = new TGraph(u_det7.size(), &u_det7[0], &v_det7[0]);
  gr7->SetTitle("Hit Position (u,v) in detID 7");
  gr7->GetXaxis()->SetTitle("u");
  gr7->GetYaxis()->SetTitle("v");
  gr7->SetMarkerStyle(20);
  gr7->SetMarkerColor(kMagenta);
  gr7->SetMarkerSize(0.5);
  gr7->Draw("AP");
  
  // detID = 9
  c1->cd(6);
  TGraph *gr9 = new TGraph(u_det9.size(), &u_det9[0], &v_det9[0]);
  gr9->SetTitle("Hit Position (u,v) in detID 9");
  gr9->GetXaxis()->SetTitle("u");
  gr9->GetYaxis()->SetTitle("v");
  gr9->SetMarkerStyle(20);
  gr9->SetMarkerColor(kOrange+1);
  gr9->SetMarkerSize(0.5);
  gr9->Draw("AP");
  
  c1->SaveAs("hit_position_scatter.svg");
  
  TFile *outFile = new TFile("hit_position_scatter.root", "RECREATE");
  gr2->Write("gr_det2");
  gr3->Write("gr_det3");
  gr32->Write("gr_det32");
  gr5->Write("gr_det5");
  gr7->Write("gr_det7");
  gr9->Write("gr_det9");
  c1->Write("canvas_individual");
  std::cout <<"OVER,The data has been saved to the ROOT file" << std::endl;
  return;
}
