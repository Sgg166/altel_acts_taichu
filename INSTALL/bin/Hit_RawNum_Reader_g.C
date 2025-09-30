void  Hit_RawNum_Reader_g(const std::string& rootFilePath){
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
 
  // 创建散点图数据容器
  std::vector<double> eventNumbers_det2, rawNums_det2;
  std::vector<double> eventNumbers_det3, rawNums_det3;
  std::vector<double> eventNumbers_det32, rawNums_det32;
  std::vector<double> eventNumbers_det5, rawNums_det5;
  std::vector<double> eventNumbers_det7, rawNums_det7;
  std::vector<double> eventNumbers_det9, rawNums_det9;
  
  std::vector<double> eventNumbers_hitCount, hitCounts_det32;
 
  size_t totalNumEvents  = ttreeReader.numEvents();
  size_t total_det32_hits = 0;
  for(size_t eventNum = 0; eventNum<totalNumEvents; eventNum++){
    std::shared_ptr<altel::TelEvent> telEvent  = ttreeReader.createTelEvent(eventNum);
    std::fprintf(stdout, "Event:  FileEvent #%zu, event #%u, clock/trigger #%lu, numTraj %zu, numMeasHit %zu \n",
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
          eventNumbers_det2.push_back(eventNum);
          rawNums_det2.push_back(numRaws);
        }
        else if(detID == 3) {
          eventNumbers_det3.push_back(eventNum);
          rawNums_det3.push_back(numRaws);
        }
        else if(detID == 32)
        {
          eventNumbers_det32.push_back(eventNum);
          rawNums_det32.push_back(numRaws);
          det32_event_hit_count++;
          total_det32_hits++;
        }
        else if(detID == 5) {
          eventNumbers_det5.push_back(eventNum);
          rawNums_det5.push_back(numRaws);
        }
        else if(detID == 7) {
          eventNumbers_det7.push_back(eventNum);
          rawNums_det7.push_back(numRaws);
        }
        else if(detID == 9) {
          eventNumbers_det9.push_back(eventNum);
          rawNums_det9.push_back(numRaws);
        }
      }
    }
    eventNumbers_hitCount.push_back(eventNum);
    hitCounts_det32.push_back(det32_event_hit_count);
  }
 
  TCanvas *c1 = new TCanvas("c1", "Raw Pixels per Hit Scatter Plot", 1200, 900);
  c1->Divide(3, 3);
 
  // 创建散点图
  TGraph *g_det2 = new TGraph(eventNumbers_det2.size(), &eventNumbers_det2[0], &rawNums_det2[0]);
  TGraph *g_det3 = new TGraph(eventNumbers_det3.size(), &eventNumbers_det3[0], &rawNums_det3[0]);
  TGraph *g_det32 = new TGraph(eventNumbers_det32.size(), &eventNumbers_det32[0], &rawNums_det32[0]);
  TGraph *g_det5 = new TGraph(eventNumbers_det5.size(), &eventNumbers_det5[0], &rawNums_det5[0]);
  TGraph *g_det7 = new TGraph(eventNumbers_det7.size(), &eventNumbers_det7[0], &rawNums_det7[0]);
  TGraph *g_det9 = new TGraph(eventNumbers_det9.size(), &eventNumbers_det9[0], &rawNums_det9[0]);
  
  TGraph *g_hitCount = new TGraph(eventNumbers_hitCount.size(), &eventNumbers_hitCount[0], &hitCounts_det32[0]);

  g_det2->SetTitle("Raw Pixels per Hit vs Event Number - detID 2");
  g_det2->SetMarkerColor(kRed);
  g_det2->SetMarkerStyle(20); // 圆形标记
  g_det2->SetMarkerSize(0.5);
  g_det2->GetXaxis()->SetTitle("Event Number");
  g_det2->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_det3->SetTitle("Raw Pixels per Hit vs Event Number - detID 3");
  g_det3->SetMarkerColor(kBlue);
  g_det3->SetMarkerStyle(20);
  g_det3->SetMarkerSize(0.5);
  g_det3->GetXaxis()->SetTitle("Event Number");
  g_det3->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_det32->SetTitle("Raw Pixels per Hit vs Event Number - detID 32");
  g_det32->SetMarkerColor(kViolet);
  g_det32->SetMarkerStyle(20);
  g_det32->SetMarkerSize(0.5);
  g_det32->GetXaxis()->SetTitle("Event Number");
  g_det32->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_det5->SetTitle("Raw Pixels per Hit vs Event Number - detID 5");
  g_det5->SetMarkerColor(kGreen);
  g_det5->SetMarkerStyle(20);
  g_det5->SetMarkerSize(0.5);
  g_det5->GetXaxis()->SetTitle("Event Number");
  g_det5->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_det7->SetTitle("Raw Pixels per Hit vs Event Number - detID 7");
  g_det7->SetMarkerColor(kMagenta);
  g_det7->SetMarkerStyle(20);
  g_det7->SetMarkerSize(0.5);
  g_det7->GetXaxis()->SetTitle("Event Number");
  g_det7->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_det9->SetTitle("Raw Pixels per Hit vs Event Number - detID 9");
  g_det9->SetMarkerColor(kOrange);
  g_det9->SetMarkerStyle(20);
  g_det9->SetMarkerSize(0.5);
  g_det9->GetXaxis()->SetTitle("Event Number");
  g_det9->GetYaxis()->SetTitle("Number of Raw Pixels per Hit");

  g_hitCount->SetTitle("Hit Count per Event vs Event Number - detID 32");
  g_hitCount->SetMarkerColor(kBlack);
  g_hitCount->SetMarkerStyle(21); // 方形标记
  g_hitCount->SetMarkerSize(0.5);
  g_hitCount->GetXaxis()->SetTitle("Event Number");
  g_hitCount->GetYaxis()->SetTitle("Number of Hits per Event");

  c1->cd(1);
  g_det2->Draw("AP");

  c1->cd(2);
  g_det3->Draw("AP");

  c1->cd(3);
  g_det32->Draw("AP");

  c1->cd(4);
  g_hitCount->Draw("AP");

  c1->cd(5);
  g_det5->Draw("AP");

  c1->cd(6);
  g_det7->Draw("AP");

//  c1->cd(7);
 // g_det9->Draw("AP");

  c1->SaveAs("Hit_RawNum_scatter.svg");

  TFile *outFile = new TFile("hits_analysis_scatter.root", "RECREATE");
  g_det2->Write();
  g_det3->Write();
  g_det32->Write();
  g_det5->Write();
  g_det7->Write();
  g_det9->Write();
  g_hitCount->Write();
  c1->Write();

  std::cout << "Total hits for detID 32: " << total_det32_hits << std::endl;
  std::cout << "Average hits per event for detID 32: " << (double)total_det32_hits/totalNumEvents << std::endl;
  std::cout << "OVER, The scatter plots have been saved to hits_analysis_scatter.root and Hit_RawNum_scatter.svg" << std::endl;

  // 清理内存
 // delete c1;
 // delete outFile;
 // delete tfile;

  return;
}
