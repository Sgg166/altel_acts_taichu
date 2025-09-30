void ITHRchange_HitNum_test () {
    ifstream file1("mean_Raw_size_32.txt");
    ifstream file2("after_mean_Raw_size_32.txt");
    ifstream file3("matchHits_32.txt");
   // ifstream file4("ITHR128.txt");

    std::vector<double> x1, y1, x2, y2, x3, y3;   //, x4, y4;
    double x, y;

    while (file1 >> x >> y) {
        x1.push_back(x);
        y1.push_back(y);
    }
    while (file2 >> x >> y) {
        x2.push_back(x);
        y2.push_back(y);
    }
    while (file3 >> x >> y) {
        x3.push_back(x);
        y3.push_back(y);
    }
  /*  while (file4 >> x >> y) {
        x4.push_back(x);
        y4.push_back(y);
    }
*/
    
/*    double globalMax = 0;
    for (double val : y1) if (val > globalMax) globalMax = val;
    for (double val : y2) if (val > globalMax) globalMax = val;
    for (double val : y3) if (val > globalMax) globalMax = val;
    for (double val : y4) if (val > globalMax) globalMax = val;

    auto normalize_max = [globalMax](std::vector<double>& y) {
      if (globalMax > 0) {
        for (double& val : y) val /= globalMax;
      }
    };

    normalize_max(y1);
    normalize_max(y2);
    normalize_max(y3);
    normalize_max(y4);
全局归一化*/



  /*  auto normalize = [](std::vector<double>& y) {
        double sum = 0;
        for (double val : y) sum += val;
        if (sum > 0) {
            for (double& val : y) val /= sum;
        }
    };
    normalize(y1);
    normalize(y2);
    normalize(y3);
    normalize(y4);
*/
    TGraph* graph1 = new TGraph(x1.size(), &x1[0], &y1[0]);
    TGraph* graph2 = new TGraph(x2.size(), &x2[0], &y2[0]);
    TGraph* graph3 = new TGraph(x3.size(), &x3[0], &y3[0]);
  //  TGraph* graph4 = new TGraph(x4.size(), &x4[0], &y4[0]);

    graph1->SetLineColor(kBlue);
   // graph1->SetMarkerColor(kBlue);
    graph1->SetLineStyle(0);
    graph1->SetLineWidth(2);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(2.0);
   

    graph2->SetLineColor(kRed);         
    graph2->SetLineStyle(1);
    graph2->SetLineWidth(2);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerSize(2.0);

    graph3->SetLineColor(kGreen+2);     
    graph3->SetLineStyle(1);
    graph3->SetLineWidth(2);
    graph3->SetMarkerStyle(22);
    graph3->SetMarkerSize(2.0);
    
   /* graph4->SetLineColor(kMagenta-3);     
    graph4->SetLineStyle(1);
    graph4->SetLineWidth(2);
    graph4->SetMarkerStyle(23);
    graph4->SetMarkerSize(1.0);*/
    graph1->SetTitle("detID=32 Matched_cluster_size  varies with the threshold; X-ITHR;Y-mean");

    TCanvas* c = new TCanvas("c", "Four Graphs on One Canvas", 1000, 700);
    graph1->Draw("APL");
    graph2->Draw("PL SAME");
    graph3->Draw("PL SAME");
  //  graph4->Draw("PL SAME");

    auto legend = new TLegend(0.6, 0.68, 0.9, 0.88);
    legend->AddEntry(graph1, "before mean cluster size", "pl");
    legend->AddEntry(graph2, "after but no screening mean cluster size", "pl");
    legend->AddEntry(graph3, "after and screen mean cluster size", "pl");
   // legend->AddEntry(graph4, "ITHR128", "pl");
    legend->Draw();
    graph1->GetYaxis()->SetRangeUser(0, 2.2);
    graph2->GetYaxis()->SetRangeUser(0, 2.2);
    graph3->GetYaxis()->SetRangeUser(0, 2.2);
   // graph4->GetYaxis()->SetRangeUser(0, 1);
    c->SaveAs("ITHR_hit_distribution.svg");
    c->Update();
}

