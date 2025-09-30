void draw_two_graphs() {
    ifstream file1("mean_Raw_size_32.txt");
    ifstream file2("after_mean_Raw_size_32.txt");

    std::vector<double> x1, y1, x2, y2;
    double x, y;

    while (file1 >> x >> y) {
        x1.push_back(x);
        y1.push_back(y);
    }

    while (file2 >> x >> y) {
        x2.push_back(x);
        y2.push_back(y);
    }

    TGraph* graph1 = new TGraph(x1.size(), &x1[0], &y1[0]);
    TGraph* graph2 = new TGraph(x2.size(), &x2[0], &y2[0]);

    graph1->SetLineColor(kBlue);
    graph1->SetLineWidth(2);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(1.0);
    graph1->SetTitle("detID=32 measHits.size() varies with the threshould; X-measHits.size();Y-Probability");

    graph2->SetLineColor(kRed);
    graph2->SetLineWidth(2);
    graph2->SetLineStyle(1);
    graph2->SetMarkerStyle(21);
    graph2->SetMarkerSize(1.0);
    TCanvas* c = new TCanvas("c", "Two Graphs on One Canvas", 1000, 700);
    graph1->Draw("APL");
    graph2->Draw("PL SAME");

    // 添加图例
    auto legend = new TLegend(0.6, 0.7, 0.9, 0.85);
    legend->AddEntry(graph1, "mean_Raw_size_32", "l");
    legend->AddEntry(graph2, "after_mean_Raw_size_32", "l");
    legend->Draw();
   // graph1->GetYaxis()->SetRangeUser(0, 1);
   // graph2->GetYaxis()->SetRangeUser(0, 1);
    c->SaveAs("detID=32 before_and_after.svg");

    c->Update();
}

