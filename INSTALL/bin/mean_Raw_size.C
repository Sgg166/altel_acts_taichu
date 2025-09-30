#include <iostream>
#include <fstream>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
void mean_Raw_size() {
    std::ifstream inputFile("Efficiency.txt", std::ios::in);
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open the file 32IDmeanRaw.txt" << std::endl;
        return;
    }
    TGraph* graph = new TGraph();
    int pointIndex = 0;
    double x, y;
    while (inputFile >> x >> y) {
        graph->SetPoint(pointIndex, x, y);
        pointIndex++;
    }
    inputFile.close();
    TCanvas* canvas = new TCanvas("canvas", "Graph from Data", 800, 600);
    graph->SetTitle("detID=32 numOriginHit and Efficiency;X-numOriginHit;Y-Efficiency");
    graph->SetMarkerStyle(20);  //点的形状
    graph->SetMarkerColor(4);   // 蓝色点
    graph->SetLineColor(2);     // 红色连线
    graph->SetMarkerSize(1.5);  // 点大小
    graph->Draw("APL");         // 绘制点、连线和坐标轴
    graph->GetXaxis()->SetRangeUser(0, 256);  // 设置 X 轴范围
    graph->GetXaxis()->SetTickLength(0.02);
   // graph->GetXaxis()->SetNdivisions(8); // 8 段 = 每段 32
    graph->GetYaxis()->SetRangeUser(0, 1);  // 设置 Y 轴范围
    graph->GetYaxis()->SetTickLength(0.02);    
    graph->GetXaxis()->SetTitleOffset(1.2);   // 设置 X 轴标题偏移
    graph->GetYaxis()->SetTitleOffset(1.2);   // 设置 Y 轴标题偏移

 //   canvas->SaveAs("Matched_claster_size32.svg");
 
    std::cout << "OVER， output_graph.png" << std::endl;

    // 释放内存
    //delete graph;
   // delete canvas;
}
