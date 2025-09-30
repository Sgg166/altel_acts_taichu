#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>

void initial_u() {
    std::ifstream inputFile("1.txt", std::ios::in);
    if (!inputFile.is_open()) {
        std::cerr << "Error" << std::endl;
        return;
    }

    TH1F* hist = new TH1F("u_err", "u_err;u_err[mm];Counts", 500, 0.003, 0.01);

    double value;
    while (inputFile >> value) {
        hist->Fill(value);  
    }
    inputFile.close();

    TCanvas* canvas = new TCanvas("canvas", "Histogram", 800, 600);
    hist->SetLineColor(kBlue+1);
    hist->SetFillColor(kAzure-9);
    hist->SetLineWidth(2);

    hist->Draw();

  
    canvas->SaveAs("output_hist.png");

    std::cout << "OVER save  output_hist.png" << std::endl;
}

