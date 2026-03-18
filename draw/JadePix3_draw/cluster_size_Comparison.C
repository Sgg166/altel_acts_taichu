#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TPad.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

void cluster_size_Comparison() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);

    TCanvas *c = new TCanvas("c", "Overlay cluster-size plots", 1100, 800);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.06);
    c->SetTicks(1, 1);

    TH1F *h1 = new TH1F("h1", "Cluster size distribution of the DUT. ", 15, 0, 15);

    double y1[15] = {
        0.000,
        0.177,
        0.307,
        0.152,
        0.117,
        0.068,
        0.059,
        0.026,
        0.020,
        0.011,
        0.009,
        0.007,
        0.006,
        0.005,
        0.004
    
    };

    for (int i = 0; i < 15; ++i) {
        h1->SetBinContent(i + 1, y1[i]);
    }

    TH1F *h2 = new TH1F("h2", "", 15, 0, 15);

    double y2[15] = {
        0.000,
        0.190,
        0.319,
        0.145,
        0.110,
        0.067,
        0.055,
        0.025,
        0.020,
        0.010,
        0.009,
        0.008,
        0.007,
        0.006,
        0.005
    };

    for (int i = 0; i < 15; ++i) {
        h2->SetBinContent(i + 1, y2[i]);
    }

    double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(ymax * 1.20);
    h1->SetMinimum(0.0);

    h1->GetXaxis()->SetTitle("Cluster size ");
    h1->GetYaxis()->SetTitle("Entries (normalized)");
    h1->GetXaxis()->SetTitleSize(0.055);
    h1->GetYaxis()->SetTitleSize(0.055);
    h1->GetXaxis()->SetLabelSize(0.045);
    h1->GetYaxis()->SetLabelSize(0.045);
    h1->GetXaxis()->SetTitleOffset(1.15);
    h1->GetYaxis()->SetTitleOffset(0.95);

    h1->SetLineColor(kGreen + 2);
    h1->SetLineWidth(2);
    h1->SetFillColor(kGreen - 9);
    h1->SetFillStyle(3354);

    h2->SetLineColor(kViolet + 1);
    h2->SetLineWidth(2);
    h2->SetFillColor(kViolet - 7);
    h2->SetFillStyle(3345);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.62, 0.73, 0.93, 0.90);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(h1, "Run316 Data", "f");
    leg->AddEntry(h2, "Paper Data" , "f");
    leg->Draw();

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(42);
    txt.SetTextSize(0.045);
    txt.DrawLatex(0.48, 0.60, "Run316 Data mean cluster size: 3.7");
    txt.DrawLatex(0.48, 0.52, "Paper Data  mean cluster size: 3.7");

    c->Update();
  //  c->SaveAs("two_cluster_plots_overlay.png");
  //  c->SaveAs("two_cluster_plots_overlay.pdf");
}
