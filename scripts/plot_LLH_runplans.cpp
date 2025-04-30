#include "TF1.h"
#include "TH1.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include <iostream>
#include <vector>

#include <fstream>
#include <sstream>
#include <string>
#include "TLatex.h"
//LLH_PRISM_100parameters_5.26percentonaxis
// Define files and labels globally LLH_PRISM_subsamples_alltimeat0m.root
std::vector<TString> files = {
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_0percentonaxis.root",
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_5.26percentonaxis.root",
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_25percentonaxis.root",
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_50percentonaxis.root",
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_75percentonaxis.root",
   "/exp/dune/app/users/abipeake/MaCh3_DUNE_binning/MaCh3_DUNE/Ruplans_LLH/LLH_PRISM_100parameters_100percentonaxis.root"
};

std::vector<TString> labels = { "0%", "5.26 %", "25%", "50%","75%", "100%"};

void SetHistogramStyle(TH1* hist, int color) {
    hist->SetLineColor(color);
    hist->SetLineWidth(12);
    hist->SetMarkerColor(color);
    hist->SetMarkerSize(1.0); // Increase marker size (default is 1.0)
    hist->SetMarkerStyle(20); // Use a filled circle marker
    hist->GetXaxis()->SetRangeUser(0.3, 1.8);
    hist->GetXaxis()->SetTitle("Nominal Value");
    hist->GetYaxis()->SetTitle("Likelihood");
}

void plot_LLH_runplans() {
    if (files.size() != labels.size()) {
        std::cerr << "Error: Number of files and labels must be the same." << std::endl;
        return;
    }

    // Open the output ROOT file for saving histograms
    TFile* outputFile = new TFile("OutputHistograms.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output ROOT file." << std::endl;
        return;
    }

    std::vector<TFile*> rootFiles;
    for (const auto& file : files) {
        TFile* rootFile = new TFile(file);
        if (!rootFile || rootFile->IsZombie()) {
            std::cerr << "Error: Could not open file " << file << std::endl;
            delete rootFile;
            continue;
        }
        rootFiles.push_back(rootFile);
    }

    if (rootFiles.empty()) {
        std::cerr << "No valid files to process." << std::endl;
        outputFile->Close();
        delete outputFile;
        return;
    }

    TList* list = rootFiles[0]->GetListOfKeys();
    if (!list) {
        std::cerr << "Error: No histograms found in the first file." << std::endl;
        outputFile->Close();
        delete outputFile;
        return;
    }

    int nHists = list->GetEntries();
    TCanvas* c0 = new TCanvas("c0", "LLH Scans", 0, 0, 700, 900);
    c0->Print("LLH_runplan3.pdf[");
    //c0->SetLogy();
    //int colors[] = {kRed, kGreen + 2, kAzure + 4, kOrange + 7, kOrange,  kCyan + 1, kAzure + 4, kPink + 5, kMagenta - 9, kViolet + 5};
    //int colors[] = {KP8Gray,KP8Azure, KP8Cyan, KP8Pink, KP8Red, KP8Orange, KP8Blue};
    int colors[] = {kGray+2, kAzure+1, kCyan+1, kMagenta-9, kRed, kOrange+7, kBlue+1};

    int nColors = sizeof(colors) / sizeof(colors[0]);

    for (int h = 0; h < nHists; h++) {
        TKey* key = (TKey*)list->At(h);
        TString keyname = key->GetName();

        if (!keyname.Contains("total_sample")) {
            continue;
        }

        TLegend* legend = new TLegend(0.4, 0.7, 0.75, 0.9);
        legend->SetHeader("Time spent at 0m");
        legend->SetFillStyle(0);  // Transparent background
        legend->SetBorderSize(0);

        bool firstHist = true;

        for (size_t f = 0; f < rootFiles.size(); f++) {
            TH1* ScanHist = dynamic_cast<TH1*>(rootFiles[f]->Get(keyname));
            if (ScanHist) {
                // Set histogram style
                SetHistogramStyle(ScanHist, colors[f % nColors]);

                // Draw the histogram
                if (firstHist) {
                    ScanHist->Draw("HIST P");
                    firstHist = false;
                } else {
                    ScanHist->Draw("HIST P SAME");
                }

                // Add histogram to legend
                legend->AddEntry(ScanHist, labels[f], "l");

                // Save histogram to output ROOT file
                outputFile->cd();
                ScanHist->Write(Form("%s_%s", keyname.Data(), labels[f].Data()));
            } else {
                std::cerr << "Warning: Histogram " << keyname << " not found in file " << files[f] << std::endl;
            }
        }

        // Draw legend and print canvas
        legend->Draw();
        c0->Print("LLH_runplan3.pdf");
        delete legend;
    }

    // === Add parameter name label under the title ===
    std::ifstream paramFile("parameter_list.txt");
    if (!paramFile.is_open()) {
        std::cerr << "Warning: Could not open parameter_list.txt for reading." << std::endl;
    } else {
        std::string line;
        std::getline(paramFile, line); // skip header

        std::string allParams;
        while (std::getline(paramFile, line)) {
            std::istringstream ss(line);
            std::string name;
            double q0_min, q0_max, q3_min, q3_max;
            char comma;

            ss >> name >> comma >> q0_min >> comma >> q0_max >> comma >> q3_min >> comma >> q3_max;

            allParams += name + " ";
        }

        // Draw text just below the plot title
        TLatex* latex = new TLatex(0.5, 0.92, allParams.c_str());
        latex->SetTextSize(0.03);
        latex->SetTextColor(kBlack);
        latex->SetNDC(true);  // Normalized Device Coordinates (so it stays fixed in position)
        latex->SetTextAlign(22); // Center alignment
        latex->Draw();
        gPad->Modified(); 
        gPad->Update();

        paramFile.close();
    }

    // Close the PDF file
    c0->Print("LLH_runplan3.pdf]");

    // Close the output ROOT file
    outputFile->Close();
    delete outputFile;

    // Clean up opened ROOT files
    for (auto file : rootFiles) {
        file->Close();
        delete file;
    }
    delete c0;
}
