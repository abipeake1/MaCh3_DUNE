
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

// Define files and labels globally LLH_PRISM_subsamples_alltimeat0m.root
std::vector<TString> files = {
   "NewLLH/LLH_PRISM_subsamples_0percentOA.root",
   "/NewLLH/LLH_PRISM_subsamples_equaltimeateachposition.root",
   "NewLLH/LLH_PRISM_subsamples_25percentOA.root",
   	"NewLLH/LLH_PRISM_subsamples_50percentOA.root",
    "NewLLH/LLH_PRISM_subsamples_75percentOA.root",
    "NewLLH/LLH_PRISM_subsamples_100percentOA.root"
   
  
};

std::vector<TString> labels = { "100%", "5.26%", "75%", "50%", "25%", "0%"};

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

void PlotLLH() {
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

    //int colors[] = {kRed, kGreen + 2, kAzure + 4, kOrange + 7, kOrange,  kCyan + 1, kAzure + 4, kPink + 5, kMagenta - 9, kViolet + 5};
    int colors[] = {KP8Gray,KP8Azure, KP8Cyan, KP8Pink, KP8Red, KP8Orange, KP8Blue};
    int nColors = sizeof(colors) / sizeof(colors[0]);

    for (int h = 0; h < nHists; h++) {
        TKey* key = (TKey*)list->At(h);
        TString keyname = key->GetName();

        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->SetHeader("Time spent at 0m", "C");

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
