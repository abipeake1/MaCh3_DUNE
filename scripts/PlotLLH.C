#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"


#include <iostream>

void PlotLLH (TString inputfile)
{
    // Open the ROOT file
    TFile* file = new TFile(inputfile);
    
    // Check if the file is successfully opened
    if (!file || file->IsZombie()) {
        std::cout << "Error: Could not open file " << inputfile << std::endl;
        return;
    }

    // Get the list of keys from the file
    TList* list = file->GetListOfKeys();
    
    // Create a canvas
    TCanvas* c0 = new TCanvas("c0", "c0", 0, 0, 700, 900);
    c0->Print("LLHScans.pdf[");

    // Loop over all keys in the file
    for (int i = 0; i < list->GetEntries(); i++) {
        TKey* key = (TKey*)list->At(i); // Cast TObject to TKey
        TString keyname = key->GetName(); // Get the name of the key
        
        // Retrieve the object as a TH1D histogram
        TH1D* ScanHist = (TH1D*)file->Get(keyname);

        // Check if the object is a valid TH1D histogram
        if (ScanHist && ScanHist->InheritsFrom("TH1")) {
            ScanHist->Draw("HIST P");  // Draw histogram without error bars
            c0->Print("LLHScans.pdf"); // Print the current canvas to the PDF
        } else {
            std::cout << "Warning: Object " << keyname << " is not a TH1D" << std::endl;
        }
    }

    // Close the PDF
    c0->Print("LLHScans.pdf]");
    
    // Cleanup
    file->Close();
    delete c0;
    delete file;
}