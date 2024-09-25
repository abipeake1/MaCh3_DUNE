#include <TROOT.h>
#include "TMath.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH2.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>


void split_write3d(){

    auto f1= new TFile("/home/abipeake/Mach3/MaCh3_DUNE/Results/offaxissamples_3D_oaposition_250924.root"); //FHC_numuosc_3D
    TH3D *h = (TH3D*)f1->Get("Abis3DHistogram");
    //TH3D *osc_3Dhis = (TH3D*)f3->Get("FHC_numu_osc_3D");

    //std::vector<TH3D*> list_of_histograms = 
    // Open a ROOT file and save the formula, function and histogram
    TFile myfile("/home/abipeake/Mach3/MaCh3_DUNE/Unfoldedhists/offaxissamples_3D_oaposition_250924.root","RECREATE");

    h->Write();
    h->SetDirectory(nullptr);
    TH1::SetDefaultSumw2(true);
    
    int number_of_ELep_bins = h->GetXaxis()->GetNbins();
    int number_of_ENuReco_bins = h->GetYaxis()->GetNbins();
    int number_of_oabins = h->GetZaxis()->GetNbins();

    auto ptpzproj = std::unique_ptr<TH2>(static_cast<TH2*>(h->Project3D("yx")));
    std::string hist_2Dname = ptpzproj->GetName();

    std::cout << "The 3D histogram " << h->GetName() << "has" << number_of_ELep_bins << "bins in pt (x axis) and " <<  number_of_ENuReco_bins << "number_of_ detectorposition  bins" <<"and" << h->GetZaxis()->GetNbins()<< "number_of_ENu_Reco"<<std::endl;
    for (int x = 0; x < h->GetXaxis()->GetNbins(); ++x) {
        for (int y = 0; y < h->GetYaxis()->GetNbins(); ++y) {
            TH1 *proj =
                h->ProjectionZ((std::string(h->GetName()) + "_x" +
                                std::to_string(x) + "_y" + std::to_string(y))
                                .c_str(),
                                x+1, x+1 , y+1 , y+1 , "e");
            std::stringstream ss;
            ss << h->GetXaxis()->GetBinLowEdge(x + 1) << " < ELep < "
            << h->GetXaxis()->GetBinUpEdge(x + 1) << " [GeV/c], "
            << h->GetYaxis()->GetBinLowEdge(y + 1) << " < theta_lep < "
            << h->GetYaxis()->GetBinUpEdge(y + 1) << " [Ge/cc]";

            double proj_cell_area = (h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1)) *
                                (h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1));
            double x_bin_width_pt = h->GetXaxis()->GetBinUpEdge(x + 1) - h->GetXaxis()->GetBinLowEdge(x + 1);
            double y_bin_width_pz = h->GetYaxis()->GetBinUpEdge(y + 1) - h->GetYaxis()->GetBinLowEdge(y + 1);
            
            proj->Scale(1.0/proj_cell_area,"WIDTH"); 
            proj->SetTitle(ss.str().c_str());
            proj->GetXaxis()->SetTitle("Reconstricted Neutrino Energy (GeV/c) ");
            std::string hist_1Dname = proj->GetName();
            proj->Write();
            proj->SetDirectory(nullptr);
            delete proj;
        }
    
    }
    ptpzproj->Scale(1.0,"WIDTH");  //get  root to divide out the bin width
    ptpzproj->Write();
    //make sure to tell ROOT that it doesn't own this histogram so that we can delete it
    ptpzproj->SetDirectory(nullptr);



}


/*
void Write(std::string drawOpt) {    
  split_write3d(Abis3DHistogram_forplotting_unosc);
  split_write3d(Abis3DHistogram_forplotting_osc);
}*/