#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TColor.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRint.h>
#include <TStyle.h>
#include "samplePDF/GenericBinningTools.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"

void Write1DHistogramsToFile(std::string OutFileName,
                             std::vector<TH1D *> Histograms) {

  // Now write out the saved hsitograms to file
  auto OutputFile =
      std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));
  OutputFile->cd();
  for (auto Hist : Histograms) {
    Hist->Write();
  }
  OutputFile->Close();

  return;
}
/*
void Write1DHistogramsToPdf(std::string OutFileName,
                            std::vector<TH1D *> Histograms) {

  // Now write out the saved hsitograms to file

  // Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName += ".pdf";

  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();
  c1->Print(std::string(OutFileName + "[").c_str());
  for (auto Hist : Histograms) {
    Hist->Draw("HIST");
    c1->Print(OutFileName.c_str());
  }
  c1->Print(std::string(OutFileName + "]").c_str());

  return;
}*/

void Write1DHistogramsToPdf(std::string OutFileName, std::vector<TH1D *> Histograms) {
  // Remove root from end of file
  OutFileName.erase(OutFileName.find('.'));
  OutFileName += ".pdf";

  // Create a canvas for drawing histograms
  auto c1 = std::unique_ptr<TCanvas>(new TCanvas("c1", "c1", 800, 600));
  c1->cd();

  // Start the multi-page PDF
  c1->Print(std::string(OutFileName + "[").c_str());
  
  // Draw each histogram on a new page
  for (auto Hist : Histograms) {
    c1->Clear();  // Clear the canvas for the next histogram
    Hist->Draw("HIST");  // Draw the histogram
    c1->Print(OutFileName.c_str());  // Print the canvas to the PDF
  }

  // End the multi-page PDF
  c1->Print(std::string(OutFileName + "]").c_str());

  return;
}


int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "Usage: bin/EventRatesDUNEBeam config.cfg" << std::endl;
    return 1;
  }

  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

  covarianceXsec *xsec = nullptr;
  covarianceOsc *osc = nullptr;

  // ####################################################################################
  // Create samplePDFFD objects

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  auto gc1 = std::unique_ptr<TCanvas>(new TCanvas("gc1", "gc1", 800, 600));
  gStyle->SetOptStat(false);
  gc1->Print("GenericBinTest_plus2sigma_3D.pdf[");
  std::vector<TH1D *> DUNEHists;
  for (auto Sample : DUNEPdfs) {
    Sample->reweight();
    Sample->addData(static_cast<TH1D*>(Sample->get1DHist()->Clone((Sample->GetName() + "_asimovdata").c_str())));
    
    xsec->setParameters();

    double nominal = xsec->getNominal(0); //get central value of parameter
    double error = xsec->getDiagonalError(0);
    
    std::cout<<"nominal  = " << nominal << std::endl; 
    std::cout<<"error  = " << error << std::endl; 
      
    Sample->reweight();

    auto myhist_nom = GetGenericBinningTH1(*Sample, "myhistnom", ";global_bin_number;rate", true);

    xsec->setParCurrProp(0, nominal+(2*error));////////// set //+(2*error)
    double current_value = xsec->getParProp(0);
    std::cout<<"current value  = " << current_value << std::endl; 
    Sample->reweight();

    auto myhist_p2 = GetGenericBinningTH1(*Sample, "myhist2", ";global_bin_number;rate", true);
    
    myhist_nom->SetLineColor(kBlack);
    myhist_nom->GetYaxis()->SetRangeUser(0,std::max(myhist_nom->GetMaximum(),myhist_p2->GetMaximum())*1.2);
    myhist_nom->Draw("EHIST");
    myhist_p2->SetLineColor(kRed);
    myhist_p2->SetLineStyle(kDashed);
    myhist_p2->Draw("EHIST SAME");

    gc1->Print("GenericBinTest_plus2sigma_3D.pdf");

       if (Sample->generic_binning.GetNDimensions() == 2) {
        
          //gc1->Divide(3,1);
          //gc1->cd(1);
          xsec->setParCurrProp(0, nominal);////////// set //+(2*error)
          Sample->reweight();
         //auto a = new THStack("a","Stacked 2D histograms");
          auto myhist2_nom = GetGenericBinningTH2(*Sample, "myhist2");
         std::cout<<"hi"<<std::endl;
         //a->Add(myhist2 );
         //a->Add(myhist2_p2);
          myhist2_nom->Draw("COLZ");
          myhist2_nom->SetMinimum(0);
          myhist2_nom->SetMaximum(90e6);
          myhist2_nom->SetTitle("Nominal");
          gc1->Print("GenericBinTest_plus2sigma_3D.pdf");
    
          //gc1->cd(2);
          xsec->setParCurrProp(0, nominal+(2*error));////////// set //+(2*error)
          Sample->reweight();
         //auto a = new THStack("a","Stacked 2D histograms");
          auto myhist2_plus2std = GetGenericBinningTH2(*Sample, "myhist2plus2std");
          myhist2_plus2std->SetMinimum(0);
          myhist2_plus2std->SetMaximum(90e6);  
          //a->Draw("NOSTACK COLZ");
          myhist2_plus2std->SetTitle("+2 #sigma");
          myhist2_plus2std->Draw("COLZ");
          gc1->Print("GenericBinTest_plus2sigma_3D.pdf");
          //gc1->cd(3);
          auto myhist2_ratio = (TH2D*)myhist2_plus2std->Clone("myhist2_ratio");
          myhist2_ratio->Divide(myhist2_nom.get());
          myhist2_ratio->SetMinimum(0.9);
          myhist2_ratio->SetMaximum(1.35);  
          
          //double maxBinContent = myhist2_ratio->GetMaximum();
          //if (maxBinContent != 0) {myhist2_ratio->Scale(1.0 / maxBinContent);}  
          myhist2_ratio->SetTitle("Ratio");
          myhist2_ratio->Draw("COLZ");



          gc1->Print("GenericBinTest_plus2sigma_3D.pdf");
        }

         /*for (auto &slice :
              GetGenericBinningTH1Slices(*Sample, 0, "myslicehist")) {
           slice->Draw();
           gc1->Print("GenericBinTest_plus2sigma_3D.pdf");
         }*/
       
       /*if (Sample->generic_binning.GetNDimensions() == 3) {
         for (auto &slice :
              GetGenericBinningTH2Slices(*Sample, {0, 1}, "myslicehist")) {
           slice->Draw("colz");
           gc1->Print("GenericBinTest_plus2sigma_3D.pdf");
         }*/
       //}

    DUNEHists.push_back(Sample->get1DHist());

    std::string EventRateString =
        fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(),
                  EventRateString);

    std::string LLHString =
        fmt::format("{:.2f}", Sample-> GetLikelihood());
    MACH3LOG_INFO("LLH for {} : {:<5}", Sample->GetName(),
                  LLHString);
       
  }

  gc1->Print("GenericBinTest_plus2sigma_3D_3D.pdf]");

   std::string OutFileName = GetFromManager<std::string>(
       fitMan->raw()["General"]["OutputFile"], "EventRates_Abi.root");

   Write1DHistogramsToFile(OutFileName, DUNEHists);
   Write1DHistogramsToPdf(OutFileName, DUNEHists);
}
