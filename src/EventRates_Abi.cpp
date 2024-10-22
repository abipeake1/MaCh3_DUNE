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
  gc1->Print("GenericBinTest_plus2sigma.pdf[");

  std::vector<TH1D *> DUNEHists;
  for (auto Sample : DUNEPdfs) {
    Sample->reweight();
    xsec->setParameters();
    double nominal =xsec->getNominal(0); //get central value of parameter
    double error = xsec->getDiagonalError(0);
    std::cout<<"nominal  = " << nominal << std::endl; 
    std::cout<<"error  = " << error << std::endl; 
    xsec->setParCurrProp(0, nominal+(2*error));////////// set //+(2*error)
    std::cout<< "nominal value = " << nominal <<std::endl;;
    double current_value = xsec->getParProp(0);
    std::cout<<"current value  = " << current_value << std::endl; 
    
    Sample->reweight();
    DUNEHists.push_back(Sample->get1DHist());

    if (Sample->generic_binning.GetNDimensions()) {

      auto myhist = GetGenericBinningTH1(*Sample, "myhist");
      myhist->Scale(1, "WIDTH");
      myhist->Draw();
      gc1->Print("GenericBinTest_plus2sigma.pdf");

      if (Sample->generic_binning.GetNDimensions() == 2) {
        auto myhist2 = GetGenericBinningTH2(*Sample, "myhist2");
        myhist2->Draw("COLZ TEXT");
        gc1->Print("GenericBinTest_plus2sigma.pdf");

        for (auto &slice :
             GetGenericBinningTH1Slices(*Sample, 0, "myslicehist")) {
          slice->Draw();
          gc1->Print("GenericBinTest_plus2sigma.pdf");
        }
      }
      if (Sample->generic_binning.GetNDimensions() == 3) {
        for (auto &slice :
             GetGenericBinningTH2Slices(*Sample, {0, 1}, "myslicehist")) {
          slice->Draw();
          gc1->Print("GenericBinTest_plus2sigma.pdf");
        }
      }
    }

    DUNEHists.push_back(Sample->get1DHist());

    std::string EventRateString =
        fmt::format("{:.2f}", Sample->get1DHist()->Integral());
    MACH3LOG_INFO("Event rate for {} : {:<5}", Sample->GetName(),
                  EventRateString);
  }

  gc1->Print("GenericBinTest_plus2sigma.pdf]");


  std::string OutFileName = GetFromManager<std::string>(
      fitMan->raw()["General"]["OutputFile"], "EventRates_Abi.root");

  Write1DHistogramsToFile(OutFileName, DUNEHists);
  Write1DHistogramsToPdf(OutFileName, DUNEHists);
}
