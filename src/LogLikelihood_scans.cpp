#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDF/GenericBinningTools.h"
#include "mcmc/mcmc.h"

//////////////// Copied over from the old repo, with some changes to how the SamplePDFObjects are read in
std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}                                                            
                                                             
void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  name = getNameNoExt(name, ".root") ;                       
  name = name + legend + ".root" ;                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".root") ;                       
  name = name + ".png" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".png") ;                        
  name = name + ".pdf" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
} 



int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }


  // Read config
  //manager *fitMan = new manager(argv[1]);

  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

  covarianceXsec *xsec = nullptr;
  covarianceOsc *osc = nullptr;

  // ####################################################################################
  // Create samplePDFFD objects

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  auto gc1 = std::unique_ptr<TCanvas>(new TCanvas("gc1", "gc1", 800, 600));
  gStyle->SetOptStat(false);
  gc1->Print("GenericBinTest_someoa.pdf[");
  

    // Setting flat priors based on XSECPARAMFLAT list in configuration file 

  std::vector<double> xsecpar = xsec->getNominalArray();  
  // Set to Asimov values
  // If wanting to do an Asimov fit to non-nominal values of xsec params set them here!
  // e.g. xsecpar[i] = ....

  xsec->setParameters();


  // Print Asimov event rates to check
  std::cout << "-------- Event rates for Asimov Data  ------------" << std::endl;
  
  std::vector<std::string> sample_names;
  std::vector<TH1D *> DUNEHists;
  for (auto Sample : DUNEPdfs) {
    Sample->reweight();
	std::string name = Sample -> GetName();
    TString NameTString = TString(name.c_str());
	sample_names.push_back(name);

	if (Sample -> GetNDim() == 1) {
      TH1D *Asimov_1D = (TH1D*)Sample->get1DHist()->Clone(NameTString+"_asimov");
      std::cout << name.c_str() << ": " << Asimov_1D->Integral() << std::endl;
      Sample -> addData(Asimov_1D); 
	}
      
	if (Sample -> GetNDim() == 2) {
      TH2D *Asimov_2D = (TH2D*)Sample->get2DHist()->Clone(NameTString+"_asimov");
      std::cout << name.c_str() << ": " << Asimov_2D->Integral() << std::endl;
      Sample -> addData(Asimov_2D); 
	}

  }

  std::string OutFileName = fitMan->raw()["General"]["OutputFile"].as<std::string>();

  TCanvas *sys_canv = new TCanvas("sys_canv","",1200,600);
  auto OutFile = std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));
  // Save pars before modifying
  std::vector<double> tpar = xsecpar;

  unsigned int n_points = fitMan->raw()["LLHScans"]["ParPoints"].as<int>();
  int n_sigma = 3;

  // Loop over parameters
  for (unsigned  par=0; par < xsecpar.size(); par++) {

    std::vector<std::string> histnames;
    std::vector<std::string> histtitles;
    std::vector<TH1D*> llh_hists;

    xsecpar = xsec->getNominalArray();
    double lowerb = xsec->getNominal(par) - (n_sigma+0.01) * sqrt((*xsec->getCovMatrix())(par,par));
    double upperb = xsec->getNominal(par) + (n_sigma+0.01) * sqrt((*xsec->getCovMatrix())(par,par));
 
	 // Create Histogram for each sample + systematic penalty + total sample	
    for(unsigned sample_i=0; sample_i < DUNEPdfs.size() ; ++sample_i) {
	  std::string histname = "xsec_" + std::to_string(par) + "_llh_" + sample_names[sample_i];
	  std::string histtitle = "xsec_" + std::to_string(par) + "_" + sample_names[sample_i];
      TH1D *hScan = new TH1D(histname.c_str(), histtitle.c_str(), n_points, lowerb, upperb);
	  llh_hists.push_back(hScan);
	}
	
    std::string histname_systpen = "xsec_" + std::to_string(par) + "_llh_syst";
    std::string histtitle_systpen = "xsec_" + std::to_string(par) + "_syst";
    TH1D *hScan_pen = new TH1D(histname_systpen.c_str(), histtitle_systpen.c_str(), n_points, lowerb, upperb);
	llh_hists.push_back(hScan_pen);

    std::string histname_total = "xsec_" + std::to_string(par) + "_llh_total_sample";
    std::string histtitle_total = "xsec_" + std::to_string(par) + "_total_sample";
    TH1D *hScan_total = new TH1D(histname_total.c_str(), histtitle_total.c_str(), n_points, lowerb, upperb);
	llh_hists.push_back(hScan_total);
	
    xsecpar[par] = xsec->getNominal(par)-n_sigma*xsec->GetError(par);
 
	// Increment in sigma units
    double dsigma = (2*n_sigma)/((double)n_points-1);
	std::cout << "dsigma = " << dsigma << std::endl;

	// Loop over parameter values
    for (unsigned val =0; val < n_points; val++) {
      double totalllh = 0;
      
	  xsec->setParameters(xsecpar);
	  //std::cout << "Val = " << xsecpar[par] << std::endl;

	  // Calc LLH for each sample
      for(unsigned sample_i=0; sample_i < DUNEPdfs.size() ; ++sample_i) {
        DUNEPdfs[sample_i]  -> reweight();
	    llh_hists[sample_i]->Fill(xsecpar[par], 2 *  DUNEPdfs[sample_i]  -> GetLikelihood());	
        totalllh += 2 *  DUNEPdfs[sample_i] -> GetLikelihood();
		//std::cout << "LLH for sample " << sample_i << " = " << 2 *  DUNEPdfs[sample_i]  -> GetLikelihood() << std::endl;
      } //end of sample loop

	  // Calc Penalty LLH
	  llh_hists[DUNEPdfs.size()]->Fill(xsecpar[par], 2 * (xsec->GetLikelihood()));
	  // Total sample LLH
	  llh_hists[DUNEPdfs.size()+1]->Fill(xsecpar[par], 2 * totalllh);

      xsecpar[par] += dsigma*sqrt((*xsec->getCovMatrix())(par,par)); // Increment parameter
    } //end of value loop

    OutFile->cd();

	for(unsigned hist=0; hist < llh_hists.size(); hist++)
	{
      llh_hists[hist]->Write();
	}

    std::cout << "Finished xsec param " << par << std::endl;
  } //end of parameter loop

  std::cout << "Finished LLH Scans!" << std::endl;
  return 0;
}

