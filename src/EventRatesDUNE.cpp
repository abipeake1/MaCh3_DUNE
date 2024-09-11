#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

//template <typename T>
#include <assert.h>
#include <stdexcept>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>

#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"
#include "manager/manager.h"

using namespace std;
/* a function to generate numpy linspace */
//template <typename T>
  
template <typename T>
  std::vector<T> linspace(T a, T b, size_t N) {
      T h = (b - a) / static_cast<T>(N);
      std::vector<T> xs(N);
      typename std::vector<T>::iterator x;
      T val;
      for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
          *x = val;
      return xs;
  }

void ExtendLinspace(std::vector<double> &in, double low, double high,size_t nbins){
  for(auto x : linspace(low, high,nbins)){
      in.push_back(x);
      
    }
    auto x =std::unique(in.begin(),in.end());  //put duplicate bins at end
    in.erase(x, in.end()); //erase duplicates
}


// Constructors for erec-binned errors

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

  std::unique_ptr<TH3D> Abis3DHistogram_forplotting_unosc;
  std::unique_ptr<TH3D> Abis3DHistogram_forplotting_osc;
  std::unique_ptr<TH2D> Abis2DHistogram_forplotting_unosc;
  std::unique_ptr<TH2D> Abis2DHistogram_forplotting_unosc_scaled;
  std::unique_ptr<TH2D> Abis2DHistogram_forplotting_osc;

int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);

  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>(); 

  double FDPOT = fitMan->raw()["General"]["FDPOT"].as<double>(); 
  double NDPOT = fitMan->raw()["General"]["NDPOT"].as<double>(); 

  // Asimov fit
  bool asimovfit = false;//fitMan->GetAsimovFitFlag();
  
 
  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save plots
  std::string OutfileName = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>(); 
  TFile *Outfile = new TFile(OutfileName.c_str() , "RECREATE");

  //initialise xsec
  covarianceXsec *xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());

  std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  std::cout << "Cross section parameters:" << std::endl;
  xsec->printNominal();

  bool XsecParsAtGen = fitMan->raw()["General"]["Systematics"]["XsecAtGen"].as<bool>();

  std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;

  covarianceOsc *osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  std::vector<double> oscpars = fitMan->raw()["General"]["OscParams"].as<std::vector<double>>();

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;
  osc->setFlipDeltaM23(true);

  // Ask config file whether to use reactor constraint
  //std::cout << "use reactor prior is : " << useRC << std::endl ;

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);
  // This line gives a crash and stack trace...
  osc->setParameters(oscpars);
  osc->acceptStep();

  bool addFD = fitMan->raw()["General"]["IncludeFD"].as<bool>();
  bool addND = false;//fitMan->raw()["General"]["IncludeND"].as<bool>();//false;//fitMan->raw()["General"]["IncludeND"].as<bool>();  // false;

  if (!addFD && !addND) {std::cerr << "[ERROR:] You've chosen NOT to include FD or ND samples in the config file... you need to add something!" << std::endl; throw;}


  std::vector<samplePDFFDBase*> SamplePDFs;
  
  if(addFD) { 
    samplePDFDUNEBase *numu_oa = new samplePDFDUNEBase(NDPOT, "configs/SamplePDFDune_FHC_offaxis-nosplines.yaml", xsec);
    SamplePDFs.push_back(numu_oa);

    //samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_numuselec-nosplines.yaml", xsec);
    //SamplePDFs.push_back(numu_pdf);

    
    //samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_nueselec.yaml", xsec);
    //SamplePDFs.push_back(nue_pdf);
    // samplePDFDUNEBase *numubar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_numuselec.yaml", xsec);
    // SamplePDFs.push_back(numubar_pdf);
    // samplePDFDUNEBase *nuebar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_nueselec.yaml", xsec);
    // SamplePDFs.push_back(nuebar_pdf);
  }
 
  if(addND) {
    //samplePDFDUNEBaseND *numu_pdf_oa = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_FHC_CCnumuselec-nosplines.yaml", xsec);
    //SamplePDFs.push_back(numu_pdf_oa);
    // * FHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_FHC_CCnumuselec.yaml", xsec);
    //SamplePDFs.push_back(FHC_numuCCND_pdf);
    //samplePDFDUNEBaseND * RHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_RHC_CCnumuselec.yaml", xsec);
    //SamplePDFs.push_back(RHC_numuCCND_pdf);
  }

  // Oscillated
  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	    << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	    << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	    << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	    << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	    << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;


  // unosc
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0;
  oscpars_un[1] = 0;
  oscpars_un[2] = 0;
  oscpars_un[3] = 0;
  oscpars_un[4] = 0;

  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  if(XsecParsAtGen){
	TFile* XsecFile = new TFile(XsecMatrixFile.c_str(), "READ");
	TVectorD* XsecGeneratedParamArray = (TVectorD*)XsecFile->Get("xsec_param_nom");
	std::cout << "Setting xsec systs to their generated values " << std::endl;
	for (unsigned param_i = 0 ; param_i < XsecParVals.size() ; ++param_i) {
	  std::cout << "Generated value for param " << param_i << " is " << (*XsecGeneratedParamArray)(param_i) << std::endl;
	  XsecParVals[param_i] = (*XsecGeneratedParamArray)(param_i);
	  std::cout << "Set parameter " << param_i << " to value " << XsecParVals[param_i] << std::endl;
	}
  }
  else{
	std::cout << "Keeping xsec parameters at their prior values" << std::endl;
  }

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());




  //////////////////////////////////////Stuff for the 3D Histogram unfolding :)
  std::vector<double> ELep_bins;
    ExtendLinspace(ELep_bins,0,40,5);
    std::vector<double> theta_bins;
    ExtendLinspace(theta_bins,0,3,10);
    ExtendLinspace(theta_bins,3,9,3);
    ExtendLinspace(theta_bins,9,100,2);
    
    //std::vector<double> ENuReco_bins;
     std::vector<double> ENuReco_bins={0.,   1.,  1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 10.};
    //ExtendLinspace(ENuReco_bins,0,4,8);

    std::vector<double> offaxis_position;
    ExtendLinspace(offaxis_position,-1,30,10);

    
    
//TH3D* Abis3DHistogram_forplotting = new TH3D("Abis3DHistogramforplotting", "", ELep_bins.size(), theta_bins.size(),  ENuReco_bins.size()); // fill this in with the binning like we used to in NUISANCE
  //std::cout<< "ELep_bins.size() = " << ELep_bins.size() <<std::endl;
  //std::cout<< "theta_bins.size() = " << theta_bins.size() <<std::endl;
  //std::cout<< "ENuReco_bins.size() = " << ENuReco_bins.size() <<std::endl;



  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
 

  std::vector<std::string> sample_names;
  
  TCanvas *canv = new TCanvas("nomcanv","", 0, 0, 700,900);
  canv->Divide(2,3);
  std::string OutPlotName = OutfileName.substr(0, OutfileName.length() - 5) + ".pdf";
  canv->Print((OutPlotName+"[").c_str());

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {

    Abis3DHistogram_forplotting_unosc  = std::make_unique<TH3D>("Abis3DHistogram_forplotting_unosc ", "", ELep_bins.size() - 1,
                                         ELep_bins.data(), theta_bins.size() - 1, theta_bins.data(),
                                         ENuReco_bins.size() - 1, ENuReco_bins.data());

    Abis3DHistogram_forplotting_osc  = std::make_unique<TH3D>("Abis3DHistogram_forplotting_osc ", "", ELep_bins.size() - 1,
                                        ELep_bins.data(), theta_bins.size() - 1, theta_bins.data(),
                                        ENuReco_bins.size() - 1, ENuReco_bins.data());


    Abis2DHistogram_forplotting_unosc  = std::make_unique<TH2D>("Abis2DHistogram_forplotting_unosc ", "", offaxis_position.size() - 1,
                                         offaxis_position.data(), ENuReco_bins.size() - 1, ENuReco_bins.data());

    Abis2DHistogram_forplotting_osc  = std::make_unique<TH2D>("Abis2DHistogram_forplotting_osc ", "", offaxis_position.size() - 1,
                                        offaxis_position.data(), ENuReco_bins.size() - 1, ENuReco_bins.data());

    Abis2DHistogram_forplotting_unosc_scaled  = std::make_unique<TH2D>("Abis2DHistogram_forplotting_unos_scaled ", "", offaxis_position.size() - 1,
                                         offaxis_position.data(), ENuReco_bins.size() - 1, ENuReco_bins.data());




    std::string name = SamplePDFs[sample_i]->GetSampleName();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    // Unoscillated
    osc -> setParameters(oscpars_un);
    osc -> acceptStep();
    SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
    SamplePDFs[sample_i] -> reweight(osc->getPropPars());
    TH1D *sample_unosc = (TH1D*)SamplePDFs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");

  //Then in your plotting script, you make the same 3D histogram and loop over the bins of the 1D histogram and just do:
    for(Int_t bin1d_i = 0; bin1d_i < sample_unosc->GetNcells(); ++bin1d_i){
      Abis3DHistogram_forplotting_unosc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
    }

	unoscillated_hists.push_back(sample_unosc);
  //unoscillated_hists.push_back(sample_unosc);

  canv->cd(1);
	sample_unosc -> SetTitle(NameTString+"_unosc");
  sample_unosc -> GetXaxis()->SetTitle("Bin Number");
	sample_unosc -> Draw("HIST");
	canv->Update();

	Outfile->cd();
	sample_unosc->Write(NameTString+"_unosc");

	osc->setParameters(oscpars);
	osc->acceptStep();

	SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH1D *sample_osc = (TH1D*)SamplePDFs[sample_i] -> get1DHist()->Clone(NameTString+"_osc");
	oscillated_hists.push_back(sample_osc);

  
	canv->cd(2);
	sample_osc -> SetTitle(NameTString+"_osc");
  sample_osc -> GetXaxis()->SetTitle("Bin Number");
	sample_osc -> Draw("HIST");
	canv->Update();


   canv->cd(3);

  //Then in your plotting script, you make the same 3D histogram and loop over the bins of the 1D histogram and just do:
    for(Int_t bin1d_i = 0; bin1d_i < sample_osc->GetNcells(); ++bin1d_i){
      Abis3DHistogram_forplotting_osc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
    }
     Abis3DHistogram_forplotting_osc->Draw("lego colz");
     canv->Update();

     canv->cd(4);
    
     Abis3DHistogram_forplotting_unosc->Draw("lego colz");
     canv->Update();



     canv->cd(5);
     /*TH2D *sample_unosc2D = (TH2D*)SamplePDFs[sample_i] -> get2DHist() -> Clone(NameTString+"_unosc");
     sample_unosc2D->Scale(1., "width");
     sample_unosc2D->Draw("colz");
     sample_unosc2D->Write(NameTString+"unosc_2D");
     canv->Update();
     */
     //TH2D *sample_osc2D = (TH2D*)SamplePDFs[sample_i] -> get2DHist() -> Clone(NameTString+"_osc");
     //sample_osc2D->Draw("colz");
     for(Int_t bin1d_i = 0; bin1d_i < sample_unosc->GetNcells(); ++bin1d_i){
      Abis2DHistogram_forplotting_unosc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
      }
    // Abis2DHistogram_forplotting_unosc->SetOptStat(0);
     Abis2DHistogram_forplotting_unosc->GetXaxis()->SetTitle("Off axis position");
     Abis2DHistogram_forplotting_unosc->GetYaxis()->SetTitle("Neutrino energy (GeV)");
     Abis2DHistogram_forplotting_unosc->Draw("colz");
     

     canv->cd(6);
     
     for(Int_t bin1d_i = 0; bin1d_i < sample_unosc->GetNcells(); ++bin1d_i){
      Abis2DHistogram_forplotting_unosc_scaled->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
      }
    // Abis2DHistogram_forplotting_unosc->SetOptStat(0);
     Abis2DHistogram_forplotting_unosc_scaled->GetXaxis()->SetTitle("Off axis position");
     Abis2DHistogram_forplotting_unosc_scaled->GetYaxis()->SetTitle("Neutrino energy (GeV)");

     
    //Divide histogram rows by number of events in each OA position------------------------------------------------------------------
    std::vector<double> bincontent_oa_position;
    for(Int_t det_i = 0; det_i < Abis2DHistogram_forplotting_unosc_scaled->GetNbinsX(); ++det_i){
      double n_events_yaxis=1.0;
      for(int x=0 ; x < Abis2DHistogram_forplotting_unosc_scaled->GetNbinsY(); ++x){

        double bin_content = Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x);
        std::cout<<"bin content" << bin_content <<std::endl;
        n_events_yaxis += bin_content;
      }
      std::cout<<"events in bin " << det_i << "=" << n_events_yaxis <<std::endl;
      bincontent_oa_position.push_back(n_events_yaxis);
      }

      for(Int_t det_i = 0; det_i < Abis2DHistogram_forplotting_unosc_scaled->GetNbinsX(); ++det_i){
      
        for(int x=0 ; x < Abis2DHistogram_forplotting_unosc_scaled->GetNbinsY(); ++x){
            //std::cout << "Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x)" << Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x) <<std::endl;
            Abis2DHistogram_forplotting_unosc_scaled->SetBinContent(det_i,x, Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x)/bincontent_oa_position[det_i]);
            //std::cout << "Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x)/bincontent_oa_position[det_i]" << Abis2DHistogram_forplotting_unosc_scaled->GetBinContent(det_i,x)/bincontent_oa_position[det_i] <<std::endl;
        
        }
      //Abis2DHistogram_forplotting_unosc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
      }
      Abis2DHistogram_forplotting_unosc_scaled->Draw("colz");

     /*
    for(Int_t bin1d_i = 0; bin1d_i < sample_osc->GetNcells(); ++bin1d_i){
      Abis2DHistogram_forplotting_osc->SetBinContent(bin1d_i,sample_osc->GetBinContent(bin1d_i));
    }

     //Abis2DHistogram_forplotting_osc->SetOptStat(0);
     Abis2DHistogram_forplotting_osc->GetXaxis()->SetTitle("Off axis position");
     Abis2DHistogram_forplotting_osc->GetYaxis()->SetTitle("Neutrino energy (GeV)");
     Abis2DHistogram_forplotting_osc->Draw("colz");
     */
	   canv->Print(OutPlotName.c_str());

	Outfile->cd();
	sample_osc->Write(NameTString+"_osc");
  Abis2DHistogram_forplotting_unosc->Write(NameTString+"unosc_2D");

    
  Abis2DHistogram_forplotting_unosc_scaled->Write(NameTString+"unosc_2D_scaled");

  ///////////////////////////////////////////////////////////
  }

  canv->Print((OutPlotName+"]").c_str());
  Outfile->Close();
// #define DEBUG

   
  //Now print out some event rates, we'll make a nice latex table at some point 
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Integrals of nominal hists: " << std::endl;
  for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	std::cout << " " << std::endl;
	std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
	std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
  }
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;

  return 0;
 }