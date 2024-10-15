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
#include "samplePDFDUNE/histograms.h"

using namespace std;

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

  // Access histograms
    Multidim_HistogramManager histograms;
    auto& Abis3DHistogram_forplotting_unosc = histograms.Abis3DHistogram;
    auto& Abis3DHistogram_forplotting_osc = histograms.Abis3DHistogram;
    //auto& oa3D_forplotting_unosc = histograms.OA3DHistogram;
    //auto& oa3D_forplotting_osc = histograms.OA3DHistogram
    auto& Abis2DHistogram_forplotting_unosc  = histograms.OA2DHistogram;
    auto& Abis2DHistogram_osc = histograms.OA2DHistogram;
    auto& Abis2DHistogram_forplotting_unosc_scaled= histograms.OA2DHistogram;
    auto& Abis2DHistogram_forplotting_osc =histograms.OA2DHistogram;


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

  //Setup the cross-section parameters
  //This should get the prior values.
  std::cout<< "xsec->getNominalArray()" <<std::endl;
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
  xsec->setParameters();
  }
  
  if(addFD) { 
    samplePDFDUNEBase *numu_oa = new samplePDFDUNEBase(NDPOT, "/home/abipeake/Mach3/MaCh3_DUNE/configs/SamplePDFDune_FHC_offaxis-nosplines.yaml", xsec); //"configs/SamplePDFDune_FHC_offaxis-nosplines.yaml"
    SamplePDFs.push_back(numu_oa);
  }
  if(addND) {
  
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
  
  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
 

  std::vector<std::string> sample_names;
  
  TCanvas *canv = new TCanvas("nomcanv","", 0, 0, 700,900);
   
  canv->Divide(2,3);
   std::cout << "OutfileName  =  "  << OutfileName <<std::endl;
  std::string OutPlotName =  OutfileName.substr(0, OutfileName.length() - 5) + ".pdf";
  std::cout << "OutPlotName  =  "  << OutPlotName<<std::endl;
  canv->Print((OutPlotName+"[").c_str());

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {
    std::string name = SamplePDFs[sample_i]->GetSampleName();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    // Unoscillated
    osc -> setParameters(oscpars_un);
    osc -> acceptStep();
    //std::cout<<"acceptStep1" <<std::endl;
    SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
    //std::cout<<"acceptStep2" <<std::endl;
    SamplePDFs[sample_i] -> reweight(osc->getPropPars());

    ////////////////////
    xsec->setParameters(XsecParVals);
    double nominal =xsec->getNominal(0); //get central value of parameter
    double error = xsec->getDiagonalError(0);
    //xsec->setParCurrProp(0, nominal+(2*error));////////// set 
    xsec->setParCurrProp(0, nominal);////////// set //+(2*error)
    std::cout<< "nominal value = " << nominal <<std::endl;;
    std::cout<< "error = " << error <<std::endl;
   // xsec->setSingleParameter(0, nominal); //nominal+(2*error)
    double current_value = xsec->getParProp(0);
    std::cout<<"current value  = " << current_value << std::endl; 
    //xsec->setParCurrProp(0, nominal);
    xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
    SamplePDFs[sample_i] -> reweight(osc->getPropPars());
    //SamplePDFs[sample_i] -> xsec->setParameters();
    
   
    TH1D *sample_unosc = (TH1D*)SamplePDFs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
	  unoscillated_hists.push_back(sample_unosc);

    canv->cd(1);
    sample_unosc -> SetTitle(NameTString+"_unosc_0sigma");
    sample_unosc -> GetXaxis()->SetTitle("Bin Number");
    sample_unosc -> Draw("HIST");
    canv->Update();

    Outfile->cd();
    sample_unosc->Write(NameTString+"_unosc");

    //osc->setParameters(oscpars);
    //osc->acceptStep();
    //SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	

    canv->cd(3);
    

    for(Int_t bin1d_i = 0; bin1d_i < sample_unosc->GetNcells(); ++bin1d_i){
      Abis3DHistogram_forplotting_unosc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
    }
    Abis3DHistogram_forplotting_unosc->Draw("colz");
    Abis3DHistogram_forplotting_unosc->Write();
    canv->Update();
      
     canv->cd(5);
     for(Int_t bin1d_i = 0; bin1d_i < sample_unosc->GetNcells(); ++bin1d_i){
      Abis2DHistogram_forplotting_unosc->SetBinContent(bin1d_i,sample_unosc->GetBinContent(bin1d_i));
      }
     Abis2DHistogram_forplotting_unosc->SetTitle("Nominal");
     Abis2DHistogram_forplotting_unosc->Draw("colz");
     Abis2DHistogram_forplotting_unosc->Write(NameTString+"unosc_2D");
     canv->Update();


     canv->cd(2);

    //Now shift the xsec parameter by +2 sigma
      xsec->setParameters(XsecParVals);
      double nominal2 =xsec->getNominal(0); //get central value of parameter
      double error2 = xsec->getDiagonalError(0);
      //xsec->setParCurrProp(0, nominal+(2*error));////////// set 
      xsec->setParCurrProp(0, nominal2+(2*error2));////////// set //+(2*error)
      std::cout<< "nominal value = " << nominal2 <<std::endl;;
      std::cout<< "error = " << error2 <<std::endl;
      double current_value2 = xsec->getParProp(0);
      std::cout<<"current value  = " << current_value2 << std::endl; 
      xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
      SamplePDFs[sample_i] -> reweight(osc->getPropPars());


      TH1D *sample_osc = (TH1D*)SamplePDFs[sample_i] -> get1DHist()->Clone(NameTString+"_osc");
      oscillated_hists.push_back(sample_osc);
      Outfile->cd();
      sample_osc->Write(NameTString+"_osc");
      sample_osc -> SetTitle(NameTString+"_osc_2sigma");
      sample_osc -> GetXaxis()->SetTitle("Bin Number");
      sample_osc -> Draw("HIST");
      canv->Update();

     canv->cd(4);
     xsec->setParameters(XsecParVals);
      //double nominal2 =xsec->getNominal(0); //get central value of parameter
      //double error2 = xsec->getDiagonalError(0);
      //xsec->setParCurrProp(0, nominal+(2*error));////////// set 
     xsec->setParCurrProp(0, nominal2+(2*error2));////////// set //+(2*error)
      //std::cout<< "nominal value = " << nominal2 <<std::endl;;
      //std::cout<< "error = " << error2 <<std::endl;
      //double current_value2 = xsec->getParProp(0);
      //std::cout<<"current value  = " << current_value2 << std::endl; 
     xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
     SamplePDFs[sample_i] -> reweight(osc->getPropPars());


     for(Int_t bin1d_i = 0; bin1d_i < sample_osc->GetNcells(); ++bin1d_i){
      Abis3DHistogram_forplotting_osc->SetBinContent(bin1d_i,sample_osc->GetBinContent(bin1d_i));
    }
     Abis3DHistogram_forplotting_osc->Draw("lego");
     Abis3DHistogram_forplotting_osc->Write();
     canv->Update();

  
     canv->cd(6);
     xsec->setParameters(XsecParVals);
      //double nominal2 =xsec->getNominal(0); //get central value of parameter
      //double error2 = xsec->getDiagonalError(0);
      //xsec->setParCurrProp(0, nominal+(2*error));////////// set 
      xsec->setParCurrProp(0, nominal2+(2*error2));////////// set //+(2*error)
      //std::cout<< "nominal value = " << nominal2 <<std::endl;;
      //std::cout<< "error = " << error2 <<std::endl;
      //double current_value2 = xsec->getParProp(0);
      //std::cout<<"current value  = " << current_value2 << std::endl; 
      xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
      SamplePDFs[sample_i] -> reweight(osc->getPropPars());


       for(Int_t bin1d_i = 0; bin1d_i < sample_osc->GetNcells(); ++bin1d_i){
        Abis2DHistogram_forplotting_osc->SetBinContent(bin1d_i,sample_osc->GetBinContent(bin1d_i));
      }
     
       Abis2DHistogram_forplotting_osc->SetTitle("+2 #sigma");
       Abis2DHistogram_forplotting_osc->Draw("colz");
      Abis2DHistogram_forplotting_osc->Write(NameTString+"osc_2D");
      canv->Update();
      
	   canv->Print(OutPlotName.c_str());

	Outfile->cd();
	
 


  ///////////////////////////////////////////////////////////
  }

  canv->Print((OutPlotName+"]").c_str());
  Outfile->Close();
// #define DEBUG

   
  //Now print out some event rates, we'll make a nice latex table at some point 
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "Integrals of nominal hists: " << std::endl;
  for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
  std::cout<< "help"<< std::endl;
	std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
  }
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;

  return 0;


  /*
    for(Int_t bin1d_i = 0; bin1d_i < sample_osc->GetNcells(); ++bin1d_i){
      Abis2DHistogram_forplotting_osc->SetBinContent(bin1d_i,sample_osc->GetBinContent(bin1d_i));
    }

     //Abis2DHistogram_forplotting_osc->SetOptStat(0);
     Abis2DHistogram_forplotting_osc->GetXaxis()->SetTitle("Off axis position");
     Abis2DHistogram_forplotting_osc->GetYaxis()->SetTitle("Neutrino energy (GeV)");
     Abis2DHistogram_forplotting_osc->Draw("colz");
     */

 }