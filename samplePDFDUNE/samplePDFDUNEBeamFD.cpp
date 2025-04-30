#include <TROOT.h>

#include "samplePDFDUNEBeamFD.h"
#include "TString.h"
#include <assert.h>
#include <stdexcept>
#include "TMath.h"
#include "manager/manager.h"

#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "TSystem.h"


#include <cstdlib>
#include <iostream>



samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_, covarianceXsec* xsec_cov_) : samplePDFFDBase(mc_version_, xsec_cov_) {
  //Call insitialise in samplePDFFD
  //std::cout << "starting samplePDFDUNEBeamFD........"<<std::endl;
  //gSystem->Setenv("XROOTD_CLIENT_VERSION", "5.6.9");

  Initialise();
  //std::cout << "sdone initialisE()........"<<std::endl;
}


samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {
}


void samplePDFDUNEBeamFD::Init() {
  //std::cout<<"in Init()"<< std::endl;
  dunemcSamples.resize(nSamples,dunemc_base());
  //double pot = 1.3628319e+23;
  //events_for_analysis = SampleManager->raw()["NEvents"].as<double>();
  pot = SampleManager->raw()["POT"].as<double>();
  //pot_cafgenwith = SampleManager->raw()["GEN_POT"].as<double>(); /////////////just as test for fermigrid
  std::cout<< "POT IS !!! ================ " << pot <<std::endl;
 
  if (CheckNodeExists(SampleManager->raw(), "DUNESampleBools", "iselike" )) {
    iselike = SampleManager->raw()["DUNESampleBools"]["iselike"].as<bool>();
  } else{
    MACH3LOG_ERROR("Did not find DUNESampleBools:iselike in {}, please add this", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  /*
  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else{
    MACH3LOG_ERROR("POT not defined in {}, please add this!", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }*/
  
  tot_escale_fd_pos = -999;
  tot_escale_sqrt_fd_pos = -999;
  tot_escale_invsqrt_fd_pos = -999;
  had_escale_fd_pos = -999;
  had_escale_sqrt_fd_pos = -999;
  had_escale_invsqrt_fd_pos = -999;
  mu_escale_fd_pos = -999;
  mu_escale_sqrt_fd_pos = -999;
  mu_escale_invsqrt_fd_pos = -999;
  n_escale_fd_pos = -999;
  n_escale_sqrt_fd_pos = -999;
  n_escale_invsqrt_fd_pos = -999;
  em_escale_fd_pos = -999;
  em_escale_sqrt_fd_pos = -999;
  em_escale_invsqrt_fd_pos = -999;
  had_res_fd_pos = -999;
  mu_res_fd_pos = -999;
  n_res_fd_pos = -999;
  em_res_fd_pos = -999;
  cvn_numu_fd_pos = -999;
  cvn_nue_fd_pos = -999;

  // create dunemc storage
  dunemcSamples.resize(nSamples);

  nFDDetectorSystPointers = funcParsIndex.size();
  std::unordered_map<std::string, const double*> FDDetectorSystPointersMap;
  FDDetectorSystPointers = std::vector<const double*>(nFDDetectorSystPointers);

  for(auto FuncPar_i  = 0 ; FuncPar_i < funcParsIndex.size() ; ++FuncPar_i){
    FDDetectorSystPointersMap.insert(std::pair<std::string, const double*>(funcParsNames.at(FuncPar_i), XsecCov->retPointer(funcParsIndex.at(FuncPar_i))));
  }

  /*
  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);
    
    if (name == "TotalEScaleFD") {
      tot_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_fd_pos);
    }
    else if (name == "TotalEScaleSqrtFD") {
      tot_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_sqrt_fd_pos);
    }
    else if (name == "TotalEScaleInvSqrtFD") {
      tot_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_invsqrt_fd_pos);
    }
    else if (name == "HadEScaleFD") {
      had_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_fd_pos);
    }
    else if (name == "HadEScaleSqrtFD") {
      had_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_sqrt_fd_pos);
    }
    else if (name == "HadEScaleInvSqrtFD") {
      had_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_invsqrt_fd_pos);
    }
    else if (name == "MuEScaleFD") {
      mu_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_fd_pos);
    }
    else if (name == "MuEScaleSqrtFD") {
      mu_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_sqrt_fd_pos);
    }
    else if (name == "MuEScaleInvSqrtFD") {
      mu_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_invsqrt_fd_pos);
    }
    else if (name == "NEScaleFD") {
      n_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_fd_pos);
    }
    else if (name == "NEScaleSqrtFD") {
      n_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_sqrt_fd_pos);
    }
    else if (name == "NEScaleInvSqrtFD") {
      n_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_invsqrt_fd_pos);
    }
    else if (name == "EMEScaleFD") {
      em_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_fd_pos);
    }
    else if (name == "EMEScaleSqrtFD") {
      em_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_sqrt_fd_pos);
    }
    else if (name == "EMEScaleInvSqrtFD") {
      em_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_invsqrt_fd_pos);
    }
    else if (name == "HadResFD") {
      had_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(had_res_fd_pos);
    }
    else if (name == "MuResFD") {
      mu_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_res_fd_pos);
    }
    else if (name == "NResFD") {
      n_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(n_res_fd_pos);
    }
    else if (name == "EMResFD") {
      em_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(em_res_fd_pos);
    }
    else if (name == "CVNNumuFD") {
      cvn_numu_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_numu_fd_pos);
    }
    else if (name == "CVNNueFD") {
      cvn_nue_fd_pos = *it;
      FDDetectorSystPointers[func_it] = XsecCov->retPointer(cvn_nue_fd_pos);
    }
    
    else { 
      std::cerr << "Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamFD:" << name << std::endl;
      throw;
    }
  }
  */
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
}

void samplePDFDUNEBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if(spline_files.size() && (XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0)){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = new splinesDUNE(XsecCov);
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }
  
  return;
}

double samplePDFDUNEBeamFD::CalculatePOT() {
  TChain calc_pot_chain("meta");  // Use correct tree name

  std::string pot_branch = "pot";  // Use the correct branch name

  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
      calc_pot_chain.AddFile((mc_files[i]).c_str());
  }

  // Check if the branch exists before proceeding
  if (!calc_pot_chain.GetBranch(pot_branch.c_str())) {
      std::cerr << "Error: Branch " << pot_branch << " not found in the tree!" << std::endl;
      return 0.0;
  }

  double pot_value = 0.0;
  calc_pot_chain.SetBranchAddress(pot_branch.c_str(), &pot_value);

  double sum_pot = 0.0;
  Long64_t nEntries = calc_pot_chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
      calc_pot_chain.GetEntry(i);
      sum_pot += pot_value;
  }

  std::cout << "Summed POT: " << sum_pot << std::endl;
  return sum_pot;
}



void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j] = new const double*[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunemcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}


  

int samplePDFDUNEBeamFD::setupExperimentMC(int iSample) {

  std::cout<< " in experiment setupMC " << std::endl; 

  auto &duneobj = dunemcSamples[iSample];

  int nutype = sample_nutype[iSample];
  int oscnutype = sample_oscnutype[iSample];
  bool signal = sample_signal[iSample];

  double newpot = CalculatePOT();
  std::cout<< "CalculatePOT() = " << newpot << " -----------------------------------------------------------------------------" << std::endl;

  
  //std::cout << "BEARER_TOKEN_FILE: " << getenv("BEARER_TOKEN_FILE") << std::endl;
  //std::cout << "XRD_SEC_TOKEN: " << getenv("XRD_SEC_TOKEN") << std::endl;
  //std::cout << "XRD_PLUGIN: " << getenv("XRD_PLUGIN") << std::endl;
  //gSystem->Setenv("XrdSecPROTOCOL", "bearer");
  //gSystem->Setenv("XRD_RUNFORKHANDLER", "1");


  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files[iSample].native());
  std::cout<< " about to do  TFile::Open " << std::endl; 
  //_sampleFile = new TFile(mc_files[iSample].c_str(), "READ");
  //_sampleFile = TFile::Open(mc_files[iSample].c_str(), "READ"); //for xrootd reading apparently
  std::string filename = mc_files[iSample].string();
  _sampleFile = TFile::Open(filename.c_str(), "READ");
  if (!_sampleFile || _sampleFile->IsZombie()) {
    std::cout << "Error: Failed to open file " << filename << std::endl;
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  std::cout << "Successfully read in sample file " << filename << std::endl;
  

  _data = (TTree*)_sampleFile->Get("cafTree");
  if (!_data) {
      MACH3LOG_ERROR("Branch 'cafTree' not found in file: {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
  }
  //_data = (TTree*)_sampleFile->Get("caf");
  _meta = (TTree*)_sampleFile->Get("meta");
  if (!_meta ){
      MACH3LOG_ERROR("Branch 'meta' not found in file: {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
  } 

  //_sampleFile = TFile::Open(mc_files[iSample].c_str(), "READ");
  

  if (!_sampleFile || _sampleFile->IsZombie()) {
      std::cerr << "Error: Failed to open file " << mc_files[iSample] << std::endl;
  }
  _data = (TTree*)_sampleFile->Get("caf");
  if (!_data) {
      std::cerr << "Error: 'caf' tree not found in " << mc_files[iSample] << std::endl;
  }

  _meta = (TTree*)_sampleFile->Get("meta");
  if (!_meta) {
      std::cerr << "Error: 'meta' tree not found in " << mc_files[iSample] << std::endl;
  }
  //if (!caf){ caf = tfile->Get<TTree>(“cafTree”); }

  //TTree *_metaTree = (TTree*)_sampleFile->Get("meta")

  
  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample].native());
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    _data = (TTree*)_sampleFile->Get("cafTree");
    
  }
  
  double _production_pot = 0.0;  // Explicitly initialize
  double gen_pot = 0.0; //set the sum of the pot from each file to be 0,befor any are read in
  
  if(_meta){
    std::cout<<"Found meta tree" << std::endl;
     
    _meta->SetBranchStatus("*", 0);
    _meta->SetBranchStatus("pot", 1);
    _meta->SetBranchAddress("pot", &_production_pot);
   
    for (int i = 0; i < _meta->GetEntries(); i++) {
      //std::cout << "Total entries in _meta: " << _meta->GetEntries() << std::endl;

      _meta->GetEntry(i);
      //std::cout << "Entry " << i << " pot: " << _production_pot << std::endl;
      if (_production_pot > 1e30) {  // Arbitrary threshold for large numbers
        std::cerr << "WARNING: Unusually large pot value at entry " << i << ": " 
                  << _production_pot << std::endl;
      }

      //gen_pot += _production_pot;
      gen_pot = gen_pot + _production_pot;
      //std::cout << "updated gen_pot = " << gen_pot << std::endl;
  }

  }
  else{
    MACH3LOG_ERROR("Could not find \"meta\" tree in {}", mc_files[iSample].native());
    gen_pot =  3.85e21;
    std::cout << "instead setting gen_pot to be = " << gen_pot << std::endl;
  }
  
  std::cout << " final gen_pot  = " << gen_pot << std::endl;

  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco_numu", 1);
  _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("Ev_reco", 1);
  _data->SetBranchAddress("Ev_reco", &_erec);
  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  _data->SetBranchStatus("RecoLepEnNumu", 1);
  _data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
  _data->SetBranchStatus("RecoLepEnNue", 1);
  _data->SetBranchAddress("RecoLepEnNue", &_erec_lep_nue);

  _data->SetBranchStatus("eRecoP", 1);
  _data->SetBranchAddress("eRecoP", &_eRecoP);
  _data->SetBranchStatus("eRecoPip", 1);
  _data->SetBranchAddress("eRecoPip", &_eRecoPip);
  _data->SetBranchStatus("eRecoPim", 1);
  _data->SetBranchAddress("eRecoPim", &_eRecoPim);
  _data->SetBranchStatus("eRecoPi0", 1);
  _data->SetBranchAddress("eRecoPi0", &_eRecoPi0);
  _data->SetBranchStatus("eRecoN", 1);
  _data->SetBranchAddress("eRecoN", &_eRecoN);

  _data->SetBranchStatus("LepE", 1);
  _data->SetBranchAddress("LepE", &_LepE);
  _data->SetBranchStatus("eP", 1);
  _data->SetBranchAddress("eP", &_eP);
  _data->SetBranchStatus("ePip", 1);
  _data->SetBranchAddress("ePip", &_ePip);
  _data->SetBranchStatus("ePim", 1);
  _data->SetBranchAddress("ePim", &_ePim);
  _data->SetBranchStatus("ePi0", 1);
  _data->SetBranchAddress("ePi0", &_ePi0);
  _data->SetBranchStatus("eN", 1);
  _data->SetBranchAddress("eN", &_eN);

  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("cvnnumu",1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue",1);
  _data->SetBranchAddress("cvnnue", &_cvnnue);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
  _data->SetBranchStatus("vtx_x", 1);
  _data->SetBranchAddress("vtx_x", &_vtx_x);
  _data->SetBranchStatus("vtx_y", 1);
  _data->SetBranchAddress("vtx_y", &_vtx_y);
  _data->SetBranchStatus("vtx_z", 1);
  _data->SetBranchAddress("vtx_z", &_vtx_z);

  _data->SetBranchStatus("NuMomX", 1);
  _data->SetBranchAddress("NuMomX", &_NuMomX);
  _data->SetBranchStatus("NuMomY", 1);
  _data->SetBranchAddress("NuMomY", &_NuMomY);
  _data->SetBranchStatus("NuMomZ", 1);
  _data->SetBranchAddress("NuMomZ", &_NuMomZ);
  _data->SetBranchStatus("LepMomX", 1);
  _data->SetBranchAddress("LepMomX", &_LepMomX);
  _data->SetBranchStatus("LepMomY", 1);
  _data->SetBranchAddress("LepMomY", &_LepMomY);
  _data->SetBranchStatus("LepMomZ", 1);
  _data->SetBranchAddress("LepMomZ", &_LepMomZ);

  _data->SetBranchStatus("LepNuAngle", 1);
  _data->SetBranchAddress("LepNuAngle", &_LepNuAngle);

  _data->SetBranchStatus("Elep_reco", 1);
  _data->SetBranchAddress("Elep_reco", &_erec_lep);

  /* TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
    throw MaCh3Exception(__FILE__, __LINE__);
    }
  */
   

  //meta->SetBranchAddress("pot", &_production_pot);

  
//double mean_pot = sum_pot / meta->GetEntries();
std::cout << "Generated POT: " << gen_pot << std::endl;
std::cout << "pot_s: " << (pot)  << std::endl; /// (newpot)

   TH1D* norm = (TH1D*)_sampleFile->Get("norm");
  if(!norm){
    //MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
    //throw MaCh3Exception(__FILE__, __LINE__);
  
    norm = new TH1D("norm","",1,0,1);
    norm->SetBinContent(1,1);
    duneobj.norm_s = 1.0; //norm->GetBinContent(1);
    //duneobj->pot_s = (pot) / norm->GetBinContent(1);
    //duneobj.pot_s = (pot) / 3.85e21;

    //duneobj->norm_s = 1.0/gen_pot;
    //duneobj->pot_s = *pot; /

    duneobj.pot_s = (pot)/(newpot) ; ////////////correct one for all CAFS
    //duneobj.pot_s = (pot)/pot_cafgenwith; //one for subsample of OA
  
     std::cout << "(pot)/(newpot)" << (pot)/(newpot) << std::endl;
     //std::cout << "(pot)/(pot_cafgenwith)" << (pot)/(pot_cafgenwith) << std::endl;

    }
  else{
    duneobj.norm_s = norm->GetBinContent(1);
    duneobj.pot_s = pot/norm->GetBinContent(2);
  }

  std::cout << "pot = " << (pot)<< std::endl;
  std::cout << "pot_s = " << duneobj.pot_s << std::endl;
  std::cout << "norm_s = " << duneobj.norm_s << std::endl;
  std::cout << "N.o of events = " <<  _data->GetEntries() << std::endl;
  
  
  //if(events_for_analysis  >  _data->GetEntries()){
    //std::cout<<"you are trying to analyse with an number of events larger than that in the sample" << std::endl;
  //}
  //duneobj.norm_s = norm->GetBinContent(1);
  //duneobj.pot_s = pot/norm->GetBinContent(2);
  duneobj.nEvents = _data->GetEntries();
  duneobj.nutype = nutype;
  duneobj.oscnutype = oscnutype;
  duneobj.signal = signal;

  // allocate memory for dunemc variables
  duneobj.rw_cvnnumu.resize(duneobj.nEvents);
  duneobj.rw_cvnnue.resize(duneobj.nEvents);
  duneobj.rw_cvnnumu_shifted.resize(duneobj.nEvents);
  duneobj.rw_cvnnue_shifted.resize(duneobj.nEvents);
  duneobj.rw_etru.resize(duneobj.nEvents);
  duneobj.rw_erec.resize(duneobj.nEvents);
  duneobj.rw_erec_shifted.resize(duneobj.nEvents);
  duneobj.rw_erec_had.resize(duneobj.nEvents);
  duneobj.rw_erec_lep.resize(duneobj.nEvents);

  duneobj.true_q0.resize(duneobj.nEvents);
  duneobj.true_q3.resize(duneobj.nEvents);

  duneobj.rw_eRecoP.resize(duneobj.nEvents);
  duneobj.rw_eRecoPip.resize(duneobj.nEvents);
  duneobj.rw_eRecoPim.resize(duneobj.nEvents);
  duneobj.rw_eRecoPi0.resize(duneobj.nEvents);
  duneobj.rw_eRecoN.resize(duneobj.nEvents);

  duneobj.rw_LepE.resize(duneobj.nEvents);
  duneobj.rw_eP.resize(duneobj.nEvents);
  duneobj.rw_ePip.resize(duneobj.nEvents);
  duneobj.rw_ePim.resize(duneobj.nEvents);
  duneobj.rw_ePi0.resize(duneobj.nEvents);
  duneobj.rw_eN.resize(duneobj.nEvents);

  duneobj.rw_theta.resize(duneobj.nEvents);
  duneobj.flux_w.resize(duneobj.nEvents);
  duneobj.rw_isCC.resize(duneobj.nEvents);
  duneobj.rw_nuPDGunosc.resize(duneobj.nEvents);
  duneobj.rw_nuPDG.resize(duneobj.nEvents);
  duneobj.rw_berpaacvwgt.resize(duneobj.nEvents); 
  duneobj.rw_vtx_x.resize(duneobj.nEvents);
  duneobj.rw_vtx_y.resize(duneobj.nEvents);
  duneobj.rw_vtx_z.resize(duneobj.nEvents);
  duneobj.lepton_momentum.resize(duneobj.nEvents);

  duneobj.global_bin_number.resize(duneobj.nEvents);

  duneobj.mode.resize(duneobj.nEvents);
  duneobj.Target.resize(duneobj.nEvents);

  duneobj.rw_pt.resize(duneobj.nEvents); //perpendicular momentum
  duneobj.rw_pz.resize(duneobj.nEvents); //parallel momentum

  _data->GetEntry(0);

  bool need_global_bin_numbers = (XVarStr == "global_bin_number");

  

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj.nEvents; ++i) { // Loop through tree
  //for (int i = 0; i < events_for_analysis ; ++i) { // Loop through just enough events...............
    _data->GetEntry(i);
    duneobj.rw_cvnnumu[i] = _cvnnumu;
    duneobj.rw_cvnnue[i] = _cvnnue;
    duneobj.rw_cvnnumu_shifted[i] = _cvnnumu; 
    duneobj.rw_cvnnue_shifted[i] = _cvnnue;

    if (iselike) {
      duneobj.rw_erec[i] = _erec_nue;
      duneobj.rw_erec_shifted[i] = _erec_nue; 
      duneobj.rw_erec_had[i] = _erec_had_nue;
      //duneobj.rw_erec_lep[i] = _erec_lep_nue;
    } else {
      duneobj.rw_erec[i] = _erec; 
      duneobj.rw_erec_shifted[i] = _erec; 
      duneobj.rw_erec_had[i] = _erec_had; 
      //duneobj.rw_erec_lep[i] = _erec_lep; 
    }
    duneobj.rw_erec_lep[i]= _Elep_reco;
    duneobj.true_q0[i] = _ev - _LepE;
    duneobj.true_q3[i] = (TVector3{_NuMomX, _NuMomY, _NuMomZ} -
                          TVector3{_LepMomX, _LepMomY, _LepMomZ})
                             .Mag();

    //duneobj.rw_pt[i] =(TVector3{_LepMomX, _LepMomY, _LepMomZ}).Dot((TVector3{TVector3{_NuMomX, _NuMomY, _NuMomZ}}.Divide(TVector3{TVector3{_NuMomX, _NuMomY, _NuMomZ}}.Mag())));
    //duneobj.rw_pz[i] =(TVector3{_LepMomX, _LepMomY, _LepMomZ}).Cross((TVector3{TVector3{_NuMomX, _NuMomY, _NuMomZ}}.Divide(TVector3{TVector3{_NuMomX, _NuMomY, _NuMomZ}}.Mag())));
    TVector3 nuMom(_NuMomX, _NuMomY, _NuMomZ);
    TVector3 nuMomNorm = nuMom.Unit(); // Normalized vector

    duneobj.rw_pt[i] = TVector3(_LepMomX, _LepMomY, _LepMomZ).Dot(nuMomNorm);
    duneobj.rw_pz[i] = (TVector3(_LepMomX, _LepMomY, _LepMomZ).Cross(nuMomNorm)).Mag();


    duneobj.lepton_momentum[i] =(TVector3{_LepMomX, _LepMomY, _LepMomZ}).Mag();
    duneobj.rw_eRecoP[i] = _eRecoP; 
    duneobj.rw_eRecoPip[i] = _eRecoPip; 
    duneobj.rw_eRecoPim[i] = _eRecoPim; 
    duneobj.rw_eRecoPi0[i] = _eRecoPi0; 
    duneobj.rw_eRecoN[i] = _eRecoN; 
    
    duneobj.rw_LepE[i] = _LepE; 
    duneobj.rw_eP[i] = _eP; 
    duneobj.rw_ePip[i] = _ePip; 
    duneobj.rw_ePim[i] = _ePim; 
    duneobj.rw_ePi0[i] = _ePi0; 
    duneobj.rw_eN[i] = _eN; 
    
    duneobj.rw_etru[i] = _ev;
    duneobj.rw_theta[i] = _LepNuAngle;
    duneobj.rw_isCC[i] = _isCC;
    duneobj.rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj.rw_nuPDG[i] = _nuPDG;
    duneobj.rw_berpaacvwgt[i] = _BeRPA_cvwgt;
    duneobj.rw_vtx_x[i] = _vtx_x;
    duneobj.rw_vtx_y[i] = _vtx_y;
    duneobj.rw_vtx_z[i] = _vtx_z;

    if(need_global_bin_numbers){
      duneobj.global_bin_number[i] = GetGenericBinningGlobalBinNumber(iSample, i);
    } 
    //Assume everything is on Argon for now....
    duneobj.Target[i] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj.mode[i]=SIMBMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj.flux_w[i] = 1.0;

    

    //std::cout<< "rw_erec_lep =  " << _erec_lep[i] << std::endl;
  }
  
  _sampleFile->Close();
  return duneobj.nEvents;
}

TH1D* samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill!=-1) {
    if (kChannelToFill>dunemcSamples.size()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}",dunemcSamples.size());
      MACH3LOG_ERROR("kChannelToFill given: {}",kChannelToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
    if (kModeToFill>kMaCh3_nModes) {
      MACH3LOG_ERROR("Required mode is not available. kModeToFill should be between 0 and {}",kMaCh3_nModes);
      MACH3LOG_ERROR("kModeToFill given: {}",kModeToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector< std::vector<double> > SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill+1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill+1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(Var1,SelectionVec,WeightStyle,Axis);
}

/*! DB New version of get1DVarHist which only fills histogram with events passing IsEventSelected
 * This works by having the Selection vector, where each component of Selection is a 2 or 3 length vector
 * If Selection[i].size()==3, Selection[i][0] is the ND280KinematicType which is being cut, and only events with ND280KinematicType values between Selection[i][1] and Selection[i][2] are accepted
 */
TH1D* samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1,std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) {

  Selection = SelectionVec;

  for (unsigned int iStoredSelection=0;iStoredSelection<StoredSelection.size();iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }

  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection].size()!=3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect size == 3, given: {}",iSelection,Selection[iSelection].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  //DB Cut on OscChannel in this function due to speed increase from considering duneSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection=0;iSelection<Selection.size();iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (fChannel && kChannelToFill>dunemcSamples.size()) {
    MACH3LOG_ERROR("Required channel is not available. kChannelToFill should be between 0 and {}",dunemcSamples.size());
    MACH3LOG_ERROR("kChannelToFill given: {}",kChannelToFill);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TH1D* _h1DVar;
  std::vector<double> xBinEdges = ReturnKinematicParameterBinning(Var1);
  _h1DVar = new TH1D("", "", xBinEdges.size()-1, xBinEdges.data());

  //This should be the same as FillArray in core basically, except that
  //events will end up in different bins
  for (int i=0;i<dunemcSamples.size();i++) {
    if (fChannel && (i!=kChannelToFill)) {
      continue;
    }
    for(int j=0;j<dunemcSamples[i].nEvents;j++) {

      //DB Determine which events pass selection
      if (!IsEventSelected(i,j)) {
		continue;
      }

      double Weight = GetEventWeight(i,j);
	  if (WeightStyle==1) {
	    Weight = *(MCSamples[i].osc_w_pointer[j]) * dunemcSamples[i].pot_s * dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[j];
	  }

	  //ETA - not sure about this
	  if (MCSamples[i].xsec_w[j] == 0.) continue;

	  double Var1_Val;

	  Var1_Val = ReturnKinematicParameter(Var1,i,j);
	  if (Var1_Val!=_DEFAULT_RETURN_VAL_) {
		_h1DVar->Fill(Var1_Val,Weight);
	  }
    }
  }

  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects the overall integral of the histogram so it should not be used anywhere we care about event rates
  // We could use Hist->Integral("width") but it would require a lot of modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

  return _h1DVar;
}

double const& samplePDFDUNEBeamFD::ReturnKinematicParameterByReference(int KinematicParameter, int iSample, int iEvent) {

  switch(KinematicParameter){
  case kTrueNeutrinoEnergy:
    return dunemcSamples[iSample].rw_etru[iEvent]; 
  case kRecoNeutrinoEnergy:
    return dunemcSamples[iSample].rw_erec_shifted[iEvent];
  case kTrueXPos:
    return dunemcSamples[iSample].rw_vtx_x[iEvent];
  case kTrueYPos:
    return dunemcSamples[iSample].rw_vtx_y[iEvent];
  case kTrueZPos:
    return dunemcSamples[iSample].rw_vtx_z[iEvent];
  case kCVNNumu:
    return dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
  case kCVNNue:
    return dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
  case kGlobalBinNumber:
    return dunemcSamples[iSample].global_bin_number[iEvent];
  case kELepRec: {
    //std::cout << "ELepRec =  " << dunemcSamples[iSample].rw_erec_lep[iEvent] << std::endl;
    return dunemcSamples[iSample].rw_LepE[iEvent];
  }
  case kq0:
    return dunemcSamples[iSample].true_q0[iEvent];
  case kq3:
    return dunemcSamples[iSample].true_q3[iEvent];
  case k_pT:
    return dunemcSamples[iSample].rw_pt[iEvent];
  case k_pz:
    return dunemcSamples[iSample].rw_pz[iEvent];
  default:
    std::stringstream ss;
    ss << "[ERROR]: " << __FILE__ << ":" << __LINE__
              << " ReturnKinematicParameterByReference Did not recognise "
                 "Kinematic Parameter type:"
              << KinematicParameter
              << ". Is it possibly only available via ReturnKinematicParameter "
                 "(not by reference?) if you need it here, give it storage in "
                 "dunemc_base and move it to here."
              ;
    throw std::runtime_error(ss.str());
  }
  }

  double samplePDFDUNEBeamFD::ReturnKinematicParameter(int KinematicParameter,
                                                       int iSample,
                                                       int iEvent) {
    switch (KinematicParameter) {
    case kERecQE: {
      constexpr double V = 0;        // 0 binding energy for now
      constexpr double mn = 939.565; // neutron mass
      constexpr double mp = 938.272; // proton mass

      double mN_eff = mn - V;
      double mN_oth = mp;

      if (dunemcSamples[iSample].rw_nuPDGunosc[iEvent] <
          0) { // if anti-neutrino, swap target/out masses
        mN_eff = mp - V;
        mN_oth = mn;
      }

      double el = dunemcSamples[iSample].rw_erec_lep[iEvent];

      // this is funky, but don't be scared, it defines an annonymous function
      // in place that grabs the lepton mass in MeV when given the neutrino PDG
      // and whether the interaction was CC or NC and then immediately calls it.
      // It's basically a generalisation of the ternary operator.
      double ml =
          [](int nupdg, bool isCC) {
            switch (std::abs(nupdg)) {
            case 12: {
              return isCC ? 0.511 : 0;
            }
            case 14: {
              return isCC ? 105.66 : 0;
            }
            case 16: {
              return isCC ? 1777.0 : 0;
            }
            }
          }(dunemcSamples[iSample].rw_nuPDGunosc[iEvent],
            dunemcSamples[iSample].rw_isCC[iEvent]);

      double pl = std::sqrt(el*el - ml*ml); // momentum of lepton

      double rEnu =
          (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
          (2 * (mN_eff - el +
                pl * std::cos(dunemcSamples[iSample].rw_theta[iEvent])));

      return rEnu;
    }
    case ktheta_lep: {  // EHadRec - EHadTrue)
      //std::cout << "theta_lep =  " << dunemcSamples[iSample].rw_theta[iEvent] <<std::endl;

      return  dunemcSamples[iSample].rw_theta[iEvent];
    }
    case kp_lep: {  // EHadRec - EHadTrue)

      //std::cout << "lepton momentum = " << dunemcSamples[iSample].lepton_momentum[iEvent] << std::endl;
      return  dunemcSamples[iSample].lepton_momentum[iEvent];
  
    }
    case kEHadRec: {

      return dunemcSamples[iSample].rw_eRecoP[iEvent] +
             dunemcSamples[iSample].rw_eRecoPip[iEvent] +
             dunemcSamples[iSample].rw_eRecoPim[iEvent] +
             dunemcSamples[iSample].rw_eRecoPi0[iEvent] +
             dunemcSamples[iSample].rw_eRecoN[iEvent];
    }
    case kERec_minus_Etrue: {  //(ERec - ETrue, EHadRec - EHadTrue)
         //std::cout << "rw_etru[iEvent] - rw_erec_shifted[iEvent] =  " << dunemcSamples[iSample].rw_etru[iEvent] - dunemcSamples[iSample].rw_erec_shifted[iEvent] << std::endl;
      return dunemcSamples[iSample].rw_etru[iEvent] - dunemcSamples[iSample].rw_erec_shifted[iEvent];
    }
    case kEHadRec_minus_EHadtrue: {  // EHadRec - EHadTrue)

      return ( dunemcSamples[iSample].rw_eRecoP[iEvent] +
             dunemcSamples[iSample].rw_eRecoPip[iEvent] +
             dunemcSamples[iSample].rw_eRecoPim[iEvent] +
             dunemcSamples[iSample].rw_eRecoPi0[iEvent] +
             dunemcSamples[iSample].rw_eRecoN[iEvent] ) - dunemcSamples[iSample].true_q0[iEvent];
    }
    
    default: {
      return ReturnKinematicParameterByReference(KinematicParameter, iSample,
                                                 iEvent);
    }
    }
  }

int samplePDFDUNEBeamFD::ReturnKinematicParameterFromString(std::string KinematicParameterStr){
  if (KinematicParameterStr.find("TrueNeutrinoEnergy") != std::string::npos) {return kTrueNeutrinoEnergy;}
  if (KinematicParameterStr.find("RecoNeutrinoEnergy") != std::string::npos) {return kRecoNeutrinoEnergy;}
  if (KinematicParameterStr.find("TrueXPos") != std::string::npos) {return kTrueXPos;}
  if (KinematicParameterStr.find("TrueYPos") != std::string::npos) {return kTrueYPos;}
  if (KinematicParameterStr.find("TrueZPos") != std::string::npos) {return kTrueZPos;}
  if (KinematicParameterStr.find("CVNNumu") != std::string::npos) {return kCVNNumu;}
  if (KinematicParameterStr.find("CVNNue") != std::string::npos) {return kCVNNue;}
  if (KinematicParameterStr.find("M3Mode") != std::string::npos) {return kM3Mode;}
  if (KinematicParameterStr.find("global_bin_number") != std::string::npos) {return kGlobalBinNumber;}
  if (KinematicParameterStr.find("q0") != std::string::npos) {return kq0;}
  if (KinematicParameterStr.find("q3") != std::string::npos) {return kq3;}
  if (KinematicParameterStr.find("ERecQE") != std::string::npos) {return kERecQE;}
  if (KinematicParameterStr.find("ELepRec") != std::string::npos) {return kELepRec;}
  if (KinematicParameterStr.find("EHadRec") != std::string::npos) {return kEHadRec;}
  if (KinematicParameterStr.find("p_lep") != std::string::npos) {return kp_lep;}
  if (KinematicParameterStr.find("pT") != std::string::npos) {return k_pT;}
  if (KinematicParameterStr.find("pz") != std::string::npos) {return k_pz;}
  if (KinematicParameterStr.find("ERec_minus_Etrue") != std::string::npos) {return kERec_minus_Etrue;}
  //if (KinematicParameterStr.find("ERec_minus_Etrue")!= std::string::npos) {return kERec_minus_Etrue;}
  if (KinematicParameterStr.find("kEHadRec_minus_EHadtrue")!= std::string::npos) {return kEHadRec_minus_EHadtrue;}
  if (KinematicParameterStr.find("theta_lep")!= std::string::npos) {return ktheta_lep;}
 
  std::stringstream ss;
  ss << "[ERROR]: " << __FILE__ << ":" << __LINE__
     << "failed to translate kinematic parameter string "
     << KinematicParameterStr << " to parameter id.";
  throw std::runtime_error(ss.str());
}

std::string samplePDFDUNEBeamFD::ReturnStringFromKinematicParameter(
    int KinematicParameter) {

  switch (KinematicParameter) {
  case kRecoNeutrinoEnergy:
    return "RecoNeutrinoEnergy";
  case kTrueNeutrinoEnergy:
    return "RecoNeutrinoEnergy";
  case kTrueXPos:
    return "TrueXPos";
  case kTrueYPos:
    return "TrueYPos";
  case kTrueZPos:
    return "TrueZPos";
  case kCVNNumu:
    return "CVNNumu";
  case kCVNNue:
    return "CVNNue";
  case kM3Mode:
    return "M3Mode";
  case kGlobalBinNumber:
    return "global_bin_number";
  case kq0:
    return "q0";
  case kq3:
    return "q3";
  case kERecQE:
    return "ERecQE";
  case kELepRec:
    return "ELepRec";
  case kEHadRec:
    return "EHadRec";
  case ktheta_lep:
    return "theta_lep";
  case kp_lep:
    return "p_lep";
  case k_pz:
    return "pz";
  case k_pT:
    return "pT";
  case kERec_minus_Etrue:
    return "ERec_minus_Etrue";
  default: {
   MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
   throw MaCh3Exception(__FILE__, __LINE__);
  }
  }
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  auto &duneobj = dunemcSamples[iSample];
  fdmc_base *fdobj = &(MCSamples[iSample]);  
  
  fdobj->nutype = duneobj.nutype;
  fdobj->oscnutype = duneobj.oscnutype;
  fdobj->signal = duneobj.signal;
  fdobj->SampleDetID = SampleDetID;
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj.rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj.mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj.Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj.rw_isCC[iEvent]);
  }
}
 
void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {
   
  //ETA - this is pretty horrific... we need to think of a nicer way to do this.
  //Don't want to add in hard checks on which systematics are defined but also don't want to hard-code
  //the order in which the systematics are specified. All of these functions should have access to the 
  //dunemc struct so they only need to have iSample and iEvent passed to them. Can probably loop over
  //a vector of std::function objects and pass each of them iSample and iEvent.
  /*
   // reset erec back to original value
  dunemcSamples[iSample].rw_erec_shifted[iEvent] = dunemcSamples[iSample].rw_erec[iEvent];

  // reset cvnnumu back to original value
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnumu[iEvent];

  // reset cvnnue back to original value
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] = dunemcSamples[iSample].rw_cvnnue[iEvent];

  //Calculate values needed
  double sqrtErecHad =  sqrt(dunemcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunemcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunemcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunemcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunemcSamples[iSample].rw_eRecoP[iEvent] + dunemcSamples[iSample].rw_eRecoPip[iEvent] + dunemcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14) && dunemcSamples[iSample].nutype==2};
  bool CCnue {dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==12) && dunemcSamples[iSample].nutype==1};
  bool NotCCnumu {!(dunemcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunemcSamples[iSample].rw_nuPDG[iEvent]==14)) && dunemcSamples[iSample].nutype==2};


  TotalEScaleFD(FDDetectorSystPointers[0], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);

  TotalEScaleSqrtFD(FDDetectorSystPointers[1], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad, sqrtErecLep, NotCCnumu);

  TotalEScaleInvSqrtFD(FDDetectorSystPointers[2], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_had[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecHad, invSqrtErecLep, NotCCnumu);

  HadEScaleFD(FDDetectorSystPointers[3], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

  HadEScaleSqrtFD(FDDetectorSystPointers[4], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, sqrtSumEhad);

  HadEScaleInvSqrtFD(FDDetectorSystPointers[5], &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, invSqrtSumEhad);

  MuEScaleFD(FDDetectorSystPointers[6], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

  MuEScaleSqrtFD(FDDetectorSystPointers[7], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, CCnumu);

  MuEScaleInvSqrtFD(FDDetectorSystPointers[8], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, CCnumu);

  NEScaleFD(FDDetectorSystPointers[9], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent]);

  NEScaleSqrtFD(FDDetectorSystPointers[10], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

  NEScaleInvSqrtFD(FDDetectorSystPointers[11], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

  EMEScaleFD(FDDetectorSystPointers[12], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], CCnue);

  EMEScaleSqrtFD(FDDetectorSystPointers[13], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, sqrteRecoPi0, CCnue);

  EMEScaleInvSqrtFD(FDDetectorSystPointers[14], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, invSqrteRecoPi0, CCnue);

  HadResFD(FDDetectorSystPointers[15], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoP[iEvent], dunemcSamples[iSample].rw_eRecoPip[iEvent], dunemcSamples[iSample].rw_eRecoPim[iEvent], dunemcSamples[iSample].rw_eP[iEvent], dunemcSamples[iSample].rw_ePip[iEvent], dunemcSamples[iSample].rw_ePim[iEvent]);

  MuResFD(FDDetectorSystPointers[16], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnumu);

  NResFD(FDDetectorSystPointers[17], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoN[iEvent], dunemcSamples[iSample].rw_eN[iEvent]);

  EMResFD(FDDetectorSystPointers[18], &dunemcSamples[iSample].rw_erec_shifted[iEvent], dunemcSamples[iSample].rw_eRecoPi0[iEvent], dunemcSamples[iSample].rw_ePi0[iEvent], dunemcSamples[iSample].rw_erec_lep[iEvent], dunemcSamples[iSample].rw_LepE[iEvent], CCnue);

  CVNNumuFD(FDDetectorSystPointers[19], &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent]);

  CVNNueFD(FDDetectorSystPointers[20], &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent]);
  */
}

std::vector<double> samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(int KinematicParameter) {
  std::vector<double> binningVector;

  int nBins = 0;
  double bin_width = 0;
  switch(KinematicParameter){
	case(kRecoNeutrinoEnergy):
	  nBins = 20; 
	  bin_width = 0.5; //GeV
	  break;
	case(kTrueNeutrinoEnergy):
	  nBins = 20; 
	  bin_width = 0.5; //GeV
	  break;
	default:
	  nBins = 10;
	  bin_width = 1.0;
	  break;
  }

  for(int bin_i = 0 ; bin_i < nBins ; bin_i++){
	binningVector.push_back(bin_i*bin_width);
  }

  return binningVector;
}
