#include <TROOT.h>

#include "TMath.h"
#include "TString.h"
#include "manager/manager.h"
#include "samplePDFDUNEBase.h"
#include <assert.h>
#include <stdexcept>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <filesystem>
#include <iostream>

// #define DEBUG

//Read in all the input CAF files




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

  std::unique_ptr<TH3D> Abis3DHistogram;
  std::unique_ptr<TH3D> OA3DHistogram; 
  std::unique_ptr<TH2D> OA2DHistogram; 

  TH1D *onedim_binnumberhisto;

// Constructors for erec-binned errors

//!!rewrite execs to give arguments in new order
samplePDFDUNEBase::samplePDFDUNEBase(double pot, std::string mc_version,
                                     covarianceXsec *xsec_cov)
    : samplePDFBase(pot) {
  std::cout << "- Using DUNE sample config in this file " << mc_version
            << std::endl;
  // ETA - safety feature so you can't pass a NULL xsec_cov
  if (xsec_cov == NULL) {
    std::cerr << "[ERROR:] You've passed me a NULL xsec covariance matrix... I "
                 "need this to setup splines!"
              << std::endl;
    throw;
  }
  init(pot, mc_version, xsec_cov);


  //Somewhere in samplePDFDUNEBase::samplePDFDUNEBase
  //int number_of_pt_bins = h->GetXaxis()->GetNbins();

    std::vector<double> ELep_bins;
    ExtendLinspace(ELep_bins,0,40,5);
    std::vector<double> theta_bins;
    ExtendLinspace(theta_bins,0,3,10);
    ExtendLinspace(theta_bins,3,9,3);
    ExtendLinspace(theta_bins,9,100,2);
    //ExtendLinspace(theta_bins,0,100,10);
    std::vector<double> ENuReco_bins={0.,   1.,  1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 10.};
    //ExtendLinspace(ENuReco_bins,0,4,8);

    std::vector<double> offaxis_position;
    ExtendLinspace(offaxis_position,-35,0,10);
    
     //TH3D* Abis3DHistogram = new TH3D("Abis3DHistogram", "", ELep_bins.size(), theta_bins.size(),  ENuReco_bins.size()); // fill this in with the binning like we used to in NUISANCE
     std::cout<< "ELep_bins.size() = " << ELep_bins.size() <<std::endl;
     std::cout<< "theta_bins.size() = " << theta_bins.size() <<std::endl;
     std::cout<< "ENuReco_bins.size() = " << ENuReco_bins.size() <<std::endl;
     
     Abis3DHistogram =
        std::make_unique<TH3D>("Abis3DHistogram", "", ELep_bins.size() - 1,
                               ELep_bins.data(), theta_bins.size() - 1, theta_bins.data(),
                               ENuReco_bins.size() - 1, ENuReco_bins.data());

      double noofbins_1Dhisto = Abis3DHistogram->GetNcells();
      std::cout<< "no of bins in 3D HISTO = " << Abis3DHistogram->GetNcells() <<std::endl;

      OA3DHistogram =
        std::make_unique<TH3D>("OA3DHistogram", "", ELep_bins.size() - 1,
                               ELep_bins.data(), theta_bins.size() - 1, theta_bins.data(),
                               offaxis_position.size() - 1, offaxis_position.data());




      OA2DHistogram =
        std::make_unique<TH2D>("OA2DHistogram", "",offaxis_position.size() - 1, offaxis_position.data(),
        ENuReco_bins.size() - 1, ENuReco_bins.data());

      //double noofbins_1Dhisto = Abis3DHistogram->GetNcells();
      //std::cout<< "no of bins in 3D HISTO = " << Abis3DHistogram->GetNcells() <<std::endl;



     
      onedim_binnumberhisto = 
      new TH1D("onedim_binnumberhisto","", noofbins_1Dhisto, -1, noofbins_1Dhisto -1 );

}

samplePDFDUNEBase::~samplePDFDUNEBase() {}

void samplePDFDUNEBase::init(double pot, std::string samplecfgfile,
                             covarianceXsec *xsec_cov) {

  Beta = 1;
  useBeta = false;
  applyBetaNue = false;
  applyBetaDiag = false;

  // doubled_angle =true ;
  useNonDoubledAngles(true);
  if (doubled_angle)
    std::cout << "- Using non doubled angles for oscillation parameters"
              << std::endl;

  osc_binned = false;
  if (osc_binned)
    std::cout << "- Using binned oscillation weights" << std::endl;

  modes = new TH1D("modes", "", 120, -60, 60);

  std::string mtupleprefix;
  std::string mtuplesuffix;
  std::string splineprefix;
  std::string splinesuffix;

  char *sample_char = (char *)samplecfgfile.c_str();
  // ETA - trying to read stuff from yaml file
  manager *SampleManager = new manager(sample_char);

  // Bools
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID =SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  // Inputs
  /*
  mtupleprefix =
      SampleManager->raw()["InputFiles"]["mtupleprefix"].as<std::string>();
  mtuplesuffix =
      SampleManager->raw()["InputFiles"]["mtuplesuffix"].as<std::string>();
  splineprefix =
      SampleManager->raw()["InputFiles"]["splineprefix"].as<std::string>();
  splinesuffix =
      SampleManager->raw()["InputFiles"]["splinesuffix"].as<std::string>();
      */

  // Binning
  BinningOpt = SampleManager->raw()["Binning"]["BinningOpt"].as<int>();
  std::vector<double> sample_erec_bins =
      SampleManager->raw()["Binning"]["XVarBins"].as<std::vector<double>>();
  std::vector<double> sample_theta_bins =
      SampleManager->raw()["Binning"]["YVarBins"].as<std::vector<double>>();

  
  
  double Number_of_uniform_bins = SampleManager->raw()["Binning"]["NTotalBins"].as<double>();
  std::vector<double> uniform_bins;
      for (int i = 0; i < Number_of_uniform_bins; ++i) {
         uniform_bins.push_back(i);
      }
      

  samplename = SampleManager->raw()["SampleName"].as<std::string>();

  std::vector<std::string> mtuple_files;
  std::vector<std::string> spline_files;
  std::vector<int> sample_vecno;
  std::vector<int> sample_oscnutype;
  std::vector<int> sample_nutype;
  std::vector<bool> sample_signal;


  //Loop over all the subsamples using a wildcard so you can loop over the whole directory

   std::string caf_file_directory = "/home/abipeake/Mach3/MaCh3_DUNE/inputs/testcafs/*";
   std::cout<< "looking in caf_file_directory :  " << caf_file_directory << std::endl;
   int number_of_subsamples=0;
        // pretends shell glob caf_file_directorys are regex any char
        size_t next_ast = caf_file_directory.find_first_of('*');
        std::cout<< "files in CAF folder =  :  " << next_ast << std::endl;

        while(next_ast != std::string::npos){
          std::cout<< "in while loop finding subsamples in CAF dir" << std::endl;
                if(!next_ast || (caf_file_directory[next_ast-1] != '.')){
                        caf_file_directory.replace(next_ast,1,".*");
                        next_ast++;
                        //number_of_subsamples = number_of_subsamples +1;
                       
                        //std::cout << "number_of_subsamples " << number_of_subsamples << std::endl;
                }
                next_ast = caf_file_directory.find_first_of('*',next_ast+1);
        }

        size_t last_fs = caf_file_directory.find_last_of('/');
        std::cout<< " found all subsamples in CAF dir" << std::endl;
        std::filesystem::path dir;
        std::regex file_re;
        if(last_fs == std::string::npos){
                dir = "./";
                file_re = std::regex(caf_file_directory);
        } else {
                dir = caf_file_directory.substr(0,last_fs);
                file_re = std::regex(caf_file_directory.substr(last_fs+1));
        }
        if(!std::filesystem::exists(dir)){
                std::cerr << "[ERROR]: caf_file_directory directory: " << dir.native() << " does not exist." << std::endl;
                throw;
        }
        for (auto const &dir_entry :
          std::filesystem::directory_iterator{dir}) {
                if (std::regex_match(dir_entry.path().filename().native(), file_re)) {
                        std::cout << "found matching file:" << dir_entry.path().native() << std::endl;
                        
                        number_of_subsamples = number_of_subsamples + 1 ;
                        mtuple_files.push_back(dir_entry.path().filename().native());
                        spline_files.push_back(dir_entry.path().filename().native());
                       
                        sample_nutype.push_back(14);
                        sample_oscnutype.push_back(14);
                        sample_signal.push_back(false);
                        sample_vecno.push_back(number_of_subsamples);
                        std::cout<< "number_of_subsamples = " << number_of_subsamples <<std::endl;
                }
        }

    /*---------------------original
    for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    std::cout << "Found sub sample" << std::endl;
    
    mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
    spline_files.push_back(osc_channel["splinefile"].as<std::string>());
    sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
    sample_nutype.push_back(
        PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
    sample_oscnutype.push_back(
        PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
    sample_signal.push_back(osc_channel["signal"].as<bool>());
  }*/


  /*
  // Loop over all the sub-samples --original version
  for (auto const &osc_channel : SampleManager->raw()["SubSamples"]) {
    std::cout << "Found sub sample" << std::endl;
    mtuple_files.push_back(osc_channel["mtuplefile"].as<std::string>());
    spline_files.push_back(osc_channel["splinefile"].as<std::string>());
    sample_vecno.push_back(osc_channel["samplevecno"].as<int>());
    sample_nutype.push_back(
        PDGToProbs(static_cast<NuPDG>(osc_channel["nutype"].as<int>())));
    sample_oscnutype.push_back(
        PDGToProbs(static_cast<NuPDG>(osc_channel["oscnutype"].as<int>())));
    sample_signal.push_back(osc_channel["signal"].as<bool>());
  }*/

  
  // Now loop over the kinematic cuts
  for (auto const &SelectionCuts : SampleManager->raw()["SelectionCuts"]) {
    std::cout << "Looping over selection cuts " << std::endl;
    SelectionStr.push_back(SelectionCuts["KinematicStr"].as<std::string>());

    SelectionBounds.push_back(
        SelectionCuts["Bounds"].as<std::vector<double>>());
    std::cout << "Found cut on string "
              << SelectionCuts["KinematicStr"].as<std::string>() << std::endl;
    std::cout << "With bounds "
              << SelectionCuts["Bounds"].as<std::vector<double>>()[0] << " to "
              << SelectionCuts["Bounds"].as<std::vector<double>>()[1]
              << std::endl;
  }
  NSelections = SelectionStr.size();

  // Make some arrays so we can initialise _hPDF1D and _hPDF2D with these
  double erec_bin_edges[sample_erec_bins.size()];
  double theta_bin_edges[sample_theta_bins.size()];
  for (unsigned erec_i = 0; erec_i < sample_erec_bins.size(); erec_i++) {
    erec_bin_edges[erec_i] = sample_erec_bins[erec_i];
  }
  for (unsigned theta_i = 0; theta_i < sample_theta_bins.size(); theta_i++) {
    theta_bin_edges[theta_i] = sample_theta_bins[theta_i];
  }

  // create dunemc storage -original 
  /*int nSamples = SampleManager->raw()["NSubSamples"].as<int>();
  for (int i = 0; i < nSamples; i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }*/
  
  
  for (int i = 0; i < number_of_subsamples; i++) {
    //std::cout << "line 315 number of subsamples = " << number_of_subsamples << std::endl;
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
    }

  // Now down with yaml file for sample - original
  delete SampleManager;
  std::cout << "Number of subsamples line 322 = : " << number_of_subsamples << endl;
  std::cout << "Oscnutype size: " << sample_oscnutype.size()
            << ", dunemcSamples size: " << dunemcSamples.size() << endl;
  if (sample_oscnutype.size() != dunemcSamples.size()) {
    std::cerr << "[ERROR:] samplePDFDUNEBase::samplePDFDUNEBase() - something "
                 "went wrong either getting information from sample config"
              << std::endl;
    throw;
  }
  

  for (unsigned iSample = 0; iSample < dunemcSamples.size(); iSample++) {
    setupDUNEMC((mtupleprefix + mtuple_files[iSample] + mtuplesuffix).c_str(),
                &dunemcSamples[sample_vecno[iSample]], pot,
                sample_nutype[iSample], sample_oscnutype[iSample],
                sample_signal[iSample]);
  }

  for (int i = 0; i < number_of_subsamples; i++) {
    // std::cout << "line 341 number of subsamples = " << number_of_subsamples << std::endl;
    struct fdmc_base obj = fdmc_base();
    MCSamples.push_back(obj);
  }

  for (unsigned iSample = 0; iSample < MCSamples.size(); iSample++) {
    setupFDMC(&dunemcSamples[sample_vecno[iSample]],
              &MCSamples[sample_vecno[iSample]],
              (splineprefix + spline_files[iSample] + splinesuffix).c_str());
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD MC   " << std::endl;
  std::cout << "################" << std::endl;

  // ETA - If xsec_cov hasn't been passed to the samplePDFDUNEBase constructor
  // then it's NULL and the old funcitonality is kept this calls this function
  // in the core code this needs to come after setupFDMC as otherwise
  // MCSamples.splinefile will be NULL
  setXsecCov(xsec_cov);

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

  nFDDetectorSystPointers = funcParsIndex.size();
  FDDetectorSystPointers = std::vector<const double *>(nFDDetectorSystPointers);

  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin();
       it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);

    if (name == "TotalEScaleFD") {
      tot_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(tot_escale_fd_pos);
    }

    else if (name == "TotalEScaleSqrtFD") {
      tot_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(tot_escale_sqrt_fd_pos);
    }

    else if (name == "TotalEScaleInvSqrtFD") {
      tot_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(tot_escale_invsqrt_fd_pos);
    }

    else if (name == "HadEScaleFD") {
      had_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_escale_fd_pos);
    }

    else if (name == "HadEScaleSqrtFD") {
      had_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(had_escale_sqrt_fd_pos);
    }

    else if (name == "HadEScaleInvSqrtFD") {
      had_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(had_escale_invsqrt_fd_pos);
    }

    else if (name == "MuEScaleFD") {
      mu_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_escale_fd_pos);
    }

    else if (name == "MuEScaleSqrtFD") {
      mu_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(mu_escale_sqrt_fd_pos);
    }

    else if (name == "MuEScaleInvSqrtFD") {
      mu_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(mu_escale_invsqrt_fd_pos);
    }

    else if (name == "NEScaleFD") {
      n_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_escale_fd_pos);
    }

    else if (name == "NEScaleSqrtFD") {
      n_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(n_escale_sqrt_fd_pos);
    }

    else if (name == "NEScaleInvSqrtFD") {
      n_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(n_escale_invsqrt_fd_pos);
    }

    else if (name == "EMEScaleFD") {
      em_escale_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_escale_fd_pos);
    }

    else if (name == "EMEScaleSqrtFD") {
      em_escale_sqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(em_escale_sqrt_fd_pos);
    }

    else if (name == "EMEScaleInvSqrtFD") {
      em_escale_invsqrt_fd_pos = *it;
      FDDetectorSystPointers[func_it] =
          xsec_cov->retPointer(em_escale_invsqrt_fd_pos);
    }

    else if (name == "HadResFD") {
      had_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(had_res_fd_pos);
    }

    else if (name == "MuResFD") {
      mu_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(mu_res_fd_pos);
    }

    else if (name == "NResFD") {
      n_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(n_res_fd_pos);
    }

    else if (name == "EMResFD") {
      em_res_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(em_res_fd_pos);
    } else if (name == "CVNNumuFD") {
      cvn_numu_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(cvn_numu_fd_pos);
    } else if (name == "CVNNueFD") {
      cvn_nue_fd_pos = *it;
      FDDetectorSystPointers[func_it] = xsec_cov->retPointer(cvn_nue_fd_pos);
    }

    else {
      std::cerr << "Found a functional parameter which wasn't specified in the "
                   "xml | samplePDFDUNEBase:"
                << name << std::endl;
      throw;
    }
  }

  std::cout << "Now setting up Splines" << std::endl;
  for (unsigned iSample = 0; iSample < MCSamples.size(); iSample++) {
    if(!splineprefix.size() || !splinesuffix.size()){
      continue;
    }
    setupSplines(&MCSamples[sample_vecno[iSample]],
                 (splineprefix + spline_files[iSample] + splinesuffix).c_str(),
                 MCSamples[iSample].nutype, MCSamples[iSample].signal);
  }

  std::cout << "################" << std::endl;
  std::cout << "Setup FD splines   " << std::endl;
  std::cout << "################" << std::endl;

  std::cout<< "about to call weight pointers" <<std::endl;
  setupWeightPointers();
  std::cout<< "finished weight pointers" <<std::endl;

  fillSplineBins();

#ifdef USE_PROB3
  std::cout << "- Setup Prob3++" << std::endl;
#else
  std::cout << "- Setup CUDAProb3" << std::endl;
#endif

  _sampleFile->Close();
  char *histname = (char *)"blah";
  char *histtitle = (char *)"blahblah";

  std::cout
      << "-------------------------------------------------------------------"
      << std::endl;

  // The binning here is arbitrary, now we get info from cfg so the
  // set1DBinning and set2Dbinning calls below will make the binning
  // to be what we actually want
  _hPDF1D = new TH1D("hErec_nue", "Reconstructed Energy", 200, 0, 50.0);
  dathist = new TH1D("dat_nue", "", 200, 0, 50.0);
  _hPDF2D = new TH2D(histname, histtitle, 15, 0, 50.0 * 1000, 15, 0, 150);
  dathist2d = new TH2D("dat2d_nue", "", 15, 0, 1500, 15, 0, 150);

  // ETA Don't forget the -1 on the size here, as it's number of bins not bin
  // edges
  set1DBinning(sample_erec_bins.size() - 1, erec_bin_edges);
  set2DBinning(sample_erec_bins.size() - 1, erec_bin_edges,
               sample_theta_bins.size() - 1, theta_bin_edges);

  return;
}

void samplePDFDUNEBase::setupWeightPointers() {
  
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    //std::cout<< "help in first for loop"<<std::endl;
	for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
    //std::cout<< "help0"<<std::endl;
	  //DB Setting total weight pointers
	  MCSamples[i].ntotal_weight_pointers[j] = 6;
   // std::cout<< "help1"<<std::endl;
	  MCSamples[i].total_weight_pointers[j] = new double*[MCSamples[i].ntotal_weight_pointers[j]];
    //std::cout<< "help2"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][0] = &(dunemcSamples[i].pot_s);
    //std::cout<< "help3"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].norm_s);
    //std::cout<< "help4"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].osc_w[j]);
    //std::cout<< "help5"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][3] = &(dunemcSamples[i].rw_berpaacvwgt[j]);
    //std::cout<< "help6"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][4] = &(MCSamples[i].flux_w[j]);
    //std::cout<< "help7"<<std::endl;
	  MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
	}
  }

  return;
}



void samplePDFDUNEBase::setupDUNEMC(const char *sampleFile,
                                    dunemc_base *duneobj, double pot,
                                    int nutype, int oscnutype, bool signal,
                                    bool hasfloats) {

  // set up splines
  std::cout
      << "-------------------------------------------------------------------"
      << std::endl;
  std::cout << "input file: " << sampleFile << std::endl;

  _sampleFile = new TFile(sampleFile, "READ");
  _data = (TTree *)_sampleFile->Get("cafTree");  //cafTree

  if (_data) {
    std::cout << "Found mtuple tree is " << sampleFile << std::endl;
    std::cout << "N of entries: " << _data->GetEntries() << std::endl;
  }

  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco_numu", 1);
  _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco", 1);
  _data->SetBranchAddress("Ev_reco", &_Ev_reco);//////for oa sampls
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  _data->SetBranchStatus("RecoLepEnNumu", 1);
  _data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
  _data->SetBranchStatus("RecoLepEnNue", 1);
  _data->SetBranchAddress("RecoLepEnNue", &_erec_lep_nue);
  _data->SetBranchStatus("RecoLepAngNumu",1);
  _data->SetBranchAddress("RecoLepAngNumu", &_erec_lep_ang_numu);

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

  _data->SetBranchStatus("mode", 1);
  _data->SetBranchAddress("mode", &_mode);
  _data->SetBranchStatus("cvnnumu", 1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue", 1);
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
  _data->SetBranchStatus("det_x", 1);
  _data->SetBranchAddress("det_x", &_det_x);
  


  TH1D *norm = (TH1D *)_sampleFile->Get("norm");
  if (!norm) {
    //std::cout << "Add a norm KEY to the root file using MakeNormHists.cxx"
             // << std::endl;
    //std::cout << "Ignoring for now" << std::endl;
    //std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    //throw;

  norm = new TH1D("norm","",1,0,1);
  norm->SetBinContent(1,1);
  }

  //pot scaling  duneobj->norm_s = (1e21/1.905e21);
  // now fill the actual variables
  duneobj->norm_s = 1.0; //norm->GetBinContent(1);
  //duneobj->pot_s = (pot) / norm->GetBinContent(1);
  duneobj->pot_s = (pot) / 3.85e21;
  //duneobj->pot_s = (pot) / norm->GetBinContent(2);
  std::cout << "pot = " << (pot)<< std::endl;
  std::cout << "pot_s = " << duneobj->pot_s << std::endl;
  std::cout << "norm_s = " << duneobj->norm_s << std::endl;
  duneobj->nEvents = _data->GetEntries();
  duneobj->nutype = nutype;
  duneobj->oscnutype = oscnutype;
  duneobj->signal = signal;

  std::cout << "signal: " << duneobj->signal << std::endl;
  std::cout << "nevents: " << duneobj->nEvents << std::endl;
  std::cout << "nutype " << duneobj->nutype << std::endl;

  // allocate memory for dunemc variables
  duneobj->rw_cvnnumu = new double[duneobj->nEvents];
  duneobj->rw_cvnnue = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu_shifted = new double[duneobj->nEvents];
  duneobj->rw_cvnnue_shifted = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];

  duneobj->rw_erec_had_nd = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];
  
  duneobj->enureco = new double[duneobj->nEvents];
  

  duneobj->rw_eRecoP = new double[duneobj->nEvents];
  duneobj->rw_eRecoPip = new double[duneobj->nEvents];
  duneobj->rw_eRecoPim = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0 = new double[duneobj->nEvents];
  duneobj->rw_eRecoN = new double[duneobj->nEvents];

  duneobj->rw_LepE = new double[duneobj->nEvents];
  duneobj->rw_eP = new double[duneobj->nEvents];
  duneobj->rw_ePip = new double[duneobj->nEvents];
  duneobj->rw_ePim = new double[duneobj->nEvents];
  duneobj->rw_ePi0 = new double[duneobj->nEvents];
  duneobj->rw_eN = new double[duneobj->nEvents];

  duneobj->rw_lep_ang_numu = new double[duneobj->nEvents];

  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->xsec_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents];
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];

  duneobj->rw_abis3dbinnumber = new double[duneobj->nEvents]; /////////////////////new bin number
  duneobj->detector_oa_position = new double[duneobj->nEvents];
  duneobj->offaxis_3dbinnumber = new double[duneobj->nEvents];
  duneobj->offaxis_2dbinnumber = new double[duneobj->nEvents];
  duneobj->uniform_bins = new double[duneobj->nEvents];
  
  //duneobk->global_binnumber = new double[duneobj->nEvents];

  duneobj->energyscale_w = new double[duneobj->nEvents];
  duneobj->mode = new int[duneobj->nEvents];
  duneobj->rw_lower_erec_1d =
      new double[duneobj->nEvents]; // lower erec bound for bin
  duneobj->rw_upper_erec_1d =
      new double[duneobj->nEvents]; // upper erec bound for bin
  duneobj->rw_lower_erec_2d =
      new double[duneobj->nEvents]; // lower erec bound for bin
  duneobj->rw_upper_erec_2d =
      new double[duneobj->nEvents]; // upper erec bound for bin

  // These spline bins get filled in fillSplineBins
  duneobj->enu_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->erec_s_bin = new unsigned int[duneobj->nEvents];
  duneobj->flux_bin = new int[duneobj->nEvents];
  duneobj->xsec_norms_bins = new std::list<int>[duneobj->nEvents];
  // change so only points to one
  duneobj->Target = new int[duneobj->nEvents];

  _data->GetEntry(0);

  // FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) // Loop through tree
  {
    _data->GetEntry(i);
    duneobj->rw_cvnnumu[i] = (double)_cvnnumu;
    duneobj->rw_cvnnue[i] = (double)_cvnnue;
    duneobj->rw_cvnnumu_shifted[i] = (double)_cvnnumu;
    duneobj->rw_cvnnue_shifted[i] = (double)_cvnnue;
    
    duneobj->enureco[i] = (double)_Ev_reco;
  
    if (iselike) {
      duneobj->rw_erec[i] = (double)_erec_nue;
      duneobj->rw_erec_shifted[i] = (double)_erec_nue;
      duneobj->rw_erec_had[i] = (double)_erec_had_nue;
      duneobj->rw_erec_lep[i] = (double)_erec_lep_nue;
    } else {
      duneobj->rw_erec[i] = (double)_erec;
      duneobj->rw_erec_shifted[i] = (double)_erec;
      duneobj->rw_erec_had[i] = (double)_erec_had;
      duneobj->rw_erec_lep[i] = (double)_erec_lep;
    }
    duneobj->rw_erec_had_nd[i] = (double)(_erec - _erec_lep);
    duneobj->detector_oa_position[i] = (((double)_det_x + (double)_vtx_x)/100.0);  ///make OA position positive and in m   +vtx  
    duneobj->rw_lep_ang_numu[i] = (double)_erec_lep_ang_numu;
    duneobj->rw_eRecoP[i] = (double)_eRecoP;
    duneobj->rw_eRecoPip[i] = (double)_eRecoPip;
    duneobj->rw_eRecoPim[i] = (double)_eRecoPim;
    duneobj->rw_eRecoPi0[i] = (double)_eRecoPi0;
    duneobj->rw_eRecoN[i] = (double)_eRecoN;

    duneobj->rw_LepE[i] = (double)_LepE;
    duneobj->rw_eP[i] = (double)_eP;
    duneobj->rw_ePip[i] = (double)_ePip;
    duneobj->rw_ePim[i] = (double)_ePim;
    duneobj->rw_ePi0[i] = (double)_ePi0;
    duneobj->rw_eN[i] = (double)_eN;

    duneobj->rw_etru[i] = (double)_ev;
    duneobj->rw_theta[i] = (double)_LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (double)_BeRPA_cvwgt;
    duneobj->rw_vtx_x[i] = (double)_vtx_x;
    duneobj->rw_vtx_y[i] = (double)_vtx_y;
    duneobj->rw_vtx_z[i] = (double)_vtx_z;

    // Assume everything is on Argon for now....
    duneobj->Target[i] = 40;

    duneobj->xsec_w[i] = 1.0;

    // fill modes
    modes->Fill(_mode);
    //!!possible cc1pi exception might need to be 11
    int mode = TMath::Abs(_mode);

    duneobj->mode[i] = SIMBMode_ToMaCh3Mode(mode, _isCC);

    duneobj->energyscale_w[i] = 1.0;

    duneobj->flux_w[i] = 1.0;

    duneobj->uniform_bins[i] = 1.0;

   //detector_oa_position[i] = "<< (-1.0*(double)_det_x)/100.0<<std::endl;
   //std::cout<< "rw_etrue = " << _ev << std::endl;
  }

  _sampleFile->Close();
  std::cout << "Sample set up OK" << std::endl;
}

double
samplePDFDUNEBase::ReturnKinematicParameter(std::string KinematicParameter,
                                            int iSample, int iEvent) {

  KinematicTypes KinPar = static_cast<KinematicTypes>(
      ReturnKinematicParameterFromString(KinematicParameter));
  double KinematicValue = -999;

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    KinematicValue = dunemcSamples[iSample].rw_etru[iEvent];
    break;
  case kTrueXPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_x[iEvent];
    break;
  case kTrueYPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_y[iEvent];
    break;
  case kTrueZPos:
    KinematicValue = dunemcSamples[iSample].rw_vtx_z[iEvent];
    break;
  case kCVNNumu:
    KinematicValue = dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
    break;
  case kCVNNue:
    KinematicValue = dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
    break;

  case kInvariantHadronicMass_W:{
       // Get W_true with assumption of initial state nucleon at rest
      // m_n =  ; //GeV//(float)PhysConst::mass_proton;
      // Q2 assuming nucleon at res
      //double W_nuc_rest = sqrt(-Q2 + 2 * m_n * q0 + m_n * m_n);
      double q0 = dunemcSamples[iSample].rw_erec_lep[iEvent]-dunemcSamples[iSample].rw_cvnnue[iEvent];
	   KinematicValue = sqrt(-dunemcSamples[iSample].rw_Q2[iEvent]+ (2* 0.938 * q0) + (0.938*0.938) );
     //std::cout<< "W = " << KinematicValue << std::endl;
     //std::cout<< "rw_Q2[iEvent] = " << _Q2 << std::endl;
	   break;
    }
    case kBinningin2D:{
        double ELep = dunemcSamples[iSample].rw_etru[iEvent];  // calculate ELep
        double offaxis_position = dunemcSamples[iSample].detector_oa_position[iEvent]; // OA position
        KinematicValue = OA2DHistogram->FindFixBin(offaxis_position,ELep);
        double global_bin_prism_2D = OA2DHistogram->FindFixBin(offaxis_position,ELep);
        //std::cout << "Kinematic Value for 2D histogram = " << KinematicValue <<std::endl;
        dunemcSamples[iSample].offaxis_2dbinnumber[iEvent] = global_bin_prism_2D;
        //std::cout<<"offaxis_2dbinnumber  = " << global_bin_prism_2D << std::endl;
        //std::cout<<"OA positon  = " << offaxis_position << std::endl;
        //std::cout<<" rw_etru  = " << ELep << std::endl;
        break; 
      }
  case kAbis3DBinning:{
        double ELep = dunemcSamples[iSample].rw_LepE[iEvent];  // calculate ELep
        double thetaLep = dunemcSamples[iSample].rw_lep_ang_numu[iEvent]; // calculate thetaLep 
        double ENuReco = dunemcSamples[iSample].rw_erec[iEvent]; // calculate Reconstructed neutrino energy
        KinematicValue = Abis3DHistogram->FindFixBin(ELep,thetaLep,ENuReco);

        /*
        std::cout << "Kinematic Value for Abis3DHistogram = " << KinematicValue <<std::endl;
        std::cout << "Elep = " << ELep <<std::endl;
        std::cout << "thetalep = " << thetaLep <<std::endl;
        std::cout << "ENuReco= " << ENuReco <<std::endl;
        */
        double global_bin = Abis3DHistogram->FindFixBin(ELep,thetaLep,ENuReco);
        //std::cout << "Kinematic Value for Abis3DHistogram = " << KinematicValue <<std::endl;
        dunemcSamples[iSample].rw_abis3dbinnumber[iEvent] = global_bin;
        onedim_binnumberhisto->Fill(global_bin);
        
        break;
      }
    

    case kDUNEPRISMBinning:{
        double ELep = dunemcSamples[iSample].rw_LepE[iEvent];  // calculate ELep
        double Ehad = dunemcSamples[iSample].rw_erec_had[iEvent]; // calculate thetaLep 
        double offaxis_position = dunemcSamples[iSample].detector_oa_position[iEvent]; // OA position
        KinematicValue = OA3DHistogram->FindFixBin(ELep,Ehad,offaxis_position);
        double global_bin_prism = OA3DHistogram->FindFixBin(ELep,Ehad,offaxis_position);
        std::cout << "Kinematic Value for kDUNEPRISMBinning = " << KinematicValue <<std::endl;
        dunemcSamples[iSample].offaxis_3dbinnumber[iEvent] = global_bin_prism;
        break;
      }

      
      //onedim_binnumberhisto->Write();
  case kM3Mode:
        KinematicValue = dunemcSamples[iSample].mode[iEvent];
        break;
  default:
    std::cout << "[ERROR]: " << __FILE__ << ":" << __LINE__
              << " Did not recognise Kinematic Parameter type..." << std::endl;
    throw;
  }

  return KinematicValue;
}


//Fill a 1D histogram with the kinematic values of Abis 3D histogram

void samplePDFDUNEBase::setupFDMC(dunemc_base *duneobj, fdmc_base *fdobj,
                                  const char *splineFile) {

  fdobj->nEvents = duneobj->nEvents;
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  fdobj->x_var = new double *[fdobj->nEvents];
  fdobj->y_var = new double *[fdobj->nEvents];
  fdobj->enu_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->xvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->yvar_s_bin = new unsigned int[fdobj->nEvents];
  fdobj->rw_etru = new double *[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];
  fdobj->NomXBin = new int[fdobj->nEvents];
  fdobj->NomYBin = new int[fdobj->nEvents];
  fdobj->XBin = new int[fdobj->nEvents];
  fdobj->YBin = new int[fdobj->nEvents];

  
  
  fdobj->rw_lower_xbinedge = new double[fdobj->nEvents];
  fdobj->rw_lower_lower_xbinedge = new double[fdobj->nEvents];
  fdobj->rw_upper_xbinedge = new double[fdobj->nEvents];
  fdobj->rw_upper_upper_xbinedge = new double[fdobj->nEvents];
  fdobj->mode = new int *[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double **[fdobj->nEvents];
  fdobj->nxsec_norm_pointers = new int[fdobj->nEvents];
  fdobj->xsec_norm_pointers = new const double **[fdobj->nEvents];
  fdobj->xsec_norms_bins = new std::list<int>[fdobj->nEvents];
  fdobj->xsec_w = new double[fdobj->nEvents];
  fdobj->flux_w = new double[fdobj->nEvents];
  fdobj->osc_w = new double[fdobj->nEvents];
  fdobj->isNC = new bool[fdobj->nEvents];
  fdobj->nxsec_spline_pointers = new int[fdobj->nEvents];
  fdobj->xsec_spline_pointers = new const double **[fdobj->nEvents];
  fdobj->ntotal_weight_pointers = new int[fdobj->nEvents];
  fdobj->total_weight_pointers = new double **[fdobj->nEvents];
  fdobj->Target = new int *[fdobj->nEvents];

  for (int iEvent = 0; iEvent < fdobj->nEvents; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
    // ETA - set these to a dummy value to begin with, these get set in
    // set1DBinning or set2DBinning
    fdobj->NomXBin[iEvent] = -1;
    fdobj->NomYBin[iEvent] = -1;
    fdobj->XBin[iEvent] = -1;
    fdobj->YBin[iEvent] = -1;
    fdobj->rw_lower_xbinedge[iEvent] = -1;
    fdobj->rw_lower_lower_xbinedge[iEvent] = -1;
    fdobj->rw_upper_xbinedge[iEvent] = -1;
    fdobj->rw_upper_upper_xbinedge[iEvent] = -1;
    fdobj->xsec_w[iEvent] = 1.0;
    fdobj->osc_w[iEvent] = 1.0;
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->flux_w[iEvent] = duneobj->flux_w[iEvent];

    // ETA - this is where the variables that you want to bin your samples in
    // are defined If you want to bin in different variables this is where you
    // put it for now
    switch (BinningOpt) {
    case 0:
    case 1:
      // Just point to xvar to the address of the variable you want to bin in
      // This way we don't have to update both fdmc and skmc when we apply
      // shifts to variables we're binning in
      fdobj->x_var[iEvent] = &(duneobj->offaxis_2dbinnumber[iEvent]); //offaxis_3dbinnumber
      fdobj->y_var[iEvent] =
          &(duneobj
                ->dummy_y); // ETA - don't think we even need this as if we have
                            // a 1D sample we never need this, just not sure I
                            // like an unitialised variable in fdmc struct?
      break;
    case 2:
      // Just point to xvar to the address of the variable you want to bin in
      // This way we don't have to update both fdmc and skmc when we apply
      // shifts to variables we're binning in
      fdobj->y_var[iEvent] = &(duneobj->detector_oa_position[iEvent]);
      fdobj->x_var[iEvent] = &(duneobj->rw_etru[iEvent]);

      
      break;
    case 3: // binning opt 3 sets itself up different but just uses the 1d binning projection
      //fdobj->x_var[iEvent] = doc["Binning"]["BinningOpt"].as<int>();
      //vector<double> uniform_binning(100);
      //double x = 0.0;
      //std::generate(uniform_binning.begin(), uniform_binning.end(), [&]{ return x++; }); 
      fdobj->x_var[iEvent] = &(duneobj->offaxis_2dbinnumber[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->dummy_y);
            //    doc["Binning"]["NTotalBins"].as<int>();

              //  SampleManager->raw()["Binning"]["BinningOpt"].as<int>();
      break;
      case 4: //case 4 useful for 3D binning, set the 'global bin number' in the yaml file from the NTotalBins options
      

      
      
      fdobj->x_var[iEvent] = &(duneobj->offaxis_2dbinnumber[iEvent]);
      fdobj->y_var[iEvent] = &(duneobj->dummy_y);
            //    doc["Binning"]["NTotalBins"].as<int>();

              //  SampleManager->raw()["Binning"]["BinningOpt"].as<int>();
      break;

    default:
      std::cout << "[ERROR:] " << __FILE__ << ":" << __LINE__
                << " unrecognised binning option" << BinningOpt << std::endl;
      throw;
      break;
    }
  }

  return;
}

void samplePDFDUNEBase::setupSplines(fdmc_base *fdobj, const char *splineFile,
                                     int nutype, int signal) {

  int nevents = fdobj->nEvents;
  std::cout << "##################" << std::endl;
  std::cout << "Initialising splines from file: " << (splineFile) << std::endl;
  std::cout << "##################" << std::endl;

  switch (BinningOpt) {
  case 0:
  case 1:
    fdobj->splineFile = new splinesDUNE((char *)splineFile, nutype, nevents,
                                        fdobj->SampleDetID, xsecCov);
    if (!(nutype == 1 || nutype == -1 || nutype == 2 || nutype == -2)) {
      std::cerr << "problem setting up splines in erec" << std::endl;
    }
    break;
  case 2:
    std::cout << "Creating splineDUNEBase" << std::endl;
    fdobj->splineFile =
        new splinesDUNE((char *)splineFile, nutype, nevents, (double)BinningOpt,
                        SampleDetID, xsecCov);
    if (!(nutype == 1 || nutype == -1 || nutype == 2 || nutype == -2)) {
      std::cerr << "problem setting up splines in erec" << std::endl;
    }
    break;
  default:
    break;
  }

  // ETA - Moved SetupSplineInfoArrays to be here
  fdobj->splineFile->SetupSplineInfoArray(xsecCov);
  fdobj->splineFile->SetSplineInfoArrays();

  return;
}

void samplePDFDUNEBase::applyShifts(int iSample, int iEvent) {

  // reset erec back to original value
  dunemcSamples[iSample].rw_erec_shifted[iEvent] =
      dunemcSamples[iSample].rw_erec[iEvent];

  // reset cvnnumu back to original value
  dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent] =
      dunemcSamples[iSample].rw_cvnnumu[iEvent];

  // reset cvnnue back to original value
  dunemcSamples[iSample].rw_cvnnue_shifted[iEvent] =
      dunemcSamples[iSample].rw_cvnnue[iEvent];

  // Calculate values needed
  double sqrtErecHad = sqrt(dunemcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep = sqrt(dunemcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunemcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunemcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunemcSamples[iSample].rw_eRecoP[iEvent] +
                   dunemcSamples[iSample].rw_eRecoPip[iEvent] +
                   dunemcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad = 1 / (sqrtErecHad + 0.1);
  double invSqrtErecLep = 1 / (sqrtErecLep + 0.1);
  double invSqrteRecoPi0 = 1 / (sqrteRecoPi0 + 0.1);
  double invSqrteRecoN = 1 / (sqrteRecoN + 0.1);
  double invSqrtSumEhad = 1 / (sqrtSumEhad + 0.1);

  bool CCnumu{dunemcSamples[iSample].rw_isCC[iEvent] == 1 &&
              abs(dunemcSamples[iSample].rw_nuPDG[iEvent] == 14) &&
              dunemcSamples[iSample].nutype == 2};
  bool CCnue{dunemcSamples[iSample].rw_isCC[iEvent] == 1 &&
             abs(dunemcSamples[iSample].rw_nuPDG[iEvent] == 12) &&
             dunemcSamples[iSample].nutype == 1};
  bool NotCCnumu{!(dunemcSamples[iSample].rw_isCC[iEvent] == 1 &&
                   abs(dunemcSamples[iSample].rw_nuPDG[iEvent] == 14)) &&
                 dunemcSamples[iSample].nutype == 2};

   /*
  TotalEScaleFD(FDDetectorSystPointers[0],
                &dunemcSamples[iSample].rw_erec_shifted[iEvent],
                dunemcSamples[iSample].rw_erec_had[iEvent],
                dunemcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);*/

//   TotalEScaleSqrtFD(FDDetectorSystPointers[1],
//                     &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                     dunemcSamples[iSample].rw_erec_had[iEvent],
//                     dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad,
//                     sqrtErecLep, NotCCnumu);

//   TotalEScaleInvSqrtFD(FDDetectorSystPointers[2],
//                        &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                        dunemcSamples[iSample].rw_erec_had[iEvent],
//                        dunemcSamples[iSample].rw_erec_lep[iEvent],
//                        invSqrtErecHad, invSqrtErecLep, NotCCnumu);

//   HadEScaleFD(FDDetectorSystPointers[3],
//               &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

//   HadEScaleSqrtFD(FDDetectorSystPointers[4],
//                   &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad,
//                   sqrtSumEhad);

//   HadEScaleInvSqrtFD(FDDetectorSystPointers[5],
//                      &dunemcSamples[iSample].rw_erec_shifted[iEvent], sumEhad,
//                      invSqrtSumEhad);

//   MuEScaleFD(FDDetectorSystPointers[6],
//              &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//              dunemcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

//   MuEScaleSqrtFD(FDDetectorSystPointers[7],
//                  &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                  dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep,
//                  CCnumu);

//   MuEScaleInvSqrtFD(FDDetectorSystPointers[8],
//                     &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                     dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep,
//                     CCnumu);

//   NEScaleFD(FDDetectorSystPointers[9],
//             &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//             dunemcSamples[iSample].rw_eRecoN[iEvent]);

//   NEScaleSqrtFD(FDDetectorSystPointers[10],
//                 &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                 dunemcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

//   NEScaleInvSqrtFD(FDDetectorSystPointers[11],
//                    &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                    dunemcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

//   EMEScaleFD(FDDetectorSystPointers[12],
//              &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//              dunemcSamples[iSample].rw_eRecoPi0[iEvent],
//              dunemcSamples[iSample].rw_erec_lep[iEvent], CCnue);

//   EMEScaleSqrtFD(FDDetectorSystPointers[13],
//                  &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                  dunemcSamples[iSample].rw_eRecoPi0[iEvent],
//                  dunemcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep,
//                  sqrteRecoPi0, CCnue);

//   EMEScaleInvSqrtFD(FDDetectorSystPointers[14],
//                     &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//                     dunemcSamples[iSample].rw_eRecoPi0[iEvent],
//                     dunemcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep,
//                     invSqrteRecoPi0, CCnue);

//   HadResFD(FDDetectorSystPointers[15],
//            &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//            dunemcSamples[iSample].rw_eRecoP[iEvent],
//            dunemcSamples[iSample].rw_eRecoPip[iEvent],
//            dunemcSamples[iSample].rw_eRecoPim[iEvent],
//            dunemcSamples[iSample].rw_eP[iEvent],
//            dunemcSamples[iSample].rw_ePip[iEvent],
//            dunemcSamples[iSample].rw_ePim[iEvent]);

//   MuResFD(FDDetectorSystPointers[16],
//           &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//           dunemcSamples[iSample].rw_erec_lep[iEvent],
//           dunemcSamples[iSample].rw_LepE[iEvent], CCnumu);

//   NResFD(FDDetectorSystPointers[17],
//          &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//          dunemcSamples[iSample].rw_eRecoN[iEvent],
//          dunemcSamples[iSample].rw_eN[iEvent]);

//   EMResFD(FDDetectorSystPointers[18],
//           &dunemcSamples[iSample].rw_erec_shifted[iEvent],
//           dunemcSamples[iSample].rw_eRecoPi0[iEvent],
//           dunemcSamples[iSample].rw_ePi0[iEvent],
//           dunemcSamples[iSample].rw_erec_lep[iEvent],
//           dunemcSamples[iSample].rw_LepE[iEvent], CCnue);

//   CVNNumuFD(FDDetectorSystPointers[19],
//             &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent]);

//   CVNNueFD(FDDetectorSystPointers[20],
//            &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent]);
}


// This is currently here just for show. We'll implement functional parameters
// soon!

int getNMCSamples() {
    //std::vector<struct dunemc_base> dunemcSamples;
    // Initialize dunemcSamples as needed
    //return dunemcSamples.size();
    return 0;
}



double samplePDFDUNEBase::CalcXsecWeightFunc(int iSample, int iEvent) {
  return 1.0;
}

std::vector<double> samplePDFDUNEBase::ReturnKinematicParameterBinning(
    std::string KinematicParameterStr) {
  std::cout << "ReturnKinematicVarBinning" << std::endl;
  std::vector<double> binningVector;
  return binningVector;
}

double samplePDFDUNEBase::getDiscVar(int iSample, int iEvent, int varindx) {
  std::cout << "getDiscVar" << std::endl;
  return 0.0;
}

double samplePDFDUNEBase::getCovLikelihood() {
  std::cout << "getCovLikelihood" << std::endl;
  return 0.0;
}

void samplePDFDUNEBase::printPosteriors() {
  std::cout << "printPosteriors" << std::endl;
}

/*
int samplePDFDUNEBase::getNMCSamples() {
  std::cout << "getNMCSamples" << std::endl;
  return 0;
}*/

int samplePDFDUNEBase::getNEventsInSample(int sample) {
  std::cout << "getNEventsInSample" << std::endl;
  return 0;
}

TH1D* get1DVarHist(std::string KinematicVar1, std::vector< std::vector<double> > SelectionVec, int kModeToFill, int kChannelToFill, int WeightStyle, TAxis* Axis) {
  bool fChannel;
  bool fMode;
  

  if (kChannelToFill!=-1) {
    if (kChannelToFill>getNMCSamples()) {
      std::cout << "Required channel is not available. kChannelToFill should be between 0 and " << getNMCSamples() << std::endl;
      std::cout << "kChannelToFill given:" << kChannelToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill!=-1) {
    if (kModeToFill>kMaCh3_nModes) {
      std::cout << "Required mode is not available. kModeToFill should be between 0 and " << kMaCh3_nModes << std::endl;
      std::cout << "kModeToFill given:" << kModeToFill << std::endl;
      std::cout << "Exitting.." << std::endl;
      throw;
    }
    fMode = true;
  } else {
    fMode = false;
  }

 // std::vector< std::vector<double> > SelectionVec;

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
  TFile my1dhisto("my1dhisto.root","RECREATE");
  onedim_binnumberhisto->Write();
  my1dhisto.Close();
  //return get1DVarHist(KinematicVar1,SelectionVec,WeightStyle,Axis);
  return get1DVarHist(KinematicVar1, SelectionVec, kModeToFill, kChannelToFill, WeightStyle, Axis);
}

//TFile my1dhisto("my1dhisto.root","RECREATE");
//std::ofstream out("1Dhisto.root".c_str()); // create a new output file or overwrite an existing one
//onedim_binnumberhisto->Write();
//my1dhisto.close();


/*
//1D hist for kinemtaticvariable 3D binning
  double noofbins = Abis3DHistogram->GetNcells();
    f1DHist_binnumber =
        std::make_unique<TH1D>("f1DHist_binnumbe", "",noofbins ,
        Erec.data()) ;*/
