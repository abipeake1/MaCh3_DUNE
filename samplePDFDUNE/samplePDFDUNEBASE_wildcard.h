#ifndef _samplePDFDUNEBASE_wildcard_h_
#define _samplePDFDUNEBASE_wildcard_h_

#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TGraph2DErrors.h>
#include <vector>
#include <omp.h>
#include <list>

#include "splines/splinesDUNE.h"
#include "covariance/covarianceXsec.h"
#include "covariance/covarianceOsc.h"
#include "samplePDF/samplePDFFDBase.h"
#include "samplePDFDUNEBase.h"
#include "StructsDUNE.h"

class samplePDFDUNEBASE_wildcard : virtual public samplePDFFDBase
{
public:
  samplePDFDUNEBASE_wildcard(double pot, std::string mc_version, covarianceXsec* xsec_cov);
  ~samplePDFDUNEBASE_wildcard();
;}

#endif

//init same for all should just take a config with file names

//think about splines

//setcov functions should be doable

//calcoscweights all same

//reweight all the same

//getcovlikelihood 
