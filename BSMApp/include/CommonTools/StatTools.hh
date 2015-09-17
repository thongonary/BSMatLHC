//-----------------------------------------------------------------
// Description:  minimal set of statistics tools
// Authors:      Maurizio Pierini, CERN
//-----------------------------------------------------------------

#ifndef StatTools_h
#define StatTools_h

#include <TH2D.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>

#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

using namespace std;

class StatTools {
public:
  
  //! constructor
  StatTools(int seed);
  //! destructor
  virtual ~StatTools();
  /// Find excluded upper value for s, given a likelihood of s*scale and a Credibility Level
  double FindUL(TH1D* Lik, double CL, double scale);
  
  /// perform hit or miss
  bool HitOrMiss(double y);

  /// Provide the posterior on a signal s, inetgrating out the nuisance parameter
  /// b (background) out of a likelihood model
  /// L = P(n|s+b)*G(b|CenB, SigB)
  /// where CenB +/- SigB is the expected background and n is the observed yield. A flat prior on b is assumed
  TH1D* LogNormPoissonConv(TString name, Double_t n, Double_t CenB, Double_t SigB,  Double_t smin = 0., Double_t smax = 100.);

  /// Convolution of a given pdf with a log-normal of a given error sigma,
  /// representing the systematic error on the signal
  /// here sigma is given in [0,1], representing the relative error
  TH1D* LogNormHistoConv(TH1D* histo, double sigma);

protected:

private:

};
#endif


#ifndef PoissonGaussian_h
#define PoissonGaussian_h

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class PoissonGaussian: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double b = x;
     double CenB = p[0];
     double SigB = p[1];
     double s = p[2];
     double n = p[3];
     return TMath::Poisson(n, s+b)*TMath::Gaus(b,CenB,SigB,true);
   }
   
   double DoEval(double x) const
   {     
     double b = x;
     double CenB = pars[0];
     double SigB = pars[1];
     double s = pars[2];
     double n = pars[3];
     return TMath::Poisson(n, s+b)*TMath::Gaus(b,CenB,SigB,true);
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new PoissonGaussian();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 4;
   }
};


#endif
