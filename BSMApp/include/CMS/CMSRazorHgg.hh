//-------------------------------------------------------
// Description:  CMS  Razor analysis 
// Authors:      Maurizio Pierini, CERN
// Authors:      Christopher Rogan, Caltech
//-------------------------------------------------------

///  The CMSRazorHgg class implements the Razor-based
///  SUSY search by the CMS collaboration.
///  http://arxiv.org/abs/1006.2727
///  MISSING ANALYSIS REFERENCE HERE

#ifndef CMSRazorHgg_h
#define CMSRazorHgg_h

#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "CMS/CMSReco.hh"
#include "CommonTools/DataAnalysis.hh"
#include "CMS/CMSHemisphere.hh"

using namespace std;

class CMSRazorHgg : public CMSReco, public DataAnalysis {
public:
  
  //! constructor
  CMSRazorHgg(TTree *tree, double Lumi, double filterEff, double xsecMax, string analysis, bool delphesFormat);
  //! destructor
  virtual ~CMSRazorHgg();
  //! loop over events
  void Loop(string outFileName);
  void SetSqrts(double sqrts);
  double DeltaPhi(TLorentzVector jet1, TLorentzVector jet2);
private:
  
  /// Return the posterior distribution for the inclusive xsec,
  /// given a pdf and efficiency for one of the boxes;
  TH1D* XsecProb(TH1D* sigPdf, double eff, string Filename, string directory, int ibin, double xmin, double xmax, bool expected, bool doubleErr);
  /// Luminosity
  double _Lumi;
  // Gen-level filter efficiency
  double _filterEff;
  // Max xsec for plots [pb]
  double _xsecMax;
  // collision energy
  double _sqrts;
  // boolean to use Delphes format
  double _delphesFormat;
};
#endif
