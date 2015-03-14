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
  CMSRazorHgg(TTree *tree, double Lumi, string analysis);
  //! destructor
  virtual ~CMSRazorHgg();
  //! loop over events
  void Loop(string outFileName);
  void SetSqrts(double sqrts);
  double DeltaPhi(TLorentzVector jet1, TLorentzVector jet2);
private:
  /// Luminosity
  double _Lumi;

  // collision energy
  double _sqrts;

};
#endif
