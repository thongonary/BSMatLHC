//-------------------------------------------------------
// Description:  CMS Razor analysis for 13 TeV
// Authors:      Javier Duarte, Caltech
//-------------------------------------------------------

#ifndef CMSRazor13TeV_h
#define CMSRazor13TeV_h

#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

#include "CMS/CMSReco.hh"
#include "CommonTools/DataAnalysis.hh"
#include "CMS/CMSHemisphere.hh"

using namespace std;



class CMSRazor13TeV : public CMSReco, public DataAnalysis {
public:
  
  //! constructor
  CMSRazor13TeV(TTree *tree, double Lumi, string analysis);
  //! destructor
  virtual ~CMSRazor13TeV();
  //! loop over events
	
  void Loop(string outFileName);
  void SetSqrts(double sqrts);

private:
  /// Luminosity
  double _Lumi;
  ///collision energy
  double _sqrts;
  int _event_counter;
	
  ///cambridge exclusive jets implementation
  vector<fastjet::PseudoJet> improper_Exclusive_Jets(const fastjet::ClusterSequence cs_E);
	
	// theses are used in gen level analysis
  //vector<int> find_all(const vector<double> objects, const item, const low_bound, const up_bound );
  double D_between_vectors( const TLorentzVector V1, const TLorentzVector V2);

};
#endif
