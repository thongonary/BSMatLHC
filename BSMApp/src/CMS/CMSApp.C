//-------------------------------------------------------
// Description:
// Authors:
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

// Supported Analyses
#include "CMS/CMSRazor.hh"
#include "CMS/CMSRazor13TeV.hh"
#include "CMS/CMSSUSYVars.hh"
//#include "CMS/CMSMonoJet.hh"
//#include "CMS/CMSSSDilepBtag.hh"
//#include "CMS/CMSDisplacedJet.hh"
// 
#include "CMS/CMSDarkMatter.hh"
#include "CMS/CMSSubstructure.hh"
#include "CMS/CMSRazorHgg.hh"
#include "CMS/CMSRazorHggHbb.hh"

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  if ( argc < 2 ){
    cout << "To run CMSApp please specify the input file" << endl; 
    cout << "Example:        ./CMSApp myBSMGenfile.root" << endl;
    cout << "OPTIONS:        " << endl;
    cout << "--verbose       Increase verbosity level for debug" << endl;
    cout << "--razor         Run Razor Analysis" << endl;
    cout << "--razor13       Run Razor 13 TeV Analysis" << endl;
    cout << "--hggrazor      Razor Hgg" << endl;
    cout << "--hgghbbrazor   Razor Hgg with Hbb converstion" << endl;
    cout << "--susyvars      Razor, alpha_T, MT2, etc" << endl;
    cout << "--monojet       Run MonoJet Analysis" << endl;
    cout << "--darkmatter    Run Dark Matter future study" << endl;
    cout << "--substructure  Run Substructure tagger analysis on jets" << endl;
    cout << "--ssdilepbtag   Run the SS dilepton btag Analysis" << endl;
    cout << "--displaced     Run DSisplaced Jet Analysis" << endl;
    cout << "--sqrts         sqrts of collisions" << endl;
    cout << "--filter        gen-level filter efficiency"  << endl;
    cout << "--lumi          luminosity in inverse picobarns" << endl;
    cout << "--output        Name of the output file (none created if not specified)" << endl;
    cout << "--delphes       Sets input tree structure to be Delphes" << endl;
    return 1;
  }

  char inputCMS[150];
  strcpy(inputCMS,argv[1]);

  char outFileName[150]= "none";

  bool verbose  = false;
  bool writeOut = false;
  bool razor = false;
  bool razor13 = false;
  bool susyvars = false;
  bool monojet = false;
  bool darkmatter = false;
  bool substructure = false;
  bool ssdilepbtag = false;
  bool displaced = false;
  bool razorhgg = false;
  bool razorhgghbb = false;
  bool delphesFormat = false;
  double sqrts = 13000.;
  double filter = 1.;
  double lumi = 19800.;

  for (int i=1;i<argc;i++){
    if (strncmp(argv[i],"--output",8)==0) {
      sscanf(argv[i],"--output=%s",outFileName);
      writeOut = true;
    }
    if (strncmp(argv[i],"--sqrts",7)==0)  {
      sscanf(argv[i],"--sqrts=%lf",&sqrts);
    }
    if (strncmp(argv[i],"--filter",8)==0)  {
      sscanf(argv[i],"--filter=%lf",&filter);
    }
    if (strncmp(argv[i],"--lumi",6)==0)  {
      sscanf(argv[i],"--lumi=%lf",&lumi);
    }
    if (strncmp(argv[i],"--verbose",9)==0)      verbose = true;
    if (strncmp(argv[i],"--monojet",9)==0)      monojet = true;
    if (strncmp(argv[i],"--delphes",9)==0)      delphesFormat = true;
    if (strncmp(argv[i],"--substructure",14)==0) substructure = true;
    if (strncmp(argv[i],"--darkmatter",12)==0)  darkmatter = true;
    if (strncmp(argv[i],"--displaced",11)==0)   displaced = true;
    if (strncmp(argv[i],"--ssdilepbtag",13)==0) ssdilepbtag = true;
    if (strncmp(argv[i],"--razor",7)==0){
      if (strncmp(argv[i],"--razor13",9)==0) razor13 = true;
      else razor = true;
    }
    if (strncmp(argv[i],"--susyvars",10)==0)        susyvars = true;
    if (strncmp(argv[i],"--hggrazor",10)==0)        razorhgg = true;
    if (strncmp(argv[i],"--hgghbbrazor",13)==0)        razorhgghbb = true;
  }
  
  if(strncmp(inputCMS,"none",4)!=0) {
    
    // RECO Tree
    TChain *cmsChain;
    if (delphesFormat) cmsChain = new TChain("Delphes");
    else cmsChain = new TChain("GenEvent");
    cmsChain->Add(inputCMS);

    // Open Output file
    if(writeOut) {
      TFile *file = new TFile(outFileName,"RECREATE");
      file->Close();
    }


    // Razor 13 TeV analysis                                                                                                                                                              
    if(razor13) {
      CMSRazor13TeV cmsRazor13(cmsChain, 0., "");
      if(!writeOut) {
        cout << "please specify output file" << endl;
        return 0;
      }
      if(verbose) cmsRazor13.SetVerbose(true);
      cmsRazor13.SetSqrts(sqrts);
      cmsRazor13.Loop(outFileName);
    }


    // SUS-12-005: Razor analysis
    if(razor) {
      //      CMSRazor cmsRazor(cmsChain, 4600., "CMSRazor_HadLik_SUS_12_005");
      CMSRazor cmsRazor(cmsChain, 0., "");
      if(!writeOut) {
	cout << "please specify output file" << endl;
	return 0;
      }
      if(verbose) cmsRazor.SetVerbose(true);
      cmsRazor.SetSqrts(sqrts);
      cmsRazor.Loop(outFileName);
    }

    //Razor + alpha_T, MT2, xE, ptOut
    if(susyvars){
        CMSSUSYVars cmssusyvars(cmsChain, 0., "");
        if(!writeOut){
            cout << "please specify output file" << endl;
            return 0;
        }
        if(verbose) cmssusyvars.SetVerbose(true);
        cmssusyvars.SetSqrts(sqrts);
        cmssusyvars.Loop(outFileName);
    }

    if(razorhgg)
      {
	//CMSRazorHgg* cmsrazorhgg = new CMSRazorHgg(cmsChain, 19800., 0.0044280, "CMSRazorHgg_Lik_SUS_14_017", delphesFormat);
	CMSRazorHgg* cmsrazorhgg = new CMSRazorHgg(cmsChain, lumi, filter, "CMSRazorHgg_Lik_SUS_14_017", delphesFormat);
	if(!writeOut){
	  cout << "please specify output file" << endl;
	  return 0;
	}
	if(verbose) cmsrazorhgg->SetVerbose(true);
	cmsrazorhgg->SetSqrts(sqrts);
	cmsrazorhgg->Loop(outFileName);
	//delete cmsrazorhgg;
      }
    
    if(razorhgghbb){
      CMSRazorHggHbb cmsrazorhgghbb(cmsChain, 0., "");
      if(!writeOut){
	cout << "please specify output file" << endl;
	return 0;
      }
      if(verbose) cmsrazorhgghbb.SetVerbose(true);
      cmsrazorhgghbb.SetSqrts(sqrts);
      cmsrazorhgghbb.Loop(outFileName);
    }
    
    /*    
    // EXO-11-059: MonoJet analysis
    if(monojet) {
    CMSMonoJet cmsMonoJet(cmsChain, inputCMS, 4980., "CMSmonojet_EXO_11_059");
    if(verbose) cmsMonoJet.SetVerbose(true);
    cmsMonoJet.Loop();
    }
    
    // SUS-11-020: SUSY Dilep SS + Btag
    if(ssdilepbtag) {
      CMSSSDilepBtag cmsSSDilepBtag(cmsChain, inputCMS, 4980., "CMSSSDilepBTag_SUS_11_020");
      if(verbose) cmsSSDilepBtag.SetVerbose(true);
      cmsSSDilepBtag.Loop();
      }
      
      // Displaced search
      if(displaced) {
      CMSDisplacedJet cmsDisplacedJet(cmsChain);
      if(!writeOut) {
      cout << "please specify output file" << endl;
      return 0;
      }
      cmsDisplacedJet.Loop(outFileName);
      }
    */
    
    // Dark Matter future studies
    if(darkmatter) {
      CMSDarkMatter cmsDarkMatter(cmsChain, inputCMS, 4980., "CMSmonojet_EXO_11_059");
      if(verbose) cmsDarkMatter.SetVerbose(true);
      cmsDarkMatter.Loop(outFileName);
    }
    
    // Dark Matter future studies
    if(substructure) {
      CMSSubstructure cmsSubstructure(cmsChain, inputCMS, 4980., "CMSsubstructure");
      if(verbose) cmsSubstructure.SetVerbose(true);
      cmsSubstructure.Loop(outFileName);
    }
    
  }

  return 0;
}
