//-------------------------------------------------------
// Description:
// Authors:
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "CommonTools/CountingExperiment.hh"
#include "CMS/CMSRazorLikelihood.hh"

using namespace std;

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  string InputFileName("none");
  for(int i=0; i< argc; i++) {
    if (string(argv[i]) == "-i") {
      if(i+1 < argc) InputFileName = argv[i+1];
    }
  }

  // read counting experiments from input file
  if(InputFileName != "none") {
    // open file
    ifstream myfile(InputFileName.c_str());
    if (myfile.is_open()) {
      string rootFileName;
      string analysisName;
      int n;
      double b, db;
      double mins, maxs;
      string line;
      string oldrootFileName("none");
      while(myfile.good()) {
	getline(myfile,line);
	if(line.find("#") == 0) continue;
	stringstream ss(line);
	ss >> rootFileName >> analysisName >> n >> b >> db >> mins >> maxs;
	CountingExperiment* cexp = new CountingExperiment(n, b, db, mins, maxs);
	if(oldrootFileName != rootFileName) {
	  cexp->CreatePosteriors(rootFileName, analysisName, "recreate");
	  oldrootFileName = rootFileName;
	} else {
	  cexp->CreatePosteriors(rootFileName, analysisName, "update");
	}
      }
    }
    myfile.close();
  }

  // Razor analysis
  CMSRazorLikelihood cmsRazorLikelihoodHighRes("data/ExpectedObserved_RazorHgg_HighRes_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodHighRes.CreatePosteriors("data/CMSRazorHgg_HighResLik_SUS_14_017.root","HighRes");
  
  CMSRazorLikelihood cmsRazorLikelihoodLowRes("data/ExpectedObserved_RazorHgg_LowRes_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodLowRes.CreatePosteriors("data/CMSRazorHgg_LowResLik_SUS_14_017.root","LowRes");
  
  CMSRazorLikelihood cmsRazorLikelihoodHighPt("data/ExpectedObserved_RazorHgg_HighPt_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodHighPt.CreatePosteriors("data/CMSRazorHgg_HighPtLik_SUS_14_017.root","HighPt");
  
  CMSRazorLikelihood cmsRazorLikelihoodHbb("data/ExpectedObserved_RazorHgg_Hbb_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodHbb.CreatePosteriors("data/CMSRazorHgg_HbbLik_SUS_14_017.root","Hbb");
  
  CMSRazorLikelihood cmsRazorLikelihoodZbb("data/ExpectedObserved_RazorHgg_Zbb_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodZbb.CreatePosteriors("data/CMSRazorHgg_ZbbLik_SUS_14_017.root","Zbb");
  
  CMSRazorLikelihood cmsRazorLikelihoodTotal("data/ExpectedObserved_RazorHgg_Total_Summer2015.root", "hOBS", "hEXP");
  cmsRazorLikelihoodTotal.CreatePosteriors("data/CMSRazorHgg_TotalLik_SUS_14_017.root","Total");
  
  return 0;
}
