#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

#define _USE_MATH_DEFINES

using namespace std;

//Simple script to for a rough evaluation of clustering performance
//Calculates a success rate, where a failure only happens if an algorithm clusters jets that "match" with two different gen-level decay chains

//Inputs: root file you are evaluating
//Currently evaluates only two algorithms: megajet and new cambridge inclusive clustering 

int main(){

    //style tips
    gStyle->SetOptStat(000000);
    gStyle->SetPalette(51, 0);
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetTitleFont(132,"xyz");
    gStyle->SetTitleSize(0.06);
    gStyle->SetTextFont(132);
    gStyle->SetTextSize(0.08);
    gStyle->SetStatFont(132);
    gStyle->SetLabelFont(132,"xyz");
    gStyle->SetLabelSize(0.06);

    //make canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->cd();

    vector<TFile*> inFiles;
    inFiles.push_back(new TFile("T1_1400_100_gen_new.root"));;


    TTree *inTree;

    double PFJets_phi[100];
    double PFJets_eta[100];
    double PFJets_pT[100];
    int Default_hem1_csts[400];
    int Default_hem2_csts[4000];
    int newjetcam_hem1_csts[1000];
    int newjetcam_hem2_csts[800];
    int gen_hem1[200];
    int gen_hem2[3000];
    int nEvents;
    int new_matches = 0;
    int def_matches = 0;
    int new_failures = 0;
    int def_failures = 0;
    int array_len;
    TFile outFile("jet_pts.root", "recreate");
    for(int i = 0; i < inFiles.size(); i++){
        inFiles[i]->cd();
	double phi_temp;
	double eta_temp;
	double pt_temp;
        inTree = (TTree*)inFiles[i]->Get("RazorInclusive");
	inTree->SetBranchAddress("PFJets_phi", &PFJets_phi);
	inTree->SetBranchAddress("PFJets_eta", &PFJets_eta);
	inTree->SetBranchAddress("PFJets_pT", &PFJets_pT);
	inTree->SetBranchAddress("Default_hem1_csts", &Default_hem1_csts);
	inTree->SetBranchAddress("Default_hem2_csts", &Default_hem2_csts);
	inTree->SetBranchAddress("newjetcam_hem1_csts", &newjetcam_hem1_csts);
	inTree->SetBranchAddress("newjetcam_hem2_csts", &newjetcam_hem2_csts);
	inTree->SetBranchAddress("gen_hem1", &gen_hem1);
	inTree->SetBranchAddress("gen_hem2", &gen_hem2);
        //loop over events
	nEvents = inTree->GetEntries();
        for (int iEntry = 0; iEntry < nEvents; iEntry++){
	  inTree->GetEntry(iEntry);
	  int count = 0;
	  bool cam1_g2_match = false;
	  bool cam1_g1_match = false;
	  bool cam2_g2_match = false;
	  bool cam2_g1_match = false;
	  int count2 = 0;

	  // loop through constituents and ask if there is any clustering with crossing decays 
	  while (newjetcam_hem1_csts[count] > -900){
	    while (gen_hem1[count2] > -900){
	      if (newjetcam_hem1_csts[count] == gen_hem1[count2]){
		cam1_g1_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  count = 0;
	  count2 = 0;
	  while (newjetcam_hem1_csts[count] > -900){
	    while (gen_hem2[count2] > -900){
	      if (newjetcam_hem1_csts[count] == gen_hem2[count2]){
		cam1_g2_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  while (newjetcam_hem2_csts[count] > -900){
	    while (gen_hem1[count2] > -900){
	      if (newjetcam_hem2_csts[count] == gen_hem1[count2]){
		cam2_g1_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  count = 0;
	  count2 = 0;
	  while (newjetcam_hem2_csts[count] > -900){
	    while (gen_hem2[count2] > -900){
	      if (newjetcam_hem2_csts[count] == gen_hem2[count2]){
		cam2_g2_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  if ((cam1_g2_match && cam1_g2_match)||(cam1_g2_match && cam1_g1_match)) {new_failures++;}
	  else { 
	    new_matches++;}


	  count = 0;
	  bool def1_g2_match = false;
	  bool def1_g1_match = false;
	  bool def2_g2_match = false;
	  bool def2_g1_match = false;
	 count2 = 0;
	  while (Default_hem1_csts[count] > -900){
	    while (gen_hem1[count2] > -900){
	      if (Default_hem1_csts[count] == gen_hem1[count2]){
		def1_g1_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  count = 0;
	  count2 = 0;
	  while (Default_hem1_csts[count] > -900){
	    while (gen_hem2[count2] > -900){
	      if (Default_hem1_csts[count] == gen_hem2[count2]){
		def1_g2_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  while (Default_hem2_csts[count] > -900){
	    while (gen_hem1[count2] > -900){
	      if (Default_hem2_csts[count] == gen_hem1[count2]){
		def2_g1_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  count = 0;
	  count2 = 0;
	  while (Default_hem2_csts[count] > -900){
	    while (gen_hem2[count2] > -900){
	      if (Default_hem2_csts[count] == gen_hem2[count2]){
		def2_g2_match = true;
		break;
	      }
	      count2++;
	    }
	    count++;
	  }
	  if ((def1_g2_match && def1_g2_match)||(def1_g2_match && def1_g1_match)) {def_failures++;}
	  else { 
	    def_matches++;}
	}
    }

    cout << "//New Algorithm//" << endl;
    cout << "Successes: " << new_matches << endl;
    cout << "Failures: " << new_failures << endl;
    cout << "Events: " << nEvents << endl;
    double fraction = float(matches)/float(nEvents);
    cout << "Fraction: " << fraction << endl;
    
    cout << "//Megajet//" << endl;
    cout << "Successes: " << def_matches << endl;
    cout << "Failures: " << def_failures << endl;
    cout << "Events: " << nEvents << endl;
    fraction = float(def_matches)/float(nEvents);
    cout << "Fraction: " << fraction << endl;
    for(int i = 0; i < inFiles.size(); i++) delete inFiles[i];
    return 0;
}

