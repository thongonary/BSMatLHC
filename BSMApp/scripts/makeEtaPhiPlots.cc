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

    // attributes of the PFAK04 jets
    double PFJets_phi[100];
    double PFJets_eta[100];
    double PFJets_pT[100];

    //plotting arrays
    vector<double> phi_points;
    vector<double> eta_points;
    vector<double> pt_points;
    vector<double> phi_points1;
    vector<double> eta_points1;

    //megajet hemisphere constituents
    int Default_hem1_csts[400];
    int Default_hem2_csts[4000];

    //new algorithm hemisphere constituents
    int newjetcam_hem1_csts[1000];
    int newjetcam_hem2_csts[800];

    //generator-level "hemisphere" constituents
    int gen_hem1[200];
    int gen_hem2[3000];
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

	// which event to look at 
	inTree->GetEntry(280); 
	int count = 0;
	// print out each jet and its phi, eta, pt values
	while (PFJets_phi[count] > -900){
	  phi_points.push_back(PFJets_phi[count]);
	  cout << "phi " << count << ": " << PFJets_phi[count] << endl;
	  eta_points.push_back(PFJets_eta[count]);
	  cout << "eta " << count << ": " <<PFJets_eta[count] << endl;
	  pt_points.push_back(PFJets_pT[count]);
	  cout << "pt " << count << ": " <<PFJets_pT[count] << endl;
	  count++;				
	}
	// print out hemispheres for the "default" = megajet
	count = 0;
	while (Default_hem1_csts[count]> -900){
	  cout << "Default number " << Default_hem1_csts[count] << endl;
	  count++;
	}
	count = 0;
	while (Default_hem2_csts[count]> -900){
	  cout << "Default number2 " << Default_hem2_csts[count] << endl;
	  count++;
	}
	// print out hemispheres for the "new algorithm"
	count = 0;
	while (newjetcam_hem1_csts[count]> -900){
	  cout << "newjetcam number " << newjetcam_hem1_csts[count] << endl;
	  count++;
	}
	count = 0;
	while (newjetcam_hem2_csts[count]> -900){
	  cout << "newjetcam number2 " << newjetcam_hem2_csts[count] << endl;
	  count++;
	}
	// print out generator-level matched info
	count = 0;
	while (gen_hem1[count]> -900){
	  cout << "genhem1 " << gen_hem1[count] << endl;
	  count++;
	  }
	count = 0;
	while (gen_hem2[count]> -900){
	  cout << "genhem2 " << gen_hem2[count] << endl;
	  count++;
	  }
 	  
    
    c1->cd();
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);

    }
    //plotting
    double phi[phi_points.size()];
    double eta[eta_points.size()];
    double pt[pt_points.size()];
    for (int i = 0; i < phi_points.size(); i++){
      phi[i] = phi_points[i];
    }
    for (int i = 0; i < eta_points.size(); i++){
      eta[i] = eta_points[i];
    }
    for (int i = 0; i < pt_points.size(); i++){
      pt[i] = pt_points[i];
    }
    TMultiGraph *mg = new TMultiGraph();
    TGraph *gr1 = new TGraph (phi_points.size(),phi,eta);
    gr1->GetXaxis()->SetTitle("#phi");
    gr1->GetYaxis()->SetTitle("#eta");
    gr1->GetYaxis()->SetTitleSize(0.06);
    gr1->GetYaxis()->SetLabelSize(0.06);
    gr1->GetYaxis()->SetTitleOffset(1.4);
    gr1->SetTitle("");
    gr1->SetMarkerStyle(12);
    gr1->SetMarkerSize(3);
    gr1->SetMarkerColor(9);
    gr1->Draw("AP");
    c1->Print("eta_phi.pdf");
    delete c1;
    for(int i = 0; i < inFiles.size(); i++) delete inFiles[i];
    return 0;
}

