//Run with root -l .x makeEtaPhi.cc
//
//Creates 2D eta phi plot.  Plot shows which 
//PFAK04 jets are clustered into which hemisphere
//with each algorithm.  Blue is newjetcam and red 
//is either megajet or Georgi.  Different shapes 
//represent different hemispheres.  Each point is
//annotated with pt, mass.  The generator-level 
//quark information is shown with a green shape.  
//Each point is labeled with pt and mass information, 
//and the MR value for each algorithm is near the title. 
//Each file has is saved as with its event number at the 
//end.


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
#include "TVectorD.h"
#define _USE_MATH_DEFINES
using namespace std;

int makeEtaPhi(){
  //beginning event
  int entry_start = 188;
  //number of events to loop through.
  int number_of_events = 1;
  //loop over events
  for (int c = 0; c < number_of_events; c ++){
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

    //vector<TFile*> inFiles;
    //TFile *New = new TFile("../smallsamples/output.root"); 
    TFile *New = new TFile("../smallsamples/T1_G_1.5.root"); 
    //TFile *New = new TFile("../smallsamples/T1_G_100.root"); 
    TTree *inTree;

    double MR_Default, MR_newjetcam;
    // attributes of the PFAK04 jets
    double PFJets_phi[100];
    double PFJets_eta[100];
    double PFJets_pT[100];
    double PFJets_m[100];
    
    //plotting arrays
    vector<double> phi_points;
    vector<double> eta_points;
    vector<double> pt_points;
    vector<double> m_points;
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
    New->cd();
	double phi_temp;
	double eta_temp;
	double pt_temp;
        inTree = (TTree*)New->Get("RazorInclusive");
	inTree->SetBranchAddress("PFJets_phi", &PFJets_phi);
	inTree->SetBranchAddress("PFJets_eta", &PFJets_eta);
	inTree->SetBranchAddress("PFJets_pT", &PFJets_pT);
	inTree->SetBranchAddress("PFJets_m", &PFJets_m);
	inTree->SetBranchAddress("Default_hem1_csts", &Default_hem1_csts);
	inTree->SetBranchAddress("Default_hem2_csts", &Default_hem2_csts);
	inTree->SetBranchAddress("newjetcam_hem1_csts", &newjetcam_hem1_csts);
	inTree->SetBranchAddress("newjetcam_hem2_csts", &newjetcam_hem2_csts);
	inTree->SetBranchAddress("gen_hem1", &gen_hem1);
	inTree->SetBranchAddress("gen_hem2", &gen_hem2);
	inTree->SetBranchAddress("MR_Default", &MR_Default);
	inTree->SetBranchAddress("MR_newjetcam", &MR_newjetcam);
	

	// which event to look at 
	inTree->GetEntry(entry_start + c); 
	int count = 0;
	//values of -999 are invalid
	while (PFJets_phi[count] > -900){
	  phi_points.push_back(PFJets_phi[count]);
	  eta_points.push_back(PFJets_eta[count]);
	  pt_points.push_back(PFJets_pT[count]);
	  m_points.push_back(PFJets_m[count]);
	  count++;				
	}
	const int phisize = phi_points.size();
	//adds the phi and eta coordinates to each hemisphere
	//and adds mass and pt information (only for Default) 
	count = 0;
        double Default1_phi[phisize], Default1_eta[phisize], Default1_pt[phisize], 
	  Default1_m[phisize], Default2_phi[phisize], Default2_eta[phisize], 
	  Default2_pt[phisize], Default2_m[phisize];
	while (Default_hem1_csts[count]> -900){
	  Default1_phi[count] = phi_points[Default_hem1_csts[count]];
	  Default1_eta[count] = eta_points[Default_hem1_csts[count]];
	  Default1_pt[count] = pt_points[Default_hem1_csts[count]];
	  Default1_m[count] = m_points[Default_hem1_csts[count]];
	  count++;
	}
	int Default_hem1_len = count;
	count = 0;
	while (Default_hem2_csts[count]> -900){
	  Default2_phi[count] = phi_points[Default_hem2_csts[count]];
	  Default2_eta[count] = eta_points[Default_hem2_csts[count]];
	  Default2_pt[count] = pt_points[Default_hem2_csts[count]];
	  Default2_m[count] = m_points[Default_hem2_csts[count]];
	  count++;
	}
	int Default_hem2_len = count;
        double newjetcam1_phi[phisize], newjetcam2_phi[phisize],
	  newjetcam1_eta[phisize], newjetcam2_eta[phisize];
	count = 0;
	while (newjetcam_hem1_csts[count]> -900){
	  newjetcam1_phi[count] = phi_points[newjetcam_hem1_csts[count]];
	  newjetcam1_eta[count] = eta_points[newjetcam_hem1_csts[count]];
	  count++;
	}
	int newjetcam_hem1_len = count;
	count = 0;
	while (newjetcam_hem2_csts[count]> -900){
	  newjetcam2_phi[count] = phi_points[newjetcam_hem2_csts[count]];
	  newjetcam2_eta[count] = eta_points[newjetcam_hem2_csts[count]];
	  count++;
	}
	int newjetcam_hem2_len = count;
	double gen_hem1_phi[phisize], gen_hem2_phi[phisize], 
	  gen_hem1_eta[phisize], gen_hem2_eta[phisize];
	count = 0;
	while (gen_hem1[count]> -900){
	  gen_hem1_phi[count] = phi_points[gen_hem1[count]];
	  gen_hem1_eta[count] = eta_points[gen_hem1[count]];
	  count++;
	  }
	int gen_hem1_len = count;
	count = 0;
	while (gen_hem2[count]> -900){
	  gen_hem2_phi[count] = phi_points[gen_hem2[count]];
	  gen_hem2_eta[count] = eta_points[gen_hem2[count]];
	  count++;
	  }
	int gen_hem2_len = count;

    c1->cd();
    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.18);

 
    //plotting

    double phi[phisize];
    double eta[phisize];
    double pt[phisize];
    for (int i = 0; i < phi_points.size(); i++){
      phi[i] = phi_points[i];
    }
    for (int i = 0; i < eta_points.size(); i++){
      eta[i] = eta_points[i];
    }
    for (int i = 0; i < pt_points.size(); i++){
      pt[i] = pt_points[i];
    }

    //Make 6 new graphs, 2 hemispheres for each algorithm, 
    //and 2 hemispheres for gen level information
    TMultiGraph *mg = new TMultiGraph();
    TGraph *gr1 = new TGraph (newjetcam_hem1_len, newjetcam1_eta, newjetcam1_phi);
    TGraph *gr2 = new TGraph (newjetcam_hem2_len, newjetcam2_eta, newjetcam2_phi);
    TGraph *gr3 = new TGraph (Default_hem1_len, Default1_eta, Default1_phi);
    TGraph *gr4 = new TGraph (Default_hem2_len, Default2_eta, Default2_phi);
    TGraph *gr5 = new TGraph (gen_hem1_len, gen_hem1_eta, gen_hem1_phi);    
    TGraph *gr6 = new TGraph (gen_hem2_len, gen_hem2_eta, gen_hem2_phi);    
    
    
    gr1->SetMarkerStyle(29);
    gr1->SetMarkerSize(3);
    gr1->SetMarkerColor(4);
    
    gr2->SetMarkerStyle(34);
    gr2->SetMarkerSize(3);
    gr2->SetMarkerColor(4);

    gr3->SetMarkerStyle(5);
    gr3->SetMarkerSize(3);
    gr3->SetMarkerColor(2);

    gr4->SetMarkerStyle(2);
    gr4->SetMarkerSize(3);
    gr4->SetMarkerColor(2);

    gr5->SetMarkerStyle(24);
    gr5->SetMarkerSize(8);
    gr5->SetMarkerColor(30);

    gr6->SetMarkerStyle(25);
    gr6->SetMarkerSize(8);
    gr6->SetMarkerColor(30);

    //add pt labels to the points.  Only done from Default.  
    for (int i = 0; i < Default_hem1_len; i++){
      //converts pt and mass values into a string
      std::ostringstream stream1;
      stream1 << int(Default1_pt[i]);// << ", " << (Default1_m[i]);
      const char* pt_out = (stream1.str()).c_str();
      //Annotates with a slight offset from the point's coordinate
      TLatex *latex = new TLatex(gr3->GetX()[i] + 0.07, gr3->GetY()[i] + 0.07, pt_out);
      latex->SetTextSize(0.04);
      gr3->GetListOfFunctions()->Add(latex);
    }
    //same thing for the other hemisphere
    for (int i = 0; i < Default_hem2_len; i++){
      std::ostringstream stream2;
      stream2 << int(Default2_pt[i]);// << ", " << (Default2_m[i]);
      const char* pt_out2 = (stream2.str()).c_str();
      TLatex *latex = new TLatex(gr4->GetX()[i] + 0.07, gr4->GetY()[i] + 0.07, pt_out2);
      latex->SetTextSize(0.04);
      gr4->GetListOfFunctions()->Add(latex);
    }

    mg->Add(gr1, "AP");
    mg->Add(gr2, "AP");
    mg->Add(gr3, "AP");
    mg->Add(gr4, "AP");
    mg->Add(gr5, "AP");
    mg->Add(gr6, "AP");


    mg->SetTitle("T1 1400 100 Georgi 1.5");

    //Adds the MR values of the two algorithms near the title
    
    std::ostringstream MRStream;
    MRStream << "MR Default: " << int(MR_Default) << "            MR new: " << int(MR_newjetcam);
    const char* MR_out = (MRStream.str()).c_str();
    TLatex *latex = new TLatex(-2.25, 3.40, MR_out);
    latex->SetTextSize(0.04);
    mg->GetListOfFunctions()->Add(latex);
    //Ends the MR adding

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("#eta");
    mg->GetYaxis()->SetTitle("#phi");
    mg->GetYaxis()->SetRangeUser(-3.3, 3.3);
    mg->GetXaxis()->SetLimits(-3, 3);
    
    //mg->SetAxisRange(0., 3.,"Y");
    gPad->Modified();

    std::ostringstream stream;
    //Creates a stream to make a string that changes with each event#
    //so that each file can have a different saved name
    stream << "8-29-2014/Eta_Phi/Georgi1.5/G1.5_"<< entry_start + c<<".png";
    const char* savename = (stream.str()).c_str();
    cout <<savename<< endl;
    c1->Print(savename);
    delete c1;
    delete New;
}
  return 0;
}
