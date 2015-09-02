//plots mean, mode, and RMS of any value calculated across many different mass points.
//median is currently set up but not implemented  
#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"

#define _USE_MATH_DEFINES

const int p = 11;
double MJ_array[p], Georgi_array[p], New_array[p], KT_array[p];
double Perfect_array[p];
bool average = true;
bool RMS = false;
bool mode = false;
bool MR = true;

//Make these not hard coded?
double M_delta[p] = {83.3333333333, 150, 210, 321.4285714286, 480, 583.3333333333, 736.6666666667, 888.8888888889, 1090.9090909091, 1242, 1392};
double gluino_mass[p] = {100, 200, 250, 350, 500, 600, 750, 900, 1100, 1250, 1400};

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 700);
    c1->cd();
    c1->SetGrid();
    gStyle->SetOptStat(000000);
    gStyle->SetPalette(51, 0);
    gStyle->SetPaintTextFormat(".2f");
    //gStyle->SetTitleFont(132,"xyz");
    //gStyle->SetTextFont(132);
    //gStyle->SetTextSize(0.05);
    //gStyle->SetStatFont(132);
    //gStyle->SetLabelFont(132,"xyz");
    
    gStyle->SetTitleFont(132,"xyz"); 
    gStyle->SetTitleFont(132," ");   
    gStyle->SetTitleSize(0.06,"xyz"); 
    gStyle->SetTitleSize(0.06,"");
    gStyle->SetTitleFontSize(0.08);
    gStyle->SetLabelFont(132,"xyz");
    gStyle->SetLabelSize(0.02,"xyz");
    gStyle->SetLabelColor(1,"xyz");
    gStyle->SetTextFont(132);
    gStyle->SetTextSize(0.08);
    gStyle->SetStatFont(132);
  
using namespace std;
int mean(){
  for(int file = 0; file < p; file ++){
      if (file == 0){const char* output_file = "T1_150_100.root";}
      if (file == 1){const char* output_file = "T1_200_100.root";}
      if (file == 2){const char* output_file = "T1_250_100.root";}
      if (file == 3){const char* output_file = "T1_350_100.root";}
      if (file == 4){const char* output_file = "T1_500_100.root";}
      if (file == 5){const char* output_file = "T1_600_100.root";}
      if (file == 6){const char* output_file = "T1_750_100.root";}
      if (file == 7){const char* output_file = "T1_900_100.root";}
      if (file == 8){const char* output_file = "T1_1100_100.root";}
      if (file == 9){const char* output_file = "T1_1250_100.root";}
      if (file == 10){const char* output_file = "T1_1400_100.root";}

      TFile *New = new TFile(output_file,"Open");

      if (MR){
      const int numMRBins = 50;
      double mrBinRange[2] = {0, 1500 + file * 250};
      double width = mrBinRange[1] / float(numMRBins);
      }

      else{
      const int numMRBins = 240 + file * 8;
      double mrBinRange[2] = {0, 1.2 + file * 0.04};
      }

      TTree *inTree;
      TH1D *hist_Default = new TH1D("hist_Default", "T1 MR", numMRBins, mrBinRange[0], mrBinRange[1]);
      TH1D *hist_Georgi = new TH1D("hist_Georgi", "T1 MR", numMRBins, mrBinRange[0], mrBinRange[1]);
      TH1D *hist_incrCAM = new TH1D("hist_incrCAM", "T1 MR", numMRBins, mrBinRange[0], mrBinRange[1]);
      TH1D *hist_gen_level = new TH1D("hist_gen_level", "T1 MR", numMRBins, mrBinRange[0], mrBinRange[1]);
      TH1D *hist_KT = new TH1D("hist_KT", "T1 MR", numMRBins, mrBinRange[0], mrBinRange[1]);
    
      Double_t MR_Default, RSQ_Default;
      Double_t MR_newjetcam, RSQ_newjetcam;
      Double_t MR_Georgi, RSQ_Georgi;
      Double_t MR_perfect_cluster, RSQ_perfect_cluster;
      Double_t MR_KT, RSQ_KT;
      
      TFile outFile("jets_study.root", "recreate");
      New->cd();
      inTree = (TTree*)New->Get("RazorInclusive");
      inTree->SetBranchAddress("MR_Default", &MR_Default);
      inTree->SetBranchAddress("RSQ_Default", &RSQ_Default);
      inTree->SetBranchAddress("MR_newjetcam", &MR_newjetcam);
      inTree->SetBranchAddress("RSQ_newjetcam", &RSQ_newjetcam);
      inTree->SetBranchAddress("MR_Georgi", &MR_Georgi);
      inTree->SetBranchAddress("RSQ_Georgi", &RSQ_Georgi);
      inTree->SetBranchAddress("MR_perfect_cluster", &MR_perfect_cluster);
      inTree->SetBranchAddress("RSQ_perfect_cluster", &RSQ_perfect_cluster);
      inTree->SetBranchAddress("MR_KT_Jets", &MR_KT);
      inTree->SetBranchAddress("RSQ_KT_Jets", &RSQ_KT);

      //loop over events
      int nEvents = inTree->GetEntries();
      for(int iEntry = 0; iEntry < nEvents; iEntry++){
	inTree->GetEntry(iEntry);
	if (MR){
	hist_Default->Fill(MR_Default);
	hist_gen_level->Fill(MR_perfect_cluster);
	hist_incrCAM->Fill(MR_newjetcam);
	hist_Georgi->Fill(MR_Georgi);
	hist_KT->Fill(MR_KT);
	}
	else{
	hist_Default->Fill(RSQ_Default);
	hist_gen_level->Fill(RSQ_perfect_cluster);
	hist_incrCAM->Fill(RSQ_newjetcam);
	hist_Georgi->Fill(RSQ_Georgi);
	hist_KT->Fill(RSQ_KT);
	}	    	    
      }

      double GetMedian(TH1D* hist){
	int k = 0;
	//find median.  Currently in units of bins.
	double integral, integral_1;
	double entries = hist->GetEntries();
	for (int i = 0; integral < 0.5; i++){
	  integral = hist->Integral(0, i) / entries;
	  integral_1 = hist->Integral(0, i-1) / entries;
	  k = i;
	}
	//approximation of median with bisection.  
	return 2*(1 - (fabs(0.5 - integral)) * k  + (1-fabs(0.5 - integral_1)) * (k-1))/ fabs(2.0-fabs(0.5 - integral) - fabs(0.5- integral_1));
      }

      int GetMode(TH1D* hist){
	//Currently in units of bin.
	int max = 0;
	int location = 0;
	double entries = hist->GetEntries();
	for (int i = 0; i < 200; i++){
	  if (max < hist->GetBinContent(i)){
	    max = hist->GetBinContent(i);
	    location = i;
	  }
	}
	return location;
      }
      string output = string(output_file);
      size_t pos1 = output.find("_100");
      size_t pos2 = output.find("T1_");
      const string T1_mass = output.substr(pos2 + 3, pos1 - pos2 - 3); 

      if (average){
	MJ_array[file] = hist_Default->GetMean();
	Georgi_array[file] = hist_Georgi->GetMean();
	New_array[file] = hist_incrCAM->GetMean();
	Perfect_array[file] = hist_gen_level->GetMean();
	KT_array[file] = hist_KT->GetMean();
	cout << MJ_array[file] << endl;
      } 
      if (RMS){
	MJ_array[file] = hist_Default->GetRMS();
	Georgi_array[file] = hist_Georgi->GetRMS();
	New_array[file] = hist_incrCAM->GetRMS();
	Perfect_array[file] = hist_gen_level->GetRMS();
	KT_array[file] = hist_KT->GetRMS();
	cout << MJ_array[file] << endl;
  }

      if (mode){
	MJ_array[file] = GetMode(hist_Default) * width;
	Georgi_array[file] = GetMode(hist_Georgi) * width;
	New_array[file] = GetMode(hist_incrCAM) * width;
	Perfect_array[file] = GetMode(hist_gen_level) * width;
	KT_array[file] = GetMode(hist_KT) * width;
	cout << MJ_array[file] << " " << width << endl;
      }
}
  TMultiGraph *mg = new TMultiGraph();
   
  TGraph *gr = new TGraph(p, gluino_mass, MJ_array);
  gr->GetXaxis()->SetTitleOffset(1);
  gr->GetXaxis()->SetTitle("Gluino mass");
  gr->GetYaxis()->SetTitle("T1 Mean");
  gr->GetYaxis()->SetTitleOffset(1.5);
  gr->SetTitle("T1 Mean");
  gr->SetLineWidth(3);
  gr->GetYaxis()->SetRangeUser(0., 2000.);
  //Used for margins??
  gr->Draw("");
  gr->SetLineColor(kPink+10);
  gr->GetXaxis()->SetRange(0, 1300);
  
  TGraph *gr2 = new TGraph(p, gluino_mass, Georgi_array);
  gr2->GetXaxis()->SetTitleOffset(1);
  gr2->SetLineWidth(3);
  gr2->SetLineColor(kViolet-1);
  
    
  TGraph *gr3 = new TGraph(p, gluino_mass, New_array);
  gr3->GetXaxis()->SetTitleOffset(1);
  gr3->SetLineWidth(3);
  gr3->SetLineColor(kCyan+1);
  
//change line colors?!
  TGraph *gr4 = new TGraph(p, gluino_mass, Perfect_array);
  gr4->GetXaxis()->SetTitleOffset(1);
  gr4->SetLineWidth(3);
  gr4->SetLineColor(kBlue+1);
    
  
  TGraph *gr5 = new TGraph(p, gluino_mass, KT_array);
  gr5->GetXaxis()->SetTitleOffset(1);
  gr5->SetLineWidth(3);
  gr5->SetLineColor(kRed+1);
    
  TGraph *gr6 = new TGraph(p, gluino_mass, M_delta);
  gr6->GetXaxis()->SetTitleOffset(1);
  gr6->SetLineWidth(3);
  gr6->SetLineColor(kGreen);
  
  //argument for legend is edge for (left, bottom, right, top)
  TLegend* legend = new TLegend(0.5, 0.14, 0.9, 0.3);
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(1);
  legend->SetTextFont(132);
  legend->SetTextSize(0.03);
  
  legend->AddEntry(gr, "MegaJet", "l");
  legend->AddEntry(gr2, "Georgi_{n=1}", "l");
  legend->AddEntry(gr3, "Incremental R Cam", "l");
  legend->AddEntry(gr4, "Gen-Level", "l");
  legend->AddEntry(gr5, "KT 1.57", "l");
  legend->AddEntry(gr6, "M_#Delta", "l");
  
  legend->Draw();
  
  mg->Add(gr, "l");
  mg->Add(gr2, "l");
  mg->Add(gr3, "l");
  mg->Add(gr4, "l");
  mg->Add(gr5, "");
  mg->Add(gr6, "");
  mg->Draw("l");
  

  c1->Print("../finalpaperplots/T1_Mean.png");
  c1->Print("../finalpaperplots/T1_Mean.pdf");
  c1->Print("../finalpaperplots/T1_Mean.root");
  c1->Print("../finalpaperplots/T1_Mean.C");
  delete c1;
  return 0;
	}
