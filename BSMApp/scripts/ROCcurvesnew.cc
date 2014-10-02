// run with root -l .x <file>
//Change input files, title, and output file each time!
//Makes a 2D hyperbolic cut along MR/RSQ.  Currently set up for MR/RSQ default.  
//Loops over all combinations of background + signal
//Currently leaks memory.  Need to delete more pointers?
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
#include "TVectorD.h"
#define _USE_MATH_DEFINES

//augments the hyperbola with these offsets


double MR_offset = 0.;
double RSQ_offset = 0.0;

// p is the array size
const int p;
//p must be greater than the number of events!!
p = 100000;
const int c;
c = 200;
Double_t mrrsqarraymj[c], backmrrsqarraymj[c];
Double_t MR_Default_temp[p], RSQ_Default_temp[p];
Double_t backMR_Default_temp[p], backRSQ_Default_temp[p];

Double_t mrrsqarrayGeorgi[c], backmrrsqarrayGeorgi[c];
Double_t MR_Georgi_temp[p], RSQ_Georgi_temp[p];
Double_t backMR_Georgi_temp[p], backRSQ_Georgi_temp[p];

Double_t mrrsqarraynew[c], backmrrsqarraynew[c];
Double_t MR_newjetcam_temp[p], RSQ_newjetcam_temp[p];
Double_t backMR_newjetcam_temp[p], backRSQ_newjetcam_temp[p];

Double_t mrrsqarrayKT[c], backmrrsqarrayKT[c];
Double_t MR_KT_temp[p], RSQ_KT_temp[p];
Double_t backMR_KT_temp[p], backRSQ_KT_temp[p];

using namespace std;
int ROCcurvesnew(){
for (int file = 2; file < 10; file++){
  if (file == 2 || file == 3){const char* signal_file = "T2_200_100.root";}
  if (file == 4 || file == 5){const char* signal_file = "T2_500_50.root";}
  if (file == 6 || file == 7){const char* signal_file = "T1_500_100.root";}
  if (file == 8 || file == 9){const char* signal_file = "T1_1400_100.root";}
  if (file == 2 || file == 4 || file == 6 || file == 8){const char* background_file = "zjets.root";}
  if (file == 3 || file == 5 || file == 7 || file == 9){const char* background_file = "ttbar.root";}
  std::ostringstream stream, stream1, stream2, stream3, stream4;
  const char* savename;
  string signal, background, output_name, title_signal, title_background, title_string;
  size_t pos_sig, pos_background;
  signal = string(signal_file);
  background = string(background_file);
  pos_sig = signal.find(".root");
  pos_background = background.find(".root");
  output_name =  "../finalpaperplots/ROC/" + signal.substr(0, pos_sig) + background.substr(0, pos_background) + "ROC";
  title_signal = signal.substr(0, pos_sig);
  title_background = background.substr(0, pos_background);
  stream << output_name + ".png";
  stream2 << output_name + ".pdf";
  stream3 << output_name + ".root";
  stream4 << output_name + ".C";
  savename = (stream.str()).c_str();
  cout << savename << " " << (stream.str()).c_str() << endl;
  //title =(stream1.str()).c_str();
  //cout << output_name<< endl;
  title_string = "M_{R}*R^{2}" + title_signal + " with " + title_background;
  stream1 << title_string;
 

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
    
    //signal input file
    TFile *New = new TFile(signal_file,"Open");
    //TFile *New = new TFile("T2_200_100.root","Open");
    //TFile *New = new TFile("T1_500_100.root","Open");
    //TFile *New = new TFile("T1_500_100.root","Open");


    const int numMRBins = 10.0 * p;
    //const int numR2Bins = 50;

    double mrBinRange[2] = {0, 3000};
    //double r2BinRange[2] = {0, 1.5};

    TTree *inTree;
 
    Double_t MR_Default, RSQ_Default, MR_newjetcam, RSQ_newjetcam, MR_Georgi, RSQ_Georgi;
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

    inTree->SetBranchAddress("MR_KT_Jets", &MR_KT);
    inTree->SetBranchAddress("RSQ_KT_Jets", &RSQ_KT);


    //loop over events
    int nEvents = inTree->GetEntries();
    for(int iEntry = 0; iEntry < nEvents; iEntry++){
      inTree->GetEntry(iEntry);
      
      //Adds MR & RSQ into a list so the hyperbolic cut can be made
      MR_Default_temp[iEntry] =  MR_Default + MR_offset;
      RSQ_Default_temp[iEntry] = RSQ_Default + RSQ_offset;

      MR_Georgi_temp[iEntry] = MR_Georgi + MR_offset;
      RSQ_Georgi_temp[iEntry] = RSQ_Georgi + RSQ_offset;

      MR_newjetcam_temp[iEntry] = MR_newjetcam + MR_offset;
      RSQ_newjetcam_temp[iEntry] = RSQ_newjetcam + RSQ_offset;

      MR_KT_temp[iEntry] = MR_KT + MR_offset;
      RSQ_KT_temp[iEntry] = RSQ_KT + RSQ_offset;


    }
    c1->cd();
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.14);
    
    //The cut
    int i;  
    i = 0;
    double cut;
    cut = 1;
    double product1[p], product2[p], product3[p], product4[p];
    int countn = 0;
    while (i < nEvents){
      //Loops through all MR and RSQ and multiplies them.
      product1[i] = (MR_Default_temp[i] * RSQ_Default_temp[i]);
      product2[i] =  (MR_newjetcam_temp[i] * RSQ_newjetcam_temp[i]);
      product3[i] =  (MR_Georgi_temp[i] * RSQ_Georgi_temp[i]);
      product4[i] =  (MR_KT_temp[i] * RSQ_KT_temp[i]);
      i += 1;
    }
    //Cut ranges from MR*RSQ = 1 to 996
    while (cut < 1000){
      i = 0;
      double count1 = 0;
      double count2 = 0;
      double count3 = 0;
      double count4 = 0;
      while (i < nEvents){
	//Checks if the product is greater than the threshold.
	if (product1[i] > cut){
	  count1 += 1;
	}
	if (product2[i] > cut){
	  count2 += 1;
	}
	if (product3[i] > cut){
	  count3 += 1;
	}
	if (product4[i] > cut){
	  count4 += 1;
	}
	i += 1;
      }
      cut += 5;
      //adds normalized counts to an array that is plotted 
      mrrsqarraymj[countn] = count1 / nEvents;
      mrrsqarraynew[countn] = count2 / nEvents;
      mrrsqarrayGeorgi[countn] = count3 / nEvents;
      mrrsqarrayKT[countn] = count4 / nEvents;
      countn += 1;
}    
    //background    
    //TFile *New = new TFile("../output/zjets.root","Open");
    //TFile *New = new TFile("zjets.root","Open");
    TFile *New = new TFile(background_file,"Open");

    Double_t back_MR_Default, back_RSQ_Default;
    Double_t back_MR_Georgi, back_RSQ_Georgi;
    Double_t back_MR_newjetcam, back_RSQ_newjetcam;
    Double_t back_MR_KT, back_RSQ_KT;
    //background
    TFile outFile("jets_study1.root", "recreate");
    New->cd();
inTree = (TTree*)New->Get("RazorInclusive");
    inTree->SetBranchAddress("MR_Default", &back_MR_Default);
    inTree->SetBranchAddress("RSQ_Default", &back_RSQ_Default);

    inTree->SetBranchAddress("MR_Georgi", &back_MR_Georgi);
    inTree->SetBranchAddress("RSQ_Georgi", &back_RSQ_Georgi);

    inTree->SetBranchAddress("MR_newjetcam", &back_MR_newjetcam);
    inTree->SetBranchAddress("RSQ_newjetcam", &back_RSQ_newjetcam);
  
    inTree->SetBranchAddress("MR_KT_Jets", &back_MR_KT);
    inTree->SetBranchAddress("RSQ_KT_Jets", &back_RSQ_KT);
     
    //loop over events
    int nEvents = inTree->GetEntries();
    for(int iEntry = 0; iEntry < nEvents; iEntry++){
      inTree->GetEntry(iEntry);

      backMR_Default_temp[iEntry] = back_MR_Default + MR_offset;
      backRSQ_Default_temp[iEntry] = back_RSQ_Default + RSQ_offset;

      backMR_Georgi_temp[iEntry] = back_MR_Georgi + MR_offset;
      backRSQ_Georgi_temp[iEntry] = back_RSQ_Georgi + RSQ_offset;

      backMR_newjetcam_temp[iEntry] = back_MR_newjetcam + MR_offset;
      backRSQ_newjetcam_temp[iEntry] = back_RSQ_newjetcam + RSQ_offset;

      backMR_KT_temp[iEntry] = back_MR_KT + MR_offset;
      backRSQ_KT_temp[iEntry] = back_RSQ_KT + RSQ_offset;

    }
    c1->cd();
    gPad->SetBottomMargin(0.14);

    int i;  
    i = 0;
    double cut;
    cut = 1;
    double product1[p], product2[p], product3[p], product4[p];
    int countn = 0;
    while (i < nEvents){
      product1[i] = (backMR_Default_temp[i] * backRSQ_Default_temp[i]);
      product2[i] = (backMR_newjetcam_temp[i] * backRSQ_newjetcam_temp[i]);
      product3[i] = (backMR_Georgi_temp[i] * backRSQ_Georgi_temp[i]);
      product4[i] = (backMR_KT_temp[i] * backRSQ_KT_temp[i]);
      i += 1;
    }
    while (cut < 1000){
      i = 0;
      double count1 = 0;
      double count2 = 0;
      double count3 = 0;
      double count4 = 0;
      while (i < nEvents){
	if (product1[i] > cut){
	  count1 += 1;
	}
	if (product2[i] > cut){
	  count2 += 1;
	}
	if (product3[i] > cut){
	  count3 += 1;
	}
	if (product4[i] > cut){
	  count4 += 1;
	}
	i += 1;
      }
      cut += 5;
      //add counts to an array that is plotted 
      //cout << count / nEvents  << endl;
      backmrrsqarraymj[countn] = 1 - (count1 / nEvents);
      backmrrsqarraynew[countn] = 1 - (count2 / nEvents);
      backmrrsqarrayGeorgi[countn] = 1 - (count3 / nEvents);
      backmrrsqarrayKT[countn] = 1 - (count4 / nEvents);
      countn += 1;
    } 

    TMultiGraph *mg = new TMultiGraph();
    
    TGraph *gr = new TGraph(c, backmrrsqarraymj, mrrsqarraymj);
    gr->GetXaxis()->SetTitleOffset(1);
    gr->GetXaxis()->SetTitle("Percent Background Rejected");
    gr->GetYaxis()->SetTitle("Percent Signal Accepted");
    gr->GetYaxis()->SetTitleOffset(0.80);
    gr->SetTitle((stream1.str()).c_str());
    gr->SetLineWidth(3);
    //gr->SetLineColor(kOrange);
    gr->SetMaximum(1.0);
    //Used for margins??
    gr->Draw("");
    gr->SetLineColor(kPink+10);
    //Measured in bins?
    gr->GetXaxis()->SetRange(0, c / 1.08);

    TGraph *gr2 = new TGraph(c, backmrrsqarraynew, mrrsqarraynew);
    gr2->GetXaxis()->SetTitleOffset(1);
    gr2->SetLineWidth(3);
    gr2->SetLineColor(kViolet-1);

    
    TGraph *gr3 = new TGraph(c,backmrrsqarrayGeorgi, mrrsqarrayGeorgi);
    gr3->GetXaxis()->SetTitleOffset(1);
    gr3->SetLineWidth(3);
    gr3->SetLineColor(kCyan+1);
  
    TGraph *gr4 = new TGraph(c,backmrrsqarrayKT, mrrsqarrayKT);
    gr4->GetXaxis()->SetTitleOffset(1);
    gr4->SetLineWidth(3);
    gr4->SetLineColor(kBlue);

    //argument for legend is edge for (left, bottom, right, top)
    TLegend* legend = new TLegend(0.1, 0.14, 0.5, 0.3);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(1);
    legend->SetTextFont(132);
    legend->SetTextSize(0.03);

    legend->AddEntry(gr, "MegaJet", "l");
    legend->AddEntry(gr2, " Incremental R CAM", "l");
    legend->AddEntry(gr3, "Georgi_{n=1}", "l");
    legend->AddEntry(gr4, "KT 1.57", "l");
    legend->Draw();
	
    mg->Add(gr, "");
    mg->Add(gr2, "");
    mg->Add(gr3, "");
    mg->Add(gr4, "");
    mg->Draw("l");
    cout << savename << endl;
    cout << (stream.str()).c_str();
    c1->Print((stream.str()).c_str());
    c1->Print((stream2.str()).c_str());
    c1->Print((stream3.str()).c_str());
    c1->Print((stream4.str()).c_str());
    delete c1;
    delete New;
    delete legend;
 }
}
