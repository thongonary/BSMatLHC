#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "RooRandom.h"
#include "TStyle.h"
#include "TCanvas.h"

using namespace std;

double smearPt(double pT);

int main(){
    gStyle->SetOptStat(000000);
    gStyle->SetPalette(51, 0);

    RooRandom::randomGenerator()->SetSeed(314159);

    vector<TFile*> inFiles;
    inFiles.push_back(new TFile("ttbar.root"));
    inFiles.push_back(new TFile("t2bb_900_50.root"));
    inFiles.push_back(new TFile("t2bb_900_250.root"));
    inFiles.push_back(new TFile("t2bb_900_500.root"));
    inFiles.push_back(new TFile("t2bb_900_800.root"));

    const int numMRBins = 7;
    const int numR2Bins = 6;
    double mrBinLowEdges[numMRBins+1] = {400, 450, 550, 700, 900, 1200, 1600, 2500};
    double r2BinLowEdges[numR2Bins+1] = {0.25, 0.3, 0.41, 0.52, 0.64, 0.8, 1.5};
    TTree *inTree;

    double MR, R2, HT, MHT, DeltaPhi, MET;
    double jetPt[50];
    double jetEta[50];
    double jetPhi[50];
    int numJets;
    TFile outFile("triggerStudy.root", "recreate");
    for(int i = 0; i < inFiles.size(); i++){
        inFiles[i]->cd();
        inTree = (TTree*)inFiles[i]->Get("RazorInclusive");
        inTree->SetBranchAddress("MRNEW", &MR);
        inTree->SetBranchAddress("RSQNEW", &R2);
        inTree->SetBranchAddress("HT", &HT);
        inTree->SetBranchAddress("MHT", &MHT);
        inTree->SetBranchAddress("DeltaPhi", &DeltaPhi);
        inTree->SetBranchAddress("jetPt", jetPt);
        inTree->SetBranchAddress("jetEta", jetEta);
        inTree->SetBranchAddress("jetPhi", jetPhi);
        inTree->SetBranchAddress("MET", &MET);
        inTree->SetBranchAddress("numJets", &numJets);
        
        TH2D *denominator = new TH2D(Form("file%ddenominator", i), Form("dijet events with pT > 40, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator1 = new TH2D(Form("file%defficiency1", i), Form("pT > 100, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator2 = new TH2D(Form("file%defficiency2", i), Form("pT > 60, DeltaPhi < 150 degrees, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator3 = new TH2D(Form("file%defficiency3", i), Form("pT > 60, DeltaPhi < 170 degrees, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator4 = new TH2D(Form("file%defficiency4", i), Form("pT > 60, DeltaPhi < 160 degrees, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator5 = new TH2D(Form("file%defficiency5", i), Form("pT > 60, MHT > 150 GeV, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator6 = new TH2D(Form("file%defficiency6", i), Form("HT > 100, MHT/HT > 0.8, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator7 = new TH2D(Form("file%defficiency7", i), Form("(pT > 60, DeltaPhi < 160 OR MHT > 150, %s", inFiles[i]->GetName()), numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);

        //loop over events
        int nEvents = inTree->GetEntries();
        for(int iEntry = 0; iEntry < nEvents; iEntry++){
            inTree->GetEntry(iEntry);

            //fill denominator histogram using unsmeared pT
            if(numJets >= 2 && numJets < 4 && jetPt[0] > 40 && jetPt[1] > 40){ //event is in dijet box
                denominator->Fill(MR, R2);

                double HTSmeared = 0;
                double pxSmeared = 0; double pySmeared = 0;
                double smearedPt1 = smearPt(jetPt[0]);
                HTSmeared += smearedPt1;
                pxSmeared += smearedPt1*cos(jetPhi[0]);
                pySmeared += smearedPt1*sin(jetPhi[0]);
                double smearedPt2 = smearPt(jetPt[1]);
                HTSmeared += smearedPt2;
                pxSmeared += smearedPt2*cos(jetPhi[1]);
                pySmeared += smearedPt2*sin(jetPhi[1]);
                for(int j = 2; j < numJets; j++){
                    double smeared = smearPt(jetPt[j]);
                    HTSmeared += smeared;
                    pxSmeared += smeared*cos(jetPhi[j]);
                    pySmeared += smeared*sin(jetPhi[j]);
                }
                double MHTSmeared = sqrt(pxSmeared*pxSmeared + pySmeared*pySmeared);

                if(smearedPt1 > 100 && smearedPt2 > 100) numerator1->Fill(MR, R2);
                if(smearedPt1 > 60 && smearedPt2 > 60){
                    if(DeltaPhi < 2.62) numerator2->Fill(MR, R2);//150 degrees
                    if(DeltaPhi < 2.96) numerator3->Fill(MR, R2);//170 degrees
                    if(DeltaPhi < 2.79) numerator4->Fill(MR, R2);//160 degrees
                    if(MHTSmeared > 150) numerator5->Fill(MR, R2);
                    if(DeltaPhi < 2.79 || MHTSmeared > 150) numerator7->Fill(MR, R2);
                }
                if(HTSmeared > 100 && MHTSmeared/HTSmeared > 0.8) numerator6->Fill(MR, R2);
            }
            /*else{
                cout << "Failed dijet cut: numJets = " << numJets << "; jetPt[0] = " << jetPt[0] << "; jetPt[1] = " << jetPt[1] << endl;
            }*/
        }

        TH2D *efficiency1 = (TH2D*)numerator1->Clone();
        TH2D *efficiency2 = (TH2D*)numerator2->Clone();
        TH2D *efficiency3 = (TH2D*)numerator3->Clone();
        TH2D *efficiency4 = (TH2D*)numerator4->Clone();
        TH2D *efficiency5 = (TH2D*)numerator5->Clone();
        TH2D *efficiency6 = (TH2D*)numerator6->Clone();
        TH2D *efficiency7 = (TH2D*)numerator7->Clone();

        efficiency1->Divide(denominator);
        efficiency2->Divide(denominator);
        efficiency3->Divide(denominator);
        efficiency4->Divide(denominator);
        efficiency5->Divide(denominator);
        efficiency6->Divide(denominator);
        efficiency7->Divide(denominator);

        TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
        efficiency1->GetXaxis()->SetTitle("M_{R}");
        efficiency1->GetYaxis()->SetTitle("R^{2}");
        efficiency2->GetXaxis()->SetTitle("M_{R}");
        efficiency2->GetYaxis()->SetTitle("R^{2}");
        efficiency3->GetXaxis()->SetTitle("M_{R}");
        efficiency3->GetYaxis()->SetTitle("R^{2}");
        efficiency4->GetXaxis()->SetTitle("M_{R}");
        efficiency4->GetYaxis()->SetTitle("R^{2}");
        efficiency5->GetXaxis()->SetTitle("M_{R}");
        efficiency5->GetYaxis()->SetTitle("R^{2}");
        efficiency6->GetXaxis()->SetTitle("M_{R}");
        efficiency6->GetYaxis()->SetTitle("R^{2}");
        efficiency7->GetXaxis()->SetTitle("M_{R}");
        efficiency7->GetYaxis()->SetTitle("R^{2}");
        denominator->GetXaxis()->SetTitle("M_{R}");
        denominator->GetYaxis()->SetTitle("R^{2}");
        efficiency1->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency2->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency3->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency4->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency5->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency6->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency7->GetZaxis()->SetRangeUser(0.0, 1);
        efficiency1->SetTitle("");
        efficiency2->SetTitle("");
        efficiency3->SetTitle("");
        efficiency4->SetTitle("");
        efficiency5->SetTitle("");
        efficiency6->SetTitle("");
        efficiency7->SetTitle("");
        denominator->SetTitle("");
        efficiency1->Draw("colz");
        c1->Print(Form("dijetpt100file%d.pdf", i));
        c1->Print(Form("dijetpt100file%d.root", i));
        efficiency2->Draw("colz");
        c1->Print(Form("dijetpt60deltaphi150file%d.pdf", i));
        c1->Print(Form("dijetpt60deltaphi150file%d.root", i));
        efficiency3->Draw("colz");
        c1->Print(Form("dijetpt60deltaphi170file%d.pdf", i));
        c1->Print(Form("dijetpt60deltaphi170file%d.root", i));
        efficiency4->Draw("colz");
        c1->Print(Form("dijetpt60deltaphi160file%d.pdf", i));
        c1->Print(Form("dijetpt60deltaphi160file%d.root", i));
        efficiency5->Draw("colz");
        c1->Print(Form("dijetpt60met150file%d.pdf", i));
        c1->Print(Form("dijetpt60met150file%d.root", i));
        efficiency6->Draw("colz");
        c1->Print(Form("dijetht100mhtht0p8file%d.pdf", i));
        c1->Print(Form("dijetht100mhtht0p8file%d.root", i));
        efficiency7->Draw("colz");
        c1->Print(Form("dijetpt60deltaphi160ORmht150file%d.root", i));
        c1->Print(Form("dijetpt60deltaphi160ORmht150file%d.pdf", i));
        denominator->Draw("colz");
        c1->Print(Form("dijetdenominatorfile%d.root", i));
        c1->Print(Form("dijetdenominatorfile%d.pdf", i));

        delete c1;

        outFile.cd();
        denominator->Write();
        efficiency1->Write();
        efficiency2->Write();
        efficiency3->Write();
        efficiency4->Write();
        efficiency5->Write();
        efficiency6->Write();
        efficiency7->Write();

        delete denominator;
        delete numerator1;
        delete numerator2;
        delete numerator3;
        delete numerator4;
        delete numerator5;
        delete numerator6;
        delete numerator7;
    }

    for(int i = 0; i < inFiles.size(); i++) delete inFiles[i];
    return 0;
}

double smearPt(double pT){
    double percentRes = 0.37474*exp(-0.0123776*pT); //this is an exponential fit to the points (36, 0.24) and (92, 0.12)
    return pT*pow(1+percentRes, RooRandom::gaussian());
}
