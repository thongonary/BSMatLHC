#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "RooRandom.h"

using namespace std;

double smearPt(double pT);

int main(){

    RooRandom::randomGenerator()->SetSeed(314159);

    vector<TFile*> inFiles;
    inFiles.push_back(new TFile("t2bb_900_50.root"));
    inFiles.push_back(new TFile("t2bb_900_250.root"));
    inFiles.push_back(new TFile("t2bb_900_500.root"));
    inFiles.push_back(new TFile("t2bb_900_800.root"));

    const int numMRBins = 10;
    const int numR2Bins = 10;
    double mrBinLowEdges[numMRBins+1] = {0.0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500};
    double r2BinLowEdges[numR2Bins+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    TTree *inTree;

    double MR, R2, HT, MHT, DeltaPhi, MET;
    double jetPt[4];
    double jetEta[4];
    double jetPhi[4];
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
        
        TH2D *denominator = new TH2D(Form("file%ddenominator", i), "dijet events with pT > 40", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator1 = new TH2D(Form("file%dnumerator1", i), "pT > 100", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator2 = new TH2D(Form("file%dnumerator2", i), "pT > 60, DeltaPhi < 150 degrees", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator3 = new TH2D(Form("file%dnumerator3", i), "pT > 60, DeltaPhi < 140 degrees", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator4 = new TH2D(Form("file%dnumerator4", i), "pT > 60, DeltaPhi < 130 degrees", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);
        TH2D *numerator5 = new TH2D(Form("file%dnumerator5", i), "pT > 60, MET > 150 GeV", numMRBins, mrBinLowEdges, numR2Bins, r2BinLowEdges);

        //loop over events
        int nEvents = inTree->GetEntries();
        for(int iEntry = 0; iEntry < nEvents; iEntry++){
            inTree->GetEntry(iEntry);

            //fill denominator histogram using unsmeared pT
            if(numJets >= 2 && numJets < 4 && jetPt[0] > 40 && jetPt[1] > 40){ //event is in dijet box
                denominator->Fill(MR, R2);

                double smearedPt1 = smearPt(jetPt[0]);
                double smearedPt2 = smearPt(jetPt[1]);

                if(smearedPt1 > 100 && smearedPt2 > 100) numerator1->Fill(MR, R2);
                if(smearedPt1 > 60 && smearedPt2 > 60){
                    if(DeltaPhi < 2.62) numerator2->Fill(MR, R2);
                    if(DeltaPhi < 2.44) numerator3->Fill(MR, R2);
                    if(DeltaPhi < 2.27) numerator4->Fill(MR, R2);
                    if(MET > 150) numerator5->Fill(MR, R2);
                }
            }
            /*else{
                cout << "Failed dijet cut: numJets = " << numJets << "; jetPt[0] = " << jetPt[0] << "; jetPt[1] = " << jetPt[1] << endl;
            }*/
        }
        outFile.cd();
        denominator->Write();
        numerator1->Write();
        numerator2->Write();
        numerator3->Write();
        numerator4->Write();
        numerator5->Write();

        delete denominator;
        delete numerator1;
        delete numerator2;
        delete numerator3;
        delete numerator4;
        delete numerator5;
    }

    for(int i = 0; i < inFiles.size(); i++) delete inFiles[i];
    return 0;
}

double smearPt(double pT){
    double percentRes = 0.37474*exp(-0.0123776*pT); //this is an exponential fit to the points (36, 0.24) and (92, 0.12)
    return pT*pow(1+percentRes, RooRandom::gaussian());
}
