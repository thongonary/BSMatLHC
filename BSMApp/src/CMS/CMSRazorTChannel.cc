//-------------------------------------------------------
// Description:
// Runs Hgg analysis
// Authors: Duarte, Pena and Wang
//-------------------------------------------------------


#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TTree.h>
#include "CMS/CMSRazorTChannel.hh"
#include <fastjet/tools/Pruner.hh>



#define _debug 0

CMSRazorTChannel::CMSRazorTChannel(TTree *tree, double Lumi, double filterEff, double xsecMax, string analysis, bool delphesFormat) : CMSReco(tree, delphesFormat) {
  _Lumi = Lumi;
  _statTools = new StatTools(-99);
  _analysis = analysis;
  _delphesFormat = delphesFormat;
}

CMSRazorTChannel::~CMSRazorTChannel(){
  std::cout << "CMSRazorTChannel destructor" << std::endl;
}

void CMSRazorTChannel::SetSqrts(double sqrts) {
  _sqrts = sqrts;
}

// loop over events - real analysis
void CMSRazorTChannel::Loop(string outFileName) {
    srand(time(0)); 
    if( _delphesFormat && DelphesTree::fChain == 0) return;
    if (!_delphesFormat && DetectorBase::fChain == 0) return;

    std::cout << "[INFO]: starting Loop" << std::endl;

    //----------------------
    //keep track of cut flow
    //----------------------
    int _noPhotonCut  = 0;
    int _twoPhotonCut = 0;
    int _higgsPtCut   = 0;
    int _twoJetCut    = 0;

    double MR, RSQ;

    double jetPt[50];
    double jetEta[50];
    double jetPhi[50];

    double pho1Pt, pho2Pt;
    double pho1Eta, pho2Eta;
    double pho1Phi, pho2Phi;

    double bjet1Pt = -1, bjet2Pt = -1;
    double bjet1Eta = 999., bjet2Eta = 999.;
    double bjet1Phi = -999., bjet2Phi = -999.;

    int numJets;
    int numJetsAbove80GeV;
    double MET, HT;

    //new variables
    int numBJets;
    int numPhotons;
    int numBox = -1;
    int jetLen = 100;
    int Hem1_csts[jetLen], Hem2_csts[jetLen];
    double dPhi;
    int leadingHemContents, subleadingHemContents;
    int numGenEvents;
    double pthat = -1;


    // Open Output file

    std::cout << "[INFO]: open output file" << std::endl;
    TFile *file = new TFile(outFileName.c_str(),"UPDATE");

    TTree* outTree = new TTree("RazorInclusive","RazorInclusive");
    outTree->Branch("MR", &MR, "MR/D");
    outTree->Branch("RSQ", &RSQ, "RSQ/D");
    outTree->Branch("numJets", &numJets, "numJets/I");
    outTree->Branch("numJetsAbove80GeV", &numJetsAbove80GeV, "numJetsAbove80GeV/I");
    outTree->Branch("numPhotons", &numPhotons, "numPhotons/I");
    outTree->Branch("jetPt", jetPt, "jetPt[numJets]/D");
    outTree->Branch("jetEta", jetEta, "jetEta[numJets]/D");
    outTree->Branch("jetPhi", jetPhi, "jetPhi[numJets]/D");
    outTree->Branch("numBox", &numBox, "numBox/I");
    outTree->Branch("MET", &MET, "MET/D");
    outTree->Branch("HT", &HT, "HT/D");
    // added branches
    outTree->Branch("dPhi", &dPhi, "dPhi/D");
    outTree->Branch("numBJets", &numBJets, "numBJets/I");
    outTree->Branch("pthat", &pthat, "pthat/D");

    double xedge[5] = {150, 250, 400, 1400, 3000};
    double yedge[10] = {0, 0.05, 0.15, 0.25, 0.40, 0.55, 0.70, 0.90, 1.2, 2.0};

    TH2D* pdfHighResRsqMR = new TH2D("pdfHighResRsqMR","pdfHighResRsqMR",4,xedge,9,yedge);


    std::cout << "[INFO]: getting number of entries" << std::endl;

    // loop over entries
    Long64_t nbytes = 0, nb = 0;
    Long64_t nentries = 0;
    if (_delphesFormat) nentries = DelphesTree::fChain->GetEntries();
    else nentries = DetectorBase::fChain->GetEntries();
    std::cout << "[INFO]: Number of entries = " << nentries << std::endl;

    // set the by-event weight
    for (Long64_t jentry=0; jentry<nentries;jentry+=1) {
        if ( verbose ) std::cout << "[VERBOSE]: new event" << std::endl;
        // clean physics-objects blocks
        CleanEvent();
        // get new event
        if ( _delphesFormat ) 
        {
            Long64_t ientry = DelphesTree::LoadTree(jentry);      
            if (ientry < 0) break;
            nb = DelphesTree::fChain->GetEntry(jentry); nbytes += nb;
        }
        else {
            Long64_t ientry = DetectorBase::LoadTree(jentry);
            nb = DetectorBase::fChain->GetEntry(jentry); nbytes += nb;
        }

        if (jentry%100 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

        //reset jet variables
        numJets = 0;
        numJetsAbove80GeV = 0;
        numBJets = 0;
        MET = 0.; 
        HT = 0.;
        for(int i = 0; i < 50; i++){
            jetPt[i] = -1;
            jetEta[i] = -999;
            jetPhi[i] = -999;
        }

        // AK4 jets    
        vector<fastjet::PseudoJet> pfAK04;
        vector<fastjet::PseudoJet> pfAK04_btag;

        std::vector<TLorentzVector> pfJets;
        std::vector<TLorentzVector> pfJetsBtag;

        if ( _delphesFormat ) 
        {
            for (int iJet = 0; iJet < Jet_size; iJet++ )
            {	
                TLorentzVector jet;
                jet.SetPtEtaPhiM(Jet_PT[iJet],Jet_Eta[iJet],Jet_Phi[iJet],Jet_Mass[iJet]);
                if ( jet.Pt() < 40. || jet.Eta() > 3. ) continue; 
                numJets++;
                if ( jet.Pt() > 80.) numJetsAbove80GeV++;
                pfAK04.push_back(fastjet::PseudoJet(jet.Px(), jet.Py(), jet.Pz(), jet.E()));
                pfJets.push_back( jet );
                if ( Jet_BTag[iJet] & (1 << 0) ) 
                {
                    numBJets++; 
                    pfAK04_btag.push_back(fastjet::PseudoJet(jet.Px(), jet.Py(), jet.Pz(), jet.E()));
                    pfJetsBtag.push_back( jet );
                }
            }
        }
        else {
            vector<fastjet::PseudoJet> empty;
            vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty); //note that this I think includes photons
            fastjet::JetDefinition AK04_def(fastjet::antikt_algorithm, 0.4);
            fastjet::ClusterSequence pfAK04ClusterSequence = JetMaker(JetsConst, AK04_def);
            pfAK04 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK04ClusterSequence.inclusive_jets()),40., 3.0); //only cluster jets with > 40 GeV, eta < 3.0
        }

        //---------------------------------------
        //I n p u t   t o   H e m i s p h e r e s
        //---------------------------------------
        std::vector<TLorentzVector> Jets4Vector;
        if( _delphesFormat )
        {
            //-----------------------------------
            //D e l p h e s  j e t s + h i g  g s
            //-----------------------------------
            int i = 0;
            for (auto tmp : pfJets)
            {
                jetPt[i]  = tmp.Pt();
                jetEta[i] = tmp.Eta();
                jetPhi[i] = tmp.Phi();
                HT += tmp.Pt();
                Jets4Vector.push_back(tmp);
                i++;
            }

            //---
            //MET
            //---
            double px = MissingET_MET[0]*cos( MissingET_Phi[0] );
            double py = MissingET_MET[0]*sin( MissingET_Phi[0] );
            PFMET = fastjet::PseudoJet(px, py, 0., MissingET_MET[0] );
            MET = PFMET.pt();
        }
        else
        {
            //-----------------------------------
            //P y t h i a  j e t s + h i g  g s
            //-----------------------------------
            for (int i = 0; i < pfAK04.size(); i++)
            {
                jetPt[i]  = pfAK04[i].pt();
                jetEta[i] = pfAK04[i].eta();
                jetPhi[i] = pfAK04[i].phi();
                HT += pfAK04[i].pt();
                Jets4Vector.push_back(ConvertTo4Vector(pfAK04[i]));
                i++;
            }
            //---
            //MET
            //---
            GenMET();
            PFMET = genMET;
            MET = PFMET.pt();
        }

        //--------------------------
        // At least 2 jets above 80 GeV
        //--------------------------
        if( numJetsAbove80GeV < 2 || Jets4Vector.size() < 2 )
        {
            if( _debug ) std::cout << "[DEBUG]: Not enough objects for razor computation" << std::endl;
            continue;
        }

        if ( _debug ) std::cout << "--> before CMSHEM" << std::endl;
        CMSHemisphere* myHem = new CMSHemisphere(Jets4Vector);
        myHem->CombineMinMass();
        std::vector<TLorentzVector> hem = myHem->GetHemispheres();
        std::vector<int> Temporary = myHem->GetHem1Constituents();
        std::vector<int> Temporary2 = myHem->GetHem2Constituents();
        if ( _debug ) std::cout << "--> pass CMSHEM" << std::endl;

        // SAVES HEM CONSTITUENTS
        int count1 = 0;
        int count2 = 0;
        for (int k=0; k < Temporary.size(); k++){
            if (Temporary[k]>-1) { 
                Hem1_csts[count1] = Temporary[k];
                count1++;
            }
        }
        for (int k=0; k < Temporary2.size(); k++){
            if (Temporary2[k]>-1) {
                Hem2_csts[count2] = Temporary2[k];
                count2++;
            }
        }

        delete myHem;

        if ( _debug ) std::cout << "delete Hem" << std::endl;
        //compute traditional RSQ and MR
        TLorentzVector j1 = hem[0];
        TLorentzVector j2 = hem[1];  
        dPhi = DeltaPhi(j1,j2); //deltaPhi
        MR = CalcMR(j1, j2);
        RSQ = pow(CalcMRT(j1, j2, PFMET),2.)/MR/MR;
        if ( _debug ) std::cout << "[DEBUG]: numBox-> "<< numBox << std::endl;
        pdfHighResRsqMR->Fill(MR,RSQ);
        outTree->Fill();

        char outname[256];
        sprintf(outname,"data/%s.root", _analysis.c_str());

    }
    // Open Output file again 
    file->cd();
    pdfHighResRsqMR->Write();
    outTree->Write();
    file->Close();
}

double CMSRazorTChannel::DeltaPhi(TLorentzVector jet1, TLorentzVector jet2) {
  double deltaPhi = jet1.Phi() - jet2.Phi();
  while (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  while (deltaPhi <= -M_PI) deltaPhi += 2*M_PI;
  return deltaPhi;
}


TH1D* CMSRazorTChannel::XsecProb(TH1D* sigPdf, double eff, string Filename, string directory, int ibin, double xmin, double xmax, bool expected, bool doubleErr) {
  
  int ibinX = sigPdf->GetXaxis()->GetNbins();
  int ibinY = sigPdf->GetYaxis()->GetNbins();
  
  char histname[256];
  if (doubleErr)
    {      
    if (expected) sprintf(histname, "%s%s%s","probVecExp", directory.c_str(), "2xErr");
    else sprintf(histname, "%s%s%s","probVec", directory.c_str(),"2xErr");
    }
  else {
    if (expected) sprintf(histname, "%s%s","probVecExp", directory.c_str());
    else sprintf(histname, "%s%s","probVec", directory.c_str());
  }
  TH1D* probVec = new TH1D(histname, histname, ibin, xmin, xmax);
  
  TFile* likFile = new TFile(TString(Filename));  
  gROOT->cd();
  // a loop over xsec should go here... 
  for(int i=0; i<ibin; i++) {
    //center of bin
    //double xsec = xmin + (i+0.5)/ibin*(xmax-xmin);
    //left edge of bin
    double xsec = xmin + (1.*i)/ibin*(xmax-xmin);
    double prob = 1;
    double logProb = 0;
    for(int ix=0; ix<ibinX; ix++) {
      for(int iy=0; iy<ibinY; iy++) {
	if (sigPdf->GetBinContent(ix+1,iy+1)<=0.) continue;
	double sBin = _Lumi*xsec*eff*sigPdf->GetBinContent(ix+1,iy+1);	
	char name[256];
	if (expected) sprintf(name, "%s/exp_%i_%i", directory.c_str(), ix, iy);
	else sprintf(name, "%s/lik_%i_%i", directory.c_str(), ix, iy);
	
	TH1D* binProb = (TH1D*) likFile->Get(name);
	double yieldmax = binProb->GetXaxis()->GetXmax();
	if (sBin >= yieldmax) {
	  cout << "Note: signal events exceed " << yieldmax << "! sBin = " << sBin << endl;
	  cout << "FindBin returns: " << binProb->FindBin(sBin) << endl;
	  cout << "bin prob returns: " << binProb->GetBinContent(binProb->FindBin(sBin)) << endl;
	}
	logProb += TMath::Log(binProb->GetBinContent(binProb->FindBin(sBin)));
	prob *= binProb->GetBinContent(binProb->FindBin(sBin)); 
	delete binProb;
      }
    }
    probVec->SetBinContent(i+1,TMath::Exp(logProb));
  }
  probVec->Scale(1./probVec->Integral());
  likFile->Close();
  return probVec;
}


TH1D* CMSRazorTChannel::LocalBayesFactor(TH1D* probVec) {  
  int ibin = probVec->GetXaxis()->GetNbins();
  double xmin = probVec->GetXaxis()->GetXmin();
  double xmax = probVec->GetXaxis()->GetXmax();  

  string name = probVec->GetName();
  name = name.substr(7,name.size());
  name = "localBayes"+name;  
  
  TH1D* localBayes = new TH1D(name.c_str(), name.c_str(), ibin, xmin, xmax);
  
  gROOT->cd();
  for(int i=0; i<ibin; i++) {
    double prob = 1;
    double logProb = TMath::Log(probVec->GetBinContent(i+1))-TMath::Log(probVec->GetBinContent(1));
    localBayes->SetBinContent(i+1,TMath::Sign(TMath::Sqrt(2.*TMath::Abs(logProb)),logProb));
  }
  return localBayes;
}
