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
#include "CMS/CMSRazorHgg.hh"
#include <fastjet/tools/Pruner.hh>


#define _debug 0

CMSRazorHgg::CMSRazorHgg(TTree *tree, double Lumi, double filterEff, string analysis, bool delphesFormat) : CMSReco(tree, delphesFormat) {
  _Lumi = Lumi;
  _filterEff = filterEff;
  _statTools = new StatTools(-99);
  _analysis = analysis;
  _delphesFormat = delphesFormat;
}

CMSRazorHgg::~CMSRazorHgg(){
  std::cout << "CMSRazorHgg destructor" << std::endl;
}

void CMSRazorHgg::SetSqrts(double sqrts) {
  _sqrts = sqrts;
}

// loop over events - real analysis
void CMSRazorHgg::Loop(string outFileName) {
  srand(time(0)); 
  if( _delphesFormat && DelphesTree::fChain == 0) return;
  if (!_delphesFormat && DetectorBase::fChain == 0) return;
  
  std::cout << "[INFO]: starting Loop" << std::endl;
  
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

  double higgsPt, higgsEta, higgsPhi, higgsMass;
  int numJets;
  double MET;
  
  //new variables
  int numBJets;
  int numPhotons;
  int numBox = -1;
  int jetLen = 100;
  int Hem1_csts[jetLen], Hem2_csts[jetLen];
  double dPhi, dPhiHiggsJ1, dPhiHiggsJ2;
  int leadingHemContents, subleadingHemContents;
  int numGenEvents;
  double pthat = -1;
  int numSbottoms;
  int numHgg=0;
  double neutMET;
  double sbottomPt[2];
  double sbottomEta[2];
  double sbottomPhi[2];
  double lspPt[2];
  double lspEta[2];
  double lspPhi[2];


  // Open Output file
  
  std::cout << "[INFO]: open output file" << std::endl;
  TFile *file = new TFile(outFileName.c_str(),"UPDATE");

  TTree* outTree = new TTree("RazorInclusive","RazorInclusive");
  outTree->Branch("MR", &MR, "MR/D");
  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("numJets", &numJets, "numJets/I");
  outTree->Branch("numPhotons", &numPhotons, "numPhotons/I");
  outTree->Branch("jetPt", jetPt, "jetPt[numJets]/D");
  outTree->Branch("jetEta", jetEta, "jetEta[numJets]/D");
  outTree->Branch("jetPhi", jetPhi, "jetPhi[numJets]/D");

  outTree->Branch("pho1Pt", &pho1Pt, "pho1Pt/D");
  outTree->Branch("pho1Eta", &pho1Eta, "pho1Eta/D");
  outTree->Branch("pho1Phi", &pho1Phi, "pho1Phi/D");
  outTree->Branch("pho2Pt", &pho2Pt, "pho2Pt/D");
  outTree->Branch("pho2Eta", &pho2Eta, "pho2Eta/D");
  outTree->Branch("pho2Phi", &pho2Phi, "pho2Phi/D");
  outTree->Branch("higgsPt", &higgsPt, "higgsPt/D");
  outTree->Branch("higgsEta", &higgsEta, "higgsEta/D");
  outTree->Branch("higgsPhi", &higgsPhi, "higgsPhi/D");
  outTree->Branch("higgsMass", &higgsMass, "higgsMass/D");

  outTree->Branch("numBox", &numBox, "numBox/I");
  outTree->Branch("MET", &MET, "MET/D");

  // added branches
  outTree->Branch("dPhi", &dPhi, "dPhi/D");
  outTree->Branch("dPhiHiggsJ1", &dPhiHiggsJ1, "dPhiHiggsJ1/D");
  outTree->Branch("dPhiHiggsJ2", &dPhiHiggsJ2, "dPhiHiggsJ2/D");
  outTree->Branch("numBJets", &numBJets, "numBJets/I");
  outTree->Branch("leadingHemContents", &leadingHemContents, "leadingHemContents/I");
  outTree->Branch("subleadingHemContents", &subleadingHemContents, "subleadingHemContents/I"); //0 if higgs only, //1 if jets only, //2 if higgs and jets
  outTree->Branch("pthat", &pthat, "pthat/D");
  outTree->Branch("sbottomPt", &sbottomPt, "sbottomPt[2]/D");
  outTree->Branch("sbottomEta", &sbottomEta, "sbottomEta[2]/D");
  outTree->Branch("sbottomPhi", &sbottomPhi, "sbottomPhi[2]/D");
  outTree->Branch("numSbottoms", &numSbottoms, "numSbottoms/I");
  outTree->Branch("numHgg", &numHgg, "numHgg/I");
  outTree->Branch("neutMET", &neutMET, "neutMET/D");
  outTree->Branch("lspPt", &lspPt, "lspPt[2]/D");
  outTree->Branch("lspEta", &lspEta, "lspEta[2]/D");
  outTree->Branch("lspPhi", &lspPhi, "lspPhi[2]/D");
  outTree->Branch("pthat", &pthat, "pthat/D");

  
  TH1D* pdfHighPt = new TH1D("pdfHighPt","pdfHighPt",15,0,15);
  TH1D* pdfHbb = new TH1D("pdfHbb","pdfHbb",3,0,3);
  TH1D* pdfZbb = new TH1D("pdfZbb","pdfZbb",3,0,3);
  TH1D* pdfHighRes = new TH1D("pdfHighRes","pdfHighRes",10,0,10);
  TH1D* pdfTotal = new TH1D("pdfTotal","pdfTotal",31,0,31);
  
  double xedge[5] = {150, 250, 400, 1400, 3000};
  double yedge[5] = {0, 0.05, 0.10, 0.15, 1.0};

  TH2D* pdfHighResRsqMR = new TH2D("pdfHighResRsqMR","pdfHighResRsqMR",4,xedge,4,yedge);


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
    numBJets = 0;
    numPhotons = 0;
    for(int i = 0; i < 50; i++){
        jetPt[i] = -1;
        jetEta[i] = -999;
        jetPhi[i] = -999;
    }
    pho1Pt = 0; pho1Eta = -999; pho1Phi = -999;
    pho2Pt = 0; pho2Eta = -999; pho2Phi = -999;
    higgsPt = 0; higgsEta = -999; higgsPhi = -999;
    higgsMass = -999;
    for(int i = 0; i < 2; i++){
      sbottomPt[i] = -999;
      sbottomEta[i] = -999;
      sbottomPhi[i] = -999;
      lspPt[i] =  -999;
      lspEta[i] = -999;
      lspPhi[i] = -999;
    }


    TLorentzVector sbottomvector1;
    TLorentzVector sbottomvector2;
    numSbottoms = 0;
    numHgg = 0;

    int numLSP = 0;
    TLorentzVector LSP1;
    TLorentzVector LSP2;

    // THIS PART SAVES SBOTTOM + LSP INFORMATION
    if (_delphesFormat) {      
      for(int iParticle = 0; iParticle < Particle_size; iParticle++){
	if (fabs(Particle_PID[iParticle])==1000005) {
	  numSbottoms++;
	  if (numSbottoms == 1){
	    sbottomvector1.SetPxPyPzE(Particle_Px[iParticle], Particle_Py[iParticle], Particle_Pz[iParticle], Particle_E[iParticle]);
	  }
	  else{
	    sbottomvector2.SetPxPyPzE(Particle_Px[iParticle], Particle_Py[iParticle], Particle_Pz[iParticle], Particle_E[iParticle]);
	  }
	}
	if (fabs(Particle_PID[iParticle])==1000022) {
	  numLSP++;
	  if (numLSP == 1){
	    LSP1.SetPxPyPzE(Particle_Px[iParticle], Particle_Py[iParticle], Particle_Pz[iParticle], Particle_E[iParticle]);
	  }
	  else{
	    LSP2.SetPxPyPzE(Particle_Px[iParticle], Particle_Py[iParticle], Particle_Pz[iParticle], Particle_E[iParticle]);
	  }
	}
      }
    }
    else {
      for(int iSUSY = 0; iSUSY < SUSY; iSUSY++){
	if (fabs(SUSYPdgId[iSUSY])==1000005 && (SUSYD1PdgId[iSUSY]==1000023 || SUSYD1PdgId[iSUSY]==1000022)) {
	  numSbottoms++;
	  if (numSbottoms == 1){
	    sbottomvector1.SetPxPyPzE(SUSYPx[iSUSY], SUSYPy[iSUSY], SUSYPz[iSUSY], SUSYE[iSUSY]);
	  }
	  else{
	    sbottomvector2.SetPxPyPzE(SUSYPx[iSUSY], SUSYPy[iSUSY], SUSYPz[iSUSY], SUSYE[iSUSY]);
	  }
	}
	if (fabs(SUSYPdgId[iSUSY])==1000022) {
	  numLSP++;
	  if (numLSP == 1){
	    LSP1.SetPxPyPzE(SUSYPx[iSUSY], SUSYPy[iSUSY], SUSYPz[iSUSY], SUSYE[iSUSY]);
	  }
	  else{
	    LSP2.SetPxPyPzE(SUSYPx[iSUSY], SUSYPy[iSUSY], SUSYPz[iSUSY], SUSYE[iSUSY]);
	  }
	}
      }
    }
    
    lspPt[0] = LSP1.Pt();
    lspPt[1] = LSP2.Pt();
    lspEta[0] = LSP1.Eta();
    lspEta[1] = LSP2.Eta();
    lspPhi[0] = LSP1.Phi();
    lspPhi[1] = LSP2.Phi();
    
    TLorentzVector LSPs = LSP1 + LSP2;
    neutMET = LSPs.Pt();
    
    TLorentzVector sbottoms = sbottomvector1 + sbottomvector2;
    sbottomPt[0] = sbottomvector1.Pt();
    sbottomPt[1] = sbottomvector2.Pt();
    sbottomEta[0] = sbottomvector1.Eta();
    sbottomEta[1] = sbottomvector2.Eta();
    sbottomPhi[0] = sbottomvector1.Phi();
    sbottomPhi[1] = sbottomvector2.Phi();
    pthat = sbottoms.Pt();
    
    // Build the event 
    PFReco();
    
    //Photons 
    //double bestSumPt = -1;
    int bestPhotIndex1 = -1;
    int bestPhotIndex2 = -1;
    double bestMass = -1;
    int numHiggs = 0;
    double secondBestSumPt = -1;
    int secondBestPhotIndex1 = -1;
    int secondBestPhotIndex2 = -1;
    double secondBestMass = -1; 
    numBox = -1;
    
    //-----------------------------------
    //S e l e c t   a l l   p h o t o n s
    //-----------------------------------
    std::vector<TLorentzVector> phoCand;
    //loop over all pairs of gen photons and find the pair with largest scalar sum Pt
    for ( auto pho : _PFPhotons )
      {
	if( pho.pt() < 25 || fabs( pho.eta() ) > 2.5 ) continue;
	if( fabs( pho.eta() ) > 1.4442 && fabs( pho.eta() ) < 1.566 ) continue;
	TLorentzVector tmp;
	tmp.SetPtEtaPhiM( pho.pt(), pho.eta(), pho.phi(), .0 );
	phoCand.push_back( tmp );
      }
    
    //---------------------------------------
    //find the "best" photon pair, highest Pt!
    //---------------------------------------
    TLorentzVector HiggsCandidate(0,0,0,0);
    double bestSumPt = -99.;
    TLorentzVector photon1, photon2;
    for( size_t i = 0; i < phoCand.size(); i++ )
      {
	for( size_t j = i+1; j < phoCand.size(); j++ )
	  {
	    TLorentzVector pho1 = phoCand[i];
	    TLorentzVector pho2 = phoCand[j];
	    if ( _debug )
	      {
		std::cout << "[DEBUG]: pho1-> " << pho1.Pt() << " [DEBUG]: pho2->" << pho2.Pt() 
			  << std::endl;
	      }
	    //need one photon in the pair to have pt > 40 GeV
	    if ( pho1.Pt() < 40.0 && pho2.Pt() < 40.0 )
	      {
		if ( _debug ) std::cout << "[DEBUG]: both photons failed PT > 40 GeV" << std::endl; 
		continue;
	      }
	    //need diphoton mass between > 50 GeV, cut offline
	    double diphotonMass = (pho1 + pho2).M();
	    if( diphotonMass < 50 )
	      {
		if ( _debug ) std::cout << "[DEBUG]: Diphoton mass < 100 GeV: mgg->" << diphotonMass << std::endl;
		if ( _debug ) std::cout << "... pho1Pt: " << pho1.Pt()  << " pho2Pt: " << pho2.Pt()  << std::endl;
		continue;
	      }

	    if ( _debug )
	      {
		std::cout << "[DEBUG] Diphoton Sum pT: " << pho1.Pt() + pho2.Pt() << std::endl;
	      }
	    //if sum PT of photons is larger than the current Higgs candidate, make this the Higgs candidate
	    if( pho1.Pt() + pho2.Pt() > bestSumPt )
	      {
		bestSumPt = pho1.Pt() + pho2.Pt();
		HiggsCandidate = pho1 + pho2;
		photon1 = pho1;
		photon2 = pho2;
	      }
	  }
      }
    
    /*
    for(int iPhot = 0; iPhot < _PFPhotons.size(); iPhot++){
      if(_PFPhotons[iPhot].pt() < 25 && abs(_PFPhotons[iPhot].eta()) > 1.44) continue; //only count photons that have pt > 25 GeV, eta <1 .44
      numPhotons++;
      for(int jPhot = iPhot+1; jPhot < _PFPhotons.size(); jPhot++){
      if(_PFPhotons[jPhot].pt() < 25 || abs(_PFPhotons[jPhot].eta()) > 1.44) continue; //only count photons that have pt > 25 GeV, eta <1 .44
      if(_PFPhotons[jPhot].pt() < 40 && _PFPhotons[iPhot].pt() < 40) continue; //at least one photon is 40 GeV
      fastjet::PseudoJet higgsCandidate = _PFPhotons[iPhot] + _PFPhotons[jPhot];
      if (higgsCandidate.pt() < 20) continue; //higgs candidate at least 20 GeV
      //check if this pair is in the correct mass range
      if(higgsCandidate.m() > 100){
      cout << "Checking!" << endl;
      numHiggs++;
	  if(_PFPhotons[iPhot].pt() + _PFPhotons[jPhot].pt() > bestSumPt){
	  cout << "Made it!" << endl;
	  //secondBestMass = bestMass;
	    //secondBestSumPt = bestSumPt;
	    //secondBestPhotIndex1 = bestPhotIndex1;
	    //secondBestPhotIndex2 = bestPhotIndex2;
	    bestMass = higgsCandidate.m();
	    bestSumPt = _PFPhotons[iPhot].pt() + _PFPhotons[jPhot].pt();
	    bestPhotIndex1 = iPhot;
	    bestPhotIndex2 = jPhot;
	    if (_PFPhotons[iPhot].pt() < _PFPhotons[jPhot].pt()){	      
	      bestPhotIndex1 = jPhot;
	      bestPhotIndex2 = iPhot;
	    }
	  }
	}
      }
    }
    
    if (numHiggs<1) continue;
    cout << "Found higgs candidate!!!!" << endl;
    */
    if ( HiggsCandidate.Pt() < 20 )
      {
	if ( _debug ) std::cout << "[DEBUG]: failed higgs pt requirement"<< std::endl;
	continue;
      }
        
    if ( _debug ) std::cout << "[DEBUG]: --> pass higgs candidate selection" << std::endl;
    //---------------------
    //fill photon variables
    //---------------------
    pho1Pt  = photon1.Pt();
    pho1Eta = photon1.Eta();
    pho1Phi = photon1.Phi();
    pho2Pt  = photon2.Pt();
    pho2Eta = photon2.Eta();
    pho2Phi = photon2.Phi();
    //---------
    //Higgs
    //---------
    higgsPt     = HiggsCandidate.Pt();
    higgsEta    = HiggsCandidate.Eta();
    higgsPhi    = HiggsCandidate.Phi();
    higgsMass   = HiggsCandidate.M();
    
    /*
    fastjet::PseudoJet pho1(pho1Pt*cos(pho1Phi), pho1Pt*sin(pho1Phi), pho1Pt*sinh(pho1Eta), pho1Pt*cosh(pho1Eta));
    fastjet::PseudoJet pho2(pho2Pt*cos(pho2Phi), pho2Pt*sin(pho2Phi), pho2Pt*sinh(pho2Eta), pho2Pt*cosh(pho2Eta));
    
    fastjet::PseudoJet higgs = pho1 + pho2;
    double higgsPhi = higgs.phi();
    double higgsEta = higgs.eta();
    double higgsPt = higgs.pt();
    double higgsEnergy = higgs.E();
    TLorentzVector higgsvector;
    higgsvector.SetPtEtaPhiE(higgsPt, higgsEta, higgsPhi, higgsEnergy);
    */
    // AK5 jets    
    vector<fastjet::PseudoJet> pfAK05;
    vector<fastjet::PseudoJet> pfAK05_btag;

    std::vector<TLorentzVector> pfJets;
    std::vector<TLorentzVector> pfJetsBtag;
    if ( _delphesFormat ) 
      {
	for (int iJet = 0; iJet < Jet_size; iJet++ )
	  {	
	    TLorentzVector jet;
	    jet.SetPtEtaPhiM(Jet_PT[iJet],Jet_Eta[iJet],Jet_Phi[iJet],Jet_Mass[iJet]);
	    if ( jet.Pt() < 30. || jet.Eta() > 3. ) continue; 
	    //check if within DR < 0.5 of a selected photon
	    double thisDR = min( jet.DeltaR( photon1 ), jet.DeltaR( photon2 ) );	  
	    if ( thisDR > 0.5 )
	      {
		numJets++;
		pfAK05.push_back(fastjet::PseudoJet(jet.Px(), jet.Py(), jet.Pz(), jet.E()));
		pfJets.push_back( jet );
		if ( Jet_BTag[iJet] ) 
		  {
		    numBJets++; 
		    pfAK05_btag.push_back(fastjet::PseudoJet(jet.Px(), jet.Py(), jet.Pz(), jet.E()));
		    pfJetsBtag.push_back( jet );
		  }
	      }
	  }
	/*
	  double minHiggsMassDiff = -1;
	double minZMassDiff = -1;
	int bestBjetIndex1 = -1;
	int bestBjetIndex2 = -1;    
	
	if (numBox!=0) {
	  for (int i = 0; i < pfAK05_btag.size(); i++){
	    for (int j = i+1; j < pfAK05_btag.size(); j++){	
	      fastjet::PseudoJet bbCandidate = pfAK05_btag[i] + pfAK05_btag[j];
	      if (bbCandidate.m()>110 && bbCandidate.m()<140) {
		numBox = 1;
		if ( fabs(bbCandidate.m()-125.) < minHiggsMassDiff || minHiggsMassDiff < 0){
		  minHiggsMassDiff = fabs(bbCandidate.m()-125.);	    
		  bestBjetIndex1 = i;
		  bestBjetIndex2 = j;
		  if (pfAK05_btag[bestBjetIndex1].pt() < pfAK05_btag[bestBjetIndex2].pt()) {
		    bestBjetIndex1 = j;
		    bestBjetIndex2 = i;
		  }
		}
	      }
	    }
	  }
	}
	if (numBox!=0 && numBox!=1){
	  for (int i = 0; i < pfAK05_btag.size(); i++){
	    for (int j = i+1; j < pfAK05_btag.size(); j++){	      
	      fastjet::PseudoJet bbCandidate = pfAK05_btag[i] + pfAK05_btag[j];
	      if ( bbCandidate.m()>76 && bbCandidate.m()<106 ){
		numBox = 2;	    
		if ( fabs(bbCandidate.m()-90.2) < minZMassDiff || minZMassDiff < 0){
		  minZMassDiff = fabs(bbCandidate.m()-91.2);	    
		  bestBjetIndex1 = i;
		  bestBjetIndex2 = j;
		  if (pfAK05_btag[bestBjetIndex1].pt() < pfAK05_btag[bestBjetIndex2].pt()) {
		    bestBjetIndex1 = j;
		    bestBjetIndex2 = i;
		  }	    
		}
	      }
	    }
	  }
	}
	if (numBox==1 || numBox==2) {
	  bjet1Pt = pfAK05_btag[bestBjetIndex1].pt();      
	  bjet1Eta = pfAK05_btag[bestBjetIndex1].eta();
	  bjet1Phi = pfAK05_btag[bestBjetIndex1].phi();
	  bjet2Pt = pfAK05_btag[bestBjetIndex2].pt();      
	  bjet2Eta = pfAK05_btag[bestBjetIndex2].eta();
	  bjet2Phi = pfAK05_btag[bestBjetIndex2].phi();
	  //fastjet::PseudoJet bjet1(bjet1Pt*cos(bjet1Phi), bjet1Pt*sin(bjet1Phi), bjet1Pt*sinh(bjet1Eta), bjet1Pt*cosh(bjet1Eta));
	  //fastjet::PseudoJet bjet2(bjet2Pt*cos(bjet2Phi), bjet2Pt*sin(bjet2Phi), bjet2Pt*sinh(bjet2Eta), bjet2Pt*cosh(bjet2Eta));
	}
	*/
      }
    else {
      vector<fastjet::PseudoJet> empty;
      vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty); //note that this I think includes photons
      fastjet::JetDefinition AK05_def(fastjet::antikt_algorithm, 0.5);
      fastjet::ClusterSequence pfAK05ClusterSequence = JetMaker(JetsConst, AK05_def);
      pfAK05 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK05ClusterSequence.inclusive_jets()),30., 3.0); //only cluster jets with > 30 GeV, eta < 3.0
    }


    //----------------------------
    //Find H->bb, Z->bb Candidates
    //----------------------------
    double mbbH = .0;
    double mbbZ = .0;
    for( int i = 0; i < numBJets; i++)
      {
	for(int j = i+1; j < numBJets; j++)
	  {
	    double mbb = (pfJetsBtag[i] + pfJetsBtag[j]).M();
	    //if mbb is closer to the higgs mass than mbbH, make mbbH = mbb
	    if( fabs(mbbH - 125.0) > fabs(mbb - 125.0) ) mbbH = mbb;
	    //same for mbbZ
	    if( fabs(mbbZ - 91.2) > fabs(mbb - 91.2) ) mbbZ = mbb;
	  }//end second jet loop
      }//end first jet loop
    
    //---------------------------------------
    //I n p u t   t o   H e m i s p h e r e s
    //---------------------------------------
    std::vector<TLorentzVector> JetsPlusHiggsCandidate;
    JetsPlusHiggsCandidate.push_back( HiggsCandidate );//Adding Higgs
    if( _delphesFormat )
      {
	//-----------------------------------
	//D e l p h e s  j e t s + h i g  g s
	//-----------------------------------
	int i = 0;
	for( auto tmp : pfJets )
	  {
	    jetPt[i]  = tmp.Pt();
	    jetEta[i] = tmp.Eta();
	    jetPhi[i] = tmp.Phi();
	    JetsPlusHiggsCandidate.push_back( tmp );
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
	for(int i = 0; i < pfAK05.size(); i++)
	  {
	    jetPt[i] = pfAK05[i].pt();
	    jetEta[i] = pfAK05[i].eta();
	    jetPhi[i] = pfAK05[i].phi();
	    JetsPlusHiggsCandidate.push_back( ConvertTo4Vector( pfAK05[i] ) );
	  }
	//---
	//MET
	//---
	GenMET();
	PFMET = genMET;
	MET = PFMET.pt();
      }
    
    //--------------------------
    //higgs+jets size at least 2
    //--------------------------
    if( JetsPlusHiggsCandidate.size() < 2 )
      {
	if( _debug ) std::cout << "[DEBUG]: Not enough objects for razor computation" << std::endl;
	continue;
      }
    if ( _debug ) std::cout << "--> before CMSHEM" << std::endl;
    CMSHemisphere* myHem = new CMSHemisphere( JetsPlusHiggsCandidate );
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
    dPhiHiggsJ1 = DeltaPhi( HiggsCandidate, j1);
    dPhiHiggsJ2 = DeltaPhi( HiggsCandidate, j2);
    MR = CalcMR(j1, j2);
    RSQ = pow(CalcMRT(j1, j2, PFMET),2.)/MR/MR;
    if ( _debug ) std::cout << "[DEBUG]: numBox-> "<< numBox << std::endl;

    //---------------------------------------
    //M a k i n g   B o x   H i e r a r c h y
    //---------------------------------------
    if ( higgsPt > 110. )
      {
	numBox = 0;
      }
    else if ( mbbH > 110. && mbbH < 140. )
      {
	numBox = 1;
      }
    else if ( mbbZ > 76.  && mbbZ < 106. )
      {
	numBox = 2;
      }
    else
      {
	numBox = 3;
      }
    
    // write event in the tree
    outTree->Fill();
    
    // fill PDF histograms
    if(numBox == 0) {
      if (MR>=150 && MR<200 && RSQ>=0.00 && RSQ<0.05) { pdfHighPt->Fill(0); pdfTotal->Fill(0); }
      else if (MR>=150 && MR<200 && RSQ>=0.05 && RSQ<0.10) { pdfHighPt->Fill(1); pdfTotal->Fill(1); }
      else if (MR>=150 && MR<200 && RSQ>=0.10 && RSQ<0.15) { pdfHighPt->Fill(2); pdfTotal->Fill(2); }
      else if (MR>=150 && MR<200 && RSQ>=0.15 && RSQ<0.20) { pdfHighPt->Fill(3); pdfTotal->Fill(3); }
      else if (MR>=150 && MR<200 && RSQ>=0.20 && RSQ<1.00) { pdfHighPt->Fill(4); pdfTotal->Fill(4); }
      else if (MR>=200 && MR<300 && RSQ>=0.00 && RSQ<0.05) { pdfHighPt->Fill(5); pdfTotal->Fill(5); }
      else if (MR>=200 && MR<300 && RSQ>=0.05 && RSQ<0.10) { pdfHighPt->Fill(6); pdfTotal->Fill(6); }
      else if (MR>=200 && MR<300 && RSQ>=0.10 && RSQ<0.15) { pdfHighPt->Fill(7); pdfTotal->Fill(7); }
      else if (MR>=200 && MR<300 && RSQ>=0.15 && RSQ<1.00) { pdfHighPt->Fill(8); pdfTotal->Fill(8); }
      else if (MR>=300 && MR<500 && RSQ>=0.00 && RSQ<0.05) { pdfHighPt->Fill(9); pdfTotal->Fill(9); }
      else if (MR>=300 && MR<500 && RSQ>=0.05 && RSQ<0.10) { pdfHighPt->Fill(10); pdfTotal->Fill(10); }
      else if (MR>=300 && MR<500 && RSQ>=0.10 && RSQ<1.00) { pdfHighPt->Fill(11); pdfTotal->Fill(11); }
      else if (MR>=500 && MR<1600 && RSQ>=0.00 && RSQ<0.05) { pdfHighPt->Fill(12); pdfTotal->Fill(12); }
      else if (MR>=500 && MR<1600 && RSQ>=0.05 && RSQ<1.00) { pdfHighPt->Fill(13); pdfTotal->Fill(13); }
      else if (MR>=1600 && MR<3000 && RSQ>=0.00 && RSQ<1.00) { pdfHighPt->Fill(14); pdfTotal->Fill(14); }
    }
    if(numBox == 1) {
      if (MR>=150 && MR<300 && RSQ>=0.00 && RSQ<0.05) { pdfHbb->Fill(0); pdfTotal->Fill(15); }
      else if (MR>=150 && MR<300 && RSQ>=0.05 && RSQ<1.00) { pdfHbb->Fill(1); pdfTotal->Fill(16); }
      else if (MR>=300 && MR<3000 && RSQ>=0.00 && RSQ<1.00) { pdfHbb->Fill(2); pdfTotal->Fill(17); }
    }
    if(numBox == 2) {
      if (MR>=150 && MR<450 && RSQ>=0.00 && RSQ<0.05) { pdfZbb->Fill(0); pdfTotal->Fill(18); }
      else if (MR>=150 && MR<450 && RSQ>=0.05 && RSQ<1.00) { pdfZbb->Fill(1); pdfTotal->Fill(19); }
      else if (MR>=450 && MR<3000 && RSQ>=0.00 && RSQ<1.00) { pdfZbb->Fill(2); pdfTotal->Fill(20); }
    }
    if(numBox == 3) {
      pdfHighResRsqMR->Fill(MR,RSQ);
      if (MR>=150 && MR<200 && RSQ>=0.00 && RSQ<0.05) { pdfHighRes->Fill(0); pdfTotal->Fill(21); }
      else if (MR>=150 && MR<250 && RSQ>=0.05 && RSQ<0.10) { pdfHighRes->Fill(1); pdfTotal->Fill(22); }
      else if (MR>=150 && MR<250 && RSQ>=0.10 && RSQ<0.15) { pdfHighRes->Fill(2); pdfTotal->Fill(23); }
      else if (MR>=150 && MR<250 && RSQ>=0.15 && RSQ<1.00) { pdfHighRes->Fill(3); pdfTotal->Fill(24); }
      else if (MR>=250 && MR<400 && RSQ>=0.00 && RSQ<0.05) { pdfHighRes->Fill(4); pdfTotal->Fill(25); }
      else if (MR>=250 && MR<400 && RSQ>=0.05 && RSQ<0.10) { pdfHighRes->Fill(5); pdfTotal->Fill(26); }
      else if (MR>=250 && MR<400 && RSQ>=0.10 && RSQ<1.00) { pdfHighRes->Fill(6); pdfTotal->Fill(27); }
      else if (MR>=400 && MR<1400 && RSQ>=0 && RSQ<0.05) { pdfHighRes->Fill(7); pdfTotal->Fill(28); }
      else if (MR>=400 && MR<1400 && RSQ>=0.05 && RSQ<1.00) { pdfHighRes->Fill(8); pdfTotal->Fill(29); }
      else if (MR>=1400 && MR<3000 && RSQ>=0 && RSQ<1.00) { pdfHighRes->Fill(9); pdfTotal->Fill(30); }
    }
    
  }

  // full event TTree
  outTree->Write();

  
  // eff TTree
  double effHighPt = _filterEff*pdfHighPt->Integral()/double(nentries);
  double effHbb = _filterEff*pdfHbb->Integral()/double(nentries);
  double effZbb = _filterEff*pdfZbb->Integral()/double(nentries);
  double effHighRes = _filterEff*pdfHighRes->Integral()/double(nentries);
  double effTotal = _filterEff*pdfTotal->Integral()/double(nentries);
  
  // normalize the PDFs
  if(pdfHighPt->Integral()>0)  pdfHighPt->Scale(1./pdfHighPt->Integral());
  if(pdfHbb->Integral()>0)  pdfHbb->Scale(1./pdfHbb->Integral());
  if(pdfZbb->Integral()>0)  pdfZbb->Scale(1./pdfZbb->Integral());
  if(pdfHighRes->Integral()>0)  pdfHighRes->Scale(1./pdfHighRes->Integral());
  if(pdfTotal->Integral()>0)  pdfTotal->Scale(1./pdfTotal->Integral());
  
  // write the PDFs
  pdfHighPt->Write();  
  pdfHbb->Write();
  pdfZbb->Write();
  pdfHighRes->Write();
  pdfTotal->Write();
  
  pdfHighResRsqMR->Write();

  char outname[256];
  sprintf(outname,"data/%s.root", _analysis.c_str());
  TH1D* xsecProbHighPt = XsecProb(pdfHighPt, effHighPt, outname, "HighPt", 100, 0., 5.0);
  TH1D* xsecProbHbb = XsecProb(pdfHbb, effHbb, outname, "Hbb", 100, 0., 5.0);
  TH1D* xsecProbZbb = XsecProb(pdfZbb, effZbb, outname, "Zbb", 100, 0., 5.0);
  TH1D* xsecProbHighRes = XsecProb(pdfHighRes, effHighRes, outname, "HighRes", 100, 0., 5.0);
  TH1D* xsecProbTotal = XsecProb(pdfTotal, effTotal, outname, "Total", 100, 0., 5.0);
  // Open Output file again 
  file->cd();
  
  double xsecULHighPt = _statTools->FindUL(xsecProbHighPt, 0.95, 1.);
  double xsecULHbb = _statTools->FindUL(xsecProbHbb, 0.95, 1.);
  double xsecULZbb = _statTools->FindUL(xsecProbZbb, 0.95, 1.);
  double xsecULHighRes = _statTools->FindUL(xsecProbHighRes, 0.95, 1.);
  double xsecULTotal = _statTools->FindUL(xsecProbTotal, 0.95, 1.);
  
  TTree* effTree = new TTree("RazorInclusiveEfficiency","RazorInclusiveEfficiency");
  effTree->Branch("effHighPt", &effHighPt, "effHighPt/D");
  effTree->Branch("effHbb", &effHbb, "effHbb/D");
  effTree->Branch("effZbb", &effZbb, "effZbb/D");
  effTree->Branch("effHighRes", &effHighRes, "effHighRes/D");
  effTree->Branch("effTotal", &effTotal, "effTotal/D");
  effTree->Branch("xsecULHighPt", &xsecULHighPt, "xsecULHighPt/D");
  effTree->Branch("xsecULHbb", &xsecULHbb, "xsecULHbb/D");
  effTree->Branch("xsecULZbb", &xsecULZbb, "xsecULZbb/D");
  effTree->Branch("xsecULHighRes", &xsecULHighRes, "xsecULHighRes/D");
  effTree->Branch("xsecULTotal", &xsecULTotal, "xsecULTotal/D");
  effTree->Fill();
  effTree->Write();

  xsecProbHighPt->Write();
  xsecProbHbb->Write();
  xsecProbZbb->Write();
  xsecProbHighRes->Write();
  xsecProbTotal->Write();
  //std::cout << "before malloc" << std::endl;
  file->Close();
  //std::cout << "after malloc" << std::endl;
}

double CMSRazorHgg::DeltaPhi(TLorentzVector jet1, TLorentzVector jet2) {
  double deltaPhi = jet1.Phi() - jet2.Phi();
  while (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  while (deltaPhi <= -M_PI) deltaPhi += 2*M_PI;
  return deltaPhi;
}


TH1D* CMSRazorHgg::XsecProb(TH1D* sigPdf, double eff, string Filename, string directory, int ibin, double xmin, double xmax) {
  
  int ibinX = sigPdf->GetXaxis()->GetNbins();
  int ibinY = sigPdf->GetYaxis()->GetNbins();
  
  char histname[256];
  sprintf(histname, "%s%s","probVec", directory.c_str());
  TH1D* probVec = new TH1D(histname, histname, ibin, xmin, xmax);
  
  TFile* likFile = new TFile(TString(Filename));  
  gROOT->cd();
  // a loop over xsec should go here... 
  for(int i=0; i<ibin; i++) {
    double xsec = xmin + (i+0.5)/ibin*(xmax-xmin);
    double prob = 1;
    double logProb = 0;
    for(int ix=0; ix<ibinX; ix++) {
      for(int iy=0; iy<ibinY; iy++) {
	//if (!((ix==8) && iy==0)) continue; //just using the bin with the excess
	if (sigPdf->GetBinContent(ix+1,iy+1)<=0.) continue;
	double sBin = _Lumi*xsec*eff*sigPdf->GetBinContent(ix+1,iy+1);	
	//cout << "xsec = " << xsec << endl;
	//cout << "sBin = " << sBin << endl;
	if (sBin >= 100) cout << "Note: signal events exceed 100! sBin = " << sBin << endl;
	char name[256];
	sprintf(name, "%s/lik_%i_%i", directory.c_str(), ix, iy);
	TH1D* binProb = (TH1D*) likFile->Get(name);
	//if(prob < 10.e-30) prob = 0.;
	//if(prob <= 0) continue;
	logProb += TMath::Log(binProb->GetBinContent(binProb->FindBin(sBin)));
	prob *= binProb->GetBinContent(binProb->FindBin(sBin)); 
	//cout << "logProb = " << logProb << endl;
	//cout << "prob = " << prob << endl;
	delete binProb;
      }
    }
    //cout << "prob = " << prob << endl;
    //cout << "exp(logprob) = " << TMath::Exp(logprob) << endl;
    //probVec->SetBinContent(i+1,prob);
    probVec->SetBinContent(i+1,TMath::Exp(logProb));
  }
  probVec->Scale(1./probVec->Integral());
  likFile->Close();
  return probVec;
}
