//-------------------------------------------------------
// Description:
// Runs Hgg analysis
// Authors: 
//-------------------------------------------------------


#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <TTree.h>
#include "CMS/CMSRazorHgg.hh"
#include <fastjet/tools/Pruner.hh>

CMSRazorHgg::CMSRazorHgg(TTree *tree, double Lumi, string analysis, bool delphesFormat) : CMSReco(tree, delphesFormat) {
  _Lumi = Lumi;
  _statTools = new StatTools(-99);
  _analysis = analysis;
  _delphesFormat = delphesFormat;
}

CMSRazorHgg::~CMSRazorHgg(){
}

void CMSRazorHgg::SetSqrts(double sqrts) {
  _sqrts = sqrts;
}

// loop over events - real analysis
void CMSRazorHgg::Loop(string outFileName) {
  srand(time(0)); 
  if( _delphesFormat && DelphesTree::fChain == 0) return;
  if (!_delphesFormat && DetectorBase::fChain == 0) return;
  
  cout << "starting Loop" << endl;

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
  
  cout << "open output file" << endl;
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

  
  //double xedge[9] = {0, 100, 200, 300, 400, 500, 700, 1000, 2500};
  //double yedge[7] = {0, 0.10, 0.20, 0.30, 0.40, 0.50, 1.0};
  TH1D* pdfHighPt = new TH1D("pdfHighPt","pdfHighPt",15,0,15);
  TH1D* pdfHbb = new TH1D("pdfHbb","pdfHbb",3,0,3);
  TH1D* pdfZbb = new TH1D("pdfZbb","pdfZbb",3,0,3);
  TH1D* pdfHighRes = new TH1D("pdfHighRes","pdfHighRes",10,0,10);


  cout << "getting number of entries" << endl;
  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = 0;
  if (_delphesFormat) nentries = DelphesTree::fChain->GetEntries();
  else nentries = DetectorBase::fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  // set the by-event weight
  for (Long64_t jentry=0; jentry<nentries;jentry+=1) {

    if(verbose) cout << "new event" << endl;

    // clean physics-objects blocks
    CleanEvent();

    // get new event
    if (_delphesFormat) {
      Long64_t ientry = DelphesTree::LoadTree(jentry);      
      if (ientry < 0) break;
      nb = DelphesTree::fChain->GetEntry(jentry); nbytes += nb;
    }
    else {
      Long64_t ientry = DetectorBase::LoadTree(jentry);
      nb = DetectorBase::fChain->GetEntry(jentry); nbytes += nb;
    }
    
    if (jentry%1 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

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
    double bestSumPt = -1;
    int bestPhotIndex1 = -1;
    int bestPhotIndex2 = -1;
    double bestMass = -1;
    int numHiggs = 0;
    double secondBestSumPt = -1;
    int secondBestPhotIndex1 = -1;
    int secondBestPhotIndex2 = -1;
    double secondBestMass = -1; 
    numBox = -1;

    //loop over all pairs of gen photons and find the pair with largest scalar sum Pt
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

    fastjet::PseudoJet bestDiphoton;
    bestDiphoton = _PFPhotons[bestPhotIndex1] + _PFPhotons[bestPhotIndex2];

    std::cout << "-->pass1" << std::endl;
    //fill photon and bjet variables
    pho1Pt = _PFPhotons[bestPhotIndex1].pt();
    pho1Eta = _PFPhotons[bestPhotIndex1].eta();
    pho1Phi = _PFPhotons[bestPhotIndex1].phi();
    pho2Pt = _PFPhotons[bestPhotIndex2].pt();
    pho2Eta = _PFPhotons[bestPhotIndex2].eta();
    pho2Phi = _PFPhotons[bestPhotIndex2].phi();
    higgsPt = bestDiphoton.pt();
    
    if (higgsPt > 110){
      //cout << "higgsPt: "<<higgsPt<<endl;
      numBox = 0;
    }
    
    higgsEta = bestDiphoton.eta();
    higgsPhi = bestDiphoton.phi();
    higgsMass = bestDiphoton.m();
    
    fastjet::PseudoJet pho1(pho1Pt*cos(pho1Phi), pho1Pt*sin(pho1Phi), pho1Pt*sinh(pho1Eta), pho1Pt*cosh(pho1Eta));
    fastjet::PseudoJet pho2(pho2Pt*cos(pho2Phi), pho2Pt*sin(pho2Phi), pho2Pt*sinh(pho2Eta), pho2Pt*cosh(pho2Eta));

    fastjet::PseudoJet higgs = pho1 + pho2;
    double higgsPhi = higgs.phi();
    double higgsEta = higgs.eta();
    double higgsPt = higgs.pt();
    double higgsEnergy = higgs.E();
    TLorentzVector higgsvector;
    higgsvector.SetPtEtaPhiE(higgsPt, higgsEta, higgsPhi, higgsEnergy);

    // AK5 jets    
    vector<fastjet::PseudoJet> pfAK05;
    vector<fastjet::PseudoJet> pfAK05_btag;

    if (_delphesFormat) {
      for (int iJet = 0; iJet<Jet_size; iJet++){	
	TLorentzVector v;
	v.SetPtEtaPhiM(Jet_PT[iJet],Jet_Eta[iJet],Jet_Phi[iJet],Jet_Mass[iJet]);
	if (v.Pt()>30. && v.Eta()<3.) {	  
	  //check if within DR < 0.5 of a selected photon
	  fastjet::PseudoJet pv = ConvertToPseudoJet(v);
	  double thisDR = min(pv.delta_R(pho1), pv.delta_R(pho2));	  
	  if (thisDR > 0.5) {
	    numJets++;
	    pfAK05.push_back(fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()));
	    if (Jet_BTag[iJet]) {
	      numBJets++; 
	      pfAK05_btag.push_back(fastjet::PseudoJet(v.Px(), v.Py(), v.Pz(), v.E()));
	    }
	  }
	}
      }
      
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
    }
    else {
      vector<fastjet::PseudoJet> empty;
      vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty); //note that this I think includes photons
      fastjet::JetDefinition AK05_def(fastjet::antikt_algorithm, 0.5);
      fastjet::ClusterSequence pfAK05ClusterSequence = JetMaker(JetsConst, AK05_def);
      pfAK05 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK05ClusterSequence.inclusive_jets()),30., 3.0); //only cluster jets with > 30 GeV, eta < 3.0
    }
    
    if(pfAK05.size()<1){
      cout << "No jets..." << endl;
      continue;
    }

    vector<fastjet::PseudoJet> jetsForHemispheres;      
    jetsForHemispheres.push_back(higgs);

    for(int i = 0; i < pfAK05.size(); i++){
      jetPt[i] = pfAK05[i].pt();
      jetEta[i] = pfAK05[i].eta();
      jetPhi[i] = pfAK05[i].phi();
      jetsForHemispheres.push_back(pfAK05[i]);
    }
       
    if (numBox<0){
      numBox = 3;
    }

    if (_delphesFormat){
      TVector3 v3;
      v3.SetPtEtaPhi(MissingET_MET[0],MissingET_Eta[0],MissingET_Phi[0]);
      PFMET = fastjet::PseudoJet(v3.Px(), v3.Py(), 0., v3.Mag());
      MET = PFMET.pt();
    }
    else{
      GenMET();
      PFMET = genMET;
      MET = PFMET.pt();
    }

    //traditional hemispheres
    if(jetsForHemispheres.size() < 2){
        cout << "Not enough objects for razor computation" << endl;
        continue;
    }
    std::cout << "--> before CMSHEM" << std::endl;
    CMSHemisphere* myHem = new CMSHemisphere(ConvertTo4Vector(jetsForHemispheres));
    myHem->CombineMinMass();
    vector<TLorentzVector> hem = myHem->GetHemispheres();
    vector<int> Temporary = myHem->GetHem1Constituents();
    vector<int> Temporary2 = myHem->GetHem2Constituents();
    std::cout << "--> pass CMSHEM" << std::endl;
    
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

    std::cout << "delete Hem" << std::endl;
    //compute traditional RSQ and MR
    TLorentzVector j1 = hem[0];
    TLorentzVector j2 = hem[1];  
    dPhi = DeltaPhi(j1,j2); //deltaPhi
    dPhiHiggsJ1 = DeltaPhi(higgsvector,j1);
    dPhiHiggsJ2 = DeltaPhi(higgsvector,j2);
    MR = CalcMR(j1, j2);
    RSQ = pow(CalcMRT(j1, j2, PFMET),2.)/MR/MR;
    cout << "numBox: "<<numBox<<endl;

    // write event in the tree
    outTree->Fill();
    
    // fill PDF histograms
    if(numBox == 0) pdfHighPt->Fill(0);
    if(numBox == 1) pdfHbb->Fill(0);
    if(numBox == 2) pdfZbb->Fill(0);
    if(numBox == 3) {
      if (MR>=150 && MR<200 && RSQ>=0.00 && RSQ<0.05) pdfHighRes->Fill(0);
      else if (MR>=150 && MR<250 && RSQ>=0.05 && RSQ<0.10) pdfHighRes->Fill(1);
      else if (MR>=150 && MR<250 && RSQ>=0.10 && RSQ<0.15) pdfHighRes->Fill(2);
      else if (MR>=150 && MR<250 && RSQ>=0.15 && RSQ<1.00) pdfHighRes->Fill(3);
      else if (MR>=250 && MR<400 && RSQ>=0.00 && RSQ<0.05) pdfHighRes->Fill(4);
      else if (MR>=250 && MR<400 && RSQ>=0.05 && RSQ<0.10) pdfHighRes->Fill(5);
      else if (MR>=250 && MR<400 && RSQ>=0.10 && RSQ<1.00) pdfHighRes->Fill(6);
      else if (MR>=400 && MR<1400 && RSQ>=0 && RSQ<0.05) pdfHighRes->Fill(7);
      else if (MR>=400 && MR<1400 && RSQ>=0.05 && RSQ<1.00) pdfHighRes->Fill(8);
      else if (MR>=1400 && MR<3000 && RSQ>=0 && RSQ<1.00) pdfHighRes->Fill(9);
    }
    
  }

  // full event TTree
  outTree->Write();

  
  // eff TTree
  double effHighPt = pdfHighPt->Integral()/double(nentries);
  double effHbb = pdfHbb->Integral()/double(nentries);
  double effZbb = pdfZbb->Integral()/double(nentries);
  double effHighRes = pdfHighRes->Integral()/double(nentries);
  
  // normalize the PDFs
  if(pdfHighPt->Integral()>0)  pdfHighPt->Scale(1./pdfHighPt->Integral());
  if(pdfHbb->Integral()>0)  pdfHbb->Scale(1./pdfHbb->Integral());
  if(pdfZbb->Integral()>0)  pdfZbb->Scale(1./pdfZbb->Integral());
  if(pdfHighRes->Integral()>0)  pdfHighRes->Scale(1./pdfHighRes->Integral());
  
  // write the PDFs
  pdfHighPt->Write();  
  pdfHbb->Write();
  pdfZbb->Write();
  pdfHighRes->Write();

  //char outname[256];
  //sprintf(outname,"data/%s.root", _analysis.c_str());
  //TH1D* xsecProb = XsecProb(pdfHighRes, effHighRes,name, 1000, 0., 1.);
  // Open Output file again 
  file->cd();
  double xsecULHighRes = 0.0;// _statTools->FindUL(xsecProb, 0.95, 1.);
  
  TTree* effTree = new TTree("RazorInclusiveEfficiency","RazorInclusiveEfficiency");
  effTree->Branch("effHighPt", &effHighPt, "effHighPt/D");
  effTree->Branch("effHbb", &effHbb, "effHbb/D");
  effTree->Branch("effZbb", &effZbb, "effZbb/D");
  effTree->Branch("effHighRes", &effHighRes, "effHighRes/D");
  effTree->Branch("xsecULHighRes", &xsecULHighRes, "xsecULHighRes/D");
  effTree->Fill();
  effTree->Write();
  
  //  xsecProb->Write();
  file->Close();
  


}

double CMSRazorHgg::DeltaPhi(TLorentzVector jet1, TLorentzVector jet2) {
  double deltaPhi = jet1.Phi() - jet2.Phi();
  while (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  while (deltaPhi <= -M_PI) deltaPhi += 2*M_PI;
  return deltaPhi;
}


TH1D* CMSRazorHgg::XsecProb(TH1D* sigPdf, double eff, TString Filename, int ibin, double xmin, double xmax) {
  
  int ibinX = sigPdf->GetXaxis()->GetNbins();
  int ibinY = sigPdf->GetYaxis()->GetNbins();
  
  TH1D* probVec = new TH1D("probVec", "probVec", ibin, xmin, xmax);
  
  TFile* likFile = new TFile(Filename);  
  gROOT->cd();
  // a loop over xsec should go here... 
  for(int i=0; i<ibin; i++) {
    double xsec = xmin + (i+0.5)/ibin*(xmax-xmin);
    double prob = 1;
    for(int ix=0; ix<ibinX; ix++) {
      for(int iy=0; iy<ibinY; iy++) {
	double sBin = _Lumi*xsec*eff*sigPdf->GetBinContent(ix,iy);
	if(sBin <= 0.) continue;
	char name[256];
	sprintf(name, "lik_%i_%i", ix, iy);
	TH1D* binProb = (TH1D*) likFile->Get(name);
	if(prob < 10.e-30) prob = 0.;
	if(prob <= 0) continue;
	prob *= binProb->GetBinContent(binProb->FindBin(sBin));
	delete binProb;
      }
    }
    probVec->SetBinContent(i+1,prob);
  }
  likFile->Close();
  return probVec;
}
