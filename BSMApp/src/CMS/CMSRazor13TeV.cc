#include <string>
#include <iostream>
#include <TTree.h>
#include "CMS/CMSRazor13TeV.hh"
#include <fastjet/tools/Pruner.hh>
#include <vector>

CMSRazor13TeV::CMSRazor13TeV(TTree *tree, double Lumi, string analysis) : CMSReco(tree) {
  _Lumi = Lumi;
  _statTools = new StatTools(-99);
  _analysis = analysis;
  _event_counter = 0;
}

CMSRazor13TeV::~CMSRazor13TeV(){
}

void CMSRazor13TeV::SetSqrts(double sqrts) {
	_sqrts = sqrts;
}


// loop over events - real analysis
void CMSRazor13TeV::Loop(string outFileName) {

	if(DetectorBase::fChain == 0) return;

	double MR_Default, RSQ_Default, MR_KT_Jets, RSQ_KT_Jets, MR_AKT_Jets, RSQ_AKT_Jets, MR_CAM_Jets, RSQ_CAM_Jets;
	double MRz_Default, MRz_KT_Jets, MRz_AKT_Jets, MRz_CAM_Jets, MRz_CAM_Cons, RSQz_Default, RSQz_KT_Jets, RSQz_AKT_Jets, RSQz_CAM_Jets;

	double sqrtsR_Default, gammaRp1_Default, pTCM_Default, gammaRp1Ana_Default;
	double W_EFF;
	int BOX_NUM;
	int SizesList[5];
	
	// Open Output file
	TFile *file = new TFile(outFileName.c_str(),"UPDATE"); 	
	
	cout << "writing to ";
	cout << outFileName << endl;
		
	TTree* outTree = new TTree("RazorInclusive","RazorInclusive");
	
	outTree->Branch("MR_Default", &MR_Default, "MR_Default/D");
	outTree->Branch("MRz_Default", &MRz_Default, "MRz_Default/D");
	outTree->Branch("RSQ_Default", &RSQ_Default, "RSQ_Default/D");
	outTree->Branch("RSQz_Default", &RSQz_Default, "RSQz_Default/D");

	outTree->Branch("sqrtsR_Default", &sqrtsR_Default, "sqrtsR_Default/D");
	outTree->Branch("gammaRp1_Default", &gammaRp1_Default, "gammaRp1_Default/D");
	outTree->Branch("gammaRp1Ana_Default", &gammaRp1Ana_Default, "gammaRp1Ana_Default/D");
	outTree->Branch("pTCM_Default", &pTCM_Default, "pTCM_Default/D");
		
	outTree->Branch("MR_KT_Jets", &MR_KT_Jets, "MR_KT_Jets/D");
	outTree->Branch("MRz_KT_Jets", &MRz_KT_Jets, "MRz_KT_Jets/D");
	outTree->Branch("RSQ_KT_Jets", &RSQ_KT_Jets, "RSQ_KT_Jets/D");
	outTree->Branch("RSQz_KT_Jets", &RSQz_KT_Jets, "RSQz_KT_Jets/D");
	
	outTree->Branch("MR_AKT_Jets", &MR_AKT_Jets, "MR_AKT_Jets/D");
	outTree->Branch("MRz_AKT_Jets", &MRz_AKT_Jets, "MRz_AKT_Jets/D");
	outTree->Branch("RSQ_AKT_Jets", &RSQ_AKT_Jets, "RSQ_AKT_Jets/D");
	outTree->Branch("RSQz_AKT_Jets", &RSQz_AKT_Jets, "RSQz_AKT_Jets/D");
		
	outTree->Branch("MR_CAM_Jets", &MR_CAM_Jets, "MR_CAM_Jets/D");
	outTree->Branch("MRz_CAM_Jets", &MRz_CAM_Jets, "MRz_CAM_Jets/D");
	outTree->Branch("RSQ_CAM_Jets", &RSQ_CAM_Jets, "RSQ_CAM_Jets/D");
	outTree->Branch("RSQz_CAM_Jets", &RSQz_CAM_Jets, "RSQz_CAM_Jets/D");

	outTree->Branch("SizesList", SizesList, "SizesList[5]/I");
	
	
	outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
	outTree->Branch("W_EFF", &W_EFF, "W_EFF/D");
	
	// loop over entries
	Long64_t nbytes = 0, nb = 0;
	Long64_t nentries = DetectorBase::fChain->GetEntries();
	std::cout << "Number of entries = " << nentries << std::endl;
	
	// set the by-event weight
	W_EFF = 1./nentries;
	//set box
	BOX_NUM = 0;
	
       

	int event_counter = _event_counter;
  for (Long64_t jentry=0; jentry<nentries;jentry+=1) {

    if(verbose) cout << "new event" << endl;

    // clean physics-objects blocks
    CleanEvent();

    // get new event
    Long64_t ientry = DetectorBase::LoadTree(jentry);
    if (ientry < 0) break;
    nb = DetectorBase::fChain->GetEntry(jentry);   nbytes += nb;
	//if (jentry > 10000) continue; //use for gen level comparison coding
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
	
    // Build the event at generator level
    PFReco();
    vector<fastjet::PseudoJet> empty;
    vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty);
	  
    // wide jets
    fastjet::JetDefinition CA08_def(fastjet::cambridge_algorithm, 0.8);
    fastjet::ClusterSequence pfCA08ClusterSequence = JetMaker(JetsConst, CA08_def);
    vector<fastjet::PseudoJet> pfCA08 = SelectByAcceptance(fastjet::sorted_by_pt(pfCA08ClusterSequence.inclusive_jets()),40., 2.4);
    fastjet::Pruner pruner(CA08_def, 0.1, 0.25);

    // narrow jets
    fastjet::JetDefinition AK04_def(fastjet::antikt_algorithm, 0.4);
    fastjet::ClusterSequence pfAK04ClusterSequence = JetMaker(JetsConst, AK04_def);
    vector<fastjet::PseudoJet> pfAK04 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK04ClusterSequence.inclusive_jets()),40., 2.4);
	//vector<fastjet::PseudoJet> pfAK04 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK04ClusterSequence.inclusive_jets()),18., 2.4);

    if(pfAK04.size()<2) continue;
    event_counter = event_counter + 1 ;
	  
    GenMET();
    PFMET = genMET;
	  
	  
    
    // 1a) cluster hemispheres using kT in exclusive mode, using jets as ingredients
    // some backward compatibility test will be needed here
    
	double cone_size = 7.0;  
	  
	fastjet::ClusterSequence cs1(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, cone_size));
    vector<TLorentzVector> hem_KT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs1.exclusive_jets(2))); //this is from original, written correctly
	//vector<TLorentzVector> hem_KT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs1.inclusive_jets()));
	  
    fastjet::ClusterSequence cs2(pfAK04, fastjet::JetDefinition(fastjet::antikt_algorithm, cone_size));
	vector<TLorentzVector> hem_AKT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(improper_Exclusive_Jets(cs2)));
	//vector<TLorentzVector> hem_AKT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs2.inclusive_jets())); 
	  
    fastjet::ClusterSequence cs3(pfAK04, fastjet::JetDefinition(fastjet::cambridge_algorithm, cone_size));
	vector<TLorentzVector> hem_CAM_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(improper_Exclusive_Jets(cs3)));
	//vector<TLorentzVector> hem_CAM_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs3.inclusive_jets()));
	  
	  
	// 1b) traditional hemispheres
    CMSHemisphere* myHem = new CMSHemisphere(ConvertTo4Vector(pfAK04));
    myHem->CombineMinMass();
    vector<TLorentzVector> hem_Default = myHem->GetHemispheres();
    delete myHem;
	  
	
	//initialize  
    MR_KT_Jets = -9999. ;
    RSQ_KT_Jets = -9999. ;
    MR_AKT_Jets = -9999. ;
    RSQ_AKT_Jets = -9999. ;
    MR_CAM_Jets = -9999. ;
    RSQ_CAM_Jets = -9999. ;
    MR_Default = -9999. ;
    RSQ_Default = -9999. ;
    MRz_KT_Jets = -9999. ;
    MRz_AKT_Jets = -9999. ;
    MRz_CAM_Jets = -9999. ;
    MRz_Default = -9999. ;
    RSQz_KT_Jets = -9999. ;
    RSQz_AKT_Jets = -9999. ;
    RSQz_CAM_Jets = -9999. ;
    RSQz_Default = -9999. ;
	  

    sqrtsR_Default = -9999. ;
    pTCM_Default = -9999. ;
    gammaRp1_Default = -9999. ;
    gammaRp1Ana_Default = -9999. ;
	  
    SizesList[0] = -9999.0;
    SizesList[1] = -9999.0;
    SizesList[2] = -9999.0;
    SizesList[3] = -9999.0;
    SizesList[4] = -9999.0;
	  
    TLorentzVector j1;
    TLorentzVector j2;
 

    

    // 2a) compute new RSQ and MR vals
    if (hem_KT_Jets.size() > 1){
      j1 = hem_KT_Jets[0];
      j2 = hem_KT_Jets[1];  
      MR_KT_Jets = CalcMR(j1, j2);
      RSQ_KT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MR_KT_Jets/MR_KT_Jets;
	  MRz_KT_Jets = CalcMR_zinvariant(j1, j2);
	  RSQz_KT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_KT_Jets/MRz_KT_Jets;
	  
    }
	  
    if (hem_AKT_Jets.size() > 1) {
      j1 = hem_AKT_Jets[0];
      j2 = hem_AKT_Jets[1];  
      MR_AKT_Jets = CalcMR(j1, j2);
      RSQ_AKT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MR_AKT_Jets/MR_AKT_Jets;
	  MRz_AKT_Jets = CalcMR_zinvariant(j1, j2);
	  RSQz_AKT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_AKT_Jets/MRz_AKT_Jets;
    }
    
    if (hem_CAM_Jets.size() > 1) {
      j1 = hem_CAM_Jets[0];
      j2 = hem_CAM_Jets[1];  
      MR_CAM_Jets = CalcMR(j1, j2);
      RSQ_CAM_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MR_CAM_Jets/MR_CAM_Jets;
	  MRz_CAM_Jets = CalcMR_zinvariant(j1, j2);
      RSQz_CAM_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_CAM_Jets/MRz_CAM_Jets;
    }

    
    // 2b) compute traditional RSQ and MR (DEFAULT)
    j1 = hem_Default[0];
    j2 = hem_Default[1];  
    MR_Default = CalcMR(j1, j2);
    RSQ_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MR_Default/MR_Default;
    MRz_Default = CalcMR_zinvariant(j1, j2);
    RSQz_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_Default/MRz_Default;

    sqrtsR_Default = CalcSqrtsR(j1, j2, PFMET);
    gammaRp1_Default = CalcGammaRp1(j1, j2, PFMET);
    gammaRp1Ana_Default = CalcGammaRp1Ana(j1, j2, PFMET);
	  
    pTCM_Default = (j1+j2+PFMET).pt();
	  
    //delete j1;
    //delete j2;
	  
	
    SizesList[0] = hem_Default.size();
    SizesList[1] = hem_KT_Jets.size();
    SizesList[2] = hem_AKT_Jets.size();
    SizesList[3] = hem_CAM_Jets.size();
    SizesList[4] = pfAK04.size(); //pre cluster pfak04 jet number
    //cout <<endl << "Number of PFJets:" << pfAK04.size()<<endl;
	 
	 

    // write event in the tree
    outTree->Fill();
	
	  
	  
    continue;  	  
	
	  
    //BEGIN GEN LEVEL ANALYSIS
    //initialize
	
	  
    int correct_clustering = -9999.0 ;
	  
    //initialize structures to store gen data pulled from cmsreco
    int genParticle;
    vector<double> genParticlePx;
    vector<double> genParticlePy;
    vector<double> genParticlePz;
    vector<double> genParticleE;
    vector<int> genParticlePdgId;
    vector<int> genParticleM1PdgId;
    vector<int> genParticleStatus;
    vector<double> genParticlem1px;
    // edited this because function changed
    //GenReturn(genParticle, genParticlePx, genParticlePy, genParticlePz, genParticleE, genParticlePdgId, genParticleM1PdgId);
    GenReturn(genParticle, genParticlePx, genParticlePy, genParticlePz, genParticleE, genParticlePdgId, genParticleM1PdgId, genParticleStatus, genParticlem1px);
	
    cout << ""<< endl;
    for (int p = 0; p < genParticle; p++) {
      cout << genParticlePdgId[p];
      cout << "           ";
      cout << genParticleM1PdgId[p]<<endl;
    } 
    cout << ""<< endl;
	  
    //begin gen level comparison
    //find indices for final visible jets ,and their pair produced parent. 
    vector<int> gen_final_ind ;
    vector<int> gen_parent_ind ; 
	
    vector<TLorentzVector> gen_particles ;

	
    for (int gp = 0; gp < genParticle; gp++) {
      TLorentzVector v_dummy ;
      gen_particles.push_back(v_dummy);
      gen_particles[gp].SetPxPyPzE(genParticlePx[gp], genParticlePy[gp], genParticlePz[gp], genParticleE[gp]);
    }
	
	
	  
    vector<int> parent(genParticle, -1);
    for (int gp = 0; gp < genParticle; gp++) {
		
      if (parent[gp] != -1) continue;      //skip if this particle parent value is filled already
		
		
      //if its one of first 2 particles, it is considered one of colliding particles
      if (gp == 0 || gp == 1){
	parent[gp] = -2.0;
	continue;
      }
		
      //finds all indices of mothers where index is less than gp index, and also where more than 1 other daughter are not already assigned to the mother
      vector<int> possible_mothers;
      for (int ind = 0; ind < gp; ind++) {
	if (genParticlePdgId[ind]==genParticleM1PdgId[gp]&&(count(parent.begin(), parent.end(), ind) < 2 )) possible_mothers.push_back(ind);
      }
		
		
      //if there arent any valid mothers, this should be one of initial 2 colliding particles
      if (possible_mothers.size() == 0) {									
	parent[gp] = -2.0;		
	continue;
      }
		
		
      //if there is one possible mother, that is assigned as the parent
      if(possible_mothers.size() == 1){
	parent[gp] = possible_mothers[0];   
	continue;
      }
		
		
      //get list of indices of other possible daughter candidates for that mother particle type, which have not been assigned to a parent already
      vector<int> others_with_same_mother; 
      for (int ind = possible_mothers.back(); ind < genParticle; ind ++) {
	if ((genParticleM1PdgId[ind] == genParticleM1PdgId[gp]) && (parent[ind] == -1)) others_with_same_mother.push_back(ind); 
      }
      others_with_same_mother.erase( remove(others_with_same_mother.begin(), others_with_same_mother.end(), gp), others_with_same_mother.end() ); //remove gp from this list
		
		
		
      //now want to go through every partner to partner with gp, minimize [(gp+partner) distance from 1 mother candidate ] + [(all other daughter candidates) distance from (all other mother candidates)]
		
		
      TLorentzVector total_with_same_mother;  //initialize then build a vector including all non-gp daughter candidates
      for (int smi =0; smi < others_with_same_mother.size(); smi++) {
	total_with_same_mother = total_with_same_mother + gen_particles[others_with_same_mother[smi]] ;
      }
      TLorentzVector total_mother_cand;
      for (int mi = 0; mi < possible_mothers.size(); mi++) {
	total_mother_cand = total_mother_cand + gen_particles[possible_mothers[mi]];
      }
		
		
      double minR = 100.0;
      vector<int> op_pairing (2, -1) ;
      for (int smi = 0; smi < others_with_same_mother.size(); smi++) {
	TLorentzVector V1 = gen_particles[gp] + gen_particles[others_with_same_mother[smi]];
	TLorentzVector V2 = total_with_same_mother - gen_particles[others_with_same_mother[smi]];
	for (int mi = 0; mi < possible_mothers.size(); mi++) {
	  TLorentzVector M1 = gen_particles[possible_mothers[mi]];
	  TLorentzVector M2 = total_mother_cand - gen_particles[possible_mothers[mi]];
	  double currentR = D_between_vectors(V1, M1) + D_between_vectors(V2, M2); 
	  if (currentR < minR){
	    minR = currentR;
	    op_pairing[0] = others_with_same_mother[smi];
	    op_pairing[1] = possible_mothers[mi] ;
	  }
	}
      }
		

		
      parent[gp] = op_pairing[1];
      parent[op_pairing[0]] = op_pairing[1];
		

		
		
    }
	
    for (int p = 0; p < parent.size(); p++) {
      cout << parent[p];
      cout << " ";
    }
	cout << "" <<endl;

	  
	 
	  
	//find hemisphere parent
	vector<int> pp_parent(parent.size(), -1);
	for (int ppi = 0; ppi < pp_parent.size(); ppi++) {
	  if (parent[ppi] == -2){
	    pp_parent[ppi] = -2;
	    continue;
	  }
		
	  if (parent[parent[ppi]] == -2){
	    pp_parent[ppi] = ppi;
	    continue;
	  }
		
	  pp_parent[ppi] = pp_parent[parent[ppi]];
		
	}  
	
	  
	
	for (int p = 0; p < parent.size(); p++) {
	  cout << pp_parent[p];
	  cout << " ";
	}
	  
	cout << "" <<endl;
	cout << "" <<endl;	  

  }

  cout << event_counter << endl; ;
  cout << "test1" << endl; 	
	
  // full event TTree
  outTree->Write();
	
  file->cd();

	
  file->Close();

}


double CMSRazor13TeV::D_between_vectors(const TLorentzVector V1, const TLorentzVector V2){
	if ((V1.Px() * V2.Px() + V1.Py() * V2.Py() + V1.Pz() * V2.Pz()) == 0) return 0; 
	if ((V1.Px() * V2.Px() + V1.Py() * V2.Py()) == 0) return 50.0;
	
	double eta1 = V1.Eta();
	double eta2 = V2.Eta();
	double phi1 = V1.Phi();
	double phi2 = V2.Phi();
	return pow(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2) , 0.5);
};


vector<fastjet::PseudoJet> CMSRazor13TeV::improper_Exclusive_Jets(const fastjet::ClusterSequence cs_E){
	//vector<fastjet::PseudoJet> input_vecs ;
	//cs_C_E.transfer_input_jets(input_vecs);
	//vector<fastjet::PseudoJet> exc_jets;
	double dcut = 0.02;
	int jetN = 0;
	vector<fastjet::PseudoJet> fin_jets;

	while ((jetN != 2) && (dcut < 7.0)) {
		fin_jets.clear();
		fin_jets = cs_E.exclusive_jets(dcut);
		jetN = cs_E.n_exclusive_jets(dcut);
		dcut = dcut + 0.1;
	}
	
	return fin_jets;
};
