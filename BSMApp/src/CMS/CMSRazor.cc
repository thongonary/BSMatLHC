#include <string>
#include <iostream>
#include <TTree.h>
#include "CMS/CMSRazor.hh"
#include "CMS/ParticleInfo.hh"
#include <fastjet/tools/Pruner.hh>
#include <vector>
//for time:
#include <time.h>

//Use clock() to return the current time.  

CMSRazor::CMSRazor(TTree *tree, double Lumi, string analysis) : CMSReco(tree) {
  _Lumi = Lumi;
  _statTools = new StatTools(-99);
  _analysis = analysis;
}

CMSRazor::~CMSRazor(){
}

//using namespace std;
int event_counter = 0;

void CMSRazor::SetSqrts(double sqrts) {
	_sqrts = sqrts;
}

// loop over events - real analysis
void CMSRazor::Loop(string outFileName) {

	if(fChain == 0) return;

	// saving info about the PFAK04 jets

	double PFJets_phi[100];
	double PFJets_eta[100];
	double PFJets_pT[100];

	double MR_Default, RSQ_Default, MR_KT_Jets, RSQ_KT_Jets;
	double MRz_Default, MRz_KT_Jets, RSQz_Default, RSQz_KT_Jets;

	//newly added

	//loop through cone sizes algorithm
	double MRz_newjetakt, MR_newjetakt, RSQ_newjetakt, RSQz_newjetakt;
	double MRz_newjetkt, MR_newjetkt, RSQ_newjetkt, RSQz_newjetkt;
	double MRz_newjetcam, MR_newjetcam, RSQ_newjetcam, RSQz_newjetcam;

	//MRT values
	double MRt_Default, MRt_KT_Jets;
        double MRt_newjetcam, MRt_newjetkt, MRt_newjetakt;

	//inclusives jet size
	double jetlenaktkt, jetlenaktakt, jetlenaktcam;
	//ofstream outfile;
	double i;
	double count1 = 1;
	
	//hemisphere properties
	double Default_mass, Default_px, Default_py, Default_pz, Default_rapidity, Default_phi;
	double KT_Jets_mass, KT_Jets_px, KT_Jets_py, KT_Jets_pz, KT_Jets_rapidity, KT_Jets_phi;
        double Default_mass_2, Default_px_2, Default_py_2, Default_pz_2, Default_rapidity_2, Default_phi_2;
        double KT_Jets_mass_2, KT_Jets_px_2, KT_Jets_py_2, KT_Jets_pz_2, KT_Jets_rapidity_2, KT_Jets_phi_2;

	double newjetcam_px, newjetcam_px_2;
	double newjetcam_py, newjetcam_py_2;
	double newjetcam_mass, newjetcam_mass_2;
        double newjetkt_px, newjetkt_px_2;
	double newjetkt_py, newjetkt_py_2;
	double newjetkt_mass, newjetkt_mass_2;
        double newjetakt_px, newjetakt_px_2;
        double newjetakt_py, newjetakt_py_2;
	double newjetakt_mass, newjetakt_mass_2;
	double Default_E, KT_Jets_E, newjetcam_E;

	int BOX_NUM;
	double W_EFF;
	int SizesList[7];
	
	// Open Output file
	TFile *file = new TFile(outFileName.c_str(),"UPDATE"); 	
	
	cout << "writing to ";
	cout << outFileName << endl;
		
	TTree* outTree = new TTree("RazorInclusive","RazorInclusive");
	
	outTree->Branch("MR_Default", &MR_Default, "MR_Default/D");
	outTree->Branch("MRz_Default", &MRz_Default, "MRz_Default/D");
	outTree->Branch("RSQ_Default", &RSQ_Default, "RSQ_Default/D");
	outTree->Branch("RSQz_Default", &RSQz_Default, "RSQz_Default/D");
	outTree->Branch("MRt_Default", &MRt_Default, "MRt_Default/D");

	
	outTree->Branch("MR_KT_Jets", &MR_KT_Jets, "MR_KT_Jets/D");
	outTree->Branch("MRz_KT_Jets", &MRz_KT_Jets, "MRz_KT_Jets/D");
	outTree->Branch("RSQ_KT_Jets", &RSQ_KT_Jets, "RSQ_KT_Jets/D");
	outTree->Branch("RSQz_KT_Jets", &RSQz_KT_Jets, "RSQz_KT_Jets/D");
	outTree->Branch("MRt_KT_Jets", &MRt_KT_Jets, "MRt_KT_Jets/D");

	
	outTree->Branch("MR_newjetakt", &MR_newjetakt, "MR_newjetakt/D");
	outTree->Branch("MRz_newjetakt", &MRz_newjetakt, "MRz_newjetakt/D");
	outTree->Branch("RSQ_newjetakt", &RSQ_newjetakt, "RSQ_newjetakt/D");
	outTree->Branch("RSQz_newjetakt", &RSQz_newjetakt, "RSQz_newjetakt/D");
	outTree->Branch("MRt_newjetakt", &MRt_newjetakt, "MRt_newjetakt/D");

	outTree->Branch("MR_newjetkt", &MR_newjetkt, "MR_newjetkt/D");
	outTree->Branch("MRz_newjetkt", &MRz_newjetkt, "MRz_newjetkt/D");
	outTree->Branch("RSQ_newjetkt", &RSQ_newjetkt, "RSQ_newjetkt/D");
	outTree->Branch("RSQz_newjetkt", &RSQz_newjetkt, "RSQz_newjetkt/D");
	outTree->Branch("MRt_newjetkt", &MRt_newjetkt, "MRt_newjetkt/D");

	outTree->Branch("MR_newjetcam", &MR_newjetcam, "MR_newjetcam/D");
	outTree->Branch("MRz_newjetcam", &MRz_newjetcam, "MRz_newjetcam/D");
	outTree->Branch("RSQ_newjetcam", &RSQ_newjetcam, "RSQ_newjetcam/D");
	outTree->Branch("RSQz_newjetcam", &RSQz_newjetcam, "RSQz_newjetcam/D");
	outTree->Branch("MRt_newjetcam", &MRt_newjetcam, "MRt_newjetcam/D");

	outTree->Branch("SizesList", SizesList, "SizesList[8]/I");
 	outTree->Branch("PFJets_phi", PFJets_phi, "PFJets_phi[100]/D");
 	outTree->Branch("PFJets_eta", PFJets_eta, "PFJets_eta[100]/D");
 	outTree->Branch("PFJets_pT", PFJets_pT, "PFJets_pT[100]/D");

	//Jet information
	
	outTree->Branch("Default_mass", &Default_mass, "Default_mass/D");
	outTree->Branch("Default_px", &Default_px, "Default_px/D");
	outTree->Branch("Default_py", &Default_py, "Default_py/D");
	outTree->Branch("Default_pz", &Default_pz, "Default_pz/D");
	outTree->Branch("Default_rapidity", &Default_rapidity, "Default_rapidity/D");
	outTree->Branch("Default_phi", &Default_phi, "Default_phi/D");

	outTree->Branch("Default_mass_2", &Default_mass_2, "Default_mass_2/D");
	outTree->Branch("Default_px_2", &Default_px_2, "Default_px_2/D");
	outTree->Branch("Default_py_2", &Default_py_2, "Default_py_2/D");
	outTree->Branch("Default_pz_2", &Default_pz_2, "Default_pz_2/D");
	outTree->Branch("Default_rapidity_2", &Default_rapidity_2, "Default_rapidity_2/D");
	outTree->Branch("Default_phi_2", &Default_phi_2, "Default_phi_2/D");
	
	outTree->Branch("KT_Jets_mass", &KT_Jets_mass, "KT_Jets_mass/D");
	outTree->Branch("KT_Jets_px", &KT_Jets_px, "KT_Jets_px/D");
	outTree->Branch("KT_Jets_py", &KT_Jets_py, "KT_Jets_py/D");
	outTree->Branch("KT_Jets_pz", &KT_Jets_pz, "KT_Jets_pz/D");
	outTree->Branch("KT_Jets_rapidity", &KT_Jets_rapidity, "KT_Jets_rapidity/D");
	outTree->Branch("KT_Jets_phi", &KT_Jets_phi, "KT_Jets_phi/D");

	outTree->Branch("KT_Jets_mass_2", &KT_Jets_mass_2, "KT_Jets_mass_2/D");
	outTree->Branch("KT_Jets_px_2", &KT_Jets_px_2, "KT_Jets_px_2/D");
	outTree->Branch("KT_Jets_py_2", &KT_Jets_py_2, "KT_Jets_py_2/D");
	outTree->Branch("KT_Jets_pz_2", &KT_Jets_pz_2, "KT_Jets_pz_2/D");
	outTree->Branch("KT_Jets_rapidity_2", &KT_Jets_rapidity_2, "KT_Jets_rapidity_2/D");
	outTree->Branch("KT_Jets_phi_2", &KT_Jets_phi_2, "KT_Jets_phi_2/D");
	
	outTree->Branch("jetlenaktkt", &jetlenaktkt, "jetlenaktkt/D");

	outTree->Branch("newjetcam_mass", &newjetcam_mass, "newjetcam_mass/D");
	outTree->Branch("newjetcam_mass_2", &newjetcam_mass_2, "newjetcam_mass_2/D");

	outTree->Branch("newjetkt_mass", &newjetkt_mass, "newjetkt_mass/D");
	outTree->Branch("newjetkt_mass_2", &newjetkt_mass_2, "newjetkt_mass_2/D");
	outTree->Branch("newjetakt_mass", &newjetakt_mass, "newjetakt_mass/D");
	outTree->Branch("newjetakt_mass_2", &newjetakt_mass_2, "newjetakt_mass_2/D");
	outTree->Branch("Default_E", &Default_E, "Default_E/D");
	outTree->Branch("KT_Jets_E", &KT_Jets_E, "KT_Jets_E/D");
	outTree->Branch("newjetcam_E", &newjetcam_E, "newjetcam_E/D");

	outTree->Branch("newjetcam_px", &newjetcam_px, "newjetcam_px/D");
	outTree->Branch("newjetcam_py", &newjetcam_py, "newjetcam_py/D");
	outTree->Branch("newjetcam_px_2", &newjetcam_px_2, "newjetcam_px_2/D");
	outTree->Branch("newjetcam_py_2", &newjetcam_py_2, "newjetcam_py_2/D");

	outTree->Branch("newjetkt_px", &newjetkt_px, "newjetkt_px/D");
	outTree->Branch("newjetkt_py", &newjetkt_py, "newjetkt_py/D");
        outTree->Branch("newjetkt_px_2", &newjetkt_px_2, "newjetkt_px_2/D");
	outTree->Branch("newjetkt_py_2", &newjetkt_py_2, "newjetkt_py_2/D");

        outTree->Branch("newjetakt_px", &newjetakt_px, "newjetakt_px/D");
        outTree->Branch("newjetakt_py", &newjetakt_py, "newjetakt_py/D");
        outTree->Branch("newjetakt_px_2", &newjetakt_px_2, "newjetakt_px_2/D");
        outTree->Branch("newjetakt_py_2", &newjetakt_py_2, "newjetakt_py_2/D");

	outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
	outTree->Branch("W_EFF", &W_EFF, "W_EFF/D");
		
	double xedge[17] = {300, 350, 400.,450.,500.,550.,600.,650.,700.,800.,900.,1000.,1200.,1600.,2000.,2800.,3500.};
	double yedge[6] = {0.11,0.18,0.20,0.30,0.40,0.50};
	TH2D* pdfHad = new TH2D("pdfHad","pdfHad",16,xedge,5,yedge);
	TH2D* pdfMuMu = new TH2D("pdfMuMu","pdfMuEle",16,xedge,5,yedge);
	TH2D* pdfMuEle = new TH2D("pdfMuEle","pdfMuEle",16,xedge,5,yedge);
	TH2D* pdfMu = new TH2D("pdfMu","pdfMu",16,xedge,5,yedge);
	TH2D* pdfEleEle = new TH2D("pdfEleEle","pdfEleEle",16,xedge,5,yedge);
	TH2D* pdfEle = new TH2D("pdfEle","pdfEle",16,xedge,5,yedge);
	
	// loop over entries
	Long64_t nbytes = 0, nb = 0;
	Long64_t nentries = fChain->GetEntries();
	std::cout << "Number of entries = " << nentries << std::endl;
	
	// set the by-event weight
	W_EFF = 1./nentries;
	
	
  for (Long64_t jentry=0; jentry<nentries;jentry+=1) {

    if(verbose) cout << "new event" << endl;

    // clean physics-objects blocks
      CleanEvent();

    // get new event
      Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
	//if (jentry > 10000) continue; //use for gen level comparison coding
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
  	
    // Build the event at generator level
    PFReco();
    vector<fastjet::PseudoJet> empty;
    vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty);
    
    bool ttbar_no_all_had = false; //true if you want to filter out all-hadronic events (no leptons)
    bool hasLepton = true;
    if (ttbar_no_all_had){
      hasLepton = false;
    }
    int particleID;
    if (ttbar_no_all_had){
      for (int k = 0; k < JetsConst.size(); k++){
       particleID = JetsConst[k].user_info<ParticleInfo>().pdg_id;
       if (abs(particleID) == 11 || abs(particleID) == 13 || abs(particleID) == 15) hasLepton = true;
      }
    }
    if (!hasLepton) continue;

    i = 0;
	  
    // wide jets

    /*
    fastjet::JetDefinition CA08_def(fastjet::cambridge_algorithm, 0.8);
    fastjet::ClusterSequence pfCA08ClusterSequence = JetMaker(JetsConst, CA08_def);
    vector<fastjet::PseudoJet> pfCA08 = SelectByAcceptance(fastjet::sorted_by_pt(pfCA08ClusterSequence.inclusive_jets()),40., 2.4);
    fastjet::Pruner pruner(CA08_def, 0.1, 0.25);
    */

    // narrow jets
    fastjet::JetDefinition AK04_def(fastjet::antikt_algorithm, 0.4);
    fastjet::ClusterSequence pfAK04ClusterSequence = JetMaker(JetsConst, AK04_def);
    vector<fastjet::PseudoJet> pfAK04 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK04ClusterSequence.inclusive_jets()),40., 2.4);
  	//vector<fastjet::PseudoJet> pfAK04 = SelectByAcceptance(fastjet::sorted_by_pt(pfAK04ClusterSequence.inclusive_jets()),18., 2.4);
    i = 0;

    if(pfAK04.size()<2) continue;
    event_counter = event_counter + 1 ;

    //cout << "Event: " << event_counter << endl;

    //saving variables for PFAK04 jets
    for(unsigned k=0; k<pfAK04.size(); k++){
      pfAK04[k].set_user_index(k); //label jet with index
      PFJets_phi[k] = pfAK04[k].phi();
      PFJets_eta[k] = pfAK04[k].eta();
      PFJets_pT[k] = pfAK04[k].pt();
    }

    GenMET();
    PFMET = genMET;

    // Ele reco: WP80 and WP95
    EleReco();

    //cout << "Ele Done" << endl;

    // Mu reco: Tight and Loose
    MuReco();
    
    // 1a) cluster hemispheres using kT in exclusive mode, using jets as ingredients
    // some backward compatibility test will be needed here

    double cone_size = 1.57; //did 0.8, 1.0  
	  
    fastjet::ClusterSequence cs1(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, cone_size));
    vector<TLorentzVector> hem_KT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs1.exclusive_jets(2))); //this is from original, written correctly
    //vector<TLorentzVector> hem_KT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs1.inclusive_jets()));

    	  
    //----------------inclusive mode used to find jets thrown out---------------------
    fastjet::ClusterSequence cs101(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, cone_size));
    vector<TLorentzVector> hem_KT_Jets_inc = ConvertTo4Vector(fastjet::sorted_by_pt(cs101.inclusive_jets()));

    jetlenaktkt = hem_KT_Jets_inc.size();


    //---------------------------new clustering algorithms that run through cone size--------------------------------//

    //New idea for clustering.  Re-run anti-kt inclusively and increase cone size to lower final jet numer
    //Do the same for kt and cam also
    double jetlen = pfAK04.size(); 
    double conesize1 = 0.2;
    double conestep = 0.05;
    bool changed = false;
    vector<TLorentzVector> newjetakt;
    vector<fastjet::PseudoJet> final_akt_jets;
    while (jetlen > 1){
      fastjet::ClusterSequence cs11(pfAK04, fastjet::JetDefinition(fastjet::antikt_algorithm, conesize1));
      newjetakt = ConvertTo4Vector(fastjet::sorted_by_pt(cs11.inclusive_jets()));
      jetlen = newjetakt.size();
      if (jetlen > 1)
	conesize1 = conesize1 + conestep;
      changed = true;
    }
    if (changed) {
      fastjet::ClusterSequence cs11(pfAK04, fastjet::JetDefinition(fastjet::antikt_algorithm, conesize1 - conestep));
      newjetakt = ConvertTo4Vector(fastjet::sorted_by_pt(cs11.inclusive_jets()));
      final_akt_jets = fastjet::sorted_by_pt(cs11.inclusive_jets());
      /*for (unsigned k=0; k<final_akt_jets.size();k++){                        
	cout << "jet "<< k << ": pt "<<final_akt_jets[k].pt() <<", eta  " << final_akt_jets[i].eta() << ", phi " << final_akt_jets[k].phi() << endl;                            
	vector<fastjet::PseudoJet> constituents = final_akt_jets[k].constituents();                                                                                             
        for (unsigned l=0; l<constituents.size();l++){                                                                                                                             
	cout <<" constituent" << l <<"'s index: "<< constituents[l].user_index() << endl;                                                                                 
	}                                                                                                                                                                 
	}*/
    }
    else {
      newjetakt = ConvertTo4Vector(pfAK04);
      final_akt_jets = pfAK04;
    }

    jetlen = pfAK04.size();
    conesize1 = 0.2;
    conestep = 0.05;
    changed = false;
    vector<TLorentzVector> newjetcam;
    vector<fastjet::PseudoJet> final_jets;

    while (jetlen > 1){
      fastjet::ClusterSequence cs1101(pfAK04, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1));
      newjetcam = ConvertTo4Vector(fastjet::sorted_by_pt(cs1101.inclusive_jets()));
      jetlen = newjetcam.size();
      if (jetlen > 1)
        conesize1 = conesize1 + conestep;
      changed = true;
    }
    if (changed) {
      fastjet::ClusterSequence cs1101(pfAK04, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1 - conestep));
      newjetcam = ConvertTo4Vector(fastjet::sorted_by_pt(cs1101.inclusive_jets()));
      final_jets = fastjet::sorted_by_pt(cs1101.inclusive_jets()); //save pseudojets
      /*for (unsigned k=0; k<final_jets.size();k++){
	cout << "jet "<< k << ": pt "<<final_jets[k].pt() <<", eta  " << final_jets[i].eta() << ", phi " << final_jets[k].phi() << endl;
	vector<fastjet::PseudoJet> constituents = final_jets[k].constituents();
	for (unsigned l=0; l<constituents.size();l++){
		  cout <<" constituent" << l <<"'s index: "<< constituents[l].user_index() << endl;
		  }
	}*/
    }
    else {
      newjetcam = ConvertTo4Vector(pfAK04);
      final_jets = pfAK04;
    }

    //kt
    jetlen = pfAK04.size(); 
    conesize1 = 0.2;
    changed = false;
    vector<TLorentzVector> newjetkt;
    while (jetlen > 1){
      fastjet::ClusterSequence cs12(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, conesize1));
      newjetkt = ConvertTo4Vector(fastjet::sorted_by_pt(cs12.inclusive_jets()));
      jetlen = newjetkt.size();
      if (jetlen > 1)
	conesize1 = conesize1 + conestep;
      changed = true;
    }
    if (changed){
      fastjet::ClusterSequence cs12(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, conesize1 - conestep));
      newjetkt = ConvertTo4Vector(fastjet::sorted_by_pt(cs12.inclusive_jets()));
    }
    else newjetkt = ConvertTo4Vector(pfAK04);

    //End of new idea of clustering
    

    //Clusters particles off of previous clustering
    /*
    jetlen = pfAK04.size();
    conesize1 = 0.2;
    changed = false;
    vector<fastjet::PseudoJet> Continualjet = pfAK04;
    while (jetlen > 2){
      fastjet::ClusterSequence cs1000(Continualjet, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1));
      Continualjet = fastjet::sorted_by_pt(cs1000.inclusive_jets()); 
      jetlen = Continualjet.size();
      conesize1 += 0.1;
      cout << jetlen << " " << conesize1 << endl;
    }
    */

    
    // 1b) traditional hemispheres
    CMSHemisphere* myHem = new CMSHemisphere(ConvertTo4Vector(pfAK04));
    myHem->CombineMinMass();
    vector<TLorentzVector> hem_Default = myHem->GetHemispheres();
    vector<int> Temporary = myHem->GetHem1Constituents();
    for (int k=0; k < Temporary.size(); k++){
      cout << "Megajet jet: " << Temporary[k] << endl;
    }
    delete myHem; 
	
    //initialize  
    MR_KT_Jets = -9999. ;
    RSQ_KT_Jets = -9999. ;
    MR_Default = -9999. ;
    RSQ_Default = -9999. ;
    MRz_KT_Jets = -9999. ;
    MRz_Default = -9999. ;
    RSQz_KT_Jets = -9999. ;
    RSQz_Default = -9999. ;

 // Edward's new algorithm

    MR_newjetakt = -9999. ;
    MRz_newjetakt = -9999. ;
    RSQ_newjetakt = -9999. ;
    RSQz_newjetakt = -9999. ;
    MR_newjetkt = -9999.;
    MRz_newjetkt = -9999. ;
    RSQ_newjetkt = -9999. ;
    RSQz_newjetkt = -9999. ;   
    MR_newjetcam = -9999.;
    MRz_newjetcam = -9999. ;
    RSQ_newjetcam = -9999. ;
    RSQz_newjetcam = -9999. ;    
 
    
    SizesList[0] = -9999.0;
    SizesList[1] = -9999.0;
    SizesList[2] = -9999.0;
    SizesList[3] = -9999.0;
    SizesList[4] = -9999.0;
    SizesList[5] = -9999.0;
    SizesList[6] = -9999.0;
    SizesList[7] = -9999.0;

    KT_Jets_mass = -9999.0;
    KT_Jets_px = -9999.0;
    KT_Jets_py = -9999.0;
    KT_Jets_pz = -9999.0;
    KT_Jets_phi = -9999.0;
    KT_Jets_rapidity = -9999.0;

    Default_mass = -9999.0;
    Default_px = -9999.0;
    Default_py = -9999.0;
    Default_pz = -9999.0;
    Default_phi = -9999.0;
    Default_rapidity = -9999.0;

    KT_Jets_mass_2 = -9999.0;
    KT_Jets_px_2 = -9999.0;
    KT_Jets_py_2 = -9999.0;
    KT_Jets_pz_2 = -9999.0;
    KT_Jets_phi_2 = -9999.0;
    KT_Jets_rapidity_2 = -9999.0;

    Default_mass_2 = -9999.0;
    Default_px_2 = -9999.0;
    Default_py_2 = -9999.0;
    Default_pz_2 = -9999.0;
    Default_phi_2 = -9999.0;
    Default_rapidity_2 = -9999.0;

    newjetcam_px = -9999.0;
    newjetcam_py = -9999.0;
    newjetcam_px_2 = -9999.0;
    newjetcam_py_2 = -9999.0;

    newjetakt_px = -9999.0;
    newjetakt_py = -9999.0;
    newjetakt_px_2 = -9999.0;
    newjetakt_py_2 = -9999.0;

    newjetkt_px = -9999.0;
    newjetkt_py = -9999.0;
    newjetkt_px_2 = -9999.0;
    newjetkt_py_2 = -9999.0;

    TLorentzVector j1;
    TLorentzVector j2;
 
    // 2a) compute new RSQ and MR vals

    if (hem_KT_Jets.size() > 1){
      j1 = hem_KT_Jets[0];
      j2 = hem_KT_Jets[1]; 
      KT_Jets_E = j1.E();
      KT_Jets_mass = j1.M();
      KT_Jets_px = j1.Px();
      KT_Jets_py = j1.Py();
      KT_Jets_pz = j1.Pz();
      KT_Jets_phi = j1.Phi();
      KT_Jets_rapidity = j1.Rapidity();
      KT_Jets_mass_2 = j2.M();
      KT_Jets_px_2 = j2.Px();
      KT_Jets_py_2 = j2.Py();
      KT_Jets_pz_2 = j2.Pz();
      KT_Jets_phi_2 = j2.Phi();
      KT_Jets_rapidity_2 = j2.Rapidity();
      MR_KT_Jets = CalcMR(j1, j2);
      MRt_KT_Jets = CalcMRT(j1, j2, PFMET);
      RSQ_KT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MR_KT_Jets/MR_KT_Jets;
      MRz_KT_Jets = CalcMR_zinvariant(j1, j2);
      RSQz_KT_Jets = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_KT_Jets/MRz_KT_Jets;
	  
    }
	  
    if (newjetakt.size() > 1) {
      j1 = newjetakt[0];
      j2 = newjetakt[1];  
      newjetakt_mass = j1.M();
      newjetakt_mass_2 = j2.M();
      newjetakt_px = j1.Px();
      newjetakt_py = j1.Py();
      newjetakt_px_2 = j2.Px();
      newjetakt_py_2 = j2.Py();
      MR_newjetakt = CalcMR(j1, j2);
      MRt_newjetakt = CalcMRT(j1, j2, PFMET);
      RSQ_newjetakt = pow(CalcMRT(j1, j2, PFMET),2.)/MR_newjetakt/MR_newjetakt;
      MRz_newjetakt = CalcMR_zinvariant(j1, j2);
      RSQz_newjetakt = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_newjetakt/MRz_newjetakt;
    }

    if (newjetkt.size() > 1) {
      j1 = newjetkt[0];
      j2 = newjetkt[1];
      newjetkt_mass = j1.M();
      newjetkt_mass_2 = j2.M();
      newjetkt_px = j1.Px();
      newjetkt_py = j1.Py();
      newjetkt_px_2 = j2.Px();
      newjetkt_py_2 = j2.Py();
      MR_newjetkt = CalcMR(j1, j2);
      MRt_newjetkt = CalcMRT(j1, j2, PFMET);
      RSQ_newjetkt = pow(CalcMRT(j1, j2, PFMET),2.)/MR_newjetkt/MR_newjetkt;
      MRz_newjetkt = CalcMR_zinvariant(j1, j2);
      RSQz_newjetkt = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_newjetkt/MRz_newjetkt;
    }

    if (newjetcam.size() > 1) {
      j1 = newjetcam[0];
      j2 = newjetcam[1];  
      newjetcam_E = j1.E();
      newjetcam_px = j1.Px();
      newjetcam_py = j1.Py();
      newjetcam_px_2 = j2.Px();
      newjetcam_py_2 = j2.Py();
      //cout << j1.Px() << " " << j1.Py() << " " << j2.Px() << " " << j2.Py() << endl;
      newjetcam_mass = j1.M();
      newjetcam_mass_2 = j2.M();
      MR_newjetcam = CalcMR(j1, j2);      
      MRt_newjetcam = CalcMRT(j1, j2, PFMET);
      RSQ_newjetcam = pow(CalcMRT(j1, j2, PFMET),2.)/MR_newjetcam/MR_newjetcam;
      MRz_newjetcam = CalcMR_zinvariant(j1, j2);
      RSQz_newjetcam = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_newjetcam/MRz_newjetcam;
    } 
    else cout<<"FAILED"<<endl;
    
    // 2b) compute traditional RSQ and MR (DEFAULT)
    j1 = hem_Default[0];
    j2 = hem_Default[1];
    Default_E = j1.E();
    Default_mass = j1.M();
    Default_px = j1.Px();
    Default_py = j1.Py();
    Default_pz = j1.Pz();
    Default_phi = j1.Phi();
    Default_rapidity = j1.Rapidity();
    Default_mass_2 = j2.M();
    Default_px_2 = j2.Px();
    Default_py_2 = j2.Py();
    Default_pz_2 = j2.Pz();
    Default_phi_2 = j2.Phi();
    Default_rapidity_2 = j2.Rapidity();
    MR_Default = CalcMR(j1, j2);
    MRt_Default = CalcMRT(j1, j2, PFMET);
    RSQ_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MR_Default/MR_Default;
    MRz_Default = CalcMR_zinvariant(j1, j2);
    RSQz_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_Default/MRz_Default;
	  

    //delete j1;
    //delete j2;
    
    
    SizesList[0] = hem_Default.size();
    SizesList[1] = hem_KT_Jets.size();
    // SizesList[2] = hem_AKT_Jets.size();
    //SizesList[3] = hem_CAM_Jets.size();
    //SizesList[4] = hem_KT_CAM.size();
    //SizesList[5] = hem_AKT_CAM.size();
    //SizesList[6] = hem_CAM_CAM.size();
    SizesList[7] = pfAK04.size(); //pre cluster pfak04 jet number
    //cout <<endl << "Number of PFJets:" << pfAK04.size()<<endl;
    
	  
    // Boxes
    BOX_NUM = 5; // Had by default
    if(MUELEBox()) BOX_NUM = 0;
    else if(MUMUBox()) BOX_NUM = 1;
    else if(ELEELEBox()) BOX_NUM = 2;
    else if(MUBox()) BOX_NUM = 3;
    else if(ELEBox()) BOX_NUM = 4;

    // write event in the tree
    outTree->Fill();
	

    // fill PDF histograms
    bool fillBox = SignalRegion(MR_Default, RSQ_Default, BOX_NUM);
    if(BOX_NUM == 0 && fillBox) pdfMuEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 1 && fillBox) pdfMuMu->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 2 && fillBox) pdfEleEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 3 && fillBox) pdfMu->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 4 && fillBox) pdfEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 5 && fillBox) pdfHad->Fill(MR_Default, RSQ_Default);
	  
	  
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
	
    GenReturn(genParticle, genParticlePx, genParticlePy, genParticlePz, genParticleE, genParticlePdgId, genParticleM1PdgId);
	
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
	
  cout << "test2" << endl; 	
  
  // eff TTree
  double effMuEle = pdfMuEle->Integral()/double(nentries);
  double effMuMu = pdfMuMu->Integral()/double(nentries);
  double effMu = pdfMu->Integral()/double(nentries);
  double effEleEle = pdfEleEle->Integral()/double(nentries);
  double effEle = pdfEle->Integral()/double(nentries);
  double effHad = pdfHad->Integral()/double(nentries);

  // normalize the PDFs
  if(pdfMuEle->Integral()>0)  pdfMuEle->Scale(1./pdfMuEle->Integral());
  if(pdfMuMu->Integral()>0)  pdfMuMu->Scale(1./pdfMuMu->Integral());
  if(pdfMu->Integral()>0)  pdfMu->Scale(1./pdfMu->Integral());
  if(pdfEle->Integral()>0)  pdfEle->Scale(1./pdfEle->Integral());
  if(pdfEleEle->Integral()>0)  pdfEleEle->Scale(1./pdfEleEle->Integral());
  if(pdfHad->Integral()>0)  pdfHad->Scale(1./pdfHad->Integral());

  // write the PDFs
  pdfMuEle->Write();  
  pdfMuMu->Write();
  pdfMu->Write();
  pdfEle->Write();
  pdfEleEle->Write();
  pdfHad->Write();

  cout << "test3" << endl; 
	
  char name[256];
  sprintf(name,"data/%s.root", _analysis.c_str());
  //  TH1D* xsecProb = XsecProb(pdfHad, effHad,name, 1000, 0., 1.);
  // Open Output file again 
  file->cd();
  double xsecULHad = 0.;//_statTools->FindUL(xsecProb, 0.95, 1.);

  TTree* effTree = new TTree("RazorInclusiveEfficiency","RazorInclusiveEfficiency");
  effTree->Branch("effMuEle", &effMuEle, "effMuEle/D");
  effTree->Branch("effMuMu", &effMuMu, "effMuMu/D");
  effTree->Branch("effEleEle", &effEleEle, "effEleEle/D");
  effTree->Branch("effMu", &effMu, "effMu/D");
  effTree->Branch("effEle", &effEle, "effEle/D");
  effTree->Branch("effHad", &effHad, "effHad/D");
  effTree->Branch("xsecULHad", &xsecULHad, "xsecULHad/D");
  effTree->Fill();
  effTree->Write();

  cout << "test4" << endl; 	
	
  //  xsecProb->Write();
  file->Close();

	
  cout << "test5";	
  //outfile.close();
}


bool CMSRazor::ELEELEBox() {
  int iEle = 0;
  for(int j=0; j< EleWP95.size(); j++) {
    if(EleWP95[j].Pt()>10.) iEle++;
  }
  return ELEBox()*(iEle>1? true : false);
}

bool CMSRazor::MUBox() {
  int iMu = 0;
  for(int i=0; i< TightMu.size(); i++) {
    if(TightMu[i].Pt()>12.) iMu++;
  }
  return (iMu>0 ? true : false);
}

bool CMSRazor::ELEBox() {
  int iEle = 0;
  for(int i=0; i< EleWP80.size(); i++) {
    if(EleWP80[i].Pt()>20.) iEle++;
  }
  return (iEle>0 ? true : false);
}

bool CMSRazor::MUELEBox() {
  return MUBox()*ELEBox();
}

bool CMSRazor::MUMUBox() {
  int iMu = 0;
  for(int i=0; i< LooseMu.size(); i++) {
    if(LooseMu[i].Pt()>10.) iMu++;
  }
  if(iMu<2) return false;
  iMu = 0;
  for(int i=0; i< TightMu.size(); i++) {
    if(TightMu[i].Pt()>15.) iMu++;
  }
  if(iMu<1) return false;
  return true;
}

bool CMSRazor::SignalRegion(double mr, double rsq, double ibox){
  if(ibox == 4 || ibox == 3) return SignalRegionLep(mr, rsq);
  else if(ibox == 5) return SignalRegionHad(mr, rsq);
  else if(ibox >=0 && ibox <=2) return SignalRegionDiLep(mr, rsq);
  cout <<"Error on CMSRazor::SignalRegion : invalid box number " << ibox << endl;
  return false;
}

bool CMSRazor::SignalRegionHad(double mr, double rsq){
  // tighter baseline cuts
  if(rsq<0.18) return false;
  if(mr<500.) return false;
  return SignalRegionLep(mr, rsq);
}

bool CMSRazor::SignalRegionLep(double mr, double rsq){
  if(rsq>0.5) return false;
  if(rsq<0.11) return false;
  if(mr>1000.) return true;
  if(rsq>0.2 && mr>650.) return true;
  if(rsq>0.3 && mr>450.) return true; 
  return false;
}

bool CMSRazor::SignalRegionDiLep(double mr, double rsq){
  if(rsq>0.5) return false;
  if(rsq<0.11) return false;
  if(mr>650.) return true;
  if(rsq>0.2 && mr>450.) return true;
  if(rsq>0.3 && mr>400.) return true; 
  return false;
}

TH1D* CMSRazor::XsecProb(TH2D* sigPdf, double eff, TString Filename, int ibin, double xmin, double xmax) {

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


//vector<int> CMSRazor::find_all(const vector<T> objects, item, low_bound, up_bound){
//	vector<int> indices;
//	for (int ind = low_bound; ind < up_bound + 1; ind ++) {
//		if (item == objects[ind]) indices.push_back(ind); 
//	}
//	return indices;
//
//};

double CMSRazor::D_between_vectors(const TLorentzVector V1, const TLorentzVector V2){
	if ((V1.Px() * V2.Px() + V1.Py() * V2.Py() + V1.Pz() * V2.Pz()) == 0) return 0; 
	if ((V1.Px() * V2.Px() + V1.Py() * V2.Py()) == 0) return 50.0;
	
	double eta1 = V1.Eta();
	double eta2 = V2.Eta();
	double phi1 = V1.Phi();
	double phi2 = V2.Phi();
	return pow(pow(eta1 - eta2, 2) + pow(phi1 - phi2, 2) , 0.5);
};


