#include <string>
#include <iostream>
#include <TTree.h>
#include "CMS/CMSRazor.hh"
#include "CMS/ParticleInfo.hh"
#include <fastjet/tools/Pruner.hh>
#include <vector>
#include <sstream>
#include <stdio.h>
//??


// SWITCHES:
// BOOLEAN ttbar_no_all_had
// BOOLEAN ttbar 

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
  int pflen = 100;
  double PFJets_phi[pflen];
  double PFJets_eta[pflen];
  double PFJets_pT[pflen];
  double PFJets_m[pflen];
  // index of PFAK04 jets in each original pair-produced particle hemisphere
  int gen_hem1[pflen];
  int gen_hem2[pflen];

  double MR_Default, RSQ_Default, MR_KT_Jets, RSQ_KT_Jets;
  double MRz_Default, MRz_KT_Jets, RSQz_Default, RSQz_KT_Jets;
  
  double MR_Georgi, RSQ_Georgi;
  double MRz_Georgi, RSQz_Georgi;

  //loop through cone sizes algorithm
  double MRz_newjetakt, MR_newjetakt, RSQ_newjetakt, RSQz_newjetakt;
  double MRz_newjetkt, MR_newjetkt, RSQ_newjetkt, RSQz_newjetkt;
  double MRz_newjetcam, MR_newjetcam, RSQ_newjetcam, RSQz_newjetcam;
  
  //MRT values
  double MRt_Default, MRt_KT_Jets, MRt_Georgi;
  double MRt_newjetcam, MRt_newjetkt, MRt_newjetakt;
  
  //inclusives jet size
  double jetlenaktkt, jetlenaktakt, jetlenaktcam;
  
  //ofstream outfile;
  double i;
  double count1 = 1;
  
  //hemisphere properties
  double Default_mass, Default_px, Default_py, Default_pz, Default_eta, Default_phi;
  double KT_Jets_mass, KT_Jets_px, KT_Jets_py, KT_Jets_pz, KT_Jets_eta, KT_Jets_phi;
  double Default_mass_2, Default_px_2, Default_py_2, Default_pz_2, Default_eta_2, Default_phi_2, Default_R;
  double KT_Jets_mass_2, KT_Jets_px_2, KT_Jets_py_2, KT_Jets_pz_2, KT_Jets_eta_2, KT_Jets_phi_2;


 double Georgi_mass, Georgi_px, Georgi_py, Georgi_pz, Georgi_eta, Georgi_phi;  
  double Georgi_mass_2, Georgi_px_2, Georgi_py_2, Georgi_pz_2, Georgi_eta_2, Georgi_phi_2, Georgi_R;
  double Georgi_eta_diff, Georgi_phi_diff;
  double Default_eta_diff, KT_Jets_eta_diff, newjetcam_eta_diff;
  double Default_phi_diff, KT_Jets_phi_diff, newjetcam_phi_diff;
  
  double newjetcam_px, newjetcam_px_2;
  double newjetcam_py, newjetcam_py_2;
  double newjetcam_mass, newjetcam_mass_2;
  double newjetkt_px, newjetkt_px_2;
  double newjetkt_py, newjetkt_py_2;
  double newjetcam_phi, newjetcam_phi_2;
  double newjetcam_eta, newjetcam_eta_2;
  
  
  double newjetkt_mass, newjetkt_mass_2;
  double newjetakt_px, newjetakt_px_2;
  double newjetakt_py, newjetakt_py_2;
  double newjetakt_mass, newjetakt_mass_2;
  double Default_E, KT_Jets_E, newjetcam_E;
  double newjetcam_R;
  
  double MR_Cont, MRz_Cont, RSQ_Cont, RSQz_Cont, MRT_Cont, Cont_px, Cont_py, Cont_pz, Cont_mass;
  double Cont_px_2, Cont_py_2, Cont_pz_2, Cont_mass_2, Cont_phi_diff, Cont_eta_diff;
  double pt_Gen_lead, pt_Gen_sub;
  double Gen_lead_phi, Gen_lead_eta, Gen_sub_phi, Gen_sub_eta;
  int MR_len = 20000;
  //double All_MR[MR_len];
  double min_MR, max_MR;
  //double All_RSQ[MR_len];
  double min_RSQ, max_RSQ;
  double parent_mass1, parent_mass2;
  double MR_estimator1, MR_estimator2, MR_estimator3,  MR_estimator4, MR_perfect_cluster, RSQ_perfect_cluster, MR_random_cluster;
  double perfect_mass_1, perfect_mass_2;
  //saving constituents of hemispheres
  int Default_hem1_csts[pflen], Default_hem2_csts[pflen];
  int newjetcam_hem1_csts[pflen], newjetcam_hem2_csts[pflen];
  int KT_Jets_hem1_csts[pflen], KT_Jets_hem2_csts[pflen];
  int Georgi_hem1_csts[pflen], Georgi_hem2_csts[pflen];
  
  //saving cone size
  double finalconesize;
  
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
  
  outTree->Branch("MR_Georgi", &MR_Georgi, "MR_Georgi/D");
  outTree->Branch("MRz_Georgi", &MRz_Georgi, "MRz_Georgi/D");
  outTree->Branch("RSQ_Georgi", &RSQ_Georgi, "RSQ_Georgi/D");
  outTree->Branch("RSQz_Georgi", &RSQz_Georgi, "RSQz_Georgi/D");
  outTree->Branch("MRt_Georgi", &MRt_Georgi, "MRt_Georgi/D");
  
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
  outTree->Branch("PFJets_m", PFJets_m, "PFJets_m[100]/D");
  
  outTree->Branch("gen_hem1", gen_hem1, "gen_hem1[100]/I");
  outTree->Branch("gen_hem2", gen_hem2, "gen_hem2[100]/I");

  //Jet information
  
  outTree->Branch("Default_mass", &Default_mass, "Default_mass/D");
  outTree->Branch("Default_px", &Default_px, "Default_px/D");
  outTree->Branch("Default_py", &Default_py, "Default_py/D");
  outTree->Branch("Default_pz", &Default_pz, "Default_pz/D");
  outTree->Branch("Default_eta", &Default_eta, "Default_eta/D");
  outTree->Branch("Default_phi", &Default_phi, "Default_phi/D");
  
  outTree->Branch("Default_mass_2", &Default_mass_2, "Default_mass_2/D");
  outTree->Branch("Default_px_2", &Default_px_2, "Default_px_2/D");
  outTree->Branch("Default_py_2", &Default_py_2, "Default_py_2/D");
  outTree->Branch("Default_pz_2", &Default_pz_2, "Default_pz_2/D");
  outTree->Branch("Default_eta_2", &Default_eta_2, "Default_eta_2/D");
  outTree->Branch("Default_phi_2", &Default_phi_2, "Default_phi_2/D");
  outTree->Branch("Default_phi_diff", &Default_phi_diff, "Default_phi_diff/D");
  outTree->Branch("Default_eta_diff", &Default_eta_diff, "Default_eta_diff/D");
  outTree->Branch("Default_R", &Default_R, "Default_R/D");

  outTree->Branch("Georgi_mass", &Georgi_mass, "Georgi_mass/D");
  outTree->Branch("Georgi_px", &Georgi_px, "Georgi_px/D");
  outTree->Branch("Georgi_py", &Georgi_py, "Georgi_py/D");
  outTree->Branch("Georgi_pz", &Georgi_pz, "Georgi_pz/D");
  outTree->Branch("Georgi_eta", &Georgi_eta, "Georgi_eta/D");
  outTree->Branch("Georgi_phi", &Georgi_phi, "Georgi_phi/D");
  
  outTree->Branch("Georgi_mass_2", &Georgi_mass_2, "Georgi_mass_2/D");
  outTree->Branch("Georgi_px_2", &Georgi_px_2, "Georgi_px_2/D");
  outTree->Branch("Georgi_py_2", &Georgi_py_2, "Georgi_py_2/D");
  outTree->Branch("Georgi_pz_2", &Georgi_pz_2, "Georgi_pz_2/D");
  outTree->Branch("Georgi_eta_2", &Georgi_eta_2, "Georgi_eta_2/D");
  outTree->Branch("Georgi_phi_2", &Georgi_phi_2, "Georgi_phi_2/D");
  outTree->Branch("Georgi_phi_diff", &Georgi_phi_diff, "Georgi_phi_diff/D");
  outTree->Branch("Georgi_eta_diff", &Georgi_eta_diff, "Georgi_eta_diff/D");
  outTree->Branch("Georgi_R", &Georgi_R, "Georgi_R/D");
  
  
  outTree->Branch("KT_Jets_mass", &KT_Jets_mass, "KT_Jets_mass/D");
  outTree->Branch("KT_Jets_px", &KT_Jets_px, "KT_Jets_px/D");
  outTree->Branch("KT_Jets_py", &KT_Jets_py, "KT_Jets_py/D");
  outTree->Branch("KT_Jets_pz", &KT_Jets_pz, "KT_Jets_pz/D");
  outTree->Branch("KT_Jets_eta", &KT_Jets_eta, "KT_Jets_eta/D");
  outTree->Branch("KT_Jets_phi", &KT_Jets_phi, "KT_Jets_phi/D");
  
  outTree->Branch("KT_Jets_mass_2", &KT_Jets_mass_2, "KT_Jets_mass_2/D");
  outTree->Branch("KT_Jets_px_2", &KT_Jets_px_2, "KT_Jets_px_2/D");
  outTree->Branch("KT_Jets_py_2", &KT_Jets_py_2, "KT_Jets_py_2/D");
  outTree->Branch("KT_Jets_pz_2", &KT_Jets_pz_2, "KT_Jets_pz_2/D");
  outTree->Branch("KT_Jets_eta_2", &KT_Jets_eta_2, "KT_Jets_eta_2/D");
  outTree->Branch("KT_Jets_phi_2", &KT_Jets_phi_2, "KT_Jets_phi_2/D");
  outTree->Branch("KT_Jets_phi_diff", &KT_Jets_phi_diff, "KT_Jets_phi_diff/D");
  outTree->Branch("KT_Jets_eta_diff", &KT_Jets_eta_diff, "KT_Jets_eta_diff/D");
  
  outTree->Branch("jetlenaktkt", &jetlenaktkt, "jetlenaktkt/D");
  
  outTree->Branch("newjetcam_mass", &newjetcam_mass, "newjetcam_mass/D");
  outTree->Branch("newjetcam_mass_2", &newjetcam_mass_2, "newjetcam_mass_2/D");
  outTree->Branch("newjetcam_phi_diff", &newjetcam_phi_diff, "newjetcam_phi_diff/D");
  outTree->Branch("newjetcam_eta_diff", &newjetcam_eta_diff, "newjetcam_eta_diff/D");
  outTree->Branch("newjetcam_phi", &newjetcam_phi, "newjetcam_phi/D");
  outTree->Branch("newjetcam_eta", &newjetcam_eta, "newjetcam_eta/D");
  outTree->Branch("newjetcam_phi_2", &newjetcam_phi_2, "newjetcam_phi_2/D");
  outTree->Branch("newjetcam_eta_2", &newjetcam_eta_2, "newjetcam_eta_2/D");

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
  outTree->Branch("newjetcam_R", &newjetcam_R, "newjetcam_R/D");
  outTree->Branch("finalconesize", &finalconesize, "finalconesize/D");
  
  outTree->Branch("MR_Cont", &MR_Cont, "MR_Cont/D");
  outTree->Branch("MRz_Cont", &MRz_Cont, "MRz_Cont/D");
  outTree->Branch("RSQ_Cont", &RSQ_Cont, "RSQ_Cont/D");
  outTree->Branch("RSQz_Cont", &RSQz_Cont, "RSQz_Cont/D");
  outTree->Branch("Cont_px", &Cont_px, "Cont_px/D");
  outTree->Branch("Cont_px_2", &Cont_px_2, "Cont_px_2/D");
  outTree->Branch("Cont_py", &Cont_py, "Cont_py/D");
  outTree->Branch("Cont_py_2", &Cont_py_2, "Cont_py_2/D");
  outTree->Branch("Cont_pz", &Cont_pz, "Cont_pz/D");
  outTree->Branch("Cont_pz_2", &Cont_pz_2, "Cont_pz_2/D");
  outTree->Branch("Cont_mass", &Cont_mass, "Cont_mass/D");
  outTree->Branch("Cont_mass_2", &Cont_mass_2, "Cont_mass_2/D");
  
  outTree->Branch("newjetkt_px", &newjetkt_px, "newjetkt_px/D");
  outTree->Branch("newjetkt_py", &newjetkt_py, "newjetkt_py/D");
  outTree->Branch("newjetkt_px_2", &newjetkt_px_2, "newjetkt_px_2/D");
  outTree->Branch("newjetkt_py_2", &newjetkt_py_2, "newjetkt_py_2/D");
  
  outTree->Branch("newjetakt_px", &newjetakt_px, "newjetakt_px/D");
  outTree->Branch("newjetakt_py", &newjetakt_py, "newjetakt_py/D");
  outTree->Branch("newjetakt_px_2", &newjetakt_px_2, "newjetakt_px_2/D");
  outTree->Branch("newjetakt_py_2", &newjetakt_py_2, "newjetakt_py_2/D");
  
  
  outTree->Branch("pt_Gen_lead", &pt_Gen_lead, "pt_Gen_lead/D");
  outTree->Branch("pt_Gen_sub", &pt_Gen_sub, "pt_Gen_sub/D");
  outTree->Branch("Gen_lead_phi", &Gen_lead_phi, "Gen_lead_phi/D");
  outTree->Branch("Gen_sub_phi", &Gen_sub_phi, "Gen_sub_phi/D");
  outTree->Branch("Gen_lead_eta", &Gen_lead_eta, "Gen_lead_eta/D");
  outTree->Branch("Gen_sub_eta", &Gen_sub_eta, "Gen_sub_eta/D");

  outTree->Branch("Default_hem1_csts", &Default_hem1_csts, "Default_hem1_csts[100]/I");
  outTree->Branch("Default_hem2_csts", &Default_hem2_csts, "Default_hem2_csts[100]/I");

  outTree->Branch("Georgi_hem1_csts", &Georgi_hem1_csts, "Georgi_hem1_csts[100]/I");
  outTree->Branch("Georgi_hem2_csts", &Georgi_hem2_csts, "Georgi_hem2_csts[100]/I");
  
  outTree->Branch("newjetcam_hem1_csts", &newjetcam_hem1_csts, "newjetcam_hem1_csts[100]/I");
  outTree->Branch("newjetcam_hem2_csts", &newjetcam_hem2_csts, "newjetcam_hem2_csts[100]/I");
  
  outTree->Branch("KT_Jets_hem1_csts", &KT_Jets_hem1_csts, "KT_Jets_hem1_csts[100]/I");
  outTree->Branch("KT_Jets_hem2_csts", &KT_Jets_hem2_csts, "KT_Jets_hem2_csts[100]/I");
  //outTree->Branch("All_MR", &All_MR, "All_MR[20000]/D");
  outTree->Branch("min_MR", &min_MR, "min_MR/D");
  outTree->Branch("max_MR", &max_MR, "max_MR/D");
  //outTree->Branch("All_RSQ", &All_RSQ, "All_RSQ[20000]/D");
  outTree->Branch("min_RSQ", &min_RSQ, "min_RSQ/D");
  outTree->Branch("max_RSQ", &max_RSQ, "max_RSQ/D");
  outTree->Branch("MR_estimator1", &MR_estimator1, "MR_estimator1/D");
  outTree->Branch("MR_estimator2", &MR_estimator2, "MR_estimator2/D");
  outTree->Branch("MR_estimator3", &MR_estimator3, "MR_estimator3/D");
  outTree->Branch("MR_estimator4", &MR_estimator4, "MR_estimator4/D");
  outTree->Branch("MR_perfect_cluster", &MR_perfect_cluster, "MR_perfect_cluster/D");
  outTree->Branch("RSQ_perfect_cluster", &RSQ_perfect_cluster, "RSQ_perfect_cluster/D");
  outTree->Branch("MR_random_cluster", &MR_random_cluster, "MR_random_cluster/D");
  outTree->Branch("parent_mass1", &parent_mass1, "parent_mass1/D");
  outTree->Branch("parent_mass2", &parent_mass2, "parent_mass2/D");
  outTree->Branch("perfect_mass_1", &perfect_mass_1, "perfect_mass_1/D");
  outTree->Branch("perfect_mass_2", &perfect_mass_2, "perfect_mass_2/D");

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
    //for (Long64_t jentry=188; jentry<189;jentry+=1) {
    
    if(verbose) cout << "new event" << endl;

    // clean physics-objects blocks
    CleanEvent();
    
    //cout << "------" << event_counter << "------" << endl;

    // get new event
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // Build the event at generator level
    PFReco();
    vector<fastjet::PseudoJet> empty;
    vector<fastjet::PseudoJet> JetsConst = PFJetConstituents(empty,empty,empty);
    
    //--------------All-hadronic events filter---------------//
    
    bool ttbar_no_all_had = true; //true if you want to filter out all-hadronic events (no leptons)
    bool ttbar = false; //controls gen-level mothers
    bool hasLepton = true;
    if (ttbar_no_all_had){
      hasLepton = false;
      ttbar = true;
    }
    int particleID;
    if (ttbar_no_all_had){
      for (int k = 0; k < JetsConst.size(); k++){
	particleID = JetsConst[k].user_info<ParticleInfo>().pdg_id;
	if (abs(particleID) == 11 || abs(particleID) == 13 || abs(particleID) == 15) hasLepton = true;
      }
    }
    if (!hasLepton) {
      continue;
    }
    //--------------End filter-----------------//
    
    
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
    //cout << "Event " << event_counter << " size " << pfAK04.size() << endl;
    i = 0;

    if(pfAK04.size()<2) continue;
    event_counter = event_counter + 1 ;
    
    //setting arrays initially to -999

    for(unsigned k=0; k<pflen; k++){
      PFJets_phi[k] = -999.;
      PFJets_eta[k] = -999.;
      PFJets_pT[k] = -999.;
      PFJets_m[k] = -999.;
    }
    
    for(unsigned k=0; k<pflen; k++){
      Default_hem1_csts[k] = -999;
      Default_hem2_csts[k] = -999;
      newjetcam_hem1_csts[k] = -999;
      newjetcam_hem2_csts[k] = -999;
      KT_Jets_hem1_csts[k] = -999;
      KT_Jets_hem2_csts[k] = -999;
      Georgi_hem1_csts[k] = -999;
      Georgi_hem2_csts[k] = -999;

    }
    
    //saving angle variables and pT for PFAK04 jets
    for(int k=0; k<pfAK04.size(); k++){
      pfAK04[k].set_user_index(k); //label jet with index
      PFJets_phi[k] = pfAK04[k].phi_std();
      PFJets_eta[k] = pfAK04[k].eta();
      PFJets_pT[k] = pfAK04[k].pt();
      PFJets_m[k] = pfAK04[k].m();
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
    
    double cone_size = 1.57; 
    
    fastjet::ClusterSequence cs1(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, cone_size));
    vector<fastjet::PseudoJet> final_KT_exc_jets = fastjet::sorted_by_pt(cs1.exclusive_jets(2));
    //saving constituents
    for (unsigned k=0; k<final_KT_exc_jets.size(); k++){
      vector<fastjet::PseudoJet> KT_constituents = final_KT_exc_jets[k].constituents();
      int count1 = 0;
      int count2 = 0;
      for (unsigned l=0; l<KT_constituents.size();l++){
	if (k==0) {
	  KT_Jets_hem1_csts[count1] = KT_constituents[l].user_index();
	  count1++;
	}
	else if (k==1) {
	  KT_Jets_hem2_csts[count2] = KT_constituents[l].user_index();
	  count2++;
	}
      }
    }
    
    vector<TLorentzVector> hem_KT_Jets = ConvertTo4Vector(fastjet::sorted_by_pt(cs1.exclusive_jets(2)));
    
    
    //----------------inclusive mode used to find jets thrown out---------------------
    fastjet::ClusterSequence cs101(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, cone_size));
    vector<TLorentzVector> hem_KT_Jets_inc = ConvertTo4Vector(fastjet::sorted_by_pt(cs101.inclusive_jets()));
    
    jetlenaktkt = hem_KT_Jets_inc.size();
    
    
    //---------------------------new clustering algorithms that run through cone size--------------//
    
    //New idea for clustering.  Re-run anti-kt inclusively and increase cone size to lower final jet numer
    //Do the same for kt and cam also
    //Need to fix incase jet length drops to 1
    double jetlen = pfAK04.size(); 
    double conesize1 = 0.4;
    double conestep = 0.05;
    bool changed = false;
    bool slowdown = true;
    vector<TLorentzVector> newjetakt;
    double max_cone_size = 20.0; //turn off by making the max size large 
    
    vector<fastjet::PseudoJet> final_akt_jets;
    while (jetlen > 1 && slowdown && conesize1 < max_cone_size){
      fastjet::ClusterSequence cs11(pfAK04, fastjet::JetDefinition(fastjet::antikt_algorithm, conesize1));
      newjetakt = ConvertTo4Vector(fastjet::sorted_by_pt(cs11.inclusive_jets()));
      jetlen = newjetakt.size();
      
      if (jetlen > 2){ //not yet reached final two so increment cone size
	conesize1 = conesize1 + conestep; 
	changed = true;
      }
      else if (jetlen == 2){ //if reached two then stop
	slowdown = false;
      }
      else if (jetlen == 1){ //if reached one then go back one step
	slowdown = false;
	conesize1 = conesize1 - conestep;
      }
    }
    if (changed) { //cluster with the final cone size 
      fastjet::ClusterSequence cs11(pfAK04, fastjet::JetDefinition(fastjet::antikt_algorithm, conesize1));
      newjetakt = ConvertTo4Vector(fastjet::sorted_by_pt(cs11.inclusive_jets()));
      final_akt_jets = fastjet::sorted_by_pt(cs11.inclusive_jets());
    }
    else { //if we already started with two 
      newjetakt = ConvertTo4Vector(pfAK04);
      //final_akt_jets = pfAK04;
    }
    
    //-------------Cambridge Form------------//
    jetlen = pfAK04.size();
    conestep = 0.05;
    conesize1 = 0.4;
    changed = false;
    slowdown = true;
    vector<TLorentzVector> newjetcam;
    vector<fastjet::PseudoJet> final_jets;
    while (jetlen > 1 && slowdown){
      fastjet::ClusterSequence cs1101(pfAK04, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1));
      newjetcam = ConvertTo4Vector(fastjet::sorted_by_pt(cs1101.inclusive_jets()));
      jetlen = newjetcam.size();

      if (jetlen > 2){ //not yet reached final two so increment cone size
	conesize1 = conesize1 + conestep;
	changed = true;
      }
      else if (jetlen == 2){ //if reached two then stop
	slowdown = false;
      }
      else if (jetlen == 1){ //if reached one step back
	slowdown = false;
	conesize1 = conesize1-conestep;
      }
      else if (max_cone_size < conesize1){slowdown = false;}
    }
    
    if (changed) {
      fastjet::ClusterSequence cs1101(pfAK04, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1));
      newjetcam = ConvertTo4Vector(fastjet::sorted_by_pt(cs1101.inclusive_jets()));
      final_jets = fastjet::sorted_by_pt(cs1101.inclusive_jets()); //save pseudojets
      for (unsigned k=0; k<final_jets.size();k++){
	vector<fastjet::PseudoJet> constituents = final_jets[k].constituents();
	int count1 = 0;
	int count2 = 0;
	for (unsigned l=0; l<constituents.size();l++){
	  if (k==0) { // first index of final_jets is leading hem
	    newjetcam_hem1_csts[count1] = constituents[l].user_index(); //saving jet number
	    //cout << l << " newjetcam1 " << constituents[l].user_index() << endl;
	    count1++;
	  }
	  else if (k==1) { //second index is subleading hem	    
	    newjetcam_hem2_csts[count2] = constituents[l].user_index(); //saving jet number
	    //cout << l << " newjetcam2 " << constituents[l].user_index() << endl;
	    count2++;
	  }
	}
      }
    }
    else {
      newjetcam = ConvertTo4Vector(pfAK04);
      final_jets = pfAK04;
      newjetcam_hem1_csts[0] = 0;
      newjetcam_hem2_csts[0] = 1;
    }
    finalconesize = conesize1;
    
    //kt
    jetlen = pfAK04.size(); 
    conesize1 = 0.4;
    changed = false;
    slowdown = true;
    vector<TLorentzVector> newjetkt;
    while (jetlen > 1 && slowdown && conesize1 < max_cone_size){
      fastjet::ClusterSequence cs12(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, conesize1));
      newjetkt = ConvertTo4Vector(fastjet::sorted_by_pt(cs12.inclusive_jets()));
      jetlen = newjetkt.size();
      if (jetlen > 2){
	conesize1 = conesize1 + conestep;
	changed = true;
      }
      else if (jetlen == 2){
	slowdown = false;
      }
      else if (jetlen == 1){
	slowdown = false;
	conesize1 = conesize1-conestep;
      }
    }
    if (changed){
      fastjet::ClusterSequence cs12(pfAK04, fastjet::JetDefinition(fastjet::kt_algorithm, conesize1));
      newjetkt = ConvertTo4Vector(fastjet::sorted_by_pt(cs12.inclusive_jets()));
    }
    else newjetkt = ConvertTo4Vector(pfAK04);
    
    //End of new idea of clustering
    
    
    //Clusters particles off of previous clustering.  Should produce the same results as CAM with varying conesize
    
    jetlen = pfAK04.size();
    conesize1 = 0.2;
    changed = false;
    vector<fastjet::PseudoJet> Continualjet = pfAK04;
    while (jetlen > 2){
      fastjet::ClusterSequence cs1000(Continualjet, fastjet::JetDefinition(fastjet::cambridge_algorithm, conesize1));
      Continualjet = fastjet::sorted_by_pt(cs1000.inclusive_jets()); 
      jetlen = Continualjet.size();
      conesize1 += conestep;
    }
    vector<TLorentzVector> Continual;
    Continual = ConvertTo4Vector(Continualjet);
    

    //MR array
    /*
    for (int k = 0; k < MR_len; k++){
      All_MR[k] = -9999.;
      All_RSQ[k] = -1.;
    }
    */
    //Used to find min/max values of MR and RSQ
    CMSHemisphere* myHemMR = new CMSHemisphere(ConvertTo4Vector(pfAK04));
    myHemMR->Find_All_MR(); 
    
    vector<TLorentzVector> All_MR_raw = myHemMR->GetHemispheres();
    min_MR = 5000; max_MR = 0; min_RSQ = 2; max_RSQ = 0;
    //permutations are divided by 2 because there are 2 jets and another 2
    //because of redundancy in finding all possible combinations
    //
    //sometimes the varying R algorithms return a more extreme value of 
    //MR or RSQ than thought possible because they will cluster to 3 jets 
    //and only use the highest 2 pt ones
    int permutations = All_MR_raw.size()/4;
    for (int k = 0; k < permutations; k++){
      double temp_MR = CalcMR(All_MR_raw[2*k], All_MR_raw[2*k + 1]);
      double temp_RSQ = pow(CalcMRT(All_MR_raw[2*k], All_MR_raw[2*k+1], PFMET),2.)/temp_MR/temp_MR;
      //All_MR[k] = temp_MR; 
      //All_RSQ[k] = temp_RSQ;
      if (temp_MR > max_MR){max_MR = temp_MR;}
      if (temp_MR < min_MR){min_MR = temp_MR;}
      if (temp_RSQ > max_RSQ){max_RSQ = temp_RSQ;}
      if (temp_RSQ < min_RSQ){min_RSQ = temp_RSQ;}
    
     
    }    
    delete myHemMR;
  
    // 1b) traditional hemispheres
    CMSHemisphere* myHem = new CMSHemisphere(ConvertTo4Vector(pfAK04));
    myHem->CombineMinMass(); //CombineMinMass() to use megajet and CombineGeorgi() to use Georgi
    vector<TLorentzVector> hem_Default = myHem->GetHemispheres();
    vector<int> Temporary = myHem->GetHem1Constituents();
    vector<int> Temporary2 = myHem->GetHem2Constituents();
    int count1 = 0;
    int count2 = 0;
    for (int k=0; k < Temporary.size(); k++){
      if (Temporary[k]>-1) { 
	Default_hem1_csts[count1] = Temporary[k];
	//cout << "Default hem: " <<Default_hem1_csts[count1] << endl;
	count1++;
      }
    }
    for (int k=0; k < Temporary2.size(); k++){
      if (Temporary2[k]>-1) {
	Default_hem2_csts[count2] = Temporary2[k];
	//cout << "Default hem2: " << Default_hem2_csts[count2] <<endl;
	count2++;
      }
    }
    delete myHem; 


    CMSHemisphere* myHemGeorgi = new CMSHemisphere(ConvertTo4Vector(pfAK04));
    myHem->CombineGeorgi(); //CombineMinMass() to use megajet and CombineGeorgi() to use Georgi
    vector<TLorentzVector> hem_Georgi = myHem->GetHemispheres();
    vector<int> Temporary3 = myHemGeorgi->GetHem1Constituents();
    vector<int> Temporary4 = myHemGeorgi->GetHem2Constituents();
    int count3 = 0;
    int count4 = 0;
    for (int k=0; k < Temporary3.size(); k++){
      if (Temporary3[k]>-1) { 
	Georgi_hem1_csts[count1] = Temporary3[k];
	//cout << "Default hem: " <<Default_hem1_csts[count1] << endl;
	count1++;
      }
    }
    for (int k=0; k < Temporary4.size(); k++){
      if (Temporary4[k]>-1) {
	Georgi_hem2_csts[count2] = Temporary4[k];
	//cout << "Default hem2: " << Default_hem2_csts[count2] <<endl;
	count2++;
      }
    }
    delete myHemGeorgi; 

    //initialize  
    MR_KT_Jets = -9999. ;
    RSQ_KT_Jets = -9999. ;
    MR_Default = -9999. ;
    RSQ_Default = -9999. ;
    MRz_KT_Jets = -9999. ;
    MRz_Default = -9999. ;
    RSQz_KT_Jets = -9999. ;
    RSQz_Default = -9999. ;
    MRt_Default = -9999.;
    MRt_KT_Jets = -9999.;
    // Edward's new algorithm
    
    MR_newjetakt = -9999. ;
    MRz_newjetakt = -9999. ;
    RSQ_newjetakt = -9999. ;
    RSQz_newjetakt = -9999. ;
    MRt_newjetakt = -9999. ;
    MR_newjetkt = -9999.;
    MRz_newjetkt = -9999. ;
    RSQ_newjetkt = -9999. ;
    RSQz_newjetkt = -9999. ;
    MRt_newjetkt = -9999.;
    MR_newjetcam = -9999.;
    MRz_newjetcam = -9999. ;
    RSQ_newjetcam = -9999. ;
    RSQz_newjetcam = -9999. ;    
    MRt_newjetakt = -9999.;
    
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
    KT_Jets_eta = -9999.0;
    
    Default_mass = -9999.0;
    Default_px = -9999.0;
    Default_py = -9999.0;
    Default_pz = -9999.0;
    Default_phi = -9999.0;
    Default_eta = -9999.0;
    
    KT_Jets_mass_2 = -9999.0;
    KT_Jets_px_2 = -9999.0;
    KT_Jets_py_2 = -9999.0;
    KT_Jets_pz_2 = -9999.0;
    KT_Jets_phi_2 = -9999.0;
    KT_Jets_eta_2 = -9999.0;
    
    Default_mass_2 = -9999.0;
    Default_px_2 = -9999.0;
    Default_py_2 = -9999.0;
    Default_pz_2 = -9999.0;
    Default_phi_2 = -9999.0;
    Default_eta_2 = -9999.0;
    
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
      KT_Jets_eta = j1.Eta();
      KT_Jets_mass_2 = j2.M();
      KT_Jets_px_2 = j2.Px();
      KT_Jets_py_2 = j2.Py();
      KT_Jets_pz_2 = j2.Pz();
      KT_Jets_phi_2 = j2.Phi();
      KT_Jets_eta_2 = j2.Eta();
      KT_Jets_eta_diff = (KT_Jets_eta - KT_Jets_eta_2);
      KT_Jets_phi_diff = (KT_Jets_phi - KT_Jets_phi_2);
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
      newjetcam_mass = j1.M();
      newjetcam_mass_2 = j2.M();
      newjetcam_phi = j1.Phi();
      newjetcam_eta = j1.Eta();
      newjetcam_phi_2 = j2.Phi();
      newjetcam_eta_2 = j2.Eta();
      newjetcam_phi_diff = fabs(j1.DeltaPhi(j2));
      newjetcam_eta_diff = fabs(newjetcam_eta - newjetcam_eta_2);
      newjetcam_R = D_between_vectors(j1, j2);
      
      MR_newjetcam = CalcMR(j1, j2);      
      MRt_newjetcam = CalcMRT(j1, j2, PFMET);
      RSQ_newjetcam = pow(CalcMRT(j1, j2, PFMET),2.)/MR_newjetcam/MR_newjetcam;
      MRz_newjetcam = CalcMR_zinvariant(j1, j2);
      RSQz_newjetcam = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_newjetcam/MRz_newjetcam;
    } 
    else cout<<"FAILED"<<endl;
    
    if (Continual.size() > 1){
      j1 = Continual[0];
      j2 = Continual[1];
      Cont_px = j1.Px();
      Cont_py = j1.Py();
      Cont_pz = j1.Pz();
      Cont_mass = j1.M();
      Cont_phi_diff = j1.Phi() - j2.Phi();
      Cont_eta_diff = j1.Eta() - j2.Eta();
      Cont_px_2 = j2.Px();
      Cont_py_2 = j2.Py();
      Cont_pz_2 = j2.Pz();
      Cont_mass_2 = j2.M();
      MR_Cont = CalcMR(j1, j2);
      MRz_Cont = CalcMR_zinvariant(j1, j2);
      RSQ_Cont = pow(CalcMRT(j1, j2, PFMET),2.)/MR_Cont/MR_Cont;
      RSQz_Cont = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_Cont/MRz_Cont;
      
      
    }
    
    // 2b) compute traditional RSQ and MR (DEFAULT)
    j1 = hem_Default[0];
    j2 = hem_Default[1];
    Default_E = j1.E();
    Default_mass = j1.M();
    Default_px = j1.Px();
    Default_py = j1.Py();
    Default_pz = j1.Pz();
    Default_phi = j1.Phi();
    Default_eta = j1.Eta();
    Default_mass_2 = j2.M();
    Default_px_2 = j2.Px();
    Default_py_2 = j2.Py();
    Default_pz_2 = j2.Pz();
    Default_phi_2 = j2.Phi();
    Default_eta_2 = j2.Eta();
    Default_eta_diff = Default_eta - Default_eta_2;
    Default_phi_diff = Default_phi - Default_phi_2;
    Default_R = pow(pow(Default_eta_diff, 2) + pow(Default_phi_diff, 2), 0.5);    
    MR_Default = CalcMR(j1, j2);
    MRt_Default = CalcMRT(j1, j2, PFMET);
    RSQ_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MR_Default/MR_Default;
    MRz_Default = CalcMR_zinvariant(j1, j2);
    RSQz_Default = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_Default/MRz_Default;
    
    //delete j1;
    //delete j2;

    j1 = hem_Georgi[0];
    j2 = hem_Georgi[1];
    Georgi_mass = j1.M();
    Georgi_px = j1.Px();
    Georgi_py = j1.Py();
    Georgi_pz = j1.Pz();
    Georgi_phi = j1.Phi();
    Georgi_eta = j1.Eta();
    Georgi_mass_2 = j2.M();
    Georgi_px_2 = j2.Px();
    Georgi_py_2 = j2.Py();
    Georgi_pz_2 = j2.Pz();
    Georgi_phi_2 = j2.Phi();
    Georgi_eta_2 = j2.Eta();
    Georgi_eta_diff = Georgi_eta - Georgi_eta_2;
    Georgi_phi_diff = Georgi_phi - Georgi_phi_2;
    Georgi_R = pow(pow(Georgi_eta_diff, 2) + pow(Georgi_phi_diff, 2), 0.5);
    
    MR_Georgi = CalcMR(j1, j2);
    MRt_Georgi = CalcMRT(j1, j2, PFMET);
    RSQ_Georgi = pow(CalcMRT(j1, j2, PFMET),2.)/MR_Georgi/MR_Georgi;
    MRz_Georgi = CalcMR_zinvariant(j1, j2);
    RSQz_Georgi = pow(CalcMRT(j1, j2, PFMET),2.)/MRz_Georgi/MRz_Georgi;
    
    
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
    //outTree->Fill();
    
    // fill PDF histograms
    bool fillBox = SignalRegion(MR_Default, RSQ_Default, BOX_NUM);
    if(BOX_NUM == 0 && fillBox) pdfMuEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 1 && fillBox) pdfMuMu->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 2 && fillBox) pdfEleEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 3 && fillBox) pdfMu->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 4 && fillBox) pdfEle->Fill(MR_Default, RSQ_Default);
    if(BOX_NUM == 5 && fillBox) pdfHad->Fill(MR_Default, RSQ_Default);
    
    
      //continue;  
    
    
    //BEGIN GEN LEVEL ANALYSIS
    //initialize
    
    int correct_clustering = -9999.0 ;
    for(unsigned k=0; k<pflen; k++){
      gen_hem1[k] = -999;
      gen_hem2[k] = -999;
    }
    //initialize structures to store gen data pulled from cmsreco
    //--------Finds parent SUSY of each LSP------//
    int SUSY;
    vector<double> SUSYPx;
    vector<double> SUSYPy;
    vector<double> SUSYPz;
    vector<double> SUSYE;
    vector<int> SUSYPdgId;
    vector<int> SUSYM1PdgId;
    vector<int> SUSYStatus;
    vector<double> SUSYm1px;
    //Vector with index of LSP location
    vector<int> SUSYLSP;
    //Vector with index of first mother location
    vector<int> SUSYLSP_Mother_Eve;
    //Function used to return information about SUSY particles
    SUSYReturn(SUSY, SUSYPx, SUSYPy, SUSYPz, SUSYE, SUSYPdgId, SUSYM1PdgId, SUSYStatus, SUSYm1px);
    //Looks for the first mother of an LSP and finds the difference in their pt
    if (SUSY > 0){
      for (int p = 0; p < SUSY; p++){
	//Status = 1 is final state (LSP) and status = 22 from the initial collision 
	//Looking to find the LSP's parent SUSY particle
	if (SUSYStatus[p] == 1){
	  SUSYLSP.push_back(p);
	  //Continue running until the parent has been found.
	  int mom_status = 1; //this line is an arbitrary initializing value 
	  int id = SUSYPdgId[p];
	  int mom_id = SUSYM1PdgId[p];
	  double mom_px = SUSYm1px[p];
	  double final_momentum = 0;
	  int finali = 0;
	  while (mom_status != 22){
	    //cout << "Number " << finali << " " << mom_px << "  " << mom_status << endl;
	    for (int i = 0; i < SUSY; i++){
	      //Finds the parent by matching px to the parent
	      if (mom_id ==  SUSYPdgId[i] && mom_px == SUSYPx[i]){
		id = SUSYPdgId[i];
		mom_id = SUSYM1PdgId[i];
		mom_px = SUSYm1px[i];
		final_momentum = SUSYPx[i];
		mom_status = SUSYStatus[i];
		finali = i;
		i = SUSY;
	      }
	    }
	  }
	  SUSYLSP_Mother_Eve.push_back(finali);
	}
      }
      double pt_diff1, pt_diff2;
      TLorentzVector Mother1 (SUSYPx[SUSYLSP_Mother_Eve[0]], SUSYPy[SUSYLSP_Mother_Eve[0]], SUSYPz[SUSYLSP_Mother_Eve[0]], SUSYE[SUSYLSP_Mother_Eve[0]]);
      TLorentzVector Mother2 (SUSYPx[SUSYLSP_Mother_Eve[1]], SUSYPy[SUSYLSP_Mother_Eve[1]], SUSYPz[SUSYLSP_Mother_Eve[1]], SUSYE[SUSYLSP_Mother_Eve[1]]);
      TLorentzVector Child1 (SUSYPx[SUSYLSP[0]], SUSYPy[SUSYLSP[0]], SUSYPz[SUSYLSP[0]], SUSYE[SUSYLSP[0]]);
      TLorentzVector Child2 (SUSYPx[SUSYLSP[1]], SUSYPy[SUSYLSP[1]], SUSYPz[SUSYLSP[1]], SUSYE[SUSYLSP[1]]);
      //V1 is the vector of the sum of the other decay products.  
      TLorentzVector V1 = Mother1 - Child1;
      TLorentzVector V2 = Mother2 - Child2;
      pt_diff1 =  pow(pow(V1.Px(), 2) + pow(V1.Py(), 2), 0.5);
      pt_diff2 =  pow(pow(V2.Px(), 2) + pow(V2.Py(), 2), 0.5);
      //These gen-level values are currently written over by using the hardest quarks.
       if (pt_diff1 > pt_diff2){
	pt_Gen_lead = pt_diff1;
	pt_Gen_sub = pt_diff2;
	Gen_lead_phi = V1.Phi();
	Gen_lead_eta = V1.Eta();
	Gen_sub_phi = V2.Phi();
	Gen_sub_eta = V2.Eta();
      }
      else{
	pt_Gen_lead = pt_diff2;
	pt_Gen_sub = pt_diff1;
	Gen_lead_phi = V2.Phi();
	Gen_lead_eta = V2.Eta();
	Gen_sub_phi = V1.Phi();
	Gen_sub_eta = V1.Eta();
      }       
    }

    int genParticle;
    int numQuarks_1=0;
    int numQuarks_2=0;
    int threshold = 4; //how many quarks we should see: set at 4 for T1, 2 for ttbar
    double mother_px = -999.;
    vector<double> genParticlePx;
    vector<double> genParticlePy;
    vector<double> genParticlePz;
    vector<double> genParticleE;
    vector<int> genParticlePdgId;
    vector<int> genParticleM1PdgId;
    vector<int> genParticleStatus;
    vector<double> genParticlem1px;
    vector<int> genParticleQuark; //Vector with index of quark from gluino location
    vector<TLorentzVector> genParticleQuarkVector_hem1;
    vector<TLorentzVector> genParticleQuarkVector_hem2;

    // to determine from which pair-produced particle PFAK04 jet came from
    GenReturn(genParticle, genParticlePx, genParticlePy, genParticlePz, genParticleE, genParticlePdgId, genParticleM1PdgId, genParticleStatus, genParticlem1px);
    if (genParticle > 0){
      for (int p = 0; p < genParticle; p++){ //loop through genParticles
	if (!ttbar){
	  if (abs(genParticlePdgId[p]) > 0 && abs(genParticlePdgId[p]) < 9){ //if it is a quark
	    if (abs(genParticleM1PdgId[p]) == 1000021 || abs((genParticleM1PdgId[p]) > 1000000 && abs(genParticleM1PdgId[p]) < 1000007) || (abs(genParticleM1PdgId[p]) > 2000000 && abs(genParticleM1PdgId[p]) < 2000007)) {// if mother is a gluino or a squark
	      
	      genParticleQuark.push_back(p); // push back index of quark
	      if (mother_px < -990.) { // set px of first gluino
		mother_px = genParticlem1px[p];
	      }
	      TLorentzVector quark;
	      quark.SetPxPyPzE(genParticlePx[p],genParticlePy[p],genParticlePz[p],genParticleE[p]);
	      if (genParticlem1px[p] == mother_px){
		genParticleQuarkVector_hem1.push_back(quark); //push back 4-vector of quark into first hemisphere
		numQuarks_1++;
	      }
	      else {
		genParticleQuarkVector_hem2.push_back(quark); //push back 4-vector of quark into second hemisphere
		numQuarks_2++; 
	      }
	    }
	  }
	}
	  
	//ttbar search
	else {
	  if (fabs(genParticlePdgId[p]) > 0 && fabs(genParticlePdgId[p]) < 9){ //if it is a quark                     
            if (fabs(genParticleM1PdgId[p]) == 6) {// if mother is a top quark
              genParticleQuark.push_back(p); // push back index of quark
              TLorentzVector quark;
              quark.SetPxPyPzE(genParticlePx[p],genParticlePy[p],genParticlePz[p],genParticleE[p]);
              if (genParticleM1PdgId[p] == 6){ // if top
                genParticleQuarkVector_hem1.push_back(quark); //push back 4-vector of quark into first hemisphere            
                numQuarks_1++;
              }
              else if (genParticleM1PdgId[p] == -6){ // if antitop
                genParticleQuarkVector_hem2.push_back(quark); //push back 4-vector of quark into second hemisphere           
                numQuarks_2++;
	      }
	    }
	    if (fabs(genParticleM1PdgId[p]) == 24) { // if mother is a W+/W- (accounts for hadronic decays)
	      genParticleQuark.push_back(p); 
	      TLorentzVector quark;
	      quark.SetPxPyPzE(genParticlePx[p], genParticlePy[p], genParticlePx[p], genParticleE[p]);
	      if (genParticleM1PdgId[p] == 24){ // if mother is W+
		genParticleQuarkVector_hem1.push_back(quark);
		numQuarks_1++;
	      }
	      else if (genParticleM1PdgId[p] == -24){ //if mother is W-
		genParticleQuarkVector_hem2.push_back(quark);
		numQuarks_2++;
	      }
	    }
	  }
	}
      }
    }

    // too many quarks 
    if (ttbar && ((numQuarks_1+numQuarks_2) >= 6)) {
      continue;
    }
    
    if ((numQuarks_1+numQuarks_2) >= threshold){ // passed threshold for num of quarks
      int hem1_count = 0;
      int hem2_count = 0;
      for (int q=0; q<numQuarks_1; q++){ //loops through quarks from first pair-produced particle
	for (int looper = 0; looper < pfAK04.size(); looper++){ //loops through all PFJets
	  double eta_diff = PFJets_eta[looper] - genParticleQuarkVector_hem1[q].Eta(); 
	  double phi_diff = PFJets_phi[looper] - genParticleQuarkVector_hem1[q].Phi();
	  double dR = sqrt(eta_diff*eta_diff + phi_diff*phi_diff); //calculate deltaR between original quark and the PFAK04 Jet
	  if (dR < 0.5) {
	    gen_hem1[hem1_count] = looper; //pushes back index of PFJet that matches deltaR
	    hem1_count++; //loops through gen_hem1
	  }
	}
      }
      for (int q=0; q<numQuarks_2; q++){
	for (int looper=0; looper < pfAK04.size(); looper++){ //loops through all PFJets
	  double eta_diff = PFJets_eta[looper] - genParticleQuarkVector_hem2[q].Eta();
	  double phi_diff = PFJets_phi[looper] - genParticleQuarkVector_hem2[q].Phi();
	  double dR = sqrt(eta_diff*eta_diff + phi_diff*phi_diff);
	  if (dR < 0.5) {
	    gen_hem2[hem2_count] = looper; //pushes back index of PFJet
	    hem2_count++; //loops through gen_hem2
	  }
	}
      }
    }
    
    //Used to find the quarks directly.  Using that as gen data now instead of
    //gluino - lsp.  Only confirmed to work for T1, T2
    TLorentzVector parent1, parent2, lsp1, lsp2, quark1, quark2, quark3, quark4, lsp;
    vector<int> parent_locations, lsp_locations, quark_locations;
    //If signal, go through this process to find gen-level info.  
    if (SUSY > 0){
    GenReturn(genParticle, genParticlePx, genParticlePy, genParticlePz, genParticleE, genParticlePdgId, genParticleM1PdgId, genParticleStatus, genParticlem1px);
    if (genParticle > 0){
      for (int p = 0; p < genParticle; p++){ //loop through genParticles
	if ((abs(genParticlePdgId[p]) > 1000000 && abs(genParticlePdgId[p]) < 1000022) || (abs(genParticlePdgId[p]) > 2000000 && abs(genParticlePdgId[p]) < 2000006)){
	  parent_locations.push_back(p);
	}
	if (abs(genParticlePdgId[p]) == 1000022){
	  lsp_locations.push_back(p);
	}
	//px requirement to ensure the incoming beam quarks are not used
	if ((abs(genParticlePdgId[p]) > 0 && abs(genParticlePdgId[p]) < 9) &&  fabs(genParticlePx[p]) > 0.1){
	    quark_locations.push_back(p);
	  }
	  }}
    parent1.SetPxPyPzE(genParticlePx[parent_locations[0]], genParticlePy[parent_locations[0]], genParticlePz[parent_locations[0]], genParticleE[parent_locations[0]]);
    parent2.SetPxPyPzE(genParticlePx[parent_locations[1]], genParticlePy[parent_locations[1]], genParticlePz[parent_locations[1]], genParticleE[parent_locations[1]]);    
    quark1.SetPxPyPzE(genParticlePx[quark_locations[0]], genParticlePy[quark_locations[0]], genParticlePz[quark_locations[0]], genParticleE[quark_locations[0]]);
    quark2.SetPxPyPzE(genParticlePx[quark_locations[1]], genParticlePy[quark_locations[1]], genParticlePz[quark_locations[1]], genParticleE[quark_locations[1]]);    
    if (lsp_locations.size() == 2){
      lsp1.SetPxPyPzE(genParticlePx[lsp_locations[0]], genParticlePy[lsp_locations[0]], genParticlePz[lsp_locations[0]], genParticleE[lsp_locations[0]]);
      lsp2.SetPxPyPzE(genParticlePx[lsp_locations[1]], genParticlePy[lsp_locations[1]], genParticlePz[lsp_locations[1]], genParticleE[lsp_locations[1]]);    
    }
    
    double temp_gen_pt_1 = pow(pow(quark1.Px(), 2) + pow(quark1.Py(), 2), 0.5);
    double temp_gen_pt_2 = pow(pow(quark2.Px(), 2) + pow(quark2.Py(), 2), 0.5);
    if (temp_gen_pt_1 > temp_gen_pt_2){
      pt_Gen_lead = temp_gen_pt_1;
      pt_Gen_sub = temp_gen_pt_2;
      Gen_lead_phi = quark1.Phi();
      Gen_lead_eta = quark1.Eta();
      Gen_sub_phi = quark2.Phi();
      Gen_sub_eta = quark2.Eta();

    }
    else {
      pt_Gen_lead = temp_gen_pt_2; 
      pt_Gen_sub = temp_gen_pt_1;
      Gen_lead_phi = quark2.Phi();
      Gen_lead_eta = quark2.Eta();
      Gen_sub_phi = quark1.Phi();
      Gen_sub_eta = quark1.Eta();
    }
    //begin MR estimation
    TVector3 beta_parent_1, beta_parent_2;
    //T1.  More ordered so there is less filtering necessary
    TLorentzVector quark1_new, quark2_new;
    bool T2 = true;
    //if more than 2 quarks are found, assume T1
    if (quark_locations.size() > 2){ 
      T2 = false;
      quark3.SetPxPyPzE(genParticlePx[quark_locations[2]], genParticlePy[quark_locations[2]], genParticlePz[quark_locations[2]], genParticleE[quark_locations[2]]);
      quark4.SetPxPyPzE(genParticlePx[quark_locations[3]], genParticlePy[quark_locations[3]], genParticlePz[quark_locations[3]], genParticleE[quark_locations[3]]);  
      //the first and second quarks always belong to the first parent for T1
      quark1_new = quark1 + quark2;
      quark2_new = quark3 + quark4;
      MR_perfect_cluster = CalcMR(quark1_new, quark2_new);
      RSQ_perfect_cluster = pow(CalcMRT(quark1_new, quark2_new, PFMET),2.)/MR_perfect_cluster/MR_perfect_cluster;
      MR_random_cluster = CalcMR(quark1 + quark3, quark2 + quark4);
      parent_mass1 = parent1.M();
      parent_mass2 = parent2.M();
      perfect_mass_1 = quark1_new.M();
      perfect_mass_2 = quark2_new.M();
      //boost into each gluino's rest frame.  
      //Verified that each gluino's 4 vector is (m, 0) after boost
      beta_parent_1 = TVector3 (parent1.Px(), parent1.Py(), parent1.Pz()) * (1/parent1.E());
      beta_parent_2 = TVector3 (parent2.Px(), parent2.Py(), parent2.Pz()) * (1/parent2.E());
      quark1_new.Boost(-beta_parent_1);
      lsp1.Boost(-beta_parent_1);
      quark2_new.Boost(-beta_parent_2);
      lsp2.Boost(-beta_parent_2);
      //MR_estimator + the mass is a constant, which is m#Delta
      //MR_estimator1 = quark1_new.P() + quark1_new.M2()/(parent1.M() * 2);
      //MR_estimator2 = quark2_new.P() + quark2_new.M2()/(parent2.M() * 2);
      if (lsp_locations.size() == 2){
	//These MR estimators are not accurate.  What should MR estimate for T1?
	MR_estimator1 = quark1_new.P();
	MR_estimator2 = lsp1.P();
	MR_estimator3 = quark2_new.P();
	MR_estimator4 = lsp2.P();
      }
      else{
	MR_estimator1 = -999;
	MR_estimator2 = -999;
	MR_estimator3 = -999;
	MR_estimator4 = -999;
	  }
    }
    //end T1
    if (T2){
      // if quark and squark are of the same sign
      //boost into CM frame
      double beta_z= (parent1.Pz() + parent2.Pz())/(parent1.E() + parent2.E());
      TVector3 beta_boostz = TVector3 (0, 0, beta_z);	 
      MR_perfect_cluster = CalcMR(quark1, quark2);
      RSQ_perfect_cluster = pow(CalcMRT(quark1, quark2, PFMET),2.)/MR_perfect_cluster/MR_perfect_cluster;
      parent_mass1 = parent1.M();
      parent_mass2 = parent2.M();
      perfect_mass_1 = quark1.M();
      perfect_mass_2 = quark2.M();
      parent1.Boost(-beta_boostz);
      parent2.Boost(-beta_boostz);
      quark1.Boost(-beta_boostz);
      quark2.Boost(-beta_boostz);
      //if the quark matches the squark, boost into each parent's rest frame
      if (genParticleM1PdgId[quark_locations[0]] == genParticlePdgId[parent_locations[0]]){
	beta_parent_1 = TVector3 (parent1.Px(), parent1.Py(), parent1.Pz()) * (1/parent1.E());
	beta_parent_2 = TVector3 (parent2.Px(), parent2.Py(), parent2.Pz()) * (1/parent2.E());
	parent_mass1 = parent1.M();
	parent_mass2 = parent2.M();
      }
      else{
	beta_parent_1 = TVector3 (parent2.Px(), parent2.Py(), parent2.Pz()) * (1/parent2.E());
	beta_parent_2 = TVector3 (parent1.Px(), parent1.Py(), parent1.Pz()) * (1/parent1.E());
	parent_mass1 = parent2.M();
	parent_mass2 = parent1.M();
	
      }
      //boost into parent's rest frame
      quark1.Boost(-beta_parent_1);
      quark2.Boost(-beta_parent_2);
      MR_estimator1 = quark1.P();
      MR_estimator2 = quark2.P();
    }}

    outTree->Fill();

    continue;
    
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
      cout <<"parent " << parent[p];
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


