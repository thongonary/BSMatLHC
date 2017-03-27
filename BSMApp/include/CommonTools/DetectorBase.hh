//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 27 11:48:04 2017 by ROOT version 6.08/00
// from TTree GenEvent/GenEvent
// found on file: test_GenTree.root
//////////////////////////////////////////////////////////

#ifndef DetectorBase_h
#define DetectorBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DetectorBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Muon;
   Float_t         MuonE[5];   //[Muon]
   Float_t         MuonPx[5];   //[Muon]
   Float_t         MuonPy[5];   //[Muon]
   Float_t         MuonPz[5];   //[Muon]
   Float_t         MuonX[5];   //[Muon]
   Float_t         MuonY[5];   //[Muon]
   Float_t         MuonMass[5];   //[Muon]
   Int_t           MuonPdgId[5];   //[Muon]
   Int_t           MuonStatus[5];   //[Muon]
   Int_t           MuonDecayLmm[5];   //[Muon]
   Int_t           MuonM1PdgId[5];   //[Muon]
   Int_t           MuonM2PdgId[5];   //[Muon]
   Int_t           MuonD1PdgId[5];   //[Muon]
   Float_t         Muonm1px[5];   //[Muon]
   Int_t           Electron;
   Float_t         ElectronE[17];   //[Electron]
   Float_t         ElectronPx[17];   //[Electron]
   Float_t         ElectronPy[17];   //[Electron]
   Float_t         ElectronPz[17];   //[Electron]
   Float_t         ElectronX[17];   //[Electron]
   Float_t         ElectronY[17];   //[Electron]
   Float_t         ElectronMass[17];   //[Electron]
   Int_t           ElectronPdgId[17];   //[Electron]
   Int_t           ElectronStatus[17];   //[Electron]
   Int_t           ElectronDecayLmm[17];   //[Electron]
   Int_t           ElectronM1PdgId[17];   //[Electron]
   Int_t           ElectronM2PdgId[17];   //[Electron]
   Int_t           ElectronD1PdgId[17];   //[Electron]
   Float_t         Electronm1px[17];   //[Electron]
   Int_t           Tau;
   Float_t         TauE[6];   //[Tau]
   Float_t         TauPx[6];   //[Tau]
   Float_t         TauPy[6];   //[Tau]
   Float_t         TauPz[6];   //[Tau]
   Float_t         TauX[6];   //[Tau]
   Float_t         TauY[6];   //[Tau]
   Float_t         TauMass[6];   //[Tau]
   Int_t           TauPdgId[6];   //[Tau]
   Int_t           TauStatus[6];   //[Tau]
   Int_t           TauDecayLmm[6];   //[Tau]
   Int_t           TauM1PdgId[6];   //[Tau]
   Int_t           TauM2PdgId[6];   //[Tau]
   Int_t           TauD1PdgId[6];   //[Tau]
   Float_t         Taum1px[6];   //[Tau]
   Int_t           b;
   Float_t         bE[57];   //[b]
   Float_t         bPx[57];   //[b]
   Float_t         bPy[57];   //[b]
   Float_t         bPz[57];   //[b]
   Float_t         bX[57];   //[b]
   Float_t         bY[57];   //[b]
   Float_t         bMass[57];   //[b]
   Int_t           bPdgId[57];   //[b]
   Int_t           bStatus[57];   //[b]
   Int_t           bDecayLmm[57];   //[b]
   Int_t           bM1PdgId[57];   //[b]
   Int_t           bM2PdgId[57];   //[b]
   Int_t           bD1PdgId[57];   //[b]
   Float_t         bm1px[57];   //[b]
   Int_t           c;
   Float_t         cE[79];   //[c]
   Float_t         cPx[79];   //[c]
   Float_t         cPy[79];   //[c]
   Float_t         cPz[79];   //[c]
   Float_t         cX[79];   //[c]
   Float_t         cY[79];   //[c]
   Float_t         cMass[79];   //[c]
   Int_t           cPdgId[79];   //[c]
   Int_t           cStatus[79];   //[c]
   Int_t           cDecayLmm[79];   //[c]
   Int_t           cM1PdgId[79];   //[c]
   Int_t           cM2PdgId[79];   //[c]
   Int_t           cD1PdgId[79];   //[c]
   Float_t         cm1px[79];   //[c]
   Int_t           Photon;
   Float_t         PhotonE[501];   //[Photon]
   Float_t         PhotonPx[501];   //[Photon]
   Float_t         PhotonPy[501];   //[Photon]
   Float_t         PhotonPz[501];   //[Photon]
   Float_t         PhotonX[501];   //[Photon]
   Float_t         PhotonY[501];   //[Photon]
   Float_t         PhotonMass[501];   //[Photon]
   Int_t           PhotonPdgId[501];   //[Photon]
   Int_t           PhotonStatus[501];   //[Photon]
   Int_t           PhotonDecayLmm[501];   //[Photon]
   Int_t           PhotonM1PdgId[501];   //[Photon]
   Int_t           PhotonM2PdgId[501];   //[Photon]
   Int_t           PhotonD1PdgId[501];   //[Photon]
   Float_t         Photonm1px[501];   //[Photon]
   Int_t           Neutrino;
   Float_t         NeutrinoE[9];   //[Neutrino]
   Float_t         NeutrinoPx[9];   //[Neutrino]
   Float_t         NeutrinoPy[9];   //[Neutrino]
   Float_t         NeutrinoPz[9];   //[Neutrino]
   Float_t         NeutrinoX[9];   //[Neutrino]
   Float_t         NeutrinoY[9];   //[Neutrino]
   Float_t         NeutrinoMass[9];   //[Neutrino]
   Int_t           NeutrinoPdgId[9];   //[Neutrino]
   Int_t           NeutrinoStatus[9];   //[Neutrino]
   Int_t           NeutrinoDecayLmm[9];   //[Neutrino]
   Int_t           NeutrinoM1PdgId[9];   //[Neutrino]
   Int_t           NeutrinoM2PdgId[9];   //[Neutrino]
   Int_t           NeutrinoD1PdgId[9];   //[Neutrino]
   Float_t         Neutrinom1px[9];   //[Neutrino]
   Int_t           SUSY;
   Float_t         SUSYE[26];   //[SUSY]
   Float_t         SUSYPx[26];   //[SUSY]
   Float_t         SUSYPy[26];   //[SUSY]
   Float_t         SUSYPz[26];   //[SUSY]
   Float_t         SUSYX[26];   //[SUSY]
   Float_t         SUSYY[26];   //[SUSY]
   Float_t         SUSYMass[26];   //[SUSY]
   Int_t           SUSYPdgId[26];   //[SUSY]
   Int_t           SUSYStatus[26];   //[SUSY]
   Int_t           SUSYDecayLmm[26];   //[SUSY]
   Int_t           SUSYM1PdgId[26];   //[SUSY]
   Int_t           SUSYM2PdgId[26];   //[SUSY]
   Int_t           SUSYD1PdgId[26];   //[SUSY]
   Float_t         SUSYm1px[26];   //[SUSY]
   Int_t           GenTreeParticle;
   Float_t         GenTreeParticleE[4];   //[GenTreeParticle]
   Float_t         GenTreeParticlePx[4];   //[GenTreeParticle]
   Float_t         GenTreeParticlePy[4];   //[GenTreeParticle]
   Float_t         GenTreeParticlePz[4];   //[GenTreeParticle]
   Float_t         GenTreeParticleX[4];   //[GenTreeParticle]
   Float_t         GenTreeParticleY[4];   //[GenTreeParticle]
   Float_t         GenTreeParticleMass[4];   //[GenTreeParticle]
   Int_t           GenTreeParticlePdgId[4];   //[GenTreeParticle]
   Int_t           GenTreeParticleStatus[4];   //[GenTreeParticle]
   Int_t           GenTreeParticleDecayLmm[4];   //[GenTreeParticle]
   Int_t           GenTreeParticleM1PdgId[4];   //[GenTreeParticle]
   Int_t           GenTreeParticleM2PdgId[4];   //[GenTreeParticle]
   Int_t           GenTreeParticleD1PdgId[4];   //[GenTreeParticle]
   Float_t         GenTreeParticlem1px[4];   //[GenTreeParticle]
   Int_t           Particle;
   Float_t         ParticleE[1019];   //[Particle]
   Float_t         ParticlePx[1019];   //[Particle]
   Float_t         ParticlePy[1019];   //[Particle]
   Float_t         ParticlePz[1019];   //[Particle]
   Float_t         ParticleX[1019];   //[Particle]
   Float_t         ParticleY[1019];   //[Particle]
   Float_t         ParticleMass[1019];   //[Particle]
   Int_t           ParticlePdgId[1019];   //[Particle]
   Int_t           ParticleStatus[1019];   //[Particle]
   Int_t           ParticleDecayLmm[1019];   //[Particle]
   Int_t           ParticleM1PdgId[1019];   //[Particle]
   Int_t           ParticleM2PdgId[1019];   //[Particle]
   Int_t           ParticleD1PdgId[1019];   //[Particle]
   Float_t         Particlem1px[1019];   //[Particle]

   // List of branches
   TBranch        *b_Muon;   //!
   TBranch        *b_MuonE;   //!
   TBranch        *b_MuonPx;   //!
   TBranch        *b_MuonPy;   //!
   TBranch        *b_MuonPz;   //!
   TBranch        *b_MuonX;   //!
   TBranch        *b_MuonY;   //!
   TBranch        *b_MuonMass;   //!
   TBranch        *b_MuonPdgId;   //!
   TBranch        *b_MuonStatus;   //!
   TBranch        *b_MuonDecayLmm;   //!
   TBranch        *b_MuonM1PdgId;   //!
   TBranch        *b_MuonM2PdgId;   //!
   TBranch        *b_MuonD1PdgId;   //!
   TBranch        *b_Muonm1px;   //!
   TBranch        *b_Electron;   //!
   TBranch        *b_ElectronE;   //!
   TBranch        *b_ElectronPx;   //!
   TBranch        *b_ElectronPy;   //!
   TBranch        *b_ElectronPz;   //!
   TBranch        *b_ElectronX;   //!
   TBranch        *b_ElectronY;   //!
   TBranch        *b_ElectronMass;   //!
   TBranch        *b_ElectronPdgId;   //!
   TBranch        *b_ElectronStatus;   //!
   TBranch        *b_ElectronDecayLmm;   //!
   TBranch        *b_ElectronM1PdgId;   //!
   TBranch        *b_ElectronM2PdgId;   //!
   TBranch        *b_ElectronD1PdgId;   //!
   TBranch        *b_Electronm1px;   //!
   TBranch        *b_Tau;   //!
   TBranch        *b_TauE;   //!
   TBranch        *b_TauPx;   //!
   TBranch        *b_TauPy;   //!
   TBranch        *b_TauPz;   //!
   TBranch        *b_TauX;   //!
   TBranch        *b_TauY;   //!
   TBranch        *b_TauMass;   //!
   TBranch        *b_TauPdgId;   //!
   TBranch        *b_TauStatus;   //!
   TBranch        *b_TauDecayLmm;   //!
   TBranch        *b_TauM1PdgId;   //!
   TBranch        *b_TauM2PdgId;   //!
   TBranch        *b_TauD1PdgId;   //!
   TBranch        *b_Taum1px;   //!
   TBranch        *b_b;   //!
   TBranch        *b_bE;   //!
   TBranch        *b_bPx;   //!
   TBranch        *b_bPy;   //!
   TBranch        *b_bPz;   //!
   TBranch        *b_bX;   //!
   TBranch        *b_bY;   //!
   TBranch        *b_bMass;   //!
   TBranch        *b_bPdgId;   //!
   TBranch        *b_bStatus;   //!
   TBranch        *b_bDecayLmm;   //!
   TBranch        *b_bM1PdgId;   //!
   TBranch        *b_bM2PdgId;   //!
   TBranch        *b_bD1PdgId;   //!
   TBranch        *b_bm1px;   //!
   TBranch        *b_c;   //!
   TBranch        *b_cE;   //!
   TBranch        *b_cPx;   //!
   TBranch        *b_cPy;   //!
   TBranch        *b_cPz;   //!
   TBranch        *b_cX;   //!
   TBranch        *b_cY;   //!
   TBranch        *b_cMass;   //!
   TBranch        *b_cPdgId;   //!
   TBranch        *b_cStatus;   //!
   TBranch        *b_cDecayLmm;   //!
   TBranch        *b_cM1PdgId;   //!
   TBranch        *b_cM2PdgId;   //!
   TBranch        *b_cD1PdgId;   //!
   TBranch        *b_cm1px;   //!
   TBranch        *b_Photon;   //!
   TBranch        *b_PhotonE;   //!
   TBranch        *b_PhotonPx;   //!
   TBranch        *b_PhotonPy;   //!
   TBranch        *b_PhotonPz;   //!
   TBranch        *b_PhotonX;   //!
   TBranch        *b_PhotonY;   //!
   TBranch        *b_PhotonMass;   //!
   TBranch        *b_PhotonPdgId;   //!
   TBranch        *b_PhotonStatus;   //!
   TBranch        *b_PhotonDecayLmm;   //!
   TBranch        *b_PhotonM1PdgId;   //!
   TBranch        *b_PhotonM2PdgId;   //!
   TBranch        *b_PhotonD1PdgId;   //!
   TBranch        *b_Photonm1px;   //!
   TBranch        *b_Neutrino;   //!
   TBranch        *b_NeutrinoE;   //!
   TBranch        *b_NeutrinoPx;   //!
   TBranch        *b_NeutrinoPy;   //!
   TBranch        *b_NeutrinoPz;   //!
   TBranch        *b_NeutrinoX;   //!
   TBranch        *b_NeutrinoY;   //!
   TBranch        *b_NeutrinoMass;   //!
   TBranch        *b_NeutrinoPdgId;   //!
   TBranch        *b_NeutrinoStatus;   //!
   TBranch        *b_NeutrinoDecayLmm;   //!
   TBranch        *b_NeutrinoM1PdgId;   //!
   TBranch        *b_NeutrinoM2PdgId;   //!
   TBranch        *b_NeutrinoD1PdgId;   //!
   TBranch        *b_Neutrinom1px;   //!
   TBranch        *b_SUSY;   //!
   TBranch        *b_SUSYE;   //!
   TBranch        *b_SUSYPx;   //!
   TBranch        *b_SUSYPy;   //!
   TBranch        *b_SUSYPz;   //!
   TBranch        *b_SUSYX;   //!
   TBranch        *b_SUSYY;   //!
   TBranch        *b_SUSYMass;   //!
   TBranch        *b_SUSYPdgId;   //!
   TBranch        *b_SUSYStatus;   //!
   TBranch        *b_SUSYDecayLmm;   //!
   TBranch        *b_SUSYM1PdgId;   //!
   TBranch        *b_SUSYM2PdgId;   //!
   TBranch        *b_SUSYD1PdgId;   //!
   TBranch        *b_SUSYm1px;   //!
   TBranch        *b_GenTreeParticle;   //!
   TBranch        *b_GenTreeParticleE;   //!
   TBranch        *b_GenTreeParticlePx;   //!
   TBranch        *b_GenTreeParticlePy;   //!
   TBranch        *b_GenTreeParticlePz;   //!
   TBranch        *b_GenTreeParticleX;   //!
   TBranch        *b_GenTreeParticleY;   //!
   TBranch        *b_GenTreeParticleMass;   //!
   TBranch        *b_GenTreeParticlePdgId;   //!
   TBranch        *b_GenTreeParticleStatus;   //!
   TBranch        *b_GenTreeParticleDecayLmm;   //!
   TBranch        *b_GenTreeParticleM1PdgId;   //!
   TBranch        *b_GenTreeParticleM2PdgId;   //!
   TBranch        *b_GenTreeParticleD1PdgId;   //!
   TBranch        *b_GenTreeParticlem1px;   //!
   TBranch        *b_Particle;   //!
   TBranch        *b_ParticleE;   //!
   TBranch        *b_ParticlePx;   //!
   TBranch        *b_ParticlePy;   //!
   TBranch        *b_ParticlePz;   //!
   TBranch        *b_ParticleX;   //!
   TBranch        *b_ParticleY;   //!
   TBranch        *b_ParticleMass;   //!
   TBranch        *b_ParticlePdgId;   //!
   TBranch        *b_ParticleStatus;   //!
   TBranch        *b_ParticleDecayLmm;   //!
   TBranch        *b_ParticleM1PdgId;   //!
   TBranch        *b_ParticleM2PdgId;   //!
   TBranch        *b_ParticleD1PdgId;   //!
   TBranch        *b_Particlem1px;   //!

   DetectorBase(TTree *tree=0);
   virtual ~DetectorBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DetectorBase_cxx
DetectorBase::DetectorBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("BSMGenTemplate.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("BSMGenTemplate.root");
      }
      f->GetObject("GenEvent",tree);

   }
   Init(tree);
}

DetectorBase::~DetectorBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DetectorBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DetectorBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DetectorBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Muon", &Muon, &b_Muon);
   fChain->SetBranchAddress("MuonE", MuonE, &b_MuonE);
   fChain->SetBranchAddress("MuonPx", MuonPx, &b_MuonPx);
   fChain->SetBranchAddress("MuonPy", MuonPy, &b_MuonPy);
   fChain->SetBranchAddress("MuonPz", MuonPz, &b_MuonPz);
   fChain->SetBranchAddress("MuonX", MuonX, &b_MuonX);
   fChain->SetBranchAddress("MuonY", MuonY, &b_MuonY);
   fChain->SetBranchAddress("MuonMass", MuonMass, &b_MuonMass);
   fChain->SetBranchAddress("MuonPdgId", MuonPdgId, &b_MuonPdgId);
   fChain->SetBranchAddress("MuonStatus", MuonStatus, &b_MuonStatus);
   fChain->SetBranchAddress("MuonDecayLmm", MuonDecayLmm, &b_MuonDecayLmm);
   fChain->SetBranchAddress("MuonM1PdgId", MuonM1PdgId, &b_MuonM1PdgId);
   fChain->SetBranchAddress("MuonM2PdgId", MuonM2PdgId, &b_MuonM2PdgId);
   fChain->SetBranchAddress("MuonD1PdgId", MuonD1PdgId, &b_MuonD1PdgId);
   fChain->SetBranchAddress("Muonm1px", Muonm1px, &b_Muonm1px);
   fChain->SetBranchAddress("Electron", &Electron, &b_Electron);
   fChain->SetBranchAddress("ElectronE", ElectronE, &b_ElectronE);
   fChain->SetBranchAddress("ElectronPx", ElectronPx, &b_ElectronPx);
   fChain->SetBranchAddress("ElectronPy", ElectronPy, &b_ElectronPy);
   fChain->SetBranchAddress("ElectronPz", ElectronPz, &b_ElectronPz);
   fChain->SetBranchAddress("ElectronX", ElectronX, &b_ElectronX);
   fChain->SetBranchAddress("ElectronY", ElectronY, &b_ElectronY);
   fChain->SetBranchAddress("ElectronMass", ElectronMass, &b_ElectronMass);
   fChain->SetBranchAddress("ElectronPdgId", ElectronPdgId, &b_ElectronPdgId);
   fChain->SetBranchAddress("ElectronStatus", ElectronStatus, &b_ElectronStatus);
   fChain->SetBranchAddress("ElectronDecayLmm", ElectronDecayLmm, &b_ElectronDecayLmm);
   fChain->SetBranchAddress("ElectronM1PdgId", ElectronM1PdgId, &b_ElectronM1PdgId);
   fChain->SetBranchAddress("ElectronM2PdgId", ElectronM2PdgId, &b_ElectronM2PdgId);
   fChain->SetBranchAddress("ElectronD1PdgId", ElectronD1PdgId, &b_ElectronD1PdgId);
   fChain->SetBranchAddress("Electronm1px", Electronm1px, &b_Electronm1px);
   fChain->SetBranchAddress("Tau", &Tau, &b_Tau);
   fChain->SetBranchAddress("TauE", TauE, &b_TauE);
   fChain->SetBranchAddress("TauPx", TauPx, &b_TauPx);
   fChain->SetBranchAddress("TauPy", TauPy, &b_TauPy);
   fChain->SetBranchAddress("TauPz", TauPz, &b_TauPz);
   fChain->SetBranchAddress("TauX", TauX, &b_TauX);
   fChain->SetBranchAddress("TauY", TauY, &b_TauY);
   fChain->SetBranchAddress("TauMass", TauMass, &b_TauMass);
   fChain->SetBranchAddress("TauPdgId", TauPdgId, &b_TauPdgId);
   fChain->SetBranchAddress("TauStatus", TauStatus, &b_TauStatus);
   fChain->SetBranchAddress("TauDecayLmm", TauDecayLmm, &b_TauDecayLmm);
   fChain->SetBranchAddress("TauM1PdgId", TauM1PdgId, &b_TauM1PdgId);
   fChain->SetBranchAddress("TauM2PdgId", TauM2PdgId, &b_TauM2PdgId);
   fChain->SetBranchAddress("TauD1PdgId", TauD1PdgId, &b_TauD1PdgId);
   fChain->SetBranchAddress("Taum1px", Taum1px, &b_Taum1px);
   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("bE", bE, &b_bE);
   fChain->SetBranchAddress("bPx", bPx, &b_bPx);
   fChain->SetBranchAddress("bPy", bPy, &b_bPy);
   fChain->SetBranchAddress("bPz", bPz, &b_bPz);
   fChain->SetBranchAddress("bX", bX, &b_bX);
   fChain->SetBranchAddress("bY", bY, &b_bY);
   fChain->SetBranchAddress("bMass", bMass, &b_bMass);
   fChain->SetBranchAddress("bPdgId", bPdgId, &b_bPdgId);
   fChain->SetBranchAddress("bStatus", bStatus, &b_bStatus);
   fChain->SetBranchAddress("bDecayLmm", bDecayLmm, &b_bDecayLmm);
   fChain->SetBranchAddress("bM1PdgId", bM1PdgId, &b_bM1PdgId);
   fChain->SetBranchAddress("bM2PdgId", bM2PdgId, &b_bM2PdgId);
   fChain->SetBranchAddress("bD1PdgId", bD1PdgId, &b_bD1PdgId);
   fChain->SetBranchAddress("bm1px", bm1px, &b_bm1px);
   fChain->SetBranchAddress("c", &c, &b_c);
   fChain->SetBranchAddress("cE", cE, &b_cE);
   fChain->SetBranchAddress("cPx", cPx, &b_cPx);
   fChain->SetBranchAddress("cPy", cPy, &b_cPy);
   fChain->SetBranchAddress("cPz", cPz, &b_cPz);
   fChain->SetBranchAddress("cX", cX, &b_cX);
   fChain->SetBranchAddress("cY", cY, &b_cY);
   fChain->SetBranchAddress("cMass", cMass, &b_cMass);
   fChain->SetBranchAddress("cPdgId", cPdgId, &b_cPdgId);
   fChain->SetBranchAddress("cStatus", cStatus, &b_cStatus);
   fChain->SetBranchAddress("cDecayLmm", cDecayLmm, &b_cDecayLmm);
   fChain->SetBranchAddress("cM1PdgId", cM1PdgId, &b_cM1PdgId);
   fChain->SetBranchAddress("cM2PdgId", cM2PdgId, &b_cM2PdgId);
   fChain->SetBranchAddress("cD1PdgId", cD1PdgId, &b_cD1PdgId);
   fChain->SetBranchAddress("cm1px", cm1px, &b_cm1px);
   fChain->SetBranchAddress("Photon", &Photon, &b_Photon);
   fChain->SetBranchAddress("PhotonE", PhotonE, &b_PhotonE);
   fChain->SetBranchAddress("PhotonPx", PhotonPx, &b_PhotonPx);
   fChain->SetBranchAddress("PhotonPy", PhotonPy, &b_PhotonPy);
   fChain->SetBranchAddress("PhotonPz", PhotonPz, &b_PhotonPz);
   fChain->SetBranchAddress("PhotonX", PhotonX, &b_PhotonX);
   fChain->SetBranchAddress("PhotonY", PhotonY, &b_PhotonY);
   fChain->SetBranchAddress("PhotonMass", PhotonMass, &b_PhotonMass);
   fChain->SetBranchAddress("PhotonPdgId", PhotonPdgId, &b_PhotonPdgId);
   fChain->SetBranchAddress("PhotonStatus", PhotonStatus, &b_PhotonStatus);
   fChain->SetBranchAddress("PhotonDecayLmm", PhotonDecayLmm, &b_PhotonDecayLmm);
   fChain->SetBranchAddress("PhotonM1PdgId", PhotonM1PdgId, &b_PhotonM1PdgId);
   fChain->SetBranchAddress("PhotonM2PdgId", PhotonM2PdgId, &b_PhotonM2PdgId);
   fChain->SetBranchAddress("PhotonD1PdgId", PhotonD1PdgId, &b_PhotonD1PdgId);
   fChain->SetBranchAddress("Photonm1px", Photonm1px, &b_Photonm1px);
   fChain->SetBranchAddress("Neutrino", &Neutrino, &b_Neutrino);
   fChain->SetBranchAddress("NeutrinoE", NeutrinoE, &b_NeutrinoE);
   fChain->SetBranchAddress("NeutrinoPx", NeutrinoPx, &b_NeutrinoPx);
   fChain->SetBranchAddress("NeutrinoPy", NeutrinoPy, &b_NeutrinoPy);
   fChain->SetBranchAddress("NeutrinoPz", NeutrinoPz, &b_NeutrinoPz);
   fChain->SetBranchAddress("NeutrinoX", NeutrinoX, &b_NeutrinoX);
   fChain->SetBranchAddress("NeutrinoY", NeutrinoY, &b_NeutrinoY);
   fChain->SetBranchAddress("NeutrinoMass", NeutrinoMass, &b_NeutrinoMass);
   fChain->SetBranchAddress("NeutrinoPdgId", NeutrinoPdgId, &b_NeutrinoPdgId);
   fChain->SetBranchAddress("NeutrinoStatus", NeutrinoStatus, &b_NeutrinoStatus);
   fChain->SetBranchAddress("NeutrinoDecayLmm", NeutrinoDecayLmm, &b_NeutrinoDecayLmm);
   fChain->SetBranchAddress("NeutrinoM1PdgId", NeutrinoM1PdgId, &b_NeutrinoM1PdgId);
   fChain->SetBranchAddress("NeutrinoM2PdgId", NeutrinoM2PdgId, &b_NeutrinoM2PdgId);
   fChain->SetBranchAddress("NeutrinoD1PdgId", NeutrinoD1PdgId, &b_NeutrinoD1PdgId);
   fChain->SetBranchAddress("Neutrinom1px", Neutrinom1px, &b_Neutrinom1px);
   fChain->SetBranchAddress("SUSY", &SUSY, &b_SUSY);
   fChain->SetBranchAddress("SUSYE", SUSYE, &b_SUSYE);
   fChain->SetBranchAddress("SUSYPx", SUSYPx, &b_SUSYPx);
   fChain->SetBranchAddress("SUSYPy", SUSYPy, &b_SUSYPy);
   fChain->SetBranchAddress("SUSYPz", SUSYPz, &b_SUSYPz);
   fChain->SetBranchAddress("SUSYX", SUSYX, &b_SUSYX);
   fChain->SetBranchAddress("SUSYY", SUSYY, &b_SUSYY);
   fChain->SetBranchAddress("SUSYMass", SUSYMass, &b_SUSYMass);
   fChain->SetBranchAddress("SUSYPdgId", SUSYPdgId, &b_SUSYPdgId);
   fChain->SetBranchAddress("SUSYStatus", SUSYStatus, &b_SUSYStatus);
   fChain->SetBranchAddress("SUSYDecayLmm", SUSYDecayLmm, &b_SUSYDecayLmm);
   fChain->SetBranchAddress("SUSYM1PdgId", SUSYM1PdgId, &b_SUSYM1PdgId);
   fChain->SetBranchAddress("SUSYM2PdgId", SUSYM2PdgId, &b_SUSYM2PdgId);
   fChain->SetBranchAddress("SUSYD1PdgId", SUSYD1PdgId, &b_SUSYD1PdgId);
   fChain->SetBranchAddress("SUSYm1px", SUSYm1px, &b_SUSYm1px);
   fChain->SetBranchAddress("GenTreeParticle", &GenTreeParticle, &b_GenTreeParticle);
   fChain->SetBranchAddress("GenTreeParticleE", GenTreeParticleE, &b_GenTreeParticleE);
   fChain->SetBranchAddress("GenTreeParticlePx", GenTreeParticlePx, &b_GenTreeParticlePx);
   fChain->SetBranchAddress("GenTreeParticlePy", GenTreeParticlePy, &b_GenTreeParticlePy);
   fChain->SetBranchAddress("GenTreeParticlePz", GenTreeParticlePz, &b_GenTreeParticlePz);
   fChain->SetBranchAddress("GenTreeParticleX", GenTreeParticleX, &b_GenTreeParticleX);
   fChain->SetBranchAddress("GenTreeParticleY", GenTreeParticleY, &b_GenTreeParticleY);
   fChain->SetBranchAddress("GenTreeParticleMass", GenTreeParticleMass, &b_GenTreeParticleMass);
   fChain->SetBranchAddress("GenTreeParticlePdgId", GenTreeParticlePdgId, &b_GenTreeParticlePdgId);
   fChain->SetBranchAddress("GenTreeParticleStatus", GenTreeParticleStatus, &b_GenTreeParticleStatus);
   fChain->SetBranchAddress("GenTreeParticleDecayLmm", GenTreeParticleDecayLmm, &b_GenTreeParticleDecayLmm);
   fChain->SetBranchAddress("GenTreeParticleM1PdgId", GenTreeParticleM1PdgId, &b_GenTreeParticleM1PdgId);
   fChain->SetBranchAddress("GenTreeParticleM2PdgId", GenTreeParticleM2PdgId, &b_GenTreeParticleM2PdgId);
   fChain->SetBranchAddress("GenTreeParticleD1PdgId", GenTreeParticleD1PdgId, &b_GenTreeParticleD1PdgId);
   fChain->SetBranchAddress("GenTreeParticlem1px", GenTreeParticlem1px, &b_GenTreeParticlem1px);
   fChain->SetBranchAddress("Particle", &Particle, &b_Particle);
   fChain->SetBranchAddress("ParticleE", ParticleE, &b_ParticleE);
   fChain->SetBranchAddress("ParticlePx", ParticlePx, &b_ParticlePx);
   fChain->SetBranchAddress("ParticlePy", ParticlePy, &b_ParticlePy);
   fChain->SetBranchAddress("ParticlePz", ParticlePz, &b_ParticlePz);
   fChain->SetBranchAddress("ParticleX", ParticleX, &b_ParticleX);
   fChain->SetBranchAddress("ParticleY", ParticleY, &b_ParticleY);
   fChain->SetBranchAddress("ParticleMass", ParticleMass, &b_ParticleMass);
   fChain->SetBranchAddress("ParticlePdgId", ParticlePdgId, &b_ParticlePdgId);
   fChain->SetBranchAddress("ParticleStatus", ParticleStatus, &b_ParticleStatus);
   fChain->SetBranchAddress("ParticleDecayLmm", ParticleDecayLmm, &b_ParticleDecayLmm);
   fChain->SetBranchAddress("ParticleM1PdgId", ParticleM1PdgId, &b_ParticleM1PdgId);
   fChain->SetBranchAddress("ParticleM2PdgId", ParticleM2PdgId, &b_ParticleM2PdgId);
   fChain->SetBranchAddress("ParticleD1PdgId", ParticleD1PdgId, &b_ParticleD1PdgId);
   fChain->SetBranchAddress("Particlem1px", Particlem1px, &b_Particlem1px);
   Notify();
}

Bool_t DetectorBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DetectorBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DetectorBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DetectorBase_cxx
