#include "CMS/CMSReco.hh"
#include "CMS/ParticleInfo.hh"

using namespace std;
using namespace stdcomb;

CMSReco::CMSReco(TTree *tree, bool delphesFormat) : DetectorReco(tree, delphesFormat) {
  dR = 0.5;
  cms = new CMSDetectorResponse(99999);
  _delphesFormat = delphesFormat;
}

CMSReco::~CMSReco() {
  delete cms;
}

void CMSReco::EleReco() {
  for(int i=0; i< Electron; i++) {
    if(sqrt(ElectronPx[i]*ElectronPx[i]+ElectronPy[i]*ElectronPy[i])<5.) continue;
    bool wp95 = false;
    if(EleSelector(i,"WP95") ) {
      EleWP95.push_back(cms->ElectronReco(TLorentzVector(ElectronPx[i],ElectronPy[i],ElectronPz[i],ElectronE[i])));
      idxEleWP95.push_back(i);
      wp95 = true;
    }
    if(EleSelector(i,"WP80") && wp95) {
      EleWP80.push_back(cms->ElectronReco(TLorentzVector(ElectronPx[i],ElectronPy[i],ElectronPz[i],ElectronE[i])));
      idxEleWP80.push_back(i);
    }
  }
}

void CMSReco::TrackReco(double minPt) {
  for(int i=0; i<Particle; i++) {
    // only charged particles
    if(IsCharged(ParticlePdgId[i])) {
      // we need to convert from mm to cm
      double d0 = sqrt(ParticleX[i]*ParticleX[i]+ParticleY[i]*ParticleY[i])/10.;
      // at least pT>500 MeV to reconstuct the track
      if(sqrt(ParticlePx[i]*ParticlePx[i]+ParticlePy[i]*ParticlePy[i])<minPt) continue;
      if(cms->TrackReco(d0)) {
	idxRecoTrack.push_back(i);
	Track.push_back(TLorentzVector(ParticlePx[i], ParticlePy[i],ParticlePz[i],ParticleE[i]));
      }
    }
  }
}

int CMSReco::MissingPixelHits(int iTrack) {
  double d0 = sqrt(ParticleX[iTrack]*ParticleX[iTrack]+ParticleY[iTrack]*ParticleY[iTrack])/10.;
  int missingLayers = 0;
  if(d0 > 4.) missingLayers++;
  if(d0 > 7.) missingLayers++;
  if(d0 > 11.) missingLayers++;
  int crossedLayers = 3-missingLayers;
  for(int i=0; i<crossedLayers; i++) if(!cms->PixelHit()) missingLayers++;
  return missingLayers;
}

void CMSReco::MuReco() {
  for(int i=0; i< Muon; i++) {
    if(sqrt(MuonPx[i]*MuonPx[i]+MuonPy[i]*MuonPy[i])<5.) continue;
    TLorentzVector myMuon(MuonPx[i],MuonPy[i],MuonPz[i],MuonE[i]);
    if(fabs(myMuon.Eta())>2.1) continue;
    bool looseMu = false;
    TLorentzVector recoMuon;
    if(MuonSelector(i,"Loose")) {
      recoMuon = cms->MuonReco(myMuon);
      LooseMu.push_back(recoMuon);
      idxLooseMu.push_back(i);
      if(MuonSelector(i,"Tight")) { 
	TightMu.push_back(recoMuon);
	idxTightMu.push_back(i);
      }
    }
  }
}

void CMSReco::PFReco() {

  if (_delphesFormat) {    
    // list of gmuons
    for(int i=0; i<Muon_size; i++) {    
      double Px = Muon_PT[i]*TMath::Cos(Muon_Phi[i]);
      double Py = Muon_PT[i]*TMath::Sin(Muon_Phi[i]);
      double Pz = Muon_PT[i]*TMath::SinH(Muon_Eta[i]);
      double E = Muon_PT[i]*TMath::CosH(Muon_Eta[i]);
      fastjet::PseudoJet p(Px, Py, Pz, E);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFMuons.push_back(p);
      //p.set_user_info(new ParticleInfo(Muon_PID[i]));
    }
  
    // list of electrons
    for(int i=0; i<Electron_size; i++) {    
      double Px = Electron_PT[i]*TMath::Cos(Electron_Phi[i]);
      double Py = Electron_PT[i]*TMath::Sin(Electron_Phi[i]);
      double Pz = Electron_PT[i]*TMath::SinH(Electron_Eta[i]);
      double E = Electron_PT[i]*TMath::CosH(Electron_Eta[i]);
      fastjet::PseudoJet p(Px, Py, Pz, E);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFElectrons.push_back(p);
      //p.set_user_info(new ParticleInfo(Electron_PID[i]));
    }

    // list of photons
    for(int i=0; i<Photon_size; i++) {    
      double Px = Photon_PT[i]*TMath::Cos(Photon_Phi[i]);
      double Py = Photon_PT[i]*TMath::Sin(Photon_Phi[i]);
      double Pz = Photon_PT[i]*TMath::SinH(Photon_Eta[i]);
      double E = Photon_PT[i]*TMath::CosH(Photon_Eta[i]);
      fastjet::PseudoJet p(Px, Py, Pz, E);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFPhotons.push_back(p);
      //p.set_user_info(new ParticleInfo(Photon_PID[i]));
    }

    // list of charged/neutral hadrons
    for(int i=0; i<Particle_size; i++) {
      fastjet::PseudoJet p(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
      if(p.pt()<0.5 || fabs(p.eta()) > 3.0) continue;
      if(abs(Particle_PID[i]) == 11) continue; // already included as electron
      if(abs(Particle_PID[i]) == 13) continue; // already included as muon
      if(abs(Particle_PID[i]) == 22) continue; // already included as photon
      if(IsCharged(Particle_PID[i])) {
	int size = _PFChHadrons.size();
	_PFChHadrons.push_back(p);
	_PFChHadrons[size].set_user_info(new ParticleInfo(Particle_PID[i]));
      }
      else {
	int size = _PFNeuHadrons.size();
	_PFNeuHadrons.push_back(p);
	_PFNeuHadrons[size].set_user_info(new ParticleInfo(Particle_PID[i]));
      }
    }  
  }
  else {
    // list of gen muons
    for(int i=0; i<Muon; i++) {
      fastjet::PseudoJet p(MuonPx[i], MuonPy[i], MuonPz[i], MuonE[i]);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFMuons.push_back(p);
      p.set_user_info(new ParticleInfo(MuonPdgId[i]));
    }
  
    // list of gen electrons
    for(int i=0; i<Electron; i++) {
      fastjet::PseudoJet p(ElectronPx[i], ElectronPy[i], ElectronPz[i], ElectronE[i]);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFElectrons.push_back(p);
      p.set_user_info(new ParticleInfo(ElectronPdgId[i]));
    }

    // list of gen photons
    for(int i=0; i<Photon; i++) {
      fastjet::PseudoJet p(PhotonPx[i], PhotonPy[i], PhotonPz[i], PhotonE[i]);
      if(p.pt()>0.5 && fabs(p.eta()) < 2.4) _PFPhotons.push_back(p);
      p.set_user_info(new ParticleInfo(PhotonPdgId[i]));
    }

    // list of gen charged/neutral hadrons
    for(int i=0; i<Particle; i++) {
      fastjet::PseudoJet p(ParticlePx[i], ParticlePy[i], ParticlePz[i], ParticleE[i]);
      if(p.pt()<0.5 || fabs(p.eta()) > 3.0) continue;
      if(abs(ParticlePdgId[i]) == 11) continue; // already included as electron
      if(abs(ParticlePdgId[i]) == 13) continue; // already included as muon
      if(abs(ParticlePdgId[i]) == 22) continue; // already included as photon
      if(IsCharged(ParticlePdgId[i])) {
	int size = _PFChHadrons.size();
	_PFChHadrons.push_back(p);
	_PFChHadrons[size].set_user_info(new ParticleInfo(ParticlePdgId[i]));
      }
      else {
	int size = _PFNeuHadrons.size();
	_PFNeuHadrons.push_back(p);
	_PFNeuHadrons[size].set_user_info(new ParticleInfo(ParticlePdgId[i]));
      }
    }
    // smear each particle category
  }

}


vector<fastjet::PseudoJet> CMSReco::PFJetConstituents(vector<fastjet::PseudoJet> muToRemove,
						      vector<fastjet::PseudoJet> eleToRemove,
						      vector<fastjet::PseudoJet> photonToRemove) {
  vector<fastjet::PseudoJet> PFcands;
  for(int i=0; i<_PFElectrons.size(); i++) {
    if(!FoundParticle(_PFElectrons[i], eleToRemove, 0.01)){
      int point = PFcands.size();
      PFcands.push_back(_PFElectrons[i]);
      PFcands[point].set_user_info(new ParticleInfo(ElectronPdgId[i]));
    }
  }
  for(int i=0; i<_PFMuons.size(); i++) {
    if(!FoundParticle(_PFMuons[i], muToRemove, 0.01)) {
      int point = PFcands.size();
      PFcands.push_back(_PFMuons[i]);
      PFcands[point].set_user_info(new ParticleInfo(MuonPdgId[i]));
    }
  }
  for(int i=0; i<_PFPhotons.size(); i++) {
    if(!FoundParticle(_PFPhotons[i], photonToRemove, 0.01)) {
      int point = PFcands.size();
      PFcands.push_back(_PFPhotons[i]);
      PFcands[point].set_user_info(new ParticleInfo(PhotonPdgId[i]));
    }
  }
  for(int i=0; i<_PFChHadrons.size(); i++) {
    int point = PFcands.size();
    PFcands.push_back(_PFChHadrons[i]);
    PFcands[point].set_user_info(new ParticleInfo(_PFChHadrons[i].user_info<ParticleInfo>().pdg_id));
  }
    
  for(int i=0; i<_PFNeuHadrons.size(); i++) {
    int point = PFcands.size();
    PFcands.push_back(_PFNeuHadrons[i]);
    PFcands[point].set_user_info(new ParticleInfo(_PFNeuHadrons[i].user_info<ParticleInfo>().pdg_id));
  }
  return PFcands;
  }

void CMSReco::SUSYReturn(int &susy1, vector<double> &susy2, vector<double> &susy3, vector<double> &susy4, vector<double> &susy5, vector<int> &susy6, vector<int> &susy7, vector<int> &susy8, vector<double> &susy9 ){
  susy1 = SUSY;
  susy2.reserve(SUSY);
  susy3.reserve(SUSY);
  susy4.reserve(SUSY);
  susy5.reserve(SUSY);
  susy6.reserve(SUSY);
  susy7.reserve(SUSY);
  susy8.reserve(SUSY);
  susy9.reserve(SUSY);
  
  for (int ni = 0; ni < SUSY; ni++){
    susy2[ni] = SUSYPx[ni];
    susy3[ni] = SUSYPy[ni];
    susy4[ni] = SUSYPz[ni];
    susy5[ni] = SUSYE[ni];
    susy6[ni] = SUSYPdgId[ni];
    susy7[ni] = SUSYM1PdgId[ni];
    susy8[ni] = SUSYStatus[ni];
    susy9[ni] = SUSYm1px[ni];

  }
}

void CMSReco::ParticleReturn(int &particle1, vector<double> &particle2, vector<double> &particle3, vector<double> &particle4, vector<double> &particle5, vector<int> &particle6, vector<int> &particle7, vector<int> &particle8, vector<double> &particle9 ){
  particle1 = Particle;
  particle2.reserve(Particle);
  particle3.reserve(Particle);
  particle4.reserve(Particle);
  particle5.reserve(Particle);
  particle6.reserve(Particle);
  particle7.reserve(Particle);
  particle8.reserve(Particle);
  particle9.reserve(Particle);
  
  for (int ni = 0; ni < Particle; ni++){
    particle2[ni] = ParticlePx[ni];
    particle3[ni] = ParticlePy[ni];
    particle4[ni] = ParticlePz[ni];
    particle5[ni] = ParticleE[ni];
    particle6[ni] = ParticlePdgId[ni];
    particle7[ni] = ParticleM1PdgId[ni];
    particle8[ni] = ParticleStatus[ni];
    particle9[ni] = Particlem1px[ni];

  }
}

void CMSReco::GenReturn(int &gen1, vector<double> &gen2, vector<double> &gen3, vector<double> &gen4, vector<double> &gen5, vector<int> &gen6, vector<int> &gen7, vector<int> &gen8, vector<double> &gen9 ){
	gen1 = GenTreeParticle; //this is fine
	
	gen2.reserve(GenTreeParticle);
	gen3.reserve(GenTreeParticle);
	gen4.reserve(GenTreeParticle);
	gen5.reserve(GenTreeParticle);
	gen6.reserve(GenTreeParticle);
	gen7.reserve(GenTreeParticle);
	gen8.reserve(GenTreeParticle);
	gen9.reserve(GenTreeParticle);
	
	for (int ni = 0; ni < GenTreeParticle; ni++) {
		gen2[ni] = GenTreeParticlePx[ni];
		gen3[ni] = GenTreeParticlePy[ni];
		gen4[ni] = GenTreeParticlePz[ni];
		gen5[ni] = GenTreeParticleE[ni];
		gen6[ni] = GenTreeParticlePdgId[ni];
		gen7[ni] = GenTreeParticleM1PdgId[ni];
		gen8[ni] = GenTreeParticleStatus[ni];
		gen9[ni] = GenTreeParticlem1px[ni];
	}
}






bool CMSReco::FoundParticle(fastjet::PseudoJet p, vector<fastjet::PseudoJet> q, double dR) {
  bool found = false;
  for(int i=0; i<q.size(); i++) {
    if(p.delta_R(p)<dR) found = true;
  }
  return found;
}

void CMSReco::CaloMETReco() {
  CaloMET = ConvertToPseudoJet(cms->CaloMETReco(ConvertTo4Vector(genMET_nomuon)));
}

void CMSReco::PFMETReco() {
  PFMET = ConvertToPseudoJet(cms->PFMETReco(ConvertTo4Vector(genMET)));
  PFMET_nomuon =  ConvertToPseudoJet(cms->PFMETReco(ConvertTo4Vector(genMET_nomuon)));
}

// The isolation is the only variable correlated to kinematic
// the rest of the selection can be described by an inefficiency factor
bool CMSReco::EleSelector(int iEle, string Selector){
  if(Selector != "WP80" && Selector != "WP95") {
    cout << "Error in CMSReco::EleSelector: selector " << Selector << " not implemented. Returning TRUE by default" << endl;
    return true;
  }
  TLorentzVector Ele(ElectronPx[iEle], ElectronPy[iEle],ElectronPz[iEle],ElectronE[iEle]);
  // ele eff vs gen-level particle-based isolation
  double genIso = RelGenIso(Ele, 0.3, true);
  return cms->EleEff(Ele, genIso, Selector);
}

// The isolation and the transverse impact parameter are the only variables correlated to kinematic
// the rest of the selection can be described by an inefficiency factor
bool CMSReco::MuonSelector(int iMu, string Selector) {
  if(Selector != "Tight" && Selector != "Loose") {
    cout << "Error in CMSReco::MuonSelector: Selector " << Selector << " not implemented. Returning TRUE by default" << endl;
    return true;
  }
  TLorentzVector Mu(MuonPx[iMu], MuonPy[iMu], MuonPz[iMu], MuonE[iMu]);  
  // apply muon efficiency vs isolation
  double genIso = RelGenIso(Mu, 0.3, true);
  return cms->MuonEff(Mu, genIso, Selector);
}

fastjet::PseudoJet CMSReco::CalcMHT(vector<fastjet::PseudoJet> jets){
  double px=0.;
  double py=0.;
  for(int i=0; i<jets.size(); i++) {
    px += jets[i].px();
    py += jets[i].py();
  }
  return fastjet::PseudoJet(-px,-py, 0., sqrt(px*px+py*py));
}

double CMSReco::CalcHT(vector<fastjet::PseudoJet> jets){
  double HT=0;
  for(int i=0; i<jets.size(); i++) {
    HT += jets[i].pt();
  }
  return HT;
}

double CMSReco::CalcAlphaT(fastjet::PseudoJet ja, fastjet::PseudoJet jb){
  fastjet::PseudoJet jT = ja+jb;
  double MT = sqrt(jT.Et()*jT.Et()-
		   jT.px()*jT.px()-
		   jT.py()*jT.py());
  return jb.Et()/MT;
}

double CMSReco::CalcMR(fastjet::PseudoJet j1, fastjet::PseudoJet j2){
  TLorentzVector ja,jb;
  ja.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), sqrt(j1.pz()*j1.pz()+j1.perp2()));
  jb.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), sqrt(j2.pz()*j2.pz()+j2.perp2()));
  return sqrt((ja.P()+jb.P())*(ja.P()+jb.P())-(ja.Pz()+jb.Pz())*(ja.Pz()+jb.Pz()));
}

double CMSReco::CalcSqrtsR(fastjet::PseudoJet j1, fastjet::PseudoJet j2, fastjet::PseudoJet met){
  TLorentzVector ja,jb;
  ja.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), sqrt(j1.pz()*j1.pz()+j1.perp2()));
  jb.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), sqrt(j2.pz()*j2.pz()+j2.perp2()));

  TLorentzVector pTcm(ja.Px()+jb.Px()+met.px(), ja.Py()+jb.Py()+met.py(), 0., 0.);
  double MR = CalcMR(j1,j2);
  double Einv = sqrt(met.perp2()+(ja+jb).M2());

  double sqrtsR = sqrt(2.)*sqrt( MR*MR + pTcm.Dot(ja+jb) + MR*sqrt(MR*MR + pTcm.Perp2() + 2.*pTcm.Dot(ja+jb)) );
  //cout << "MR+Einv = " << MR+Einv << endl;
  //cout << "sqrtsR = " << sqrtsR << endl;
  return sqrtsR;
}

double CMSReco::CalcGammaRp1Ana(fastjet::PseudoJet j1, fastjet::PseudoJet j2, fastjet::PseudoJet myMet){

  //Reconstructed leptons and missing transverse energy
  TLorentzVector L1,L2;
  L1.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), sqrt(j1.pz()*j1.pz()+j1.perp2()));
  L2.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), sqrt(j2.pz()*j2.pz()+j2.perp2()));
  TVector3 MET;
  MET.SetXYZ(myMet.px(),myMet.py(),0.);

  // MR

  double MR = CalcMR(L1,L2);
  double gamma_Rp1_Ana = 1./sqrt(1. - (L1-L2).Perp2()/MR/MR - 4.*(L2.E()*L1.Pz() - L1.E()*L2.Pz())*(L2.E()*L1.Pz() - L1.E()*L2.Pz())/MR/MR/MR/MR );

  return gamma_Rp1_Ana;

}

double CMSReco::CalcGammaRp1(fastjet::PseudoJet j1, fastjet::PseudoJet j2, fastjet::PseudoJet myMet){

  //Reconstructed leptons and missing transverse energy
  TLorentzVector L1,L2;
  L1.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), sqrt(j1.pz()*j1.pz()+j1.perp2()));
  L2.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), sqrt(j2.pz()*j2.pz()+j2.perp2()));
  TVector3 MET;
  MET.SetXYZ(myMet.px(),myMet.py(),0.);

  TVector3 vBETA_z = (1./(L1.E()+L2.E()))*(L1+L2).Vect(); 
  vBETA_z.SetX(0.0);         
  vBETA_z.SetY(0.0);

  //transformation from lab frame to approximate rest frame along beam-axis
  L1.Boost(-vBETA_z);
  L2.Boost(-vBETA_z);

  TVector3 pT_CM = (L1+L2).Vect() + MET;
  pT_CM.SetZ(0.0);     

  TLorentzVector LL = L1+L2;
  double SHATR = sqrt( 2.*(LL.E()*LL.E() - LL.Vect().Dot(pT_CM) + LL.E()*sqrt( LL.E()*LL.E() + pT_CM.Mag2() - 2.*LL.Vect().Dot(pT_CM) )));

  TVector3 vBETA_R = (1./sqrt(pT_CM.Mag2() + SHATR*SHATR))*pT_CM;

  double gamma_R = 1./sqrt(1.-vBETA_R.Mag2());

  //transformation from lab frame to R frame
  L1.Boost(-vBETA_R);
  L2.Boost(-vBETA_R);
  LL.Boost(-vBETA_R);  

  /////////////
  //
  // R-frame
  //
  /////////////

  double dphi_BETA_R = fabs((LL.Vect()).DeltaPhi(vBETA_R));

  double dphi_L1_L2 = fabs(L1.Vect().DeltaPhi(L2.Vect()));

  TVector3 vBETA_Rp1 = (1./(L1.E()+L2.E()))*(L1.Vect() - L2.Vect());

  double gamma_Rp1 = 1./sqrt(1.-vBETA_Rp1.Mag2());

  return gamma_Rp1;
}

double CMSReco::CalcMR_zinvariant(fastjet::PseudoJet j1, fastjet::PseudoJet j2){
	TLorentzVector ja,jb;
	ja.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), j1.E());//  E was originally set to sqrt(j1.pz()*j1.pz()+j1.perp2()));
	jb.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), j2.E());//  E was originally set to sqrt(j2.pz()*j2.pz()+j2.perp2()));
	return sqrt((ja.E()+jb.E())*(ja.E()+jb.E())-(ja.Pz()+jb.Pz())*(ja.Pz()+jb.Pz()));
}


double CMSReco::CalcMRT(fastjet::PseudoJet j1, fastjet::PseudoJet j2, fastjet::PseudoJet met){
  TLorentzVector ja,jb;
  ja.SetPtEtaPhiE(j1.pt(), j1.eta(), j1.phi(), sqrt(j1.pz()*j1.pz()+j1.perp2()));
  jb.SetPtEtaPhiE(j2.pt(), j2.eta(), j2.phi(), sqrt(j2.pz()*j2.pz()+j2.perp2()));
  TLorentzVector tSum(ja.Px()+jb.Px(), ja.Py()+jb.Py(), 0., 0.);
  return sqrt((met.pt()*(ja.Pt()+jb.Pt()) - met.px()*tSum.Px()-met.py()*tSum.Py())/2.);
}

// apply btag efficiency  
bool CMSReco::BTagCSVM(fastjet::PseudoJet j, double dR) {
  bool btag = true;
  if(IsBJet(j, dR, 20.)) {
      btag = cms->BTagCSVM_b(ConvertTo4Vector(j));
  } else if(IsCJet(j, dR, 20.)) {
    btag = cms->BTagCSVM_c(ConvertTo4Vector(j));
  } else {
    btag = cms->BTagCSVM_udsg(ConvertTo4Vector(j));
  }	
  return btag;
}

// apply btag efficiency  
bool CMSReco::BTagHiggs(fastjet::PseudoJet j, double dR) {
  bool btag = true;
  if(IsHbbJet(j, dR)) {
    btag = cms->BTagHbb_Hbb(ConvertTo4Vector(j));
  } else if(IsZbbJet(j, dR)) {
    btag = cms->BTagHbb_Zbb(ConvertTo4Vector(j));
  } else if(IsTopJet(j, dR)) {
    btag = cms->BTagHbb_Top(ConvertTo4Vector(j));
  } else if(IsWcsJet(j, dR)) {
    btag = cms->BTagHbb_Wcs(ConvertTo4Vector(j));
  } else {
    btag = cms->BTagHbb_udcsg(ConvertTo4Vector(j));
  }	
  return btag;
}

void CMSReco::CleanEvent() {
  idxEleWP80.clear();
  EleWP80.clear();
  idxEleWP95.clear();
  EleWP95.clear();
  idxTightMu.clear();
  TightMu.clear();
  idxLooseMu.clear();
  LooseMu.clear();
  idxRecoTrack.clear();
  Track.clear();
  _PFElectrons.clear();
  _PFMuons.clear();
  _PFPhotons.clear();
  _PFChHadrons.clear();
  _PFNeuHadrons.clear();
}
