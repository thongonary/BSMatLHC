#include "Pythia8/Pythia.h"
//#include "Pythia8/Pythia8ToHepMC.h"
#include "Pythia8Plugins/HepMC2.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"

#include <TFile.h>
#include <TTree.h>
#include <GenTree.hh>
#include <GenCandidateFiller.hh>
#include <EventFilter.hh>

#include <stdio.h>
#include <string.h>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc < 3) {
    cerr << " To run the code provide the name of the input pythia card and the output LHE file. \n"
	 << " example: ./GenPythia data/pythiaCards/EXO/RSGraviton_gg_EXAMPLE.pythia outFile [--filter] \n" << endl;
    return 1;
  }

  // Check that the provided input name corresponds to an existing file.
  ifstream is(argv[1]);  
  if (!is) {
    cerr << " Command-line file " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::Pythia8ToHepMC ToHepMC;
  

  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_warnings(false);

  // Confirm that external files will be used for input and output.
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;
  
  std::string cfg = argv[1];
  std::string outfilename = argv[2];
  bool boolFilter = false;
  for (int i = 3; i < argc; i++){    
    if (strncmp(argv[i],"--filter",8)==0) {
      boolFilter = true;
    }
  }
  cout << " boolFilter = " << boolFilter << endl;
  
  // Specify file where HepMC events will be stored.
  char hepmcname[256];
  sprintf(hepmcname,"%s.hepmc", outfilename.c_str());  
  HepMC::IO_GenEvent ascii_io(hepmcname, std::ios::out);
  
  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  pythia.readFile(cfg);  
  
  // Initialize Les Houches Event File run. List initialization information.
  LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  
  // Open a file on which LHEF events should be stored, and write header.  
  char lhename[256];
  sprintf(lhename,"%s.lhe", outfilename.c_str());
  myLHA.openLHEF(lhename);

  // Initialize. Beam parameters set in .pythia file.
  pythia.init();
  
  //pythia.particleData.doForceWidth(1000005,true);
  //pythia.particleData.mayDecay(1000005, false);
  
  // Values for filter.
  int    pdgId   = 22; //ask for photons
  int    pdgMothId   = 25; //ask for Higgs mother
  double etaMax   = -1; //no requirement on eta
  double pTmin = -1; //no requirement on pt

  // Declare Event Filter according to specification.
  EventFilter filter( pdgId, pdgMothId, etaMax, pTmin );
  
  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  //int nList    = pythia.mode("Main:numberToList");
  int nShow    = pythia.mode("Next:numberCount");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  //bool showCS  = pythia.flag("Main:showChangedSettings");
  bool showAS  = pythia.flag("Init:showAllSettings");
  bool showCPD = pythia.flag("Init:showChangedParticleData");
  bool showAPD = pythia.flag("Init:showAllParticleData");  

  // List settings.
  //if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // the output file 
  char name[256];
  sprintf(name,"%s_GenTree.root", outfilename.c_str());
  TFile* treeOut = new TFile(name,"recreate");

  // the output TTree with information on the model
  double xsec, filtereff;
  TTree* infoTree = new TTree("infoTree", "infoTree");
  infoTree->Branch("xsec", &xsec, "xsec/D");
  infoTree->Branch("filtereff", &filtereff, "filtereff/D");

  // the event TTree
  GenTree* myTree = new GenTree("GenEvent","GenEvent");
  GenCandidateFiller* muonFiller = new GenCandidateFiller(myTree,"Muon");
  GenCandidateFiller* electronFiller = new GenCandidateFiller(myTree,"Electron");
  GenCandidateFiller* tauFiller = new GenCandidateFiller(myTree,"Tau");
  GenCandidateFiller* bFiller = new GenCandidateFiller(myTree,"b");
  GenCandidateFiller* cFiller = new GenCandidateFiller(myTree,"c");
  GenCandidateFiller* photonFiller = new GenCandidateFiller(myTree,"Photon");
  GenCandidateFiller* neutrinoFiller = new GenCandidateFiller(myTree,"Neutrino");
  GenCandidateFiller* susyFiller = new GenCandidateFiller(myTree,"SUSY");
  GenCandidateFiller* gentreeparticleFiller = new GenCandidateFiller(myTree,"GenTreeParticle");
  GenCandidateFiller* particleFiller = new GenCandidateFiller(myTree,"Particle");

  int iFilteredEvent = 0;
  int nUnfilteredEvent = 0;
  
  // Begin event loop.
  int nPace = max(1, nShow); 
  int iAbort = 0;
  for (int iEvent = 0; ; ++iEvent) {
    
    if (iFilteredEvent >= nEvent) break;

    
    if (nShow > 0 && iEvent%nPace == 0)      
    //cout << " Now begin event " << iEvent << endl;
      cout << " Now begin event " << iFilteredEvent << endl;

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }

    nUnfilteredEvent++;    

    // Find final state photons
    filter.filter( pythia.process);
    
    //cout << "number of final state photons from Higgs = " << filter.size() << endl;
    if (filter.size()<2) {
      //cout << "< 2 final state photons from Higgs; not saving" << endl;
      //filter.list();
      //pythia.process.list();
      if (boolFilter) continue;
    }
    
    iFilteredEvent++;
			
			 
    // List first few events.
    //if (iEvent < nList) { 
    //  pythia.info.list();
    //  pythia.event.list();
    //}

    // Construct new empty HepMC event. 
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );
    
    // write event to root file
    muonFiller->FillEvent(hepmcevt);
    electronFiller->FillEvent(hepmcevt);
    tauFiller->FillEvent(hepmcevt);
    bFiller->FillEvent(hepmcevt);
    cFiller->FillEvent(hepmcevt);
    photonFiller->FillEvent(hepmcevt);
    neutrinoFiller->FillEvent(hepmcevt);
    susyFiller->FillEvent(hepmcevt);
    gentreeparticleFiller->FillEvent(hepmcevt);
    particleFiller->FillEvent(hepmcevt);

    // write data in TTree
    myTree->dumpData();

    // Clear the event from memory
    muonFiller->ClearEvent();
    electronFiller->ClearEvent();
    tauFiller->ClearEvent();
    bFiller->ClearEvent();
    cFiller->ClearEvent();
    photonFiller->ClearEvent();
    neutrinoFiller->ClearEvent();
    susyFiller->ClearEvent();
    gentreeparticleFiller->ClearEvent();
    particleFiller->ClearEvent();

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;
    
    // Store event info in the LHAup object.
    myLHA.setEvent();

    // Write out this event info on the file.
    // With optional argument (verbose =) false the file is smaller.
    myLHA.eventLHEF();
    
    // End of event loop.
  }
  cout << "nEvent = " << nEvent << endl;
  cout << "nUnfilteredEvent = " << nUnfilteredEvent << endl;
  filtereff = nEvent*1.0/nUnfilteredEvent;
  
  // Write output file
  infoTree->Fill();

  TTree* roottree = myTree->getTree();
  roottree->Write();
  infoTree->Write();
  treeOut->Close();

  // Give statistics. 
  pythia.stat();
  
  // Update the cross section info based on Monte Carlo integration during run.
  myLHA.updateSigma();

  // Write endtag. Overwrite initialization info with new cross sections.
  myLHA.closeLHEF(true);
    
  // Done.
  return 0;
}
