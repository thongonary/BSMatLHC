// Generate a sample of events, starting from 
// a sample of LHE and write the output in a ROOT file

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

#include <stdio.h>
#include <string.h>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments. \n You are"
         << " expected to provide the name of an input card, an input LHE file, and an output file. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Check that the provided input names corresponds to existing files
  ifstream is1(argv[1]);  
  if (!is1) {
    cerr << " The requested PYTHIA card " << argv[1] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }
  ifstream is2(argv[2]);  
  if (!is2) {
    cerr << " The requested LHE event file " << argv[2] << " was not found. \n"
         << " Program stopped! " << endl;
    return 1;
  }

  // Confirm that external files will be used for input and output.
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;
  
  std::string cfg = argv[1];
  std::string infilename = argv[2]; 
  std::string outfilename = argv[3]; 

  std::cout << "Using config " << cfg << std::endl;
  std::cout << "Will write decayed file to " << outfilename << std::endl;

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::Pythia8ToHepMC ToHepMC;
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_warnings(false);

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
  
  pythia.init();
  
  //pythia.init(infilename);

  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nList    = pythia.mode("Main:numberToList");
  int nShow    = pythia.mode("Main:timesToShow");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  bool showCS  = pythia.flag("Main:showChangedSettings");
  bool showAS  = pythia.flag("Main:showAllSettings");
  bool showCPD = pythia.flag("Main:showChangedParticleData");
  bool showAPD = pythia.flag("Main:showAllParticleData");
  

  // List settings.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // the output file 
  char name[256];
  sprintf(name,"%s_GenTree.root", outfilename.c_str());
  TFile* treeOut = new TFile(name,"recreate");

  // the output TTree with information on the model
  double xsec;
  TTree* infoTree = new TTree("infoTree", "infoTree");
  infoTree->Branch("xsec", &xsec, "xsec/D");
  
// the event TTree
  GenTree* myTree = new GenTree("GenEvent","GenEvent");
  GenCandidateFiller* muonFiller = new GenCandidateFiller(myTree,"Muon");
  GenCandidateFiller* electronFiller = new GenCandidateFiller(myTree,"Electron");
  GenCandidateFiller* tauFiller = new GenCandidateFiller(myTree,"Tau");
  GenCandidateFiller* bFiller = new GenCandidateFiller(myTree,"b");
  GenCandidateFiller* photonFiller = new GenCandidateFiller(myTree,"Photon");
  GenCandidateFiller* neutrinoFiller = new GenCandidateFiller(myTree,"Neutrino");
  GenCandidateFiller* susyFiller = new GenCandidateFiller(myTree,"SUSY");
  GenCandidateFiller* particleFiller = new GenCandidateFiller(myTree,"Particle");

  // Begin event loop.
  int nPace = max(1, nEvent / max(1, nShow) ); 
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (nShow > 0 && iEvent%nPace == 0) 
      cout << " Now begin event " << iEvent << endl;
    
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break; 
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }
 
    // List first few events.
    if (iEvent < nList) { 
      pythia.LHAeventList();               
      pythia.info.list();
      pythia.process.list();          
      pythia.event.list();
    }

    // Construct new empty HepMC event. 
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // write event to root file
    muonFiller->FillEvent(hepmcevt);
    electronFiller->FillEvent(hepmcevt);
    tauFiller->FillEvent(hepmcevt);
    bFiller->FillEvent(hepmcevt);
    photonFiller->FillEvent(hepmcevt);
    neutrinoFiller->FillEvent(hepmcevt);
    susyFiller->FillEvent(hepmcevt);
    particleFiller->FillEvent(hepmcevt);

    // write data in TTree
    myTree->dumpData();

    // Clear the event from memory
    muonFiller->ClearEvent();
    electronFiller->ClearEvent();
    tauFiller->ClearEvent();
    bFiller->ClearEvent();
    photonFiller->ClearEvent();
    neutrinoFiller->ClearEvent();
    susyFiller->ClearEvent();
    particleFiller->ClearEvent();
    delete hepmcevt;
    
    // Store event info in the LHAup object.
    myLHA.setEvent();

    // Write out this event info on the file.
    // With optional argument (verbose =) false the file is smaller.
    myLHA.eventLHEF();

    // End of event loop.
  }

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
