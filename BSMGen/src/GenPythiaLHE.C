#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include <unistd.h>
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"

#include <TFile.h>
#include <TTree.h>
#include <TDatime.h>
#include <GenTree.hh>
#include <GenCandidateFiller.hh>
#include <EventFilter.hh>

#include <stdio.h>
#include <string.h>

#define _debug 0

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings (.cmnd)" << endl
         << " 2. Full name of the input LHE file (with path) (.lhe)" << endl
         << " 3. Outut files" << endl
         << " Program stopped. " << endl;
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
 // ToHepMC.set_free_parton_exception(false);

  // Confirm that external files will be used for input and output.
  cout << " PYTHIA settings will be read from file " << argv[1] << endl;
  
  std::string cfg = string(argv[1]);
  std::string iPath = string(argv[2]);
  std::string outfilename = string(argv[3]);
  bool boolFilter = false;
  cout << " boolFilter = " << boolFilter << endl;
  
  // Specify file where HepMC events will be stored.
  char hepmcname[256];
  sprintf(hepmcname,"%s.hepmc", outfilename.c_str());  
  HepMC::IO_GenEvent ascii_io(hepmcname, std::ios::out);
 
  if (_debug) cout << "[DEBUG] Prepare the generator...\n";
  // Generator. 
  Pythia pythia;

  // Read in commands from external file.
  if (_debug) cout << "[DEBUG] Reading from " << cfg << endl;
  pythia.readFile(cfg); 
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + iPath); 
  
  // set seed                                                                                                             
  int jobpid = getpid();
  TDatime *now = new TDatime();
  int today = now->GetDate();
  int clock = now->GetTime();
  int myseed = today+clock+jobpid*1000;
  if(myseed>900000000) myseed = myseed - 900000000;
  pythia.readString("Random:setSeed=on");
  char command[512];
  sprintf(command,"Random:seed=%i",myseed);
  pythia.readString(command);
  pythia.readString("Next:numberShowInfo = 1");
  pythia.readString("Next:numberShowProcess = 1");
  pythia.readString("Next:numberShowEvent = 1");
//
  // Initialize. Beam parameters set in .pythia file.
  if (_debug) cout << "[DEBUG] Initializing pythia... " << cfg << endl;
  pythia.init();
  
  //pythia.particleData.doForceWidth(1000005,true);
  //pythia.particleData.mayDecay(1000005, false);
  
  // Values for filter.
  int    pdgId   = 1000022; // ask for neutralinos
  int    pdgMothId   = 100001; // ask for squarks mother. (Disabled)
  double etaMax   = -1; // no requirement on eta
  double pTmin = 50; // MET > 100 GeV

  // Declare Event Filter according to specification.
  EventFilter filter( pdgId, pdgMothId, etaMax, pTmin );
  
  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  //int nList    = pythia.mode("Next:numberShowEvent");
  int nShow    = pythia.mode("Next:numberCount");
  int nAbort   = pythia.mode("Main:timesAllowErrors"); 
  //bool showCS  = pythia.flag("Init:showChangedSettings");
  bool showAS  = pythia.flag("Init:showAllSettings");
  bool showCPD = pythia.flag("Init:showChangedParticleData");
  bool showAPD = pythia.flag("Init:showAllParticleData");  

  // List settings.
  //if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  
// Begin event loop.
  int nPace = max(1, nShow); 
  int iAbort = 0;
  
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    
    if (nShow > 0 && iEvent%nPace == 0)      
    cout << " Now begin event " << iEvent << endl;
    //  cout << " Now begin event " << iFilteredEvent << endl;

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (_debug) cout << "iAbort = " << iAbort << endl;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      cerr << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
	
    if (_debug) cout << "Doing fine after !pythia.next()" << endl; 
    
    // List first few events.
    //if (iEvent < nList) { 
    //  pythia.info.list();
    //  pythia.event.list();
    //}

    // Construct new empty HepMC event. 
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );
    

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;
    
  }
  cout << "nEvent = " << nEvent << endl;

  // Give statistics. 
  pythia.stat();
  
  // Done.
  return 0;
}
