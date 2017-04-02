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

#include <unistd.h>

// Include UserHooks for Jet Matching.
#include "Pythia8Plugins/CombineMatchingInput.h"
// Include UserHooks for randomly choosing between integrated and
// non-integrated treatment for unitarised merging.
#include "Pythia8Plugins/aMCatNLOHooks.h"

using namespace Pythia8;

//==========================================================================

// Example main programm to illustrate merging.

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc != 4) {
    cerr << " Unexpected number of command-line arguments ("<<argc<<"). \n"
         << " You are expected to provide the arguments" << endl
         << " 1. Input file for settings" << endl
         << " 2. Input LHE file" << endl
         << " 3. Output file for HepMC events" << endl
         << " Program stopped. " << endl;
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
  
  Pythia pythia;

  // New setting to allow processing of multiple input LHEFs.
  pythia.settings.addMode("LHEFInputs:nSubruns",0,true,false,0,100);

  // Input parameters:
  pythia.readFile(cfg,0);

  // Interface for conversion from Pythia8::Event to HepMC one.
  HepMC::Pythia8ToHepMC ToHepMC;
  // Specify file where HepMC events will be stored.  
  char hepmcname[256];
  sprintf(hepmcname,"%s.hepmc", outfilename.c_str());  
  HepMC::IO_GenEvent ascii_io(hepmcname, std::ios::out);
  // Switch off warnings for parton-level events.
  ToHepMC.set_print_inconsistency(false);
  ToHepMC.set_free_parton_exception(false);
  // Do not store cross section information, as this will be done manually.
  ToHepMC.set_store_pdf(false);
  ToHepMC.set_store_proc(false);
  ToHepMC.set_store_xsec(false);

  // Check if jet matching should be applied.
  bool doMatch   = pythia.settings.flag("JetMatching:merge");

  
  if (doMatch) {
    cout << "\n Do MATCHING! \n" << endl;
  }
  
  // Check if internal merging should be applied.
  bool doMerge   = !(pythia.settings.word("Merging:Process").compare("void")
    == 0);

  // Currently, only one scheme at a time is allowed.
  if (doMatch && doMerge) {
    cerr << " Jet matching and merging cannot be used simultaneously.\n"
         << " Program stopped.";
  }

  // Get number of subruns.
  int nMerge = pythia.mode("LHEFInputs:nSubruns");

  // Number of events. Negative numbers mean all events in the LHEF will be
  // used.
  long nEvent = pythia.settings.mode("Main:numberOfEvents");
  if (nEvent < 1) nEvent = 1000000000000000;

  // For jet matching, initialise the respective user hooks code.
  CombineMatchingInput* combined = NULL;
  UserHooks* matching            = NULL;

  // Allow to set the number of addtional partons dynamically.
  amcnlo_unitarised_interface* setting = NULL;
  if ( doMerge ) {
    // Store merging scheme.
    int scheme = ( pythia.settings.flag("Merging:doUMEPSTree")
                || pythia.settings.flag("Merging:doUMEPSSubt")) ?
                1 :
                 ( ( pythia.settings.flag("Merging:doUNLOPSTree")
                || pythia.settings.flag("Merging:doUNLOPSSubt")
                || pythia.settings.flag("Merging:doUNLOPSLoop")
                || pythia.settings.flag("Merging:doUNLOPSSubtNLO")) ?
                2 :
                0 );
    setting = new amcnlo_unitarised_interface(scheme);
    pythia.setUserHooksPtr(setting);
  }

  // For jet matching, initialise the respective user hooks code.
  if (doMatch) {
    matching = combined->getHook(pythia);
    if (!matching) {
      cerr << " Failed to initialise jet matching structures.\n"
           << " Program stopped.";
      return 1;
    }
    pythia.setUserHooksPtr(matching);
  }

  // Cross section and error.
  double sigmaTotal  = 0.;
  double errorTotal  = 0.;

  // Allow abort of run if many errors.
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  int  iAbort  = 0;
  bool doAbort = false;

  cout << endl << endl << endl;
  cout << "Start generating events" << endl;

  // Loop over subruns with varying number of jets.
  for (int iMerge = 0; iMerge < nMerge; ++iMerge) {

    double sigmaSample = 0., errorSample = 0.;

    // Read in name of LHE file for current subrun and initialize.
    pythia.readFile(cfg, iMerge);
    // Initialise.
    pythia.init();

    // Get the inclusive x-section by summing over all process x-sections.
    double xs = 0.;
    for (int i=0; i < pythia.info.nProcessesLHEF(); ++i)
      xs += pythia.info.sigmaLHEF(i);

    // Start generation loop
    while( pythia.info.nSelected() < nEvent ){

      // Generate next event
      if( !pythia.next() ) {
        if ( pythia.info.atEndOfFile() ) break;
        else if (++iAbort > nAbort) {doAbort = true; break;}
        else continue;
      }

      // Get event weight(s).
      double evtweight         = pythia.info.weight();
      // Additional PDF/alphaS weight for internal merging.
      if (doMerge) evtweight  *= pythia.info.mergingWeightNLO()
      // Additional weight due to random choice of reclustered/non-reclustered
      // treatment. Also contains additional sign for subtractive samples.
                                *setting->getNormFactor();

      // Do not print zero-weight events.
      if ( evtweight == 0. ) continue;
      // Construct new empty HepMC event.
      HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();

      // Work with weighted (LHA strategy=-4) events.
      double normhepmc = 1.;
      if (abs(pythia.info.lhaStrategy()) == 4)
        normhepmc = 1. / double(1e9*nEvent);
      // Work with unweighted events.
      else
        normhepmc = xs / double(1e9*nEvent);

      // Set event weight
      hepmcevt->weights().push_back(evtweight*normhepmc);
      // Fill HepMC event
      ToHepMC.fill_next_event( pythia, hepmcevt );
      // Add the weight of the current event to the cross section.
      sigmaTotal  += evtweight*normhepmc;
      sigmaSample += evtweight*normhepmc;
      errorTotal  += pow2(evtweight*normhepmc);
      errorSample += pow2(evtweight*normhepmc);
      // Report cross section to hepmc
      HepMC::GenCrossSection xsec;
      xsec.set_cross_section( sigmaTotal*1e9, pythia.info.sigmaErr()*1e9 );
      hepmcevt->set_cross_section( xsec );
      // Write the HepMC event to file. Done with it.
      //      cout << "writing HepMC event" << endl;
      ascii_io << hepmcevt;
      delete hepmcevt;

    } // end loop over events to generate.
    if (doAbort) break;

    // print cross section, errors
    pythia.stat();

    cout << endl << " Contribution of sample " << iMerge
         << " to the inclusive cross section : "
         << scientific << setprecision(8)
         << sigmaSample << "  +-  " << sqrt(errorSample)  << endl;

  }

  cout << endl << endl << endl;
  if (doAbort)
    cout << " Run was not completed owing to too many aborted events" << endl;
  else
    cout << "Inclusive cross section: " << scientific << setprecision(8)
         << sigmaTotal << "  +-  " << sqrt(errorTotal) << " mb " << endl;
  cout << endl << endl << endl;

  // Done
  return 0;

}
