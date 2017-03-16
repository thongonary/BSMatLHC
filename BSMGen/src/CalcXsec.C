#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"

#include <TDatime.h>
#include <TH1D.h>
#include <TFile.h>

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fstream>

using namespace Pythia8; 

int main(int argc, char* argv[]) {

  // Check that correct number of command-line arguments
  if (argc != 6) {
    cerr << " To run the code provide the name of the input pythia card, the output ROOT file, the string containing the parameter to scan, the min and max values \n"
	 << " example: ./CalcXsec data/pythiaCards/EXO/RSGraviton_gg_EXAMPLE.pythia outFile.root 1000001:m0 10. 100. \n" << endl;
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

  float xMin = atof(argv[4]);
  float xMax = atof(argv[5]);
  double xsecVal[100];
  double errXsec[100];

  for(int i=0; i<100; i++) {
 
    // Generator. 
    Pythia pythia;

    // Read in commands from external file.
    pythia.readFile(argv[1]);    

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

    // Initialize. Beam parameters set in .pythia file.
    pythia.init();

    // Extract settings to be used in the main program.
    int nEvent   = 300;
    int nList    = pythia.mode("Next:numberShowEvent");
    int nShow    = pythia.mode("Next:numberCount");
    int nAbort   = pythia.mode("Main:timesAllowErrors"); 
    bool showCS  = pythia.flag("Init:showChangedSettings");
    bool showAS  = pythia.flag("Init:showAllSettings");
    bool showCPD = pythia.flag("Init:showChangedParticleData");
    bool showAPD = pythia.flag("Init:showAllParticleData");
  
    // List settings.
    //  if (showCS) pythia.settings.listChanged();
    //  if (showAS) pythia.settings.listAll();
    
    // List particle data.  
    //  if (showCPD) pythia.particleData.listChanged();
    //  if (showAPD) pythia.particleData.listAll();
  
    // set new mass
    double thisValue = xMin + (i+0.5)/100.*(xMax-xMin);
    sprintf(command,"%s=%f",argv[3], thisValue);
    //    cout << "READTHIS " << command << endl;
    pythia.readString(command);
    // re-initialize pythia
    //    pythia.init();

    // Begin event loop.
    int nPace = max(1, nShow ); 
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (nShow > 0 && iEvent%nPace == 0) 
	cout << " Now begin event " << iEvent << endl;
      
      // Generate events. Quit if many failures.
      if (!pythia.next()) {
	// First few failures write off as "acceptable" errors, then quit.
	if (++iAbort < nAbort) continue;
	cout << " Event generation aborted prematurely, owing to error!\n"; 
	break;
      }
    }
    xsecVal[i] = pythia.info.sigmaGen();
    errXsec[i] = pythia.info.sigmaErr();
    ParticleData& pdt = pythia.particleData;
//    cout << "READTHIS " << xsecVal[i] << " " << errXsec[i] << " " << pdt.m0(1000022) << endl;
  }
  
  TH1D* xsec = new TH1D("xsec", "", 100,  xMin, xMax);
  for(int i=0; i<100; i++) {
    xsec->SetBinContent(i+1, xsecVal[i]);
    xsec->SetBinError(i+1, errXsec[i]);
  }

  TFile* out = new TFile(argv[2], "RECREATE");
  xsec->Write();
  out->Close();

//  std::string textout = argv[2].replace(argv[2].find(".root"), 5, ".txt")
//  ofstream textfile;
//  textfile.open(textout);
//  textfile << "Mass       xs\n";
//  for (int i = 0; i<100; i++)
//  {
//      textfile << xsec->GetBinLowEdge(i+1) << "    " << xsecVal[i] << "\n";
//  }
//  textfile.close();
  // Done.
  return 0;
}
