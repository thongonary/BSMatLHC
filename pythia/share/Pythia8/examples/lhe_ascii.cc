// main11.cc is a part of the PYTHIA event generator.
// Copyright (C) 2016 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how Les Houches Event File input can be used in Pythia8.
// It uses the ttsample.lhe(.gz) input file, the latter only with 100 events.
#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include "p4.h"

using namespace std;
using namespace Pythia8;

int main(int argc, char** argv) {

  Pythia pythia;

  cout<<"file: "<<string(argv[1])<<", "<<argv[2]<<endl;

  //open file
  ofstream myfile(argv[2]);
  myfile <<"evt, m, deta, dtheta, dphi, e_ratio, ptgg"<<endl;

  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + string(argv[1]));
  pythia.readString("ProcessLevel:all = 0");
  pythia.init();

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; ; ++iEvent) {

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;

      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    myfile<<iEvent<<",";
    
    int nphoton=0;
    P4 p1,p2;
    

    // Sum up final charged multiplicity and fill in histogram.
    for (int i = 0; i < pythia.event.size() && nphoton <2; ++i)
      if(pythia.event[i].isFinal())
	if(pythia.event[i].id()==22)
	{
	  if(nphoton==0)
	    p1.PtEtaPhi(pythia.event[i].pT(),
			pythia.event[i].eta(),
			pythia.event[i].phi()
			);

	  else
	    p2.PtEtaPhi(pythia.event[i].pT(),
			pythia.event[i].eta(),
			pythia.event[i].phi()
			);
	  nphoton++;
	}

    myfile<< (p1+p2).m() << ","
	  << fabs(p1.eta()-p2.eta()) << ","
	  << p1.dtheta(p2) <<","      
	  << p1.dphi(p2) <<","
	  << p1.eratio(p2) <<","
	  << (p1+p2).pt()
	  <<endl;
  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();
  
  // Done.
  return 0;
}
