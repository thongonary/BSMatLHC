// main18.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to write an event filter.
// No new functionality is involved - all could be done in the main program
// - but the division of tasks may be more convenient for recurrent cuts.

#ifndef EventFilter_h
#define EventFilter_h

#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================


class EventFilter {

public:

  // Constructor sets properties of filter.
  EventFilter( int pdgIdIn = 22, int pdgMothIdIn = 25, double etaMaxIn = -1,
    double pTminIn = 0.)
    : pdgId(pdgIdIn), pdgMothId(pdgMothIdIn), etaMax(etaMaxIn), pTmin(pTminIn) {}

  // Analysis of each new event to find acceptable particles.
  void filter(Event& event);

  // Return size of array, and index of a particle.
  int size()       const {return keptPtrs.size();}
  int index(int i) const {return keptIndx[i];}

  // Return pointer or reference to a particle.
  Particle* particlePtr(int i) {return  keptPtrs[i];}
  Particle& particleRef(int i) {return *keptPtrs[i];}

  // List kept particles only.
  void list(ostream& os = cout);

private:

  // Filter properties, set by constructor.
  int    pdgId, pdgMothId;
  double etaMax, pTmin;

  // Kept particle indices and pointers, referring to original event.
  vector<int>       keptIndx;
  vector<Particle*> keptPtrs;

};

#endif // EventFilter_h

