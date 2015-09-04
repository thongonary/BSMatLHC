// main18.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It illustrates how to write an event filter.
// No new functionality is involved - all could be done in the main program
// - but the division of tasks may be more convenient for recurrent cuts.

#include <EventFilter.hh>
#include "Pythia8/Pythia.h"

using namespace Pythia8;

//==========================================================================

// The EventFilter class.

// The constructor takes the following arguments
// select = 1 : keep only final particles.
//        = 2 : keep only final visible particles (i.e. not neutrinos).
//        = 3 : keep only final charged particles.
// etaMax (default = 50) : keep only particles with pseudorapidity
//        |eta| < etaMax.
// pTminCharged (default = 0) : keep a charged particle only if
//        its transverse momentum pT < pTminCharged.
// pTminNeutral (default = 0) : keep a neutral particle only if
//        its transverse momentum pT < pTminNeutral.

// Main methods:
// filter( event) takes an event record as input and analyzes it.
// size() returns the number of particles kept.
// index(i) returns the index in the full event of the i'th kept particle.
// particlePtr(i) returns a pointer to the i'th kept particle.
// particleRef(i) returns a reference to the i'th kept particle.
// list() gives a listing of the kept particles only.


void EventFilter::filter(Event& event) {

  // Reset arrays in preparation for new event.
  keptIndx.resize(0);
  keptPtrs.resize(0);

  // Loop over all particles in the event record.
  for (int i = 0; i < event.size(); ++i) {

    // Skip if particle kind selection criteria not fulfilled.
    if (!event[i].isFinal()) continue;
    if (select == 2 && !event[i].isVisible()) continue;
    bool isCharged = event[i].isCharged();
    if (select == 3 && !isCharged) continue;

    // Skip if too large pseudorapidity.
    if (abs(event[i].eta()) > etaMax) continue;

    // Skip if too small pT.
    if       (isCharged && event[i].pT() < pTminCharged) continue;
    else if (!isCharged && event[i].pT() < pTminNeutral) continue;

    // Add particle to vectors of indices and pointers.
    keptIndx.push_back( i );
    keptPtrs.push_back( &event[i] );

  // End of particle loop. Done.
  }

}

//--------------------------------------------------------------------------

// The list method: downscaled version of Event::list.

void EventFilter::list(ostream& os) {

  // Header.
  os << "\n --------  PYTHIA Event Listing  (filtered)  ------------------"
     << "-----------------------------------------------------------------"
     << "----\n \n    no        id   name            status     mothers  "
     << " daughters     colours      p_x        p_y        p_z         e  "
     << "        m \n";

  // At high energy switch to scientific format for momenta.
  double eSum = 0.;
  for (int iKept = 0; iKept < size(); ++iKept) eSum += keptPtrs[iKept]->e();
  bool useFixed = (eSum < 1e5);

  // Listing of kept particles in event.
  for (int iKept = 0; iKept < size(); ++iKept) {
    int i = keptIndx[iKept];
    Particle& pt = *keptPtrs[iKept];

    // Basic line for a particle, always printed.
    os << setw(6) << i << setw(10) << pt.id() << "   " << left
       << setw(18) << pt.nameWithStatus(18) << right << setw(4)
       << pt.status() << setw(6) << pt.mother1() << setw(6)
       << pt.mother2() << setw(6) << pt.daughter1() << setw(6)
       << pt.daughter2() << setw(6) << pt.col() << setw(6) << pt.acol()
       << ( (useFixed) ? fixed : scientific ) << setprecision(3)
       << setw(11) << pt.px() << setw(11) << pt.py() << setw(11)
       << pt.pz() << setw(11) << pt.e() << setw(11) << pt.m() << "\n";
  }

  // Listing finished.
  os << "\n --------  End PYTHIA Event Listing  ----------------------------"
     << "-------------------------------------------------------------------"
     << endl;
}
