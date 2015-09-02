//Enables user info for pseudojets
//Here we save PDG ID


#ifndef ParticleInfo_h
#define ParticleInfo_h

// std includes                                                                                                                                
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <algorithm>
#include <numeric>
#include <list>
#include "CommonTools/combination.hh"

// ParticleInfoApp includes                                                                                                                         
#include <CommonTools/DetectorReco.hh>
#include <CMS/CMSDetectorResponse.hh>
#include <fastjet/SharedPtr.hh>
#include "fastjet/PseudoJet.hh"

using namespace std;

class ParticleInfo: public fastjet::PseudoJet::UserInfoBase {
public:
  ParticleInfo(int id) : pdg_id(id){}

  int pdg_id;
};

#endif
