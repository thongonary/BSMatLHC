# PYTHIA configuration file.
# Generated on Wed May 17 10:05:23 CEST 2017 with the user supplied options:
# --enable-64bit
# --with-hepmc2=/afs/cern.ch/user/q/qnguyen/BSMatLHC/hepmc/build/

# Install directory prefixes.
PREFIX_BIN=/afs/cern.ch/user/q/qnguyen/BSMatLHC/pythia/bin
PREFIX_INCLUDE=/afs/cern.ch/user/q/qnguyen/BSMatLHC/pythia/include
PREFIX_LIB=/afs/cern.ch/user/q/qnguyen/BSMatLHC/pythia/lib
PREFIX_SHARE=/afs/cern.ch/user/q/qnguyen/BSMatLHC/pythia/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2 -m64 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=false
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

# FASTJET3 configuration.
FASTJET3_USE=false
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/afs/cern.ch/user/q/qnguyen/BSMatLHC/hepmc/build//bin/
HEPMC2_INCLUDE=/afs/cern.ch/user/q/qnguyen/BSMatLHC/hepmc/build//include
HEPMC2_LIB=/afs/cern.ch/user/q/qnguyen/BSMatLHC/hepmc/build//lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=
LHAPDF5_PLUGIN=LHAPDF5.h

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=
LHAPDF6_PLUGIN=LHAPDF5.h

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=
PROMC_LIB=

# ROOT configuration.
ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=
BOOST_LIB=
