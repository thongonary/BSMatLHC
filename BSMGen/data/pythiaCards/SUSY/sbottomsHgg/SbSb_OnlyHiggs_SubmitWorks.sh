#!/bin/sh
export SCRAM_ARCH=slc5_amd64_gcc472
source /cvmfs/cms.cern.ch/cmsset_default.sh
mkdir SbSb_OnlyHiggs_Prod_800_650
cd SbSb_OnlyHiggs_Prod_800_650

cmsrel CMSSW_6_2_0;cd CMSSW_6_2_0/src;
cp -r /home/cmorgoth/git/BSMatLHC .
cd BSMatLHC/
echo "BUILDING HEPMC"
export PATH=($PATH):/usr/bin/:/home/cmorgoth/cmake/bin
rm -rf hepmc;mkdir hepmc;cd hepmc;cp ../extraCode/HepMC-2.06.08.tar.gz .; tar -xzf HepMC-2.06.08.tar.gz; rm HepMC-2.06.08.tar.gz
cmsenv
mkdir build;cmake -m64 -DCMAKE_INSTALL_PREFIX=$PWD/build/ $PWD/HepMC-2.06.08/ -Dmomentum:STRING=GEV -Dlength:STRING=MM;
make;make install
export LD_LIBRARY_PATH=/share/apps/root_v5.34.19//lib/root:/opt/gridengine/lib/lx26-amd64:/cvmfs/cms.cern.ch/slc5_amd64_gcc481/external/gcc/4.8.1/lib:$PWD/lib/
cd ../BSMGen/
echo "BUILDING BSMAPP"
cmsenv;source /share/apps/root_v5.34.19/bin/thisroot.sh
make clean;make;
cp -r /home/cmorgoth/git/BSMatLHC/BSMGen/data/pythiaCards/SUSY/sbottomsHgg/SbSb_OnlyHiggs_Prod .
source setup.sh
./GenPythia SbSb_OnlyHiggs_Prod/SbSb_OnlyHiggs_800_650.pythia SbSb_OnlyHiggs_Prod_800_650 
cp SbSb_OnlyHiggs_Prod_800_650_GenTree.root /mnt/hadoop/store/user/cmorgoth/SbSb_OnlyHiggs_Prod/. 
cp SbSb_OnlyHiggs_Prod_800_650.lhe /mnt/hadoop/store/user/cmorgoth/SbSb_OnlyHiggs_Prod/.
cp SbSb_OnlyHiggs_Prod_800_650.hepmc /mnt/hadoop/store/user/cmorgoth/SbSb_OnlyHiggs_Prod/.
echo FINISHING!!!

