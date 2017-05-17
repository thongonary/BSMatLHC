#! /usr/bin/env python
import os
import sys
# IMPORTANT: SET ROOTSYS FIRST!!!!
# FOR MACOSX gfortran is needed to compile softsusy, and for MadGraph5
# C++ boost libraries are also needed for LHAPDF-6, can be gotten with MacPorts:
# sudo port install boost

# CURRENT DIRECTORY
BSMDIR =  os.environ['PWD']

# 64bit?
m64 = True
for i in range(len(sys.argv)):
    if sys.argv[i].find("--m32") != -1: m64 = False
    continue

# external code: softsusy-3.4.0
os.system("cp extraCode/softsusy-3.4.0.tar.gz . ; tar -xzf softsusy-3.4.0.tar.gz; rm softsusy-3.4.0.tar.gz")
os.system("mv softsusy-3.4.0 softsusy")
os.chdir("softsusy")
options = "CXXFLAGS=-m32"
if m64: options = "CXXFLAGS=-m64"
os.system("./configure %s; make" %options)
os.chdir(BSMDIR)

# get external code: HepMC 2.06.08
os.system("mkdir hepmc")
os.chdir("hepmc")
os.system("cp ../extraCode/HepMC-2.06.08.tar.gz .; tar -xzf HepMC-2.06.08.tar.gz; rm HepMC-2.06.08.tar.gz")
os.system("mkdir build")
sourcedir = os.environ['PWD']+"/hepmc/HepMC-2.06.08"
builddir = os.environ['PWD']+"/hepmc/build"
os.system("cmake -m64 -DCMAKE_INSTALL_PREFIX=%s %s -Dmomentum:STRING=GEV -Dlength:STRING=MM" %(builddir, sourcedir))
os.system("make; make install")
BSMDIR = os.environ['PWD']
#print "BSMDIR: ", BSMDIR
os.chdir(BSMDIR)

#get external code: boost_1_63_0
os.system("curl -L https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz > boost_1_63_0.tar.gz; tar -zxf boost_1_63_0.tar.gz; rm boost_1_63_0.tar.gz")
os.system("mv boost_1_63_0 boost")
#get external code: LHAPDF 6.1
os.system("cp extraCode/LHAPDF-6.1.5.tar.gz .; tar -xzf LHAPDF-6.1.5.tar.gz; rm LHAPDF-6.1.5.tar.gz")
os.system("cd LHAPDF-6.1.5; ./configure   CC=clang CXX=clang++ --prefix="+BSMDIR+"/lhapdf/ --with-boost="+BSMDIR+"/boost/; make; make install; cd ..; rm -r LHAPDF-6.1.5")

get external code: PYTHIA 8.223
#os.system("cp extraCode/pythia8223.tar .; tar -xvf pythia8223.tar; rm pythia8223.tar")
#os.system("mv pythia8223 pythia")
#os.system("cd pythia;  ./configure --enable-64bit --with-hepmc2="+BSMDIR+"/hepmc/build/ --with-lhapdf6="+BSMDIR+"/lhapdf/; make")
os.system("cd pythia;  ./configure --enable-64bit --with-hepmc2="+BSMDIR+"/hepmc/build/; make")  #try ignoring lhapdf

# get external code: MADGRAPH 5
os.system("cp extraCode/MG5_aMC_v2.3.2.tar.gz .; tar xvzf MG5_aMC_v2.3.2.tar.gz; rm MG5_aMC_v2.3.2.tar.gz")
os.system("mv MG5_aMC_v2_3_2 madgraph; cp -r BSMGen/data/madgraphModels/* madgraph/models/")

# get external code: DELPHES 3
#os.system("cp extraCode/delphes-master.tar.gz .; tar xvzf delphes-master.tar.gz; rm delphes-master.tar.gz; mv delphes-master delphes")
os.system("git clone git@github.com:thongonary/delphes.git")
os.system("cd delphes; make")

#Compile BSMGen
os.system("cd BSMGen; source setup.sh; make")

# fastjet 3.1.3
os.chdir("BSMApp")
os.system("cp ../extraCode/fastjet-3.1.3.tar.gz .")
os.system("tar -xzf fastjet-3.1.3.tar.gz; mkdir fastjet; cd fastjet-3.1.3; ./configure CXXFLAGS=-m64 --prefix=%s/BSMApp/fastjet/; make; make check; make install" %BSMDIR)
os.system("rm -r fastjet-3.1.3.tar.gz fastjet-3.1.3")

#Compile BSMApp
os.system("make")
