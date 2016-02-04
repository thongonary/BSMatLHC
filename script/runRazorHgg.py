#!/usr/bin/env python
##########################################################################################################
# AUTHOR: Maurizio Pierini (CERN)                                                                        #
# DESCRIPTION:                                                                                           #
# Script to generate PYTHIA events and process them through BSMApp                                       #
##########################################################################################################

import os, sys
from optparse import OptionParser
import ROOT as rt

def exec_me(command,dryRun=True):
    print command
    if not dryRun: os.system(command)

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('--pythia-card',dest="pythiaCardTemplate",type="string",default="config/run2.config",
                  help="Template pythia card to use")
    parser.add_option('--slha',dest="slhaTemplate",type="string",default="config/run2.config",
                  help="Template pythia card to use")
    parser.add_option('-d','--dir',dest="outDir",default="output",type="string",
                  help="Output directory to store output")
    parser.add_option('-n','--nevents',dest="nevents", default=10000,type="int",
                  help="number of events")
    parser.add_option('-s','--sqrts',dest="sqrts", default=8000,type="float",
                  help="CM energy")
    parser.add_option('-l','--lumi',dest="lumi", default=19800,type="float",
                  help="integrated luminosity")
    parser.add_option('--mLSP',dest="mLSP", default=1,type="int",
                  help="value to replace LSP in all templates")
    parser.add_option('--mParent',dest="mParent", default=130,type="int",
                  help="value to replace parent mass in all templates")
    parser.add_option('--mOtherParent',dest="mOtherParent", default=130,type="int",
                  help="value to replace other parent mass in all templates")
    parser.add_option('--xsec',dest="xsecMax", default=20,type="float",
                  help="max xsec for priors, etc.")
    parser.add_option('-m','--model',dest="model",type="string",default="TChiwh",
                  help="model name")
    parser.add_option('--dry-run',dest='dryRun',default=False,action='store_true',
                  help="dry run")
    parser.add_option('--no-gen',dest='noGen',default=False,action='store_true',
                  help="no generation (GenPythia) step")
    parser.add_option('--no-sim',dest='noSim',default=False,action='store_true',
                  help="no simulation (Delphes) step")


    (options,args) = parser.parse_args()
    

    # Get the slha name and number of events from the arguments:
    pwd = os.environ['PWD']

    # Set the directories of the different packages:
    susygendir = pwd+'/BSMGen/'
    susyappdir = pwd+'/BSMApp/'
    delphesdir = pwd+'/delphes/'

    outDir = options.outDir
    tmpDir = pwd+'/tmp'
    exec_me('mkdir -p %s'%outDir,options.dryRun)
    exec_me('mkdir -p %s'%tmpDir,options.dryRun)
    if outDir[0]!='/':
        outDir = pwd+'/'+outDir
        
    delphesCard = 'cards/delphes_card_CMS_Cpena.tcl'
    
    massPoint = '%i.%i'%(options.mParent,options.mLSP)
    if options.model=="T21bH":
        massPoint = '%i.%i.%i'%(options.mParent,options.mOtherParent,options.mLSP)
        
    
    #create the SLHA card from template    
    slha = options.slhaTemplate.replace(options.model,'%s.%s'%(options.model,massPoint))
    slhaFile = open(slha, 'w')
    slhaTemplateFile = open(options.slhaTemplate,'r')
    for line in slhaTemplateFile:
        if 'CHARGINO' in line and options.model=="TChiwh":
            newline = line.replace('CHARGINO','%i'%options.mParent)
            slhaFile.write(newline)
        elif 'SBOTTOM' in line and options.model=="T2bH":
            newline = line.replace('SBOTTOM','%i'%options.mParent)
            slhaFile.write(newline)
        elif 'SBOTTOM2' in line and options.model=="T21bH":
            newline = line.replace('SBOTTOM2','%i'%options.mParent)
            slhaFile.write(newline)
        elif 'SBOTTOM1' in line and options.model=="T21bH":
            newline = line.replace('SBOTTOM1','%i'%(options.mOtherParent))
            slhaFile.write(newline)
        elif 'NLSP' in line and (options.model=="T2bH" or options.model=="T21bH"):
            newline = line.replace('NLSP','%i'%(options.mLSP+130))
            slhaFile.write(newline)
        elif 'LSP' in line:
            newline = line.replace('LSP','%i'%options.mLSP)
            slhaFile.write(newline)
        else:
            slhaFile.write(line)
    slhaFile.close()
    slhaTemplateFile.close()

    #create the PYTHIA8 card from template
    pythiaCard = options.pythiaCardTemplate.replace(options.model,'%s.%s'%(options.model,massPoint))
    pythiaCardFile = open(pythiaCard, "w")
    pythiaCardTemplateFile = open(options.pythiaCardTemplate,"r")
    for line in pythiaCardTemplateFile:
        if "SLHA:file = " in line:
            pythiaCardFile.write("SLHA:file = %s\n" %slha.replace('BSMGen/',''))
        elif "Beams:eCM =" in line:
            pythiaCardFile.write("Beams:eCM = %f\n" %options.sqrts)
        elif "Main:numberOfEvents =" in line:
            pythiaCardFile.write("Main:numberOfEvents = %i\n" %options.nevents)
        elif "SLHA:useDecayTable =" in line:
            pythiaCardFile.write("SLHA:useDecayTable = on\n")
        else:
            pythiaCardFile.write(line)
    pythiaCardFile.close()
    pythiaCardTemplateFile.close()

    pythiaCard = pwd+'/'+pythiaCard
    slha = pwd+'/'+slha

    if options.noGen:
        # copy GenPythia output to tmpdir:
        pythiaOut = outDir+'/'+pythiaCard.split('/')[-1].replace(".pythia","")
    else:            
        # Run GenPythia to generate events:
        pythiaOut = tmpDir+'/'+pythiaCard.split('/')[-1].replace(".pythia","")
        command = 'source ../script/ToBuild/setupHepmc.sh; source setup.sh; ./GenPythia %s %s --filter'%(pythiaCard,pythiaOut)
        os.chdir(susygendir)
        exec_me(command,options.dryRun)
        exec_me('cp %s/simplifiedModel.%s.%s_GenTree.root %s'%(tmpDir, options.model, massPoint, outDir),options.dryRun)
        exec_me('cp %s/simplifiedModel.%s.%s.hepmc %s'%(tmpDir, options.model, massPoint, outDir),options.dryRun)
        exec_me('cp %s/simplifiedModel.%s.%s.lhe %s'%(tmpDir, options.model, massPoint, outDir),options.dryRun)

    # Get gen-level filter information
    if options.dryRun:
        filtereff = 1
    else:
        genTreeFile = rt.TFile.Open(pythiaOut+'_GenTree.root')
        infoTree = genTreeFile.Get('infoTree')
        infoTree.GetEntry(0)
        filtereff = infoTree.filtereff
    if options.model=='TChiwh' or options.model=="T21bH":
        filtereff *= 2.28E-03 # for 125/126 GeV from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR3#Higgs_2_gauge_bosons

    if options.noSim:        
        delphesOut = pythiaOut.replace(tmpDir,outDir)+'_delphes.root'
    else:
        # Run Delphes to simulate the CMS detector
        delphesOut = pythiaOut.replace(outDir,tmpDir)+'_delphes.root'
        command = './DelphesHepMC %s %s %s'%(delphesCard,delphesOut,pythiaOut+'.hepmc')
        os.chdir(delphesdir)
        exec_me(command,options.dryRun)
        exec_me('cp %s/simplifiedModel.%s.%s_delphes.root %s'%(tmpDir, options.model, massPoint, outDir),options.dryRun)

    # Run CMSApp to do apply the CMS analyses to the generated sample
    cmsOut = delphesOut.replace('_delphes.root','_cmsapp.root')
    command = "./CMSApp %s --hggrazor --delphes --output=%s --filter=%f --lumi=%f --xsec=%f --sqrts=%f" %(delphesOut,cmsOut,filtereff,options.lumi,options.xsecMax,options.sqrts)    
    os.chdir(susyappdir)
    exec_me(command,options.dryRun)
    exec_me('cp %s/simplifiedModel.%s.%s_cmsapp.root %s'%(tmpDir, options.model, massPoint, outDir),options.dryRun)
    
    exec_me('rm -f %s/simplifiedModel.%s.%s.*'%(tmpDir, options.model, massPoint),options.dryRun)
    exec_me('rm -f %s/simplifiedModel.%s.%s_*'%(tmpDir, options.model, massPoint),options.dryRun)
 
