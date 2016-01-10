import sys
import os 
import ROOT as rt
from runRazorHgg import exec_me

def writeScript(outputName,pwd,model,mParent,mLSP,xsec,nevents):

    #Executable = sleep.sh 
    #Universe = vanilla 
    #Output = sleep.out.$(Process)
    #Log = sleep.log
    #Error = sleep.err
    #getenv = True

    user = os.environ['USER'] # This gives the current user.
    
    submitCards=outputName+".sub"
    submitScript=outputName+".ssh"
    outFile=outputName+".out"
    logFile=outputName+".log"
    errFile=outputName+".err"
    

    # script
    outputFile = open(submitScript,'w')
    outputFile.write('#!/bin/sh\n')
    outputFile.write('hostname -f\n')
    outputFile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    outputFile.write('cd /home/%s/CMSSW_7_1_0/src\n'%user)
    outputFile.write('eval `scramv1 runtime -sh`\n')
    #outputFile.write('mkdir -p /wntmp/%s\n'%(user))
    #outputFile.write('mkdir -p /wntmp/%s/%s\n'%(user,outputName)) 
    #outputFile.write('cd /wntmp/%s/%s\n'%(user,outputName))
    outputFile.write('cd %s\n'%(pwd))
    outputFile.write('time python script/runRazorHgg.py -m %s --mParent %i --mLSP %i --xsec %f -d /mnt/hadoop/store/user/%s/%s/ -n %i --pythia-card BSMGen/data/pythiaCards/SUSY/SimplifiedModels/simplifiedModel.%s.pythia --slha BSMGen/data/pythiaCards/SUSY/SimplifiedModels/simplifiedModel.%s.slha\n'%(model, mParent, mLSP, xsec, user, model,nevents,model,model))
    outputFile.close()

    # cards
    outputFile = open(submitCards,'w')
    outputFile.write('Executable = %s\n'%submitScript)
    outputFile.write('Universe = vanilla\n')
    outputFile.write('Output = %s\n'%outFile)
    outputFile.write('Log = %s\n'%logFile)
    outputFile.write('Error = %s\n'%errFile)
    outputFile.write('getenv = True\n')
    outputFile.write('Queue')
    outputFile.close()
    

if __name__ == '__main__':
    pwd = os.environ['PWD']

    model = 'T2bH'
    xsec = 8
    nevents = 10
    
    for mLSP in [100]:
       for mParent in [800]: 
            outputName = "%s_%i_%i"%(model,mParent,mLSP)
            writeScript(outputName,pwd,model,mLSP,mParent,xsec,nevents)
            exec_me('condor_submit %s.sub'%outputName,True)
