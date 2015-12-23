from optparse import OptionParser
import ROOT as rt
import sys
import glob
from math import *
import os
from array import *

def getFileName(mg, mchi, model, directory):
    modelPoint = "%i.%i"%(mg,mchi)
    fileName = "%s/simplifiedModel.%s.%s_cmsapp.root"%(directory,model,modelPoint)
    return fileName


def writeXsecTree(directory, model, mg, mchi, xsecULTotal, xsecULHighPt, xsecULHighRes, xsecULHbb, xsecULZbb):
    outputFileName = "%s/xsecUL_%s_mg_%s_mchi_%s.root" %(directory, model, mg, mchi)
    print "INFO: xsec UL values being written to %s"%outputFileName
    fileOut = rt.TFile.Open(outputFileName, "recreate")
    
    xsecTree = rt.TTree("xsecTree", "xsecTree")
    myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi;"
    for ixsecUL in range(0,5):
        myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL)
    myStructCmd += "}"
    rt.gROOT.ProcessLine(myStructCmd)
    from ROOT import MyStruct

    s = MyStruct()
    xsecTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
    xsecTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')    
    s.mg = mg
    s.mchi = mchi
    
    ixsecUL = 0
    for box in ['Total','HighPt','HighRes','Hbb','Zbb']:
        xsecTree.Branch("xsecUL%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL)),'xsecUL%i/D'%(ixsecUL))
        exec 's.xsecUL%i = xsecUL%s[0]'%(ixsecUL,box)
        ixsecUL += 1

    xsecTree.Fill()

    fileOut.cd()
    xsecTree.Write()
    
    fileOut.Close()
    
    return outputFileName

if __name__ == '__main__':

    
    parser = OptionParser()
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store output")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--xsec-file',dest="refXsecFile",default="./data/gluino13TeV.txt",type="string",
                  help="Input directory")

    (options,args) = parser.parse_args()

    model = options.model
    outDir = options.outDir
    inDir = options.inDir
    

    haddOutputs = []

    gchipairs = []
    for mg in range(300,850,50):
        for mchi in [1]+range(50,mg,50):
            gchipairs.append((mg,mchi))
            
    for mg, mchi in gchipairs:
        
        if not glob.glob(getFileName(mg,mchi,model,inDir)): continue
        print "INFO: opening %s"%(getFileName(mg,mchi,model,inDir))
        tFile = rt.TFile.Open(getFileName(mg,mchi,model,inDir))

        try:
            if tFile.InheritsFrom("TFile") is False:
                continue
        except:
            continue
            
        limit = tFile.Get("RazorInclusiveEfficiency")
        try:
            if limit.InheritsFrom("TTree") is False: 
                tFile.cd()
                tFile.Close()
                continue
        except:
            tFile.cd()
            tFile.Close()
            continue
        if limit.GetEntries() < 1: 
            tFile.cd()
            tFile.Close()
            continue
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limit.GetEntry(entry)
        limits = []
        limit.GetEntry(0)
        limits.append(limit.xsecULTotal)
        limits.append(limit.xsecULHighPt)
        limits.append(limit.xsecULHighRes)
        limits.append(limit.xsecULHbb)
        limits.append(limit.xsecULZbb)
        tFile.cd()
        tFile.Close()
            
        print mg, mchi
        print limits

        haddOutput = writeXsecTree(outDir, model, mg, mchi, [limits[0]],[limits[1]],[limits[2]],[limits[3]],[limits[4]])
        haddOutputs.append(haddOutput)

    os.system("hadd -f %s/xsecUL_Bayesian_%s.root %s"%(outDir,model," ".join(haddOutputs)))
