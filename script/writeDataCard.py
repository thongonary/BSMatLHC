from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import re
import glob



def writeDataCard(box,model,txtfileName,obs_yields,bkg_yields,bkg_errors,sig_yields):
        nChan = len(bkg_yields)
        nBkgd = 1
        divider = "------------------------------------------------------------\n"        
        binString = "bin"
        obsString = "observation"
        for i in range(0,nChan):
            binString +="\tb%i"%i
            obsString +="\t%.3f"%obs_yields[i]
        binString+="\n"; obsString +="\n"
        
        datacard = "imax %i number of channels\n"%nChan + \
                   "jmax %i number of backgrounds\n"%nBkgd + \
                   "kmax * number of nuisance parameters\n" + \
                   divider + \
                   binString + \
                   obsString + \
                   divider
        binString = "bin"
        processString = "process"
        processNumberString = "process"
        rateString = "rate"
        proc = ['sig','bg']
        yields = [sig_yields, bkg_yields]
        for i in range(0,nChan):
            for j in range(0,len(proc)):                
                binString +="\tb%i"%i
                processString += "\t%s"%proc[j]
                processNumberString += "\t%i"%j
                rateString += "\t%.3f" %yields[j][i]
        binString+="\n"; processString+="\n"; processNumberString+="\n"; rateString +="\n"
        datacard+=binString+processString+processNumberString+rateString+divider
        # now nuisances
        for i in range(0,nChan):
            bgLnString='bgb%i\tlnN'%i
            for j in range(0,nChan):
                bgLnString+='\t%.3f\t%.3f'%(1.0,1.0+(i==j)*bkg_errors[i])
            bgLnString +="\n"
            datacard+=bgLnString
        txtfile = open(txtfileName,"w")
        txtfile.write(datacard)
        txtfile.close()


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="Total",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="TChiwh",type="string",
                  help="signal model name")
    parser.add_option('--mLSP',dest="mLSP", default=1,type="float",
                  help="mass of LSP")
    parser.add_option('--mParent',dest="mParent", default=130,type="float",
                  help="mass of Parent")
    parser.add_option('-l','--lumi',dest="lumi", default=19800,type="float",
                  help="luminosity in pb^-1")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")
    
    
    (options,args) = parser.parse_args()

    xsec = 1.
     
    if not glob.glob('%s/simplifiedModel.%s.%i.%i_cmsapp.root'%(options.inDir,options.model,options.mParent,options.mLSP)):
        sys.exit()
    tfile = rt.TFile.Open('%s/simplifiedModel.%s.%i.%i_cmsapp.root'%(options.inDir,options.model,options.mParent,options.mLSP))
    tfile.Print()
    datfile = rt.TFile.Open('BSMApp/data/ExpectedObserved_RazorHgg_%s_Summer2015.root'%options.box)
    obsHisto = datfile.Get("hOBS")
    expHisto = datfile.Get("hEXP")

    obs_yields = []
    for i in range(1,obsHisto.GetNbinsX()+1):
        obs_yields.append(obsHisto.GetBinContent(i))
        
    bkg_yields = []
    bkg_errors = []
    for i in range(1,expHisto.GetNbinsX()+1):
        bkg_yields.append(expHisto.GetBinContent(i))
        bkg_errors.append(expHisto.GetBinError(i)/expHisto.GetBinContent(i))
    
    sig_yields = []
    tree = tfile.Get('RazorInclusiveEfficiency')
    tree.GetEntry(0)
    exec 'eff = tree.eff%s'%options.box
    pdfHisto = tfile.Get('pdf%s'%options.box)
    for i in range(1,pdfHisto.GetNbinsX()+1):
        sig_yields.append(pdfHisto.GetBinContent(i)*eff*options.lumi*xsec)

    txtfileName = '%s/datacard_%s_%i_%i_%s.txt'%(options.outDir,options.model,options.mParent,options.mLSP,options.box)
    writeDataCard(options.box,options.model,txtfileName,obs_yields,bkg_yields,bkg_errors,sig_yields)
