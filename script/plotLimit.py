from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import re
import glob

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="Total",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="TChiwh",type="string",
                  help="signal model name")
    parser.add_option('--mLSP',dest="mLSP", default=1,type="float",
                  help="mass of LSP")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")
    
    (options,args) = parser.parse_args()


    inDir = options.inDir
    box = options.box
    model = options.model
    mLSP = options.mLSP
    
    obsArray = array('d')
    cmsArray = array('d')
    expp1sigmaArray = array('d')
    expm1sigmaArray = array('d')
    expp2sigmaArray = array('d')
    expm2sigmaArray = array('d')
    xsecArray = array('d')
    xsec1sigmaArray = array('d')

    massArray = array('d')
    massArray2 = array('d')
    massArray3 = array('d')
    
    for mNLSP in range(120,215,5):
        for line in open('chipmchi08TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mNLSP))==line.split(',')[0]:
                massArray.append(mNLSP)
                xsecArray.append(float(line.split(',')[1]))
                xsec1sigmaArray.append(float(line.split(',')[2]))
                
    zeroArray = array('d',[0 for mNLSP in massArray])
                
    for mNLSP in range(120,215,5):
        for line in open('xsec.txt','r'):
            line = line.replace('\n','')
            if str(int(mNLSP))==line.split(',')[0]:
                massArray3.append(mNLSP)
                cmsArray.append(float(line.split(',')[1]))
                
           
    for mNLSP in range(120,215,5):
        print '%s/simplifiedModel.%s.%i.%i_cmsapp.root'%(inDir,model,mNLSP,mLSP)
        if not glob.glob('%s/simplifiedModel.%s.%i.%i_cmsapp.root'%(inDir,model,mNLSP,mLSP)):
            print 'no file %s/simplifiedModel.%s.%i.%i_cmsapp.root'%(inDir,model,mNLSP,mLSP)
            continue
        tfile = rt.TFile.Open('%s/simplifiedModel.%s.%i.%i_cmsapp.root'%(inDir,model,mNLSP,mLSP))
        try:
            tfile.Print('v')
        except:
            continue
        limit = tfile.Get('RazorInclusiveEfficiency')
        limit.Draw('>>elist','','entrylist')
        limit.GetEntry(0)
        exec('mylimit = limit.xsecUL%s'%box)
        obsArray.append(mylimit)
        massArray2.append(mNLSP)

    obsGraph = rt.TGraph(len(massArray2),massArray2,obsArray)
    obsGraph.SetLineStyle(2)
    obsGraph.SetLineWidth(2)
    
    cmsGraph = rt.TGraph(len(massArray3),massArray3,cmsArray)
    cmsGraph.SetLineWidth(2)
    
    xsecGraphErr = rt.TGraphErrors(len(massArray),massArray,xsecArray,zeroArray,xsec1sigmaArray)
    xsecGraphErr.SetFillStyle(1001)
    xsecGraphErr.SetLineColor(rt.kOrange)
    xsecGraphErr.SetFillColor(rt.kBlue-7)
        
    xsecGraph = rt.TGraph(len(massArray),massArray,xsecArray)
    xsecGraph.SetMarkerSize(0)
    xsecGraph.SetLineStyle(1)
    xsecGraph.SetLineWidth(1)
    xsecGraph.SetLineColor(rt.kOrange)
    
    h_limit = rt.TMultiGraph()
    c = rt.TCanvas('c','c',500,360)
    c.SetLogy(1)
    h_limit.Add(xsecGraphErr)
    h_limit.Add(obsGraph)
    h_limit.Add(cmsGraph)
    h_limit.Draw("a3")
    h_limit.GetXaxis().SetLimits(123,207)
    h_limit.GetXaxis().SetTitle("m_{#chi^{#pm}_{1}} [GeV]")
    h_limit.GetYaxis().SetTitle("95% C.L. Upper Limit Cross Section [pb]")
    h_limit.SetMaximum(100)
    h_limit.SetMinimum(0.01)
    h_limit.Draw("a3")
    
    obsGraph.Draw("l same")
    cmsGraph.Draw("l same")
    xsecGraph.Draw("l same")
        
    rt.gPad.Update()


    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.84,"Pythia+Delphes (8 TeV)")
    l.SetTextFont(52)
    l.DrawLatex(0.15,0.77,"razor %s"%box)
    l.SetTextFont(42)
    if model=="TChiwh":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}, #tilde{#chi}_{1}^{#pm}#rightarrowW^{#pm}#tilde{#chi}^{0}_{1},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
    leg = rt.TLegend(0.5,0.65,0.85,0.8)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)

        
    if model=="TChiwh":
        leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}) #pm 1 #sigma (theory)","lf")
    leg.AddEntry(obsGraph, "observed (emulation)","l")
    leg.AddEntry(cmsGraph, "observed (CMS)","l")

    leg.Draw()

    c.Print("%s/xsecUL_%s_%s.pdf"%(options.outDir,model,box))
    c.Print("%s/xsecUL_%s_%s.C"%(options.outDir,model,box))


