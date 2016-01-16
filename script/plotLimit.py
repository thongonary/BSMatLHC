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
    parser.add_option('--mOtherParent',dest="mOtherParent", default=130,type="float",
                  help="mass of other parent")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")    
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal systematic (for bayesian result)")
    parser.add_option('--freq',dest="freq",default=False,action='store_true',
                  help="Use frequentist (asymptotic CLs) result")
    
    (options,args) = parser.parse_args()


    inDir = options.inDir
    box = options.box
    model = options.model
    mLSP = options.mLSP
    
    obsArray = array('d')
    expArray = array('d')
    cmsObsArray = array('d')
    cmsExpArray = array('d')
    expp1sigmaArray = array('d')
    expm1sigmaArray = array('d')
    expp2sigmaArray = array('d')
    expm2sigmaArray = array('d')
    xsecArray = array('d')
    xsec1sigmaArray = array('d')

    massArray = array('d')
    massArray2 = array('d')
    massArray3 = array('d')

    if model=="TChiwh":
        mParentRange = range(120,215,5)
    if model=="T2bH" or model=="T21bH":
        mParentRange = range(250,805,5)
        
    if model=="TChiwh":
        for mParent in mParentRange:
            for line in open('BSMApp/data/chipmchi08TeV.txt','r'):
                line = line.replace('\n','')
                if str(int(mParent))==line.split(',')[0]:
                    massArray.append(mParent)
                    xsecArray.append(float(line.split(',')[1]))
                    xsec1sigmaArray.append(float(line.split(',')[2]))
                    
        for mParent in mParentRange:
            for line in open('BSMApp/data/SUS14017xsec.txt','r'):
                line = line.replace('\n','')
                if str(int(mParent))==line.split(',')[0]:
                    massArray3.append(mParent)
                    cmsExpArray.append(float(line.split(',')[1]))
                    cmsObsArray.append(float(line.split(',')[2]))
                
    if model=="T2bH":
        rootFile = rt.TFile.Open('BSMApp/data/stop.root','r')
        h_xsec = rootFile.Get("stop")
        N_xsec = h_xsec.GetNbinsX()
        for i in range(0,N_xsec):
            if h_xsec.GetBinContent(i+1) <= 0.0: continue
            massArray.append(h_xsec.GetXaxis().GetBinCenter(i+1))
            xsecArray.append(h_xsec.GetBinContent(i+1))
            xsec1sigmaArray.append(h_xsec.GetBinError(i+1))
            
    if model=="T21bH":
        for mParent in mParentRange:
            for line in open('BSMApp/data/sb1sb28TeV.txt','r'):
                line = line.replace('\n','')
                if str(int(mParent))==line.split(',')[0]:
                    massArray.append(mParent)
                    xsecArray.append(float(line.split(',')[2]))
                    xsec1sigmaArray.append(float(line.split(',')[3]))
                
                
           
    for mParent in mParentRange:
        
        massPoint = '%i.%i'%(mParent,options.mLSP)
        if model=="T21bH":
            massPoint = '%i.%i.%i'%(mParent,options.mOtherParent,options.mLSP)
    
        if options.freq:
            if not glob.glob('%s/higgsCombine%s_%i_%i_%s.Asymptotic.mH120.root'%(inDir,model,mParent,mLSP,box)):
                continue
            tfile = rt.TFile.Open('%s/higgsCombine%s_%i_%i_%s.Asymptotic.mH120.root'%(inDir,model,mParent,mLSP,box))
            try:
                tfile.Print('v')
            except:
                continue
            limit = tfile.Get('limit')
            limit.GetEntry(5)
            obslimit = limit.limit
            obsArray.append(obslimit)
            limit.GetEntry(2)
            explimit = limit.limit
            expArray.append(explimit)
            massArray2.append(mParent)            
        else:
            if not glob.glob('%s/simplifiedModel.%s.%s_cmsapp.root'%(inDir,model,massPoint)):
                continue
            tfile = rt.TFile.Open('%s/simplifiedModel.%s.%s_cmsapp.root'%(inDir,model,massPoint))
            try:
                tfile.Print('v')
            except:
                continue
            limit = tfile.Get('RazorInclusiveEfficiency')
            limit.GetEntry(0)
            signalSys = 'WithSignalSYS'
            if options.noSignalSys:
                signalSys = ''
            exec('obslimit = limit.xsecUL%s%s'%(box,signalSys))
            obsArray.append(obslimit)
            exec('explimit = limit.xsecULExp%s%s'%(box,signalSys))
            expArray.append(explimit)
            massArray2.append(mParent)

    expGraph = rt.TGraph(len(massArray2),massArray2,expArray)
    expGraph.SetLineStyle(2)
    expGraph.SetLineWidth(2)
    
    obsGraph = rt.TGraph(len(massArray2),massArray2,obsArray)
    obsGraph.SetLineWidth(2)
    
    if model=='TChiwh':
        cmsObsGraph = rt.TGraph(len(massArray3),massArray3,cmsObsArray)
        cmsObsGraph.SetLineWidth(2)
        cmsObsGraph.SetLineColor(rt.kRed+1)
        
        cmsExpGraph = rt.TGraph(len(massArray3),massArray3,cmsExpArray)
        cmsExpGraph.SetLineWidth(2)
        cmsExpGraph.SetLineStyle(2)
        cmsExpGraph.SetLineColor(rt.kRed+1)    
    
    zeroArray = array('d',[0 for mParent in massArray])
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
    c = rt.TCanvas('c','c',500,364)
    c.SetLogy(1)
    h_limit.Add(xsecGraphErr)
    h_limit.Add(obsGraph)
    if model=="TChiwh":
        h_limit.Add(cmsObsGraph)
        h_limit.Add(cmsExpGraph)
    h_limit.Draw("a3")   
    if model=="TChiwh":
        h_limit.GetXaxis().SetTitle("m_{#tilde{#chi}^{#pm}_{1}} [GeV]")
        #h_limit.GetXaxis().SetLimits(123,207)
        h_limit.GetXaxis().SetLimits(130,200) 
    if model=="T2bH":
        h_limit.GetXaxis().SetTitle("m_{#tilde{b}} [GeV]")
        h_limit.GetXaxis().SetLimits(250,800) 
    if model=="T21bH":
        h_limit.GetXaxis().SetTitle("m_{#tilde{b}_{2}} [GeV]")
        h_limit.GetXaxis().SetLimits(250,800) 
    h_limit.GetYaxis().SetTitle("95% C.L. Upper Limit Cross Section [pb]")
    h_limit.SetMaximum(1000)
    h_limit.SetMinimum(0.1)
    h_limit.Draw("a3")
    
    obsGraph.Draw("l same")
    expGraph.Draw("l same")
    if model=="TChiwh":
        cmsObsGraph.Draw("l same")
        cmsExpGraph.Draw("l same")
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
        l.DrawLatex(0.52,0.78,"m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mLSP))
        leg = rt.TLegend(0.5,0.43,0.85,0.68)
        h_limit.SetMinimum(0.5)
        h_limit.SetMaximum(100)
    if model=="T2bH":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{b}#tilde{b}, #tilde{b}#rightarrowb#tilde{#chi}^{0}_{2},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.78,"m_{#tilde{#chi}^{0}_{2}} = %i GeV, m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mLSP+130,options.mLSP))
        leg = rt.TLegend(0.5,0.53,0.85,0.68)
        h_limit.SetMinimum(0.005)
        #h_limit.SetMaximum(100)
    if model=="T21bH":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{b}_{1}#tilde{b}_{2}, #tilde{b}_{2}#rightarrowb#tilde{#chi}^{0}_{2},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.78,"m_{#tilde{b}_{1}} = %i GeV"%(options.mOtherParent))
        l.DrawLatex(0.52,0.72,"m_{#tilde{#chi}^{0}_{2}} = %i GeV, m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mLSP+130,options.mLSP))
        leg = rt.TLegend(0.5,0.53,0.85,0.68)
        h_limit.SetMinimum(0.005)
        #h_limit.SetMaximum(100)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)

        
    if model=="TChiwh":
        leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}) #pm 1 #sigma (theory)","lf")
    if model=="T2bH":
        leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{b}#tilde{b}) #pm 1 #sigma (theory)","lf")
    if model=="T21bH":
        leg.AddEntry(xsecGraphErr, "#sigma_{LO} (#tilde{b}_{1}#tilde{b}_{2}) #pm 1 #sigma (theory)","lf")
    leg.AddEntry(obsGraph, "observed (emulation)","l")
    leg.AddEntry(expGraph, "expected (emulation)","l")
    
    if model=="TChiwh":
        leg.AddEntry(cmsObsGraph, "observed (CMS)","l")
        leg.AddEntry(cmsExpGraph, "expected (CMS)","l")

    leg.Draw()

    massPoint = '%i'%(mLSP)
    if model=="T21bH":
        massPoint = '%i_%i'%(options.mOtherParent,mLSP)
        
    c.Print("%s/xsecUL_%s_%s_%s.pdf"%(options.outDir,model,massPoint,box))
    c.Print("%s/xsecUL_%s_%s_%s.C"%(options.outDir,model,massPoint,box))

    print 'masses  ', massArray2
    print 'expArray', expArray
    print 'obsArray', obsArray

    
    if model!="TChiwh": sys.exit()
    
    expRatio = []
    commMasses = [] 
    obsRatio = []
    for i in range(0,len(massArray3)):       
        for j in range(0,len(massArray2)):        
            if massArray3[i]==massArray2[j]:
                commMasses.append(massArray2[j])
                expRatio.append(expArray[j]/cmsExpArray[i])
                obsRatio.append(obsArray[j]/cmsObsArray[i])
    print 'masses  ', commMasses
    print 'expRatio', expRatio
    print 'obsRatio', obsRatio
    print 'expArray', expArray
    print 'obsArray', obsArray
