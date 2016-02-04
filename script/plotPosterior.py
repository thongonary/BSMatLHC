from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import re
import glob

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-m','--model',dest="model", default="TChiwh",type="string",
                  help="signal model name")
    parser.add_option('--mLSP',dest="mLSP", default=1,type="float",
                  help="mass of LSP")
    parser.add_option('--mParent',dest="mParent", default=1,type="float",
                  help="mass of parent")
    parser.add_option('--mOtherParent',dest="mOtherParent", default=130,type="float",
                  help="mass of other parent")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")
    parser.add_option('-i','--indir',dest="inDir",default="./",type="string",
                  help="Input directory")
    parser.add_option('--expected',dest='expected',default=False,action='store_true',
                  help="expected")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal systematic (for bayesian result)")
    
    (options,args) = parser.parse_args()

    massPoint = '%i.%i'%(options.mParent,options.mLSP)
    if options.model=="T21bH":
        massPoint = '%i.%i.%i'%(options.mParent,options.mOtherParent,options.mLSP)
        
    tfile = rt.TFile.Open('%s/simplifiedModel.%s.%s_cmsapp.root'%(options.inDir,options.model,massPoint))
    
    #boxes = ['HighPt','HighRes','Hbb','Zbb','Total']
    boxes = ['HighPt','HighRes','LowRes','Total']
    colors = {'HighPt': rt.kBlue, 'HighRes': rt.kRed, 'LowRes': rt.kViolet, 'Hbb': rt.kGreen, 'Zbb': rt.kMagenta, 'Total': rt.kBlack}
    hists = {}
    histsFill = {}
    xsecUL = {}
    maxBox = boxes[0]
    maxVal = 0

    signalSys = "WithSignalSYS"
    if options.noSignalSys:
        signalSys = ""
    limit = tfile.Get('RazorInclusiveEfficiency')
    limit.GetEntry(0)
    for box in boxes:
        if options.expected:            
            hist = tfile.Get("probVecExp"+box+signalSys)
        else:
            hist = tfile.Get("probVec"+box+signalSys)
        hist.SetLineColor(colors[box])
        hists[box] = hist
        exec('mylimit = float(limit.xsecUL%s%s)'%(box,signalSys))
        xsecUL[box] = mylimit
        if hist.GetBinContent(hist.GetMaximumBin()) > maxVal:
            maxVal = hist.GetBinContent(hist.GetMaximumBin())
            maxBox = box

    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    
    
    c = rt.TCanvas('c','c',500,360)
    #c.SetLogy()
    if options.model=="T2bH":
        hists[maxBox].GetXaxis().SetTitle("Cross Section #sigma(#tilde{b}#tilde{b}) [pb]")
    if options.model=="T21bH":
        hists[maxBox].GetXaxis().SetTitle("Cross Section #sigma(#tilde{b}_{1}#tilde{b}_{2}) [pb]")
    if options.model=="TChiwh":
        hists[maxBox].GetXaxis().SetTitle("Cross Section #sigma(#tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}) [pb]")
    #hists[maxBox].GetXaxis().SetRangeUser(0,10)
    if options.expected:
        hists[maxBox].GetYaxis().SetTitle("Expected Posterior Probability")
    else:
        hists[maxBox].GetYaxis().SetTitle("Posterior Probability")
    hists[maxBox].GetYaxis().SetTitleOffset(1.3)
    hists[maxBox].Draw()
    graphs = []
    #lines = []
    for box in boxes:
        hists[box].Draw("same")
        graph = rt.TGraph(4)
        graph.SetPoint(0,xsecUL[box],0)
        graph.SetPoint(1,hists[box].GetXaxis().GetBinUpEdge(hists[box].FindBin(xsecUL[box])),0)
        graph.SetPoint(2,hists[box].GetXaxis().GetBinUpEdge(hists[box].FindBin(xsecUL[box])),hists[box].GetBinContent(hists[box].FindBin(xsecUL[box])))
        graph.SetPoint(3,xsecUL[box],hists[box].GetBinContent(hists[box].FindBin(xsecUL[box])))
        graph.SetFillColor(colors[box])
        graphs.append(graph)        
        #lines.append(rt.TLine(xsecUL[box],0,xsecUL[box],hists[box].GetBinContent(hists[box].FindBin(xsecUL[box]))))
        histsFill[box] = hists[box].Clone('hist%sFill'%box)
        histsFill[box].GetXaxis().SetRangeUser(xsecUL[box]+hists[box].GetBinWidth(1),hists[box].GetXaxis().GetXmax())
        histsFill[box].SetFillColor(colors[box])
        histsFill[box].Draw("same")
        graph.Draw("fill")
        #lines[-1].SetLineColor(colors[box])
        #lines[-1].SetLineWidth(2)

        

        
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.84,"Pythia+Delphes (8 TeV)")
    if options.model=="TChiwh":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{#chi}^{#pm}_{1}#tilde{#chi}^{0}_{2}, #tilde{#chi}_{1}^{#pm}#rightarrowW^{#pm}#tilde{#chi}^{0}_{1},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.78,"m_{#tilde{#chi}^{#pm}_{1}} = %i GeV, m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mParent,options.mLSP))
    if options.model=="T2bH":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{b}#tilde{b}, #tilde{b}#rightarrowb#tilde{#chi}^{0}_{2},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.78,"m_{#tilde{b}} = %i GeV, m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mParent,options.mLSP))
    if options.model=="T21bH":
        l.DrawLatex(0.52,0.84,"pp #rightarrow #tilde{b}_{1}#tilde{b}_{2}, #tilde{b}_{2}#rightarrowb#tilde{#chi}^{0}_{2},  #tilde{#chi}_{2}^{0}#rightarrowH#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.52,0.78,"m_{#tilde{b}_{2}} = %i GeV, m_{#tilde{b}_{1}} = %i GeV"%(options.mParent,options.mOtherParent))
        l.DrawLatex(0.52,0.72,"m_{#tilde{#chi}^{0}_{2}} = %i GeV, m_{#tilde{#chi}^{0}_{1}} = %i GeV"%(options.mLSP+130, options.mLSP))
    if len(boxes)==5:
        leg = rt.TLegend(0.7,0.4,0.89,0.65)
    elif len(boxes)==4:
        leg = rt.TLegend(0.7,0.45,0.89,0.65)
    elif len(boxes)==3:
        leg = rt.TLegend(0.7,0.5,0.89,0.65)
    elif len(boxes)==2:
        leg = rt.TLegend(0.7,0.55,0.89,0.65)
    elif len(boxes)==1:
        leg = rt.TLegend(0.7,0.6,0.89,0.65)
    
        
    #leg = rt.TLegend(0.7,0.5,0.89,0.7)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)

    for box in boxes:
        leg.AddEntry(hists[box], box,"l")

    leg.Draw()

    massPoint = '%i_%i'%(options.mParent,options.mLSP)
    if options.model=="T21bH":
        massPoint = '%i_%i_%i'%(options.mParent,options.mOtherParent,options.mLSP)

    if options.expected:        
        c.Print("%s/posterior_expected_%s_%s.pdf"%(options.outDir,options.model,massPoint))
        c.Print("%s/posterior_expected_%s_%s.C"%(options.outDir,options.model,massPoint))
    else:
        c.Print("%s/posterior_%s_%s.pdf"%(options.outDir,options.model,massPoint))
        c.Print("%s/posterior_%s_%s.C"%(options.outDir,options.model,massPoint))
