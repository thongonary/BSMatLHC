#! /usr/bin/env python
import ROOT as rt
import os.path
import sys, glob, re
from array import *
import math

def setstyle():
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(600) #Width of canvas
    rt.gStyle.SetCanvasDefX(0)   #POsition on screen
    rt.gStyle.SetCanvasDefY(0)
    
    # For the Pad:
    rt.gStyle.SetPadBorderMode(0)
    # rt.gStyle.SetPadBorderSize(Width_t size = 1)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
    
    # For the frame:
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    # set the paper & margin sizes
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.09)
    #rt.gStyle.SetPadRightMargin(0.065)
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.065,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.055," ")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.065,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    #rt.gStyle.SetMarkerStyle(20)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    #rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    #rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(1111)
    rt.gStyle.SetStatY(0.85)        
    rt.gStyle.SetStatX(0.92)                
    rt.gStyle.SetStatW(0.15)                
    rt.gStyle.SetStatH(0.15)                
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    
    ncontours = 999
    
    stops = [ 0.00, 0.1, 0.25, 0.65, 1.00 ]
    #stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
    green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
    blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
        
    npoints = len(s)
    #rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)
   
    rt.gStyle.cd()

if __name__ == '__main__':
    setstyle()
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetNumberContours(999)
    fileName = sys.argv[1]
    outdir = sys.argv[2]

    #rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    c  = rt.TCanvas("c","c",1,1,650,376)
    #c  = rt.TCanvas("c","c",1,1,600,500)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    c.Update()
    #c.SetLogy()
    #c.SetLogx()
    emptyEtaPhi = rt.TH2D("emptyEtaPhi","emptyEtaPhi",1,-4,4,1,-rt.TMath.Pi(),rt.TMath.Pi())
    emptyEtaPhi.SetXTitle("#eta")
    emptyEtaPhi.SetYTitle("#phi")
    
    tFile = rt.TFile.Open(fileName)
    tree = tFile.Get("outTree")
    tree.Draw('>>elist','','entrylist')
    elist = rt.gDirectory.Get('elist')
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)

        pTJetHem1 = []
        etaJetHem1 = []
        phiJetHem1 = []

        pTGGHem1 = []
        etaGGHem1 = []
        phiGGHem1 = []
        
        pTJetHem2 = []
        etaJetHem2 = []
        phiJetHem2 = []

        pTGGHem2 = []
        etaGGHem2 = []
        phiGGHem2 = []


        tlatexList = []
        
        for i in range(0,tree.NPFJets):
            if tree.constPFHem1[i]>-1 and tree.constPFHem1[i]<tree.NPFJets:
                pTJetHem1.append(tree.pTPFJet[tree.constPFHem1[i]])
                etaJetHem1.append(tree.etaPFJet[tree.constPFHem1[i]])
                phiJetHem1.append(tree.phiPFJet[tree.constPFHem1[i]])
                tlatex = rt.TLatex(tree.etaPFJet[tree.constPFHem1[i]]-0.6,tree.phiPFJet[tree.constPFHem1[i]]+0.25,"%2.1f"%tree.pTPFJet[tree.constPFHem1[i]])
                tlatex.SetTextSize(0.05)
                tlatex.SetTextFont(42)
                tlatex.SetTextColor(rt.kBlue)
                tlatexList.append(tlatex)
            elif tree.constPFHem1[i]>-1 and tree.constPFHem1[i]==tree.NPFJets:
                pTGGHem1.append(tree.pTGG)
                etaGGHem1.append(tree.etaGG)
                phiGGHem1.append(tree.phiGG)
                tlatex = rt.TLatex(tree.etaGG-0.4,tree.phiGG+0.25,"%2.1f"%tree.pTGG)
                tlatex.SetTextSize(0.05)
                tlatex.SetTextFont(42)
                tlatex.SetTextColor(rt.kBlue)
                tlatexList.append(tlatex)

        for i in range(0,tree.NPFJets):
            if tree.constPFHem2[i]>-1 and tree.constPFHem2[i]<tree.NPFJets:
                pTJetHem2.append(tree.pTPFJet[tree.constPFHem2[i]])
                etaJetHem2.append(tree.etaPFJet[tree.constPFHem2[i]])
                phiJetHem2.append(tree.phiPFJet[tree.constPFHem2[i]])
                tlatex = rt.TLatex(tree.etaPFJet[tree.constPFHem2[i]]-0.4,tree.phiPFJet[tree.constPFHem2[i]]+0.25,"%2.1f"%tree.pTPFJet[tree.constPFHem2[i]])
                tlatex.SetTextSize(0.05)
                tlatex.SetTextFont(42)
                tlatex.SetTextColor(rt.kRed)
                tlatexList.append(tlatex)
            elif tree.constPFHem2[i]>-1 and tree.constPFHem2[i]==tree.NPFJets:
                pTGGHem2.append(tree.pTGG)
                etaGGHem2.append(tree.etaGG)
                phiGGHem2.append(tree.phiGG)
                tlatex = rt.TLatex(tree.etaGG-0.4,tree.phiGG+0.25,"%2.1f"%tree.pTGG)
                tlatex.SetTextSize(0.05)
                tlatex.SetTextFont(42)
                tlatex.SetTextColor(rt.kRed)
                tlatexList.append(tlatex)
                
        if len(pTJetHem1)>0:
            graphJetHem1 = rt.TGraph(len(pTJetHem1))
            graphJetHem1.SetMarkerColor(rt.kBlue)
            graphJetHem1.SetMarkerSize(2)
            graphJetHem1.SetMarkerStyle(8)
            for i in range(0,len(pTJetHem1)):
                graphJetHem1.SetPoint(i,etaJetHem1[i],phiJetHem1[i])

        if len(pTGGHem1)>0:
            graphGGHem1 = rt.TGraph(len(pTGGHem1))
            graphGGHem1.SetMarkerColor(rt.kBlue)
            graphGGHem1.SetMarkerSize(2)
            graphGGHem1.SetMarkerStyle(21)
            graphGGHem1.SetPoint(0,etaGGHem1[0],phiGGHem1[0])

        if len(pTJetHem2)>0:
            graphJetHem2 = rt.TGraph(len(pTJetHem2))
            graphJetHem2.SetMarkerColor(rt.kRed)
            graphJetHem2.SetMarkerSize(2)
            graphJetHem2.SetMarkerStyle(8)
            for i in range(0,len(pTJetHem2)):
                graphJetHem2.SetPoint(i,etaJetHem2[i],phiJetHem2[i])

        if len(pTGGHem2)>0:
            graphGGHem2 = rt.TGraph(len(pTGGHem2))
            graphGGHem2.SetMarkerColor(rt.kRed)
            graphGGHem2.SetMarkerSize(2)
            graphGGHem2.SetMarkerStyle(21)
            graphGGHem2.SetPoint(0,etaGGHem2[0],phiGGHem2[0])

        graphG = rt.TGraph(2)
        graphG.SetMarkerColor(rt.kGreen)
        graphG.SetMarkerSize(2)
        graphG.SetMarkerStyle(22)
        graphG.SetPoint(0,tree.etaPh1,tree.phiPh1)
        tlatex = rt.TLatex(tree.etaPh1+0.1,tree.phiPh1,"%2.1f"%tree.pTPh1)
        tlatex.SetTextSize(0.05)
        tlatex.SetTextFont(42)
        tlatex.SetTextColor(rt.kGreen)
        tlatexList.append(tlatex)
        graphG.SetPoint(1,tree.etaPh2,tree.phiPh2)
        tlatex = rt.TLatex(tree.etaPh2+0.1,tree.phiPh2,"%2.1f"%tree.pTPh2)
        tlatex.SetTextSize(0.05)
        tlatex.SetTextFont(42)
        tlatex.SetTextColor(rt.kGreen)
        tlatexList.append(tlatex)
                
                
            
        emptyEtaPhi.Draw()
        #emptyEtaPhi.SetTitle("%i : %i : %i, M_{R} = %2.1f, R^{2} = %2.3f, m_{#gamma#gamma} = %2.1f"%(tree.run,tree.lumi,tree.evNum,tree.MR,tree.RSQ,tree.mGG))
        emptyEtaPhi.SetTitle("Run %i, M_{R} = %2.1f, R^{2} = %2.3f, m_{#gamma#gamma} = %2.1f"%(tree.run,tree.MR,tree.RSQ,tree.mGG))
        if len(pTJetHem1)>0:
            graphJetHem1.Draw("psame")
        if len(pTGGHem1)>0:
            graphGGHem1.Draw("psame")
        if len(pTJetHem2)>0:
            graphJetHem2.Draw("psame")
        if len(pTGGHem2)>0:
            graphGGHem2.Draw("psame")
        graphG.Draw("psame")
        for tlatex in tlatexList: tlatex.Draw("same")

        c.Print("%s/etaPhi_%i.pdf"%(outdir,tree.run))
        c.Print("%s/etaPhi_%i.C"%(outdir,tree.run))
    
    
