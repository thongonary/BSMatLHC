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

def makeHistos(fileName):
    
    #rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    c  = rt.TCanvas("c","c",1,1,650,376)
    #c  = rt.TCanvas("c","c",1,1,600,500)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    #c.SetLogy()
    #c.SetLogx()
    c.Update()
    tleg = rt.TLegend(0.68, 0.75, 0.9, 0.9)
    
    tfile = rt.TFile(fileName)
    tree = tfile.Get('RazorInclusive')

    histoData = {('MR_Default',''):
                 ['MR','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,4000,"M_{R} [GeV]","Prob.",""],
              ('RSQ_Default',''):
              ['RSQ','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,1,"R^{2}","Prob.",""],
              ('RSQ_Default','SizesList[4]<4'):
              ['RSQlt4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, <4 jets',100,0,1,"R^{2}","Prob.",""],
              ('RSQ_Default','SizesList[4]>=4'):
              ['RSQgeq4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, #geq4 jets',100,0,1,"R^{2}","Prob.",""],
              ('pTCM_Default',''):
              ['pTCM','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,300,"p_{T}^{CM} [GeV]","Prob.",""],
              ('sqrtsR_Default/2.',''):
              ['sqrtsR','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,4000,"#sqrt{#hat{s}_{R}}/2 [GeV]","Prob.",""],
              ('sqrtsR_Default/2./MR_Default',''):
              ['sqrtsRoverMR','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0.9,1.1,"#sqrt{#hat{s}_{R}}/2/M_{R}","Prob.",""],
              ('1./gammaRp1Ana_Default',''):
              ['gammaRp1Ana','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,1,"1/#gamma_{R+1} (Approx.)","Prob.",""],
              ('gammaRp1_Default/gammaRp1Ana_Default',''):
              ['gammaRp1AnaovergammaRp1','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0.75,1.25,"1/#gamma_{R+1} (Approx.) / 1/#gamma_{R+1}","Prob.",""],
              ('1./gammaRp1_Default',''):
              ['gammaRp1','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,1,"1/#gamma_{R+1}","Prob.",""],
              ('1./gammaRp1_Default','SizesList[4]<4'):
              ['gammaRp1lt4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, <4 jets',100,0,1,"1/#gamma_{R+1}","Prob.",""],
              ('1./gammaRp1_Default','SizesList[4]>=4'):
              ['gammaRp1geq4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, #geq4 jets',100,0,1,"1/#gamma_{R+1}","Prob.",""],
              ('1./gammaRp1_Default/RSQ_Default',''):
              ['gammaRp1overRSQ','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,5,"1/#gamma_{R+1}/R^{2}","Prob.",""],
              ('RSQ_Default:MR_Default',''):
              ['MRRSQ','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,4000,100,0,1,"M_{R} [GeV]","R^{2}","colz"],
              ('RSQ_Default:MR_Default','SizesList[4]<4'):
              ['MRRSQlt4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, <4 jets',100,0,4000,100,0,1,"M_{R} [GeV]","R^{2}","colz"],
              ('RSQ_Default:MR_Default','SizesList[4]>=4'):
              ['MRRSQgeq4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, #geq4 jets',100,0,4000,100,0,1,"M_{R} [GeV]","R^{2}","colz"],
              ('1./gammaRp1_Default:sqrtsR_Default/2.',''):
              ['sqrtsRgammaRp1','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,4000,100,0,1,"#sqrt{#hat{s}_{R}}/2 [GeV]","1/#gamma_{R+1}","colz"],
              ('1./gammaRp1_Default:sqrtsR_Default/2.','SizesList[4]<4'):
              ['sqrtsRgammaRp1lt4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, <4 jets',100,0,4000,100,0,1,"#sqrt{#hat{s}_{R}}/2","1/#gamma_{R+1}","colz"],
              ('1./gammaRp1_Default:sqrtsR_Default/2.','SizesList[4]>=4'):
              ['sqrtsRgammaRp1geq4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, #geq4 jets',100,0,4000,100,0,1,"#sqrt{#hat{s}_{R}}/2","1/#gamma_{R+1}","colz"],
              ('sqrtsR_Default/2.:MR_Default',''):
              ['MRsqrtsR','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,4000,100,0,4000,"M_{R}","#sqrt{#hat{s}_{R}}/2 [GeV]","colz"],
              ('1./gammaRp1_Default:RSQ_Default',''):
              ['RSQgammaRp1','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,1,100,0,1,"R^{2}","1/#gamma_{R+1}","colz"],
              ('1./gammaRp1_Default:RSQ_Default','SizesList[4]<4'):
              ['RSQgammaRp1lt4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, <4 jets',100,0,1,100,0,1,"R^{2}","1/#gamma_{R+1}","colz"],
              ('1./gammaRp1_Default:RSQ_Default','SizesList[4]>=4'):
              ['RSQgammaRp1geq4jets','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV, #geq4 jets',100,0,1,100,0,1,"R^{2}","1/#gamma_{R+1}","colz"],
              ('1./gammaRp1Ana_Default:1./gammaRp1_Default',''):
              ['gammaRp1gammaRp1Ana','pp#rightarrow#tilde{b}#tilde{b}, m(#tilde{b}) = 900 GeV, m(#tilde{#chi}^{0}) = 100 GeV',100,0,1,100,0,1,"1/#gamma_{R+1}","1/#gamma_{R+1} (Approx.)","colz"]
              }
    histo = {}
    for (var,wgt), histoDatum in histoData.iteritems():
        if len(var.split(":"))==1:
            histo[var,wgt] = rt.TH1D(histoDatum[0],histoDatum[1],histoDatum[2],histoDatum[3],histoDatum[4])
            histo[var,wgt].SetXTitle(histoDatum[5])
            histo[var,wgt].SetYTitle(histoDatum[6])
            tree.Project(histo[var,wgt].GetName(),var,wgt)
            histo[var,wgt].Scale(1./histo[var,wgt].Integral())
            histo[var,wgt].Draw(histoDatum[7])
        elif len(var.split(":"))==2:
            histo[var,wgt] = rt.TH2D(histoDatum[0],histoDatum[1],histoDatum[2],histoDatum[3],histoDatum[4],histoDatum[5],histoDatum[6],histoDatum[7])
            histo[var,wgt].SetXTitle(histoDatum[8])
            histo[var,wgt].SetYTitle(histoDatum[9])
            tree.Project(histo[var,wgt].GetName(),var,wgt)
            histo[var,wgt].Scale(1./histo[var,wgt].Integral())
            histo[var,wgt].Draw(histoDatum[10])
        

        c.Print(histo[var,wgt].GetName()+".pdf")
        c.Print(histo[var,wgt].GetName()+".C")
        
    outputFile = rt.TFile.Open("%s_plots.root"%fileName.split('.root')[0],"recreate")
    for var,histogram in histo.iteritems():
        histogram.Write()
def printPlots(plotDict, c, tleg, colors, varName, pdfName):
    c.Clear()
    tleg.Clear()
    c.cd()
    if varName=="R^{2}":
        c.SetLogy()
    c.Print(pdfName)


    
if __name__ == '__main__':
    setstyle()
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetNumberContours(999)
    fileName = sys.argv[1]
    makeHistos(fileName)
