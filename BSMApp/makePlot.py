import ROOT as rt
import argparse
import os
rt.gROOT.SetBatch()

parser = argparse.ArgumentParser()
parser.add_argument('infile', help = "List of root files to plot")
args = parser.parse_args()

with open(args.infile) as f:
    lines = f.readlines()
    for i in range(len(lines)):
        if ".root" not in lines[i]: continue
        rName = lines[i].replace('\n','')
        print rName
        name = rName.replace('.root','')
        rt.gStyle.SetOptStat(10)

        hMR = rt.TH1D("hMR",name+";M_{R}",60,0,10000)
        hMR.SetLineWidth(3)
        hRsq = rt.TH1D("hRsq",name+";R^{2}",60,0,2.)
        hRsq.SetLineWidth(3)
        hMET = rt.TH1D("hMET",name+";MET",60,0,10000.)
        hMET.SetLineWidth(3)
        hRsqMR = rt.TH2D("hRsqMR",name+";M_{R};R^{2}",60,0,10000,60,0,2.)
        
        rFile = rt.TFile.Open(rName,'r')
        rTree = rFile.Get("RazorInclusive")
        nev = rTree.GetEntries()

        for iev in range(nev):
            rTree.GetEntry(iev)
            hMR.Fill(rTree.MR)
            hRsq.Fill(rTree.RSQ)
            hMET.Fill(rTree.MET)
            hRsqMR.Fill(rTree.MR,rTree.RSQ)

        canvas = rt.TCanvas(name,name,900,900)
        padRsqMR = rt.TPad("padRsqMR",name,0.0,0.5,0.5,1)
        padMET = rt.TPad("padMET","",0.5,0.5,1,1)
        padRsq = rt.TPad("padRsq","",0,0,0.5,0.5)
        padMR = rt.TPad("padMR","",0.5,0,1,0.5)
        padRsqMR.Draw()
        padRsq.Draw()
        padMET.Draw()
        padMR.Draw()

        padRsqMR.cd()
        padRsqMR.SetLogz()
        hRsqMR.Draw("COLZ")

        padRsq.cd()
        hRsq.Draw()

        padMR.cd()
        hMR.Draw()

        padMET.cd()
        hMET.Draw()

        canvas.SaveAs("/eos/user/q/qnguyen/tchannelPlots/"+name+".png")
        del hMR, hRsq, hRsqMR, hMET, rFile, rName, name, rTree, canvas, padRsqMR, padRsq, padMR, padMET
