import ROOT as rt
import re
import argparse
import os
from subprocess import call
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
        mass_point = [int(s) for s in re.findall(r'\d+',name)]
        hist_name = "m_{#tilde{q}}="+str(mass_point[0])+". m_{#chi0}="+str(mass_point[1]) + " (GeV)"

        hMR = rt.TH1D("hMR",hist_name+";M_{R}",60,0,4000)
        hMR.SetLineWidth(3)
        hRsq = rt.TH1D("hRsq",hist_name+";R^{2}",60,0,1.)
        hRsq.SetLineWidth(3)
        hMET = rt.TH1D("hMET",hist_name+";MET",60,0,2000.)
        hMET.SetLineWidth(3)
        hDiJetInvMass = rt.TH1D("hDiJetInvMass",hist_name+";DiJetInvMass",60,0,4000.)
        hDiJetInvMass.SetLineWidth(3)
        hScale = rt.TH1D("hScale",hist_name+";Scale",60,0,2000.)
        hScale.SetLineWidth(3)
        hLeadingJetPt = rt.TH1D("hLeadingJetPt",hist_name+";LeadingJetPt",60,0,2000.)
        hLeadingJetPt.SetLineWidth(3)
        hLeadingJetEta = rt.TH1D("hLeadingJetEta",hist_name+";LeadingJetEta",60,-3,3.)
        hLeadingJetEta.SetLineWidth(3)
        hSecondaryJetPt = rt.TH1D("hSecondaryJetPt",hist_name+";SecondaryJetPt",60,0,2000.)
        hSecondaryJetPt.SetLineWidth(3)
        hSecondaryJetEta = rt.TH1D("hSecondaryJetEta",hist_name+";SecondaryJetEta",60,-3,3.)
        hSecondaryJetEta.SetLineWidth(3)
        hnJets = rt.TH1D("hnJets",hist_name+";nJets",10,0,10.)
        hnJets.SetLineWidth(3)
        hnJetsAbove80GeV = rt.TH1D("hnJetsAbove80GeV",hist_name+";nJetsAbove80GeV",10,0,10.)
        hnJetsAbove80GeV.SetLineWidth(3)
        hdeltaRLeadingJetSquark = rt.TH1D("hdeltaRLeadingJetSquark",hist_name+";min(#DeltaR(LeadingJet,#tilde{q}))",60,0,5.)
        hdeltaRLeadingJetSquark.SetLineWidth(3)
        hdeltaRSecondaryJetSquark = rt.TH1D("hdeltaRSecondaryJetSquark",hist_name+";min(#DeltaR(SecondaryJet,#tilde{q}))",60,0,5.)
        hdeltaRSecondaryJetSquark.SetLineWidth(3)

        hRsqMR = rt.TH2D("hRsqMR",hist_name+";M_{R};R^{2}",60,0,4000,60,0,1.)
        hScaleLeadingJetPt = rt.TH2D("hScaleLeadingJetPt",hist_name+";Scale;Leading Jet Pt",60,0,2000,60,0,2000.)
        hScaleMR = rt.TH2D("hScaleMR",hist_name+";Scale;M_{R}",60,0,2000,60,0,4000.)
        hScaleRsq = rt.TH2D("hScaleRsq",hist_name+";Scale;R^{2}",60,0,2000,60,0,1.)
        hDiJetMassMR = rt.TH2D("hDiJetMassMR",hist_name+";M_{R};Leading Dijet Invariant Mass",60,0,4000,60,0,4000)       
        
        hdeltaRLeadingPt = rt.TH2D("hdeltaRLeadingPt",hist_name+";Leading Jet Pt; min(#DeltaR(LeadingJet,#tilde{q}))",60,0,2000,60,0,5.)
        hdeltaRSecondaryPt = rt.TH2D("hdeltaRSecondaryPt",hist_name+";Secondary Jet Pt; min(#DeltaR(SecondaryJet,#tilde{q}))",60,0,2000,60,0,5.)
        
        rFile = rt.TFile.Open(rName,'r')
        rTree = rFile.Get("RazorInclusive")
        nev = rTree.GetEntries()

        for iev in range(nev):
            rTree.GetEntry(iev)
            hMR.Fill(rTree.MR)
            hRsq.Fill(rTree.RSQ)
            hMET.Fill(rTree.MET)
            hScale.Fill(rTree.Scale)
          
            myJetPtList = []
            for i in range(len(rTree.jetPt)): myJetPtList.append(rTree.jetPt[i])
            leadingIndex = myJetPtList.index(sorted(myJetPtList)[-1])
            secondaryIndex = myJetPtList.index(sorted(myJetPtList)[-2])
            leadingJetPt = rTree.jetPt[leadingIndex]
            leadingJetEta = rTree.jetEta[leadingIndex]
            leadingJetPhi = rTree.jetPhi[leadingIndex]
            leadingJetM = rTree.jetM[leadingIndex]
            secondaryJetPt = rTree.jetPt[secondaryIndex]
            secondaryJetEta = rTree.jetEta[secondaryIndex]
            secondaryJetPhi = rTree.jetPhi[secondaryIndex]
            secondaryJetM = rTree.jetM[secondaryIndex]
           
            leadingJet = rt.TLorentzVector()
            leadingJet.SetPtEtaPhiM(leadingJetPt, leadingJetEta, leadingJetPhi, leadingJetM)
            secondaryJet = rt.TLorentzVector()
            secondaryJet.SetPtEtaPhiM(secondaryJetPt, secondaryJetEta, secondaryJetPhi, secondaryJetM)
            TwoLeadingJets = rt.TLorentzVector()
            TwoLeadingJets = leadingJet + secondaryJet

            if (rTree.squarkPt[0] > 0 and rTree.squarkPt[1] > 0):
                squark1 = rt.TLorentzVector()
                squark1.SetPtEtaPhiM(rTree.squarkPt[0], rTree.squarkEta[0], rTree.squarkPhi[0], rTree.squarkM[0])
                squark2 = rt.TLorentzVector()
                squark2.SetPtEtaPhiM(rTree.squarkPt[1], rTree.squarkEta[0], rTree.squarkPhi[0], rTree.squarkM[0])
        
                deltaR_leadingJet_squark = min(leadingJet.DeltaR(squark1), leadingJet.DeltaR(squark2))
                hdeltaRLeadingJetSquark.Fill(deltaR_leadingJet_squark)
                hdeltaRLeadingPt.Fill(leadingJetPt, deltaR_leadingJet_squark)
        
                deltaR_secondaryJet_squark = min(secondaryJet.DeltaR(squark1), secondaryJet.DeltaR(squark2))
                hdeltaRSecondaryJetSquark.Fill(deltaR_secondaryJet_squark)
                hdeltaRSecondaryPt.Fill(secondaryJetPt, deltaR_secondaryJet_squark)
            
            hLeadingJetPt.Fill(leadingJetPt)
            hSecondaryJetPt.Fill(secondaryJetPt)
            hLeadingJetEta.Fill(leadingJetEta)
            hSecondaryJetEta.Fill(secondaryJetEta)
            hnJets.Fill(rTree.numJets)
            hnJetsAbove80GeV.Fill(rTree.numJetsAbove80GeV)
            hDiJetInvMass.Fill(TwoLeadingJets.M())

            hRsqMR.Fill(rTree.MR,rTree.RSQ)
            hScaleLeadingJetPt.Fill(rTree.Scale, leadingJetPt)
            hScaleMR.Fill(rTree.Scale, rTree.MR)
            hScaleRsq.Fill(rTree.Scale, rTree.RSQ)
            hDiJetMassMR.Fill(rTree.MR, TwoLeadingJets.M())

        if not os.path.isdir('/eos/user/q/qnguyen/www/BSMGen/SUSYall/'+name): os.mkdir('/eos/user/q/qnguyen/www/BSMGen/SUSYall/'+name)

        c1 = rt.TCanvas("c1","",650,600)
        rt.gStyle.SetOptLogz()
        hMR.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/MR.png")
        hRsq.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/Rsq.png")
        hDiJetInvMass.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/DiJetInvMass.png")
        hMET.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/MET.png")
        hScale.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/Scale.png")
        hLeadingJetPt.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/LeadingJetPt.png")
        hLeadingJetEta.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/LeadingJetEta.png")
        hSecondaryJetPt.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/SecondaryJetPt.png")
        hSecondaryJetEta.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/SecondaryJetEta.png")
        hnJets.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/nJets.png")
        hnJetsAbove80GeV.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/nJetsAbove80GeV.png")
        hRsqMR.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/RsqMR.png")
        hScaleLeadingJetPt.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/ScaleLeadingJetPt.png")
        hScaleMR.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/ScaleMR.png")
        hScaleRsq.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/ScaleRsq.png")
        hDiJetMassMR.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/DiJetMassMR.png")
        hdeltaRLeadingJetSquark.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/deltaRLeadingJetSquark.png")
        hdeltaRSecondaryJetSquark.Draw()
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/deltaRSecondaryJetSquark.png")
        hdeltaRLeadingPt.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/deltaRLeadingPt.png")
        hdeltaRSecondaryPt.Draw("COLZ")
        c1.SaveAs("/eos/user/q/qnguyen/www/BSMGen/SUSYall/"+name+"/deltaRSecondaryPt.png")
        del hMR, hRsq, hRsqMR, hMET, rFile, rName, name, rTree
        del hLeadingJetPt, hSecondaryJetPt, hLeadingJetEta, hSecondaryJetEta
        del hnJets, hnJetsAbove80GeV, hDiJetInvMass
        del hScale, hScaleLeadingJetPt, hScaleMR, hScaleRsq, hDiJetMassMR, c1
        del hdeltaRLeadingJetSquark, hdeltaRLeadingPt
        del hdeltaRSecondaryJetSquark, hdeltaRSecondaryPt
