from optparse import OptionParser
import ROOT as rt
from array import *
import os
import sys
from glob import glob
import csv

            
if __name__ == '__main__':
    
    parser = OptionParser()

    (options,args) = parser.parse_args()

    # format is
    # mrMin,mrMax,rsqMin,rsqMax,nObs,nExp,nExp+1s,nExp-1s,pvalue,nSigma
    # 150,250,0.00,0.05,363,357.6,+9.6,-9.4,0.40,0.3

    rows = []
    for f in args:
        rootFile = rt.TFile.Open(f.split('.csv')[0] + '.root','recreate')
        with open(f,'rb') as csvfile:
            reader = csv.reader(csvfile,delimiter=',')
            for row in reader:
                rows.append(row)

        hEXP = rt.TH1D("hEXP","hEXP",len(rows),0,len(rows))
        hOBS = rt.TH1D("hOBS","hOBS",len(rows),0,len(rows))

        for i in range(0,hEXP.GetNbinsX()):
            hOBS.SetBinContent(i+1,float(rows[i][4]))
            hEXP.SetBinContent(i+1,float(rows[i][5]))
            hEXP.SetBinError(i+1,(float(rows[i][7])-float(rows[i][6]))/2.)
        hEXP.Write()
        hOBS.Write()
        rootFile.Close()
        
    
