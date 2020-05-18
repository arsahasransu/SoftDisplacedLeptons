import os, sys
import ROOT as rt
import math


file1= rt.TFile("HEPData-ins1317640-v1-Table_5.root","READ")
#file1.ls()
histoele= file1.Get("Table 5/Hist1D_y1")
file2= rt.TFile("HEPData-ins1317640-v1-Table_6.root","READ")
#file2.ls()
histomu= file2.Get("Table 6/Hist1D_y1")

SR1low=0.2
SR2low=0.5
SR3low=1.0
SR3high=100.

filelist = ("../../SoftDisplacedLeptons/BP_200_20_02.root",    "../../SoftDisplacedLeptons/BP_200_20_DM.root",
    "../../SoftDisplacedLeptons/BP_200_20_200.root",
    "../../SoftDisplacedLeptons/BP_200_40_20.root",
    "../../SoftDisplacedLeptons/BP_200_20_20.root",
    "../../SoftDisplacedLeptons/BP_324_20_DM.root",
    "../../SoftDisplacedLeptons/BP_200_20_2.root"
)
treelist = {"varTree"}#,"varTree_SR1","varTree_SR2","varTree_SR3"}
for filename in filelist:
#    print(filename)
    for treename in treelist:
#        print(treename)
        ch = rt.TChain(treename,treename)
        ch.Add(filename)
#        ch.Print()
        
        sumnoweightSR1=0.0
        sumweightSR1=0.0
        sumnoweightSR2=0.0
        sumweightSR2=0.0
        sumnoweightSR3=0.0
        sumweightSR3=0.0
        
        for ev in ch: # event loop
            
            weightel = histoele.GetBinContent(min(histoele.GetXaxis().FindBin(0.1*ev.D0El),20))
            weightmu = histomu.GetBinContent(min(histomu.GetXaxis().FindBin(0.1*ev.D0Mu),20))
            weight = weightel * weightmu
#            print(ev.D0El,ev.D0Mu,weight)
            if ev.D0El > SR3high :
                continue
            if ev.D0Mu > SR3high :
                continue
            if ev.D0El > SR3low and ev.D0Mu > SR3low : # SR3!
                sumweightSR3+=weight
                sumnoweightSR3+=1.0
            elif ev.D0El > SR2low and ev.D0Mu > SR2low : # SR2!
                sumweightSR2+=weight
                sumnoweightSR2+=1.0
            elif ev.D0El > SR1low and ev.D0Mu > SR1low : # SR2!
                sumweightSR1+=weight
                sumnoweightSR1+=1.0
        
        print(" model ",filename,
        #" events total:",
#            sumnoweightSR1+sumnoweightSR2+sumnoweightSR3,
            " correction per SR: SR1: ",round(sumweightSR1/sumnoweightSR1,5),
            " SR2: ",round(sumweightSR2/sumnoweightSR2,5),
            " SR3: ",round(sumweightSR3/sumnoweightSR3,5) )
        
