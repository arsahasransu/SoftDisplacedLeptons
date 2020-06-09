import os, sys
import ROOT as rt
import math as m
import numpy as np

file1= rt.TFile("Analysis/HEPData-ins1317640-v1-Table_5.root","READ")
#file1.ls()
histoeleCMS= file1.Get("Table 5/Hist1D_y1") # in cm!!
file2= rt.TFile("Analysis/HEPData-ins1317640-v1-Table_6.root","READ")
#file2.ls()
histomuoCMS= file2.Get("Table 6/Hist1D_y1") # in cm!!

histoeleFreya = rt.TH1D("d0EffEleFreya","",100,0,10) # in cm
histomuoFreya = rt.TH1D("d0EffMuoFreya","",100,0,10) # in cm

for d0bin in range(histoeleCMS.GetNbinsX()):

    d0bin = d0bin+1

    histoeleFreya.SetBinContent(d0bin, histoeleCMS.GetBinContent(d0bin))
    histomuoFreya.SetBinContent(d0bin, histomuoCMS.GetBinContent(d0bin))

transitionele = 1.0 # in cm
transitionmuo = 2.0 # in cm
maxvald0ele = 0.05
maxvald0muo = 0.5
startvald0ele = histoeleCMS.GetBinContent(histoeleCMS.GetXaxis().FindBin(transitionele))
startvald0muo = histomuoCMS.GetBinContent(histomuoCMS.GetXaxis().FindBin(transitionmuo-0.00001))
ricoele = (maxvald0ele-startvald0ele)/(10.0-transitionele)
ricomuo = (maxvald0muo-startvald0muo)/(10.0-transitionmuo)

lowedge = np.zeros(100)
histoeleFreya.GetLowEdge(lowedge)
for d0binlow in lowedge:

    if(d0binlow<transitionele):
        continue

    histoeleFreya.SetBinContent(histoeleFreya.FindBin(d0binlow+0.05), startvald0ele+(ricoele*(d0binlow+0.1-transitionele)))
    print(d0binlow,"\t",startvald0ele+(ricoele*(d0binlow+0.1-transitionele)))

lowedge = np.zeros(100)
histomuoFreya.GetLowEdge(lowedge)
for d0binlow in lowedge:

    if(d0binlow<transitionmuo):
        continue

    histomuoFreya.SetBinContent(histomuoFreya.FindBin(d0binlow+0.05), startvald0muo+(ricomuo*(d0binlow+0.1-transitionmuo)))
    print(d0binlow,"\t",startvald0muo+(ricomuo*(d0binlow+0.1-transitionmuo)))

outfile = rt.TFile("Analysis/d0weightFreya.root","RECREATE")
histoeleFreya.Write()
histomuoFreya.Write()
outfile.Close()
