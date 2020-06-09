import ROOT as rt
import numpy as np

outfile = rt.TFile("Analysis/d0weightNishita.root","RECREATE")

histoeleNishita = rt.TH1D("d0EffEleNishita","",100,0,10) # in cm
histomuoNishita = rt.TH1D("d0EffMuoNishita","",100,0,10) # in cm


weightelelast = 0.0
with open("Analysis/elec-d0.txt") as elefile:
    for lines in elefile:
        line = lines.split()
        xmin = float(line[1])
        weight = float(line[3])
        weightelelast = weight
        histoeleNishita.Fill(xmin+0.05,weight)

weightmuolast = 0.0
with open("Analysis/muon-d0.txt") as muofile:
    for lines in muofile:
        line = lines.split()
        xmin = float(line[1])
        weight = float(line[3])
        weightmuolast = weight
        histomuoNishita.Fill(xmin+0.05,weight)

for d0val in np.arange(2.1,10,0.1):
    d0val = round(d0val,3)
    histoeleNishita.Fill(d0val+0.05, weightelelast)
    histomuoNishita.Fill(d0val+0.05, weightmuolast)
        
histoeleNishita.Write()
histomuoNishita.Write()
outfile.Close()
