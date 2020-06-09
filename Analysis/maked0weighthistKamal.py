import ROOT as rt
import math as m
import numpy as np

binsele = np.array([0,0.05,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,3.0,3.5,4.0,4.5,5.0,6.0,10.0], dtype=float)
histoeleKamal = rt.TH1D("d0EffEleKamal","",20,binsele) # in cm
binsmuo = np.array([0,0.2,0.4,0.6,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,6.5,7.0,7.5,8.0,10.0], dtype = float)
histomuoKamal = rt.TH1D("d0EffMuoKamal","",16,binsmuo) # in cm

histoeleKamal.SetBinContent(1,0.945)
histoeleKamal.SetBinContent(2,0.851)
histoeleKamal.SetBinContent(3,0.772)
histoeleKamal.SetBinContent(4,0.765)
histoeleKamal.SetBinContent(5,0.746)
histoeleKamal.SetBinContent(6,0.719)
histoeleKamal.SetBinContent(7,0.688)
histoeleKamal.SetBinContent(8,0.679)
histoeleKamal.SetBinContent(9,0.572)
histoeleKamal.SetBinContent(10,0.457)
histoeleKamal.SetBinContent(11,0.373)
histoeleKamal.SetBinContent(12,0.363)
histoeleKamal.SetBinContent(13,0.339)
histoeleKamal.SetBinContent(14,0.327)
histoeleKamal.SetBinContent(15,0.291)
histoeleKamal.SetBinContent(16,0.233)
histoeleKamal.SetBinContent(17,0.184)
histoeleKamal.SetBinContent(18,0.087)
histoeleKamal.SetBinContent(19,0.063)
histoeleKamal.SetBinContent(20,0.051)
''' 0 TO 0.05; 0.945 +0.002,-0.002;
 0.05 TO 0.2; 0.851 +0.003,-0.003;
 0.2 TO 0.4; 0.772 +0.004,-0.004;
 0.4 TO 0.6; 0.765 +0.006,-0.006;
 0.6 TO 0.8; 0.746 +0.008,-0.008;
 0.8 TO 1.0; 0.719 +0.01,-0.01;
 1.0 TO 1.2; 0.688 +0.011,-0.011;
 1.2 TO 1.4; 0.679 +0.013,-0.012;
 1.4 TO 1.6; 0.572 +0.014,-0.014;
 1.6 TO 1.8; 0.457 +0.016,-0.016;
 1.8 TO 2.0; 0.373 +0.017,-0.017;
 2.0 TO 2.25; 0.363 +0.016,-0.016;
 2.25 TO 2.5; 0.339 +0.017,-0.018;
 2.5 TO 3.0; 0.327 +0.013,-0.013;
 3.0 TO 3.5; 0.291 +0.014,-0.015;
 3.5 TO 4.0; 0.233 +0.015,-0.016;
 4.0 TO 4.5; 0.184 +0.021,-0.022;
 4.5 TO 5.0; 0.087 +0.012,-0.014;
 5.0 TO 6.0; 0.063 +0.009,-0.01;
 6.0 TO 10.0; 0.051 +0.006,-0.007;
'''

c1e = rt.TCanvas()
histoeleKamal.Draw()
c1e.SaveAs("kamald0eleweight.pdf")
c1e.Close()

histomuoKamal.SetBinContent(1,0.959)
histomuoKamal.SetBinContent(2,0.957)
histomuoKamal.SetBinContent(3,0.951)
histomuoKamal.SetBinContent(4,0.953)
histomuoKamal.SetBinContent(5,0.949)
histomuoKamal.SetBinContent(6,0.938)
histomuoKamal.SetBinContent(7,0.913)
histomuoKamal.SetBinContent(8,0.921)
histomuoKamal.SetBinContent(9,0.902)
histomuoKamal.SetBinContent(10,0.869)
histomuoKamal.SetBinContent(11,0.867)
histomuoKamal.SetBinContent(12,0.819)
histomuoKamal.SetBinContent(13,0.784)
histomuoKamal.SetBinContent(14,0.796)
histomuoKamal.SetBinContent(15,0.683)
histomuoKamal.SetBinContent(16,0.614)

'''
 0 TO 0.2; 0.959 +0.001,-0.001;
 0.2 TO 0.4; 0.957 +0.002,-0.002;
 0.4 TO 0.6; 0.951 +0.003,-0.003;
 0.6 TO 1.0; 0.953 +0.003,-0.003;
 1.0 TO 1.5; 0.949 +0.004,-0.003;
 1.5 TO 2.0; 0.938 +0.005,-0.004;
 2.0 TO 2.5; 0.913 +0.007,-0.006;
 2.5 TO 3.0; 0.921 +0.007,-0.007;
 3.0 TO 4.0; 0.902 +0.007,-0.006;
 4.0 TO 5.0; 0.869 +0.01,-0.01;
 5.0 TO 6.0; 0.867 +0.013,-0.012;
 6.0 TO 6.5; 0.819 +0.026,-0.023;
 6.5 TO 7.0; 0.784 +0.028,-0.025;
 7.0 TO 7.5; 0.796 +0.031,-0.028;
 7.5 TO 8.0; 0.683 +0.038,-0.036;
 8.0 TO 10.0; 0.614 +0.025,-0.024;
'''

c1 = rt.TCanvas()
histomuoKamal.Draw()
c1.SaveAs("kamald0muoweight.pdf")
c1.Close()

outfile = rt.TFile("Analysis/d0weightKamal.root","RECREATE")
histoeleKamal.Write()
histomuoKamal.Write()
outfile.Close()
