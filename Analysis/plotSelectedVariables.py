import os, sys
import ROOT as rt
import math

# Read the variables to chain
sigTree=rt.TChain("varTree","varTree")
bgTree=rt.TChain("varTree","varTree")
sigTree.Add("signal_SR3.root")
bgTree.Add("background_SR3.root")
print("Now filled trees:",sigTree.GetEntries(),bgTree.GetEntries())

# Create a canvas and list of variables
canv = rt.TCanvas("canv","canv")
varList = ["YDelpObj", "dRLL", "dPhiLepMET",
           "Sphericity", "Spherocity", "YUserObj",
           "dPhiLepMETSelObj", "alphaT", "HtDiffHtLepJet",
           "Eta_El", "Eta_Mu"]

# Loop on variables
for varName in varList:
    print("Plotting for variable: ", varName)

    sigTree.Draw(varName+">>sigHisto")
    sigHisto = rt.gDirectory.Get("sigHisto")
    minX = sigHisto.GetXaxis().GetXmin()
    maxX = sigHisto.GetXaxis().GetXmax()
    if minX > sigHisto.GetXaxis().GetXmin() :
        minX = sigHisto.GetXaxis().GetXmin()
    if maxX < sigHisto.GetXaxis().GetXmax() :
        maxX = sigHisto.GetXaxis().GetXmax()
    nBins = sigHisto.GetXaxis().GetNbins()
    sigHisto.Delete()

    histSig = rt.TH1D("histSig"+varName, "", nBins, minX, maxX)
    histSig.SetXTitle(varName)
    histBkg = rt.TH1D("histBkg"+varName, "", nBins, minX, maxX)
    histBkg.SetXTitle(varName)
    drawCommand = varName + ">>"+histSig.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
    sigTree.Draw(drawCommand)
    drawCommand = varName + ">>"+histBkg.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
    bgTree.Draw(drawCommand)
    histSig = rt.gDirectory.Get("histSig"+varName)
    histBkg = rt.gDirectory.Get("histBkg"+varName)
    
    # Normalisation to Luminosity
    # N_bkg = 1937.09 # CR2
    # N_bkg = 3646.28 # SR1
    # N_bkg = 569.73 # SR2
    N_bkg = 22.79 # SR3
    # N_sig = 31.2*271/(5*pow(10,6)) # CR2
    # N_sig = 31.2*28976/(5*pow(10,6)) # SR1
    # N_sig = 31.2*19077/(5*pow(10,6)) # SR2
    N_sig = 31.2*10860/(5*pow(10,6)) # SR3

    histSig.Scale(N_sig/histSig.Integral())
    histBkg.Scale(N_bkg/histBkg.Integral())
    #histSig.Scale(1./histSig.GetSum())
    #histBkg.Scale(1./histBkg.GetSum())

    # Cosmetics for the plots
    histSig.SetLineColor(2)
    histSig.SetLineWidth(3)
    histSig.GetXaxis().SetRange(0, histSig.GetNbinsX()+1)
    histBkg.SetLineColor(4)
    histBkg.SetLineWidth(3)
    histBkg.GetXaxis().SetRange(0, histBkg.GetNbinsX()+1)

    histSig.GetXaxis().SetTitle(varName)
    histSig.GetYaxis().SetTitle("No. of events (Scaled to L=2.6fb-1)")
    canv.SetLogy()
    rt.gStyle.SetOptStat(0)

    # Fix the y axis range
    # axisRangeYmin = histSig.GetMaximum()
    # for hist in [histSig,histBkg]:
    #     if(hist.GetMaximum()<axisRangeYmin):
    #         axisRangeYmin = hist.GetMaximum()
    # histSig.SetMinimum(axisRangeYmin*0.0001)
    # histBkg.SetMinimum(axisRangeYmin*0.0001)
    histSig.SetMinimum(0.000001)
    histBkg.SetMinimum(0.000001)

    histSig.Draw("same hist E")
    histBkg.Draw("same hist E")

    legc1 = rt.TLegend(0.7, 0.9, 0.89, 1.0, "", "brNDC")
    legc1.AddEntry(histSig, "Signal", "l")
    legc1.AddEntry(histBkg, "Background", "l")
    legc1.SetTextSize(0.03)
    legc1.SetBorderSize(0)
    legc1.Draw()

    nx = histSig.GetNbinsX()+1
    bw1 = histSig.GetBinWidth(0)
    bw2 = histSig.GetBinWidth(nx)
    x1 = histSig.GetBinLowEdge(0)+bw1
    x2 = histSig.GetBinLowEdge(nx)+bw1
    y1 = histSig.GetBinContent(0)/histSig.Integral()
    y2 = histSig.GetBinContent(nx)/histSig.Integral()
    for hist in [histSig, histBkg]:
        y1 = max(y1,hist.GetBinContent(0)/hist.Integral())
        y2 = max(y2,hist.GetBinContent(nx)/hist.Integral())

    y1 *= 1.1
    y2 *= 1.1

    tUFlw = rt.TText(x1-0.5*bw1, y1, "Underflow")
    tUFlw.SetTextAngle(90)
    tUFlw.SetTextAlign(12)
    tUFlw.SetTextSize(0.03)
    tUFlw.Draw()
    
    tOFlw = rt.TText(x2-0.5*bw2, y1, "Overflow")
    tOFlw.SetTextAngle(90)
    tOFlw.SetTextAlign(12)
    tOFlw.SetTextSize(0.03)
    tOFlw.Draw()

    canv.SaveAs(varName+"_SR3.pdf")
