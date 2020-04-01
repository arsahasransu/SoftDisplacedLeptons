import os, sys
import ROOT as rt
import math

def plotSingleRegion(d0_region):
    # Read the variables to chain
    sigTree=rt.TChain("varTree","varTree")
    sigTreeBP_100_200=rt.TChain("varTree","varTree")
    bgTree=rt.TChain("varTree","varTree")
    if d0_region==0:
        sigTree.Add("signal_CR1.root")
        sigTreeBP_100_200.Add("signal_BP_100_200_CR1.root")
        bgTree.Add("background_CR1.root")
    elif d0_region==1:
        sigTree.Add("signal_CR2.root")
        sigTreeBP_100_200.Add("signal_BP_100_200_CR2.root")
        bgTree.Add("background_CR2.root")
    elif d0_region==2:
        sigTree.Add("signal_SR1.root")
        sigTreeBP_100_200.Add("signal_BP_100_200_SR1.root")
        bgTree.Add("background_SR1.root")
    elif d0_region==3:
        sigTree.Add("signal_SR2.root")
        sigTreeBP_100_200.Add("signal_BP_100_200_SR2.root")
        bgTree.Add("background_SR2.root")
    elif d0_region==4:
        sigTree.Add("signal_SR3.root")
        sigTreeBP_100_200.Add("signal_BP_100_200_SR3.root")
        bgTree.Add("background_SR3.root")
    else:
        sigTree.Add("signal.root")
        sigTreeBP_100_200.Add("signal_BP_100_200.root")
        bgTree.Add("background.root")
        
    print("Now filled trees:",sigTree.GetEntries(),sigTreeBP_100_200.GetEntries(),bgTree.GetEntries())
    
    # Create a canvas and list of variables
    canv = rt.TCanvas("canv","canv")
    varList = ["YDelpObj", "dRLL",
               "Sphericity", "Spherocity", "YUserObj",
               "dPhiLepMETSelObj", "alphaT",
               "MtLeadLepMET", "HtJet"]
    nBinsList = [51, 51, 50, 50, 51, 51, 51, 51, 51]
    xMinList = [-1, -0.1, 0, 0, -1, -0.05, -0.02, -1, -1]
    xMaxList = [40, 6, 1, 1, 30, 3.5, 2.02, 101, 101]
    
    # Loop on variables
    ctr = 0
    for varName in varList:
        print("Plotting for variable: ", varName)
        
        sigTree.Draw(varName+">>sigHisto")
        sigHisto = rt.gDirectory.Get("sigHisto")
        sigTreeBP_100_200.Draw(varName+">>sigHisto_BP_100_200")
        sigHisto_BP_100_200 = rt.gDirectory.Get("sigHisto_BP_100_200")
        bgTree.Draw(varName+">>bgHisto")
        bgHisto = rt.gDirectory.Get("bgHisto")
        minX = sigHisto.GetXaxis().GetXmin()
        maxX = sigHisto.GetXaxis().GetXmax()
#        if minX > sigHisto.GetXaxis().GetXmin() :
#            minX = sigHisto.GetXaxis().GetXmin()
#            if maxX < sigHisto.GetXaxis().GetXmax() :
#                maxX = sigHisto.GetXaxis().GetXmax()
        nBins = sigHisto.GetXaxis().GetNbins()
        for hist in sigHisto, sigHisto_BP_100_200, bgHisto:
            minX = minX if minX<sigHisto.GetXaxis().GetXmin() else sigHisto.GetXaxis().GetXmin()
            maxX = maxX if maxX>sigHisto.GetXaxis().GetXmax() else sigHisto.GetXaxis().GetXmax()
            nBins = nBins if nBins<sigHisto.GetXaxis().GetNbins() else sigHisto.GetXaxis().GetNbins() 

        nBins = nBinsList[ctr]
        minX = xMinList[ctr]
        maxX = xMaxList[ctr]
            
        sigHisto.Delete()
        sigHisto_BP_100_200.Delete()
        bgHisto.Delete()

        histSig = rt.TH1D("histSig"+varName, "", nBins, minX, maxX)
        histSig.SetXTitle(varName)
        histSigBP_100_200 = rt.TH1D("histSigBP_100_200"+varName, "", nBins, minX, maxX)
        histSigBP_100_200.SetXTitle(varName)
        histBkg = rt.TH1D("histBkg"+varName, "", nBins, minX, maxX)
        histBkg.SetXTitle(varName)
        drawCommand = varName + ">>"+histSig.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        sigTree.Draw(drawCommand)
        drawCommand = varName + ">>"+histSigBP_100_200.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        sigTreeBP_100_200.Draw(drawCommand)
        drawCommand = varName + ">>"+histBkg.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        bgTree.Draw(drawCommand)
        histSig = rt.gDirectory.Get("histSig"+varName)
        histSigBP_100_200 = rt.gDirectory.Get("histSigBP_100_200"+varName)
        histBkg = rt.gDirectory.Get("histBkg"+varName)
        
        # Normalisation to Luminosity
        if d0_region==1:
            N_bkg = 1937.09 # CR2
            N_sig = 31.2*271/(5*pow(10,6)) # CR2
            N_sig_BP_100_200 = 202.8*164/(1*pow(10,6)) # CR2
            histSig.Scale(N_sig/histSig.Integral())
            histBkg.Scale(N_bkg/histBkg.Integral())
        elif d0_region==2:
            N_bkg = 3646.28 # SR1
            N_sig = 31.2*28976/(5*pow(10,6)) # SR1
            N_sig_BP_100_200 = 202.8*164/(1*pow(10,6)) # SR1
            histSig.Scale(N_sig/histSig.Integral())
            histSigBP_100_200.Scale(N_sig_BP_100_200/histSigBP_100_200.Integral())
            histBkg.Scale(N_bkg/histBkg.Integral())
        elif d0_region==3:
            N_bkg = 569.73 # SR2
            N_sig = 31.2*19077/(5*pow(10,6)) # SR2
            histSig.Scale(N_sig/histSig.Integral())
            histBkg.Scale(N_bkg/histBkg.Integral())
        elif d0_region==4:
            N_bkg = 22.79 # SR3
            N_sig = 31.2*10860/(5*pow(10,6)) # SR3
            histSig.Scale(N_sig/histSig.Integral())
            histBkg.Scale(N_bkg/histBkg.Integral())
        else:
            N_bkg = 10000 # All
            N_sig = 31.2*44440/(5*pow(10,6)) # All
            N_sig_BP_100_200 = 202.8*228099/(1*pow(10,6)) # All
            histSig.Scale(1.0/histSig.GetSum())
            histSigBP_100_200.Scale(1.0/histSigBP_100_200.GetSum())
            histBkg.Scale(1.0/histBkg.GetSum())
        
        # Cosmetics for the plots
        histSig.SetLineColor(2)
        histSig.SetLineWidth(3)
        histSig.GetXaxis().SetRange(0, histSig.GetNbinsX()+1)
        histBkg.SetLineColor(4)
        histBkg.SetLineWidth(3)
        histBkg.GetXaxis().SetRange(0, histBkg.GetNbinsX()+1)
        histSigBP_100_200.SetLineColor(5)
        histSigBP_100_200.SetLineWidth(3)
        histSigBP_100_200.GetXaxis().SetRange(0, histSigBP_100_200.GetNbinsX()+1)
        
        histBkg.GetXaxis().SetTitle(varName)
        if d0_region==-1:
            histBkg.GetYaxis().SetTitle("No. of events (Scaled to 1)")
        else:
            histBkg.GetYaxis().SetTitle("No. of events (Scaled to L=2.6fb-1)")
            
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
        histSigBP_100_200.SetMinimum(0.000001)
        histBkg.SetMinimum(0.000001)
        
        histSig.Draw("same hist E")
        histSigBP_100_200.Draw("same hist E")
        histBkg.Draw("same hist E")
        
        legc1 = rt.TLegend(0.7, 0.9, 0.89, 1.0, "", "brNDC")
        legc1.AddEntry(histSig, "Signal BP=(304,324)", "l")
        legc1.AddEntry(histSigBP_100_200, "Signal BP=(100,200)", "l")
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
        
        if d0_region==0:
            canv.SaveAs("VarPlots_CR1/"+varName+"_CR1.pdf")
        elif d0_region==1:
            canv.SaveAs("VarPlots_CR2/"+varName+"_CR2.pdf")
        elif d0_region==2:
            canv.SaveAs("VarPlots_SR1/"+varName+"_SR1.pdf")
        elif d0_region==3:
            canv.SaveAs("VarPlots_SR2/"+varName+"_SR2.pdf")
        elif d0_region==4:
            canv.SaveAs("VarPlots_SR3/"+varName+"_SR3.pdf")
        else:
            canv.SaveAs("VarPlots/"+varName+".pdf")

        ctr = ctr+1
            
plotSingleRegion(-1)
plotSingleRegion(2)
#plotSingleRegion(3)
#plotSingleRegion(4)
