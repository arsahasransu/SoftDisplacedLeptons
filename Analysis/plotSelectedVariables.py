import os, sys
import ROOT as rt
import math

def plotSingleRegion(d0_region):
    # Read the variables to chain
    sigTree=rt.TChain("varTree","varTree")
    sigTree2mm=rt.TChain("varTree","varTree")
    sigTree20cm=rt.TChain("varTree","varTree")
    bgTree=rt.TChain("varTree","varTree")
    if d0_region==0:
        sigTree.Add("signal_BP_200_220_2cm_CR1.root")
        sigTree2mm.Add("signal_BP_200_220_2mm_CR1.root")
        sigTree20cm.Add("signal_BP_200_220_20cm_CR1.root")
        bgTree.Add("background.root")
    elif d0_region==1:
        sigTree.Add("signal_BP_200_220_2cm_CR2.root")
        sigTree2mm.Add("signal_BP_200_220_2mm_CR2.root")
        sigTree20cm.Add("signal_BP_200_220_20cm_CR2.root")
        bgTree.Add("background.root")
    elif d0_region==2:
        sigTree.Add("signal_BP_200_220_2cm_SR1.root")
        sigTree2mm.Add("signal_BP_200_220_2mm_SR1.root")
        sigTree20cm.Add("signal_BP_200_220_20cm_SR1.root")
        bgTree.Add("background.root")
    elif d0_region==3:
        sigTree.Add("signal_BP_200_220_2cm_SR2.root")
        sigTree2mm.Add("signal_BP_200_220_2mm_SR2.root")
        sigTree20cm.Add("signal_BP_200_220_20cm_SR2.root")
        bgTree.Add("background.root")
    elif d0_region==4:
        sigTree.Add("signal_BP_200_220_2cm_SR3.root")
        sigTree2mm.Add("signal_BP_200_220_2mm_SR3.root")
        sigTree20cm.Add("signal_BP_200_220_20cm_SR3.root")
        bgTree.Add("background.root")
    else:
        sigTree.Add("signal_BP_200_220_2cm.root")
        sigTree2mm.Add("signal_BP_180_220_2cm.root")
        sigTree20cm.Add("signal_BP_304_324_2cm.root")
        bgTree.Add("background.root")

    print("Now filled trees:",sigTree.GetEntries(),sigTree2mm.GetEntries(),sigTree20cm.GetEntries(),bgTree.GetEntries())
    
    # Create a canvas and list of variables
    canv = rt.TCanvas("canv","canv")
    varList = ["YDelpObj", "dRLL",
               "Sphericity", "Spherocity", "YUserObj",
               "dPhiLepMETSelObj", "alphaT",
               "MtLeadLepMET", "HtJet"]
    # varList = ["YDelpObj", "dRLL",
    #            "Sphericity", "Spherocity", "YUserObj",
    #            "dPhiLepMETSelObj", "alphaT"]
    nBinsList = [51, 30, 50, 50, 51, 51, 51, 51, 51]
    xMinList = [-1, -0.1, 0, 0, -1, -0.05, -0.02, -1, -1]
    xMaxList = [40, 4.5, 1, 1, 30, 3.5, 2.02, 101, 101]
    
    # Loop on variables
    ctr = 0
    for varName in varList:
        print("Plotting for variable: ", varName)
        
        nBins = nBinsList[ctr]
        minX = xMinList[ctr]
        maxX = xMaxList[ctr]
            
        histSig = rt.TH1D("histSig"+varName, "", nBins, minX, maxX)
        histSig.SetXTitle(varName)
        histSig2mm = rt.TH1D("histSig2mm"+varName, "", nBins, minX, maxX)
        histSig2mm.SetXTitle(varName)
        histSig20cm = rt.TH1D("histSig20cm"+varName, "", nBins, minX, maxX)
        histSig20cm.SetXTitle(varName)
        histBkg = rt.TH1D("histBkg"+varName, "", nBins, minX, maxX)
        histBkg.SetXTitle(varName)

        drawCommand = varName + ">>"+histSig.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        sigTree.Draw(drawCommand)
        drawCommand = varName + ">>"+histSig2mm.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        sigTree2mm.Draw(drawCommand)
        drawCommand = varName + ">>"+histSig20cm.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        sigTree20cm.Draw(drawCommand)
        drawCommand = varName + ">>"+histBkg.GetName()+"("+str(nBins)+","+str(minX)+","+str(maxX)+")"
        bgTree.Draw(drawCommand)
        histSig = rt.gDirectory.Get("histSig"+varName)
        histSig2mm = rt.gDirectory.Get("histSig2mm"+varName)
        histSig20cm = rt.gDirectory.Get("histSig20cm"+varName)
        histBkg = rt.gDirectory.Get("histBkg"+varName)
        
        # Normalisation to Luminosity
        if d0_region==1:
            N_bkg = 1937.09 # CR2
            N_sig = 2.6*3.8*math.pow(10,-1)*7351/(2*pow(10,6)) # SR3
            N_sig2mm = 2.6*3.8*math.pow(10,-3)*145/(2*pow(10,6)) # SR3
            N_sig20cm = 2.6*3.8*math.pow(10,1)*16128/(2*pow(10,6)) # SR3
            histBkg.Scale(N_bkg/histBkg.Integral())
            histSig.Scale(N_sig/histSig.Integral())
            histSig2mm.Scale(N_sig2mm/histSig2mm.Integral())
            histSig20cm.Scale(N_sig20cm/histSig20cm.Integral())
        elif d0_region==2:
            N_bkg = 3646.28 # SR1
            N_sig = 2.6*3.8*math.pow(10,-1)*14706/(2*pow(10,6)) # SR3
            N_sig2mm = 2.6*3.8*math.pow(10,-3)*3803/(2*pow(10,6)) # SR3
            N_sig20cm = 2.6*3.8*math.pow(10,1)*18315/(2*pow(10,6)) # SR3
            histBkg.Scale(N_bkg/histBkg.Integral())
            histSig.Scale(N_sig/histSig.Integral())
            histSig2mm.Scale(N_sig2mm/histSig2mm.Integral())
            histSig20cm.Scale(N_sig20cm/histSig20cm.Integral())
        elif d0_region==3:
            N_bkg = 569.73 # SR2
            N_sig = 2.6*3.8*math.pow(10,-1)*7351/(2*pow(10,6)) # SR3
            N_sig2mm = 2.6*3.8*math.pow(10,-3)*145/(2*pow(10,6)) # SR3
            N_sig20cm = 2.6*3.8*math.pow(10,1)*16128/(2*pow(10,6)) # SR3
            histBkg.Scale(N_bkg/histBkg.Integral())
            histSig.Scale(N_sig/histSig.Integral())
            histSig2mm.Scale(N_sig2mm/histSig2mm.Integral())
            histSig20cm.Scale(N_sig20cm/histSig20cm.Integral())
        elif d0_region==4:
            N_bkg = 22.79 # SR3
            N_sig = 2.6*3.8*math.pow(10,-1)*7351/(2*pow(10,6)) # SR3
            N_sig2mm = 2.6*3.8*math.pow(10,-3)*145/(2*pow(10,6)) # SR3
            N_sig20cm = 2.6*3.8*math.pow(10,1)*16128/(2*pow(10,6)) # SR3
            histBkg.Scale(N_bkg/histBkg.Integral())
            histSig.Scale(N_sig/histSig.Integral())
            histSig2mm.Scale(N_sig2mm/histSig2mm.Integral())
            histSig20cm.Scale(N_sig20cm/histSig20cm.Integral())
        else:
            N_bkg = 10000 # All
            N_sig = 2.6*3.8*math.pow(10,-1)*634178/(2*pow(10,6)) # SR3
            N_sig2mm = 2.6*3.8*math.pow(10,-3)*633891/(2*pow(10,6)) # SR3
            N_sig20cm = 2.6*3.8*math.pow(10,1)*627644/(2*pow(10,6)) # SR3
            histSig.Scale(1.0/histSig.GetSum())
            histSig2mm.Scale(1.0/histSig2mm.Integral())
            histSig20cm.Scale(1.0/histSig20cm.Integral())
            histBkg.Scale(1.0/histBkg.GetSum())
        
        # Cosmetics for the plots
        histSig.SetLineColor(2)
        histSig.SetLineWidth(3)
        histSig.GetXaxis().SetRange(0, histSig.GetNbinsX()+1)
        histSig2mm.SetLineColor(3)
        histSig2mm.SetLineWidth(3)
        histSig2mm.GetXaxis().SetRange(0, histSig2mm.GetNbinsX()+1)
        histSig20cm.SetLineColor(4)
        histSig20cm.SetLineWidth(3)
        histSig20cm.GetXaxis().SetRange(0, histSig20cm.GetNbinsX()+1)
        histBkg.SetLineColor(5)
        histBkg.SetLineWidth(3)
        histBkg.GetXaxis().SetRange(0, histBkg.GetNbinsX()+1)
        
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
        histSig.SetMinimum(0.00000001)
        histSig2mm.SetMinimum(0.00000001)
        histSig20cm.SetMinimum(0.00000001)
        histBkg.SetMinimum(0.00000001)
        
        histSig.Draw("same hist E")
        histSig2mm.Draw("same hist E")
        histSig20cm.Draw("same hist E")
        histBkg.Draw("same hist E")
        
        legc1 = rt.TLegend(0.7, 0.9, 0.89, 1.0, "", "brNDC")
        legc1.AddEntry(histSig20cm, "Signal BP=(304,324,2cm)", "l")
        legc1.AddEntry(histSig, "Signal BP=(200,220,2cm)", "l")
        legc1.AddEntry(histSig2mm, "Signal BP=(180,220,2cm)", "l")
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
        for hist in [histSig, histSig2mm, histSig20cm, histBkg]:
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
            canv.SaveAs("./Analysis/VarPlots_CR1/"+varName+"_CR1.pdf")
        elif d0_region==1:
            canv.SaveAs("./Analysis/VarPlots_CR2/"+varName+"_CR2.pdf")
        elif d0_region==2:
            canv.SaveAs("./Analysis/VarPlots_SR1/"+varName+"_SR1.pdf")
        elif d0_region==3:
            canv.SaveAs("./Analysis/VarPlots_SR2/"+varName+"_SR2.pdf")
        elif d0_region==4:
            canv.SaveAs("./Analysis/VarPlots_SR3/"+varName+"_SR3.pdf")
        else:
            canv.SaveAs("./Analysis/VarPlots/"+varName+".pdf")


        ctr = ctr+1
            
plotSingleRegion(-1)
#plotSingleRegion(2)
#plotSingleRegion(3)
#plotSingleRegion(4)
