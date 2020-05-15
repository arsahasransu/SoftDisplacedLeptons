import os, sys
import ROOT as rt
import math

NNcut = 0.950010
normBGSR1 = 4122.8
normBGSR2 = 644.2
normBGSR3 = 24.479

filelist = ("../../SoftDisplacedLeptons/BP_200_20_02_Disc.root",    "../../SoftDisplacedLeptons/BP_200_20_DM_Disc.root",
    "../../SoftDisplacedLeptons/BP_200_20_200_Disc.root",
    "../../SoftDisplacedLeptons/BP_200_40_20_Disc.root",
    "../../SoftDisplacedLeptons/BP_200_20_20_Disc.root",
    "../../SoftDisplacedLeptons/BP_324_20_DM_Disc.root",
    "../../SoftDisplacedLeptons/BP_200_20_2_Disc.root"
)
lumivals = {2.9, 36., 140., 300.}
for filename in filelist:
    print(filename)
    file = rt.TFile(filename,"READ")
#    file.ls()
    
    histosigSR1 = file.Get("SR1")
    histosigSR2 = file.Get("SR2")
    histosigSR3 = file.Get("SR3")
    histobg = file.Get("background")
    histobgSR1 = histobg.Clone("histobgSR1")
    histobgSR1.Scale(normBGSR1)
    histobgSR2 = histobg.Clone("histobgSR2")
    histobgSR2.Scale(normBGSR2)
    histobgSR3 = histobg.Clone("histobgSR3")
    histobgSR3.Scale(normBGSR3)
    


    sumSR1 = histosigSR1.Clone("sumSR1")
    sumSR2 = histosigSR2.Clone("sumSR2")
    sumSR3 = histosigSR3.Clone("sumSR3")
    sumBGSR1 = histobgSR1.Clone("sumBGSR1")
    sumBGSR2 = histobgSR2.Clone("sumBGSR2")
    sumBGSR3 = histobgSR3.Clone("sumBGSR3")
    sqrtsumBGSR1 = histobgSR1.Clone("sqrtsumBGSR1")
    sqrtsumBGSR2 = histobgSR2.Clone("sqrtsumBGSR2")
    sqrtsumBGSR3 = histobgSR3.Clone("sqrtsumBGSR3")


    for lumival in lumivals:

        lumicorr = lumival/2.9
#    print(round(histobg.Integral(),3),round(histosigSR1.Integral(),3),round(histosigSR2.Integral(),3),round(histosigSR3.Integral(),3))
        maxbins =histobg.GetNbinsX()+2


        for ibin in range(0,maxbins):

            sumSR1.SetBinContent(ibin,lumicorr*histosigSR1.Integral(ibin,maxbins))
            sumSR1.SetBinError(ibin,0)
            sumSR2.SetBinContent(ibin,lumicorr*histosigSR2.Integral(ibin,maxbins))
            sumSR2.SetBinError(ibin,0)
            sumSR3.SetBinContent(ibin,lumicorr*histosigSR3.Integral(ibin,maxbins))
            sumSR3.SetBinError(ibin,0)
            
            sumBGSR1.SetBinContent(ibin,lumicorr*histobgSR1.Integral(ibin,maxbins))
            sumBGSR1.SetBinError(ibin,0)
            sumBGSR2.SetBinContent(ibin,lumicorr*histobgSR3.Integral(ibin,maxbins))
            sumBGSR2.SetBinError(ibin,0)
            sumBGSR3.SetBinContent(ibin,lumicorr*histobgSR3.Integral(ibin,maxbins))
            sumBGSR3.SetBinError(ibin,0)
            
            sqrtsumBGSR1.SetBinContent(ibin,math.sqrt(lumicorr*histobgSR1.Integral(ibin,maxbins)))
            sqrtsumBGSR1.SetBinError(ibin,0)
            sqrtsumBGSR2.SetBinContent(ibin,math.sqrt(lumicorr*histobgSR3.Integral(ibin,maxbins)))
            sqrtsumBGSR2.SetBinError(ibin,0)
            sqrtsumBGSR3.SetBinContent(ibin,math.sqrt(lumicorr*histobgSR3.Integral(ibin,maxbins)))
            sqrtsumBGSR3.SetBinError(ibin,0)
            
    #        if ibin %10 == 0 or ibin< 10 or ibin>maxbins-5:
    #            print(ibin,histobg.GetBinLowEdge(ibin),histobg.GetBinWidth(ibin), histobg.GetBinContent(ibin),histosigSR1.GetBinContent(ibin),histosigSR2.GetBinContent(ibin),histosigSR3.GetBinContent(ibin))
    #            print("bin:",ibin,"low edge", histobg.GetBinLowEdge(ibin),"bin width", histobg.GetBinWidth(ibin),
    #               "bg:", round(histobg.Integral(ibin,maxbins),3),
    #                "SR1bg:",round(histobgSR1.Integral(ibin,maxbins),3),
    #                "SR2bg:",round(histobgSR2.Integral(ibin,maxbins),3),
    #                "SR3bg:",round(histobgSR3.Integral(ibin,maxbins),3),
    #                "SR1:",round(histosigSR1.Integral(ibin,maxbins),3),
    #                "SR2:",round(histosigSR2.Integral(ibin,maxbins),3),
    #                "SR3:",round(histosigSR3.Integral(ibin,maxbins),3))
        

        binval095 = sumSR1.GetXaxis().FindBin(NNcut)
        print("cutting at ", NNcut," : true bin cut value: ",sumSR1.GetBinLowEdge(binval095) )
        signifSR1 =sumSR1.Clone("signifSR1")
        signifSR1.Divide(sqrtsumBGSR1)
        signifSR2 =sumSR2.Clone("signifSR2")
        signifSR2.Divide(sqrtsumBGSR2)
        signifSR3 =sumSR3.Clone("signifSR3")
        signifSR3.Divide(sqrtsumBGSR3)
        canv1 = rt.TCanvas()
        #sumBGSR1.Draw()
        #sumSR1.SetLineColor(rt.kPink)
        #sumSR1.Draw("same")
        signifSR1.SetTitle(filename+" SR1")
        signifSR1.Draw()
        canv1.Update()

        canv2 = rt.TCanvas()
    #    sumBGSR2.Draw()
    #    sumSR2.SetLineColor(rt.kPink)
    #    sumSR2.Draw("same")
        signifSR2.SetTitle(filename+" SR2")
        signifSR2.Draw()
        canv2.Update()



        canv3 = rt.TCanvas()
    #    sumBGSR3.Draw()
    #    sumSR3.SetLineColor(rt.kPink)
    #    sumSR3.Draw("same")
        signifSR3.SetTitle(filename+" SR3")
        signifSR3.Draw()
        canv3.Update()
        
        print("for luminosity: ", lumicorr*2.9," fb-1, model" , filename)
        print("significance SR1 with NN>",round(sumSR1.GetBinLowEdge(binval095),2),": ",round(signifSR1.GetBinContent(binval095),3)," signal: ",round(sumSR1.GetBinContent(binval095),3),"background: ",round(sumBGSR1.GetBinContent(binval095),3) )
        print("significance SR2 with NN>",round(sumSR2.GetBinLowEdge(binval095),2),": ",round(signifSR2.GetBinContent(binval095),3)," signal: ",round(sumSR2.GetBinContent(binval095),3),"background: ",round(sumBGSR2.GetBinContent(binval095),3) )
        print("significance SR3 with NN>",round(sumSR3.GetBinLowEdge(binval095),2),": ",round(signifSR3.GetBinContent(binval095),3)," signal: ",round(sumSR3.GetBinContent(binval095),3),"background: ",round(sumBGSR3.GetBinContent(binval095),3) )

