import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l

import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import joblib
from sklearn.metrics import confusion_matrix, fbeta_score
from scikitplot.metrics import plot_roc, plot_confusion_matrix

import seaborn as sn
import pandas as pd

import ROOT as rt
import copy

def plotCosmetics(histlist,
                  labellist,
                  xaxistitle,
                  yaxistitle,
                  outplotname,
                  logy):

    c1 = rt.TCanvas("","",10,32,782,552);

    ctr = 0
    for hist in histlist:
        hist.GetXaxis().SetTitle(xaxistitle)
        hist.GetYaxis().SetTitle(yaxistitle)
        if ctr<4:
            hist.SetLineColorAlpha(ctr+1,0.5)
        if ctr>=4:
            hist.SetLineColorAlpha(ctr+2,0.5)
        ctr = ctr+1
        hist.SetLineWidth(2)

    c1.SetLogy(logy)
    rt.gStyle.SetOptStat(0)

    for hist in histlist:
        hist.Draw("hist same")

    legc1 = rt.TLegend(0.5, 0.9, 0.89, 1.0, "", "brNDC")
    ctr = 0
    for hist in histlist:
        legc1.AddEntry(hist, labellist[ctr], "l")
        ctr = ctr+1
    legc1.SetTextSize(0.03)
    legc1.SetBorderSize(0)
    legc1.Draw()

    c1.SaveAs(outplotname+".pdf")
    return

def plotCosmeticswithDivide(histlist,
                            labellist,
                            xaxistitle,
                            yaxistitle,
                            outplotname,
                            logy):

    c1 = rt.TCanvas("","",1000,1000)
    p1 = rt.TPad("","",0,0.25,1,1)
    p2 = rt.TPad("","",0,0,1,0.25)
    p1.Draw()
    p2.Draw()

    p1.cd()
    ctr = 0
    for hist in histlist:
        hist.GetXaxis().SetTitle(xaxistitle)
        hist.GetYaxis().SetTitle(yaxistitle)
        if ctr<4:
            hist.SetLineColorAlpha(ctr+1,0.5)
        if ctr>=4:
            hist.SetLineColorAlpha(ctr+2,0.5)
        ctr = ctr+1
        hist.SetLineWidth(2)

    p1.SetLogy(logy)
    rt.gStyle.SetOptStat(0)

    for hist in histlist:
        hist.Draw("hist same")

    legc1 = rt.TLegend(0.5, 0.9, 0.89, 1.0, "", "brNDC")
    ctr = 0
    for hist in histlist:
        legc1.AddEntry(hist, labellist[ctr], "l")
        ctr = ctr+1
    legc1.SetTextSize(0.03)
    legc1.SetBorderSize(0)
    legc1.Draw()

    p2.cd()
    histratiolist = copy.deepcopy(histlist)
    for ctr in range(len(histratiolist)):
        histratiolist[ctr].Divide(histlist[0])
        histratiolist[ctr].GetYaxis().SetTitle("Ratio with unweighted histogram")
        histratiolist[ctr].GetYaxis().SetTitleSize(0.055)
        #histratiolist[ctr].SetMinimum(0.5)
        #histratiolist[ctr].SetMaximum(1)
        histratiolist[ctr].Draw("hist same")

    c1.SaveAs(outplotname+".pdf")
    return

# Load the signal and variable
sigchain = rt.TChain("varTree")
sigchain.Add("BP_324_20_DM.root")
sigfull = sigchain.AsMatrix()
sigfull = np.transpose(sigfull)
d0el = sigfull[0]
d0mu = sigfull[1]
var = sigfull[2:11]
# [eleCMS-11, muoCMS-12, eleF-13, muoF-14, eleN-15, muoN-16, eleK-17, muoK-18]
weight = sigfull[11:19]

# Plot the d0 electron
d0elhist = rt.TH1D("d0_el","",100,0,100)
d0elhist.FillN(d0el.shape[0],d0el.astype(float),np.ones(d0el.shape[0]))
d0elhist_CMS = rt.TH1D("d0_el_CMS","",100,0,100)
d0elhist_CMS.FillN(d0el.shape[0],d0el.astype(float),weight[0].astype(float))
d0elhist_F = rt.TH1D("d0_el_F","",100,0,100)
d0elhist_F.FillN(d0el.shape[0],d0el.astype(float),weight[2].astype(float))
d0elhist_N = rt.TH1D("d0_el_N","",100,0,100)
d0elhist_N.FillN(d0el.shape[0],d0el.astype(float),weight[4].astype(float))
d0elhist_K = rt.TH1D("d0_el_K","",100,0,100)
d0elhist_K.FillN(d0el.shape[0],d0el.astype(float),weight[6].astype(float))
histlist = [d0elhist, d0elhist_CMS, d0elhist_F, d0elhist_N, d0elhist_K]
labellist = ["Unweighted", "CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
plotCosmetics(histlist, labellist, "d0el", "events", "d0el", True)

d0elhist_CMS.Divide(d0elhist)
d0elhist_F.Divide(d0elhist)
d0elhist_N.Divide(d0elhist)
d0elhist_K.Divide(d0elhist)
histlist = [d0elhist_CMS, d0elhist_F, d0elhist_N, d0elhist_K]
labellist = ["CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
plotCosmetics(histlist, labellist, "d0el", "weight", "weightel", False)

# Plot the d0 muon
d0muhist = rt.TH1D("d0_mu","",100,0,100)
d0muhist.FillN(d0mu.shape[0],d0mu.astype(float),np.ones(d0mu.shape[0]))
d0muhist_CMS = rt.TH1D("d0_mu_CMS","",100,0,100)
d0muhist_CMS.FillN(d0mu.shape[0],d0mu.astype(float),weight[1].astype(float))
d0muhist_F = rt.TH1D("d0_mu_F","",100,0,100)
d0muhist_F.FillN(d0mu.shape[0],d0mu.astype(float),weight[3].astype(float))
d0muhist_N = rt.TH1D("d0_mu_N","",100,0,100)
d0muhist_N.FillN(d0mu.shape[0],d0mu.astype(float),weight[5].astype(float))
d0muhist_K = rt.TH1D("d0_mu_K","",100,0,100)
d0muhist_K.FillN(d0mu.shape[0],d0mu.astype(float),weight[7].astype(float))
histlist = [d0muhist, d0muhist_CMS, d0muhist_F, d0muhist_N, d0muhist_K]
labellist = ["Unweighted", "CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
plotCosmetics(histlist, labellist, "d0mu", "events", "d0mu", True)

d0muhist_CMS.Divide(d0muhist)
d0muhist_F.Divide(d0muhist)
d0muhist_N.Divide(d0muhist)
d0muhist_K.Divide(d0muhist)
histlist = [d0muhist_CMS, d0muhist_F, d0muhist_N, d0muhist_K]
labellist = ["CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
plotCosmetics(histlist, labellist, "d0mu", "weight", "weightmu", False)

# Plot the variables
weight_CMS = weight[0]*weight[1]
weight_F = weight[2]*weight[3]
weight_N = weight[4]*weight[5]
weight_K = weight[6]*weight[7]
varname = np.array(["ht","drll","dphilepmet","metsig","metsigselobj","alphat","sphericity","spherocity","mt"])
binmin = [20, 0, 0, 0, 0, 0, 0, 0, 0]
nbin = [16, 30, 50, 25, 20, 25, 50, 50, 25]
binmax = np.array([100, 5, rt.TMath.Pi(), 25, 10, 1.25, 1, 1, 200], dtype=float)
for ctr in range(varname.shape[0]):
    varhist = rt.TH1D("","",nbin[ctr],binmin[ctr],binmax[ctr])
    varhist_wtCMS = rt.TH1D("","",nbin[ctr],binmin[ctr],binmax[ctr])
    varhist_wtF = rt.TH1D("","",nbin[ctr],binmin[ctr],binmax[ctr])
    varhist_wtN = rt.TH1D("","",nbin[ctr],binmin[ctr],binmax[ctr])
    varhist_wtK = rt.TH1D("","",nbin[ctr],binmin[ctr],binmax[ctr])
    varhist.FillN(var[ctr].shape[0],var[ctr].astype(float),np.ones(var[ctr].shape[0]))
    varhist_wtCMS.FillN(var[ctr].shape[0],var[ctr].astype(float),weight_CMS.astype(float))
    varhist_wtF.FillN(var[ctr].shape[0],var[ctr].astype(float),weight_F.astype(float))
    varhist_wtN.FillN(var[ctr].shape[0],var[ctr].astype(float),weight_N.astype(float))
    varhist_wtK.FillN(var[ctr].shape[0],var[ctr].astype(float),weight_K.astype(float))

    histlist = [varhist, varhist_wtCMS, varhist_wtF, varhist_wtN, varhist_wtK]
    labellist = ["Unweighted", "CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
    plotCosmeticswithDivide(histlist, labellist, varname[ctr], "events", "VarPlots/"+varname[ctr], True)


# Load the input data scaler
scaler = joblib.load("./Classifier/scaler.save")

# Load the model
loaded_model = m.load_model("./Classifier/simplePer.h5")

# Plot the discriminator
var = np.transpose(var)
varscaled = scaler.transform(var)
modelpredict = loaded_model.predict(varscaled)
modelsigprob = np.array(modelpredict)[:,0]

nBins = 20
dischisto = rt.TH1D("","",nBins,0,1)
dischisto_wtCMS = rt.TH1D("","",nBins,0,1)
dischisto_wtF = rt.TH1D("","",nBins,0,1)
dischisto_wtN = rt.TH1D("","",nBins,0,1)
dischisto_wtK = rt.TH1D("","",nBins,0,1)

dischisto.FillN(modelsigprob.shape[0],modelsigprob.astype(float),np.ones(modelsigprob.shape[0]))
dischisto_wtCMS.FillN(modelsigprob.shape[0],modelsigprob.astype(float),weight_CMS.astype(float))
dischisto_wtF.FillN(modelsigprob.shape[0],modelsigprob.astype(float),weight_F.astype(float))
dischisto_wtN.FillN(modelsigprob.shape[0],modelsigprob.astype(float),weight_N.astype(float))
dischisto_wtK.FillN(modelsigprob.shape[0],modelsigprob.astype(float),weight_K.astype(float))
histlist = [dischisto, dischisto_wtCMS, dischisto_wtF, dischisto_wtN, dischisto_wtK]
labellist = ["Unweighted", "CMS HEP Database", "Freya weight scheme", "Nishita weight scheme", "Kamal weight scheme"]
plotCosmeticswithDivide(histlist, labellist, "signal probability", "events", "discriminator", True)
