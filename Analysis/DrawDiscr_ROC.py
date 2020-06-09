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

from ROOT import TFile, TTree, TChain
from ROOT import TH1D, TH2D, TCanvas
from ROOT import gDirectory, TGraph, TMultiGraph
import ROOT as rt
import os
import sys
sys.path.insert(0, os.path.abspath('/home/arsahasransu/Documents/SoftDisplacedLeptons/Classifier/'))

print("All classes initialized successfully!!!")


import plotBeautifier as pB
pB.trial_func("AR")

bkg_chain = TChain("varTree")
bkg_chain.Add("../background.root")
bkg_Full = bkg_chain.AsMatrix()
    
# Load the input data scaler
scaler = joblib.load("../Classifier/scaler.save")

# Load the model
loaded_model = m.load_model("../Classifier/simplePer.h5")

bkg_scaled = scaler.transform(bkg_Full)
bkg_predict = loaded_model.predict(bkg_scaled)

bkg_sigprob = np.array(bkg_predict)[:,0]

def discr_modify(discr, val):

    modDiscr = discr
    for ctr in range(discr.shape[0]):

        modDiscr[ctr] = val if val<modDiscr[ctr] else modDiscr[ctr]

    return modDiscr

def discr_ROC_maker(rootFileName):

    sig_chain = TChain("varTree")
    
    sig_chain.Add("../"+rootFileName+".root")

    sig_Full = sig_chain.AsMatrix()
    sig_Full = np.transpose(sig_Full)
    sig_Full = sig_Full[2:11]
    sig_Full = np.transpose(sig_Full)

    sig_scaled = scaler.transform(sig_Full)

    sig_predict = loaded_model.predict(sig_scaled)

    sig_sigprob = np.array(sig_predict)[:,0]

    tpr = []
    fpr = []

    sigProb = np.arange(0, 1.01, 0.01)

    for x in sigProb:

        sig_class = sig_sigprob>=x
        bkg_class = bkg_sigprob>=x

        tp1 = sig_class.sum()
        fn1 = (1-sig_class).sum()
        tn1 = (1-bkg_class).sum()
        fp1 = bkg_class.sum()
        tpr.append(tp1/(tp1+fn1))
        fpr.append(fp1/(fp1+tn1))

    print(rootFileName+" completed.")
    return [sig_predict, tpr, fpr]

[sigprob_1, tpr_1, fpr_1] = discr_ROC_maker("BP_200_20_DM")
[sigprob_2, tpr_2, fpr_2] = discr_ROC_maker("BP_324_20_DM")
#[sigprob_3, tpr_3, fpr_3] = discr_ROC_maker("BP_200_20_2")
[sigprob_4, tpr_4, fpr_4] = discr_ROC_maker("BP_200_40_20")

nBins = 20
bkghisto = TH1D("","",nBins,0,1)
sig1histo = TH1D("","",nBins,0,1)
sig2histo = TH1D("","",nBins,0,1)
sig3histo = TH1D("","",nBins,0,1)
sig4histo = TH1D("","",nBins,0,1)

#discr_modify(bkg_predict[:,0], 0.999999)
bkghisto.FillN(bkg_predict.shape[0],(bkg_predict[:,0]).astype(float),np.ones(bkg_predict.shape[0]))
#discr_modify(sigprob_1[:,0], 0.999999)
sig1histo.FillN(sigprob_1.shape[0],(sigprob_1[:,0]).astype(float),np.ones(sigprob_1.shape[0]))
#discr_modify(sigprob_2[:,0], 0.999999)
sig2histo.FillN(sigprob_2.shape[0],(sigprob_2[:,0]).astype(float),np.ones(sigprob_2.shape[0]))
#discr_modify(sigprob_3[:,0], 0.999999)
#sig3histo.FillN(sigprob_3.shape[0],(sigprob_3[:,0]).astype(float),np.ones(sigprob_3.shape[0]))
#discr_modify(sigprob_4[:,0], 0.999999)
sig4histo.FillN(sigprob_4.shape[0],(sigprob_4[:,0]).astype(float),np.ones(sigprob_4.shape[0]))

histList = [bkghisto, sig1histo, sig2histo, sig4histo]
labelList = ["HF background", "DM: (200, 20)", "DM: (324, 20)", "(200, 40)"]
xAxisTitle = "signal probability"
yAxisTitle = "normalized number of events"
outPlotName = "discriminator"
pB.plotBeautifier(histList, labelList, xAxisTitle, yAxisTitle, outPlotName)

c1 = TCanvas("","",10,32,782,600)
rocS = TMultiGraph()
roc1 = TGraph(len(tpr_1), np.ones(len(fpr_1))-np.array(fpr_1), np.array(tpr_1))
roc1.SetLineColor(rt.kBlue)
roc1.SetLineWidth(3)
roc2 = TGraph(len(tpr_2), np.ones(len(fpr_2))-np.array(fpr_2), np.array(tpr_2))
roc2.SetLineColor(rt.kBlue-7)
roc2.SetLineWidth(3)
roc3 = TGraph(len(tpr_4), np.ones(len(fpr_4))-np.array(fpr_4), np.array(tpr_4))
roc3.SetLineColor(rt.kPink+7)
roc3.SetLineWidth(3)
rocS.Add(roc1)
rocS.Add(roc2)
rocS.Add(roc3)
rocS.Draw("AL")
c1.SaveAs("roc_ROOT.C")

plt.clf()
plt.plot(fpr_1, tpr_1, c="#0000ff", label="DM: (200, 20)")
plt.plot(fpr_2, tpr_2, c="#6666ff", label="DM: (324, 20)")
plt.plot(fpr_4, tpr_4, c="#ff0000", label="(200, 20)")
#plt.plot(fpr_4, tpr_4, c="#ff0099", label="(200, 40)")
#plt.xscale('log')
#plt.yscale('log')
plt.xlabel("background rejection")
plt.ylabel("signal efficiency")
plt.legend()
plt.savefig("roc.pdf")
