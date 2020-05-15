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
from ROOT import TH1F, TCanvas

import os
import sys
sys.path.insert(0, os.path.abspath('/home/arsahasransu/Documents/SoftDisplacedLeptons/Classifier/'))

print("All classes initialized successfully!!!")


import plotBeautifier as pB
pB.trial_func("AR")


def yieldCalc(rootFileName, crossSec, nSimu):

    # Signal Significance Calculation
    luminosity = 2.6 # in fb^{-1}
    nSimu = 20*100000
    nEvent = luminosity*crossSec
    wt = nEvent/nSimu

    signal_SR1_Chain = TChain("varTree_SR1")
    signal_SR2_Chain = TChain("varTree_SR2")
    signal_SR3_Chain = TChain("varTree_SR3")
    signal_SR1_Chain.Add("../"+rootFileName+".root")
    signal_SR2_Chain.Add("../"+rootFileName+".root")
    signal_SR3_Chain.Add("../"+rootFileName+".root")
    background_Chain = TChain("varTree")
    background_Chain.Add("../background.root")

    signal_SR1_SampleSize = signal_SR1_Chain.GetEntries()
    signal_SR2_SampleSize = signal_SR2_Chain.GetEntries()
    signal_SR3_SampleSize = signal_SR3_Chain.GetEntries()
    background_SampleSize = background_Chain.GetEntries()

    signal_SR1_Full = signal_SR1_Chain.AsMatrix()
    #print(signal_SR1_Full.shape)
    signal_SR2_Full = signal_SR2_Chain.AsMatrix()
    #print(signal_SR2_Full.shape)
    signal_SR3_Full = signal_SR3_Chain.AsMatrix()
    #print(signal_SR3_Full.shape)
    background_Full = background_Chain.AsMatrix()
    #print(background_Full.shape)

    # Load the input data scaler
    scaler = joblib.load("../Classifier/scaler.save")

    # Load the model
    loaded_model = m.load_model("../Classifier/simplePer.h5")
    #loaded_model.summary()

    signal_SR1_Scaled = scaler.transform(signal_SR1_Full)
    signal_SR2_Scaled = scaler.transform(signal_SR2_Full)
    signal_SR3_Scaled = scaler.transform(signal_SR3_Full)
    background_Scaled = scaler.transform(background_Full)

    signal_SR1_Predict = loaded_model.predict(signal_SR1_Scaled)
    signal_SR2_Predict = loaded_model.predict(signal_SR2_Scaled)
    signal_SR3_Predict = loaded_model.predict(signal_SR3_Scaled)
    background_Predict = loaded_model.predict(background_Scaled)

    '''
    # In[10]:
    
    
    # Using plotBeautifier
    
    histList = [background_histo, signal_histo]
    labelList = ["HF_background", "(200,20,02)"]
    xAxisTitle = "Sig_Prob"
    yAxisTitle = "normalized number of events"
    outPlotName = "signal_Disc"
    pB.plotBeautifier(histList, labelList, xAxisTitle, yAxisTitle, outPlotName)
    
    
    # In[11]:
    
    
    # For plots based on discriminator cut
    
    sigProb = np.arange(0, 1, 0.001) 
    '''
    signal_SR1_SigProb = np.array(signal_SR1_Predict)[:,0] 
    signal_SR2_SigProb = np.array(signal_SR2_Predict)[:,0] 
    signal_SR3_SigProb = np.array(signal_SR3_Predict)[:,0] 

    '''
    background_SigProb = np.array(background_Predict)[:,0]
    
    tpr = [] 
    fpr = []
    
    for x in sigProb: 
    signal_Class = signal_SigProb>x
    background_Class = background_SigProb>x
    
    tp1 = signal_Class.sum()
    fn1 = (1-signal_Class).sum()
    tn1 = (1-background_Class).sum()
    fp1 = background_Class.sum()
    tpr.append(tp1/(tp1+fn1))
    fpr.append(fp1/(fp1+tn1))


    # In[13]:
    
    
    # ROC Curve
    
    plt.clf()
    plt.plot(fpr,tpr, label="(200,220,02)")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    #plt.xscale('log')
    plt.legend()
    plt.savefig("roc.pdf")
    
    
    # In[14]:

    '''
    discCut = 0
    signalYield_SR1 = wt*(signal_SR1_SigProb>discCut).sum()
    signalYield_SR2 = wt*(signal_SR2_SigProb>discCut).sum()
    signalYield_SR3 = wt*(signal_SR3_SigProb>discCut).sum()
    print(rootFileName,"\t",round(signalYield_SR1,5),"\t",round(signalYield_SR2,5),"\t",round(signalYield_SR3,5))

    # Discriminator Shape in ROOT plotting
    discFile = TFile("../"+rootFileName+"_Disc.root","RECREATE")
    nBins = 103
    signal_SR1_histo = TH1F("SR1","",nBins,-0.015,1.015)
    signal_SR2_histo = TH1F("SR2","",nBins,-0.015,1.015)
    signal_SR3_histo = TH1F("SR3","",nBins,-0.015,1.015)
    background_histo = TH1F("background","",nBins,0,1.01)
    signal_SR1_histo.Sumw2()
    signal_SR2_histo.Sumw2()
    signal_SR3_histo.Sumw2()
    signal_SR1_histo.FillN(signal_SR1_Predict.shape[0],(signal_SR1_Predict[:,0]).astype(float),np.ones(signal_SR1_Predict.shape[0]))
    signal_SR2_histo.FillN(signal_SR2_Predict.shape[0],(signal_SR2_Predict[:,0]).astype(float),np.ones(signal_SR2_Predict.shape[0]))
    signal_SR3_histo.FillN(signal_SR3_Predict.shape[0],(signal_SR3_Predict[:,0]).astype(float),np.ones(signal_SR3_Predict.shape[0]))
    background_histo.FillN(background_Predict.shape[0],(background_Predict[:,0]).astype(float),np.ones(background_Predict.shape[0]))
    signal_SR1_histo.Scale(signalYield_SR1/signal_SR1_histo.Integral())
    signal_SR2_histo.Scale(signalYield_SR2/signal_SR2_histo.Integral())
    signal_SR3_histo.Scale(signalYield_SR3/signal_SR3_histo.Integral())
    background_histo.Scale(1.0/background_histo.Integral())
    signal_SR1_histo.Write()
    signal_SR2_histo.Write()
    signal_SR3_histo.Write()
    background_histo.Write()
    discFile.Close()

    #return [signalYield_SR1,signalYield_SR2,signalYield_SR3]

yieldCalc("BP_200_20_DM",903*0.014,20*100000)
yieldCalc("BP_324_20_DM",128*0.025,20*100000)
yieldCalc("BP_200_20_02",903,20*100000)
yieldCalc("BP_200_20_2",903,20*100000)
yieldCalc("BP_200_20_20",903,20*100000)
yieldCalc("BP_200_20_200",903,20*100000)
yieldCalc("BP_200_40_20",903,20*100000)
