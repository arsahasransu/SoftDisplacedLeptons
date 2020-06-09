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
from ROOT import gDirectory

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

    # Read from corresponding signal and background chain
    signal_Chain = TChain("varTree")
    signal_SR1_Chain = TChain("varTree_SR1")
    signal_SR2_Chain = TChain("varTree_SR2")
    signal_SR3_Chain = TChain("varTree_SR3")
    signal_Chain.Add("../"+rootFileName+".root")
    signal_SR1_Chain.Add("../"+rootFileName+".root")
    signal_SR2_Chain.Add("../"+rootFileName+".root")
    signal_SR3_Chain.Add("../"+rootFileName+".root")
    background_Chain = TChain("varTree")
    background_Chain.Add("../background.root")

    # Get the sample size
    signal_SR1_SampleSize = signal_SR1_Chain.GetEntries()
    signal_SR2_SampleSize = signal_SR2_Chain.GetEntries()
    signal_SR3_SampleSize = signal_SR3_Chain.GetEntries()
    background_SampleSize = background_Chain.GetEntries()

    # Include 2d D0 plots
    bins2D = np.array([0.01,0.1,0.2,0.5,1.0,100.0,1000.0])
    nbins2D = 6
    d02d = TH2D("d02d","d02d",nbins2D,bins2D,nbins2D,bins2D)
    d02d_SR1 = TH2D("d02d_SR1","d02d_SR1",nbins2D,bins2D,nbins2D,bins2D)
    d02d_SR2 = TH2D("d02d_SR2","d02d_SR2",nbins2D,bins2D,nbins2D,bins2D)
    d02d_SR3 = TH2D("d02d_SR3","d02d_SR3",nbins2D,bins2D,nbins2D,bins2D)
    signal_Chain.Draw("D0El:D0Mu>>d02d")
    signal_SR1_Chain.Draw("D0El_SR1:D0Mu_SR1>>d02d_SR1")
    signal_SR2_Chain.Draw("D0El_SR2:D0Mu_SR2>>d02d_SR2")
    signal_SR3_Chain.Draw("D0El_SR3:D0Mu_SR3>>d02d_SR3")
    d02d = gDirectory.Get("d02d")
    d02d_SR1 = gDirectory.Get("d02d_SR1")
    d02d_SR2 = gDirectory.Get("d02d_SR2")
    d02d_SR3 = gDirectory.Get("d02d_SR3")

    # Extract the variables to a numpy array
    signal_SR1_Full = signal_SR1_Chain.AsMatrix()
    signal_SR1_Full = np.transpose(signal_SR1_Full)
    weight_SR1 = signal_SR1_Full[11:19]
    signal_SR1_Full = signal_SR1_Full[2:11]
    signal_SR1_Full = np.transpose(signal_SR1_Full)
    #print(signal_SR1_Full.shape)
    signal_SR2_Full = signal_SR2_Chain.AsMatrix()
    signal_SR2_Full = np.transpose(signal_SR2_Full)
    weight_SR2 = signal_SR2_Full[11:19]
    signal_SR2_Full = signal_SR2_Full[2:11]
    signal_SR2_Full = np.transpose(signal_SR2_Full)
    #print(signal_SR2_Full.shape)
    signal_SR3_Full = signal_SR3_Chain.AsMatrix()
    signal_SR3_Full = np.transpose(signal_SR3_Full)
    weight_SR3 = signal_SR3_Full[11:19]
    signal_SR3_Full = signal_SR3_Full[2:11]
    signal_SR3_Full = np.transpose(signal_SR3_Full)
    #print(signal_SR3_Full.shape)
    background_Full = background_Chain.AsMatrix()
    #print(background_Full.shape)

    # Decide the weight scheme
    #wt_SR1 = weight_SR1[0]*weight_SR1[1] # CMS HEP Scheme
    #wt_SR1 = weight_SR1[2]*weight_SR1[3] # Freya Scheme
    #wt_SR1 = weight_SR1[4]*weight_SR1[5] # Nishita Scheme
    wt_SR1 = weight_SR1[6]*weight_SR1[7] # Kamal Scheme
    #wt_SR2 = weight_SR2[0]*weight_SR2[1]
    #wt_SR2 = weight_SR2[2]*weight_SR2[3]
    #wt_SR2 = weight_SR2[4]*weight_SR2[5]
    wt_SR2 = weight_SR2[6]*weight_SR2[7]
    #wt_SR3 = weight_SR3[0]*weight_SR3[1]
    #wt_SR3 = weight_SR3[2]*weight_SR3[3]
    #wt_SR3 = weight_SR3[4]*weight_SR3[5]
    wt_SR3 = weight_SR3[6]*weight_SR3[7]
    
    # Load the input data scaler
    scaler = joblib.load("../Classifier/scaler.save")

    # Load the model
    loaded_model = m.load_model("../Classifier/simplePer.h5")
    #loaded_model.summary()

    # Scale the variables
    signal_SR1_Scaled = scaler.transform(signal_SR1_Full)
    signal_SR2_Scaled = scaler.transform(signal_SR2_Full)
    signal_SR3_Scaled = scaler.transform(signal_SR3_Full)
    background_Scaled = scaler.transform(background_Full)

    # Prdict on the variables
    signal_SR1_Predict = loaded_model.predict(signal_SR1_Scaled)
    signal_SR2_Predict = loaded_model.predict(signal_SR2_Scaled)
    signal_SR3_Predict = loaded_model.predict(signal_SR3_Scaled)
    background_Predict = loaded_model.predict(background_Scaled)

    # Obtain the signal probability
    signal_SR1_SigProb = np.array(signal_SR1_Predict)[:,0] 
    signal_SR2_SigProb = np.array(signal_SR2_Predict)[:,0] 
    signal_SR3_SigProb = np.array(signal_SR3_Predict)[:,0] 

    # Cut on the discriminator to calculate the yield
    discCut = 0.0
    signalYield_SR1 = wt*(signal_SR1_SigProb>=discCut).sum()
    signalYield_SR2 = wt*(signal_SR2_SigProb>=discCut).sum()
    signalYield_SR3 = wt*(signal_SR3_SigProb>=discCut).sum()
    #print(rootFileName,"\t",round(signalYield_SR1,5),"\t",round(signalYield_SR2,5),"\t",round(signalYield_SR3,5))

    # Discriminator Shape in ROOT plotting with proper weight
    discFile = TFile("../"+rootFileName+"_wtK_Disc.root","RECREATE")
    nBins = 101
    signal_SR1_histo = TH1D("SR1","",nBins,0,1.01)
    signal_SR2_histo = TH1D("SR2","",nBins,0,1.01)
    signal_SR3_histo = TH1D("SR3","",nBins,0,1.01)
    signal_SR1_wt_histo = TH1D("SR1_wt","",nBins,0,1.01)
    signal_SR2_wt_histo = TH1D("SR2_wt","",nBins,0,1.01)
    signal_SR3_wt_histo = TH1D("SR3_wt","",nBins,0,1.01)
    background_histo = TH1D("background","",nBins,0,1.01)
    signal_SR1_histo.Sumw2()
    signal_SR2_histo.Sumw2()
    signal_SR3_histo.Sumw2()
    signal_SR1_wt_histo.Sumw2()
    signal_SR2_wt_histo.Sumw2()
    signal_SR3_wt_histo.Sumw2()
    signal_SR1_histo.FillN(signal_SR1_Predict.shape[0],(signal_SR1_Predict[:,0]).astype(float),np.ones(signal_SR1_Predict.shape[0]))
    signal_SR2_histo.FillN(signal_SR2_Predict.shape[0],(signal_SR2_Predict[:,0]).astype(float),np.ones(signal_SR2_Predict.shape[0]))
    signal_SR3_histo.FillN(signal_SR3_Predict.shape[0],(signal_SR3_Predict[:,0]).astype(float),np.ones(signal_SR3_Predict.shape[0]))
    signal_SR1_wt_histo.FillN(signal_SR1_Predict.shape[0],(signal_SR1_Predict[:,0]).astype(float),wt_SR1.astype(float))
    signal_SR2_wt_histo.FillN(signal_SR2_Predict.shape[0],(signal_SR2_Predict[:,0]).astype(float),wt_SR2.astype(float))
    signal_SR3_wt_histo.FillN(signal_SR3_Predict.shape[0],(signal_SR3_Predict[:,0]).astype(float),wt_SR3.astype(float))
    background_histo.FillN(background_Predict.shape[0],(background_Predict[:,0]).astype(float),np.ones(background_Predict.shape[0]))
    # Print the yield
    #print(rootFileName,
    #      "\t",round(signal_SR1_histo.Integral(),5),"\t",round(signal_SR1_wt_histo.Integral(),5),
    #      "\t",round(signal_SR2_histo.Integral(),5),"\t",round(signal_SR2_wt_histo.Integral(),5),
    #      "\t",round(signal_SR3_histo.Integral(),5),"\t",round(signal_SR3_wt_histo.Integral(),5))

    signal_SR1_wt_histo.Scale(signalYield_SR1/signal_SR1_histo.Integral())
    signal_SR2_wt_histo.Scale(signalYield_SR2/signal_SR2_histo.Integral())
    signal_SR3_wt_histo.Scale(signalYield_SR3/signal_SR3_histo.Integral())
    signal_SR1_histo.Scale(signalYield_SR1/signal_SR1_histo.Integral())
    signal_SR2_histo.Scale(signalYield_SR2/signal_SR2_histo.Integral())
    signal_SR3_histo.Scale(signalYield_SR3/signal_SR3_histo.Integral())
    background_histo.Scale(1.0/background_histo.Integral())

    # Print the yield
    #print(rootFileName,
    #      "\t",round(signal_SR1_histo.Integral(),5),"\t",round(signal_SR1_wt_histo.Integral(),5),
    #      "\t",round(signal_SR2_histo.Integral(),5),"\t",round(signal_SR2_wt_histo.Integral(),5),
    #      "\t",round(signal_SR3_histo.Integral(),5),"\t",round(signal_SR3_wt_histo.Integral(),5))

    print(rootFileName,
          "\t",round(signal_SR1_wt_histo.Integral(),5),
          "\t",round(signal_SR2_wt_histo.Integral(),5),
          "\t",round(signal_SR3_wt_histo.Integral(),5))

    signal_SR1_histo.Write()
    signal_SR2_histo.Write()
    signal_SR3_histo.Write()
    signal_SR1_wt_histo.Write()
    signal_SR2_wt_histo.Write()
    signal_SR3_wt_histo.Write()
    background_histo.Write()
    d02d.Write()
    d02d_SR1.Write()
    d02d_SR2.Write()
    d02d_SR3.Write()
    discFile.Close()

          

    #return [signalYield_SR1,signalYield_SR2,signalYield_SR3]

yieldCalc("BP_324_20_DM",128*0.025,20*100000)
yieldCalc("BP_200_20_DM",903*0.014,20*100000)
yieldCalc("BP_200_20_02",903,20*100000)
yieldCalc("BP_200_20_2",903,20*100000)
yieldCalc("BP_200_20_20",903,20*100000)
yieldCalc("BP_200_20_200",903,20*100000)
yieldCalc("BP_200_40_20",903,20*100000)
