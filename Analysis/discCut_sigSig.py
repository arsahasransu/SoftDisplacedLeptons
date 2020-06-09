import ROOT as rt
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

# Load the input data scaler
scaler = joblib.load("./Classifier/scaler.save")

# Load the model
loaded_model = m.load_model("./Classifier/simplePer.h5")
#loaded_model.summary()

normBGSR1 = 4122.8
normBGSR2 = 644.2
normBGSR3 = 24.479

plt.clf()

def sigsig(filename, sigsr1, sigsr2, sigsr3, bkgsr1, bkgsr2, bkgsr3):

    lumivals = [2.6, 5, 10, 30, 90, 100, 300, 900, 1000, 3000, 9000]
    limitvals = []
    
    nsig = sigsr1+sigsr2+sigsr3
    nbkg = bkgsr1+bkgsr2+bkgsr3

    for lumi in lumivals:
        lumi = lumi/lumivals[0]
        sigfrombkgonly = bkgsr1*lumi*rt.TMath.Log(1+sigsr1/bkgsr1)
        sigfrombkgonly = sigfrombkgonly + bkgsr2*lumi*rt.TMath.Log(1+sigsr2/bkgsr2)
        sigfrombkgonly = sigfrombkgonly + bkgsr3*lumi*rt.TMath.Log(1+sigsr3/bkgsr3)

        qbkgonly = -2*(sigfrombkgonly-nsig*lumi)
        qlim = 5.99

        limitvals.append(qbkgonly/qlim)

    plt.plot(lumivals, limitvals, label=filename)

def yieldcalc(rootFileName, crossSec, nSimu):

    signal_file = rt.TFile("./"+rootFileName+".root","READ")

    # Signal Significance Calculation
    luminosity = 2.6 # in fb^{-1}
    nSimu = 20*100000
    nEvent = luminosity*crossSec
    wt = nEvent/nSimu
    
    # Read from corresponding signal and background chain
    signal_Chain = rt.TChain("varTree")
    signal_SR1_Chain = rt.TChain("varTree_SR1")
    signal_SR2_Chain = rt.TChain("varTree_SR2")
    signal_SR3_Chain = rt.TChain("varTree_SR3")
    signal_Chain.Add("./"+rootFileName+".root")
    signal_SR1_Chain.Add("./"+rootFileName+".root")
    signal_SR2_Chain.Add("./"+rootFileName+".root")
    signal_SR3_Chain.Add("./"+rootFileName+".root")
    background_Chain = rt.TChain("varTree")
    background_Chain.Add("./background.root")

    # Get the sample size
    signal_SR1_SampleSize = signal_SR1_Chain.GetEntries()
    signal_SR2_SampleSize = signal_SR2_Chain.GetEntries()
    signal_SR3_SampleSize = signal_SR3_Chain.GetEntries()
    background_SampleSize = background_Chain.GetEntries()

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
    wt_SR1 = weight_SR1[0]*weight_SR1[1] # CMS HEP Scheme
    #wt_SR1 = weight_SR1[2]*weight_SR1[3] # Freya Scheme
    #wt_SR1 = weight_SR1[4]*weight_SR1[5] # Nishita Scheme
    #wt_SR1 = weight_SR1[6]*weight_SR1[7] # Kamal Scheme
    wt_SR2 = weight_SR2[0]*weight_SR2[1]
    #wt_SR2 = weight_SR2[2]*weight_SR2[3]
    #wt_SR2 = weight_SR2[4]*weight_SR2[5]
    #wt_SR2 = weight_SR2[6]*weight_SR2[7]
    wt_SR3 = weight_SR3[0]*weight_SR3[1]
    #wt_SR3 = weight_SR3[2]*weight_SR3[3]
    #wt_SR3 = weight_SR3[4]*weight_SR3[5]
    #wt_SR3 = weight_SR3[6]*weight_SR3[7]
    
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
    background_SigProb = np.array(background_Predict)[:,0]

    # Cut on the discriminator to calculate the yield
    discCut = 0.9
    signalYield_SR1 = wt*((signal_SR1_SigProb>=discCut)*(wt_SR1)).sum()
    signalYield_SR2 = wt*((signal_SR2_SigProb>=discCut)*(wt_SR2)).sum()
    signalYield_SR3 = wt*((signal_SR3_SigProb>=discCut)*(wt_SR3)).sum()
    wtbkg = (background_SigProb>=0.0).sum()
    backgroundYield = (1.0/wtbkg)*(background_SigProb>=discCut).sum()
    print(rootFileName,"\t",round(signalYield_SR1,5),"\t",round(signalYield_SR2,5),"\t",round(signalYield_SR3,5),"\t",backgroundYield*normBGSR1,"\t",backgroundYield*normBGSR2,"\t",backgroundYield*normBGSR3)

    # Calculate the signal significance
    sigsig(rootFileName,signalYield_SR1,signalYield_SR2,signalYield_SR3,backgroundYield*normBGSR1,backgroundYield*normBGSR2,backgroundYield*normBGSR3)
    
yieldcalc("BP_324_20_DM",128*0.025,20*100000)
yieldcalc("BP_200_20_DM",903*0.014,20*100000)
yieldcalc("BP_200_20_02",903,20*100000)
yieldcalc("BP_200_20_2",903,20*100000)
yieldcalc("BP_200_20_20",903,20*100000)
yieldcalc("BP_200_20_200",903,20*100000)
yieldcalc("BP_200_40_20",903,20*100000)

lumivals = [2.6, 5, 10, 30, 90, 100, 300, 900, 1000, 3000, 9000]
limitvals = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
plt.plot(lumivals, limitvals, linestyle='dashed')
plt.xlabel('luminosity')
plt.ylabel('Q/5.99')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001,1000)
#plt.legend()
plt.savefig("limitCMSHEP.pdf")
