#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l

import math as mt
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import joblib
from sklearn.metrics import confusion_matrix, fbeta_score
from scikitplot.metrics import plot_roc, plot_confusion_matrix

import seaborn as sn
import pandas as pd

from ROOT import TFile, TTree, TChain

print("All classes initialized successfully!!!")


# In[2]:


sigChan3 = TChain("varTree")
sigChan3.Add("signal_SR3.root")
bkgChan3 = TChain("varTree")
bkgChan3.Add("background_SR3.root")
print("Data read from the trees. Printing out the contents.")


# In[3]:


sigChan3.Print()
bkgChan3.Print()


# In[4]:


sig3SampleSize = sigChan3.GetEntries()
bkg3SampleSize = bkgChan3.GetEntries()

sig3Full = sigChan3.AsMatrix()
bkg3Full = bkgChan3.AsMatrix()


# In[5]:


# Load the input data scaler
scaler = joblib.load("../scaler.save")

# Load the model
loaded_model = m.load_model("../simplePer.h5")
loaded_model.summary()


# In[6]:


sig3FullScaled = scaler.transform(sig3Full)
bkg3FullScaled = scaler.transform(bkg3Full)

sig3Predict = loaded_model.predict(sig3FullScaled)
bkg3Predict = loaded_model.predict(bkg3FullScaled)

print(sig3FullScaled.shape)
print(bkg3FullScaled.shape)
print(sig3Predict.shape)
print(bkg3Predict.shape)


# In[12]:


# Weights for normalisation to Luminosity  

# N_bkg = 1937.09 # CR2                                                                                                                                                                                 
# N_bkg = 3646.28 # SR1                                                                                                                                                                                 
# N_bkg = 569.73 # SR2                                                                                                                                                                                  
N_bkg = 22.79 # SR3                                                                                                                                                                                     
# N_sig = 31.2*271/(5*pow(10,6)) # CR2                                                                                                                                                                  
# N_sig = 31.2*28976/(5*pow(10,6)) # SR1                                                                                                                                                                
# N_sig = 31.2*19077/(5*pow(10,6)) # SR2                                                                                                                                                                
N_sig = 31.2*11619/(5*pow(10,6)) # SR3  

w_sig = N_sig/sig3Predict.shape[0]
w_bkg = N_bkg/bkg3Predict.shape[0]
#w_sig=0.35
#w_bkg=1
print(N_sig,w_sig,N_bkg,w_bkg)


# In[28]:


# Discriminator shape
plt.clf()
plt.yscale('log')
plt.hist(np.array(sig3Predict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='signal SR3')
plt.hist(np.array(bkg3Predict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='background SR3')
plt.xlabel("Sig_Prob")
plt.ylabel("#Events (Area scaled to 1)")
plt.legend()
plt.savefig("Discriminator_SR3.pdf")
print("Discriminator plotted!!!")

plt.clf()
plt.yscale('log')
plt.hist(np.array(sig3Predict)[:,0], 
         bins=20, 
         range=(0,1), 
         density=False, 
         weights=(N_sig/sig3Predict.shape[0])*np.ones(sig3Predict.shape[0]), 
         color=None, 
         histtype='step', 
         label='signal SR3')
plt.hist(np.array(bkg3Predict)[:,0], 
         bins=20, 
         range=(0,1), 
         density=False, 
         weights=(N_bkg/bkg3Predict.shape[0])*np.ones(bkg3Predict.shape[0]), 
         color=None, 
         histtype='step', 
         label='background SR3')
plt.xlabel("Sig_Prob")
plt.ylabel("#Events (Area scaled to L=$2.6fb^{-1}$)")
plt.legend()
plt.savefig("WeightedDiscriminator_SR3.pdf")
print("Weighted Discriminator plotted!!!")


# In[22]:


# Significance Plot

sigProb = np.arange(0, 1, 0.001) 
sig3PredictSigProb = np.array(sig3Predict)[:,0]
bkg3PredictSigProb = np.array(bkg3Predict)[:,0]

SigWithCut = []

for x in sigProb: 
    sig3PredictClass = sig3PredictSigProb>x
    bkg3PredictClass = bkg3PredictSigProb>x
    nSig3TP = sig3PredictClass.sum()
    nSig3FP = bkg3PredictClass.sum()
    
    nSig3AfterCut = nSig3TP*w_sig
    nBkg3AfterCut = nSig3FP*w_bkg
    
    SigWithCut.append(nSig3AfterCut/mt.sqrt(nSig3AfterCut+nBkg3AfterCut))
    
    print(x, nSig3TP, nSig3FP, nSig3AfterCut, nBkg3AfterCut, nSig3AfterCut/mt.sqrt(nSig3AfterCut+nBkg3AfterCut))


# In[20]:


print(np.where(SigWithCut==0.13656728475323998)[0])


# In[18]:


# Significance with Discriminant

plt.clf()
plt.plot(sigProb,SigWithCut, label="NN")
plt.xlabel("Discriminant")
plt.ylabel("Significance")
plt.legend()
plt.savefig("significance_SR3.pdf")


# In[27]:


# Significance with Luminosity

disc = 0.971
lumisName = np.array([2015, "Im1", "Im2", "Run2", "Run2+Run3", "Sig1", "Sig2"])
lumis = np.array([2.6, 10, 30, 160, 300, 1000, 3000])
sigLumi = []
for lumi in lumis:
    
    sig3PredictClass = sig3PredictSigProb>disc
    bkg3PredictClass = bkg3PredictSigProb>disc
    nSig3TP = sig3PredictClass.sum()
    nSig3FP = bkg3PredictClass.sum()
    
    nSig3AfterCut = nSig3TP*w_sig*lumi/2.6
    nBkg3AfterCut = nSig3FP*w_bkg*lumi/2.6

    sigLumi.append(nSig3AfterCut/mt.sqrt(nSig3AfterCut+nBkg3AfterCut))
    
plt.clf()
plt.plot(lumisName,sigLumi, label="NN, disc=0.971")
plt.xlabel("Luminosity")
plt.ylabel("Significance")
plt.legend()
plt.savefig("significance_Lumi_SR3.pdf")


# In[ ]:




