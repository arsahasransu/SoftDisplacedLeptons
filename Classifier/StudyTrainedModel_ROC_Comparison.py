#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l

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


# In[2]:


import plotBeautifier as pB


# In[3]:


pB.trial_func("AR")


# In[4]:


sigChan = TChain("varTree")
sigChan.Add("signal.root")
bkgChan = TChain("varTree")
bkgChan.Add("background.root")
print("Data read from the trees. Printing out the contents.")


# In[33]:


brNameList = []
for br in sigChan.GetListOfBranches():
    brNameList.append(br.GetName())
    
print(brNameList)
print(len(brNameList))


# In[6]:


# Read input data from root files
sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()

# Convert the input data to matrices
sigFull = sigChan.AsMatrix()
bkgFull = bkgChan.AsMatrix()

print(sigFull.shape)
print(bkgFull.shape)


# In[7]:


# Load the input data scaler
scaler = joblib.load("../scaler.save")

# Load the model
loaded_model = m.load_model("../simplePer.h5")
loaded_model.summary()


# In[9]:


# Predict on the samples

sigFullScaled = scaler.transform(sigFull)
bkgFullScaled = scaler.transform(bkgFull)

sigFullPredict = loaded_model.predict(sigFullScaled)
bkgFullPredict = loaded_model.predict(bkgFullScaled)

print(sigFullScaled.shape)


# In[40]:


for i in np.arange(0,9,1):
    plt.clf()
    plt.yscale('log')
    #plt.xscale('log')
    plt.hist(np.array(np.transpose(sigFull)[i]), bins=20, density=True, color=None, histtype='step', label='sig BP(304,324) Unscaled')
    plt.hist(np.array(np.transpose(bkgFull)[i]), bins=20, density=True, color=None, histtype='step', label='bkg Unscaled')
    plt.hist(np.array(np.transpose(sigFullScaled)[i]), bins=20, density=True, color=None, histtype='step', label='sig BP(304,324)')
    plt.hist(np.array(np.transpose(bkgFullScaled)[i]), bins=20, density=True, color=None, histtype='step', label='bkg')
    plt.xlabel(brNameList[i])
    plt.legend()
    plt.show()
    


# In[13]:


np.transpose(sigFullScaled)[0]


# In[74]:


tpr = []
fpr = []

sigProb = np.arange(0, 1, 0.01) 
sigPredictSigProb = np.array(sigFullPredict)[:,0]
bkgPredictSigProb = np.array(bkgFullPredict)[:,0]

for x in sigProb: 
     
    tp = (sigPredictSigProb>x).sum() 
    fp = (bkgPredictSigProb>x).sum() 
    fn = (sigPredictSigProb<x).sum() 
    tn = (bkgPredictSigProb<x).sum()

    tpr.append(tp/(tp+fn))
    fpr.append(fp/(fp+tn))
    
plt.clf()
plt.plot(fpr,tpr, label="NN")

minValList = [0, 0, 0]
maxValList = [200, 5, 3.14]
stepValList = [1, 0.01, 0.01]

tpr = []
fpr = []

for i in range(3):

    for varDisc in np.arange(minValList[i], maxValList[i], stepValList[i]):
    
        tp = (np.transpose(sigFull)[i]>varDisc).sum()
        fp = (np.transpose(bkgFull)[i]>varDisc).sum()
        tn = (np.transpose(bkgFull)[i]<varDisc).sum()
        fn = (np.transpose(sigFull)[i]<varDisc).sum()
    
        tpr.append(tp/(tp+fn))
        fpr.append(fp/(fp+tn))

    plt.plot(fpr,tpr, label=brNameList[i].format(i=i))

plt.legend(loc='best')
plt.show()


# In[ ]:




