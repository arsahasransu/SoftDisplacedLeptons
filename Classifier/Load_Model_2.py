#!/usr/bin/env python
# coding: utf-8

# In[44]:


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

import plotBeautifier as pB

print("All classes initialized successfully!!!")


# In[46]:


pB.trial_func("AR")


# In[2]:


sigChan = TChain("varTree")
sigChan.Add("signal.root")
bkgChan = TChain("varTree")
bkgChan.Add("background.root")
print("Data read from the trees. Printing out the contents.")


# In[3]:


sigChan.Print()
bkgChan.Print()


# In[4]:


# Read input data from root files
sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()

# Convert the input data to matrices
sigFull = sigChan.AsMatrix()
bkgFull = bkgChan.AsMatrix()

# Choose 20% of the data for testing and 80% of the data for training
sigTrain = sigFull[0:int(0.8*sigFull.shape[0])][:]
bkgTrain = bkgFull[0:int(0.8*bkgFull.shape[0])][:]
sigTest = sigFull[int(0.8*sigFull.shape[0]):][:]
bkgTest = bkgFull[int(0.8*bkgFull.shape[0]):][:]
print(sigTrain.shape)
print(bkgTrain.shape)
print(sigTest.shape)
print(bkgTest.shape)


# In[5]:


# Load the input data scaler
scaler = joblib.load("../scaler.save")

# Load the model
loaded_model = m.load_model("../simplePer.h5")
loaded_model.summary()


# In[6]:


xTrain = np.concatenate((sigTrain,bkgTrain))
xTest = np.concatenate((sigTest,bkgTest))
yTrain = np.matrix([[1,0]]*sigTrain.shape[0]+[[0,1]]*bkgTrain.shape[0])
yTest = np.matrix([[1,0]]*sigTest.shape[0]+[[0,1]]*bkgTest.shape[0])

print(xTrain.shape)
print(xTest.shape)
print(yTrain.shape)
print(yTest.shape)


# In[7]:


# Randomize the training and testing samples
arr = np.arange(xTrain.shape[0])
np.random.shuffle(arr)
xTrain = xTrain[arr,:]
yTrain = yTrain[arr,:]

arr = np.arange(xTest.shape[0])
np.random.shuffle(arr)
xTest = xTest[arr,:]
yTest = yTest[arr,:]


# In[8]:


# Predict on the samples
sigTrainScaled = scaler.transform(sigTrain)
bkgTrainScaled = scaler.transform(bkgTrain)
sigTestScaled = scaler.transform(sigTest)
bkgTestScaled = scaler.transform(bkgTest)
xTestScaled = scaler.transform(xTest)
sigTrainPredict = loaded_model.predict(sigTrainScaled)
bkgTrainPredict = loaded_model.predict(bkgTrainScaled)
sigTestPredict = loaded_model.predict(sigTestScaled)
bkgTestPredict = loaded_model.predict(bkgTestScaled)
xTestPredict = loaded_model.predict(xTestScaled)


# In[9]:


print(sigTrainPredict[0:5])


# In[10]:


print(np.array(sigTrainPredict)[0:5,0])


# In[11]:


# Discriminator shape
plt.clf()
plt.yscale('log')
plt.hist(np.array(sigTrainPredict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='sigTrain')
plt.hist(np.array(bkgTrainPredict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='bkgTrain')
plt.hist(np.array(sigTestPredict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='sigTest')
plt.hist(np.array(bkgTestPredict)[:,0], bins=20, range=(0,1), density=True, color=None, histtype='step', label='bkgTest')
plt.xlabel("Sig_Prob")
plt.ylabel("#Events (Area scaled to 1)")
plt.legend()
plt.savefig("Discriminator.pdf")
print("Discriminator plotted!!!")


# In[33]:


# Discriminator Shape in ROOT plotting
c1 = TCanvas()
nBins = 10
sigTrainHisto = TH1F("","",nBins,0,1)
bkgTrainHisto = TH1F("","",nBins,0,1)
sigTestHisto = TH1F("","",nBins,0,1)
bkgTestHisto = TH1F("","",nBins,0,1)
sigTrainHisto.FillN(sigTrainPredict.shape[0],(sigTrainPredict[:,0]).astype(float),np.ones(sigTrainPredict.shape[0]))
bkgTrainHisto.FillN(bkgTrainPredict.shape[0],(bkgTrainPredict[:,0]).astype(float),np.ones(bkgTrainPredict.shape[0]))
sigTestHisto.FillN(sigTestPredict.shape[0],(sigTestPredict[:,0]).astype(float),np.ones(sigTestPredict.shape[0]))
bkgTestHisto.FillN(bkgTestPredict.shape[0],(bkgTestPredict[:,0]).astype(float),np.ones(bkgTestPredict.shape[0]))
sigTrainHisto.DrawNormalized("SAME E")
bkgTrainHisto.DrawNormalized("SAME E")
sigTestHisto.DrawNormalized("SAME E")
bkgTestHisto.DrawNormalized("SAME E")
c1.SetLogy()
c1.SaveAs("Discriminator.pdf")


# In[39]:


# Using plotBeautifier

histList = [sigTrainHisto, bkgTrainHisto, sigTestHisto, bkgTestHisto]
labelList = ["sigTrainHisto", "bkgTrainHisto", "sigTestHisto", "bkgTestHisto"]
xAxisTitle = "Sig_Prob"
yAxisTitle = "# Events (scaled to 1)"
outPlotName = "Discriminator_Beautified.pdf"
pB.plotBeautifier(histList, labelList, xAxisTitle, yAxisTitle, outPlotName)


# In[27]:


type((sigTrainPredict[:,0]).astype(float))


# In[12]:


# Confusion Matrix
plt.clf()
plt.figure(figsize = (10,7))
plot_confusion_matrix(np.array(yTest)[:,0], np.array(xTestPredict)[:,0]>0.5)
plt.savefig("ConfusionMatrix.pdf")
print("Confusion Matrix printed!!!")


# In[15]:


# For plots based on discriminator cut

sigProb = np.arange(0, 1, 0.01) 
xTestPredictSigProb = np.array(xTestPredict)[:,0] 
yTestClass = np.array(yTest)[:,0]

tpr = [] 
fpr = []

for x in sigProb: 
    xTestPredictClass = xTestPredictSigProb>x 
    tp = (xTestPredictClass&yTestClass).sum() 
    fp = (xTestPredictClass&(1-yTestClass)).sum() 
    fn = ((1-xTestPredictClass)&yTestClass).sum() 
    tn = ((1-xTestPredictClass)&(1-yTestClass)).sum()

    tpr.append(tp/(tp+fn))
    fpr.append(fp/(fp+tn))


# In[16]:


# ROC Curve

plt.clf()
plt.plot(fpr,tpr, label="NN")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend()
plt.savefig("roc.pdf")


# In[17]:


# Plots for scaled input variables above

for varCtr in np.arange(sigTest.shape[1]):
    plt.clf()
    plt.yscale('log')
    plt.hist(np.transpose(sigTest)[varCtr], bins=20, range=(-5.3,5.3), density=True, color=None, histtype='step', label='sigTest')
    plt.hist(np.transpose(bkgTest)[varCtr], bins=20, range=(-5.3,5.3), density=True, color=None, histtype='step', label='bkgTest')
    plt.show()


# In[ ]:




