import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.externals import joblib
from sklearn.metrics import confusion_matrix, fbeta_score
from scikitplot.metrics import plot_roc, plot_confusion_matrix

import seaborn as sn
import pandas as pd

from ROOT import TFile, TTree, TChain

print("All classes initialized successfully!!!")

sigChan = TChain("varTree")
sigChan.Add("sigVar.root")
bkgChan = TChain("varTree")
bkgChan.Add("bkgVar.root")
print("Data read from the trees.")

sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()
sampleSize = sigSampleSize if sigSampleSize<bkgSampleSize else bkgSampleSize

sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()
sampleSize = sigSampleSize if sigSampleSize<bkgSampleSize else bkgSampleSize

# Loading the input model
sigFull = sigChan.AsMatrix()
bkgFull = bkgChan.AsMatrix()

if sampleSize==bkgSampleSize:

    bkg = bkgFull
    sigIndx = np.random.choice(range(sigSampleSize), size=sampleSize, replace=False)
    sig = sigFull[sigIndx]

else:

    sig = sigFull
    bkgIndx = np.random.choice(range(bkgSampleSize), size=sampleSize, replace=False)
    bkg = bkgFull[bkgIndx]

# Load the input data scaler
scaler = joblib.load("scaler.save")

# Load the model
loaded_model = m.load_model("simplePer.h5")
loaded_model.summary()

# Prep the data for model evaluation
sigT = scaler.transform(sig)
sigLabel = np.matrix([1]*sigT.shape[0])
sigLabel = np.transpose(sigLabel)

bkgT = scaler.transform(bkg)
bkgLabel = np.matrix([0]*bkgT.shape[0])
bkgLabel = np.transpose(bkgLabel)

# Predict on the data
sig_pred = loaded_model.predict(sigT)
bkg_pred = loaded_model.predict(bkgT)

# Select for a single output discriminator variable
sig_pred_one = np.transpose(sig_pred)[0]
bkg_pred_one = np.transpose(bkg_pred)[0]

# Combine signal and background
pred = np.concatenate((sig_pred_one, bkg_pred_one))
labl = np.concatenate((sigLabel, bkgLabel))
print(pred.shape)
print(labl.shape)

# Discriminator shape
plt.clf()
plt.hist(sig_pred_one, bins=100, range=(0,1), density=True, color=None, histtype='step', label='sig')
plt.hist(bkg_pred_one, bins=100, range=(0,1), density=True, color=None, histtype='step', label='bkg')
plt.xlabel("Sig_Prob")
plt.ylabel("#Events (Area scaled to 1)")
plt.legend()
plt.savefig("Discriminator.pdf")

# Confusion matrix
threshold = 0.5
pred_class = pred>threshold
plt.clf()
plt.figure(figsize = (10,7))
plot_confusion_matrix(labl, pred_class)
plt.savefig("ConfusionMatrix.pdf")

# F1 Score

nofvalues = 100
thresh = np.empty(nofvalues)
f1score = np.empty(nofvalues)

plt.clf()
for beta in [0.1,0.25,0.5,1,1.25,1.5,2,5,10]:
    for i in range(100):
        thresh[i] = i*1.0/100
        pred_class = pred>thresh[i]
        f1score[i] = fbeta_score(labl, pred_class, beta)

    plt.plot(thresh,f1score,label=beta)

plt.xlabel("threshold")
plt.ylabel("F1 Score")
plt.legend()
plt.savefig('f1score.pdf')

# ROC Curve

threshold = 0.5
pred_class = pred>threshold
labl_class = labl>threshold
plt.clf()
fig, ax = plt.subplots()
print(labl_class.shape)
print(pred_class.shape)
plot_roc(labl_class, pred_class, ax=ax)
plt.savefig('roc.pdf')

###### END OF STUDY #######
