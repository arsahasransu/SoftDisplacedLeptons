#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l


# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import joblib


# In[3]:


from ROOT import TFile, TTree, TChain
print("All classes initialized succesfully.")


# In[4]:


sigChan = TChain("varTree")
sigChan.Add("signal.root")
bkgChan = TChain("varTree")
bkgChan.Add("background.root")
print("Data read from the trees. Printing out the contents.")


# In[5]:


sigChan.Print()
bkgChan.Print()


# In[6]:


# Read input data from root files
sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()
print(sigSampleSize)
print(bkgSampleSize)


# In[7]:


# Convert the input data to matrices
sigFull = sigChan.AsMatrix()
bkgFull = bkgChan.AsMatrix()
print(sigFull[0])
print(sigFull.shape)
print(bkgFull.shape)


# In[9]:


# Choose 20% of the data for testing and 80% of the data for training
sigTrain = sigFull[0:int(0.8*sigFull.shape[0])][:]
bkgTrain = bkgFull[0:int(0.8*bkgFull.shape[0])][:]
sigTest = sigFull[int(0.8*sigFull.shape[0]):][:]
bkgTest = bkgFull[int(0.8*bkgFull.shape[0]):][:]
print(sigTrain.shape)
print(bkgTrain.shape)
print(sigTest.shape)
print(bkgTest.shape)


# In[10]:


# Set a scaler for input features
scaler = StandardScaler()
scaler_input = np.concatenate((sigTrain,bkgTrain))
print(scaler.fit(scaler_input))
print(scaler.mean_)
joblib.dump(scaler, "scaler.save")


# In[11]:


# Build the model
per = m.Sequential()
per.add(l.Dense(18, input_dim=sigTrain.shape[1], activation='relu'))
per.add(l.Dropout(rate=0.2))
per.add(l.Dense(18, activation='relu'))
per.add(l.Dropout(rate=0.2))
per.add(l.Dense(2, activation='softmax'))
print("Model building complete!!!")

per.compile(optimizer='adam',
            loss='categorical_crossentropy',
            metrics=['accuracy'])


# In[12]:


# Loop to change the training sample every time
# by randomly choosing from the avalaible sample space.
# Make sure to run atleast so that each and every event has been used once.

nShuffleRun = 10
trainingSampleSize = sigTrain.shape[0] if sigTrain.shape[0]<bkgTrain.shape[0] else bkgTrain.shape[0]
testSampleSize = sigTest.shape[0] if sigTest.shape[0]<bkgTest.shape[0] else bkgTest.shape[0]

print(trainingSampleSize)
print(testSampleSize)


# In[13]:


nEpochs = 100

lossTrain = []
lossTest = []
accTrain = []
accTest = []

for iterTrain in np.arange(nShuffleRun):
    
    sigTrainSampleRange = np.arange(sigTrain.shape[0])
    sigTrainSampleIndex = np.random.choice(sigTrainSampleRange, trainingSampleSize, replace=False)
    sigTrainChosen = sigTrain[sigTrainSampleIndex,:]
    bkgTrainSampleRange = np.arange(bkgTrain.shape[0])
    bkgTrainSampleIndex = np.random.choice(bkgTrainSampleRange, trainingSampleSize, replace=False)
    bkgTrainChosen = bkgTrain[bkgTrainSampleIndex,:]
    
    sigTestSampleRange = np.arange(sigTest.shape[0])
    sigTestSampleIndex = np.random.choice(sigTestSampleRange, testSampleSize, replace=False)
    sigTestChosen = sigTest[sigTestSampleIndex,:]
    bkgTestSampleRange = np.arange(bkgTest.shape[0])
    bkgTestSampleIndex = np.random.choice(bkgTestSampleRange, testSampleSize, replace=False)
    bkgTestChosen = bkgTest[bkgTestSampleIndex,:]
    
    #print(sigTrainChosen.shape)
    #print(bkgTrainChosen.shape)
    
    # Concatenate the signal and background with proper labels
    x = np.concatenate((sigTrainChosen,bkgTrainChosen))
    y = np.matrix([[1,0]]*sigTrainChosen.shape[0]+[[0,1]]*bkgTrainChosen.shape[0])
    x_test = np.concatenate((sigTestChosen,bkgTestChosen))
    y_test = np.matrix([[1,0]]*sigTestChosen.shape[0]+[[0,1]]*bkgTestChosen.shape[0])
    print("Feature Space: ",x.shape)
    
    # Randomize the training samples
    arr = np.arange(x.shape[0])
    np.random.shuffle(arr)
    x = x[arr,:]
    y = y[arr,:]
    
    # Scale the input features
    x = scaler.transform(x)
    x_test = scaler.transform(x_test)
    
    history = per.fit(x, 
                      y, 
                      validation_data=(x_test, y_test), 
                      epochs=nEpochs, 
                      batch_size=int(0.1*x.shape[0]))
    
    lossTrain.extend(history.history['loss'])
    lossTest.extend(history.history['val_loss'])
    accTrain.extend(history.history['accuracy'])
    accTest.extend(history.history['val_accuracy'])


# In[14]:


# Training Curve
plt.clf()
plt.plot(lossTrain, label='training')
plt.plot(lossTest, label='test')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.savefig("loss.pdf")


# In[15]:


plt.clf()
plt.plot(accTrain, label='training')
plt.plot(accTest, label='test')
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.legend()
plt.savefig("accuracy.pdf")


# In[ ]:




