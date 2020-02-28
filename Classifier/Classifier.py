import tensorflow as tf
from tensorflow.python.keras import models as m
from tensorflow.python.keras import layers as l

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.externals import joblib

from ROOT import TFile, TTree, TChain

print("All classes initialized successfully!!!")

sigChan = TChain("varTree")
sigChan.Add("signal.root")
bkgChan = TChain("varTree")
bkgChan.Add("background.root")
print("Data read from the trees. Printing out the contents.")

sigChan.Print()
bkgChan.Print()

sigSampleSize = sigChan.GetEntries()
bkgSampleSize = bkgChan.GetEntries()
sampleSize = sigSampleSize if sigSampleSize<bkgSampleSize else bkgSampleSize

sigFull = sigChan.AsMatrix()
print(sigFull[0])
print(sigFull.shape)
bkgFull = bkgChan.AsMatrix()
print(bkgFull.shape)

if sampleSize==bkgSampleSize:

    bkg = bkgFull
    sigIndx = np.random.choice(range(sigSampleSize), size=sampleSize, replace=False)
    sig = sigFull[sigIndx]

else:

    sig = sigFull
    bkgIndx = np.random.choice(range(bkgSampleSize), size=sampleSize, replace=False)
    bkg = bkgFull[bkgIndx]

x = np.concatenate((sig,bkg))
y = np.matrix([[1,0]]*sig.shape[0]+[[0,1]]*bkg.shape[0])
print("Feature Space: ",x.shape)

print((x.transpose())[0])

# Visualise the data
#plt.hist((x.transpose())[0])
#plt.savefig("dEtaLL.png")
#print("dEtaLL printed!")

# Randomize the data
arr = np.arange(x.shape[0])
np.random.shuffle(arr)
x = x[arr,:]
y = y[arr,:]

# Scale the input data
scaler = StandardScaler()
print(scaler.fit(x))
print(scaler.mean_)
x = scaler.transform(x)
joblib.dump(scaler, "scaler.save")

# Split the data to training and testing samples
x_train = x[0:int(0.8*x.shape[0])][:]
x_test = x[int(0.8*x.shape[0]):][:]
y_train = y[0:int(0.8*y.shape[0])][:]
y_test = y[int(0.8*y.shape[0]):][:]
print("Pre-modification of data successfull!!!")
print("Shape of training sample: ",x_train.shape)
print("Shape of testing sample: ",x_test.shape)

# Build the model
per = m.Sequential()
per.add(l.Dense(18, input_dim=x_train.shape[1], activation='relu'))
per.add(l.Dropout(rate=0.2))
per.add(l.Dense(18, activation='relu'))
per.add(l.Dropout(rate=0.2))
per.add(l.Dense(2, activation='softmax'))
print("Model building complete!!!")

per.compile(optimizer='adam',
            loss='categorical_crossentropy',
            metrics=['accuracy'])

history = per.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=1000, batch_size=int(0.1*x_train.shape[0]))
per.evaluate(x_test, y_test, verbose=2)

print(history.history.keys())

# Training Curve
plt.clf()
plt.plot(history.history['loss'], label='training')
plt.plot(history.history['val_loss'], label='test')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.legend()
plt.savefig("loss.pdf")

plt.clf()
plt.plot(history.history['accuracy'], label='training')
plt.plot(history.history['val_accuracy'], label='test')
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.legend()
plt.savefig("accuracy.pdf")

print("Predicting for signal samples.")
sigT = scaler.transform(sigFull)
sigLabel = np.matrix([[1,0]]*sigT.shape[0])
per.evaluate(sigT, sigLabel, verbose=2)

print("Predicting for background samples.")
bkgT = scaler.transform(bkgFull)
bkgLabel = np.matrix([[1,0]]*bkgT.shape[0])
per.evaluate(bkgT, bkgLabel, verbose=2)

# Save the model
per.save("simplePer.h5")
