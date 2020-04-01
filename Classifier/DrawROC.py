#!/usr/bin/env python
# coding: utf-8

# In[4]:


import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt


# In[1]:


# From histogram as lists

def from_histo(discSig=[], discBkg=[], name=""):
    
    sigEff = []
    bkgEff = []
    
    maxVal = max(discSig) if(max(discSig)>max(discBkg)) else max(discBkg)
    minVal = min(discSig) if(min(discSig)<min(discBkg)) else min(discBkg)
    
    discStep = np.arange(minVal, maxVal, (maxVal-minVal)*1.0/1000.0)
    
    for discCut in discStep:
        
        tp = (discSig>discCut).sum()*1.0
        fp = (discBkg>discCut).sum()*1.0
        tn = (discBkg<discCut).sum()*1.0
        fn = (discSig<discCut).sum()*1.0
        
        sigEff.append(tp/(tp+fn))
        bkgEff.append(fp/(fp+tn))
        
    roc = plt.plot(bkgEff,sigEff,label="roc_"+name)
    
    return roc


# In[ ]:


# Overload for numpy arrays

