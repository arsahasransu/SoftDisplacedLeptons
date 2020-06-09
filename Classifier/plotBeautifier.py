#!/usr/bin/env python
# coding: utf-8

# In[2]:


import ROOT as rt


# In[3]:


def plotBeautifier(histList, 
                   labelList, 
                   xAxisTitle, 
                   yAxisTitle, 
                   outPlotName):
    
    # Create a canvas 
    c1 = rt.TCanvas("c1", "c1", 10, 32, 782, 552)
    print("Created a canvas")
    
    # Add cosmetics, underflow bin and overflow bin
    ctr=0;
    for hist in histList:
        hist.SetLineColor(ctr+1)
        ctr = ctr+1
        hist.SetLineWidth(3)
        hist.GetXaxis().SetRange(0, hist.GetNbinsX()+1)
        print(ctr)
        
    histList[0].GetXaxis().SetTitle(xAxisTitle)
    histList[0].GetYaxis().SetTitle(yAxisTitle)
    
    c1.SetLogy()
    rt.gStyle.SetOptStat(0)
    
    for hist in histList:
        hist.DrawNormalized("hist SAME E")
        
    legc1 = rt.TLegend(0.5, 0.9, 0.89, 1.0, "", "brNDC")
    ctr=0
    for hist in histList:
        legc1.AddEntry(hist, labelList[ctr], "l")
        ctr = ctr+1
    legc1.SetTextSize(0.03)
    legc1.SetBorderSize(0)
    legc1.Draw()
    
    nx = histList[0].GetNbinsX()+1
    bw1 = histList[0].GetBinWidth(0)
    bw2 = histList[0].GetBinWidth(nx)
    x1 = histList[0].GetBinLowEdge(0)+bw1
    x2 = histList[0].GetBinLowEdge(nx)+bw2
    y1 = histList[0].GetBinContent(0)/histList[0].Integral()
    y2 = histList[0].GetBinContent(nx)/histList[0].Integral()
    for hist in histList:
        y1 = y1 if y1>(hist.GetBinContent(0)/hist.Integral()) else hist.GetBinContent(0)/hist.Integral()
        
        y2 = y2 if y2>(hist.GetBinContent(nx)/hist.Integral()) else hist.GetBinContent(nx)/hist.Integral()
        
    y1 = y1*1.1
    y2 = y2*1.1
    
    tUFlw = rt.TText(x1-0.5*bw1, y1, "Underflow")
    tUFlw.SetTextAngle(90)
    tUFlw.SetTextAlign(12)
    tUFlw.SetTextSize(0.03)
    tUFlw.Draw()
    
    tOFlw = rt.TText(x2-0.5*bw2, y2, "Overflow")
    tOFlw.SetTextAngle(90)
    tOFlw.SetTextAlign(12)
    tOFlw.SetTextSize(0.03)
    tOFlw.Draw()
    
    c1.SaveAs(outPlotName+".C")
    
    return


# In[4]:


def trial_func(name):
    
    print("Trial works fine. Hello,"+name)


# In[ ]:




