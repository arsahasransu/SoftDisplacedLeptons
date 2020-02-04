import os, sys
import ROOT as rt
import math
from itertools import permutations

# small macro that just takes all branches in a flat tree and calculates the ranking

noplots=False
extracut=" ";
sigtree=rt.TChain("varTree","varTree")
bgtree=rt.TChain("varTree","varTree")
sigtree.Add("signal.root")
bgtree.Add("background.root")
print("now filled trees:",sigtree.GetEntries(),bgtree.GetEntries())
canv = rt.TCanvas("canv","canv")
listofbranches = sigtree.GetListOfBranches()
ranking=[[100,"variable"]]
for branchname in listofbranches:
#    print(branchname)
    workname=branchname.GetName()
#    print(workname)
    min=0
    max=1
    sigtree.Draw(workname+">>htemp"+workname)
    htemp = rt.gDirectory.Get("htemp"+workname)
    min =htemp.GetXaxis().GetXmin()
    max = htemp.GetXaxis().GetXmax()
    print(workname, htemp.GetXaxis().GetXmin() , htemp.GetXaxis().GetXmax(), min, max )
    if min > htemp.GetXaxis().GetXmin() :
        min = htemp.GetXaxis().GetXmin()
    if max < htemp.GetXaxis().GetXmax() :
        max = htemp.GetXaxis().GetXmax()
    print("checking out axis ranges etc...",min,max)
    histsig = rt.TH1D("histsig"+workname,"",20,min,max)
    histsig.Sumw2()
    histsig.SetXTitle(workname)
    histbg= rt.TH1D("histbg"+workname,"",20,min,max)
    histbg.Sumw2()
    histbg.SetXTitle(workname)
    drawcommand = workname + ">>"+histsig.GetName()+"(20,"+str(min)+","+str(max)+")"
    sigtree.Draw(drawcommand)
    print(drawcommand)
    drawcommand = workname + ">>"+histbg.GetName()+"(20,"+str(min)+","+str(max)+")"
    print(drawcommand)
    bgtree.Draw(drawcommand)
    histbg=rt.gDirectory.Get("histbg"+workname)
    histsig=rt.gDirectory.Get("histsig"+workname)
#    print("now histos are filled")
    if histbg.GetSum()*histsig.GetSum()==0 :
        print("histograms of ",workname, " " ,histsig.GetName()," and ", histbg.GetName()," with ", histsig.GetEntries()," and ",histbg.GetEntries()," events, and integrals ",histsig.GetSum()," ",histbg.GetSum()," exiting")
        continue
    
#    print("histograms ",histsig.GetName()," and ", histbg.GetName()," with ", histsig.GetEntries()," and ",histbg.GetEntries()," events")
  
    print(workname," comparison, probability overlap:")
    histsig.Scale(1./histsig.GetSum())
    histbg.Scale(1./histbg.GetSum())
    probhist1=histsig.Clone("probhist1"+workname)
    probhist2=histbg.Clone("probhist2"+workname)
    probhist1.Multiply(probhist2)
    print("overlap probability: ",probhist1.GetSum())
    rankworker=[round(probhist1.GetSum(),4),workname]
    ranking.append(rankworker)
    if noplots == False :
          canv.cd()
          histsig.SetLineColor(rt.kRed)
          histbg.SetLineColor(rt.kAzure)
          if histsig.GetMaximum()< histbg.GetMaximum():
              histsig.SetMaximum(1.05*histbg.GetMaximum())
          histsig.Draw()
          histbg.Draw("same")
          
          canv.Update()
          canv.BuildLegend()
          canv.Update()
#          canv.Print("overviewplot_"+workname+".root")
#          canv.Print("overviewplot_"+workname+".pdf")
          canv.Print("overviewplot_"+workname+".png")
    
print("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-")
print("variable ranking:")
ranking.sort()
for rank in ranking:
    if rank[0]<100 :
        print(rank[1],rank[0])

# now do this in 2D:
print("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-")

print("2D correlations (calculating)")
correlationranking=[[]]

perm = permutations(listofbranches,2)
for p1 in list(perm):
    workname1=p1[0].GetName()
    workname2=p1[1].GetName()
    min1=sigtree.GetMinimum(workname1)
    max1=sigtree.GetMaximum(workname1)
    min2=sigtree.GetMinimum(workname2)
    max2=sigtree.GetMaximum(workname2)
    workhist = rt.TH2D("workhist"+workname1+workname2,"",40,min1,max1,40,min2,min2)
    workhistbg=workhist.Clone("workhistbg"+workname1+workname2)
    sigtree.Draw(workname1+":"+workname2+">>+workhist"+workname1+workname2)
    bgtree.Draw(workname1+":"+workname2+">>+workhistbg"+workname1+workname2)
   
    print(workname1,workname2," signal ",round(workhist.GetCorrelationFactor(),4)," background ",round(workhistbg.GetCorrelationFactor(),4))


input("Press Enter to continue...")



