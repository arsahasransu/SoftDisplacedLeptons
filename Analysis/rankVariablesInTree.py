import os, sys
import ROOT as rt
import math
from itertools import combinations

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
# create directory to dump all the plots in
dirname="rankVariablesInTree_plots"
rt.gSystem.Exec("mkdir "+dirname)
for branchname in listofbranches:
    print("next variable is:",branchname)
    workname=branchname.GetName()
#    print(workname)
    min=0
    max=1
    sigtree.Draw(workname+">>htemp"+workname)
    htemp = rt.gDirectory.Get("htemp"+workname)
    min =htemp.GetXaxis().GetXmin()
    max = htemp.GetXaxis().GetXmax()
#    print(workname, htemp.GetXaxis().GetXmin() , htemp.GetXaxis().GetXmax(), min, max )
    if min > htemp.GetXaxis().GetXmin() :
        min = htemp.GetXaxis().GetXmin()
    if max < htemp.GetXaxis().GetXmax() :
        max = htemp.GetXaxis().GetXmax()
    NBINS = htemp.GetXaxis().GetNbins()
    htemp.Delete()
#    print("checking out axis ranges etc...",min,max)
    histsig = rt.TH1D("histsig"+workname,"",NBINS,min,max)
    histsig.Sumw2()
    histsig.SetXTitle(workname)
    histbg= rt.TH1D("histbg"+workname,"",NBINS,min,max)
    histbg.Sumw2()
    histbg.SetXTitle(workname)
    probhist1=histsig.Clone("probhist1"+workname)
    drawcommand = workname + ">>"+histsig.GetName()+"("+str(NBINS)+","+str(min)+","+str(max)+")"
    sigtree.Draw(drawcommand)
#    print(drawcommand)
    drawcommand = workname + ">>"+histbg.GetName()+"("+str(NBINS)+","+str(min)+","+str(max)+")"
#    print(drawcommand)
    bgtree.Draw(drawcommand)
    histbg=rt.gDirectory.Get("histbg"+workname)
    histsig=rt.gDirectory.Get("histsig"+workname)
#    print("now histos are filled")
    if histbg.GetSum()*histsig.GetSum()==0 :
        print("histograms of ",workname, " " ,histsig.GetName()," and ", histbg.GetName()," with ", histsig.GetEntries()," and ",histbg.GetEntries()," events, and integrals ",histsig.GetSum()," ",histbg.GetSum()," exiting")
        continue
    
#    print("histograms ",histsig.GetName()," and ", histbg.GetName()," with ", histsig.GetEntries()," and ",histbg.GetEntries()," events")
  
    # Normalisation to Luminosity
    # N_bkg = 1937.09 # CR2
    N_bkg = 3646.28 # SR1
    # N_bkg = 569.73 # SR2
    # N_bkg = 22.79 # SR3
    # N_sig = 6.24*28976*0.000001 # CR2
    N_sig = 6.24*28976*0.000001 # SR1
    # N_sig = 6.24*19077*0.000001 # SR2
    # N_sig =  6.24*10860*0.000001 # SR3

    print(workname," comparison, probability overlap:")
    #histsig.Scale(N_sig)
    #histbg.Scale(N_bkg)
    histsig.Scale(1./histsig.GetSum())
    histbg.Scale(1./histbg.GetSum())

    for ibin in range(NBINS+1):
        probhist1.SetBinContent(ibin,rt.TMath.Min(histsig.GetBinContent(ibin),histbg.GetBinContent(ibin)))
        probhist1.SetBinError(ibin,rt.TMath.Max(histsig.GetBinError(ibin),histbg.GetBinError(ibin)))
    print("overlap probability: ",probhist1.GetSum())
    rankworker=[round(probhist1.GetSum(),4),workname]
    ranking.append(rankworker)
    if noplots == False :
          canv.cd()
          histsig.SetLineColor(rt.kRed)
          histbg.SetLineColor(rt.kAzure)
          histbg.SetTitle("background")
          histsig.SetTitle("signal")
          histsig.SetXTitle(workname)
          if histsig.GetMaximum()< histbg.GetMaximum():
              histsig.SetMaximum(1.05*histbg.GetMaximum())
          histsig.Draw()
          histbg.Draw("same")
          
          canv.Update()
          canv.BuildLegend()
          canv.Update()
          canv.Print(dirname+"/overviewplot_"+workname+".png")
    probhist1.Delete()
    histsig.Delete()
    histbg.Delete()
    
print("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-")
print("variable ranking:")
ranking.sort()
for rank in ranking:
    if rank[0]<100 :
        print(rank[1],rank[0])

# now do this in 2D:
print("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-")

print("2D correlations (press enter to continue)")
correlationranking=[[]]

input("Press Enter to continue...")

perm = combinations(listofbranches,2)
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

    sigcor = round(workhist.GetCorrelationFactor(),4)
    bkgcor= round(workhistbg.GetCorrelationFactor(),4)
    #if abs(bkgcor)>0.3: # Added control for printing with condition on value of co-relation
    print(workname1,workname2," signal ",sigcor," background ",bkgcor)


input("Press Enter to continue...")



