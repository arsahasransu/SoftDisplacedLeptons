import ROOT as rt

def sigSig(rootFileName):
    signal_file = rt.TFile("../"+rootFileName+"_Disc.root","READ")
    
    print("====================================")
    print(rootFileName)
    print("====================================")
    
    normBGSR1 = 4122.8
    normBGSR2 = 644.2
    normBGSR3 = 24.479
    
    SR1_hist = signal_file.Get("SR1")
    SR2_hist = signal_file.Get("SR2")
    SR3_hist = signal_file.Get("SR3")
    bkg_hist = signal_file.Get("background")
    
    discCut = 0.95
    minBin = SR1_hist.FindBin(discCut)
    
    sr1_AC = SR1_hist.Integral(minBin, 101)
    sr2_AC = SR2_hist.Integral(minBin, 101)
    sr3_AC = SR3_hist.Integral(minBin, 101)
    bkg_SR1_AC = bkg_hist.Integral(minBin, 101)*normBGSR1
    bkg_SR2_AC = bkg_hist.Integral(minBin, 101)*normBGSR2
    bkg_SR3_AC = bkg_hist.Integral(minBin, 101)*normBGSR3

    lumivals = [2.6, 36, 140, 300]
    
    for lumi in lumivals:
        
        sigEvntSR1 = sr1_AC*lumi/2.6
        sigEvntSR2 = sr2_AC*lumi/2.6
        sigEvntSR3 = sr3_AC*lumi/2.6
        bkgEvntSR1 = bkg_SR1_AC*lumi/2.6
        bkgEvntSR2 = bkg_SR2_AC*lumi/2.6
        bkgEvntSR3 = bkg_SR3_AC*lumi/2.6
        
        signalSigSR1 = sigEvntSR1/rt.TMath.Sqrt(sigEvntSR1+bkgEvntSR1)
        signalSigSR2 = sigEvntSR2/rt.TMath.Sqrt(sigEvntSR2+bkgEvntSR2)
        signalSigSR3 = sigEvntSR3/rt.TMath.Sqrt(sigEvntSR3+bkgEvntSR3)
        
        print(round(sigEvntSR1,4),
              round(bkgEvntSR1,4),
              round(signalSigSR1,4),
              " SR1 || ",
              round(sigEvntSR2,4),
              round(bkgEvntSR2,4),
              round(signalSigSR2,4),
              " SR2 || ",
              round(sigEvntSR3,4),
              round(bkgEvntSR3,4),
              round(signalSigSR3,4),
              " SR3 || ",
              "Luminosity: ",lumi)
        
    signal_file.Close()

sigSig("BP_200_20_DM")
sigSig("BP_324_20_DM")
sigSig("BP_200_20_02")
sigSig("BP_200_20_2")
sigSig("BP_200_20_20")
sigSig("BP_200_20_200")
sigSig("BP_200_40_20")

