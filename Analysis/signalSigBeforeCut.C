#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVectorD.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <iostream>

class Calculator {

public:
  Calculator();
  ~Calculator();
  void yieldCalc(TChain* signal, int*, int*, int*, int*);

  void d0_weight(double eled0, double muod0,
		 std::vector<double> &);

  // Declare histogram suffix
  std::vector<TString> histSuff;
};

Calculator::Calculator() {
}

Calculator::~Calculator() {
}

void Calculator::yieldCalc(TChain* signal, int*tot, int* SR1, int* SR2, int*SR3) {

  bool isSignal = true;
  
  // Declare variables to read from the trees
  int numElec, numMuon, numJet;
  double MET, MET_Eta, MET_Phi, HT;
  std::vector<int> *PID = 0;
  std::vector<double> *E = 0;
  std::vector<double> *PT = 0;
  std::vector<double> *Eta = 0;
  std::vector<double> *Phi = 0;
  std::vector<double> *Iso = 0;
  std::vector<double> *D0 = 0;
  std::vector<double> *Dz = 0;
  std::vector<double> *JetPT = 0;
  std::vector<double> *JetM = 0;
  std::vector<double> *JetEta = 0;
  std::vector<double> *JetPhi = 0;

  // Link the defined variables to the leaf on signal tree
  signal->SetBranchAddress("NumElec",&numElec);
  signal->SetBranchAddress("NumMuon",&numMuon);
  signal->SetBranchAddress("MET", &MET);
  signal->SetBranchAddress("MET_Eta", &MET_Eta);
  signal->SetBranchAddress("MET_Phi", &MET_Phi);
  signal->SetBranchAddress("PID", &PID);
  signal->SetBranchAddress("simulatedE", &E);
  signal->SetBranchAddress("PT", &PT);
  signal->SetBranchAddress("Eta", &Eta);
  signal->SetBranchAddress("Phi", &Phi);
  signal->SetBranchAddress("Iso", &Iso);
  signal->SetBranchAddress("D0", &D0);
  signal->SetBranchAddress("Dz", &Dz);
  signal->SetBranchAddress("HT", &HT);
  signal->SetBranchAddress("NumJet", &numJet);
  signal->SetBranchAddress("JetEta", &JetEta);
  signal->SetBranchAddress("JetPhi", &JetPhi);
  signal->SetBranchAddress("JetPT", &JetPT);
  signal->SetBranchAddress("JetM", &JetM);

  // Selected Events for the total yield
  int SelectedEvents = 0;
  int SelEvntSR1 = 0;
  int SelEvntSR2 = 0;
  int SelEvntSR3 = 0;

  // Open Event Loop
  for(int evtCtr=0; evtCtr<signal->GetEntries(); evtCtr++) {

    if(evtCtr%100000==0) cout<<signal->GetEntries()<<" total. Ongoing event: "<<evtCtr<<endl; 
    
    signal->GetEntry(evtCtr);
    double lenObj = PID->size();

    // Select for good leptons
    bool foundGoodLep = false;
    int firstPos=-1;
    int secondPos=-1;
    int signalRegion = -1;
    int numGoodLep = 0;
    int numGoodEl = 0;
    int numGoodMu = 0;
    std::vector<int> lepPos;

    // Loop to select detector acceptable leptons
    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;
      if(TMath::Abs(D0->at(objCtr))>100) continue; // D0 acceptance to 10 cm
      if(isSignal && TMath::Abs(PID->at(objCtr))==11 && Iso->at(objCtr)>0.12) continue;
      if(isSignal && TMath::Abs(PID->at(objCtr))==13 && Iso->at(objCtr)>0.15) continue;
	
      lepPos.push_back(objCtr);
    }

    if(lepPos.size()<2) continue; // Not enough leptons qualified for event selection

    // Loop to select the leptons and the signal region
    for(int lep1Ctr=0; lep1Ctr<lepPos.size()-1; lep1Ctr++) {
      for(int lep2Ctr=lep1Ctr+1; lep2Ctr<lepPos.size(); lep2Ctr++) {
	if(PID->at(lepPos[lep1Ctr])*PID->at(lepPos[lep2Ctr]) != -11*13) continue;

	if(isSignal && TMath::Abs(D0->at(lepPos[lep1Ctr]))>=1 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=1) { // break if event is signal region 3
	  signalRegion = 4;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	  break;
	}
	else if(isSignal && TMath::Abs(D0->at(lepPos[lep1Ctr]))>=0.5 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=0.5) {
	  if(signalRegion>=3) continue;
	  signalRegion = 3;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	}
	else if(isSignal && TMath::Abs(D0->at(lepPos[lep1Ctr]))>=0.2 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=0.2) {
	  if(signalRegion>=2) continue;
	  signalRegion = 2;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	}
	else {
	  if(signalRegion!=-1) continue;
	  signalRegion = 0;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	}
      }
      if(signalRegion==4) break; // break if event is signal region 3
    }

    
    if(foundGoodLep==false) continue;

    SelectedEvents++;
    if(isSignal && signalRegion==4) SelEvntSR3++;
    if(isSignal && signalRegion==3) SelEvntSR2++;
    if(isSignal && signalRegion==2) SelEvntSR1++;

    int NLep = numGoodLep;
    int NEl = numGoodEl;
    int NMu = numGoodMu;      

    // For the selected (e,mu) pair of leptons
    std::vector<TLorentzVector*> lepvec;
    TLorentzVector lep1;
    lep1.SetPtEtaPhiE(PT->at(firstPos), Eta->at(firstPos), Phi->at(firstPos), E->at(firstPos));
    TLorentzVector lep2;
    lep2.SetPtEtaPhiE(PT->at(secondPos), Eta->at(secondPos), Phi->at(secondPos), E->at(secondPos));
    TLorentzVector lepSum, jetSum, objSum, metVec, diffObjMet;
    metVec.SetPtEtaPhiE(MET, MET_Eta, MET_Phi, MET);
    lepvec.push_back(&lep1);
    lepvec.push_back(&lep2);

    // Fill Histogram for good jet when a good event is found based on Leptonic selection
    double htlep=0.0, htjet=0.0, htlepjet=0.0;
    int numGoodJet = 0;
    int jetFirstPos = -1;
    int jetSecondPos = -1;

    // Loop to fill variables dependent on all leptons in the event
    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;
      if(TMath::Abs(D0->at(objCtr))>100) continue;
      if(TMath::Abs(PID->at(objCtr))==11 && Iso->at(objCtr)>0.12) continue;
      if(TMath::Abs(PID->at(objCtr))==13 && Iso->at(objCtr)>0.15) continue;
      
      TLorentzVector lepSingle;
      lepSingle.SetPtEtaPhiE(PT->at(objCtr), Eta->at(objCtr), Phi->at(objCtr), E->at(objCtr));
      lepSum += lepSingle;
      htlep += PT->at(objCtr);
      
    } // End lepton loop

    // Loop to fill variables dependent on all leptons in the event
    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;
      
      numGoodJet++;

      TLorentzVector jetSingle;
      jetSingle.SetPtEtaPhiM(JetPT->at(jetCtr), JetEta->at(jetCtr), JetPhi->at(jetCtr), JetM->at(jetCtr));
      jetSum += jetSingle;
      htjet += JetPT->at(jetCtr);
      
      if(jetFirstPos==-1) {
	jetFirstPos = jetCtr;
	continue;
      }
      if(jetFirstPos!=-1 && jetSecondPos==-1) {
	jetSecondPos = jetCtr;
	continue;
      }
    } // End jet loop

    objSum = lepSum+jetSum;
    diffObjMet = objSum-metVec;
    htlepjet = htlep+htjet;
    
    int NJet = numGoodJet;
                
  } // End of Event Loop

  std::cout<<SelectedEvents<<"\t"
	   <<SelEvntSR1<<"\t"
	   <<SelEvntSR2<<"\t"
	   <<SelEvntSR3<<std::endl;

  *tot = SelectedEvents;
  *SR1 = SelEvntSR1;
  *SR2 = SelEvntSR2;
  *SR3 = SelEvntSR3;
}


void signalSigBeforeCut() {

  std::vector<TString> histLabel;
  std::vector<TString> outFile;
  std::vector<TChain*> t;
  std::vector<TString> dataPath;
  std::vector<double> crossSecVec; // in fb
  std::vector<double> sigmaCrossSecVec; // in fb
  std::vector<long> nSimuVec;

  //histLabel.push_back("HF Background");
  histLabel.push_back("(220, 20, DM)");
  histLabel.push_back("(324, 20, DM)");
  histLabel.push_back("(220, 20, 0.2)");
  histLabel.push_back("(220, 20, 2)");
  histLabel.push_back("(220, 20, 20)");
  histLabel.push_back("(220, 20, 200)");
  histLabel.push_back("(220, 40, 20)");
  
  //dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_DM_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_304_324_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_304_324_DM_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2mm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2mm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_20cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_20cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2m/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2m_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_180_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_180_220_2cm_Batch*.root");

  //crossSecVec.push_back(0);
  crossSecVec.push_back(903*0.014);
  crossSecVec.push_back(128*0.025);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  
  //sigmaCrossSecVec.push_back(0);
  sigmaCrossSecVec.push_back(54*0.014);
  sigmaCrossSecVec.push_back(10*0.025);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);

  //nSimuVec.push_back(1);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  nSimuVec.push_back(20*100000);
  
  for(int ctr=0; ctr<histLabel.size(); ctr++) {
    t.push_back(new TChain("SelectedObjects"));
    t[ctr]->Add(dataPath[ctr]);
    cout<<"Initialized "<<histLabel[ctr]<<" chain"<<endl;
    std::cout<<"No.of  "<<histLabel[ctr]<<" Entries: "<<t[ctr]->GetEntries()<<std::endl;
  }

  Calculator *C = new Calculator();
  C->histSuff = histLabel;
  cout<<"Initialized instance of plotter class"<<endl;

  for(int ctr=0; ctr<histLabel.size(); ctr++) {
    double lumi = 2.6; // In fb^{-1}
    double nEvent = lumi*crossSecVec[ctr];
    double wt = nEvent/nSimuVec[ctr];
    double err_xSec = TMath::Power(sigmaCrossSecVec[ctr]/crossSecVec[ctr],2);
    double err_Stat = TMath::Power(TMath::Sqrt(nSimuVec[ctr])/nSimuVec[ctr],2);
    double err_tot = TMath::Sqrt(err_xSec+err_Stat);

    double normBGSR1 = 4122.8;
    double normBGSR2 = 644.2;
    double normBGSR3 = 24.479;

    int selNumEvnt;
    int SR1;
    int SR2;
    int SR3;
    C->yieldCalc(t[ctr], &selNumEvnt, &SR1, &SR2, &SR3);
    
    std::cout<<histLabel[ctr]<<" & "<<crossSecVec[ctr]<<" & \\begin{tabular}[c]{@{}c@{}} "
	     <<setprecision(2)<<SR1*wt<<" $\\pm$ \\\\ "<<setprecision(2)<<SR1*wt*err_tot<<" \\end{tabular} & "<<setprecision(2)<<SR1*wt/TMath::Sqrt(SR1*wt+normBGSR1)<<" & \\begin{tabular}[c]{@{}c@{}} "
	     <<setprecision(2)<<SR2*wt<<" $\\pm$ \\\\ "<<setprecision(2)<<SR2*wt*err_tot<<" \\end{tabular} & "<<setprecision(2)<<SR2*wt/TMath::Sqrt(SR2*wt+normBGSR2)<<" & \\begin{tabular}[c]{@{}c@{}} "
	     <<setprecision(2)<<SR3*wt<<" $\\pm$ \\\\ "<<setprecision(2)<<SR3*wt*err_tot<<" \\end{tabular} & "<<setprecision(2)<<SR3*wt/TMath::Sqrt(SR3*wt+normBGSR3)
	     <<std::endl;
    
  }

}

