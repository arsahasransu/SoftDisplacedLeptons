#include <iostream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMath.h"

// Number code for d0 based control and signal regions.
// 0 - CR1
// 1 - CR2
// 2 - SR1
// 3 - SR2
// 4 - SR3

// void function essential for bug-free root compilation
void VariableMaker(void){
  return;
}

// Functions required for computing sphericity and spherocity
TMatrixD makeMomentumTensor3D(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total(3,3);
  double normaliser = 0.0;
  TVector3 vect;

  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      normaliser=0.0;
      total[ii][jj]=0.0;
      for(unsigned int vec=0; vec<lvarray.size(); vec++) {
	vect = lvarray[vec]->Vect();
	total[ii][jj] += vect[ii]*vect[jj];
	normaliser += vect.Mag()*vect.Mag();
      }
      if(normaliser>0) {
	total[ii][jj] /= normaliser;
      }
    }
  }
  return total;
}

// Functions required for computing sphericity and spherocity
double sphericity(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total = makeMomentumTensor3D(lvarray);

  TVectorD eigenvals;
  auto eigenvectors = total.EigenVectors(eigenvals);

  double sphericity = 1.5*(eigenvals[1]+eigenvals[2]);
  return sphericity;
}

double transversespherocity(std::vector<TLorentzVector*> lvarray) {

  Double_t lowestval = 100000000.0;

  TVector3 workvec;
  double workerNumerator = 0.0;
  double workerDenominator = 0.0;
  TVector3 unitvec;

  for(unsigned int vecCtr=0; vecCtr<lvarray.size(); vecCtr++) {
    TLorentzVector* candvec = lvarray[vecCtr];

    unitvec.SetXYZ(candvec->Px()/candvec->Pt(), candvec->Py()/candvec->Pt(), 0.0);
    workvec.SetXYZ(0, 0, 0);
    workerNumerator = 0.0;
    workerDenominator = 0.0;

    for(unsigned int vecCtr2=0; vecCtr2<lvarray.size(); vecCtr2++) {
      TLorentzVector* workLVvec = lvarray[vecCtr2];

      workerDenominator += workLVvec->Pt();
      workvec.SetXYZ(workLVvec->Px(), workLVvec->Py(), 0.0);
      workerNumerator += workvec.Cross(unitvec).Mag();
    }
    lowestval = TMath::Min(lowestval, workerNumerator/workerDenominator);
  }

  double spherocity = TMath::Power(lowestval*TMath::Pi()*0.5, 2);
  return spherocity;
}

// Function to write read from objects and create variables in a Tree for one file
int createVarOutTree(TChain* tree, TString outFileName, bool signal, int d0_choice){
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

  // Link the defined variables to the leaf on trees
  tree->SetBranchAddress("NumElec",&numElec);
  tree->SetBranchAddress("NumMuon",&numMuon);
  tree->SetBranchAddress("MET", &MET);
  tree->SetBranchAddress("MET_Eta", &MET_Eta);
  tree->SetBranchAddress("MET_Phi", &MET_Phi);
  tree->SetBranchAddress("PID", &PID);
  tree->SetBranchAddress("simulatedE", &E);
  tree->SetBranchAddress("PT", &PT);
  tree->SetBranchAddress("Eta", &Eta);
  tree->SetBranchAddress("Phi", &Phi);
  tree->SetBranchAddress("Iso", &Iso);
  tree->SetBranchAddress("D0", &D0);
  tree->SetBranchAddress("Dz", &Dz);
  tree->SetBranchAddress("HT", &HT);
  tree->SetBranchAddress("NumJet", &numJet);
  tree->SetBranchAddress("JetEta", &JetEta);
  tree->SetBranchAddress("JetPhi", &JetPhi);
  tree->SetBranchAddress("JetPT", &JetPT);
  tree->SetBranchAddress("JetM", &JetM);

  // Define the output file and corresponding branches
  // From plotVarDisb_Objects_2
  double NLep, NEl, NMu, NJet;
  double HtJet, dRLL, dPhiLepMETSelObj, YDelpObj, YUserObj, alphaT, Sphericity, Spherocity, MtLeadLepMET;
  
  auto varOutFile = new TFile(outFileName, "recreate");
  auto varTree = new TTree("varTree", "Input Variables List for Algorithms");
  varTree->Branch("HtJet", &HtJet);
  varTree->Branch("dRLL", &dRLL);
  varTree->Branch("dPhiLepMETSelObj", &dPhiLepMETSelObj);
  varTree->Branch("YDelpObj", &YDelpObj);
  varTree->Branch("YUserObj", &YUserObj);
  varTree->Branch("alphaT", &alphaT);
  varTree->Branch("Sphericity", &Sphericity);
  varTree->Branch("Spherocity", &Spherocity);
  varTree->Branch("MtLeadLepMET", &MtLeadLepMET);
  
  // Statistic variable
  int SelectedEvents = 0;

  // Open Event Loop
  for(int evtCtr=0; evtCtr<tree->GetEntries() /*&& evtCtr<10*/; evtCtr++) {

    if(evtCtr%100000==0) cout<<tree->GetEntries()<<" total. Ongoing event: "<<evtCtr<<endl; 
    
    tree->GetEntry(evtCtr);
    double lenObj = PID->size();

    // Select for event based on oppositely charged pair of (e,mu)
    bool foundGoodLep = false;
    int firstPos=-1;
    int secondPos=-1;
    int numGoodLep = 0;
    int numGoodEl = 0;
    int numGoodMu = 0;

    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;

      // Signal control for d0 Region. Uncomment one of them only.
      if(signal && d0_choice==0) if(TMath::Abs(D0->at(objCtr))>0.1) continue; // CR1
      if(signal && d0_choice==1) if(TMath::Abs(D0->at(objCtr))<0.1 || TMath::Abs(D0->at(objCtr))>0.2) continue; // CR2
      if(signal && d0_choice==2) if(TMath::Abs(D0->at(objCtr))<0.2) continue; // SR1
      if(signal && d0_choice==3) if(TMath::Abs(D0->at(objCtr))<0.5) continue; // SR2
      if(signal && d0_choice==4) if(TMath::Abs(D0->at(objCtr))>1) continue; // SR3

      numGoodLep++;
      if(TMath::Abs(PID->at(objCtr))==11) numGoodEl++;
      if(TMath::Abs(PID->at(objCtr))==13) numGoodMu++;
      
      if(firstPos==-1) {
	firstPos = objCtr;
	continue;
      }
      if(secondPos==-1) {
	if(PID->at(firstPos)*PID->at(objCtr) != -11*13) continue;
	secondPos = objCtr;
	foundGoodLep = true;
	break;
      }
      cout<<"neverprint"<<endl;
    }
    
    if(foundGoodLep==false) continue;
    
    
    SelectedEvents++;
    
    NLep = numGoodLep;
    NEl = numGoodEl;
    NMu = numGoodMu;      

    // For the selected (e,mu) pair of leptons
    std::vector<TLorentzVector*> lepvec;
    TLorentzVector lep1;
    lep1.SetPtEtaPhiE(PT->at(firstPos), Eta->at(firstPos), Phi->at(firstPos), E->at(firstPos));
    TLorentzVector lep2;
    lep2.SetPtEtaPhiE(PT->at(secondPos), Eta->at(secondPos), Phi->at(secondPos), E->at(secondPos));
    TLorentzVector lepSum, jetSum, objSum, metVec, diffObjMet;
    metVec.SetPtEtaPhiM(MET, 0.0, MET_Phi, 0.0);
    lepvec.push_back(&lep1);
    lepvec.push_back(&lep2);

    dRLL = TMath::Abs(lep1.DeltaR(lep2));
		   
    double htlep=0.0, htjet=0.0, htlepjet=0.0;
    int numGoodJet=0;
    int jetFirstPos = -1;
    int jetSecondPos = -1;

    // Loop to fill variables dependent on all leptons in the event
    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;

      TLorentzVector lepSingle;
      lepSingle.SetPtEtaPhiE(PT->at(objCtr), Eta->at(objCtr), Phi->at(objCtr), E->at(objCtr));
      lepSum += lepSingle;
      htlep += PT->at(objCtr);
      
    } // End lepton loop
    
    MtLeadLepMET = TMath::Sqrt(2*MET*PT->at(firstPos)*(1-TMath::Cos(lep1.DeltaPhi(metVec))));

    // Loop for jet
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
    
    NJet = numGoodJet;
    
    dPhiLepMETSelObj = TMath::Abs(lep1.DeltaPhi(objSum));
    if(numGoodJet!=0) {
      HtJet = htjet;
    }
    else {
      HtJet = 0;
    }
    
    YDelpObj = MET/TMath::Sqrt(HT);
        
    YUserObj = TMath::Abs(objSum.Pt())/TMath::Sqrt(htlep);
    double mt = TMath::Sqrt(htlep*htlep-lepSum.Pt()*lepSum.Pt());
    alphaT = PT->at(secondPos)/mt;
    Sphericity = sphericity(lepvec);
    Spherocity = transversespherocity(lepvec);
    
    // Fill the tree with variables from the selected event
    varTree->Fill();
    
  } // End of Event Loop
  
  varOutFile->Write();
  
  return SelectedEvents;
}

// Function to run on one file
void runononesample(TString chainpath, TString outputname, bool isSignal, int d0_choice=-1){

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add(chainpath);
  
    string typeofdata= "background";
    if (isSignal)
        typeofdata="signal";
    
    cout<<createVarOutTree(sigChain, outputname, isSignal, d0_choice)
      <<" events selected from "
      <<sigChain->GetEntries()
    
      <<" "<<typeofdata<< " events'"<<endl;
    
  return;
}

// Function to run on multiple files
void execute(int d0_choice=-1,
	     TString sigOutFile="signal.root",
	     TString bkgOutFile="background.root",
	     TString sigchainpath="./Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root",
	     TString backgroundchainpath="./Data/ppTobb_Cuts2/Objects_sorted_ppTobb_*.root") {

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add(sigchainpath);

  TChain *bkgChain = new TChain("SelectedObjects");
  bkgChain->Add(backgroundchainpath);

  cout<<createVarOutTree(sigChain, sigOutFile, true, d0_choice)
      <<" events selected from "
      <<sigChain->GetEntries()
      <<" SIGNAL events'"<<endl;
  
  cout<<createVarOutTree(bkgChain, bkgOutFile, false, d0_choice)
      <<" events selected from "
      <<bkgChain->GetEntries()
      <<" BACKGROUND events'"<<endl;

    return;
}
