#include <iostream>
#include <cmath>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TLorentzVector.h"

int createVarOutTree(TChain* tree, TString outFileName){
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
  double PT_Lep, PT_El, PT_Mu, PT_Jet, leadPT_Lep, leadPT_El, leadPT_Mu, leadPT_Jet, subLeadPT_Lep, subLeadPT_El, subLeadPT_Mu, subLeadPT_Jet, diffPT_Lep, diffPT_Jet;

  auto varOutFile = new TFile(outFileName, "recreate");
  auto varTree = new TTree("varTree", "Input Variables List for Algorithms");
  varTree->Branch("NLep", &NLep);
  varTree->Branch("PT_Lep", &PT_Lep);
  varTree->Branch("leadPT_Lep", &leadPT_Lep);
  varTree->Branch("subLeadPT_Lep", &subLeadPT_Lep);
  varTree->Branch("NEl", &NEl);
  varTree->Branch("PT_El", &PT_El);
  varTree->Branch("leadPT_El", &leadPT_El);
  varTree->Branch("subLeadPT_El", &subLeadPT_El);
  varTree->Branch("NMu", &NMu);
  varTree->Branch("PT_Mu", &PT_Mu);
  varTree->Branch("leadPT_Mu", &leadPT_Mu);
  varTree->Branch("subLeadPT_Mu", &subLeadPT_Mu);
  varTree->Branch("NJet", &NJet);
  varTree->Branch("PT_Jet", &PT_Jet);
  varTree->Branch("leadPT_Jet", &leadPT_Jet);
  varTree->Branch("subLeadPT_Jet", &subLeadPT_Jet);
  varTree->Branch("diffPT_Lep", &diffPT_Lep);
  varTree->Branch("diffPT_Jet", &diffPT_Jet);
  
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

    // For the selected (e,mu) pair of leptons
    std::vector<TLorentzVector*> lepvec;
    TLorentzVector lep1;
    lep1.SetPtEtaPhiE(PT->at(firstPos), Eta->at(firstPos), Phi->at(firstPos), E->at(firstPos));
    TLorentzVector lep2;
    lep2.SetPtEtaPhiE(PT->at(secondPos), Eta->at(secondPos), Phi->at(secondPos), E->at(secondPos));
    TLorentzVector lepSum, jetSum, objSum, metVec;
    metVec.SetPtEtaPhiM(MET, 0.0, MET_Phi, 0.0);
    lepvec.push_back(&lep1);
    lepvec.push_back(&lep2);

    double st=0.0, htlep=0.0;
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
      st += TMath::Abs(PT->at(objCtr));
      htlep += TMath::Abs(PT->at(objCtr));
    } // End lepton loop

    // Loop for jet
    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;
      
      numGoodJet++;

      TLorentzVector jetSingle;
      jetSingle.SetPtEtaPhiM(JetPT->at(jetCtr), JetEta->at(jetCtr), JetPhi->at(jetCtr), JetM->at(jetCtr));
      jetSum += jetSingle;
      st += TMath::Abs(JetPT->at(jetCtr));

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

    // Fill the tree with variables from the selected event
    varTree->Fill();
    
  } // End of Event Loop

  varOutFile->Write();
    
  return SelectedEvents;
}

void runononesample(TString chainpath, TString outputname){

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add(chainpath);
  
  
  cout<<createVarOutTree(sigChain, outputname)
      <<" events selected from "
      <<sigChain->GetEntries()
      <<" SIGNAL events'"<<endl;
    
  return;
}

void execute(TString sigchainpath="./Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root",
	     TString backgroundchainpath="./Data/ppTobb_Cuts2/Objects_sorted_ppTobb_*.root") {

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add(sigchainpath);

  TChain *bkgChain = new TChain("SelectedObjects");
  bkgChain->Add(backgroundchainpath);

  cout<<createVarOutTree(sigChain, "sigVar.root")
      <<" events selected from "
      <<sigChain->GetEntries()
      <<" SIGNAL events'"<<endl;
  
  cout<<createVarOutTree(bkgChain, "bkgVar.root")
      <<" events selected from "
      <<bkgChain->GetEntries()
      <<" BACKGROUND events'"<<endl;

    return;
}
