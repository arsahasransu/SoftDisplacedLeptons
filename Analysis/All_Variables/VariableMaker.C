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
  std::vector<double> PT_Lep, PT_El, PT_Mu, PT_Jet, leadPT_Lep, leadPT_El, leadPT_Mu, leadPT_Jet, subLeadPT_Lep, subLeadPT_El, subLeadPT_Mu, subLeadPT_Jet, diffPT_Lep, diffPT_Jet;
  // From plotVarDisb_Objects_3
  std::vector<double> Eta_Lep, Eta_El, Eta_Mu, Eta_Jet, leadEta_Lep, leadEta_El, leadEta_Mu, leadEta_Jet, subLeadEta_Lep, subLeadEta_El, subLeadEta_Mu, subLeadEta_Jet, diffEta_Lep, diffEta_Jet;
  // From plotVarDisb_Objects_4
  std::vector<double> Phi_Lep, Phi_El, Phi_Mu, Phi_Jet, leadPhi_Lep, leadPhi_El, leadPhi_Mu, leadPhi_Jet, subLeadPhi_Lep, subLeadPhi_El, subLeadPhi_Mu, subLeadPhi_Jet, diffPhi_Lep, diffPhi_Jet;
  // From plotVarDisb_Objects_5
  std::vector<double> D0_Lep, D0_El, D0_Mu, leadD0_Lep, leadD0_El, leadD0_Mu, subLeadD0_Lep, subLeadD0_El, subLeadD0_Mu, diffD0_Lep;
  // From plotVarDisb_Objects_6
  std::vector<double> Iso_Lep, Iso_El, Iso_Mu, leadIso_Lep, leadIso_El, leadIso_Mu, subLeadIso_Lep, subLeadIso_El, subLeadIso_Mu, diffIso_Lep;
  // From plotVarDisb_Objects_7
  std::vector<double> Met, MetLepJet, MetLep, MetJet, MetDiffMetLepJet, Met_Phi, MetLepJet_Phi, MetLep_Phi, MetJet_Phi, MetDiffMetLepJet_Phi, Ht, HtLepJet, HtLep, HtJet, HtDiffHtLepJet, DiffMetHt, DiffMetHtLepJet, DiffMetHtLep, DiffMetHtJet;
  // From plotVarDisb_Objects_8
  std::vector<double> dRLL, dEtaLepJet, dPhiLepJet, dRLepJet, dPhiLepMET, dPhiLepMETSelObj;
  // From plotVarDisb_Objects_Misc
  std::vector<double> YDelpObj, YUserObj, ST, alphaT, MLL, MLLMET, M, Sphericity, Spherocity, MtLepMET, MtLepMET_El, MtLepMET_Mu, MtLeadLepMET;
      
  auto varOutFile = new TFile(outFileName, "recreate");
  auto varTree = new TTree("varTree", "Input Variables List for Algorithms");
  //varTree->Branch("NLep", &NLep);
  //varTree->Branch("PT_Lep", &PT_Lep);
  //varTree->Branch("leadPT_Lep", &leadPT_Lep);
  //varTree->Branch("subLeadPT_Lep", &subLeadPT_Lep);
  //varTree->Branch("NEl", &NEl);
  //varTree->Branch("PT_El", &PT_El);
  //varTree->Branch("leadPT_El", &leadPT_El);
  //varTree->Branch("subLeadPT_El", &subLeadPT_El);
  //varTree->Branch("NMu", &NMu);
  //varTree->Branch("PT_Mu", &PT_Mu);
  //varTree->Branch("leadPT_Mu", &leadPT_Mu);
  //varTree->Branch("subLeadPT_Mu", &subLeadPT_Mu);
  //varTree->Branch("NJet", &NJet);
  //varTree->Branch("PT_Jet", &PT_Jet);
  //varTree->Branch("leadPT_Jet", &leadPT_Jet);
  //varTree->Branch("subLeadPT_Jet", &subLeadPT_Jet);
  //varTree->Branch("diffPT_Lep", &diffPT_Lep);
  //varTree->Branch("diffPT_Jet", &diffPT_Jet);
  ////////////////////////////////////////////////
  //varTree->Branch("Eta_Lep", &Eta_Lep);
  //varTree->Branch("leadEta_Lep", &leadEta_Lep);
  //varTree->Branch("subLeadEta_Lep", &subLeadEta_Lep);
  varTree->Branch("Eta_El", &Eta_El);
  //varTree->Branch("leadEta_El", &leadEta_El);
  //varTree->Branch("subLeadEta_El", &subLeadEta_El);
  varTree->Branch("Eta_Mu", &Eta_Mu);
  //varTree->Branch("leadEta_Mu", &leadEta_Mu);
  //varTree->Branch("subLeadEta_Mu", &subLeadEta_Mu);
  //varTree->Branch("Eta_Jet", &Eta_Jet);
  //varTree->Branch("leadEta_Jet", &leadEta_Jet);
  //varTree->Branch("subLeadEta_Jet", &subLeadEta_Jet);
  //varTree->Branch("diffEta_Lep", &diffEta_Lep);
  //varTree->Branch("diffEta_Jet", &diffEta_Jet);
  ////////////////////////////////////////////////
  //varTree->Branch("Phi_Lep", &Phi_Lep);
  //varTree->Branch("leadPhi_Lep", &leadPhi_Lep);
  //varTree->Branch("subLeadPhi_Lep", &subLeadPhi_Lep);
  //varTree->Branch("Phi_El", &Phi_El);
  //varTree->Branch("leadPhi_El", &leadPhi_El);
  //varTree->Branch("subLeadPhi_El", &subLeadPhi_El);
  //varTree->Branch("Phi_Mu", &Phi_Mu);
  //varTree->Branch("leadPhi_Mu", &leadPhi_Mu);
  //varTree->Branch("subLeadPhi_Mu", &subLeadPhi_Mu);
  //varTree->Branch("Phi_Jet", &Phi_Jet);
  //varTree->Branch("leadPhi_Jet", &leadPhi_Jet);
  //varTree->Branch("subLeadPhi_Jet", &subLeadPhi_Jet);
  //varTree->Branch("diffPhi_Lep", &diffPhi_Lep);
  //varTree->Branch("diffPhi_Jet", &diffPhi_Jet);
  ////////////////////////////////////////////////
  //varTree->Branch("D0_Lep", &D0_Lep);
  //varTree->Branch("leadD0_Lep", &leadD0_Lep);
  //varTree->Branch("subLeadD0_Lep", &subLeadD0_Lep);
  //varTree->Branch("D0_El", &D0_El);
  //varTree->Branch("leadD0_El", &leadD0_El);
  //varTree->Branch("subLeadD0_El", &subLeadD0_El);
  //varTree->Branch("D0_Mu", &D0_Mu);
  //varTree->Branch("leadD0_Mu", &leadD0_Mu);
  //varTree->Branch("subLeadD0_Mu", &subLeadD0_Mu);
  //varTree->Branch("diffD0_Lep", &diffD0_Lep);
  ////////////////////////////////////////////////
  //varTree->Branch("Iso_Lep", &Iso_Lep);
  //varTree->Branch("leadIso_Lep", &leadIso_Lep);
  //varTree->Branch("subLeadIso_Lep", &subLeadIso_Lep);
  //varTree->Branch("Iso_El", &Iso_El);
  //varTree->Branch("leadIso_El", &leadIso_El);
  //varTree->Branch("subLeadIso_El", &subLeadIso_El);
  //varTree->Branch("Iso_Mu", &Iso_Mu);
  //varTree->Branch("leadIso_Mu", &leadIso_Mu);
  //varTree->Branch("subLeadIso_Mu", &subLeadIso_Mu);
  //varTree->Branch("diffIso_Lep", &diffIso_Lep);
  ////////////////////////////////////////////////  
  //varTree->Branch("Met", &Met);
  //varTree->Branch("MetLepJet", &MetLepJet);
  //varTree->Branch("MetLep", &MetLep);
  //varTree->Branch("MetJet", &MetJet);
  //varTree->Branch("MetDiffMetLepJet", &MetDiffMetLepJet);
  //varTree->Branch("Met_Phi", &Met_Phi);
  //varTree->Branch("MetLepJet_Phi", &MetLepJet_Phi);
  //varTree->Branch("MetLep_Phi", &MetLep_Phi);
  //varTree->Branch("MetJet_Phi", &MetJet_Phi);
  //varTree->Branch("MetDiffMetLepJet_Phi", &MetDiffMetLepJet_Phi);
  //varTree->Branch("Ht", &Ht);
  //varTree->Branch("HtLepJet", &HtLepJet);
  //varTree->Branch("HtLep", &HtLep);
  varTree->Branch("HtJet", &HtJet);
  //varTree->Branch("HtDiffHtLepJet", &HtDiffHtLepJet);
  //varTree->Branch("DiffMetHt", &DiffMetHt);
  //varTree->Branch("DiffMetHtLepJet", &DiffMetHtLepJet);
  //varTree->Branch("DiffMetHtLep", &DiffMetHtLep);
  //varTree->Branch("DiffMetHtJet", &DiffMetHtJet);
  ////////////////////////////////////////////////  
  varTree->Branch("dRLL", &dRLL);
  //varTree->Branch("dEtaLepJet", &dEtaLepJet);
  //varTree->Branch("dPhiLepJet", &dPhiLepJet);
  //varTree->Branch("dRLepJet", &dRLepJet);
  //varTree->Branch("dPhiLepMET", &dPhiLepMET);
  varTree->Branch("dPhiLepMETSelObj", &dPhiLepMETSelObj);
  ////////////////////////////////////////////////  
  varTree->Branch("YDelpObj", &YDelpObj);
  varTree->Branch("YUserObj", &YUserObj);
  //varTree->Branch("ST", &ST);
  varTree->Branch("alphaT", &alphaT);
  //varTree->Branch("MLL", &MLL);
  //varTree->Branch("MLLMET", &MLLMET);
  //varTree->Branch("M", &M);
  varTree->Branch("Sphericity", &Sphericity);
  varTree->Branch("Spherocity", &Spherocity);
  varTree->Branch("MtLepMET", &MtLepMET);
  //varTree->Branch("MtLepMET_El", &MtLepMET_El);
  //varTree->Branch("MtLepMET_Mu", &MtLepMET_Mu);
  varTree->Branch("MtLeadLepMET", &MtLeadLepMET);
  
  // Statistic variable
  int SelectedEvents = 0;

  // Open Event Loop
  for(int evtCtr=0; evtCtr<tree->GetEntries() /*&& evtCtr<10*/; evtCtr++) {

    PT_Lep.clear();
    PT_El.clear();
    PT_Mu.clear();
    PT_Jet.clear();
    leadPT_Lep.clear();
    leadPT_El.clear();
    leadPT_Mu.clear();
    leadPT_Jet.clear();
    subLeadPT_Lep.clear();
    subLeadPT_El.clear();
    subLeadPT_Mu.clear();
    subLeadPT_Jet.clear();
    diffPT_Lep.clear();
    diffPT_Jet.clear();
    ////////////////////////////
    Eta_Lep.clear();
    Eta_El.clear();
    Eta_Mu.clear();
    Eta_Jet.clear();
    leadEta_Lep.clear();
    leadEta_El.clear();
    leadEta_Mu.clear();
    leadEta_Jet.clear();
    subLeadEta_Lep.clear();
    subLeadEta_El.clear();
    subLeadEta_Mu.clear();
    subLeadEta_Jet.clear();
    diffEta_Lep.clear();
    diffEta_Jet.clear();
    diffD0_Lep.clear();
    ////////////////////////////
    Phi_Lep.clear();
    Phi_El.clear();
    Phi_Mu.clear();
    Phi_Jet.clear();
    leadPhi_Lep.clear();
    leadPhi_El.clear();
    leadPhi_Mu.clear();
    leadPhi_Jet.clear();
    subLeadPhi_Lep.clear();
    subLeadPhi_El.clear();
    subLeadPhi_Mu.clear();
    subLeadPhi_Jet.clear();
    diffPhi_Lep.clear();
    diffPhi_Jet.clear();
    ////////////////////////////
    D0_Lep.clear();
    D0_El.clear();
    D0_Mu.clear();
    leadD0_Lep.clear();
    leadD0_El.clear();
    leadD0_Mu.clear();
    subLeadD0_Lep.clear();
    subLeadD0_El.clear();
    subLeadD0_Mu.clear();
    diffPhi_Lep.clear();
    ////////////////////////////
    Iso_Lep.clear();
    Iso_El.clear();
    Iso_Mu.clear();
    leadIso_Lep.clear();
    leadIso_El.clear();
    leadIso_Mu.clear();
    subLeadIso_Lep.clear();
    subLeadIso_El.clear();
    subLeadIso_Mu.clear();
    diffPhi_Lep.clear();
    diffIso_Lep.clear();
    ////////////////////////////////////////////////  
    Met.clear();
    MetLepJet.clear();
    MetLep.clear();
    MetJet.clear();
    MetDiffMetLepJet.clear();
    Met_Phi.clear();
    MetLepJet_Phi.clear();
    MetLep_Phi.clear();
    MetJet_Phi.clear();
    MetDiffMetLepJet_Phi.clear();
    Ht.clear();
    HtLepJet.clear();
    HtLep.clear();
    HtJet.clear();
    HtDiffHtLepJet.clear();
    DiffMetHt.clear();
    DiffMetHtLepJet.clear();
    DiffMetHtLep.clear();
    DiffMetHtJet.clear();
    ////////////////////////////////////////////////  
    dRLL.clear();
    dEtaLepJet.clear();
    dPhiLepJet.clear();
    dRLepJet.clear();
    dPhiLepMET.clear();
    dPhiLepMETSelObj.clear();
    ////////////////////////////////////////////////  
    YDelpObj.clear();
    YUserObj.clear();
    ST.clear();
    alphaT.clear();
    MLL.clear();
    MLLMET.clear();
    M.clear();
    Sphericity.clear();
    Spherocity.clear();
    MtLepMET.clear();
    MtLepMET_El.clear();
    MtLepMET_Mu.clear();
    MtLeadLepMET.clear();

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
      if(signal && d0_choice==4) if(TMath::Abs(D0->at(objCtr))<1) continue; // SR3

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

    dRLL.push_back(TMath::Abs(lep1.DeltaR(lep2)));
		   
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
      
      PT_Lep.push_back(PT->at(objCtr));
      Eta_Lep.push_back(Eta->at(objCtr));
      Phi_Lep.push_back(Phi->at(objCtr));
      D0_Lep.push_back(D0->at(objCtr));
      Iso_Lep.push_back(Iso->at(objCtr));
      if(TMath::Abs(PID->at(objCtr))==11) {
	PT_Mu.push_back(PT->at(objCtr));
	Eta_Mu.push_back(Eta->at(objCtr));
	Phi_Mu.push_back(Phi->at(objCtr));
	D0_Mu.push_back(D0->at(objCtr));
	Iso_Mu.push_back(Iso->at(objCtr));
	MtLepMET_El.push_back(TMath::Sqrt(2*MET*PT->at(objCtr)*(1-TMath::Cos(lepSingle.DeltaPhi(metVec)))));
      }
      if(TMath::Abs(PID->at(objCtr))==13) {
	PT_El.push_back(PT->at(objCtr));
	Eta_El.push_back(Eta->at(objCtr));
	Phi_El.push_back(Phi->at(objCtr));
	D0_El.push_back(D0->at(objCtr));
	Iso_El.push_back(Iso->at(objCtr));
	MtLepMET_Mu.push_back(TMath::Sqrt(2*MET*PT->at(objCtr)*(1-TMath::Cos(lepSingle.DeltaPhi(metVec)))));
      }
      
    } // End lepton loop
    
    leadPT_Lep.push_back(PT->at(firstPos));
    subLeadPT_Lep.push_back(PT->at(secondPos));
    diffPT_Lep.push_back(TMath::Abs(PT->at(firstPos)-PT->at(secondPos)));
    leadEta_Lep.push_back(Eta->at(firstPos));
    subLeadEta_Lep.push_back(Eta->at(secondPos));
    diffEta_Lep.push_back(TMath::Abs(Eta->at(firstPos)-Eta->at(secondPos)));
    leadPhi_Lep.push_back(Phi->at(firstPos));
    subLeadPhi_Lep.push_back(Phi->at(secondPos));
    diffPhi_Lep.push_back(TMath::Abs(Phi->at(firstPos)-Phi->at(secondPos)));
    leadD0_Lep.push_back(D0->at(firstPos));
    subLeadD0_Lep.push_back(D0->at(secondPos));
    diffD0_Lep.push_back(TMath::Abs(D0->at(firstPos)-D0->at(secondPos)));
    leadIso_Lep.push_back(Iso->at(firstPos));
    subLeadIso_Lep.push_back(Iso->at(secondPos));
    diffIso_Lep.push_back(TMath::Abs(Iso->at(firstPos)-Iso->at(secondPos)));
    if(TMath::Abs(PID->at(firstPos))==11) {
      leadPT_El.push_back(PT->at(firstPos));
      subLeadPT_Mu.push_back(PT->at(secondPos));      
      leadEta_El.push_back(Eta->at(firstPos));
      subLeadEta_Mu.push_back(Eta->at(secondPos));      
      leadPhi_El.push_back(Phi->at(firstPos));
      subLeadPhi_Mu.push_back(Phi->at(secondPos));      
      leadD0_El.push_back(D0->at(firstPos));
      subLeadD0_Mu.push_back(D0->at(secondPos));      
      leadIso_El.push_back(Iso->at(firstPos));
      subLeadIso_Mu.push_back(Iso->at(secondPos));      
    }
    if(TMath::Abs(PID->at(firstPos))==13) {
      leadPT_Mu.push_back(PT->at(firstPos));
      subLeadPT_El.push_back(PT->at(secondPos));      
      leadEta_Mu.push_back(Eta->at(firstPos));
      subLeadEta_El.push_back(Eta->at(secondPos));      
      leadPhi_Mu.push_back(Phi->at(firstPos));
      subLeadPhi_El.push_back(Phi->at(secondPos));      
      leadD0_Mu.push_back(D0->at(firstPos));
      subLeadD0_El.push_back(D0->at(secondPos));      
      leadIso_Mu.push_back(Iso->at(firstPos));
      subLeadIso_El.push_back(Iso->at(secondPos));      
    }
    MtLeadLepMET.push_back(TMath::Sqrt(2*MET*PT->at(firstPos)*(1-TMath::Cos(lep1.DeltaPhi(metVec)))));
    MtLepMET.push_back(TMath::Sqrt(2*MET*lepSum.Pt()*(1-TMath::Cos(lepSum.DeltaPhi(metVec)))));

    // Loop for jet
    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;

      PT_Jet.push_back(JetPT->at(jetCtr));
      Eta_Jet.push_back(JetEta->at(jetCtr));
      Phi_Jet.push_back(JetPhi->at(jetCtr));
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
    if(numGoodJet>=1) {
      TLorentzVector leadJet;
      leadJet.SetPtEtaPhiM(JetPT->at(jetFirstPos), JetEta->at(jetFirstPos), JetPhi->at(jetFirstPos), JetM->at(jetFirstPos));
      leadPT_Jet.push_back(JetPT->at(jetFirstPos));
      leadEta_Jet.push_back(JetEta->at(jetFirstPos));
      leadPhi_Jet.push_back(JetPhi->at(jetFirstPos));
      dEtaLepJet.push_back(TMath::Abs(Eta->at(firstPos)-JetEta->at(jetFirstPos)));
      dPhiLepJet.push_back(TMath::Abs(lep1.DeltaPhi(leadJet)));
      dRLepJet.push_back(lep1.DeltaR(leadJet));
    }
    if(numGoodJet>=2) {
      subLeadPT_Jet.push_back(JetPT->at(jetSecondPos));
      diffPT_Jet.push_back(JetPT->at(jetFirstPos)-JetPT->at(jetSecondPos));
      subLeadEta_Jet.push_back(JetEta->at(jetSecondPos));
      diffEta_Jet.push_back(JetEta->at(jetFirstPos)-JetEta->at(jetSecondPos));
      subLeadPhi_Jet.push_back(JetPhi->at(jetSecondPos));
      diffPhi_Jet.push_back(JetPhi->at(jetFirstPos)-JetPhi->at(jetSecondPos));
    }
    
    Met.push_back(MET);
    MetDiffMetLepJet.push_back(TMath::Abs(diffObjMet.Pt()));
    Met_Phi.push_back(MET_Phi);
    MetDiffMetLepJet_Phi.push_back(diffObjMet.Phi());
    Ht.push_back(HT);
    HtDiffHtLepJet.push_back(TMath::Abs(HT-htlepjet));
    dPhiLepMET.push_back(TMath::Abs(lep1.DeltaPhi(metVec)));
    
    MetLepJet.push_back(objSum.Pt());
    MetLep.push_back(lepSum.Pt());
    MetLepJet_Phi.push_back(objSum.Phi());
    MetLep_Phi.push_back(lepSum.Phi());
    HtLepJet.push_back(htlepjet);
    HtLep.push_back(htlep);
    DiffMetHtLepJet.push_back(TMath::Abs(objSum.Pt()-htlepjet));
    DiffMetHtLep.push_back(TMath::Abs(lepSum.Pt()-htlep));
    dPhiLepMETSelObj.push_back(TMath::Abs(lep1.DeltaPhi(objSum)));
    
    if(numGoodJet!=0) {
      MetJet.push_back(jetSum.Pt());
      MetJet_Phi.push_back(jetSum.Phi());
      HtJet.push_back(htjet);
      DiffMetHtJet.push_back(TMath::Abs(jetSum.Pt()-htjet));
    }
    
    DiffMetHt.push_back(TMath::Abs(MET-HT));
    YDelpObj.push_back(MET/TMath::Sqrt(HT));
        
    YUserObj.push_back(TMath::Abs(objSum.Pt())/TMath::Sqrt(htlep));
    ST.push_back(htlepjet+MET);
    double mt = TMath::Sqrt(htlep*htlep-lepSum.Pt()*lepSum.Pt());
    alphaT.push_back(PT->at(secondPos)/mt);
    MLL.push_back(TMath::Abs((lep1+lep2).M()));
    MLLMET.push_back(TMath::Abs((lep1+lep2+metVec).M()));
    M.push_back(TMath::Abs((objSum+metVec).M()));
    Sphericity.push_back(sphericity(lepvec));
    Spherocity.push_back(transversespherocity(lepvec));
    
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
	     TString sigOutFileBP_100_200="signal_BP_100_200.root",
	     TString bkgOutFile="background.root",
	     TString sigchainpath="/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root",
	     TString sigchainpathBP_100_200="/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DislacedLepton_BP_100_200/Objects_sorted_DisplacedLepton_BP_100_200_*.root",
	     TString backgroundchainpath="/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root") {

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add(sigchainpath);

  TChain *sigChainBP_100_200 = new TChain("SelectedObjects");
  sigChainBP_100_200->Add(sigchainpathBP_100_200);

  TChain *bkgChain = new TChain("SelectedObjects");
  bkgChain->Add(backgroundchainpath);

  cout<<createVarOutTree(sigChain, sigOutFile, true, d0_choice)
      <<" events selected from "
      <<sigChain->GetEntries()
      <<" SIGNAL events'"<<endl;
  
  cout<<createVarOutTree(sigChainBP_100_200, sigOutFileBP_100_200, true, d0_choice)
      <<" events selected from "
      <<sigChainBP_100_200->GetEntries()
      <<" SIGNAL_BP_100_200 events'"<<endl;

  cout<<createVarOutTree(bkgChain, bkgOutFile, false, d0_choice)
      <<" events selected from "
      <<bkgChain->GetEntries()
      <<" BACKGROUND events'"<<endl;

    return;
}
