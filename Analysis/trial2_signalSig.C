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

class Calculator {

public:
  Calculator();
  ~Calculator();
  int yieldCalc(TChain* signal,
		TString outFileName);

  // Declare histogram suffix
  std::vector<TString> histSuff;
};

Calculator::Calculator() {
}

Calculator::~Calculator() {
}

int Calculator::yieldCalc(TChain* signal, TString outFileName) {

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

  // Define the output file and corresponding branches
  double HtJet, dRLL, dPhiLepMETSelObj, YDelpObj, YUserObj, alphaT, Sphericity, Spherocity, MtLeadLepMET;
  double HtJet_SR1, dRLL_SR1, dPhiLepMETSelObj_SR1, YDelpObj_SR1, YUserObj_SR1, alphaT_SR1, Sphericity_SR1, Spherocity_SR1, MtLeadLepMET_SR1;
  double HtJet_SR2, dRLL_SR2, dPhiLepMETSelObj_SR2, YDelpObj_SR2, YUserObj_SR2, alphaT_SR2, Sphericity_SR2, Spherocity_SR2, MtLeadLepMET_SR2;
  double HtJet_SR3, dRLL_SR3, dPhiLepMETSelObj_SR3, YDelpObj_SR3, YUserObj_SR3, alphaT_SR3, Sphericity_SR3, Spherocity_SR3, MtLeadLepMET_SR3;
  
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
  auto varTree_SR1 = new TTree("varTree_SR1", "Input Variables List for Algorithms SR1");
  varTree_SR1->Branch("HtJet_SR1", &HtJet_SR1);
  varTree_SR1->Branch("dRLL_SR1", &dRLL_SR1);
  varTree_SR1->Branch("dPhiLepMETSelObj_SR1", &dPhiLepMETSelObj_SR1);
  varTree_SR1->Branch("YDelpObj_SR1", &YDelpObj_SR1);
  varTree_SR1->Branch("YUserObj_SR1", &YUserObj_SR1);
  varTree_SR1->Branch("alphaT_SR1", &alphaT_SR1);
  varTree_SR1->Branch("Sphericity_SR1", &Sphericity_SR1);
  varTree_SR1->Branch("Spherocity_SR1", &Spherocity_SR1);
  varTree_SR1->Branch("MtLeadLepMET_SR1", &MtLeadLepMET_SR1);
  auto varTree_SR2 = new TTree("varTree_SR2", "Input Variables List for Algorithms SR2");
  varTree_SR2->Branch("HtJet_SR2", &HtJet_SR2);
  varTree_SR2->Branch("dRLL_SR2", &dRLL_SR2);
  varTree_SR2->Branch("dPhiLepMETSelObj_SR2", &dPhiLepMETSelObj_SR2);
  varTree_SR2->Branch("YDelpObj_SR2", &YDelpObj_SR2);
  varTree_SR2->Branch("YUserObj_SR2", &YUserObj_SR2);
  varTree_SR2->Branch("alphaT_SR2", &alphaT_SR2);
  varTree_SR2->Branch("Sphericity_SR2", &Sphericity_SR2);
  varTree_SR2->Branch("Spherocity_SR2", &Spherocity_SR2);
  varTree_SR2->Branch("MtLeadLepMET_SR2", &MtLeadLepMET_SR2);
  auto varTree_SR3 = new TTree("varTree_SR3", "Input Variables List for Algorithms SR3");
  varTree_SR3->Branch("HtJet_SR3", &HtJet_SR3);
  varTree_SR3->Branch("dRLL_SR3", &dRLL_SR3);
  varTree_SR3->Branch("dPhiLepMETSelObj_SR3", &dPhiLepMETSelObj_SR3);
  varTree_SR3->Branch("YDelpObj_SR3", &YDelpObj_SR3);
  varTree_SR3->Branch("YUserObj_SR3", &YUserObj_SR3);
  varTree_SR3->Branch("alphaT_SR3", &alphaT_SR3);
  varTree_SR3->Branch("Sphericity_SR3", &Sphericity_SR3);
  varTree_SR3->Branch("Spherocity_SR3", &Spherocity_SR3);
  varTree_SR3->Branch("MtLeadLepMET_SR3", &MtLeadLepMET_SR3);

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

      lepPos.push_back(objCtr);
    }

    if(lepPos.size()<2) continue; // Not enough leptons qualified for event selection

    // Loop to select the leptons and the signal region
    for(int lep1Ctr=0; lep1Ctr<lepPos.size()-1; lep1Ctr++) {
      for(int lep2Ctr=lep1Ctr+1; lep2Ctr<lepPos.size(); lep2Ctr++) {
	if(PID->at(lepPos[lep1Ctr])*PID->at(lepPos[lep2Ctr]) != -11*13) continue;

	if(TMath::Abs(D0->at(lepPos[lep1Ctr]))>=1 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=1) { // break if event is signal region 3
	  signalRegion = 4;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	  break;
	}
	else if(TMath::Abs(D0->at(lepPos[lep1Ctr]))>=0.5 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=0.5) {
	  if(signalRegion>=3) continue;
	  signalRegion = 3;
	  firstPos = lepPos[lep1Ctr];
	  secondPos = lepPos[lep2Ctr];
	  foundGoodLep = true;
	}
	else if(TMath::Abs(D0->at(lepPos[lep1Ctr]))>=0.2 && TMath::Abs(D0->at(lepPos[lep2Ctr]))>=0.2) {
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
    if(signalRegion==4) SelEvntSR3++;
    if(signalRegion==3) SelEvntSR2++;
    if(signalRegion==2) SelEvntSR1++;

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

    dRLL = TMath::Abs(lep1.DeltaR(lep2));

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
      
      TLorentzVector lepSingle;
      lepSingle.SetPtEtaPhiE(PT->at(objCtr), Eta->at(objCtr), Phi->at(objCtr), E->at(objCtr));
      lepSum += lepSingle;
      htlep += PT->at(objCtr);
      
    } // End lepton loop

    MtLeadLepMET = TMath::Sqrt(2*MET*PT->at(firstPos)*(1-TMath::Cos(lep1.DeltaPhi(metVec))));

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

    // Fill the corresponding signal tree
    if(signalRegion==4) {
      HtJet_SR3 = HtJet;
      dRLL_SR3 = dRLL;
      dPhiLepMETSelObj_SR3 = dPhiLepMETSelObj;
      YDelpObj_SR3 = YDelpObj;
      YUserObj_SR3 = YUserObj;
      alphaT_SR3 = alphaT;
      Sphericity_SR3 = Sphericity;
      Spherocity_SR3 = Spherocity;
      MtLeadLepMET_SR3 = MtLeadLepMET;      
      varTree_SR3->Fill();
    }
    if(signalRegion==3) {
      HtJet_SR2 = HtJet;
      dRLL_SR2 = dRLL;
      dPhiLepMETSelObj_SR2 = dPhiLepMETSelObj;
      YDelpObj_SR2 = YDelpObj;
      YUserObj_SR2 = YUserObj;
      alphaT_SR2 = alphaT;
      Sphericity_SR2 = Sphericity;
      Spherocity_SR2 = Spherocity;
      MtLeadLepMET_SR2 = MtLeadLepMET;      
      varTree_SR2->Fill();
    }
    if(signalRegion==2) {
      HtJet_SR1 = HtJet;
      dRLL_SR1 = dRLL;
      dPhiLepMETSelObj_SR1 = dPhiLepMETSelObj;
      YDelpObj_SR1 = YDelpObj;
      YUserObj_SR1 = YUserObj;
      alphaT_SR1 = alphaT;
      Sphericity_SR1 = Sphericity;
      Spherocity_SR1 = Spherocity;
      MtLeadLepMET_SR1 = MtLeadLepMET;      
      varTree_SR1->Fill();
    }
    
  } // End of Event Loop

  varOutFile->Write();
  std::cout<<SelEvntSR1<<"\t"<<SelEvntSR2<<"\t"<<SelEvntSR3<<std::endl;
  return SelectedEvents;
}


void trial2_signalSig() {

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

  //outFile.push_back(TODO);
  outFile.push_back("BP_200_20_DM.root");
  outFile.push_back("BP_324_20_DM.root");
  outFile.push_back("BP_200_20_02.root");
  outFile.push_back("BP_200_20_2.root");
  outFile.push_back("BP_200_20_20.root");
  outFile.push_back("BP_200_20_200.root");
  outFile.push_back("BP_200_40_20.root");
  
  //dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_DM_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_304_324_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_304_324_DM_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2mm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2mm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_20cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_20cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2m/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2m_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_180_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_180_220_2cm_Batch*.root");

  //crossSecVec.push_back(TODO);
  crossSecVec.push_back(903*0.014);
  crossSecVec.push_back(128*0.025);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  crossSecVec.push_back(903);
  
  //sigmaCrossSecVec.push_back(TODO);
  sigmaCrossSecVec.push_back(54*0.014);
  sigmaCrossSecVec.push_back(10*0.025);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);
  sigmaCrossSecVec.push_back(54);

  //nSimuVec.push_back(TODO);
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
    // Compute the error in wt
    int selNumEvnt = C->yieldCalc(t[ctr],outFile[ctr]);
    /*
      cout<<"For signal BP: "<<histLabel[ctr]
      <<" - Selected:"<<selNumEvnt<<" events from "<<t[ctr]->GetEntries()
      <<". Yield = "<<wt*selNumEvnt<<endl;
    */
  }

}

