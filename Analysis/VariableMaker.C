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
  double dEtaLL, dPhiLL, dRLL, MLL;
  auto varOutFile = new TFile(outFileName, "recreate");
  auto varTree = new TTree("varTree", "Input Variables List for Algorithms");
  varTree->Branch("dEtaLL", &dEtaLL);
  varTree->Branch("dPhiLL", &dPhiLL);
  varTree->Branch("dRLL", &dEtaLL);
  varTree->Branch("MLL", &MLL);
  
  // Statistic variable
  int SelectedEvents = 0;

  // Open Event Loop
  for(int evtCtr=0; evtCtr<tree->GetEntries() /*&& evtCtr<10*/; evtCtr++) {

    if(evtCtr%100000==0) cout<<tree->GetEntries()<<" total. Ongoing event: "<<evtCtr<<endl; 
    
    tree->GetEntry(evtCtr);
    double lenObj = PID->size();

    // Select for good leptons
    bool foundGoodLep = false;
    int firstPos=-1;
    int secondPos=-1;

    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;
      //if(TMath::Abs(D0->at(objCtr))<0.1 || TMath::Abs(D0->at(objCtr))>0.2) continue;
      //if(TMath::Abs(D0->at(objCtr))<0.1) continue;
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
    }
    
    if(foundGoodLep==false) continue;

    SelectedEvents++;

    TLorentzVector lep1;
    lep1.SetPtEtaPhiE(PT->at(firstPos), Eta->at(firstPos), Phi->at(firstPos), E->at(firstPos));
    TLorentzVector lep2;
    lep2.SetPtEtaPhiE(PT->at(secondPos), Eta->at(secondPos), Phi->at(secondPos), E->at(secondPos));

    dEtaLL = TMath::Abs(lep1.Eta()-lep2.Eta());
    dPhiLL = lep1.DeltaPhi(lep2);
    dRLL = lep1.DeltaR(lep2);
    MLL = TMath::Abs((lep1+lep2).M());
    //cout<<SelectedEvents<<"\t"<<dEtaLL<<"\t"<<dPhiLL<<"\t"<<dRLL<<"\t"<<MLL<<endl;

    varTree->Fill();
  }
  
  varOutFile->Write();

  return SelectedEvents;
}

void execute() {

  TChain *sigChain = new TChain("SelectedObjects");
  sigChain->Add("./Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root");

  TChain *bkgChain = new TChain("SelectedObjects");
  bkgChain->Add("./Data/ppTobb_Cuts2/Objects_sorted_ppTobb_*.root");

  cout<<createVarOutTree(sigChain, "sigVar.root")<<" events selected from "<<sigChain->GetEntries()<<" SIGNAL events'"<<endl;
  cout<<createVarOutTree(bkgChain, "bkgVar.root")<<" events selected from "<<bkgChain->GetEntries()<<" SIGNAL events'"<<endl;

}
