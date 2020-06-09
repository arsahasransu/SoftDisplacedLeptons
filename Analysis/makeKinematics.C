class plotVarDisb_Objects {
public:
  // Define Methods
  plotVarDisb_Objects(int);
  ~plotVarDisb_Objects();
  void readTreeFillHist(TChain*, int);
  
  // Regular Lepton Kinematics
  std::vector<TH1F*> hist_D0;
  std::vector<TH1F*> hist_D0_El;
  std::vector<TH1F*> hist_D0_Mu;
  std::vector<TH1F*> hist_PT;
  std::vector<TH1F*> hist_PT_El;
  std::vector<TH1F*> hist_PT_Mu;
  std::vector<TH1F*> hist_Eta;
  std::vector<TH1F*> hist_Eta_El;
  std::vector<TH1F*> hist_Eta_Mu;
  std::vector<TH1F*> hist_Phi;
  std::vector<TH1F*> hist_Phi_El;
  std::vector<TH1F*> hist_Phi_Mu;
  std::vector<TH1F*> hist_Iso;
  std::vector<TH1F*> hist_Iso_El;
  std::vector<TH1F*> hist_Iso_Mu;
  std::vector<TH1F*> hist_PT_Jet;
  std::vector<TH1F*> hist_MET;
    
};

// Define the constructor
plotVarDisb_Objects::plotVarDisb_Objects(int fileTypeNum) {

  // Define Histograms
  for(int fileTypeCtr=0; fileTypeCtr<fileTypeNum; fileTypeCtr++) {

    hist_D0.push_back(new TH1F("D0(cm)","",10000,0,10));
    hist_D0_El.push_back(new TH1F("D0_El(cm)","",10000,0,10));
    hist_D0_Mu.push_back(new TH1F("D0_Mu(cm)","",10000,0,10));
    hist_PT.push_back(new TH1F("PT(GeV)","",10000,0,1000));
    hist_PT_El.push_back(new TH1F("PT_El(GeV)","",10000,0,10));
    hist_PT_Mu.push_back(new TH1F("PT_Mu(GeV)","",25,20,70));
    hist_Eta.push_back(new TH1F("Eta","",51,-2.6,2.6));
    hist_Eta_El.push_back(new TH1F("Eta_El","",51,-2.6,2.6));
    hist_Eta_Mu.push_back(new TH1F("Eta_Mu","",51,-2.6,2.6));
    hist_Phi.push_back(new TH1F("Phi","",51,-3.2,3.2));
    hist_Phi_El.push_back(new TH1F("Phi_El","",51,-3.2,3.2));
    hist_Phi_Mu.push_back(new TH1F("Phi_Mu","",51,-3.2,3.2));
    hist_Iso.push_back(new TH1F("Iso","",13,0,0.25));
    hist_Iso_El.push_back(new TH1F("Iso_El","",13,0,0.25));
    hist_Iso_Mu.push_back(new TH1F("Iso_Mu","",13,0,0.25));
    hist_PT_Jet.push_back(new TH1F("PT_Jet(GeV)","",10,20,40));
    hist_MET.push_back(new TH1F("MET(GeV)","",36,0,180));
  }
}

// Read from Tree and Fill Histogram
void plotVarDisb_Objects::readTreeFillHist(TChain *tree, int fileTypeCtr) {

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

  // Statistic variable
  int SelectedEvents = 0;

  // Open Event Loop
  for(int evtCtr=0; evtCtr<tree->GetEntries(); evtCtr++) {

    if(evtCtr%100000==0) cout<<tree->GetEntries()<<" total. Ongoing event: "<<evtCtr<<endl; 
    
    tree->GetEntry(evtCtr);
    double lenObj = PID->size();

    // Select for good leptons
    bool foundGoodLep = false;
    int firstPos=-1;
    int secondPos=-1;
    int numGoodLep = 0;
    int numGoodEl = 0;
    int numGoodMu = 0;

    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;
      if(TMath::Abs(D0->at(objCtr))>100) continue;

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
    TLorentzVector METvec;
    METvec.SetPtEtaPhiE(MET, MET_Eta, MET_Phi, MET);
    lepvec.push_back(&lep1);
    lepvec.push_back(&lep2);

    hist_PT[fileTypeCtr]->Fill(PT->at(firstPos));
    hist_Eta[fileTypeCtr]->Fill(Eta->at(firstPos));
    hist_Phi[fileTypeCtr]->Fill(Phi->at(firstPos));
    hist_Iso[fileTypeCtr]->Fill(Iso->at(firstPos));
    hist_PT[fileTypeCtr]->Fill(PT->at(secondPos));
    hist_Eta[fileTypeCtr]->Fill(Eta->at(secondPos));
    hist_Phi[fileTypeCtr]->Fill(Phi->at(secondPos));
    hist_Iso[fileTypeCtr]->Fill(Iso->at(secondPos));

    if(TMath::Abs(PID->at(firstPos))==11) {
      hist_PT_El[fileTypeCtr]->Fill(PT->at(firstPos));
      hist_Eta_El[fileTypeCtr]->Fill(Eta->at(firstPos));
      hist_Phi_El[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_Iso_El[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_PT_Mu[fileTypeCtr]->Fill(PT->at(secondPos));
      hist_Eta_Mu[fileTypeCtr]->Fill(Eta->at(secondPos));
      hist_Phi_Mu[fileTypeCtr]->Fill(Phi->at(secondPos));
      hist_Iso_Mu[fileTypeCtr]->Fill(Iso->at(secondPos));
    }

    if(TMath::Abs(PID->at(firstPos))==13) {
      hist_PT_Mu[fileTypeCtr]->Fill(PT->at(firstPos));
      hist_Eta_Mu[fileTypeCtr]->Fill(Eta->at(firstPos));
      hist_Phi_Mu[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_Iso_Mu[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_PT_El[fileTypeCtr]->Fill(PT->at(secondPos));
      hist_Eta_El[fileTypeCtr]->Fill(Eta->at(secondPos));
      hist_Phi_El[fileTypeCtr]->Fill(Phi->at(secondPos));
      hist_Iso_El[fileTypeCtr]->Fill(Iso->at(secondPos));
    }
    hist_MET[fileTypeCtr]->Fill(MET);

    // Fill Histogram for good jet when a good event is found based on Leptonic selection
    int numGoodJet = 0;
    int jetFirstPos = -1;
    int jetSecondPos = -1;
    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;
      
      hist_PT_Jet[fileTypeCtr]->Fill(JetPT->at(jetCtr));
      numGoodJet++;

      if(jetFirstPos==-1) {
	jetFirstPos = jetCtr;
	continue;
      }
      if(jetFirstPos!=-1 && jetSecondPos==-1) {
	jetSecondPos = jetCtr;
	continue;
      }

    }

  } // End of Event Loop
  std::cout<<"Selected Events: "<<SelectedEvents<<" : out of Total Entries: "<<tree->GetEntries()<<std::endl;
}

void makeKinematics() {

  std::vector<TString> histLabel;
  std::vector<TChain*> t;
  std::vector<TString> dataPath;
  std::vector<int> histColor;
  std::vector<int> histLineStyle;

  histLabel.push_back("HF Background");
  histLabel.push_back("DM: (220, 20)");
  histLabel.push_back("DM: (324, 20)");
  histLabel.push_back("(220, 20, 0.2)");
  histLabel.push_back("(220, 20)");
  histLabel.push_back("(220, 20, 20)");
  histLabel.push_back("(220, 20, 200)");
  histLabel.push_back("(220, 40)");

  histColor.push_back(1);
  histColor.push_back(kBlue);
  histColor.push_back(kBlue-7);
  //histColor.push_back(6);
  histColor.push_back(kRed);
  //histColor.push_back(28);
  //histColor.push_back(36);
  histColor.push_back(kPink+5);

  histLineStyle.push_back(1);
  histLineStyle.push_back(1);
  histLineStyle.push_back(1);
  //histLineStyle.push_back(1);
  histLineStyle.push_back(1);
  //histLineStyle.push_back(1);
  //histLineStyle.push_back(1);
  histLineStyle.push_back(1);

  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_DM_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_304_324_DM/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_304_324_DM_Batch*.root");
  //dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2mm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2mm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2cm_Batch*.root");
  //dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_20cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_20cm_Batch*.root");
  //dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2m/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2m_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_180_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_180_220_2cm_Batch*.root");

  for(int ctr=0; ctr<histLabel.size(); ctr++) {
    t.push_back(new TChain("SelectedObjects"));
    t[ctr]->Add(dataPath[ctr]);
    cout<<"Initialized "<<histLabel[ctr]<<" chain"<<endl;
    std::cout<<"No.of  "<<histLabel[ctr]<<" Entries: "<<t[ctr]->GetEntries()<<std::endl;
  }

  plotVarDisb_Objects *pVDO = new plotVarDisb_Objects(histLabel.size());
  pVDO->histSuff = histLabel;
  cout<<"Initialized instance of plotter class"<<endl;

  for(int ctr=0; ctr<histLabel.size(); ctr++) {
    pVDO->readTreeFillHist(t[ctr], ctr);
    cout<<"Filled corresponding "<<histLabel[ctr]<<" histograms"<<endl;
  }
}
