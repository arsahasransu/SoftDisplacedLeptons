class plotVarDisb_Objects {
public:
  // Define Methods
  plotVarDisb_Objects(int);
  void readTreeFillHist(TChain*, int);
  void plotBeautifier(std::vector<TH1F*>, std::vector<TString>, TString, TString, TString);

  // Declare histogram suffix
  std::vector<TString> histSuff;
  
  // Regular Lepton Kinematics
  std::vector<TH1F*> hist_PT_El;
  std::vector<TH1F*> hist_PT_Mu;
  std::vector<TH1F*> hist_Eta_El;
  std::vector<TH1F*> hist_Eta_Mu;
  std::vector<TH1F*> hist_Phi_El;
  std::vector<TH1F*> hist_Phi_Mu;
  std::vector<TH1F*> hist_D0_El;
  std::vector<TH1F*> hist_D0_Mu;
  std::vector<TH1F*> hist_Log10D0_El;
  std::vector<TH1F*> hist_Log10D0_Mu;
  std::vector<TH1F*> hist_Iso_El;
  std::vector<TH1F*> hist_Iso_Mu;
  std::vector<TH1F*> hist_PT_Jet;
  std::vector<TH1F*> hist_MET;
    
};


// Define the constructor
plotVarDisb_Objects::plotVarDisb_Objects(int fileTypeNum) {

  // Define Histograms
  for(int fileTypeCtr=0; fileTypeCtr<fileTypeNum; fileTypeCtr++) {
    TString temp = "";
    histSuff.push_back(temp);
    
    hist_PT_El.push_back(new TH1F("PT_El"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_PT_Mu.push_back(new TH1F("PT_Mu"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_Eta_El.push_back(new TH1F("Eta_El"+histSuff[fileTypeCtr],"",51,-2.6,2.6));
    hist_Eta_Mu.push_back(new TH1F("Eta_Mu"+histSuff[fileTypeCtr],"",51,-2.6,2.6));
    hist_Phi_El.push_back(new TH1F("Phi_El"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_Phi_Mu.push_back(new TH1F("Phi_Mu"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_D0_El.push_back(new TH1F("D0_El"+histSuff[fileTypeCtr],"",51,0,100));
    hist_D0_Mu.push_back(new TH1F("D0_Mu"+histSuff[fileTypeCtr],"",51,0,100));
    hist_Log10D0_El.push_back(new TH1F("Log10D0_El"+histSuff[fileTypeCtr],"",51,-4,3));
    hist_Log10D0_Mu.push_back(new TH1F("Log10D0_Mu"+histSuff[fileTypeCtr],"",51,-4,3));
    hist_Iso_El.push_back(new TH1F("Iso_El"+histSuff[fileTypeCtr],"",51,0,1));
    hist_Iso_Mu.push_back(new TH1F("Iso_Mu"+histSuff[fileTypeCtr],"",51,0,1));
    hist_PT_Jet.push_back(new TH1F("PT_Jet"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_MET.push_back(new TH1F("MET"+histSuff[fileTypeCtr],"",51,-1,101));
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
    if(TMath::Abs(PID->at(firstPos))==11 || TMath::Abs(PID->at(firstPos))==13) {
      hist_PT_El[fileTypeCtr]->Fill(PT->at(firstPos));
      hist_Eta_El[fileTypeCtr]->Fill(Eta->at(firstPos));
      hist_Phi_El[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_D0_El[fileTypeCtr]->Fill(TMath::Abs(D0->at(firstPos)));
      hist_Log10D0_El[fileTypeCtr]->Fill(TMath::Log10(TMath::Abs(D0->at(firstPos))));
      hist_Iso_El[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_PT_Mu[fileTypeCtr]->Fill(PT->at(secondPos));
      hist_Eta_Mu[fileTypeCtr]->Fill(Eta->at(secondPos));
      hist_Phi_Mu[fileTypeCtr]->Fill(Phi->at(secondPos));
      hist_D0_Mu[fileTypeCtr]->Fill(TMath::Abs(D0->at(secondPos)));
      hist_Log10D0_Mu[fileTypeCtr]->Fill(TMath::Log10(TMath::Abs(D0->at(secondPos))));
      hist_Iso_Mu[fileTypeCtr]->Fill(Iso->at(secondPos));
    }

    if(TMath::Abs(PID->at(firstPos))==11 || TMath::Abs(PID->at(firstPos))==13) {
      hist_PT_Mu[fileTypeCtr]->Fill(PT->at(firstPos));
      hist_Eta_Mu[fileTypeCtr]->Fill(Eta->at(firstPos));
      hist_Phi_Mu[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_D0_Mu[fileTypeCtr]->Fill(TMath::Abs(D0->at(firstPos)));
      hist_Log10D0_Mu[fileTypeCtr]->Fill(TMath::Log10(TMath::Abs(D0->at(firstPos))));
      hist_Iso_Mu[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_PT_El[fileTypeCtr]->Fill(PT->at(secondPos));
      hist_Eta_El[fileTypeCtr]->Fill(Eta->at(secondPos));
      hist_Phi_El[fileTypeCtr]->Fill(Phi->at(secondPos));
      hist_D0_El[fileTypeCtr]->Fill(TMath::Abs(D0->at(secondPos)));
      hist_Log10D0_El[fileTypeCtr]->Fill(TMath::Log10(TMath::Abs(D0->at(secondPos))));
      hist_Iso_El[fileTypeCtr]->Fill(Iso->at(secondPos));
    }
    hist_MET[fileTypeCtr]->Fill(MET);
    //cout<<D0->at(firstPos)<<"\t"<<D0->at(secondPos)<<endl;
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

void plotVarDisb_Objects::plotBeautifier(std::vector<TH1F*> hist, std::vector<TString> label, TString XaxisTitle, TString YaxisTitle, TString saveName) {

  TCanvas *c1 = new TCanvas("c1", "c1", 10,32,782,552);

  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    hist[histCtr]->SetLineColor(histCtr+1);
    hist[histCtr]->SetLineWidth(3);
    hist[histCtr]->GetXaxis()->SetRange(0, hist[histCtr]->GetNbinsX()+1);
  }
  hist[0]->GetXaxis()->SetTitle(XaxisTitle);
  hist[0]->GetYaxis()->SetTitle(YaxisTitle);

  c1->SetLogy();
  gStyle->SetOptStat(0);

  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    hist[histCtr]->DrawNormalized("hist SAME E");
  }

  TLegend* legc1;
  legc1 = new TLegend(0.1, 0.9, 0.89, 1.0, NULL, "brNDC");
  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    legc1->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc1->SetTextSize(0.03);
  legc1->SetBorderSize(0);
  //legc1->SetMargin(false);
  legc1->Draw();

  
  Int_t nx    = hist[0]->GetNbinsX()+1;
  Double_t bw1 = hist[0]->GetBinWidth(0);
  Double_t bw2 = hist[0]->GetBinWidth(nx);
  Double_t x1 = hist[0]->GetBinLowEdge(0)+bw1;
  Double_t x2 = hist[0]->GetBinLowEdge(nx)+bw2;
  Double_t y1 = hist[0]->GetBinContent(0)/hist[0]->Integral();
  Double_t y2 = hist[0]->GetBinContent(nx)/hist[0]->Integral();
  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    y1 = y1>(hist[histCtr]->GetBinContent(0)/hist[histCtr]->Integral())
      ?y1
      :(hist[histCtr]->GetBinContent(0)/hist[histCtr]->Integral());
    y2 = y2>(hist[histCtr]->GetBinContent(nx)/hist[histCtr]->Integral())
      ?y2
      :(hist[histCtr]->GetBinContent(nx)/hist[histCtr]->Integral());
  }
  y1 = y1*1.1;
  y2 = y2*1.1;
  TText *tUFlw = new TText(x1-bw1/2,y1,"Underflow");
  tUFlw->SetTextAngle(90);
  tUFlw->SetTextAlign(12);
  tUFlw->SetTextSize(0.03);
  tUFlw->Draw();

  TText *t = new TText(x2-bw2/2,y2,"Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);
  t->Draw();
  
  c1->SaveAs("./Analysis/baseKin/"+saveName+".pdf");
}

void execute() {

  std::vector<TString> histLabel;
  std::vector<TChain*> t;
  std::vector<TString> dataPath;

  
  histLabel.push_back("pp > bb~ HF Background");
  histLabel.push_back("m_{c}=324 GeV, #Deltam=20 GeV, c#tau_{c}=2 cm");
  histLabel.push_back("m_{c}=200 GeV, #Deltam=20 GeV, c#tau_{c}=20 cm");
  histLabel.push_back("m_{c}=200 GeV, #Deltam=20 GeV, c#tau_{c}=2 cm");
  histLabel.push_back("m_{c}=200 GeV, #Deltam=20 GeV, c#tau_{c}=2 mm");
  histLabel.push_back("m_{c}=180 GeV, #Deltam=40 GeV, c#tau_{c}=2 cm");

  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_20cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_20cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2cm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2cm_Batch*.root");
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DisplacedModel_BP_200_220_2mm/ObjectSorted_Delp341_Mg5v266_PY8243_DisplacedModel_BP_200_220_2mm_Batch*.root");
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

  // Regular Lepton Kinematics
  pVDO->plotBeautifier(pVDO->hist_PT_El, histLabel, "pT_{e}", "Number of Events scaled to unity", "PT_El");
  pVDO->plotBeautifier(pVDO->hist_PT_Mu, histLabel, "pT_{#mu}", "Number of Events scaled to unity", "PT_Mu");
  pVDO->plotBeautifier(pVDO->hist_Eta_El, histLabel, "#eta_{e}", "Number of Events scaled to unity", "Eta_El");
  pVDO->plotBeautifier(pVDO->hist_Eta_Mu, histLabel, "#eta_{#mu}", "Number of Events scaled to unity", "Eta_Mu");
  pVDO->plotBeautifier(pVDO->hist_Phi_El, histLabel, "#phi_{e}", "Number of Events scaled to unity", "Phi_El");
  pVDO->plotBeautifier(pVDO->hist_Phi_Mu, histLabel, "#phi_{#mu}", "Number of Events scaled to unity", "Phi_Mu");
  pVDO->plotBeautifier(pVDO->hist_Iso_El, histLabel, "Iso_{e}", "Number of Events scaled to unity", "Iso_El");
  pVDO->plotBeautifier(pVDO->hist_Iso_Mu, histLabel, "Iso_{#mu}", "Number of Events scaled to unity", "Iso_Mu");
  pVDO->plotBeautifier(pVDO->hist_D0_El, histLabel, "d0_{e}", "Number of Events scaled to unity", "D0_El");
  pVDO->plotBeautifier(pVDO->hist_D0_Mu, histLabel, "d0_{#mu}", "Number of Events scaled to unity", "D0_Mu");
  pVDO->plotBeautifier(pVDO->hist_Log10D0_El, histLabel, "log10d0_{e}", "Number of Events scaled to unity", "Log10D0_El");
  pVDO->plotBeautifier(pVDO->hist_Log10D0_Mu, histLabel, "log10d0_{#mu}", "Number of Events scaled to unity", "Log10D0_Mu");
  pVDO->plotBeautifier(pVDO->hist_PT_Jet, histLabel, "pT_(jet)", "Number of Events scaled to unity", "PT_Jet");
  pVDO->plotBeautifier(pVDO->hist_MET, histLabel, "#slash{E}_{T}", "Number of Events scaled to unity", "MET");

}
