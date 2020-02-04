class plotVarDisb_Objects {
public:
  // Define Methods
  plotVarDisb_Objects(int);
  void readTreeFillHist(TChain*, int);
  void plotBeautifier(std::vector<TH1F*>, std::vector<TString>, TString, TString, TString);
  void plotBeautifier2D(std::vector<TH2F*>, std::vector<TString>, TString, TString, TString);
  static TMatrixD makeMomentumTensor3D(std::vector<TLorentzVector*>);
  static double sphericity(std::vector<TLorentzVector*>);
  static double spherocity(std::vector<TLorentzVector*>);

  // Declare histogram suffix
  std::vector<TString> histSuff;
  
  // Declare histograms
  
  // Regular Lepton Kinematics
  std::vector<TH1F*> hist_Phi_Lep;
  std::vector<TH1F*> hist_Phi_El;
  std::vector<TH1F*> hist_Phi_Mu;
  std::vector<TH1F*> hist_Phi_Jet;
  std::vector<TH1F*> hist_leadPhi_Lep;
  std::vector<TH1F*> hist_leadPhi_El;
  std::vector<TH1F*> hist_leadPhi_Mu;
  std::vector<TH1F*> hist_leadPhi_Jet;
  std::vector<TH1F*> hist_subLeadPhi_Lep;
  std::vector<TH1F*> hist_subLeadPhi_El;
  std::vector<TH1F*> hist_subLeadPhi_Mu;
  std::vector<TH1F*> hist_subLeadPhi_Jet;
  std::vector<TH1F*> hist_diffPhi_Lep;
  std::vector<TH1F*> hist_diffPhi_Jet;
    
};

TMatrixD plotVarDisb_Objects::makeMomentumTensor3D(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total(3,3);
  double normaliser = 0.0;
  TVector3 vect;

  for(int ii=0; ii<3; ii++) {
    for(int jj=0; jj<3; jj++) {
      normaliser=0.0;
      total[ii][jj]=0.0;
      for(int vec=0; vec<lvarray.size(); vec++) {
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

double plotVarDisb_Objects::sphericity(std::vector<TLorentzVector*> lvarray) {
  TMatrixD total = plotVarDisb_Objects::makeMomentumTensor3D(lvarray);

  TVectorD eigenvals;
  auto eigenvectors = total.EigenVectors(eigenvals);

  double sphericity = 1.5*(eigenvals[1]+eigenvals[2]);
  return sphericity;
}

// Define the constructor
plotVarDisb_Objects::plotVarDisb_Objects(int fileTypeNum) {

  // Define Histograms
  for(int fileTypeCtr=0; fileTypeCtr<fileTypeNum; fileTypeCtr++) {
    TString temp = "";
    histSuff.push_back(temp);
        
    // Regular Lepton Kinematics
    hist_Phi_Lep.push_back(new TH1F("Phi_Lepton"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_Phi_El.push_back(new TH1F("Phi_El"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_Phi_Mu.push_back(new TH1F("Phi_Mu"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_Phi_Jet.push_back(new TH1F("Phi_Jet"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_leadPhi_Lep.push_back(new TH1F("leadPhi_Lepton"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_leadPhi_El.push_back(new TH1F("leadPhi_El"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_leadPhi_Mu.push_back(new TH1F("leadPhi_Mu"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_leadPhi_Jet.push_back(new TH1F("leadPhi_Jet"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_subLeadPhi_Lep.push_back(new TH1F("subLeadPhi_Lepton"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_subLeadPhi_El.push_back(new TH1F("subLeadPhi_El"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_subLeadPhi_Mu.push_back(new TH1F("subLeadPhi_Mu"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_subLeadPhi_Jet.push_back(new TH1F("subLeadPhi_Jet"+histSuff[fileTypeCtr],"",51,-3.4,3.4));
    hist_diffPhi_Lep.push_back(new TH1F("diffPhi_Lept"+histSuff[fileTypeCtr],"",51,0,6.5));
    hist_diffPhi_Jet.push_back(new TH1F("diffPhi_Jet"+histSuff[fileTypeCtr],"",51,0,6.5));
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
  for(int evtCtr=0; evtCtr<tree->GetEntries() /*&& evtCtr<10*/; evtCtr++) {

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

    // Fill histograms for leptons selected above

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
    
    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;

      hist_Phi_Lep[fileTypeCtr]->Fill(Phi->at(objCtr));
      if(TMath::Abs(PID->at(objCtr))==11) hist_Phi_Mu[fileTypeCtr]->Fill(Phi->at(objCtr));
      if(TMath::Abs(PID->at(objCtr))==13) hist_Phi_El[fileTypeCtr]->Fill(Phi->at(objCtr));
    }

    hist_leadPhi_Lep[fileTypeCtr]->Fill(Phi->at(firstPos));
    hist_subLeadPhi_Lep[fileTypeCtr]->Fill(Phi->at(secondPos));
    hist_diffPhi_Lep[fileTypeCtr]->Fill(TMath::Abs(Phi->at(firstPos)-Phi->at(secondPos)));
    if(TMath::Abs(PID->at(firstPos))==11) {
      hist_leadPhi_El[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_subLeadPhi_Mu[fileTypeCtr]->Fill(Phi->at(secondPos));      
    }
    if(TMath::Abs(PID->at(firstPos))==13) {
      hist_leadPhi_Mu[fileTypeCtr]->Fill(Phi->at(firstPos));
      hist_subLeadPhi_El[fileTypeCtr]->Fill(Phi->at(secondPos));      
    }

    // Fill Histogram for good jet when a good event is found based on Leptonic selection
    int numGoodJet = 0;
    int jetFirstPos = -1;
    int jetSecondPos = -1;
    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;
      
      hist_Phi_Jet[fileTypeCtr]->Fill(JetPhi->at(jetCtr));
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

    if(numGoodJet>=1) hist_leadPhi_Jet[fileTypeCtr]->Fill(JetPhi->at(jetFirstPos));
    if(numGoodJet>=2) {
      hist_subLeadPhi_Jet[fileTypeCtr]->Fill(JetPhi->at(jetSecondPos));
      hist_diffPhi_Jet[fileTypeCtr]->Fill(JetPhi->at(jetFirstPos)-JetPhi->at(jetSecondPos));
      std::cout<<JetPhi->at(jetFirstPos)<<"\t"<<JetPhi->at(jetSecondPos)<<std::endl;
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
  legc1 = new TLegend(0.5, 0.9, 0.89, 1.0, NULL, "brNDC");
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
  
  c1->SaveAs(saveName+".pdf");
}

void plotVarDisb_Objects::plotBeautifier2D(std::vector<TH2F*> hist, std::vector<TString> label, TString XaxisTitle, TString YaxisTitle, TString saveName) {

  TCanvas *c1 = new TCanvas("c1", "c1", 10,10,782,782);

  hist[0]->GetXaxis()->SetTitle(XaxisTitle);
  hist[0]->GetYaxis()->SetTitle(YaxisTitle);

  c1->SetLogz();
  gStyle->SetOptStat(0);

  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    //hist[histCtr]->SetLineColor
    TExec *ex = new TExec("ex",("gStyle->SetPalette(53+"+std::to_string(histCtr)+");").c_str());
    hist[histCtr]->DrawNormalized("COL SAME");
    ex->Draw();
  }

  TLegend* legc1;
  legc1 = new TLegend(0.5, 0.9, 0.89, 1.0, NULL, "brNDC");
  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    legc1->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc1->SetTextSize(0.03);
  legc1->SetBorderSize(0);
  //legc1->SetMargin(false);
  legc1->Draw();

  c1->SaveAs(saveName+".pdf");
}


void execute() {

  std::vector<TString> histLabel;

  histLabel.push_back("Displaced Lepton Signal");
  TChain *t0 = new TChain("SelectedObjects");
  t0->Add("../Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root");
  cout<<"Initialized "<<histLabel[0]<<" chain"<<endl;
  
  histLabel.push_back("pp > bb~ HF Background");
  TChain *t1 = new TChain("SelectedObjects");
  t1->Add("../Data/ppTobb_Cuts2/Objects_sorted_ppTobb_Cuts2_*.root");
  cout<<"Initialized "<<histLabel[1]<<" chain"<<endl;

  std::cout<<"No.of  "<<histLabel[0]<<" Entries: "<<t0->GetEntries()<<std::endl;
  std::cout<<"No.of  "<<histLabel[1]<<" Entries: "<<t1->GetEntries()<<std::endl;

  plotVarDisb_Objects *pVDO = new plotVarDisb_Objects(histLabel.size());
  pVDO->histSuff = histLabel;
  cout<<"Initialized instance of plotter class"<<endl;

  pVDO->readTreeFillHist(t0, 0);
  cout<<"Filled corresponding "<<histLabel[0]<<" histograms"<<endl;
  pVDO->readTreeFillHist(t1, 1);
  cout<<"Filled corresponding "<<histLabel[1]<<" histograms"<<endl;
  
  // Regular Lepton Kinematics
  pVDO->plotBeautifier(pVDO->hist_Phi_Lep, histLabel, "#phi_{lep}", "Number of Events scaled to unity", "Phi_Lep");
  pVDO->plotBeautifier(pVDO->hist_Phi_El, histLabel, "#phi_{e}", "Number of Events scaled to unity", "Phi_El");
  pVDO->plotBeautifier(pVDO->hist_Phi_Mu, histLabel, "#phi_{#mu}", "Number of Events scaled to unity", "Phi_Mu");
  pVDO->plotBeautifier(pVDO->hist_Phi_Jet, histLabel, "#phi_(jet)", "Number of Events scaled to unity", "Phi_Jet");
  pVDO->plotBeautifier(pVDO->hist_leadPhi_Lep, histLabel, "#phi^{lead}_{lep}", "Number of Events scaled to unity", "Phi_LeadLep");
  pVDO->plotBeautifier(pVDO->hist_leadPhi_El, histLabel, "#phi^{lead}_{e}", "Number of Events scaled to unity", "Phi_LeadEl");
  pVDO->plotBeautifier(pVDO->hist_leadPhi_Mu, histLabel, "#phi^{lead}_{#mu}", "Number of Events scaled to unity", "Phi_LeadMu");
  pVDO->plotBeautifier(pVDO->hist_leadPhi_Jet, histLabel, "#phi^{lead}_{jet}", "Number of Events scaled to unity", "Phi_LeadJet");
  pVDO->plotBeautifier(pVDO->hist_subLeadPhi_Lep, histLabel, "#phi^{sub-lead}_{lep}", "Number of Events scaled to unity", "Phi_SubLeadLep");
  pVDO->plotBeautifier(pVDO->hist_subLeadPhi_El, histLabel, "#phi^{sub-lead}_{e}", "Number of Events scaled to unity", "Phi_SubLeadEl");
  pVDO->plotBeautifier(pVDO->hist_subLeadPhi_Mu, histLabel, "#phi^{sub-lead}_{#mu}", "Number of Events scaled to unity", "Phi_SubLeadMu");
  pVDO->plotBeautifier(pVDO->hist_subLeadPhi_Jet, histLabel, "#phi^{sub-lead}_{jet}", "Number of Events scaled to unity", "Phi_SubLeadJet");
  pVDO->plotBeautifier(pVDO->hist_diffPhi_Lep, histLabel, "#Delta#phi_{lep}", "Number of Events scaled to unity", "DeltaPhi_Lep");
  pVDO->plotBeautifier(pVDO->hist_diffPhi_Jet, histLabel, "#Delta#phi_{jet}", "Number of Events scaled to unity", "DeltaPhi_Jet");

}
