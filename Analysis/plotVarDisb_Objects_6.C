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
  std::vector<TH1F*> hist_Iso_Lep;
  std::vector<TH1F*> hist_Iso_El;
  std::vector<TH1F*> hist_Iso_Mu;
  std::vector<TH1F*> hist_leadIso_Lep;
  std::vector<TH1F*> hist_leadIso_El;
  std::vector<TH1F*> hist_leadIso_Mu;
  std::vector<TH1F*> hist_subLeadIso_Lep;
  std::vector<TH1F*> hist_subLeadIso_El;
  std::vector<TH1F*> hist_subLeadIso_Mu;
  std::vector<TH1F*> hist_diffIso_Lep;
    
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
    hist_Iso_Lep.push_back(new TH1F("Iso_Lepton"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_Iso_El.push_back(new TH1F("Iso_El"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_Iso_Mu.push_back(new TH1F("Iso_Mu"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_leadIso_Lep.push_back(new TH1F("leadIso_Lepton"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_leadIso_El.push_back(new TH1F("leadIso_El"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_leadIso_Mu.push_back(new TH1F("leadIso_Mu"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_subLeadIso_Lep.push_back(new TH1F("subLeadIso_Lepton"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_subLeadIso_El.push_back(new TH1F("subLeadIso_El"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_subLeadIso_Mu.push_back(new TH1F("subLeadIso_Mu"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
    hist_diffIso_Lep.push_back(new TH1F("diffIso_Lept"+histSuff[fileTypeCtr],"",51,-0.01,1.01));
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

      hist_Iso_Lep[fileTypeCtr]->Fill(Iso->at(objCtr));
      if(TMath::Abs(PID->at(objCtr))==11) hist_Iso_Mu[fileTypeCtr]->Fill(Iso->at(objCtr));
      if(TMath::Abs(PID->at(objCtr))==13) hist_Iso_El[fileTypeCtr]->Fill(Iso->at(objCtr));
    }

    hist_leadIso_Lep[fileTypeCtr]->Fill(Iso->at(firstPos));
    hist_subLeadIso_Lep[fileTypeCtr]->Fill(Iso->at(secondPos));
    hist_diffIso_Lep[fileTypeCtr]->Fill(TMath::Abs(Iso->at(firstPos)-Iso->at(secondPos)));
    if(TMath::Abs(PID->at(firstPos))==11) {
      hist_leadIso_El[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_subLeadIso_Mu[fileTypeCtr]->Fill(Iso->at(secondPos));      
    }
    if(TMath::Abs(PID->at(firstPos))==13) {
      hist_leadIso_Mu[fileTypeCtr]->Fill(Iso->at(firstPos));
      hist_subLeadIso_El[fileTypeCtr]->Fill(Iso->at(secondPos));      
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
  pVDO->plotBeautifier(pVDO->hist_Iso_Lep, histLabel, "#Iso_{lep}", "Number of Events scaled to unity", "Iso_Lep");
  pVDO->plotBeautifier(pVDO->hist_Iso_El, histLabel, "#Iso_{e}", "Number of Events scaled to unity", "Iso_El");
  pVDO->plotBeautifier(pVDO->hist_Iso_Mu, histLabel, "#Iso_{#mu}", "Number of Events scaled to unity", "Iso_Mu");
  pVDO->plotBeautifier(pVDO->hist_leadIso_Lep, histLabel, "#Iso^{lead}_{lep}", "Number of Events scaled to unity", "Iso_LeadLep");
  pVDO->plotBeautifier(pVDO->hist_leadIso_El, histLabel, "#Iso^{lead}_{e}", "Number of Events scaled to unity", "Iso_LeadEl");
  pVDO->plotBeautifier(pVDO->hist_leadIso_Mu, histLabel, "#Iso^{lead}_{#mu}", "Number of Events scaled to unity", "Iso_LeadMu");
  pVDO->plotBeautifier(pVDO->hist_subLeadIso_Lep, histLabel, "#Iso^{sub-lead}_{lep}", "Number of Events scaled to unity", "Iso_SubLeadLep");
  pVDO->plotBeautifier(pVDO->hist_subLeadIso_El, histLabel, "#Iso^{sub-lead}_{e}", "Number of Events scaled to unity", "Iso_SubLeadEl");
  pVDO->plotBeautifier(pVDO->hist_subLeadIso_Mu, histLabel, "#Iso^{sub-lead}_{#mu}", "Number of Events scaled to unity", "Iso_SubLeadMu");
  pVDO->plotBeautifier(pVDO->hist_diffIso_Lep, histLabel, "#Delta#Iso_{lep}", "Number of Events scaled to unity", "DeltaIso_Lep");

}
