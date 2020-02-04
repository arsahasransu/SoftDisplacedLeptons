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
  std::vector<TH1F*> hist_MET;
  std::vector<TH1F*> hist_METLepJet;
  std::vector<TH1F*> hist_METLep;
  std::vector<TH1F*> hist_METJet;
  std::vector<TH1F*> hist_METDiffMETLepJet;
    
  std::vector<TH1F*> hist_MET_Phi;
  std::vector<TH1F*> hist_METLepJet_Phi;
  std::vector<TH1F*> hist_METLep_Phi;
  std::vector<TH1F*> hist_METJet_Phi;
  std::vector<TH1F*> hist_METDiffMETLepJet_Phi;

  std::vector<TH1F*> hist_HT;
  std::vector<TH1F*> hist_HTLepJet;
  std::vector<TH1F*> hist_HTLep;
  std::vector<TH1F*> hist_HTJet;
  std::vector<TH1F*> hist_HTDiffHTLepJet;

  std::vector<TH1F*> hist_DiffMETHT;
  std::vector<TH1F*> hist_DiffMETHTLepJet;
  std::vector<TH1F*> hist_DiffMETHTLep;
  std::vector<TH1F*> hist_DiffMETHTJet;
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
    hist_MET.push_back(new TH1F("MET"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_METLepJet.push_back(new TH1F("METLepJet"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_METLep.push_back(new TH1F("MET"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_METJet.push_back(new TH1F("MET"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_METDiffMETLepJet.push_back(new TH1F("MET"+histSuff[fileTypeCtr],"",51,-1,101));

    hist_MET_Phi.push_back(new TH1F("MET_Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_METLepJet_Phi.push_back(new TH1F("METLepJet_Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_METLep_Phi.push_back(new TH1F("MET_Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_METJet_Phi.push_back(new TH1F("MET_Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_METDiffMETLepJet_Phi.push_back(new TH1F("MET_Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));

    hist_HT.push_back(new TH1F("HT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_HTLepJet.push_back(new TH1F("HTLepJet"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_HTLep.push_back(new TH1F("HT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_HTJet.push_back(new TH1F("HT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_HTDiffHTLepJet.push_back(new TH1F("HT"+histSuff[fileTypeCtr],"",51,-1,101));

    hist_DiffMETHT.push_back(new TH1F("DiffMETHT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_DiffMETHTLepJet.push_back(new TH1F("DiffMETHTLepJet"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_DiffMETHTLep.push_back(new TH1F("DiffMETHT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_DiffMETHTJet.push_back(new TH1F("DiffMETHT"+histSuff[fileTypeCtr],"",51,-1,101));
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
    TLorentzVector lepSum, jetSum, sum, diff, metvec;
    int numGoodJet=0;
    double htlepjet=0, htlep=0, htjet=0;

    metvec.SetPtEtaPhiE(MET, 0.0, MET_Phi, MET);
    
    for(int objCtr=0; objCtr<lenObj; objCtr++) {
      if(TMath::Abs(Eta->at(objCtr))>2.4) continue;
      if(PT->at(objCtr)<20) continue;

      htlep += TMath::Abs(PT->at(objCtr));
      TLorentzVector lepSingle;
      lepSingle.SetPtEtaPhiE(PT->at(objCtr), Eta->at(objCtr), Phi->at(objCtr), E->at(objCtr));
      lepSum += lepSingle;

    }

    for(int jetCtr=0; jetCtr<numJet; jetCtr++) {
      if(JetPT->at(jetCtr)<20) continue;
      if(TMath::Abs(JetEta->at(jetCtr))>2.5) continue;
      
      numGoodJet++;
      htjet += TMath::Abs(JetPT->at(jetCtr));
      TLorentzVector jetSingle;
      jetSingle.SetPtEtaPhiM(JetPT->at(jetCtr), JetEta->at(jetCtr), JetPhi->at(jetCtr), JetM->at(jetCtr));
      jetSum += jetSingle;

    }

    sum = lepSum+jetSum;
    diff = sum-metvec;
    htlepjet = htlep+htjet;
 
    if(MET>20.0) hist_MET[fileTypeCtr]->Fill(MET);
    hist_METLepJet[fileTypeCtr]->Fill(sum.Pt());
    if(numGoodJet!=0) hist_METJet[fileTypeCtr]->Fill(jetSum.Pt());
    hist_METLep[fileTypeCtr]->Fill(lepSum.Pt());
    if(MET>20.0) hist_METDiffMETLepJet[fileTypeCtr]->Fill(TMath::Abs(diff.Pt()));

    if(MET>20.0) hist_MET_Phi[fileTypeCtr]->Fill(MET_Phi);
    hist_METLepJet_Phi[fileTypeCtr]->Fill(sum.Phi());
    if(numGoodJet!=0) hist_METJet_Phi[fileTypeCtr]->Fill(jetSum.Phi());
    hist_METLep_Phi[fileTypeCtr]->Fill(lepSum.Phi());
    if(MET>20.0) hist_METDiffMETLepJet_Phi[fileTypeCtr]->Fill(diff.Phi());

    if(HT>20.0) hist_HT[fileTypeCtr]->Fill(HT);
    hist_HTLepJet[fileTypeCtr]->Fill(htlepjet);
    if(numGoodJet!=0) hist_HTJet[fileTypeCtr]->Fill(htjet);
    hist_HTLep[fileTypeCtr]->Fill(htlep);
    if(HT>20.0) hist_HTDiffHTLepJet[fileTypeCtr]->Fill(TMath::Abs(HT-htlepjet));

    if(MET>20.0 && HT>20.0) hist_DiffMETHT[fileTypeCtr]->Fill(TMath::Abs(MET-HT));
    hist_DiffMETHTLepJet[fileTypeCtr]->Fill(TMath::Abs(sum.Pt()-htlepjet));
    if(numGoodJet!=0) hist_DiffMETHTJet[fileTypeCtr]->Fill(TMath::Abs(jetSum.Pt()-htjet));
    hist_DiffMETHTLep[fileTypeCtr]->Fill(TMath::Abs(lepSum.Pt()-htlep));

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
  pVDO->plotBeautifier(pVDO->hist_MET, histLabel, "#slash{E}_{T} (GeV)", "Number of Events scaled to unity", "MET");
  pVDO->plotBeautifier(pVDO->hist_METLepJet, histLabel, "#Sigma_{lep+jet} #vec{p_{T}} (GeV)", "Number of Events scaled to unity", "METLepJet");
  pVDO->plotBeautifier(pVDO->hist_METLep, histLabel, "#Sigma_{lep} #vec{p_{T}} (GeV)", "Number of Events scaled to unity", "METLep");
  pVDO->plotBeautifier(pVDO->hist_METJet, histLabel, "#Sigma_{jet} #vec{p_{T}} (GeV)", "Number of Events scaled to unity", "METJet");
  pVDO->plotBeautifier(pVDO->hist_METDiffMETLepJet, histLabel, "#cbar#slash{E}_{T}-#Sigma_{lep+jet} #vec{p_{T}}#cbar (GeV)", "Number of Events scaled to unity", "METDiffMETLepJet");

  pVDO->plotBeautifier(pVDO->hist_MET_Phi, histLabel, "#phi(#slash{E}_{T})", "Number of Events scaled to unity", "MET_Phi");
  pVDO->plotBeautifier(pVDO->hist_METLepJet_Phi, histLabel, "#phi(#Sigma_{lep+jet} #vec{p_{T}})", "Number of Events scaled to unity", "METLepJet_Phi");
  pVDO->plotBeautifier(pVDO->hist_METLep_Phi, histLabel, "#phi(#Sigma_{lep} #vec{p_{T}})", "Number of Events scaled to unity", "METLep_Phi");
  pVDO->plotBeautifier(pVDO->hist_METJet_Phi, histLabel, "#phi(#Sigma_{jet} #vec{p_{T}})", "Number of Events scaled to unity", "METJet_Phi");
  pVDO->plotBeautifier(pVDO->hist_METDiffMETLepJet_Phi, histLabel, "#phi(#cbar#slash{E}_{T}-#Sigma_{lep+jet} #vec{p_{T}}#cbar)", "Number of Events scaled to unity", "METDiffMETLepJet_Phi");

  pVDO->plotBeautifier(pVDO->hist_HT, histLabel, "HT (GeV)", "Number of Events scaled to unity", "HT");
  pVDO->plotBeautifier(pVDO->hist_HTLepJet, histLabel, "#Sigma_{lep+jet}#cbar p_{T}#cbar (GeV)", "Number of Events scaled to unity", "HTLepJet");
  pVDO->plotBeautifier(pVDO->hist_HTLep, histLabel, "#Sigma_{lep}#cbar p_{T}#cbar (GeV)", "Number of Events scaled to unity", "HTLep");
  pVDO->plotBeautifier(pVDO->hist_HTJet, histLabel, "#Sigma_{jet}#cbar p_{T}#cbar (GeV)", "Number of Events scaled to unity", "HTJet");
  pVDO->plotBeautifier(pVDO->hist_HTDiffHTLepJet, histLabel, "#cbar#slash{E}_{T}-(#Sigma_{lep+jet}#cbar p_{T}#cbar)#cbar (GeV)", "Number of Events scaled to unity", "HTDiffHTLepJet");

  pVDO->plotBeautifier(pVDO->hist_DiffMETHT, histLabel, "#cbar#slash{E}_{T}-HT#cbar (GeV)", "Number of Events scaled to unity", "DiffMETHT");
  pVDO->plotBeautifier(pVDO->hist_DiffMETHTLepJet, histLabel, "#cbar#Sigma_{lep+jet} #vec{p_{T}}-(#Sigma_{lep+jet}#cbar p_{T}#cbar)#cbar (GeV)", "Number of Events scaled to unity", "DiffMETHTLepJet");
  pVDO->plotBeautifier(pVDO->hist_DiffMETHTLep, histLabel, "#cbar#Sigma_{lep} #vec{p_{T}}-(#Sigma_{lep}#cbar p_{T}#cbar)#cbar (GeV)", "Number of Events scaled to unity", "DiffMETHTLep");
  pVDO->plotBeautifier(pVDO->hist_DiffMETHTJet, histLabel, "#cbar#Sigma_{jet} #vec{p_{T}}-(#Sigma_{jet}#cbar p_{T}#cbar)#cbar (GeV)", "Number of Events scaled to unity", "DiffMETHTJet");
}
