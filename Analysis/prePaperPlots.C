class plotVarDisb_Objects {
public:
  // Define Methods
  plotVarDisb_Objects(int);
  ~plotVarDisb_Objects();
  void readTreeFillHist(TChain*, int);
  void plotBeautifier(std::vector<TH1F*>, std::vector<TString>,
		      TString, TString, TString,
		      std::vector<int>, std::vector<int>, 
		      bool, bool);

  // Declare histogram suffix
  std::vector<TString> histSuff;
  
  // Regular Lepton Kinematics
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
    TString temp = std::to_string(fileTypeCtr);
    histSuff.push_back(temp);

    int nD0bins = 50;
    double D0xlow = 0.001; // in cm
    double D0xup = 100; // in cm
    double D0xdiff = (TMath::Log10(D0xup)-TMath::Log10(D0xlow))/nD0bins;
    double D0xbins[nD0bins+1];
    for(int ctr=0; ctr<=nD0bins; ctr++) {
      D0xbins[ctr] = TMath::Power(10,TMath::Log10(D0xlow)+D0xdiff*ctr);
    }
    hist_PT.push_back(new TH1F("PT"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_PT_El.push_back(new TH1F("PT_El"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_PT_Mu.push_back(new TH1F("PT_Mu"+histSuff[fileTypeCtr],"",51,-1,101));
    hist_Eta.push_back(new TH1F("Eta"+histSuff[fileTypeCtr],"",51,-2.6,2.6));
    hist_Eta_El.push_back(new TH1F("Eta_El"+histSuff[fileTypeCtr],"",51,-2.6,2.6));
    hist_Eta_Mu.push_back(new TH1F("Eta_Mu"+histSuff[fileTypeCtr],"",51,-2.6,2.6));
    hist_Phi.push_back(new TH1F("Phi"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_Phi_El.push_back(new TH1F("Phi_El"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_Phi_Mu.push_back(new TH1F("Phi_Mu"+histSuff[fileTypeCtr],"",51,-3.2,3.2));
    hist_Iso.push_back(new TH1F("Iso"+histSuff[fileTypeCtr],"",51,0,1));
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
  for(int evtCtr=0; evtCtr<10000/*tree->GetEntries()*/; evtCtr++) {

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

void plotVarDisb_Objects::plotBeautifier(std::vector<TH1F*> hist, std::vector<TString> label,
					 TString XaxisTitle, TString YaxisTitle, TString saveName,
					 std::vector<int> histColor,
					 std::vector<int> histLineStyle, 
					 bool ylog=true, bool xlog=false) {

  TCanvas *c1 = new TCanvas("c1", "c1", 10,32,782,600);
  TPad *pMain = new TPad("","",0.09,0.08,1.0,0.95);
  TPad *pLeg = new TPad("","",0.001,0.9,1.0,1.0);
  TPad *pyAxisTitle = new TPad("","",0.001,0.08,0.08,0.85);
  TPad *pxAxisTitle = new TPad("","",0.08,0.001,1.0,0.08);
  pMain->Draw();
  pyAxisTitle->Draw();
  pxAxisTitle->Draw();
  pLeg->Draw();

  pMain->cd();
  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    hist[histCtr]->SetLineColor(histColor[histCtr]);
    hist[histCtr]->SetLineWidth(3);
    hist[histCtr]->SetLineStyle(histLineStyle[histCtr]);
    hist[histCtr]->Scale(1.0/hist[histCtr]->Integral());
    hist[histCtr]->GetYaxis()->SetRangeUser(0.002,0.2);
  }
  hist[0]->GetXaxis()->SetTitle("");
  hist[0]->GetYaxis()->SetTitle("");
  hist[0]->GetYaxis()->SetTicks("+");
  hist[0]->GetYaxis()->SetLabelOffset(-0.05);
  hist[0]->GetYaxis()->SetLabelFont(42);
  hist[0]->GetYaxis()->SetLabelSize(0.055);
  hist[0]->GetXaxis()->SetTicks("-");
  hist[0]->GetXaxis()->SetLabelOffset(-0.07);
  hist[0]->GetXaxis()->SetLabelFont(42);
  hist[0]->GetXaxis()->SetLabelSize(0.055);
  
  if(ylog) pMain->SetLogy();
  if(xlog) pMain->SetLogx();
  gStyle->SetOptStat(0);

  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    hist[histCtr]->Draw("hist SAME E");
  }

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
  tUFlw->SetTextSize(0.045);
  tUFlw->Draw();

  TText *t = new TText(x2-bw2/2,y2,"Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.045);
  t->Draw();
  
  pxAxisTitle->cd();
  auto xaxisTitle = new TLatex(0.7,0.35,XaxisTitle);
  xaxisTitle->SetTextAngle(0);
  xaxisTitle->SetTextSize(0.5);
  xaxisTitle->SetTextFont(42);
  xaxisTitle->Draw();

  pyAxisTitle->cd();
  auto yaxisTitle = new TLatex(0.9,0.4,YaxisTitle);
  yaxisTitle->SetTextAngle(90);
  yaxisTitle->SetTextSize(0.38);
  yaxisTitle->SetTextFont(42);
  yaxisTitle->Draw();
  
  pLeg->cd();
  
  TLegend* legc1;
  legc1 = new TLegend(0, 0.5, 0.25, 1.0, NULL, "brNDC");
  for(int histCtr=0; histCtr<1/*hist.size()/2.0*/; histCtr++) {
    legc1->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc1->SetTextSize(0.4);
  legc1->SetTextFont(42);
  legc1->SetBorderSize(0);
  //legc1->SetMargin(false);
  legc1->Draw();

  /*
  TLatex lt;
  lt.SetTextFont(42);
  lt.SetTextSize(0.045);
  lt.DrawLatex(0.01,0.101,"BP syntax (m_{c},#Deltam,c#tau_{c})");
  
  TPaveText *pt2 = new TPaveText(0.003,0.08,0.007,0.11);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->Draw();
  */
  auto abbrvFormat = new TLatex(0.01,0.15,"Signal (m_{c}, #Deltam, c#scale[1.2]{#tau}_{c})");
  abbrvFormat->SetTextAngle(0);
  abbrvFormat->SetTextSize(0.4);
  abbrvFormat->SetTextFont(42);
  abbrvFormat->Draw();

  float legStart = 0.3;
  float legDiff = 0.22;
  
  TLegend* legc2;
  legc2 = new TLegend(legStart, 0, legStart+legDiff, 1.0, NULL, "brNDC");
  for(int histCtr=1/*(int)(hist.size()/2.0)*/; histCtr<3/*hist.size()*/; histCtr++) {
    legc2->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc2->SetTextSize(0.4);
  legc2->SetTextFont(42);
  legc2->SetBorderSize(0);
  //legc2->SetMargin(false);
  legc2->Draw();
  
  TLegend* legc3;
  legc3 = new TLegend(legStart+legDiff, 0, legStart+2*legDiff, 1.0, NULL, "brNDC");
  for(int histCtr=3/*(int)(hist.size()/2.0)*/; histCtr<5/*hist.size()*/; histCtr++) {
    legc3->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc3->SetTextSize(0.4);
  legc3->SetTextFont(42);
  legc3->SetBorderSize(0);
  //legc3->SetMargin(false);
  legc3->Draw();

  TLegend* legc4;
  legc4 = new TLegend(legStart+2*legDiff, 0, legStart+3*legDiff, 1.0, NULL, "brNDC");
  for(int histCtr=5/*(int)(hist.size()/2.0)*/; histCtr<hist.size(); histCtr++) {
    legc4->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc4->SetTextSize(0.4);
  legc4->SetTextFont(42);
  legc4->SetBorderSize(0);
  //legc4->SetMargin(false);
  legc4->Draw();

  c1->SaveAs("./Analysis/PaperPlots/"+saveName+".pdf");
  delete c1;
}

void prePaperPlots() {

  std::vector<TString> histLabel;
  std::vector<TChain*> t;
  std::vector<TString> dataPath;
  std::vector<int> histColor;
  std::vector<int> histLineStyle;

  histLabel.push_back("HF Background");
  histLabel.push_back("(220, 20, DM)");
  histLabel.push_back("(324, 20, DM)");
  //histLabel.push_back("(220, 20, 0.2)");
  histLabel.push_back("(220, 20, x)");
  //histLabel.push_back("(220, 20, 20)");
  //histLabel.push_back("(220, 20, 200)");
  histLabel.push_back("(220, 40, 2)");

  histColor.push_back(1);
  histColor.push_back(2);
  histColor.push_back(4);
  //histColor.push_back(6);
  histColor.push_back(8);
  //histColor.push_back(28);
  //histColor.push_back(36);
  histColor.push_back(28);

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
  dataPath.push_back("/home/arsahasransu/Documents/SoftDisplacedLeptons/Data/DislacedLepton/Objects_sorted_DisplacedLepton_*.root");
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

  // Regular Lepton Kinematics
  pVDO->plotBeautifier(pVDO->hist_PT, histLabel, "pT (GeV)", "# Events (arbitrary units)", "PT", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_PT_El, histLabel, "Electron pT (GeV)", "# Events (arbitrary units)", "PT_El", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_PT_Mu, histLabel, "Muon pT (GeV)", "# Events (arbitrary units)", "PT_Mu", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Eta, histLabel, "#eta", "# Events (arbitrary units)", "Eta", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Eta_El, histLabel, "Electron #eta", "# Events (arbitrary units)", "Eta_El", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Eta_Mu, histLabel, "Muon #eta", "# Events (arbitrary units)", "Eta_Mu", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Phi, histLabel, "#phi", "# Events (arbitrary units)", "Phi", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Phi_El, histLabel, "Electron #phi", "# Events (arbitrary units)", "Phi_El", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Phi_Mu, histLabel, "Muon #phi", "# Events (arbitrary units)", "Phi_Mu", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Iso, histLabel, "Iso", "# Events (arbitrary units)", "Iso", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Iso_El, histLabel, "Electron Iso", "# Events (arbitrary units)", "Iso_El", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_Iso_Mu, histLabel, "Muon Iso", "# Events (arbitrary units)", "Iso_Mu", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_PT_Jet, histLabel, "Jet pT (GeV)", "# Events (arbitrary units)", "PT_Jet", histColor, histLineStyle);
  pVDO->plotBeautifier(pVDO->hist_MET, histLabel, "#slash{E}_{T} (GeV)", "# Events (arbitrary units)", "MET", histColor, histLineStyle);

}
