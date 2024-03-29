void stackHistos(std::vector<TH1D*> sig,
		 TH1D* bkg,
		 std::vector<int> histColor,
		 TString saveName,
		 bool weight_d0=true) {

  TCanvas *c1 = new TCanvas("c1", "c1", 10,32,782,600);
  TPad *pMain = new TPad("","",0.09,0.08,1.0,0.95);
  TPad *pLeg = new TPad("","",0.001,0.88,1.0,0.98);
  TPad *pyAxisTitle = new TPad("","",0.001,0.08,0.08,0.85);
  TPad *pxAxisTitle = new TPad("","",0.08,0.001,1.0,0.08);
  pMain->Draw();
  pyAxisTitle->Draw();
  pxAxisTitle->Draw();
  pLeg->Draw();
  gStyle->SetOptStat(0);
  
  pMain->cd();
  pMain->SetLogy();

  int bkgnbinsInit = bkg->GetNbinsX();
  bkg->SetBinContent(bkgnbinsInit-1, bkg->GetBinContent(bkgnbinsInit-1)+bkg->GetBinContent(bkgnbinsInit));
  bkg->SetBinError(bkgnbinsInit-1, bkg->GetBinError(bkgnbinsInit-1)+bkg->GetBinError(bkgnbinsInit));
  bkg->SetBins(100,0.0,1.0);
  int bkgnbin = 10;
  bkg->Rebin(bkgnbin);
  bkg->SetBinContent(bkgnbin+1, 0);
  bkg->SetBinError(bkgnbin+1, 0);
  bkg->SetFillColor(kGray);
  bkg->SetLineColor(kBlack);
  bkg->GetXaxis()->SetTicks("-");
  bkg->GetXaxis()->SetLabelOffset(-0.08);
  bkg->GetXaxis()->SetLabelFont(42);
  bkg->GetXaxis()->SetLabelSize(0.055);
  bkg->GetYaxis()->SetTicks("+");
  bkg->GetYaxis()->SetLabelOffset(-0.05);
  bkg->GetYaxis()->SetLabelFont(42);
  bkg->GetYaxis()->SetLabelSize(0.055);
  bkg->Draw("hist");

  for(int sigCtr=0; sigCtr<sig.size(); sigCtr++){
    int nbinsInit = sig[sigCtr]->GetNbinsX();
    sig[sigCtr]->SetBinContent(nbinsInit-1, sig[sigCtr]->GetBinContent(nbinsInit-1)+sig[sigCtr]->GetBinContent(nbinsInit));
    sig[sigCtr]->SetBinError(nbinsInit-1, sig[sigCtr]->GetBinError(nbinsInit-1)+sig[sigCtr]->GetBinError(nbinsInit));
    sig[sigCtr]->SetBins(100,0.0,1.0);
    int nBin = 10;
    sig[sigCtr]->Rebin(nBin);
    sig[sigCtr]->Add(bkg);
    sig[sigCtr]->SetBinContent(nBin+1, 0);
    sig[sigCtr]->SetBinError(nBin+1, 0);
    std::cout<<bkg->GetNbinsX()<<"\t"<<sig[sigCtr]->GetNbinsX()<<std::endl;
    sig[sigCtr]->SetFillStyle(0);
    sig[sigCtr]->SetLineWidth(3);
    sig[sigCtr]->SetLineColor(histColor[sigCtr]);
    sig[sigCtr]->Draw("same hist");
  }

  pxAxisTitle->cd();
  auto xaxisTitle = new TLatex(0.6,0.3,"signal probability");
  xaxisTitle->SetTextAngle(0);
  xaxisTitle->SetTextSize(0.6);
  xaxisTitle->SetTextFont(42);
  xaxisTitle->Draw();

  pyAxisTitle->cd();
  auto yaxisTitle = new TLatex(0.7,0.5,"Events");
  yaxisTitle->SetTextAngle(90);
  yaxisTitle->SetTextSize(0.45);
  yaxisTitle->SetTextFont(42);
  yaxisTitle->Draw();

  pLeg->cd();
  TLegend *leg = new TLegend(0.02,0,0.3,1.0,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TLegendEntry *entry=leg->AddEntry("D00","HF background","f");
  entry->SetFillColor(kGray);
  entry->SetFillStyle(1001);
  entry->SetLineColor(kBlack);
  entry->SetLineWidth(1);
  entry=leg->AddEntry("D03","(200, 20, 0.2)","l");
  entry->SetLineColor(kRed-9);
  entry->SetLineWidth(3);
  leg->Draw();
  
  float legStart = 0.365;
  float legDiff = 0.19;
  
  leg = new TLegend(legStart,0,legStart+legDiff,1.0,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  entry=leg->AddEntry("D01","DM: (220, 20, 3)","l");
  entry->SetLineColor(kBlue);
  entry->SetLineWidth(3);
  entry=leg->AddEntry("D04","(200, 20, 2)","l");
  entry->SetLineColor(kRed-7);
  entry->SetLineWidth(3);
  leg->Draw();

  leg = new TLegend(legStart+legDiff,0,legStart+2*legDiff,1.0,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  entry=leg->AddEntry("D02","DM: (324, 20, 2)","l");
  entry->SetLineColor(kBlue-7);
  entry->SetLineWidth(3);
  entry=leg->AddEntry("D05","(200, 20, 20)","l");
  entry->SetLineColor(kRed-4);
  entry->SetLineWidth(3);
  leg->Draw();

  leg = new TLegend(legStart+2*legDiff,0,legStart+3*legDiff,1.0,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  entry=leg->AddEntry("D07","(220, 40, 2)","l");
  entry->SetLineColor(kPink+7);
  entry->SetLineWidth(3);
  entry=leg->AddEntry("D06","(200, 20, 200)","l");
  entry->SetLineColor(kRed+2);
  entry->SetLineWidth(3);
  leg->Draw();

  c1->SaveAs(saveName+".C");
}

void stackByYield_trial2() {

  double normBGSR1 = 4122.8;
  double normBGSR2 = 644.2;
  double normBGSR3 = 24.479;

  //std::vector<TH1D*>
  std::vector<TH1D*> sig_hist;
  TH1D* bkg_hist = new TH1D();
  std::vector<TString> fileName;
  std::vector<int> histColor;

  histColor.push_back(kBlue);
  histColor.push_back(kBlue-7);
  histColor.push_back(kRed-9);
  histColor.push_back(kRed-7);
  histColor.push_back(kRed-4);
  histColor.push_back(kRed+2);
  histColor.push_back(kPink+7);

  fileName.push_back("BP_200_20_DM_Disc.root");
  fileName.push_back("BP_324_20_DM_Disc.root");
  fileName.push_back("BP_200_20_02_Disc.root");
  fileName.push_back("BP_200_20_2_Disc.root");
  fileName.push_back("BP_200_20_20_Disc.root");
  fileName.push_back("BP_200_20_200_Disc.root");
  fileName.push_back("BP_200_40_20_Disc.root");

  for(int sigCtr=0; sigCtr<fileName.size(); sigCtr++) {

    TFile *sig_file = TFile::Open(fileName[sigCtr],"READ");

    sig_hist.push_back((TH1D*) sig_file->Get("SR1"));
    bkg_hist = (TH1D*) sig_file->Get("background");

    bkg_hist->Scale(normBGSR1/bkg_hist->Integral());

    std::cout<<sig_hist[sigCtr]->Integral()<<std::endl;

    //delete sig_file;
  }

  stackHistos(sig_hist, bkg_hist, histColor, "SR1", false);
  stackHistos(sig_hist, bkg_hist, histColor, "SR1_d0weight", true);
}
