#include <cstring>
#include <vector>

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/IteratorRange.h"
#include "HepMC/GenRanges.h"

#include "TH1F.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h" 
#include "TLegend.h"
#include "TText.h"

// Define class
class CompareHepmc {
public:
  // Default Constructor
  CompareHepmc(int);

  // Function to process a single event
  void processEvent(std::string, int);

  // plot Designer
  void plotDesigner(std::vector<TH1F*>, std::vector<TString>, TString, TString, TString);

  // Class Variables
  std::vector<TH1F*> hist_E;
  std::vector<TH1F*> hist_PT;
  std::vector<TH1F*> hist_Eta;
  std::vector<TH1F*> hist_Phi;
  std::vector<TH1F*> hist_D0;
  std::vector<TH1F*> hist_Log10D0;
};

// Constructor Definition
CompareHepmc::CompareHepmc(int fileTypeNum) {

  for(int fileTypeCtr=0; fileTypeCtr<fileTypeNum; fileTypeCtr++) {

    hist_E.push_back(new TH1F("E", "", 51, -1, 101));
    hist_PT.push_back(new TH1F("PT", "", 51, -1, 101));
    hist_Eta.push_back(new TH1F("Eta", "", 51, -5, 5));
    hist_Phi.push_back(new TH1F("Phi", "", 51, -3.2, 3.2));
    hist_D0.push_back(new TH1F("D0", "", 51, -0.1, 10.1));
    hist_Log10D0.push_back(new TH1F("Log10D0", "", 100, -5, 2));
  }
  
}

// Plot Designer Definition
void CompareHepmc::plotDesigner(std::vector<TH1F*> hist, std::vector<TString> label, TString XaxisTitle, TString YaxisTitle, TString saveName) {

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
    hist[histCtr]->DrawNormalized("hist SAME");
  }
  
  TLegend* legc1;
  legc1 = new TLegend(0.5, 0.9, 0.89, 1.0, NULL, "brNDC");
  for(int histCtr=0; histCtr<hist.size(); histCtr++) {
    legc1->AddEntry(hist[histCtr], label[histCtr], "l");
  }
  legc1->SetTextSize(0.03);
  legc1->SetBorderSize(0);
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

// Process Event Definition
void CompareHepmc::processEvent(std::string fName, int num) {

  int icount = 0;
  HepMC::IO_GenEvent ascii_in(fName, std::ios::in);

  HepMC::GenEvent* evt = ascii_in.read_next_event();
  while ( evt ) {

    icount++;
    if ( icount%500==1 ) std::cout << "Processing Event Number " << icount
				   << " its # " << evt->event_number() 
				   << std::endl;

    // Loop on all the particles in the event
    for ( HepMC::GenEvent::particle_const_iterator p 
	    = evt->particles_begin(); p != evt->particles_end(); ++p ){

      // Select Leptons in the event
      if ( (*p)->status() == 1
	   && (TMath::Abs((*p)->pdg_id()) == 11 || TMath::Abs((*p)->pdg_id()) == 13)
	   && TMath::Abs((*p)->momentum().pseudoRapidity()) < 2.4) {
	
	double e = (*p)->momentum().e();
	double pt = (*p)->momentum().perp();
	double eta = (*p)->momentum().pseudoRapidity();
	double phi = (*p)->momentum().phi();
	
	// Calculating the d0
	HepMC::GenVertex* v = (*p)->production_vertex();
	double vx = v->point3d().x();
	double vy = v->point3d().y();
	double px = (*p)->momentum().x();
	double py = (*p)->momentum().y();
	double d0 = vx*py-vy*px;
	d0 = d0/((*p)->momentum().perp());
	double log10d0 = TMath::Log10(TMath::Abs(d0));
	
	hist_E[num]->Fill(e);
	hist_PT[num]->Fill(pt);
	hist_Eta[num]->Fill(eta);
	hist_Phi[num]->Fill(phi);
	hist_D0[num]->Fill(TMath::Abs(d0));
	hist_Log10D0[num]->Fill(log10d0);
      }
    }
    
    delete evt;
    ascii_in >> evt;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// MAIN FUNCTION /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {

  std::vector<TString> histLabel;

  histLabel.push_back("AR Sample");
  histLabel.push_back("Nastya Sample");

  CompareHepmc *CH = new CompareHepmc(histLabel.size());

  std::string filePath1 = "soft_llp_ST_BP2_AR.hepmc";
  CH->processEvent(filePath1, 0);
  std::cout<<"Filled corresponding "<<histLabel[0]<<" histograms"<<std::endl;
  std::string filePath2 = "soft_llp_ST_BP2_Nastya.hepmc";
  CH->processEvent(filePath2, 1);
  std::cout<<"Filled corresponding "<<histLabel[1]<<" histograms"<<std::endl;

  // Regular Jet Kinematics
  CH->plotDesigner(CH->hist_E, histLabel, "E (GeV)", "Events (Scaled to 1)", "E");
  CH->plotDesigner(CH->hist_PT, histLabel, "PT (GeV)", "Number of Events scaled to unity", "PT");
  CH->plotDesigner(CH->hist_Eta, histLabel, "#eta", "Events (Scaled to 1)", "Eta");
  CH->plotDesigner(CH->hist_Phi, histLabel, "#phi", "Events (Scaled to 1)", "Phi");
  CH->plotDesigner(CH->hist_D0, histLabel, "D0 (mm)", "Events (Scaled to 1)", "D0");
  CH->plotDesigner(CH->hist_Log10D0, histLabel, "lepton Log10D0 (mm)", "Events (Scaled to 1)", "Log10D0");

  return -1;
}
