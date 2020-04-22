#include <string>
#include <cstring>
#include <vector>


#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void delphesObjectSelector(int n)
{
  gSystem->Load("libDelphes");
  
  // Define ROOT File and Tree
  TFile outFile("ObjectSelected_Delp342_MG270_Py8244_DisplacedModel_BP_200_220_2cm_Nishita.root","RECREATE");
  TTree event("SelectedObjects","Selected Objects in the Event");

  // Define the variables and associated tree branches
  std::vector<double> simuE, PT, Eta, Phi;
  std::vector<double> Iso, D0, Dz;
  std::vector<int> PID;
  double MET, MET_Eta, MET_Phi, HT;
  std::vector<double> Jet_PT, Jet_Eta, Jet_Phi, Jet_M, Jet_Flavor;
  int nEl, nMu, nJet;

  auto branchE = event.Branch("simulatedE",&simuE);
  auto branchPT = event.Branch("PT",&PT);
  auto branchEta = event.Branch("Eta",&Eta);
  auto branchPhi = event.Branch("Phi",&Phi);
  auto branchIso = event.Branch("Iso",&Iso);
  auto branchD0 = event.Branch("D0",&D0);
  auto branchDz = event.Branch("Dz",&Dz);
  auto branchPID = event.Branch("PID",&PID);
  auto branchNEl = event.Branch("NumElec",&nEl,"nEl/I");
  auto branchNMu = event.Branch("NumMuon",&nMu,"nMu/I");
  auto branchMet = event.Branch("MET",&MET,"MET/D");
  auto branchMET_Eta = event.Branch("MET_Eta",&MET_Eta,"MET_Eta/D");
  auto branchMET_Phi = event.Branch("MET_Phi",&MET_Phi,"MET_Phi/D");
  auto branchJet_PT = event.Branch("JetPT",&Jet_PT);
  auto branchJet_Eta = event.Branch("JetEta",&Jet_Eta);
  auto branchJet_Phi = event.Branch("JetPhi",&Jet_Phi);
  auto branchJet_M = event.Branch("JetM",&Jet_M);
  auto branchJet_Flavor = event.Branch("JetFlavor",&Jet_Flavor);
  auto branchHt = event.Branch("HT",&HT,"HT/D");
  auto branchNJet = event.Branch("NumJet",&nJet,"NJet/I");

  // Create chain of root trees
  int batchCtr, fileCtr;
  TChain chain("Delphes");
  chain.Add("../../Data/NewSimu_TestData/Delp342_MG270_Py8244_DisplacedModel_BP_200_220_2cm_Nishita.root");
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  TClonesArray *branchHT = treeReader->UseBranch("ScalarHT");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  std::vector<Electron *> elec_sel;
  std::vector<Muon *> muon_sel;
  std::vector<Jet *> jet_sel;
  std::vector<GenParticle *> part_ElSel;
  std::vector<Track *> trk_ElSel;
  std::vector<GenParticle *> part_MuSel;
  std::vector<Track *> trk_MuSel;

  std::vector<pair <int, char>> sortOrdrLep;
  std::vector<pair <double, bool>> LepPtPair;
  
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Clear the arrays before the next event
    simuE.clear();
    PT.clear();
    Eta.clear();
    Phi.clear();
    Iso.clear();
    D0.clear();
    Dz.clear();
    PID.clear();
    Jet_PT.clear();
    Jet_Eta.clear();
    Jet_Phi.clear();
    Jet_M.clear();
    Jet_Flavor.clear();

    elec_sel.clear();
    muon_sel.clear();
    jet_sel.clear();
    part_ElSel.clear();
    part_MuSel.clear();
    trk_ElSel.clear();
    trk_MuSel.clear();
    sortOrdrLep.clear();
    LepPtPair.clear();

    if(entry%10000==9999) std::cout<<entry<<" done"<<std::endl;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    int NumElec = branchElectron->GetEntries(), elecCtr=0;
    int NumMuon = branchMuon->GetEntries(), muonCtr=0;
    int NumGenPart = branchParticle->GetEntries(), partCtr=0;
    int NumTrkPart = branchTrack->GetEntries(), trkCtr=0;
    int NumJet = branchJet->GetEntries(), jetCtr=0;

    // Select Leptons
    for(elecCtr=0; elecCtr<NumElec; elecCtr++) {
      Electron *elec = (Electron*) branchElectron->At(elecCtr);
      if(elec->PT>15) elec_sel.push_back(elec);
    }

    for(muonCtr=0; muonCtr<NumMuon; muonCtr++) {
      Muon *muon = (Muon*) branchMuon->At(muonCtr);
      if(muon->PT>15)	muon_sel.push_back(muon);
    }

    nEl = elec_sel.size();
    nMu = muon_sel.size();
    
    // Select Jets                                                                                                                 
    for(jetCtr=0; jetCtr<NumJet; jetCtr++) {
      Jet *jet = (Jet*) branchJet->At(jetCtr);
      if(jet->PT>15) jet_sel.push_back(jet);
    }

    nJet = jet_sel.size();

    // Sort Jets by PT                                                                                                             
    for(jetCtr=0; jetCtr<nJet; jetCtr++) {
      for(int jetCtr2=jetCtr; jetCtr2<nJet; jetCtr2++) {
        if(jet_sel[jetCtr]->PT<jet_sel[jetCtr2]->PT) {
          Jet *temp_jet = jet_sel[jetCtr];
          jet_sel[jetCtr] = jet_sel[jetCtr2];
          jet_sel[jetCtr2] = temp_jet;
        }
      }
    }

    if(nEl+nMu < 1) continue;

    // Sort Leptons by PT
    for(elecCtr=0; elecCtr<elec_sel.size(); elecCtr++) {
      LepPtPair.push_back(make_pair(elec_sel[elecCtr]->PT,true));
    }
    for(muonCtr=0; muonCtr<muon_sel.size(); muonCtr++) {
      LepPtPair.push_back(make_pair(muon_sel[muonCtr]->PT,true));
    }
    for(int lepCtr=0; lepCtr<LepPtPair.size(); lepCtr++) {
      double maxPt = -1;
      int nPos;
      char part;

      for(int lepCtr2=0; lepCtr2<LepPtPair.size(); lepCtr2++) {
	if(LepPtPair[lepCtr2].second==false) continue;
	if(maxPt<0) {
	  maxPt = LepPtPair[lepCtr2].first;
	  nPos = lepCtr2>=nEl?lepCtr2-nEl:lepCtr2;
	  part = lepCtr2>=nEl?'u':'e';
	  continue;
	}
	if(maxPt < LepPtPair[lepCtr2].first) {
	  maxPt = LepPtPair[lepCtr2].first;
	  nPos = lepCtr2>=nEl?lepCtr2-nEl:lepCtr2;
	  part = lepCtr2>=nEl?'u':'e';	  
	}
      }
      if(part=='u') LepPtPair[nPos+nEl].second = false;
      if(part=='e') LepPtPair[nPos].second = false;

      sortOrdrLep.push_back(make_pair(nPos,part));
    }

    if(nEl+nMu>2) std::cout<<sortOrdrLep.size()<<LepPtPair.size()<<nEl<<nMu<<"\t";
    // Fill Leptons
    for(int lepCtr=0; lepCtr<sortOrdrLep.size(); lepCtr++) {
      if(sortOrdrLep[lepCtr].second=='e') {
	auto gen = (GenParticle*) (elec_sel[sortOrdrLep[lepCtr].first]->Particle).GetObject();
	//      if(gen->D0 > 10) continue;
	D0.push_back(gen->D0);
	Dz.push_back(gen->DZ);
	PID.push_back(gen->PID);
	simuE.push_back(gen->E);
	PT.push_back(elec_sel[sortOrdrLep[lepCtr].first]->PT);
	if(nEl+nMu>2) std::cout<<elec_sel[sortOrdrLep[lepCtr].first]->PT<<"\t";
	Eta.push_back(elec_sel[sortOrdrLep[lepCtr].first]->Eta);
	Phi.push_back(elec_sel[sortOrdrLep[lepCtr].first]->Phi);
	Iso.push_back(elec_sel[sortOrdrLep[lepCtr].first]->IsolationVar);
      }
      if(sortOrdrLep[lepCtr].second=='u') {
	PT.push_back(muon_sel[sortOrdrLep[lepCtr].first]->PT);
	if(nEl+nMu>2) std::cout<<muon_sel[sortOrdrLep[lepCtr].first]->PT<<"\t";
	Eta.push_back(muon_sel[sortOrdrLep[lepCtr].first]->Eta);
	Phi.push_back(muon_sel[sortOrdrLep[lepCtr].first]->Phi);
	Iso.push_back(muon_sel[sortOrdrLep[lepCtr].first]->IsolationVar);
	
	auto gen = (GenParticle*) (muon_sel[sortOrdrLep[lepCtr].first]->Particle).GetObject();
	D0.push_back(gen->D0);
	Dz.push_back(gen->DZ);
	PID.push_back(gen->PID);
	simuE.push_back(gen->E);
      }
    }
    if(nEl+nMu>2) std::cout<<std::endl;
      
    // Fill MET data of the Event
    MissingET* eventMET = (MissingET*) branchMET->At(0);
    MET = eventMET->MET;
    MET_Eta = eventMET->Eta;
    MET_Phi = eventMET->Phi;

    ScalarHT* scalarHT = (ScalarHT*) branchHT->At(0);
    HT = scalarHT->HT;

    // File Jet Vectors                                                                                                            
    for(jetCtr=0; jetCtr<nJet; jetCtr++) {
      Jet_PT.push_back(jet_sel[jetCtr]->PT);
      Jet_Eta.push_back(jet_sel[jetCtr]->Eta);
      Jet_Phi.push_back(jet_sel[jetCtr]->Phi);
      Jet_M.push_back(jet_sel[jetCtr]->Mass);
      Jet_Flavor.push_back(jet_sel[jetCtr]->Flavor);
    }

    // Fill one event in the output tree
    event.Fill();
    
  } // End of the event loop

  outFile.Write();
  //std::cout<<numberOfEntries<<std::endl;

}
/*
void delphesObjectSelector() {

  for(int i=101; i<=101; i++) {
    std::cout<<"Processing File: "<<i<<std::endl;
    delphesObjectSelectorSingleFile(i);
  }
}
*/
