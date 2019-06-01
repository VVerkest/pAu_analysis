//  Plotting macro for out/HTjets/pAu_HTjets.root
//  Veronica Verkest     May 31, 2019

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main() {

  TFile* inFile = new TFile( "out/HTjets/pAu_HTjets.root", "READ" );

  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH1D *hTowerMult = new TH1D("hTowerMult","Tower Multiplicity by ID;Tower ID", 2400,0,4800);
  TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 180,0,6.3, 100,-1.0,1.0);
  TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 180,0,6.3, 100,-1.0,1.0);
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 300,0,300, 40,0,200, 3500,0,3500000);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 80,0,35, 100,0,1000);
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 80,0,35, 140,0,70);
  
  TTree *jt = HTjetTree;
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, sigma;
  
  jt->SetBranchAddress("EventID", &EventID);   jt->SetBranchAddress("RunID", &RunID);   jt->SetBranchAddress("Vz", &Vz);
  jt->SetBranchAddress("nTowers", &nTowers);    jt->SetBranchAddress("nPrimary", &nPrimary);   jt->SetBranchAddress("nGlobal", &nGlobal);
  jt->SetBranchAddress("nVertices", &nVertices);     jt->SetBranchAddress("refMult", &refMult);     jt->SetBranchAddress("gRefMult", &gRefMult);
  jt->SetBranchAddress("BbcCoincidenceRate", &BbcCoincidenceRate);     jt->SetBranchAddress("BbcEastRate", &BbcEastRate);
  jt->SetBranchAddress("BbcWestRate", &BbcWestRate);    jt->SetBranchAddress("vpdVz", &vpdVz);    jt->SetBranchAddress("nJets", &nJets);
  jt->SetBranchAddress("leadPt", &leadPt);    jt->SetBranchAddress("leadEta", &leadEta);    jt->SetBranchAddress("leadPhi", &leadPhi);
  jt->SetBranchAddress("leadEt", &leadEt);    jt->SetBranchAddress("leadNcons", &leadNcons);     jt->SetBranchAddress("subPt", &subPt);
  jt->SetBranchAddress("subEta", &subEta);    jt->SetBranchAddress("subPhi", &subPhi);     jt->SetBranchAddress("subEt", &subEt);
  jt->SetBranchAddress("subNcons", &subNcons);    HTjetTree->SetBranchAddress("rho", &rho);     HTjetTree->SetBranchAddress("sigma", &sigma);
  nEntries=jt->GetEntries();

  for (int i=0; i<nEntries; ++i){
    jt->GetEntry(i);
    
    hPrimaryVsBBCE->Fill(nPrimary,BbcEastRate);
    hGlobalVsBBCE->Fill(nGlobal,BbcEastRate);
    hTowerMult->Fill(towID);
    hLeadEtaPhi->Fill(leadPhi,leadEta);
    hSubEtaPhi->Fill(subPhi,subEta);
    hPt_UE_BBCE->Fill(leadPt,rho,BbcEastRate);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadPtVsRho->Fill(rho,leadPt);
  }
  
  hPrimaryVsBBCE->Write();
  hGlobalVsBBCE->Write();
  hTowerMult->Write();
  hLeadEtaPhi->Write();
  hSubEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hTowersVsRho->Write();
  hLeadPtVsRho->Write();
    
  inFile->Close();

  TFile *outFile = new TFile( "plots/pAu_HT_jetPlot.root" ,"RECREATE");

  outFile->Write();
  
  return 0;
}
