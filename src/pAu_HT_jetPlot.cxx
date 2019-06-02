//  Plotting macro for out/HTjets/pAu_HTjets.root
//  Veronica Verkest     June 2, 2019

#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
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
  TString lpf = "lpf";

  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 128,-3.2,3.2, 100,-1.0,1.0);
  TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 128,-3.2,3.2, 100,-1.0,1.0);
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 300,0,300, 40,0,200, 3500,0,3500000);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 80,0,35, 100,0,1000);
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 70,0.05,35, 70,0,70);

  TH1D *hscale = new TH1D( "hscale", "Underlying Event by Lead Jet p_{T};#rho (GeV)", 70,0,35);
  
  TTree *jt = (TTree*) inFile->Get("HTjetTree");
  
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
  jt->SetBranchAddress("subNcons", &subNcons);    jt->SetBranchAddress("rho", &rho);     jt->SetBranchAddress("sigma", &sigma);
  int nEntries=jt->GetEntries();

  for (int i=0; i<nEntries; ++i){
  // for (int i=0; i<1000; ++i){
    jt->GetEntry(i);
    
    hPrimaryVsBBCE->Fill(nPrimary,BbcEastRate);
    hGlobalVsBBCE->Fill(nGlobal,BbcEastRate);
    hLeadEtaPhi->Fill(leadPhi,leadEta);
    hSubEtaPhi->Fill(subPhi,subEta);
    hPt_UE_BBCE->Fill(leadPt,rho,BbcEastRate);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadPtVsRho->Fill(rho,leadPt);
  }
      
  TFile *outFile = new TFile( "plots/pAu_HT_jetPlot.root" ,"RECREATE");

  hPrimaryVsBBCE->Write();
  hGlobalVsBBCE->Write();
  hLeadEtaPhi->Write();
  hSubEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hTowersVsRho->Write();
  hLeadPtVsRho->Write();

  const int nPtBins = 7;

  TH1D * hRho[nPtBins];
  
  int ptBinLo[nPtBins] = { 0, 5, 10, 15, 25, 35, 50 };
  int ptBinHi[nPtBins] = { 5, 10, 15, 25, 35, 50, 70 };
  TString ptBinString[nPtBins] = { "<5 GeV", "5-10 GeV", "10-15 GeV", "15-25 GeV",  "25-35 GeV", "35-50 GeV", "<50 GeV" };
  TString ptBinName[nPtBins] = { "_5", "_5_10", "_10_15", "_15_25", "_25_35", "_35_50", "_50" };
  int color[nPtBins] = { 1, 633, 613, 596, 839, 414, 797 };
  TString name, title;

  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );       gStyle->SetOptStat(0);   // CANVAS
  hscale->Draw();

  TLegend *leg = new TLegend(0.75, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg->SetBorderSize(1);   leg->SetLineColor(1);   leg->SetLineStyle(1);   leg->SetLineWidth(1);   leg->SetFillColor(0);   leg->SetFillStyle(1001);

  
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadPtVsRho" + ptBinName[i];
    title = ptBinString[i];
    hRho[i] = hLeadPtVsRho->ProjectionX( name, ptBinLo[i], ptBinHi[i] );       // PROJECT
    hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
    hRho[i]->SetLineWidth(2);
    hRho[i]->SetLineColor( color[i] );
    hRho[i]->SetMarkerColor( color[i] );
    hRho[i]->Draw("SAME");                                                    // DRAW
    TLegendEntry *entry = leg->AddEntry( name, title, lpf );                            // ADD TO LEGEND
    entry->SetLineColor( color[i] );   entry->SetMarkerColor( color[i] );
    entry->SetFillStyle(1001);   entry->SetTextFont(42);
    entry->SetLineStyle(1);   entry->SetLineWidth(2);
    entry->SetMarkerStyle(1);   entry->SetMarkerSize(1);
    lpf += "lpf";
  }

  leg->Draw();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  c1->SaveAs("pdfs/RhoByLeadPt.pdf","PDF");
  
  outFile->Write();
  outFile->Close();
  inFile->Close();

  return 0;
}
