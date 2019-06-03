//  Plotting macro for out/HTjets/pAu_HTjets.root
//  Veronica Verkest     June 2, 2019

#include "TPad.h"
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
#include <stdlib.h>

using namespace std;

int main() {

  const double twopi = 2*3.14159265358979;
  TString Ndj, avg;
  
  TFile* inFile = new TFile( "out/HTjets/pAu_HTjets.root", "READ" );
  TString lpf = "lpf";

  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hPrimaryVsGlobal = new TH2D("hPrimaryVsGlobal","# Primary Tracks vs. # Global Tracks;# Global Tracks;# Primary Tracks", 150,0,3000, 150,0,150 );
  TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 50,0.0,6.3, 50,-1.0,1.0);
  TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 50,0.0,6.3, 50,-1.0,1.0);
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);#rho (GeV);BBC East Rate", 50,0.0,25, 10,0,10, 150,0,15000000);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 80,0,35, 100,0,1000);
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 70,0.05,35, 70,0,70);

  TH2D *hscale = new TH2D( "hscale", "Underlying Event by Lead Jet p_{T};#rho (GeV);", 60,0,30, 10,0.0001, 1.0 );
  // hscale->GetYaxis()->SetRangeUser( 0.000001, 1 );
  
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
    jt->GetEntry(i);

    hPrimaryVsGlobal->Fill(nGlobal,nPrimary);
    hPrimaryVsBBCE->Fill(BbcEastRate,nPrimary);
    hGlobalVsBBCE->Fill(BbcEastRate,nGlobal);

    if ( (abs(leadPhi) > twopi) || (abs(subPhi) > twopi) ) { cout<<"phi range is funky   "<<leadPhi<<"  "<<subPhi<<endl; }
    
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
  TH2D * hUE_BBCE[nPtBins];
  //hPt_UE_BBCE->Scale(1./nEntries);
  
  int ptBinLo[nPtBins] = { 0, 5, 10, 15, 25, 35, 50 };
  int ptBinHi[nPtBins] = { 5, 10, 15, 25, 35, 50, 70 };
  TString ptBinString[nPtBins] = { "<5 GeV", "5-10 GeV", "10-15 GeV", "15-25 GeV",  "25-35 GeV", "35-50 GeV", ">50 GeV" };
  TString ptBinName[nPtBins] = { "_5", "_5_10", "_10_15", "_15_25", "_25_35", "_35_50", "_50" };
  int color[nPtBins] = { 1, 633, 613, 596, 839, 414, 797 };
  int marker[nPtBins] = { 29, 33, 34, 23, 22, 21, 20 };
  TString name, title;

  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c1->SetLogy();
  hscale->SetStats(0);
  hscale->Draw();

  TLegend *leg1 = new TLegend(0.65, 0.68, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  TLegendEntry *entry;
  leg1->SetNColumns(3);
  leg1->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg1->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg1->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadPtVsRho" + ptBinName[i];
    title = ptBinString[i];
    hRho[i] = hLeadPtVsRho->ProjectionX( name, ptBinLo[i], ptBinHi[i] );       // PROJECT
    hRho[i]->SetStats(0);
    hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
    hRho[i]->SetLineColor( color[i] );
    hRho[i]->SetMarkerStyle( marker[i] );
    hRho[i]->SetMarkerColor( color[i] );
    hRho[i]->Draw("SAME");                                                    // DRAW
    /*entry=*/
    Ndj = ""; avg = "";
    Ndj += hRho[i]->GetEntries();
    avg += hRho[i]->GetMean(1);                                           // 1 denotes x-axis
    leg1->AddEntry( name, title, lpf );                            // ADD TO LEGEND
    leg1->AddEntry((TObject*)0,Ndj, "");
    leg1->AddEntry((TObject*)0,avg, "");
    lpf += "lpf";
  }

  leg1->Draw();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
  c1->SaveAs("plots/RhoByLeadPt.pdf","PDF");


  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c0->SetLogz();
  
  for ( int i=0; i<nPtBins; ++i ) {
    name = "plots/UEvsBBCE" + ptBinName[i] + ".pdf";
    title = "Underlying Event vs. BBC East Rate - p_{T}^{lead}:" + ptBinString[i];
    hPt_UE_BBCE->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hUE_BBCE[i] = (TH2D*)hPt_UE_BBCE->Project3D( "yz" );       // PROJECT
    hUE_BBCE[i]->Scale( 1./hUE_BBCE[i]->GetEntries() );                     // NORMALIZE
    hUE_BBCE[i]->SetTitle(title);
    hUE_BBCE[i]->Draw("COLZ");
    c0->SaveAs( name,"PDF");
  }
  
  gStyle->SetOptStat(1);
  TCanvas * c2 = new TCanvas( "c2" , "" ,0 ,23 ,1280 ,700 );              // CANVAS

  double norm = 1./nEntries;
  
  hLeadEtaPhi->Scale( 1./(double)hLeadEtaPhi->GetEntries() );
  hLeadEtaPhi->Draw("COLZ");
  c2->SaveAs("plots/LeadEtaPhi.pdf","PDF");

  hSubEtaPhi->Scale( 1./(double)hSubEtaPhi->GetEntries() );
  hSubEtaPhi->Draw("COLZ");
  c2->SaveAs("plots/SubEtaPhi.pdf","PDF");

  hTowersVsRho->Scale( 1./(double)hTowersVsRho->GetEntries() );
  hTowersVsRho->Draw("COLZ");
  c2->SaveAs("plots/TowersVsRho.pdf","PDF");
  
  c2->SetLogz();
  hPrimaryVsGlobal->Scale( 1./(double)hPrimaryVsGlobal->GetEntries() );
  hPrimaryVsGlobal->Draw("COLZ");
  c2->SaveAs("plots/PrimaryVsGlobal.pdf","PDF");

  outFile->Write();
  outFile->Close();
  inFile->Close();

  return 0;
}
