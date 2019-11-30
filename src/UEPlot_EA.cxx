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

void UEPlot_EA(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  gStyle->SetErrorX(0.0001);
  
  const double pi = 3.14159265;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const int ptLoVal[nPtBins] = { 10, 15, 20 };
  const int ptHiVal[nPtBins] = { 15, 20, 30 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  int eval, pval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, midRho, eastRho, westRho, rho, stdev;
  
  TString fileName = "out/UE/pAuHTjetUE.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  jetTree->SetBranchAddress( "RunID", &RunID );      	       		jetTree->SetBranchAddress( "EventID", &EventID );					jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );       		jetTree->SetBranchAddress( "nGlobal", &nGlobal );					jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );			jetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );	       			jetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );	       		jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );			jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );

  int nEntries = jetTree->GetEntries();

  const int nEAbins = 2;
  TString EAbinName[nEAbins] = { "Lo", "Hi" };
  TString BBCselection[nEAbins] = { "BbcAdcEastSum>4107 && BbcAdcEastSum<11503", "BbcAdcEastSum>28537" };

  TH1D *hRho_pt[nPtBins][nEtaBins][nEAbins];
  TH1D *hRho_eta[nPtBins][nEtaBins][nEAbins];

  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho+neuEastRho)", "(chgMidRho+neuMidRho)", "(chgWestRho+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta<-0.3", "leadEta>=-0.3 && leadEta<=0.3", "leadEta>0.3" };

  int EAcolor[nEAbins] = { 810, 884 };
  int EAmarker[nEAbins] = { 23, 22 };

  int MARKER[nEAbins][nPtBins] = { { 24, 25, 30 }, { 20, 21, 29 } };
  int COLOR[nEtaBins] = { 877, 596, 814 };

  TH1D *hRhoByEta_pt[nPtBins][nEAbins];
  TH1D *hRho[nPtBins][nEtaBins][nEtaBins][nEAbins];
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {

      name = "hRhoByEta" + EAbinName[a] + ptBinName[p];
      hRhoByEta_pt[p][a] = new TH1D( name, "", 3,-1.5,1.5 );

    }
  }
  
  TCanvas *cx = new TCanvas("cx");

  TH2D *sRhoByEta = new TH2D("sRhoByEta","", 3,-1.5,1.5, 10,0.5,1.5);
  sRhoByEta->GetXaxis()->SetLabelSize(0);
  sRhoByEta->GetYaxis()->SetLabelSize(0.06);
  sRhoByEta->GetYaxis()->SetNdivisions(10);
  sRhoByEta->SetStats(0);
  
  TString hname[nPtBins][nEtaBins][nEAbins];
  TCanvas *cpt = new TCanvas( "cpt", "", 0, 23, 300, 900 );
  cpt->Divide(1,3,0,0);
  for ( int e=0; e<nEtaBins; ++e ) { cpt->cd(e+1);  sRhoByEta->Draw(); }
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {
      for ( int e=0; e<nEtaBins; ++e ) {
	hname[p][e][a] = "h" + eastmidwest[e] + EAbinName[a] + ptBinName[p];
	TString drawString = rhoVal[e] + ">>" + hname[p][e][a];
	TString drawCuts = ptSelection[p] + " && " + BBCselection[a];
	cx->cd();
	jetTree->Draw( drawString, drawCuts, "" );
	name = hname[p][e][a];
	hRho_pt[p][e][a] = (TH1D*)gDirectory->Get( name );
	hRho_pt[p][e][a]->SetStats(0);
	hRho_pt[p][e][a]->SetStats(0);

	rho = hRho_pt[p][e][a]->GetMean();
	stdev = hRho_pt[p][e][a]->GetMeanError();
      
	hRhoByEta_pt[p][a]->SetBinContent( e+1, rho );
	hRhoByEta_pt[p][a]->SetBinError( e+1, stdev );
	hRhoByEta_pt[p][a]->SetMarkerStyle( EAmarker[a] );
	hRhoByEta_pt[p][a]->SetMarkerSize( 2 );
	hRhoByEta_pt[p][a]->SetMarkerColor( EAcolor[a] );
	hRhoByEta_pt[p][a]->SetLineColor( EAcolor[a] );	
      }
      cpt->cd();
      cpt->cd(p+1);
      gPad->SetTickx();
      gPad->SetTicky();
      gPad->SetGridy();
      hRhoByEta_pt[p][a]->Draw("PSAME");
    }
  }


  cpt->SaveAs("plots/UE/pt3plot.pdf","PDF");


  


  TH1D *hRhoByEta_eta[nPtBins][nEAbins];
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {

      name = "hRhoByEta" + EAbinName[a] + etaBinName[p];
      hRhoByEta_eta[p][a] = new TH1D( name, "", 3,-1.5,1.5 );
    }
  }

  TCanvas *ceta = new TCanvas( "ceta", "", 0, 23, 900, 300 );
  ceta->Divide(3,1,0,0);
  for ( int e=0; e<nEtaBins; ++e ) { ceta->cd(e+1);  sRhoByEta->Draw(); }
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int je=0; je<nPtBins; ++je ) {
      for ( int e=0; e<nEtaBins; ++e ) {
	hname[je][e][a] = "h" + eastmidwest[e] + EAbinName[a] + etaBinName[je];
	TString drawString = rhoVal[e] + ">>" + hname[je][e][a];
	TString drawCuts = etaSelection[je] + " && " + BBCselection[a];
	cx->cd();
	jetTree->Draw( drawString, drawCuts, "" );
	name = hname[je][e][a];
	hRho_eta[je][e][a] = (TH1D*)gDirectory->Get( name );
	hRho_eta[je][e][a]->SetStats(0);
	hRho_eta[je][e][a]->SetStats(0);

	rho = hRho_eta[je][e][a]->GetMean();
	stdev = hRho_eta[je][e][a]->GetMeanError();
      
	hRhoByEta_eta[je][a]->SetBinContent( e+1, rho );
	hRhoByEta_eta[je][a]->SetBinError( e+1, stdev );
	hRhoByEta_eta[je][a]->SetMarkerStyle( EAmarker[a] );
	hRhoByEta_eta[je][a]->SetMarkerSize( 2 );
	hRhoByEta_eta[je][a]->SetMarkerColor( EAcolor[a] );
	hRhoByEta_eta[je][a]->SetLineColor( EAcolor[a] );	
      }
      ceta->cd();
      ceta->cd(je+1);
      gPad->SetTickx();
      gPad->SetTicky();
      gPad->SetGridy();
      hRhoByEta_eta[je][a]->Draw("PSAME");
    }
  }


  ceta->SaveAs("plots/UE/eta3plot.pdf","PDF");
  



  //    WITH HI/LO EA AND JET ETA
  sRhoByEta->GetYaxis()->SetRangeUser( 0.5, 1.5 );
  
  TString hname2[nPtBins][nEtaBins][nEtaBins][nEAbins];
  TH1D *hRhoByEta_jetEta__pt[nPtBins][nEAbins][nEtaBins];
  
  TCanvas *cpt2 = new TCanvas( "cpt2", "", 0, 23, 300, 900 );
  cpt2->Divide(1,3,0,0);
  for ( int e=0; e<nEtaBins; ++e ) { cpt2->cd(e+1);  sRhoByEta->Draw(); }

  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {
      for ( int je=0; je<nEtaBins; ++je ) {
	name = "hRho" + EAbinName[a] + jetEtaBinName[je] + ptBinName[p];
	hRhoByEta_jetEta__pt[p][a][je] = new TH1D( name, "", 3,-1.5,1.5 );
	hRhoByEta_jetEta__pt[p][a][je]->SetMarkerStyle( MARKER[a][p] );
	hRhoByEta_jetEta__pt[p][a][je]->SetMarkerColor( COLOR[je] );
	hRhoByEta_jetEta__pt[p][a][je]->SetLineColor( COLOR[je] );
      }
    }
  }
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {
      for ( int je=0; je<nEtaBins; ++je ) {
	for ( int e=0; e<nEtaBins; ++e ) {

	  
	  hname2[p][je][e][a] = "h" + eastmidwest[e] + EAbinName[a] + jetEtaBinName[je] + ptBinName[p];
	  name = hname2[p][je][e][a];
	  TString drawString = rhoVal[e] + ">>" + name;
	  TString drawCuts = ptSelection[p] + " && " + BBCselection[a] + "&&" + etaSelection[je];
	  cx->cd();
	  jetTree->Draw( drawString, drawCuts, "" );
	  hRho[p][je][e][a] = (TH1D*)gDirectory->Get( name );
	  hRho[p][je][e][a]->SetStats(0);

	  rho = hRho[p][je][e][a]->GetMean();
	  stdev = hRho[p][je][e][a]->GetMeanError();
      
	  hRhoByEta_jetEta__pt[p][a][je]->SetBinContent( e+1, rho );
	  hRhoByEta_jetEta__pt[p][a][je]->SetBinError( e+1, stdev );
	  hRhoByEta_jetEta__pt[p][a][je]->SetMarkerSize( 1.2 );
	}
	cpt2->cd();
	cpt2->cd(p+1);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridy();
	hRhoByEta_jetEta__pt[p][a][je]->Draw("PSAME");
      }
    }
  }

  cpt2->SaveAs("plots/UE/ptJetEta3plot.pdf","PDF");


  



  TCanvas *ceta2 = new TCanvas( "ceta2", "", 0, 23, 900, 300 );
  ceta2->Divide(3,1,0,0);
  for ( int e=0; e<nEtaBins; ++e ) { ceta2->cd(e+1);  sRhoByEta->Draw(); }
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {
      for ( int je=0; je<nEtaBins; ++je ) {
	ceta2->cd();
	ceta2->cd(je+1);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridy();
	hRhoByEta_jetEta__pt[p][a][je]->Draw("PSAME");
      }
    }
  }


  ceta2->SaveAs("plots/UE/etaJetEta3plot.pdf","PDF");
  

}
