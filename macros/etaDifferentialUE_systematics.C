// Veronica Verkest
// May 23, 2020


//#include <macros/plot.h>
// Standard library
#include <math.h>
#include <iostream>
#include <fstream>

// ROOT Library
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TCut.h>
#include <TPad.h>
#include <TLatex.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>


void etaDifferentialUE_systematics(){

  //gStyle->SetErrorX(0.0001);
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  //gSystem->Load("macros/plot.h");
  
  const double pi = 3.14159265;
  const double AREA = 4.*(pi/3.);

  const double eastArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )
  const double midArea = 2.*(0.6)*(pi/3.);   // (  0.6 in eta  ) X (  2pi/3 in phi  )
  const double westArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };
  const TString ptCorrectedBinString[nPtBins] = { "10 < p_{T} < 15", "15 < p_{T} < 20",  "20 < p_{T} < 30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };
  const int ptMarker[nPtBins] = { 20, 21, 29 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6 < #eta_{jet}^{lead} < -0.3", "-0.3 < #eta_{jet}^{lead} < 0.3", "0.3 < #eta_{jet}^{lead} < 0.6" };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgUE", "_neuUE", "_allUE" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

  const int nEAbins = 2;
  TString EAbinName[nEAbins] = { "Lo", "Hi" };
  TString EAbinString[nEAbins] = { "Low EA", "High EA" };
  TString BBCselection[nEAbins] = { "BbcAdcSumEast>3559.12 && BbcAdcSumEast<11503", "BbcAdcSumEast>26718.1" };
 
  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho_te+neuEastRho)", "(chgMidRho_te+neuMidRho)", "(chgWestRho_te+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  double area[nEtaBins] = { eastArea, midArea, westArea };
  
  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 24, 20 };
  
  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  string fileSuffix = "10cmVzCut_6cmVzDiff";
  
  TString fileName[nEAbins] = { ("out/UE/pAuHTjetUE_"+fileSuffix+"_loEA_diffEta.root").c_str(),
				("out/UE/pAuHTjetUE_"+fileSuffix+"_hiEA_diffEta.root").c_str() };
  TString efficFileName[nEAbins] = { "src/trackeffic_loEA.root", "src/trackeffic_hiEA.root" };

  TFile* inFile[nEAbins];
  TFile* efficFile[nEAbins];
    
  for (int i=0; i<nEAbins; ++i) {
    efficFile[i] = new TFile( efficFileName[i], "READ" );
  }

  const int eta_bins = 10;
  TH1D *hEff[eta_bins][nEAbins];
  
  for (int i=0; i<eta_bins; ++i) {
    for (int j=0; j<nEAbins; ++j){    
      int binno = i+1;
    
      if (j==0 ) { name = "eff_s_bin_1_3_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
      else if (j==1) {name = "eff_s_bin_7_10_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
    
      hEff[i][j] = (TH1D*)efficFile[j]->Get(name);

    }
  }

  const int ybins = 10;
  double ybinEdge[ybins+1] = { -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1 };
  
  //  Tree variables
  int RunID[nEAbins], EventID[nEAbins], nTowers[nEAbins], nPrimary[nEAbins], nGlobal[nEAbins], nVertices[nEAbins], refMult[nEAbins],
    gRefMult[nEAbins], nUEpart_chg[nEAbins], nUEpart_neu[nEAbins];
  double Vz[nEAbins], BbcAdcSumEast[nEAbins], leadPt[nEAbins], leadEta[nEAbins], leadPhi[nEAbins], chgEastRho[nEAbins], chgMidRho[nEAbins],
    chgWestRho[nEAbins], neuEastRho[nEAbins], neuMidRho[nEAbins], neuWestRho[nEAbins], leadArea[nEAbins], eastRho[nEAbins], midRho[nEAbins],
    westRho[nEAbins], leadPtCorrected[nEAbins], chgEastRho_te[nEAbins], chgMidRho_te[nEAbins], chgWestRho_te[nEAbins], rho_te[nEAbins], rho[nEAbins];

  TH2D *hChgUEpt[nEAbins][nEtaBins]; TH2D *hNeuUEpt[nEAbins][nEtaBins];
  TH1D *hRho[nEAbins][nEtaBins];  TH1D *hNchg[nEAbins][nEtaBins];  TH1D *hNchg_te[nEAbins][nEtaBins];  TH1D *hNneu[nEAbins][nEtaBins];
  TTree *jetTree[nEAbins]; int nEntries[nEAbins];

  for (int i=0; i<nEAbins;++i){
    inFile[i] = new TFile( fileName[i], "READ" );

    jetTree[i] = (TTree*) inFile[i]->Get("HTjetTree");

    jetTree[i]->SetBranchAddress( "RunID", &RunID[i] );
    jetTree[i]->SetBranchAddress( "EventID", &EventID[i] );
    jetTree[i]->SetBranchAddress( "nTowers", &nTowers[i] );
    jetTree[i]->SetBranchAddress( "nPrimary", &nPrimary[i] );
    jetTree[i]->SetBranchAddress( "nGlobal", &nGlobal[i] );
    jetTree[i]->SetBranchAddress( "nVertices", &nVertices[i] );
    jetTree[i]->SetBranchAddress( "refMult", &refMult[i] );
    jetTree[i]->SetBranchAddress( "gRefMult", &gRefMult[i] );
    jetTree[i]->SetBranchAddress( "Vz", &Vz[i] );
    jetTree[i]->SetBranchAddress( "leadPt", &leadPt[i] );
    jetTree[i]->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast[i] );
    jetTree[i]->SetBranchAddress( "leadEta", &leadEta[i] );
    jetTree[i]->SetBranchAddress( "leadPhi", &leadPhi[i] );
    jetTree[i]->SetBranchAddress( "chgEastRho", &chgEastRho[i] );
    jetTree[i]->SetBranchAddress( "chgMidRho", &chgMidRho[i] );
    jetTree[i]->SetBranchAddress( "chgWestRho", &chgWestRho[i] );
    jetTree[i]->SetBranchAddress( "chgEastRho_te", &chgEastRho_te[i] );
    jetTree[i]->SetBranchAddress( "chgMidRho_te", &chgMidRho_te[i] );
    jetTree[i]->SetBranchAddress( "chgWestRho_te", &chgWestRho_te[i] );
    jetTree[i]->SetBranchAddress( "neuEastRho", &neuEastRho[i] );
    jetTree[i]->SetBranchAddress( "neuMidRho", &neuMidRho[i] );
    jetTree[i]->SetBranchAddress( "neuWestRho", &neuWestRho[i] );
    jetTree[i]->SetBranchAddress( "leadArea", &leadArea[i] );
    jetTree[i]->SetBranchAddress( "leadPtCorrected", &leadPtCorrected[i] );
    jetTree[i]->SetBranchAddress( "nUEpart_chg", &nUEpart_chg[i] );
    jetTree[i]->SetBranchAddress( "nUEpart_neu", &nUEpart_neu[i] );
    jetTree[i]->SetBranchAddress( "rho", &rho[i] );
    jetTree[i]->SetBranchAddress( "rho_te", &rho_te[i] );

    nEntries[i] = jetTree[i]->GetEntries();

    for (int e=0; e<nEtaBins; ++e) {
      name = "hChgUePtEta"; name += etaBinName[e];
      hChgUEpt[i][e] = (TH2D*)inFile[i]->Get( name );
      name = "hNeuUePtEta"; name += etaBinName[e];
      hNeuUEpt[i][e] = (TH2D*)inFile[i]->Get( name );
    
      name = "hRho_"; name += EAbinName[i]; name += etaBinName[e];
      hRho[i][e] = new TH1D(name,"UE;<#rho> [GeV/#it{c}]", 30,0.0,15.0);
      
      name = "hNchg_"; name += EAbinName[i]; name += etaBinName[e];
      hNchg[i][e] = new TH1D( name,"# charged particles in UE", 50,0.0,50.0);

      name = "hNchg_te_"; name += EAbinName[i]; name += etaBinName[e];
      hNchg_te[i][e] = new TH1D( name,"# track effic. corr. charged particles in UE", 50,0.0,50.0);
      
      name = "hNneu_"; name += EAbinName[i]; name += etaBinName[e];
      hNneu[i][e] = new TH1D( name,"# neutral particles in UE", 50,0.0,50.0);
    }


  }


  for (int a=0; a<nEAbins; ++a) {
    for (int i=0; i<nEntries[a]; ++i) {
      jetTree[a]->GetEntry(i);
      if ( !(leadPtCorrected[a]>ptLo[0] && leadPtCorrected[a]<ptHi[2]) ) { continue; }
      // if ( !(leadPtCorrected[a]>ptLo[1] && leadPtCorrected[a]<ptHi[1]) ) { continue; }
      else {

	jeval = 99;

	for ( int e=0; e<nEtaBins; ++e ) {
	  if ( leadEta[a] >= etaLo[e]  &&  leadEta[a] <= etaHi[e] ) { jeval = e; }
	}
	if ( jeval==99 /*|| jeval==99*/ ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<leadEta[a]<<endl<<endl; }
	else {
	  hRho[a][jeval]->Fill(rho[a]);
	  hNchg[a][jeval]->Fill(nUEpart_chg[a]);
	  hNneu[a][jeval]->Fill(nUEpart_neu[a]);
	}
      }
      
    }
  }

  
  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {

      hChgUEpt[a][e]->Scale(1./hRho[a][e]->GetEntries());

      new TCanvas;
      hChgUEpt[a][e]->Draw();
      
      hNeuUEpt[a][e]->Scale(1./hRho[a][e]->GetEntries());

    }
  }
  
  TH1D *hChgPtDist[ybins][nEAbins][nEtaBins];
  
  int ecolor[ybins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  TString ybinEdgeString[ybins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };
  
  TCanvas * ceta[nEtaBins];

  TH1D *hTreeChgRho[nEAbins][nEtaBins];
  TH1D *hChgDist[nEAbins][nEtaBins];

  TLegend *legeta[nEtaBins];
  
  for (int e=0; e<nEtaBins; ++e) {
    name = "ceta_"; name += e;
    ceta[e] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 0
    ceta[e]->SetLogy();
    for (int a=0; a<nEAbins; ++a) {
      for (int i=0; i<ybins; ++i) {
	int binno = i+1;
	name = "hChgpTdist"; name += etaBinName[e]; name += EAbinName[a]; name += "EA_etabin"; name += binno;
	hChgUEpt[a][e]->GetYaxis()->SetRange( binno, binno+1 );
	hChgUEpt[a][e]->SetLineColor( ecolor[i] );
	hChgUEpt[a][e]->SetMarkerColor( ecolor[i] );
	hChgUEpt[a][e]->SetMarkerStyle( etaMarker[e] );	
	hChgPtDist[i][a][e] = (TH1D*) hChgUEpt[a][e]->ProjectionX( name, binno, binno+1 );
	title = "Charged UE p_{T} dist. for "; title += etaBinString[e]; title += " (";
	title += ybinEdgeString[i]; title +="<#eta<"; title += ybinEdgeString[binno]; title += ");p_{T} [GeV/#it{c}]";
	hChgPtDist[i][a][e]->SetTitle( title );
	hChgPtDist[i][a][e]->GetYaxis()->SetRangeUser(0.000001,0.3);
	if (i==0) { hChgPtDist[i][a][e]->Draw(""); }
	else { hChgPtDist[i][a][e]->Draw("SAME"); }
      }
    }

    legeta[e] = (TLegend*) ceta[e]->BuildLegend(0.4,0.45,0.75,0.8);
    legeta[e]->SetFillColorAlpha(0,0);  legeta[e]->SetLineWidth(0);
    //ceta[e]->Close();
  }


  //  TRACKING EFFICIENCY CORRECTIONS!
  TH2D *hChgUEpt_te[nEAbins][nEtaBins];
  TH2D *hChgUEpt_te_sys1[nEAbins][nEtaBins];
  TH2D *hChgUEpt_te_sys2[nEAbins][nEtaBins];

  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      name = "hChgUEpt"; name += etaBinName[e]; name += EAbinName[a];  // rename these!
      hChgUEpt[a][e]->SetName( name );

      name = "hNeuUEpt"; name += etaBinName[e]; name += EAbinName[a];
      hNeuUEpt[a][e]->SetName( name );
      
      hChgUEpt_te[a][e] = (TH2D*) hChgUEpt[a][e]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += etaBinName[e]; name += EAbinName[a]; name += "_te";
      hChgUEpt_te[a][e]->SetName( name );

      hChgUEpt_te_sys1[a][e] = (TH2D*) hChgUEpt[a][e]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += etaBinName[e]; name += EAbinName[a]; name += "_te_sys1";
      hChgUEpt_te_sys1[a][e]->SetName( name );

      hChgUEpt_te_sys2[a][e] = (TH2D*) hChgUEpt[a][e]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += etaBinName[e]; name += EAbinName[a]; name += "_te_sys2";
      hChgUEpt_te_sys2[a][e]->SetName( name );


      for (int iy=0; iy<hChgUEpt_te[a][e]->GetNbinsY(); ++iy) {  // loop over 20 eta bins

	int fileNo = 99;
	int ybin = iy+1;
	
	double binCenter = hChgUEpt_te[a][e]->GetYaxis()->GetBinCenter( ybin );

	for (int j=0; j<ybins; ++j) {
	  if ( binCenter>ybinEdge[j] && binCenter<ybinEdge[j+1] ) { fileNo = j; };
	}
	if ( fileNo==99 ) { cerr<<"CANNOT FIND EFFICIENCY HISTOGRAM"<<endl; }


	for(int ix = 1; ix <= hChgUEpt[a][e]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][e]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin);
	  double old_value = hChgUEpt[a][e]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][e]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te[a][e]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te[a][e]->SetBinError( xbin, ybin, corr_err );
	}


	// sys1: trackeffic+0.05
	for(int ix = 1; ix <= hChgUEpt[a][e]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][e]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin) + 0.05;
	  double old_value = hChgUEpt[a][e]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][e]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te_sys1[a][e]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te_sys1[a][e]->SetBinError( xbin, ybin, corr_err );
	}


	// sys2: trackeffic-0.05
	for(int ix = 1; ix <= hChgUEpt[a][e]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][e]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin) - 0.05;
	  double old_value = hChgUEpt[a][e]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][e]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te_sys2[a][e]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te_sys2[a][e]->SetBinError( xbin, ybin, corr_err );
	}

	
      }
	
    }
  }



  double nChg[nEtaBins][nEAbins][nEtaBins];
  double nChg_err[nEtaBins][nEAbins][nEtaBins];
  double meanChgPt[nEtaBins][nEAbins][nEtaBins];
  double te_meanChgPt[nEtaBins][nEAbins][nEtaBins];
  
  double nChg_te[nEtaBins][nEAbins][nEtaBins];
  double nChg_te_err[nEtaBins][nEAbins][nEtaBins];
  double meanChgPt_err[nEtaBins][nEAbins][nEtaBins];
  double te_meanChgPt_err[nEtaBins][nEAbins][nEtaBins];

  double nChg_te_sys1[nEtaBins][nEAbins][nEtaBins];
  double te_meanChgPt_sys1[nEtaBins][nEAbins][nEtaBins];
  double nChg_te_sys2[nEtaBins][nEAbins][nEtaBins];
  double te_meanChgPt_sys2[nEtaBins][nEAbins][nEtaBins];

  double chgRho[nEtaBins][nEAbins][nEtaBins];
  
  for ( int a=0; a<nEAbins; ++a) {
    for ( int je=0; je<nEtaBins; ++je) {

      for ( int e=0; e<nEtaBins; ++e) {

	hChgUEpt[a][je]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );	
	chgRho[e][a][je] = ( hChgUEpt[a][je]->GetMean(1) * hChgUEpt[a][je]->Integral() ) / area[e];


	//  WITH TRACKING EFFICIENCY
	hChgUEpt_te[a][je]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt[e][a][je] = hChgUEpt_te[a][je]->GetMean(1);
	te_meanChgPt_err[e][a][je] = hChgUEpt_te[a][je]->GetMeanError(1);

	nChg_te[e][a][je] = hChgUEpt_te[a][je]->Integral();
	nChg_te_err[e][a][je] = hChgUEpt_te[a][je]->GetMeanError();

	// TRACKING EFFICIENCY SYSTEMATICS
	hChgUEpt_te_sys1[a][je]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt_sys1[e][a][je] = hChgUEpt_te_sys1[a][je]->GetMean(1);
	nChg_te_sys1[e][a][je] = hChgUEpt_te_sys1[a][je]->Integral();

	hChgUEpt_te_sys2[a][je]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt_sys2[e][a][je] = hChgUEpt_te_sys2[a][je]->GetMean(1);
	nChg_te_sys2[e][a][je] = hChgUEpt_te_sys1[a][je]->Integral();

      }      
    }
  }

  double nChg_te_sys[nEtaBins][nEAbins][nEtaBins];
  double te_meanChgPt_sys[nEtaBins][nEAbins][nEtaBins];

  for ( int a=0; a<nEAbins; ++a) {
    for ( int je=0; je<nEtaBins; ++je) {
      for ( int e=0; e<nEtaBins; ++e) {

	//nChg_te
	double sys1 = fabs( 1 - ( nChg_te_sys1[e][a][je]/nChg_te[e][a][je] ) );
	double sys2 = fabs( 1 - ( nChg_te_sys2[e][a][je]/nChg_te[e][a][je] ) );
	nChg_te_sys[e][a][je] = max( sys1, sys2 );

	sys1 = fabs( 1 - ( te_meanChgPt_sys1[e][a][je] / te_meanChgPt[e][a][je] ) );
	sys2 = fabs( 1 - ( te_meanChgPt_sys2[e][a][je] / te_meanChgPt[e][a][je] ) );
	te_meanChgPt_sys[e][a][je] = max( sys1, sys2 );
	
      }
    }
  }

  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


  double etaBinEdgeCenter[nEtaBins+1] = { -0.5, 0.5, 1.5, 2.5 };
  double etaBinEdge[nEtaBins+1][nEtaBins]; //  {  {...}, {...}, {...}, {...}  }

  double shift[nEtaBins] = { -0.25, 0.0, 0.25 };

  // for (int i=0; i<nEtaBins+1; ++i){
  //   for (int j=0; j<nEtaBins; ++j){
  //     etaBinEdge[i][j] = 
  //   }
  // }

  double yEdge[2] = {0.5,1.8};


  
  //East <#rho>          Mid <#rho>          West <#rho>
  TH2D *hscale0 = new TH2D("hscale0",";;#LT #frac{d#it{N}_{ch}}{d#eta d#phi} #GT", 3,0.5,2.5,11,0.5,1.8);
  hscale0->GetYaxis()->SetTitleOffset(1);
  hscale0->GetXaxis()->SetTitleOffset(1.0);
  hscale0->GetYaxis()->SetTitleSize(0.05);
  hscale0->GetXaxis()->SetLabelSize(0.07);
  hscale0->GetYaxis()->SetLabelSize(0.04);
  hscale0->GetXaxis()->SetTickSize(0.0);
  hscale0->GetXaxis()->CenterTitle();
  hscale0->SetName("");
  //hscale0->GetYaxis()->SetTitleOffset(1.25);
  hscale0->GetXaxis()->SetBinLabel(1,UEetaBinString[0]);
  hscale0->GetXaxis()->SetBinLabel(2,UEetaBinString[1]);
  hscale0->GetXaxis()->SetBinLabel(3,UEetaBinString[2]);
  //hscale0->GetXaxis()->SetTitleOffset(1.25);
  TCanvas *c0 = new TCanvas;  
  c0->SetMargin(0.125,0.05,0.125,0.05);
  
  TH1D *hCHGdNdetadphi_te[nEAbins][nEtaBins];
  TH1D *hCHGdNdetadphi_te_sys[nEAbins][nEtaBins];
  hscale0->Draw();
  hscale0->SetName(" ");
  hscale0->SetLineWidth(0);
  hscale0->SetStats(0);

  TLegend *leg0;
  
  leg0 = new TLegend(0.45, 0.67, 0.95, 0.95,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(0);  leg0->SetLineColor(1);  leg0->SetLineStyle(1);
  leg0->SetLineWidth(1);  leg0->SetFillColorAlpha(0,0.0);  leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);

  leg0->AddEntry((TObject*)0,"0-30% EA", "");
  leg0->AddEntry((TObject*)0,"70-90% EA", "");

  //errorFill[nEAbins] = {  };
  int binno = 0;

  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hCHGdNdetadphi_te_" + EAbinName[a] + etaBinName[e];
      title = etaBinString[e] + "   ";
      hCHGdNdetadphi_te[a][e] = new TH1D( name, title , 11,0.5,2.5 );

      hCHGdNdetadphi_te[a][e]->SetStats(0);
      hCHGdNdetadphi_te[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi_te[a][e]->SetMarkerSize( 2 );

      leg0->AddEntry(hCHGdNdetadphi_te[a][e],title,"lpf");
      
      name = "hCHGdNdetadphi_te_sys_" + EAbinName[a] + etaBinName[e];
      hCHGdNdetadphi_te_sys[a][e] = new TH1D( name, "" , 11,0.5,2.5 );
      
      hCHGdNdetadphi_te_sys[a][e]->SetStats(0);
      hCHGdNdetadphi_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerSize( 0 );

        if (a==0) {
	  hCHGdNdetadphi_te_sys[a][e]->SetFillStyle(3002);
	  hCHGdNdetadphi_te_sys[a][e]->SetFillColor( etaColor[e] );
	}
	else { hCHGdNdetadphi_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 ); }

      
      //hCHGdNdetadphi_te[a][e];

      
      for ( int je=0; je<nEtaBins; ++je) {
	
	binno = (3*je) + e + je + 1; 
	
	double value = (double) ( nChg_te[je][a][e] / area[je] );
	hCHGdNdetadphi_te[a][e]->SetBinContent( binno, value );
	hCHGdNdetadphi_te[a][e]->SetBinError( binno, nChg_te_err[je][a][e] );

	hCHGdNdetadphi_te_sys[a][e]->SetBinContent( binno, value );
	hCHGdNdetadphi_te_sys[a][e]->SetBinError( binno, value*nChg_te_sys[je][a][e] );

      }
    }
  }

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hCHGdNdetadphi_te_sys[a][e]->Draw("e2SAME");
    }
  }

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hCHGdNdetadphi_te[a][e]->Draw("lpfX0ESAME");
    }
  }

  
  TPaveText *sp1 = new TPaveText(0.175,0.75,0.5,1.0,"NDC");
  //sp1->AddText("STAR Preliminary");
  sp1->AddText("");
  sp1->SetTextSize(0.04);
  sp1->SetTextColor(kRed);
  sp1->SetFillColorAlpha(0,0.0);
  sp1->SetLineColorAlpha(0,0.0);
  sp1->Draw("NB");

  
  TLatex *tex0 = new TLatex(0.175,0.8,"p+Au #sqrt{#it{s}_{NN}} = 200 GeV");
  tex0->SetTextFont(63);
  tex0->SetTextSize(16);
  tex0->SetTextColor(kBlack);
  tex0->SetLineWidth(1);
  tex0->SetNDC();
  tex0->Draw();
  TLatex *tex1 = new TLatex(0.175,0.75,"E^{trig}_{T} > 5.4 GeV");
  tex1->SetTextFont(63);
  tex1->SetTextSize(16);
  tex1->SetTextColor(kBlack);
  tex1->SetLineWidth(1);
  tex1->SetNDC();
  tex1->Draw();
  TLatex *tex2 = new TLatex(0.175,0.7,"anti-k_{T} R=0.4 jets");
  tex2->SetTextFont(63);
  tex2->SetTextSize(16);
  tex2->SetTextColor(kBlack);
  tex2->SetLineWidth(1);
  tex2->SetNDC();
  tex2->Draw();
  TLatex *tex3 = new TLatex(0.175,0.65,"#left|#eta^{jets}#right| < 0.6");
  tex3->SetTextFont(63);
  tex3->SetTextSize(16);
  tex3->SetTextColor(kBlack);
  tex3->SetLineWidth(1);
  tex3->SetNDC();
  tex3->Draw();
  TLatex *tex4 = new TLatex(0.175,0.6,"10 < p_{T,lead}^{reco} < 30 [GeV/#it{c}]");
  tex4->SetTextFont(63);
  tex4->SetTextSize(16);
  tex4->SetTextColor(kBlack);
  tex4->SetLineWidth(1);
  tex4->SetNDC();
  tex4->Draw();
  
  
  leg0->Draw();
  c0->SaveAs(("plots/UE/"+fileSuffix+"/CHGdNdetadphi_te_eta_systematics.pdf").c_str(),"PDF");
  //c0->Close();
  

  TCanvas *c1 = new TCanvas;
  c1->SetMargin(0.125,0.05,0.125,0.05);

  double yEdge1[2] = { 0.55, 0.85 };
  TH2D *hscale3 = new TH2D("hscale3",";;#LT p_{T}^{ch} #GT [GeV/#it{c}]", 3,0.5,2.5,10,0.55,0.85);
  hscale3->GetYaxis()->SetTitleOffset(1.25);
  hscale3->GetYaxis()->SetTitleOffset(1.25);
  hscale3->GetYaxis()->SetTitleOffset(1.0);
  hscale3->GetYaxis()->SetTitleSize(0.05);
  hscale3->GetXaxis()->SetLabelSize(0.07);
  hscale3->GetYaxis()->SetLabelSize(0.04);
  hscale3->GetXaxis()->SetTickSize(0.0);
  hscale3->GetXaxis()->CenterTitle();
  hscale3->SetName("");
  //hscale3->GetYaxis()->SetTitleOffset(1.25);
  hscale3->GetXaxis()->SetBinLabel(1,UEetaBinString[0]);
  hscale3->GetXaxis()->SetBinLabel(2,UEetaBinString[1]);
  hscale3->GetXaxis()->SetBinLabel(3,UEetaBinString[2]);

  
  TH1D *hChg_pt_te[nEAbins][nEtaBins];
  TH1D *hChg_pt_te_sys[nEAbins][nEtaBins];
  hscale3->SetName(" ");
  hscale3->SetLineWidth(0);
  hscale3->SetStats(0);
  hscale3->Draw();

  // TLegend *leg1;
  
  // leg1 = new TLegend(0.5, 0.65, 0.92, 0.92,NULL,"brNDC");    // LEGEND 0
  // leg1->SetBorderSize(0);  leg1->SetLineColor(1);  leg1->SetLineStyle(1);
  // leg1->SetLineWidth(1);  leg1->SetFillColor(0);  leg1->SetFillStyle(1001);
  // leg1->SetNColumns(2);

  for ( int e=0; e<nEtaBins; ++e) {
    for ( int a=0; a<nEAbins; ++a) {

      name = "hChg_MeanPt_te_" + EAbinName[a] + etaBinName[e];
      title = EAbinString[a] + "   " + etaBinString[e] + "   ";
      hChg_pt_te[a][e] = new TH1D( name, title , 11,0.5,2.5 );

      hChg_pt_te[a][e]->SetStats(0);
      hChg_pt_te[a][e]->SetLineColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt_te[a][e]->SetMarkerSize( 2 );

      //leg1->AddEntry( hChg_pt_te[a][e], title, "p" );
      
      name = "hChg_MeanPt_te_sys_" + EAbinName[a] + etaBinName[e];
      hChg_pt_te_sys[a][e] = new TH1D( name, "" , 11,0.5,2.5 );

      hChg_pt_te_sys[a][e]->SetStats(0);
      hChg_pt_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      hChg_pt_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      hChg_pt_te_sys[a][e]->SetMarkerSize( 0.0 );

      if (a==0) {
	hChg_pt_te_sys[a][e]->SetFillStyle(3002);
	hChg_pt_te_sys[a][e]->SetFillColor( etaColor[e]);
      }
      else { hChg_pt_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 ); }

      //hCHGdNdetadphi_te[a][e];
	
      for ( int je=0; je<nEtaBins; ++je) {
	
	binno = (3*je) + e + je + 1; 
	
	double value = te_meanChgPt[je][a][e];
	hChg_pt_te[a][e]->SetBinContent( binno, value );
	hChg_pt_te[a][e]->SetBinError( binno, te_meanChgPt_err[je][a][e] );

	hChg_pt_te_sys[a][e]->SetBinContent( binno, value );
	hChg_pt_te_sys[a][e]->SetBinError( binno, value*te_meanChgPt_sys[je][a][e] );
	
      }
    }
  }

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hChg_pt_te_sys[a][e]->Draw("E2SAME");
    }
  }
  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hChg_pt_te[a][e]->Draw("lpfX0ESAME");
    }
  }
  
  leg0->Draw();


  sp1->Draw("NB");
  tex0->Draw();
  tex1->Draw();
  tex2->Draw();
  tex3->Draw();
  tex4->Draw();

 
  c1->SaveAs(("plots/UE/"+fileSuffix+"/Chg_MeanPt_te_eta_systematics.pdf").c_str(),"PDF");



}

