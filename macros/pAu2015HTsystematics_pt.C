// Veronica Verkest
// May 25, 2020

void pAu2015HTsystematics_pt(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  const double AREA = 4.*(pi/3.);

  const double eastArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )
  const double midArea = 2.*(0.6)*(pi/3.);   // (  0.6 in eta  ) X (  2pi/3 in phi  )
  const double westArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  
)
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const TString ptCorrectedBinString[nPtBins] = { "10<p_{T}<15", "15<p_{T}<20",  "20<p_{T}<30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };
  const int ptMarker[nPtBins] = { 20, 21, 29 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
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
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  TString fileName[nEAbins] = { "out/UE/pAuHTjetUE_loEA_diffPt.root", "out/UE/pAuHTjetUE_hiEA_diffPt.root" };
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
    gRefMult[nEAbins], nBGpart_chg[nEAbins], nBGpart_neu[nEAbins];
  double Vz[nEAbins], BbcAdcSumEast[nEAbins], leadPt[nEAbins], leadEta[nEAbins], leadPhi[nEAbins], chgEastRho[nEAbins], chgMidRho[nEAbins],
    chgWestRho[nEAbins], neuEastRho[nEAbins], neuMidRho[nEAbins], neuWestRho[nEAbins], leadArea[nEAbins], eastRho[nEAbins], midRho[nEAbins],
    westRho[nEAbins], leadPtCorrected[nEAbins], chgEastRho_te[nEAbins], chgMidRho_te[nEAbins], chgWestRho_te[nEAbins], rho_te[nEAbins], rho[nEAbins];

  TH2D *hChgUEpt[nEAbins][nPtBins]; TH2D *hNeuUEpt[nEAbins][nPtBins];
  TH1D *hRho[nEAbins][nPtBins];  TH1D *hNchg[nEAbins][nPtBins];  TH1D *hNchg_te[nEAbins][nPtBins];  TH1D *hNneu[nEAbins][nPtBins];
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
    jetTree[i]->SetBranchAddress( "nBGpart_chg", &nBGpart_chg[i] );
    jetTree[i]->SetBranchAddress( "nBGpart_neu", &nBGpart_neu[i] );
    jetTree[i]->SetBranchAddress( "rho", &rho[i] );
    jetTree[i]->SetBranchAddress( "rho_te", &rho_te[i] );

    nEntries[i] = jetTree[i]->GetEntries();

    for (int p=0; p<nPtBins; ++p) {
      name = "hChgBgPtEta"; name += ptBinName[p];
      hChgUEpt[i][p] = (TH2D*)inFile[i]->Get( name );
      name = "hNeuBgPtEta"; name += ptBinName[p];
      hNeuUEpt[i][p] = (TH2D*)inFile[i]->Get( name );
    
      name = "hRho_"; name += EAbinName[i]; name += ptBinName[p];
      hRho[i][p] = new TH1D(name,"UE;<#rho> [GeV]", 30,0.0,15.0);
      
      name = "hNchg_"; name += EAbinName[i]; name += ptBinName[p];
      hNchg[i][p] = new TH1D( name,"# charged particles in UE", 50,0.0,50.0);

      name = "hNchg_te_"; name += EAbinName[i]; name += ptBinName[p];
      hNchg_te[i][p] = new TH1D( name,"# track effic. corr. charged particles in UE", 50,0.0,50.0);
      
      name = "hNneu_"; name += EAbinName[i]; name += ptBinName[p];
      hNneu[i][p] = new TH1D( name,"# neutral particles in UE", 50,0.0,50.0);
    }


  }


  for (int a=nEAbins-1; a>=0; --a) {
    for (int i=0; i<nEntries[a]; ++i) {
      jetTree[a]->GetEntry(i);
      if ( !(leadPtCorrected[a]>ptLo[0] && leadPtCorrected[a]<ptHi[2]) ) { continue; }
      // if ( !(leadPtCorrected[a]>ptLo[1] && leadPtCorrected[a]<ptHi[1]) ) { continue; }
      else {

	pval = 99;

	for ( int p=0; p<3; ++p ) {
	  if ( leadPtCorrected[a] >= ptLo[p]  &&  leadPtCorrected[a] <= ptHi[p] ) { pval = p; }
	}
	if ( pval==99 /*|| jeval==99*/ ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<leadPtCorrected<<endl<<endl; }

	hRho[a][pval]->Fill(rho[a]);
	hNchg[a][pval]->Fill(nBGpart_chg[a]);
	hNneu[a][pval]->Fill(nBGpart_neu[a]);
      }
      
    }
  }

  
  for (int a=nEAbins-1; a>=0; --a) {
    for (int p=0; p<nPtBins; ++p) {

      hChgUEpt[a][p]->Scale(1./hRho[a][p]->GetEntries());      
      hNeuUEpt[a][p]->Scale(1./hRho[a][p]->GetEntries());

    }
  }
  
  TH1D *hChgPtDist[ybins][nEAbins][nPtBins];
  
  int ecolor[ybins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  TString ybinEdgeString[ybins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };
  
  TCanvas * cpt[nPtBins];

  TH1D *hTreeChgRho[nEAbins][nPtBins];
  TH1D *hChgDist[nEAbins][nPtBins];

  TLegend *legpt[nPtBins];
  
  for (int p=0; p<nPtBins; ++p) {
    name = "cpt_"; name += p;
    cpt[p] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 0
    cpt[p]->SetLogy();
    for (int a=nEAbins-1; a>=0; --a) {
      for (int i=0; i<ybins; ++i) {
	int binno = i+1;
	name = "hChgpTdist"; name += ptBinName[p]; name += EAbinName[a]; name += "EA_etabin"; name += binno;
	hChgUEpt[a][p]->GetYaxis()->SetRange( binno, binno+1 );
	hChgUEpt[a][p]->SetLineColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerStyle( ptMarker[p] );	
	hChgPtDist[i][a][p] = (TH1D*) hChgUEpt[a][p]->ProjectionX( name, binno, binno+1 );
	title = "Charged UE p_{T} dist. for "; title += ptBinString[p]; title += " (";
	title += ybinEdgeString[i]; title +="<#eta<"; title += ybinEdgeString[binno]; title += ");p_{T} [GeV]";
	hChgPtDist[i][a][p]->SetTitle( title );
	hChgPtDist[i][a][p]->GetYaxis()->SetRangeUser(0.000001,0.3);
	if (i==0) { hChgPtDist[i][a][p]->Draw(""); }
	else { hChgPtDist[i][a][p]->Draw("SAME"); }
      }
    }

    legpt[p] = (TLegend*) cpt[p]->BuildLegend(0.4,0.45,0.75,0.8);
    legpt[p]->SetFillColorAlpha(0,0);  legpt[p]->SetLineWidth(0);
    cpt[p]->Close();
  }


  //  TRACKING EFFICIENCY CORRECTIONS!
  TH2D *hChgUEpt_te[nEAbins][nPtBins];
  TH2D *hChgUEpt_te_sys1[nEAbins][nPtBins];
  TH2D *hChgUEpt_te_sys2[nEAbins][nPtBins];

  for (int a=nEAbins-1; a>=0; --a) {
    for (int p=0; p<nPtBins; ++p) {
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinName[a];  // rename these!
      hChgUEpt[a][p]->SetName( name );

      name = "hNeuUEpt"; name += ptBinName[p]; name += EAbinName[a];
      hNeuUEpt[a][p]->SetName( name );
      
      hChgUEpt_te[a][p] = (TH2D*) hChgUEpt[a][p]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinName[a]; name += "_te";
      hChgUEpt_te[a][p]->SetName( name );

      hChgUEpt_te_sys1[a][p] = (TH2D*) hChgUEpt[a][p]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinName[a]; name += "_te_sys1";
      hChgUEpt_te_sys1[a][p]->SetName( name );

      hChgUEpt_te_sys2[a][p] = (TH2D*) hChgUEpt[a][p]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinName[a]; name += "_te_sys2";
      hChgUEpt_te_sys2[a][p]->SetName( name );


      for (int iy=0; iy<hChgUEpt_te[a][p]->GetNbinsY(); ++iy) {  // loop over 20 eta bins

	int fileNo = 99;
	int ybin = iy+1;
	
	double binCenter = hChgUEpt_te[a][p]->GetYaxis()->GetBinCenter( ybin );

	for (int j=0; j<ybins; ++j) {
	  if ( binCenter>ybinEdge[j] && binCenter<ybinEdge[j+1] ) { fileNo = j; };
	}
	if ( fileNo==99 ) { cerr<<"CANNOT FIND EFFICIENCY HISTOGRAM"<<endl; }


	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin);
	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te[a][p]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te[a][p]->SetBinError( xbin, ybin, corr_err );
	}


	// sys1: trackeffic+0.05
	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin) + 0.05;
	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te_sys1[a][p]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te_sys1[a][p]->SetBinError( xbin, ybin, corr_err );
	}


	// sys2: trackeffic-0.05
	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
	  
	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin) - 0.05;
	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te_sys2[a][p]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te_sys2[a][p]->SetBinError( xbin, ybin, corr_err );
	}

	
      }
	
    }
  }



  double nChg[nEtaBins][nEAbins][nPtBins];
  double nChg_err[nEtaBins][nEAbins][nPtBins];
  double meanChgPt[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt[nEtaBins][nEAbins][nPtBins];
  
  double nChg_te[nEtaBins][nEAbins][nPtBins];
  double nChg_te_err[nEtaBins][nEAbins][nPtBins];
  double meanChgPt_err[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt_err[nEtaBins][nEAbins][nPtBins];

  double nChg_te_sys1[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt_sys1[nEtaBins][nEAbins][nPtBins];
  double nChg_te_sys2[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt_sys2[nEtaBins][nEAbins][nPtBins];

  double chgRho[nEtaBins][nEAbins][nPtBins];
  
  for (int a=nEAbins-1; a>=0; --a) {
    for ( int p=0; p<nPtBins; ++p) {

      for ( int e=0; e<nEtaBins; ++e) {

	hChgUEpt[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );	
	chgRho[e][a][p] = ( hChgUEpt[a][p]->GetMean(1) * hChgUEpt[a][p]->Integral() ) / area[e];


	//  WITH TRACKING EFFICIENCY
	hChgUEpt_te[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt[e][a][p] = hChgUEpt_te[a][p]->GetMean(1);
	te_meanChgPt_err[e][a][p] = hChgUEpt_te[a][p]->GetMeanError(1);

	nChg_te[e][a][p] = hChgUEpt_te[a][p]->Integral();
	nChg_te_err[e][a][p] = hChgUEpt_te[a][p]->GetMeanError();

	// TRACKING EFFICIENCY SYSTEMATICS
	hChgUEpt_te_sys1[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt_sys1[e][a][p] = hChgUEpt_te_sys1[a][p]->GetMean(1);
	nChg_te_sys1[e][a][p] = hChgUEpt_te_sys1[a][p]->Integral();

	hChgUEpt_te_sys2[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt_sys2[e][a][p] = hChgUEpt_te_sys2[a][p]->GetMean(1);
	nChg_te_sys2[e][a][p] = hChgUEpt_te_sys1[a][p]->Integral();

      }      
    }
  }

  double nChg_te_sys[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt_sys[nEtaBins][nEAbins][nPtBins];

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int p=0; p<nPtBins; ++p) {
      for ( int e=0; e<nEtaBins; ++e) {

	//nChg_te
	double sys1 = fabs( 1 - ( nChg_te_sys1[e][a][p]/nChg_te[e][a][p] ) );
	double sys2 = fabs( 1 - ( nChg_te_sys2[e][a][p]/nChg_te[e][a][p] ) );
	nChg_te_sys[e][a][p] = max( sys1, sys2 );

	sys1 = fabs( 1 - ( te_meanChgPt_sys1[e][a][p] / te_meanChgPt[e][a][p] ) );
	sys2 = fabs( 1 - ( te_meanChgPt_sys2[e][a][p] / te_meanChgPt[e][a][p] ) );
	te_meanChgPt_sys[e][a][p] = max( sys1, sys2 );
	
      }
    }
  }

  double dNx[nEAbins][nPtBins], dNy[nEAbins][nPtBins], dNx_err[nEAbins][nPtBins], dNy_err[nEAbins][nPtBins];

  double stager[nEtaBins] = { -0.25, 0.0, 0.25 };
  
  for (int a=0, a<nEAbins; ++a) {
    for ( int e=0; e<nEtaBins; ++e) {
      for ( int p=0; p<nPtBins; ++p) {
	dNx[a][p] = e+1+stagger[e];
	dNy[a][p] = nChg_te[nEtaBins][nEAbins][nPtBins];
	dNx_err[a][p] = 0.3;
	dNy_err[a][p] = nChg_te_err[nEtaBins][nEAbins][nPtBins];
	ptx[a][p] = ;
	pty[a][p] = ;
	ptx_err[a][p] = ;
	pty_err[a][p] = ;
      }
    }
  }
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


  double binEdge[nEtaBins+1] = { -0.5, 0.5, 1.5, 2.5 };

  double ptBinEdge[nPtBins+1] = { 10.0, 15.0, 20.0, 30.0 };
  double yEdge[2] = {0.5,1.5};



  auto Ngraph[nEAbins][nPtBins];

  for (int a=0, a<nEtaBins; ++a) {
    for ( int p=0; p<nPtBins; ++p) {
      Ngraph[a][p] = new TGraphErrors(3, dNx[a][p], dNy[a][p], dNx_err[a][p], dNy_err[a][p]);
      new TCanvas;
      Ngraph[a][p]->Draw("5");
    }
  }
  

  TH2D *hscale0 = new TH2D("hscale0",";UE-subtracted lead jet p_{T} [GeV];<#frac{dN_{ch}}{d#eta d#phi}>", nPtBins,ptBinEdge,1,yEdge);
  hscale0->SetName("");
  hscale0->GetYaxis()->SetTitleOffset(1.25);
  hscale0->GetXaxis()->SetTitleOffset(1.25);
  TCanvas *c0 = new TCanvas;  
  c0->SetMargin(0.125,0.05,0.125,0.05);
  
  TH1D *hCHGdNdetadphi_te[nEAbins][nEtaBins];
  TH1D *hCHGdNdetadphi_te_sys[nEAbins][nEtaBins];
  hscale0->Draw();
  hscale0->SetName(" ");
  hscale0->SetLineWidth(0);
  hscale0->SetStats(0);

  TLegend *leg0;
  
  leg0 = new TLegend(0.5, 0.65, 0.92, 0.92,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(0);  leg0->SetLineColor(1);  leg0->SetLineStyle(1);
  leg0->SetLineWidth(1);  leg0->SetFillColor(0);  leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);

  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hCHGdNdetadphi_te_" + EAbinName[a] + etaBinName[e];
      title = EAbinString[a] + "   " + eastmidwest[e] + "   ";
      hCHGdNdetadphi_te[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hCHGdNdetadphi_te[a][e]->SetStats(0);
      hCHGdNdetadphi_te[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi_te[a][e]->SetMarkerSize( 2 );

      leg0->AddEntry(hCHGdNdetadphi_te[a][e],title,"lpf");
      
      name = "hCHGdNdetadphi_te_sys_" + EAbinName[a] + etaBinName[e];
      hCHGdNdetadphi_te_sys[a][e] = new TH1D( name, "" , nPtBins,ptBinEdge );
      
      hCHGdNdetadphi_te_sys[a][e]->SetStats(0);
      hCHGdNdetadphi_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 );
      hCHGdNdetadphi_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerSize( 0 );
      
      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	
	double value = (double) ( nChg_te[e][a][p] / area[e] );
	hCHGdNdetadphi_te[a][e]->SetBinContent( pbin, value );
	hCHGdNdetadphi_te[a][e]->SetBinError( pbin, nChg_te_err[e][a][p] );

	hCHGdNdetadphi_te_sys[a][e]->SetBinContent( pbin, value );
	hCHGdNdetadphi_te_sys[a][e]->SetBinError( pbin, value*nChg_te_sys[e][a][p] );

      }
    }
  }

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hCHGdNdetadphi_te[a][e]->Draw("lpfX0ESAME");
      hCHGdNdetadphi_te_sys[a][e]->Draw("e2SAME");
    }
  }
  
  TPaveText *sp1 = new TPaveText(0.15,0.7,0.5,1.0,"NDC");
  sp1->AddText("STAR Preliminary");
  sp1->SetTextSize(0.04);
  sp1->SetTextColor(kRed);
  sp1->SetFillColorAlpha(0,0.0);
  sp1->SetLineColorAlpha(0,0.0);
  sp1->Draw("NB");
  
  leg0->Draw();
  c0->SaveAs("plots/UE/rho_debug/CHGdNdetadphi_te_pt_systematics.pdf","PDF");
  c0->Close();
  

  TCanvas *c1 = new TCanvas;
  c1->SetMargin(0.125,0.05,0.125,0.05);

  double yEdge1[2] = { 0.5, 0.8 };
  TH2D *hscale3 = new TH2D("hscale3",";UE-subtracted lead jet p_{T} [GeV];<p_{T}^{ch}> [GeV]", nPtBins,ptBinEdge,1,yEdge1);
  hscale3->GetYaxis()->SetTitleOffset(1.25);
  hscale3->GetXaxis()->SetTitleOffset(1.25);

  TH1D *hChg_pt_te[nEAbins][nEtaBins];
  TH1D *hChg_pt_te_sys[nEAbins][nEtaBins];
  hscale3->SetName(" ");
  hscale3->SetLineWidth(0);
  hscale3->SetStats(0);
  hscale3->Draw();

  TLegend *leg1;
  
  leg1 = new TLegend(0.5, 0.65, 0.92, 0.92,NULL,"brNDC");    // LEGEND 0
  leg1->SetBorderSize(0);  leg1->SetLineColor(1);  leg1->SetLineStyle(1);
  leg1->SetLineWidth(1);  leg1->SetFillColor(0);  leg1->SetFillStyle(1001);
  leg1->SetNColumns(2);

  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hChg_MeanPt_te_" + EAbinName[a] + etaBinName[e];
      title = EAbinString[a] + "   "  + eastmidwest[e] + "   ";
      hChg_pt_te[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hChg_pt_te[a][e]->SetStats(0);
      hChg_pt_te[a][e]->SetLineColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt_te[a][e]->SetMarkerSize( 2 );

      leg1->AddEntry( hChg_pt_te[a][e], title, "lpf" );
      
      name = "hChg_MeanPt_te_sys_" + EAbinName[a] + etaBinName[e];
      hChg_pt_te_sys[a][e] = new TH1D( name, "" , nPtBins,ptBinEdge );

      hChg_pt_te_sys[a][e]->SetStats(0);
      hChg_pt_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 );
      hChg_pt_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      hChg_pt_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      hChg_pt_te_sys[a][e]->SetMarkerSize( 0.0 );

      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	double value = te_meanChgPt[e][a][p];
	hChg_pt_te[a][e]->SetBinContent( pbin, value );
	hChg_pt_te[a][e]->SetBinError( pbin, te_meanChgPt_err[e][a][p] );

	hChg_pt_te_sys[a][e]->SetBinContent( pbin, value );
	hChg_pt_te_sys[a][e]->SetBinError( pbin, value*te_meanChgPt_sys[e][a][p] );
	
      }
    }
  }

  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hChg_pt_te[a][e]->Draw("lpfX0ESAME");
      hChg_pt_te_sys[a][e]->Draw("E2SAME");
    }
  }
  
  leg1->Draw();


  sp1->Draw("NB");
 
  c1->SaveAs("plots/UE/rho_debug/Chg_MeanPt_te_pt_systematics.pdf","PDF");

 

}
