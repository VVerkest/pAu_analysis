// Veronica Verkest
// May 8, 2020

void pAu2015HTsystematics_bad(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const TString ptCorrectedBinString[nPtBins] = { "10<p_{T}^{corrected}<15", "15<p_{T}^{corrected}<20",  "20<p_{T}^{corrected}<30" };
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
  
  const int nEA = 3;
  const TString EAstring[nEA] = { "low", "mid", "high" };
  const double BBCEsumLo[nEA] = { 3559.12, 10126.1, 26718.1 };   // LO: 3559.12-10126.1;  HI: 26718.1+
  const double BBCEsumHi[nEA] = { 10126.1, 26718.1, 64000.0 };  
  // { 0, 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000 }

  const int nFiles = 6;
  TFile* inFile[nFiles];
  TString fileName[nFiles]= { "out/UE/pAuHTjetUE_inLead.root", "out/UE/pAuHTjetUE_trackeffic_sys1.root", "out/UE/pAuHTjetUE_trackeffic_sys2.root", "out/UE/pAuHTjetUE_towereffic_sys1.root", "out/UE/pAuHTjetUE_towereffic_sys2.root", "out/UE/pAuHTjetUE_hadcorr_sys.root" };
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma;

  //  Tree variables
  int RunID[nFiles], EventID[nFiles], nTowers[nFiles], nPrimary[nFiles], nGlobal[nFiles], nVertices[nFiles], refMult[nFiles],
    gRefMult[nFiles], nBGpart_chg[nFiles], nBGpart_neu;
  double Vz[nFiles], BbcAdcSumEast[nFiles], leadPt[nFiles], leadEta[nFiles], leadPhi[nFiles], chgEastRho[nFiles], chgMidRho[nFiles], chgWestRho[nFiles], neuEastRho[nFiles], neuMidRho[nFiles], neuWestRho[nFiles], leadArea[nFiles], eastRho[nFiles], midRho[nFiles], westRho[nFiles], leadPtCorrected[nFiles], chgEastRho_te[nFiles], chgMidRho_te[nFiles], chgWestRho_te[nFiles], rho_te[nFiles], rho[nFiles];

  TTree *jetTree[nFiles];

  
  for ( int i=0; i<nFiles; ++i ) {
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
    jetTree[i]->SetBranchAddress( "rho", &rho[i] );
    jetTree[i]->SetBranchAddress( "rho_te", &rho_te[i] );
  }

  TH1D *hRho[nFiles];
  TH1D *hRho_rat[nFiles];
  for ( int j=0; j<nFiles; ++j ) {
    name = "hRho_"; name += j;
    hRho[j]= new TH1D( name ,"Underlying Event;#rho (GeV)",20,0,10);
    name = "hRhoRat_"; name += j;
    hRho_rat[j]= new TH1D( name ,"",20,0,10);
  }    

  const int nSys = 4; // total, track, tower, hadcorr (unfolding soon...)
  TString sysName[nSys] = { "total", "track", "tower", "hadcorr" };

  int nEntries = jetTree[0]->GetEntries();

  
  for ( int i=0; i<nEntries; ++i ) {
    
    for ( int j=0; j<nFiles; ++j ) {
      jetTree[j]->GetEntry(i);
      hRho[j]->Fill( rho_te[j] );
      hRho_rat[j]->Fill( rho_te[j] );

      if ( j==3 ) {
	int k=i;
	jetTree[j]->GetEntry(k);
	while ( !( RunID[j]==RunID[0] && EventID[j]==EventID[0] ) ) {
	  k++;
	  jetTree[j]->GetEntry(k);
	}
	if ( RunID[j]==RunID[0] && EventID[j]==EventID[0] ) {
	  hRho[j]->Fill( rho_te[j] );
	  hRho_rat[j]->Fill( rho_te[j] );
	}
      }
    }  
  }

  for ( int j=0; j<nFiles; ++j ) {
    hRho[j]->Scale(1./hRho[j]->GetEntries());
  }

  for ( int j=1; j<nFiles; ++j ) {
    name = "hRho_rat_"; name += j;
    title = "#rho ratio "; title += j;
    hRho[j]->Divide( hRho[0] );
  }


  TH1D* hRho_sys[nSys];
  TH1D* hRho_sys_ratio[nSys];

  for ( int i=0; i<nSys; ++i ) {
    name = "hRho_"; name += sysName[i];
    hRho_sys[i] = new TH1D( name ,"",20,0,10);
    name = "hRho_"; name += sysName[i]; name += "_ratio";
    hRho_sys_ratio[i] = new TH1D( name ,"",20,0,10);
  }


  for ( int i=1; i<(hRho[0]->GetNbinsX()+1); ++i ) {

    if ( hRho[0]->GetBinContent(i) == 0 ) { continue; }

    int fileNo = 1;
    for (int j=1; j<nSys; ++j) {

      if ( sysName[j]=="hadcorr" ) {
  	double maxrat = fabs( 1 - hRho[fileNo]->GetBinContent(i) );
  	hRho_sys_ratio[j]->SetBinContent(i,maxrat);
	
  	hRho_sys[j]->SetBinContent(i,hRho[0]->GetBinContent(i));    
  	hRho_sys[j]->SetBinError(i,maxrat*hRho[0]->GetBinContent(i));
  	cout<<maxrat*hRho[0]->GetBinContent(i)<<endl;
  	//j = 1;
  	fileNo = 1;
      }
      else {
  	double sys1 = fabs( 1 - hRho[fileNo]->GetBinContent(i) );
  	double sys2 = fabs( 1- hRho[fileNo+1]->GetBinContent(i) );
  	double maxrat = max( sys1, sys2 );
  	hRho_sys_ratio[j]->SetBinContent(i,maxrat);
	
  	hRho_sys[j]->SetBinContent(i,hRho[0]->GetBinContent(i));    
  	hRho_sys[j]->SetBinError(i,maxrat*hRho[0]->GetBinContent(i));
  	++fileNo;
  	++fileNo;
      }
      
    }
    
  }

  gStyle->SetOptStat(0);

  TCanvas * c_rho_sys = new TCanvas( "c_rho_sys" , "" ,700 ,500 );              // CANVAS
  c_rho_sys->SetLogy();

  TH2D *sRho_sys = new TH2D( "sRho_sys", "", 10,0,10, 10,0,2 );
  //sRho_sys->Draw();
  
  int sysColor[nSys] = { 875, 856, 835, 798/*, 896*/ };
  for (int i=0; i<nSys; ++i) {
    
    hRho_sys[i]->SetFillColor(sysColor[i]);
    hRho_sys[i]->Draw("e3same");
  }

  hRho[0]->SetLineColor(kBlack);
  hRho[0]->SetMarkerColor(kBlack);
  hRho[0]->SetLineWidth(2);
  hRho[0]->Draw("same");
  
  c_rho_sys->SaveAs("plots/UE/systematics/rhoDist.pdf","PDF");


  TCanvas * c_rho_sys_ratio = new TCanvas( "c_rho_sys_ratio" , "" ,700 ,500 );              // CANVAS

  c_rho_sys_ratio->SetLogy(0);
  //sRho_sys->Draw();
  for (int i=0; i<nSys; ++i) {
    hRho_sys_ratio[i]->SetFillColor(sysColor[i]);
    hRho_sys_ratio[i]->Draw("e3same");
  }
  c_rho_sys_ratio->SaveAs("plots/UE/systematics/rhoDistSys.pdf","PDF");

}
