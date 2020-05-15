// Veronica Verkest
// May 8, 2020

void pAu2015HTsystematics(){

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
  TH1D *hBBCEastSum_byEta[nEtaBins][nFiles];
  TH1D *hBBCEastSum_byEta_rat[nEtaBins][nFiles];
  
  for ( int j=0; j<nFiles; ++j ) {
    
    name = "hRho_"; name += j;
    hRho[j]= new TH1D( name ,"Underlying Event;#rho (GeV)",20,0,10);
    name = "hRhoRat_"; name += j;
    hRho_rat[j]= new TH1D( name ,"",20,0,10);
    
    for ( int e=0; e<nEtaBins; ++e ) {

      name = "hBBCEastSum_byEta" + etaBinName[e] + "_"; name+= j;
      title = "BBC ADC East Sum:  " + etaBinString[e] + ";BBC East Sum";
      hBBCEastSum_byEta[e][j] = new TH1D( name, title, 18,0,70000 );
      hBBCEastSum_byEta[e][j]->SetMarkerStyle( etaMarker[e] );
      hBBCEastSum_byEta[e][j]->SetMarkerColor( etaColor[e] );
      hBBCEastSum_byEta[e][j]->SetLineColor( etaColor[e] );
      //hBBCEastSum_byEta[e][j]->SetFillColor( etaColor[e] );

      name = "hBBCEastSum_byEta_ratio" + etaBinName[e] + "_"; name+= j;
      title = "BBC ADC East Sum Ratio:  " + etaBinString[e] + ";BBC East Sum";
      hBBCEastSum_byEta_rat[e][j] = new TH1D( name, title, 18,0,70000 );
      hBBCEastSum_byEta_rat[e][j]->SetMarkerColor( etaColor[e] );
      hBBCEastSum_byEta_rat[e][j]->SetLineColor( etaColor[e] );
      //hBBCEastSum_byEta_rat[e][j]->SetFillColor( etaColor[e] );
    }
  }

 
  const int nSys = 4; // total, track, tower, hadcorr (unfolding soon...)
  TString sysName[nSys] = { "total", "track", "tower", "hadcorr" };

  int nEntries = jetTree[0]->GetEntries();

  
  for ( int i=0; i<nEntries; ++i ) {
    for ( int j=0; j<nFiles; ++j ) {

      jetTree[j]->GetEntry(i);

      pval = 99;    jeval = 99;    eaval = 99;
      if ( leadPt[j] >= 10.0 && leadPt[j] <= 30.0 ) {
	for ( int ea=0; ea<3; ++ea ) {
	  if ( BbcAdcSumEast[j] > BBCEsumLo[ea]  &&  BbcAdcSumEast[j] < BBCEsumHi[ea] ) { eaval = ea; }
	}    
	for ( int p=0; p<3; ++p ) {
	  if ( leadPt[j] >= ptLo[p]  &&  leadPt[j] <= ptHi[p] ) { pval = p; }
	}
	for ( int je=0; je<3; ++je ) {
	  if ( leadEta[j] >= etaLo[je]  &&  leadEta[j] <= etaHi[je] ) { jeval = je; }
	}
	if ( eaval==99 ) { continue; }
	if ( pval==99 || jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<leadPt[j]<<endl<<endl; }
      }
      else {continue;}

    

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
	  hBBCEastSum_byEta[jeval][j]->Fill( BbcAdcSumEast[j] );
	  hBBCEastSum_byEta_rat[jeval][j]->Fill( BbcAdcSumEast[j] );
	}
      }
      else {
	hRho[j]->Fill( rho_te[j] );
	hRho_rat[j]->Fill( rho_te[j] );
	hBBCEastSum_byEta[jeval][j]->Fill( BbcAdcSumEast[j] );
	hBBCEastSum_byEta_rat[jeval][j]->Fill( BbcAdcSumEast[j] );
      }


      
    }  
  }

  for ( int j=0; j<nFiles; ++j ) {
    hRho[j]->Scale(1./hRho[j]->GetEntries());
    for ( int je=0; je<nEtaBins; ++je ) {
      hBBCEastSum_byEta[je][j]->Scale(1./hBBCEastSum_byEta[je][j]->GetEntries());
    }
  }

  for ( int j=1; j<nFiles; ++j ) {
    hRho[j]->Divide( hRho[0] );
    for ( int je=0; je<nEtaBins; ++je ) {
      hBBCEastSum_byEta[je][j]->Divide( hBBCEastSum_byEta[je][0] );
    }
  }


  TH1D* hRho_sys[nSys];
  TH1D* hRho_sys_ratio[nSys];
  TH1D* hEAeta_sys[nSys][nEtaBins];
  TH1D* hEAeta_sys_rat[nSys][nEtaBins];

  for ( int i=0; i<nSys; ++i ) {
    name = "hRho_"; name += sysName[i];
    hRho_sys[i] = new TH1D( name ,"",20,0,10);
    name = "hRho_"; name += sysName[i]; name += "_ratio";
    hRho_sys_ratio[i] = new TH1D( name ,"",20,0,10);

    for ( int je=0; je<nEtaBins; ++je ) {

      name = "hEAeta_"; name += sysName[i] + etaBinName[je];
      hEAeta_sys[i][je] = new TH1D( name, "", 18, 0, 70000 ); 
      name = "hEAeta_"; name += sysName[i]; name += "_ratio" + etaBinName[je];
      hEAeta_sys_rat[i][je] = new TH1D( name, "", 18, 0, 70000 ); 
    }
    
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

	for (int je=0; je<nEtaBins; ++je) {
	  maxrat = fabs( 1 - hBBCEastSum_byEta[je][fileNo]->GetBinContent(i) );

	  hEAeta_sys_rat[j][je]->SetBinContent(i,maxrat);
	  
	  hEAeta_sys[j][je]->SetBinContent(i,hBBCEastSum_byEta[je][0]->GetBinContent(i));
	  hEAeta_sys[j][je]->SetBinError(i,maxrat*hBBCEastSum_byEta[je][0]->GetBinContent(i));
	}
	
	j = 0;
  	fileNo = 1;
      }
      else {
  	double sys1 = fabs( 1 - hRho[fileNo]->GetBinContent(i) );
  	double sys2 = fabs( 1- hRho[fileNo+1]->GetBinContent(i) );
  	double maxrat = max( sys1, sys2 );
  	hRho_sys_ratio[j]->SetBinContent(i,maxrat);
	
  	hRho_sys[j]->SetBinContent(i,hRho[0]->GetBinContent(i));    
  	hRho_sys[j]->SetBinError(i,maxrat*hRho[0]->GetBinContent(i));

	for ( int je=0; je<nEtaBins; ++je ) {
	  sys1 = fabs( 1 - hBBCEastSum_byEta[je][fileNo]->GetBinContent(i) );
	  sys2 = fabs( 1- hBBCEastSum_byEta[je][fileNo+1]->GetBinContent(i) );
	  maxrat = max( sys1, sys2 );
	  hEAeta_sys_rat[j][je]->SetBinContent(i,maxrat);
	  
	  hEAeta_sys[j][je]->SetBinContent(i,hBBCEastSum_byEta[je][0]->GetBinContent(i));
	  hEAeta_sys[j][je]->SetBinError(i,maxrat*hBBCEastSum_byEta[je][0]->GetBinContent(i));
	}
	
  	++fileNo;
  	++fileNo;
      }


      if (j==0) {
      	double  quadsum = 0;
      	for (int k=1; k<nSys; ++k){
      	  quadsum += (hRho_sys_ratio[k]->GetBinContent(i))*(hRho_sys_ratio[k]->GetBinContent(i));
      	}
      	hRho_sys_ratio[j]->SetBinContent( i, sqrt(quadsum) );
	
      	hRho_sys[j]->SetBinContent(i,hRho[0]->GetBinContent(i));    
      	hRho_sys[j]->SetBinError(i, sqrt(quadsum)*hRho[0]->GetBinContent(i));

	for (int je=0; je<nEtaBins; ++je) {
	  quadsum = 0;
	  for (int k=1; k<nSys; ++k){
	    quadsum += (hEAeta_sys_rat[k][je]->GetBinContent(i))*(hEAeta_sys_rat[k][je]->GetBinContent(i));
	  }
	  
	  hEAeta_sys_rat[0][je]->SetBinContent( i, sqrt(quadsum) );
	  
	  hEAeta_sys[0][je]->SetBinContent(i,hBBCEastSum_byEta[je][0]->GetBinContent(i));    
	  hEAeta_sys[0][je]->SetBinError(i, sqrt(quadsum)*hBBCEastSum_byEta[je][0]->GetBinContent(i));
	}
	
	j=999;
      }
      
    }
    
  }

  gStyle->SetOptStat(0);

  
  int sysColor[nSys] = { 875, 856, 835, 798/*, 896*/ };
  int sysStyle[nSys] = { 3003, 3315, 3351, 3007/*, 3007*/ };

  TCanvas * c_rho_sys = new TCanvas( "c_rho_sys" , "" ,700 ,500 );              // CANVAS
  //c_rho_sys->SetLogy();

  TH2D *sRho_sys = new TH2D( "sRho_sys", "", 10,0.0,10.0, 50,0.0,2.0 );
  c_rho_sys->SetMargin(0.15,0.1,0.1,0.1);
  sRho_sys->GetYaxis()->SetRangeUser(0,0.35);
  sRho_sys->SetTitle(";<#rho> (GeV);#frac{1}{N_{jets}} #frac{dN_{jets}}{d#rho}");
  sRho_sys->GetYaxis()->SetTitleOffset(2);
  sRho_sys->Draw();

  // for (int i=0; i<nSys; ++i) {
    
  //   //hRho_sys[i]->SetFillColorAlpha(sysColor[i],0.3);
  //   hRho_sys[i]->SetFillColor(sysColor[i]);
  //   //hRho_sys[i]->SetFillStyle(sysStyle[i]);
  //   hRho_sys[i]->Draw("e2same");
  // }

  hRho_sys[0]->SetFillColor(430);
  hRho_sys[0]->SetLineColor(430);
  hRho_sys[0]->SetMarkerColor(430);
  hRho_sys[0]->Draw("e2same");

  hRho[0]->SetLineColor(kBlack);
  hRho[0]->SetMarkerColor(kBlack);
  //hRho[0]->SetLineWidth(2);
  //hRho[0]->SetMarkerStyle(20);
  hRho[0]->Draw("same");

  c_rho_sys->BuildLegend();
  c_rho_sys->SaveAs("plots/UE/systematics/rhoDist.pdf","PDF");


  TCanvas * c_rho_sys_ratio = new TCanvas( "c_rho_sys_ratio" , "" ,700 ,300 );              // CANVAS
  c_rho_sys_ratio->SetTopMargin(0);
  c_rho_sys_ratio->SetLeftMargin(0.15);
  c_rho_sys_ratio->SetBottomMargin(0.20);

  //c_rho_sys_ratio->SetLogy(0);
  sRho_sys->SetTitle(";<#rho> (GeV);Systematic Uncertainty       ");
  sRho_sys->GetYaxis()->SetTitleOffset(0.6);
  sRho_sys->GetYaxis()->SetLabelSize(0.05);
  sRho_sys->GetYaxis()->SetTitleSize(0.06);
  sRho_sys->GetXaxis()->SetTitleOffset(1);
  sRho_sys->GetXaxis()->SetLabelSize(0.05);
  sRho_sys->GetXaxis()->SetTitleSize(0.06);
  sRho_sys->Draw();
  sRho_sys->GetYaxis()->SetRangeUser(0,1.1);
  for (int i=0; i<nSys; ++i) {
    // hRho_sys_ratio[i]->SetFillColorAlpha(sysColor[i],0.3);
    hRho_sys_ratio[i]->SetFillColor(sysColor[i]);
    hRho_sys_ratio[i]->SetLineColor(sysColor[i]);
    hRho_sys_ratio[i]->SetLineWidth(2);
    hRho_sys_ratio[i]->SetMarkerColor(sysColor[i]);
    hRho_sys_ratio[i]->SetFillStyle(sysStyle[i]);
    hRho_sys_ratio[i]->Draw("e3same");
  }

  c_rho_sys_ratio->BuildLegend(0.19,0.50,0.49,0.89);
  c_rho_sys_ratio->SaveAs("plots/UE/systematics/rhoDistSys.pdf","PDF");




  // TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 1
  // c1->SetMargin(0.15,0.1,0.12,0.1);
  // TH2D *sBBCbyEta = new TH2D("sBBCbyEta", "BBC ADC East Sum by Lead Jet #eta;BBC East Sum;#frac{1}{N_{jets}} #frac{dN_{jets}}{d<BBCEsum>}", 35,0,70000, 10,0,0.15);
  // sBBCbyEta->SetStats(0);
  // sBBCbyEta->GetXaxis()->SetLabelSize(.04);
  // sBBCbyEta->GetYaxis()->SetLabelSize(.04);
  // sBBCbyEta->GetXaxis()->SetTitleOffset(1.2);
  // sBBCbyEta->GetYaxis()->SetTitleOffset(1.6);
  // sBBCbyEta->Draw();

  
  // TLegend *leg0 = new TLegend(0.55, 0.55, 0.89, 0.89,NULL,"brNDC");    // LEGEND 0
  // leg0->SetBorderSize(0);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  // leg0->SetNColumns(2);
  // leg0->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
  // leg0->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");
  
  // for ( int e=0; e<nEtaBins; ++e ) {

  //   hEAeta_sys[0][e]->SetFillColorAlpha(etaColor[e],0.50);
  //   hEAeta_sys[0][e]->Draw("e2same");
    
  //   hBBCEastSum_byEta[e][0]->SetStats(0);
  //   //hBBCEastSum_byEta[e][0]->Scale(1./hBBCEastSum_byEta[e]->Integral());
  //   hBBCEastSum_byEta[e][0]->Draw("SAME");
  //   avg = "";
  //   avg += hBBCEastSum_byEta[e][0]->GetMean(1);                                           // 1 denotes x-axis
  //   avg = avg(0,5);
  //   //name = "hBBCEastSum_byEta" + etaBinName[e];
  //   name = "hBBCEastSum_byEta" + etaBinName[e] + "_"; name+= 0;
  //   title = etaBinString[e];
  //   leg0->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
  //   leg0->AddEntry((TObject*)0,avg, "");

  // }
  // leg0->Draw();
  // c1->SaveAs( "plots/UE/systematics/BBCEastSum_by_eta.pdf" , "PDF" );
  // c1->SetLogy(0);
  
}
