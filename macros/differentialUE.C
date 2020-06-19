// Veronica Verkest
// May 14, 2020

void differentialUE(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

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
  TString EAbinString[nEAbins] = { "Lo", "Hi" };
  TString BBCselection[nEAbins] = { "BbcAdcSumEast>3559.12 && BbcAdcSumEast<11503", "BbcAdcSumEast>26718.1" };
 
  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho_te+neuEastRho)", "(chgMidRho_te+neuMidRho)", "(chgWestRho_te+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  const double eastArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  double area[nEtaBins] = { eastArea, midArea, westArea };
  
  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 23, 22 };
  
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
    
      name = "hRho_"; name += EAbinString[i]; name += ptBinName[p];
      hRho[i][p] = new TH1D(name,"UE;<#rho> (GeV)", 30,0.0,15.0);
      
      name = "hNchg_"; name += EAbinString[i]; name += ptBinName[p];
      hNchg[i][p] = new TH1D( name,"# charged particles in UE", 50,0.0,50.0);

      name = "hNchg_te_"; name += EAbinString[i]; name += ptBinName[p];
      hNchg_te[i][p] = new TH1D( name,"# track effic. corr. charged particles in UE", 50,0.0,50.0);
      
      name = "hNneu_"; name += EAbinString[i]; name += ptBinName[p];
      hNneu[i][p] = new TH1D( name,"# neutral particles in UE", 50,0.0,50.0);
    }


  }


  for (int a=0; a<nEAbins; ++a) {
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
	//hNchg_te[a][pval]->Fill(nBGpart_chg[a]);
	hNneu[a][pval]->Fill(nBGpart_neu[a]);


      }
      
    }
  }

  
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {

      hChgUEpt[a][p]->Scale(1./hRho[a][p]->GetEntries());
      // cout<<hChgUEpt[a][p]->Integral()<<endl;
      // cout<<hRho[a][p]->GetEntries()<<endl<<endl;
      
      hNeuUEpt[a][p]->Scale(1./hRho[a][p]->GetEntries());

    }
  }
  
  TH1D *hChgPtDist[ybins][nEAbins][nPtBins];
  
  int ecolor[ybins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  TString ybinEdgeString[ybins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };
  
  TCanvas * c0[nPtBins];

  TH1D *hTreeChgRho[nEAbins][nPtBins];
  TH1D *hChgDist[nEAbins][nPtBins];

  TLegend *leg0[nPtBins];
  
  for (int p=0; p<nPtBins; ++p) {
    name = "c0_"; name += p;
    c0[p] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 0
    c0[p]->SetLogy();
    for (int a=0; a<nEAbins; ++a) {
      for (int i=0; i<ybins; ++i) {
	int binno = i+1;
	name = "hChgpTdist"; name += ptBinName[p]; name += EAbinString[a]; name += "EA_etabin"; name += binno;
	hChgUEpt[a][p]->GetYaxis()->SetRange( binno, binno+1 );
	hChgUEpt[a][p]->SetLineColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerStyle( ptMarker[p] );	
	hChgPtDist[i][a][p] = (TH1D*) hChgUEpt[a][p]->ProjectionX( name, binno, binno+1 );
	title = "Charged UE p_{T} dist. for "; title += ptBinString[p]; title += " (";
	title += ybinEdgeString[i]; title +="<#eta<"; title += ybinEdgeString[binno]; title += ");p_{T} (GeV)";
	hChgPtDist[i][a][p]->SetTitle( title );
	hChgPtDist[i][a][p]->GetYaxis()->SetRangeUser(0.000001,0.3);
	if (i==0) { hChgPtDist[i][a][p]->Draw(""); }
	else { hChgPtDist[i][a][p]->Draw("SAME"); }
      }
    }

    leg0[p] = (TLegend*) c0[p]->BuildLegend(0.4,0.45,0.75,0.8);
    leg0[p]->SetFillColorAlpha(0,0);  leg0[p]->SetLineWidth(0);
    //c0[p]->Close();
  }


  //  TRACKING EFFICIENCY CORRECTIONS!
  TH2D *hChgUEpt_te[nEAbins][nPtBins];

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinString[a];  // rename these!
      hChgUEpt[a][p]->SetName( name );

      name = "hNeuUEpt"; name += ptBinName[p]; name += EAbinString[a];
      hNeuUEpt[a][p]->SetName( name );
      
      hChgUEpt_te[a][p] = (TH2D*) hChgUEpt[a][p]->Clone();       // create clone of 2D pt distributions
      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinString[a]; name += "_te";
      hChgUEpt_te[a][p]->SetName( name );

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
	  //if ( pt > 3.0 ) { pt = 3.0; }
	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  double eff = hEff[a][fileNo]->GetBinContent(ptbin);
	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double corr_value = (double) old_value/eff;
	  double corr_err = (double) old_err/eff;
	  hChgUEpt_te[a][p]->SetBinContent( xbin, ybin, corr_value);
	  hChgUEpt_te[a][p]->SetBinError( xbin, ybin, corr_err );
	}
	
      }
	
    }
  }



  TCanvas *c1[nPtBins];
  TLegend *leg[nPtBins];
  
  for (int j=0; j<nEAbins; ++j) {

    for (int p=0; p<nPtBins; ++p) {

      name = "c1_"; name += p;
      c1[p] = new TCanvas( name , "" ,700 ,500 );
      c1[p]->SetLogy();

      name = "hTreeChgRho_"; name += EAbinString[j];  name += ptBinName[p];
      hTreeChgRho[j][p] = new TH1D(name,"",30,0,15);
      drawString = "((0.35*chgEastRho)+(0.3*chgMidRho)+(0.35*chgWestRho))>>"; drawString += name;
      jetTree[j]->Draw( drawString , ptSelection[p] );
      //cout<<jetTree[j]->GetName()<<"->Draw( "<<drawString<<" , "<<ptSelection[p]<<" )"<<endl;;
      hTreeChgRho[j][p]->Scale(1./hTreeChgRho[j][p]->GetEntries());

      hChgDist[j][p] = (TH1D*)hChgPtDist[0][j][p]->DrawCopy();
      name = "hChgDist_"; name += EAbinString[j];  name += ptBinName[p];
      title = EAbinString[j]; title += "EA: "; title += "Charged #rho (p_{T}) distribution for ";
      title += ptBinString[p]; title += ";#rho or p_{T} (GeV)";
      hChgDist[j][p]->SetNameTitle( name , title );
      for (int i=1; i<ybins; ++i) {
	hChgDist[j][p]->Add( hChgPtDist[i][j][p] );
      }
      hChgDist[j][p]->SetStats(0);


      leg[p] = new TLegend(0.4, 0.6, 0.89, 0.89,NULL,"brNDC");    // LEGEND 0
      leg[p]->SetBorderSize(0); leg[p]->SetLineColor(1); leg[p]->SetLineStyle(1); leg[p]->SetLineWidth(1); leg[p]->SetFillColor(0); leg[p]->SetFillStyle(1001);
      leg[p]->SetNColumns(3);
      leg[p]->AddEntry((TObject*)0,"name", "");
      leg[p]->AddEntry((TObject*)0,"#bf{<#rho_{chg}>}", "");
      leg[p]->AddEntry((TObject*)0,"#bf{<#frac{dN_{chg}}{dp_{T}d#eta}>}", "");

      hChgDist[j][p]->SetAxisRange(0.000005,3.0, "Y");
      hChgDist[j][p]->Draw();
      hTreeChgRho[j][p]->SetStats(0);
      hTreeChgRho[j][p]->SetAxisRange(0.000005,3.0, "Y");
      hTreeChgRho[j][p]->Draw("SAME");
      //c1[p]->BuildLegend();
	
      title = "Avg. UE p_{T} dist.";  //title += ptBinString[p];
      leg[p]->AddEntry( hChgDist[j][p], title , "lpf" );                            // ADD TO LEGEND
      avg = "";
      avg += hChgDist[j][p]->GetMean(1)*(hChgDist[j][p]->Integral()/AREA);                       // 1 denotes x-axis
      avg = avg(0,5);
      leg[p]->AddEntry((TObject*)0,avg, "");
      avg = "";
      avg += hChgDist[j][p]->Integral()/AREA;//hRho[j]->Integral();                       // 1 denotes x-axis
      avg = avg(0,5);
      leg[p]->AddEntry((TObject*)0,avg, "");

      title = "Event-by-event #rho dist.";  //title += ptBinString[p];
      leg[p]->AddEntry( hTreeChgRho[j][p], title, "lpf" );                            // ADD TO LEGEND
      avg = "";
      avg += hTreeChgRho[j][p]->GetMean(1);                                           // 1 denotes x-axis
      avg = avg(0,5);
      leg[p]->AddEntry((TObject*)0,avg, "");
      avg = "";
      avg += hNchg[j][p]->GetMean()/AREA;                       // 1 denotes x-axis
      avg = avg(0,5);
      leg[p]->AddEntry((TObject*)0,avg, "");
    
      leg[p]->Draw();
      saveName = "plots/UE/rho_debug/UEptMeanCompare_"; saveName += ptBinName[p]; saveName += EAbinString[j]; saveName += "EA.pdf";
      c1[p]->SaveAs( saveName ,"PDF");
      c1[p]->Close();
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

  double chgRho[nEtaBins][nEAbins][nPtBins];
  
  for ( int a=0; a<nEAbins; ++a) {
    for ( int p=0; p<nPtBins; ++p) {

      for ( int e=0; e<nEtaBins; ++e) {

	hChgUEpt[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	meanChgPt[e][a][p] = hChgUEpt[a][p]->GetMean(1);
	meanChgPt_err[e][a][p] = hChgUEpt[a][p]->GetMeanError(1);

	nChg[e][a][p] = hChgUEpt[a][p]->Integral();
	nChg_err[e][a][p] = hChgUEpt[a][p]->GetMeanError();
	
	chgRho[e][a][p] = ( hChgUEpt[a][p]->GetMean(1) * hChgUEpt[a][p]->Integral() ) / area[e];


	//  WITH TRACKING EFFICIENCY
	hChgUEpt_te[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	te_meanChgPt[e][a][p] = hChgUEpt_te[a][p]->GetMean(1);
	te_meanChgPt_err[e][a][p] = hChgUEpt_te[a][p]->GetMeanError(1);

	nChg_te[e][a][p] = hChgUEpt_te[a][p]->Integral();
	nChg_te_err[e][a][p] = hChgUEpt_te[a][p]->GetMeanError();


      }      
    }
  }

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


  double binEdge[nEtaBins+1] = { -0.5, 0.5, 1.5, 2.5 };
    
  TH2D *hscale0 = new TH2D("hscale0","", 10.0,-0.5,2.5, 10,0.2,0.8);
  hscale0->SetLineColorAlpha(1,0.0);
  hscale0->SetName(" ");
  hscale0->GetYaxis()->SetTitleSize(0.06);
  hscale0->GetYaxis()->SetTitleOffset(0.6);
  hscale0->GetXaxis()->CenterTitle();
  hscale0->GetXaxis()->SetLabelSize(0.0);
  hscale0->GetXaxis()->SetTitleSize(0.06);
  hscale0->GetXaxis()->SetTitleOffset(0.4);
  TCanvas *c2 = new TCanvas;
  hscale0->SetStats(0);
  hscale0->SetTitle(";east #rho       mid #rho       west #rho;<#rho_{chg}> (GeV);lead p_{T} (GeV)");
  hscale0->Draw(); 

  TH1D *hUE[nEAbins][nPtBins];
  
  for ( int a=0; a<nEAbins; ++a) {
    for ( int p=0; p<nPtBins; ++p) {

      name = "hUE_" + EAbinString[a] + ptBinName[p];
      title = EAbinString[a] + " EA   " + ptBinString[p];
      hUE[a][p] = new TH1D( name , title ,nEtaBins,binEdge);
      hUE[a][p]->SetStats(0);

      hUE[a][p]->SetLineColor( ptColor[p] );
      hUE[a][p]->SetMarkerColor( ptColor[p] );
      hUE[a][p]->SetMarkerStyle( EAmarker[a] );
      hUE[a][p]->SetMarkerSize( 2 );
      
      for ( int e=0; e<nEtaBins; ++e) {
	int ebin = e+1;
	double value = (double) ( nChg[e][a][p] * meanChgPt[e][a][p] ) / area[e];
	hUE[a][p]->SetBinContent( ebin, value );
	double errval = sqrt( (meanChgPt_err[e][a][p]*meanChgPt_err[e][a][p]) + ( nChg_err[e][a][p]*nChg_err[e][a][p]) );
	hUE[a][p]->SetBinError( ebin, errval );
	//cout<<value<<"\t";
      }
      hUE[a][p]->Draw("lpfSAME");
    }
  }
  c2->BuildLegend(0.65,0.7,1.0,1.0);
  c2->SaveAs("plots/UE/rho_debug/hiloptUE.pdf","PDF");
  
  hscale0->SetTitle("Track Effic. Corr.;east #rho       mid #rho       west #rho;track effic. corr. <#rho_{chg}> (GeV);lead p_{T} (GeV)");
  hscale0->Draw();
  hscale0->SetName(" ");

  TH1D *hUE_te[nEAbins][nPtBins];
  
  for ( int a=0; a<nEAbins; ++a) {
    for ( int p=0; p<nPtBins; ++p) {

      name = "hUE_te_" + EAbinString[a] + ptBinName[p];
      title = EAbinString[a] + " EA   " + ptBinString[p];
      hUE_te[a][p] = new TH1D( name ,title,nPtBins,binEdge);
      hUE_te[a][p]->SetStats(0);

      //hUE_te[a][p]->SetFillColorAlpha( ptColor[p], 0.20 );
      hUE_te[a][p]->SetLineColor( ptColor[p] );
      hUE_te[a][p]->SetMarkerColor( ptColor[p] );
      hUE_te[a][p]->SetMarkerStyle( EAmarker[a] );
      hUE_te[a][p]->SetMarkerSize( 2 );
      
      for ( int e=0; e<nEtaBins; ++e) {
	int ebin = e+1;
	double value = (double) ( nChg_te[e][a][p] * te_meanChgPt[e][a][p] ) / area[e];
	hUE_te[a][p]->SetBinContent( ebin, value );
	double errval = sqrt( (te_meanChgPt_err[e][a][p]*te_meanChgPt_err[e][a][p]) + ( nChg_te_err[e][a][p]*nChg_te_err[e][a][p]) );
	hUE_te[a][p]->SetBinError( ebin, errval );
	//cout<<value<<"\t";
      }
      hUE_te[a][p]->Draw("lpfSAME");
    }
  }
  c2->BuildLegend(0.65,0.7,1.0,1.0);
  c2->SaveAs("plots/UE/rho_debug/hiloptUE_te.pdf","PDF");


  double ptBinEdge[nPtBins+1] = { 10.0, 15.0, 20.0, 30.0 };
  double yEdge[2] = {0.4,1.4};
  
  TH2D *hscale1 = new TH2D("hscale1",";lead p_{T} (GeV);<#frac{dN_{chg}}{d#eta d#phi}>", nPtBins,ptBinEdge,1,yEdge);
  hscale1->SetName("");
  TCanvas *c3 = new TCanvas;  

  TH1D *hCHGdNdetadphi[nEAbins][nEtaBins];
  hscale1->Draw();
  hscale1->SetName(" ");
  hscale1->SetLineWidth(0);
  hscale1->SetStats(0);

  for ( int a=0; a<nEAbins; ++a) {
    for ( int e=0; e<nEtaBins; ++e) {

      name = "hCHGdNdetadphi_" + EAbinString[a] + etaBinName[e];
      title = "hCHGdNdetadphi_" + EAbinString[a] + etaBinName[e];
      hCHGdNdetadphi[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hCHGdNdetadphi[a][e]->SetStats(0);
      hCHGdNdetadphi[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi[a][e]->SetMarkerSize( 2 );
	    
      //hCHGdNdetadphi[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	double value = (double) ( nChg[e][a][p] / area[e] );
	hCHGdNdetadphi[a][e]->SetBinContent( pbin, value );
	hCHGdNdetadphi[a][e]->SetBinError( pbin, nChg_err[e][a][p] );
	
      }
      hCHGdNdetadphi[a][e]->Draw("lpfSAME");
    }
  }
    
  c3->BuildLegend();
  c3->SaveAs("plots/UE/rho_debug/CHGdNdetadphi.pdf","PDF");  


  double yEdge2[2] = { 0.5, 0.8 };
  TH2D *hscale2 = new TH2D("hscale2",";lead p_{T} (GeV);Charged <p_{T}> (GeV)", nPtBins,ptBinEdge,1,yEdge2);
  hscale2->SetName("");



  TH1D *hChg_pt[nEAbins][nEtaBins];
  hscale2->Draw();
  hscale2->SetName(" ");
  hscale2->SetLineWidth(0);
  hscale2->SetStats(0);

  for ( int a=0; a<nEAbins; ++a) {
    for ( int e=0; e<nEtaBins; ++e) {

      name = "hChg_MeanPt_" + EAbinString[a] + etaBinName[e];
      title = "hChg_MeanPt_" + EAbinString[a] + etaBinName[e];
      hChg_pt[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hChg_pt[a][e]->SetStats(0);
      hChg_pt[a][e]->SetLineColor( etaColor[e] );
      hChg_pt[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt[a][e]->SetMarkerSize( 2 );
	    
      //hCHGdNdetadphi[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	double value = meanChgPt[e][a][p];
	hChg_pt[a][e]->SetBinContent( pbin, value );
	hChg_pt[a][e]->SetBinError( pbin, meanChgPt_err[e][a][p] );
	
      }
      hChg_pt[a][e]->Draw("lpfSAME");
    }
  }
    
  c3->BuildLegend(0.5,0.75,1.0,1.0);
  c3->SaveAs("plots/UE/rho_debug/Chg_MeanPt.pdf","PDF");








  TH2D *hscale1a = new TH2D("hscale1a","Track Effic. Corr.;lead p_{T} (GeV);track effic. corr. <#frac{dN_{chg}}{d#eta d#phi}>", nPtBins,ptBinEdge,1,yEdge);
  hscale1a->SetName("");
  TCanvas *c3a = new TCanvas;  

  TH1D *hCHGdNdetadphi_te[nEAbins][nEtaBins];
  hscale1a->Draw();
  hscale1a->SetName(" ");
  hscale1a->SetLineWidth(0);
  hscale1a->SetStats(0);

  for ( int a=0; a<nEAbins; ++a) {
    for ( int e=0; e<nEtaBins; ++e) {

      name = "hCHGdNdetadphi_te_te_" + EAbinString[a] + etaBinName[e];
      title = "hCHGdNdetadphi_te_te_" + EAbinString[a] + etaBinName[e];
      hCHGdNdetadphi_te[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hCHGdNdetadphi_te[a][e]->SetStats(0);
      hCHGdNdetadphi_te[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi_te[a][e]->SetMarkerSize( 2 );
	    
      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	double value = (double) ( nChg_te[e][a][p] / area[e] );
	hCHGdNdetadphi_te[a][e]->SetBinContent( pbin, value );
	hCHGdNdetadphi_te[a][e]->SetBinError( pbin, nChg_te_err[e][a][p] );
	
      }
      hCHGdNdetadphi_te[a][e]->Draw("lpfSAME");
    }
  }
    
  c3a->BuildLegend();
  c3a->SaveAs("plots/UE/rho_debug/CHGdNdetadphi_te.pdf","PDF");

  

  TCanvas *c3b = new TCanvas;  

  double yEdge3[2] = { 0.5, 0.8 };
  TH2D *hscale3 = new TH2D("hscale3",";lead p_{T} (GeV);Track Effic. Corr. Charged <p_{T}> (GeV)", nPtBins,ptBinEdge,1,yEdge3);


  TH1D *hChg_pt_te[nEAbins][nEtaBins];
  hscale3->SetName(" ");
  hscale3->SetLineWidth(0);
  hscale3->SetStats(0);
  hscale3->Draw();

  for ( int a=0; a<nEAbins; ++a) {
    for ( int e=0; e<nEtaBins; ++e) {

      name = "hChg_MeanPt_te_" + EAbinString[a] + etaBinName[e];
      title = "hChg_MeanPt_te_" + EAbinString[a] + etaBinName[e];
      hChg_pt_te[a][e] = new TH1D( name, title , nPtBins,ptBinEdge );

      hChg_pt_te[a][e]->SetStats(0);
      hChg_pt_te[a][e]->SetLineColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt_te[a][e]->SetMarkerSize( 2 );
	    
      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int pbin = p+1;
	double value = te_meanChgPt[e][a][p];
	hChg_pt_te[a][e]->SetBinContent( pbin, value );
	hChg_pt_te[a][e]->SetBinError( pbin, te_meanChgPt_err[e][a][p] );
	
      }
      hChg_pt_te[a][e]->Draw("lpfSAME");
    }
  }
    
  c3b->BuildLegend(0.5,0.75,1.0,1.0);
  c3b->SaveAs("plots/UE/rho_debug/Chg_MeanPt_te.pdf","PDF");



  
  TFile *outFile = new TFile("src/differentialUE.root","RECREATE");
  
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {

      hChgUEpt[a][p]->Write();
      hChgUEpt_te[a][p]->Write();
      hNeuUEpt[a][p]->Write();
	
    }
  }
    
  outFile->Write();
  //outFile->Close();

}
