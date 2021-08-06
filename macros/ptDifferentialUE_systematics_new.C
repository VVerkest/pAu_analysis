// Veronica Verkest
// September 13, 2020

void ptDifferentialUE_systematics_new(){

  //gStyle->SetErrorX(0.0001);  // gStyle->SetOptFit(1);
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

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
  const TString etaBinString[nEtaBins] = { "-0.6<#eta^{lead}_{jet}<-0.3", "-0.3<#eta^{lead}_{jet}<0.3", "0.3<#eta^{lead}_{jet}<0.6" };
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
  TString emw[nEtaBins] = { "east", "mid", "west" };
  TString rhoVal[nEtaBins] = { "(chgEastRho_te+neuEastRho)", "(chgMidRho_te+neuMidRho)", "(chgWestRho_te+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  const double eastArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  double area[nEtaBins] = { eastArea, midArea, westArea };
  
  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 24, 20 };
  
  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  string dirName = "prelim_new";
  // TString fileName[nEAbins] = { "out/UE/pAuHTjetUE_30cmVzCut_3cmVzDiff_loEA_variableBins.root", "out/UE/pAuHTjetUE_30cmVzCut_3cmVzDiff_hiEA_variableBins.root" };
  TString fileName[nEAbins] = { "out/UE/pAuHTjetUE_halfGeVbins_loEA_uncorrected.root", "out/UE/pAuHTjetUE_halfGeVbins_hiEA_uncorrected.root" };
  TString efficFileName[nEAbins] = { "src/trackeffic_loEA.root", "src/trackeffic_hiEA.root" };

  const int eta_bins = 10;
  double etabinEdge[eta_bins+1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0 };

  TFile* inFile[nEAbins];
  TFile* efficFile[nEAbins];
  // TF1* funcpp[nEAbins][eta_bins];
    
  for (int i=0; i<nEAbins; ++i) {
    efficFile[i] = new TFile( efficFileName[i], "READ" );
  }

  TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  TH1D *hEff[eta_bins][nEAbins];
  TF1* efficFit[eta_bins][nEAbins];

  TCanvas *can = new TCanvas;  

  for (int i=0; i<eta_bins; ++i) {
    for (int j=0; j<nEAbins; ++j){    
      int binno = i+1;
    
      if (j==0 ) { name = "eff_s_bin_1_3_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
      else if (j==1) {name = "eff_s_bin_7_10_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
    
      hEff[i][j] = (TH1D*)efficFile[j]->Get(name);

      name = EAbinName[j] + "EA_effic_etaBin"; name += i;

      hEff[i][j]->SetAxisRange(0.0,1.3,"Y");
      hEff[i][j]->SetAxisRange(0.0,15.0,"X");
      hEff[i][j]->Fit( "eff", "EM" );

      efficFit[i][j] = hEff[i][j]->GetFunction("eff");
      efficFit[i][j]->SetName(name);
      
      hEff[i][j]->Draw();
      
      saveName = "plots/UE/" + dirName + "/" + name + ".pdf";
      can->SaveAs(saveName,"PDF");

    }
  }

  const int xbins = 55;
  double xbinEdge[xbins+1];
  const int zbins = 20;
  double zbinEdge[zbins+1];  
  for (int i=0; i<xbins+1; ++i) {
    xbinEdge[i] = (double) i+4.0;
    if (i<zbins+1) { zbinEdge[i] = (double)(-10.0 + i)/10.0; }
  }
  // const int ybins = 14;
  // double ybinEdge[ybins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0 };
  const int ybins = 30;
  const double ybinEdge[ybins+1] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0 };
  
  TH3D *h_UE3D[nEAbins][nEtaBins];
  TH1D *h_leadPt[nEAbins][nEtaBins];
  
  TH2D *hChgUEpt[nEAbins][nPtBins];
  TH1D *hLeadPt[nEAbins][nPtBins];
  TH1D *hNchg[nEAbins][nPtBins];
  TH1D *hNchg_te[nEAbins][nPtBins];
  int nEntries[nEAbins];

  TH3D *hChgUE3D[nEAbins]; 
  TH1D *hALLleadPt[nEAbins];

  for (int i=0; i<nEAbins;++i){
    inFile[i] = new TFile( fileName[i], "READ" );

    name = "hChgUE3D_" + EAbinName[i]; //55,4.0,59.0, 30,0.0,15.0, 20,-1.0,1.0);//
    hChgUE3D[i] = new TH3D(name,";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);

    name = "hALLleadPt_" + EAbinName[i];
    hALLleadPt[i] = new TH1D(name,";leading jet p_{T} (GeV)", 55,4.0,59.0);
      
    for (int e=0; e<nEtaBins; ++e) {
      name = "hChgUE_" + emw[e] + "EtaJet";
      h_UE3D[i][e] = (TH3D*)inFile[i]->Get(name);
      hChgUE3D[i]->Add( h_UE3D[i][e] );
      name = "hLeadPt_" + emw[e] + "EtaJet";
      h_leadPt[i][e] = (TH1D*)inFile[i]->Get(name);
      hALLleadPt[i]->Add( h_leadPt[i][e] );
    }

    for (int p=0; p<nPtBins; ++p) {

      hChgUE3D[i]->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]); // X=leadPt, Y=UEpt, Z=UEeta
      hALLleadPt[i]->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]);
	
      name = "hChgUE" + EAbinName[i] + ptBinString[p];
      hChgUEpt[i][p] = (TH2D*)hChgUE3D[i]->Project3D( "ZY" );
      hChgUEpt[i][p]->SetName(name);
      
      name = "hLeadPt_"; name += EAbinName[i]; name += ptBinName[p];
      hLeadPt[i][p] = (TH1D*)hALLleadPt[i]->DrawClone(name);
      hLeadPt[i][p]->SetEntries(hALLleadPt[i]->Integral());
    }


  }


  
  for (int a=nEAbins-1; a>=0; --a) {
    for (int p=0; p<nPtBins; ++p) {

      hChgUEpt[a][p]->Scale(1./hLeadPt[a][p]->GetEntries());
      cout<<EAbinString[a]<<"  "<<ptBinString[p]<<"  "<<hLeadPt[a][p]->GetEntries()<<endl;

    }
  }
  
  TH1D *hChgPtDist[eta_bins][nEAbins][nPtBins];
  
  int ecolor[eta_bins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  TString etabinEdgeString[eta_bins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };
  
  TCanvas * cpt[nPtBins];

  TLegend *legpt[nPtBins];
  
  for (int p=0; p<nPtBins; ++p) {
    name = "cpt_"; name += p;
    cpt[p] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 0
    cpt[p]->SetLogy();
    for (int a=nEAbins-1; a>=0; --a) {
      for (int i=0; i<eta_bins; ++i) {
	int binno = i+1;
	name = "hChgpTdist"; name += ptBinName[p]; name += EAbinName[a]; name += "EA_etabin"; name += binno;
	hChgUEpt[a][p]->GetYaxis()->SetRange( binno, binno+1 );
	hChgUEpt[a][p]->SetLineColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerColor( ecolor[i] );
	hChgUEpt[a][p]->SetMarkerStyle( ptMarker[p] );	
	hChgPtDist[i][a][p] = (TH1D*) hChgUEpt[a][p]->ProjectionX( name, binno, binno+1 );
	title = "Charged UE p_{T} dist. for "; title += ptBinString[p]; title += " (";
	title += etabinEdgeString[i]; title +="<#eta<"; title += etabinEdgeString[binno]; title += ");p_{T} [GeV/#it{c}]";
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

	for (int j=0; j<eta_bins; ++j) {
	  if ( binCenter>etabinEdge[j] && binCenter<etabinEdge[j+1] ) { fileNo = j; };
	}
  	if ( fileNo==99 ) { cerr<<"CANNOT FIND EFFICIENCY HISTOGRAM"<<endl; }

	// cout<<fileNo<<"  "<<binCenter<<"  "<<hEff[fileNo][a]->GetTitle()<<endl;

  	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
  	  int xbin = ix;
	  
  	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
  	  // if ( pt > 3.0 ) { pt = 3.0; }
  	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  
  	  // double eff = hEff[a][fileNo]->GetBinContent(ptbin);
	  double eff = efficFit[a][fileNo]->Eval( hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin) );
	  // double eff_err = efficFit[a][fileNo]->Eval( hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin) );
	  
	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
  	  double corr_value = (double) old_value/eff;
  	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double relativeError = sqrt( (old_err*old_err) + ( hEff[a][fileNo]->GetBinError(ptbin)*hEff[a][fileNo]->GetBinError(ptbin)) );
  	  double corr_err = corr_value*relativeError;
  	  hChgUEpt_te[a][p]->SetBinContent( xbin, ybin, corr_value);
  	  hChgUEpt_te[a][p]->SetBinError( xbin, ybin, corr_err );
  	}


  	// sys1: trackeffic+0.05
  	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
  	  int xbin = ix;
	  
  	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
  	  // if ( pt > 3.0 ) { pt = 3.0; }
  	  double ptbin = hEff[a][fileNo]->FindBin(pt);
  	  // double eff = hEff[a][fileNo]->GetBinContent(ptbin) + 0.05;
 	  double eff = efficFit[a][fileNo]->Eval( hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin) ) + 0.05;
	  
 	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
  	  double corr_value = (double) old_value/eff;
  	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double relativeError = sqrt( (old_err*old_err) + ( hEff[a][fileNo]->GetBinError(ptbin)*hEff[a][fileNo]->GetBinError(ptbin)) );
  	  double corr_err = corr_value*relativeError;
  	  hChgUEpt_te_sys1[a][p]->SetBinContent( xbin, ybin, corr_value);
  	  hChgUEpt_te_sys1[a][p]->SetBinError( xbin, ybin, corr_err );
  	}


  	// sys2: trackeffic-0.05
  	for(int ix = 1; ix <= hChgUEpt[a][p]->GetNbinsX(); ++ix){
  	  int xbin = ix;
	  
  	  double pt = hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin);
  	  // if ( pt > 3.0 ) { pt = 3.0; }
  	  double ptbin = hEff[a][fileNo]->FindBin(pt);
  	  // double eff = hEff[a][fileNo]->GetBinContent(ptbin) - 0.05;
	  double eff = efficFit[a][fileNo]->Eval( hChgUEpt[a][p]->GetXaxis()->GetBinCenter(xbin) ) - 0.05;
  	  double old_value = hChgUEpt[a][p]->GetBinContent(xbin,ybin);
  	  double corr_value = (double) old_value/eff;
  	  double old_err = hChgUEpt[a][p]->GetBinError(xbin,ybin);
	  double relativeError = sqrt( (old_err*old_err) + ( hEff[a][fileNo]->GetBinError(ptbin)*hEff[a][fileNo]->GetBinError(ptbin)) );
  	  double corr_err = corr_value*relativeError;
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
  	nChg_te_sys2[e][a][p] = hChgUEpt_te_sys2[a][p]->Integral();

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

  	cout<<te_meanChgPt[e][a][p]<<"\t"; 
      }
      cout<<endl;
    }
    cout<<endl;
  }
  cout<<endl;
  
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


  double binEdge[nEtaBins+1] = { -0.5, 0.5, 1.5, 2.5 };

  //  double ptBinEdge[nPtBins+1] = { 9.5, 14.5, 19.5, 29.5 };
  double ptBinEdge[nPtBins+1] = { 10, 15, 20, 20 };
  double yEdge[2] = {0.5,1.8};
  int startBin[nPtBins] = { 11, 30, 59 };
  

  TH2D *hscale0 = new TH2D("hscale0",";p_{T,lead}^{reco} [GeV/#it{c}];#LT #frac{d#it{N}_{ch}}{d#eta d#phi} #GT", 10,10,30,10,0.6,1.9);
  hscale0->GetYaxis()->SetTitleOffset(1);
  hscale0->GetXaxis()->SetTitleOffset(1.0);
  hscale0->GetYaxis()->SetTitleSize(0.05);
  hscale0->GetXaxis()->SetLabelSize(0.04);
  hscale0->GetYaxis()->SetLabelSize(0.04);
  hscale0->GetXaxis()->SetTitleOffset(1.0);
  hscale0->GetXaxis()->SetTitleSize(0.05);

  hscale0->GetXaxis()->SetNdivisions(6);
  
  hscale0->SetName("");
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

  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hCHGdNdetadphi_te_" + EAbinName[a] + etaBinName[e];
      title = UEetaBinString[e];
      hCHGdNdetadphi_te[a][e] = new TH1D( name, title , 79,9.25,29.75 );

      hCHGdNdetadphi_te[a][e]->SetStats(0);
      hCHGdNdetadphi_te[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi_te[a][e]->SetMarkerSize( 2 );

      leg0->AddEntry(hCHGdNdetadphi_te[a][e],title,"lpf");
      
      name = "hCHGdNdetadphi_te_sys_" + EAbinName[a] + etaBinName[e];
      hCHGdNdetadphi_te_sys[a][e] = new TH1D( name, "" , 79,9.25,29.75 );
      
      hCHGdNdetadphi_te_sys[a][e]->SetStats(0);
      hCHGdNdetadphi_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      hCHGdNdetadphi_te_sys[a][e]->SetMarkerSize( 0 );

      if (a==0) {
  	hCHGdNdetadphi_te_sys[a][e]->SetFillStyle(3002);
  	hCHGdNdetadphi_te_sys[a][e]->SetFillColor( etaColor[e]);
      }
      else { hCHGdNdetadphi_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 ); }

      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
  	int binno = startBin[p] + (2*e);
	
  	double value = (double) ( nChg_te[e][a][p] / area[e] );
  	hCHGdNdetadphi_te[a][e]->SetBinContent( binno, value );
  	hCHGdNdetadphi_te[a][e]->SetBinError( binno, nChg_te_err[e][a][p] );

  	for (int i=-3; i<4; ++i) {
  	  hCHGdNdetadphi_te_sys[a][e]->SetBinContent( binno+i, value );
  	  hCHGdNdetadphi_te_sys[a][e]->SetBinError( binno+i, value*nChg_te_sys[e][a][p] );
  	}
	
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
  //tex4->Draw();

  //c0->SetGridx();
  
  leg0->Draw();
  name = "plots/UE/"+dirName+"/CHGdNdetadphi_te_pt_systematics.pdf";
  c0->SaveAs(name,"PDF");
  //c0->Close();

 
  TCanvas *c1 = new TCanvas;
  c1->SetMargin(0.125,0.05,0.125,0.05);

  double yEdge1[2] = { 0.55, 0.85 };
  TH2D *hscale3 = new TH2D("hscale3",";p_{T,lead}^{reco} [GeV/#it{c}];#LT p_{T}^{ch} #GT [GeV/#it{c}]", 10,10,30,10,0.45,0.95);
  hscale3->GetYaxis()->SetTitleOffset(1.0);
  hscale3->GetYaxis()->SetTitleSize(0.05);
  hscale3->GetXaxis()->SetTitleOffset(1.0);
  hscale3->GetXaxis()->SetTitleSize(0.05);
  hscale3->GetXaxis()->SetLabelSize(0.04);
  hscale3->GetYaxis()->SetLabelSize(0.04);
  hscale3->SetName("");

  hscale3->GetXaxis()->SetNdivisions(6);

  TH1D *hChg_pt_te[nEAbins][nEtaBins];
  TH1D *hChg_pt_te_sys[nEAbins][nEtaBins];
  hscale3->SetName(" ");
  hscale3->SetLineWidth(0);
  hscale3->SetStats(0);
  hscale3->Draw();


  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hChg_MeanPt_te_" + EAbinName[a] + etaBinName[e];
      title = EAbinString[a] + "   "  + eastmidwest[e] + "   ";
      hChg_pt_te[a][e] = new TH1D( name, title , 79,9.25,29.75 );

      hChg_pt_te[a][e]->SetStats(0);
      hChg_pt_te[a][e]->SetLineColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt_te[a][e]->SetMarkerSize( 2 );
      
      name = "hChg_MeanPt_te_sys_" + EAbinName[a] + etaBinName[e];
      hChg_pt_te_sys[a][e] = new TH1D( name, "" , 79,9.25,29.75 );

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
	
      for ( int p=0; p<nPtBins; ++p) {
	
  	int binno = startBin[p] + (2*e);

  	double value = te_meanChgPt[e][a][p];
  	hChg_pt_te[a][e]->SetBinContent( binno, value );
  	hChg_pt_te[a][e]->SetBinError( binno, te_meanChgPt_err[e][a][p] );

  	for (int i=-2; i<3; ++i) {
  	  hChg_pt_te_sys[a][e]->SetBinContent( binno+i, value );
  	  hChg_pt_te_sys[a][e]->SetBinError( binno+i, value*nChg_te_sys[e][a][p] );
  	}
		
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
  //tex4->Draw();

  name = "plots/UE/"+dirName+"/Chg_MeanPt_te_pt_systematics.pdf";
  c1->SaveAs(name,"PDF");
  //c1->SetGridx();


  TFile *outFile = new TFile("test.root","RECREATE");
  
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {

      hChgUEpt[a][p]->Write();
      hChgUEpt_te[a][p]->Write();
      // hNeuUEpt[a][p]->Write();
	
    }
  }
  // c0->Write();
  // c1->Write();

  for (int i=0; i<eta_bins; ++i) {
    for (int j=0; j<nEAbins; ++j){
      efficFit[i][j]->Write();
      hEff[i][j]->Write();
    }
  }

  
  outFile->Write();
  
}
