// Veronica Verkest
// May 14, 2020

void rho_debug(){

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
  TString ptSelection[nPtBins] = { "leadPtCorrected>10.0 && leadPtCorrected<15.0", "leadPtCorrected>=15.0 && leadPtCorrected<=20.0", "leadPtCorrected>20.0 && leadPtCorrected<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 23, 22 };
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;
  double chgRho, neuRho;

  TString fileName[nEAbins] = { "out/UE/pAuHTjetUE_loEA.root", "out/UE/pAuHTjetUE_hiEA.root" };
  TString efficFileName[nEAbins] = { "src/trackeffic_loEA.root", "src/trackeffic_hiEA.root" };

  TFile* inFile[nEAbins];
  TFile* efficFile[nEAbins];
    
  for (int i=0; i<nEAbins; ++i) {
    efficFile[i] = new TFile( efficFileName[i], "READ" );
  }

  const int eta_bins = 10;
  TH1D *hEff[eta_bins][nEAbins];
  TF1 *hEff_fit1[eta_bins][nEAbins];
  TF1 *hEff_fit2[eta_bins][nEAbins];

  
  for (int i=0; i<eta_bins; ++i) {
    for (int j=0; j<nEAbins; ++j){    
    int binno = i+1;
    
    if (j==0 ) { name = "eff_s_bin_1_3_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
    else if (j==1) {name = "eff_s_bin_7_10_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
    
    hEff[i][j] = (TH1D*)efficFile[j]->Get(name);
    name += "_fit";
    hEff_fit1[i][j] = new TF1(name,"pol3",0,10);
    hEff[i][j]->Fit( hEff_fit1[i][j], "R" );
    name += "2";
    hEff_fit2[i][j] = new TF1(name,"pol0",10,30);
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

  TH3D *hChgUEpt[nEAbins]; TH3D *hNeuUEpt[nEAbins];  TH1D *hRho[nEAbins];  TH1D *hNchg[nEAbins];  TH1D *hNneu[nEAbins];
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

    hChgUEpt[i] = (TH3D*)inFile[i]->Get("hChgBgPtEta_leadPt");
    hNeuUEpt[i] = (TH3D*)inFile[i]->Get("hNeuBgPtEta_leadPt");

    name = "hRho_"; name += EAbinString[i];
    hRho[i] = new TH1D(name,"UE;<#rho> (GeV)", 30,0.0,15.0);

    name = "hNchg_"; name += EAbinString[i];
    hNchg[i] = new TH1D("hNchg","# charged particles in UE", 50,0.0,50.0);

    name = "hNneu_"; name += EAbinString[i];
    hNneu[i] = new TH1D("hNneu","# neutral particles in UE", 50,0.0,50.0);

  }


  TH2D *hChgUEpt2D[nEAbins];  TH2D *hNeuUEpt2D[nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    for (int i=0; i<nEntries[a]; ++i) {
      jetTree[a]->GetEntry(i);
      if ( !(leadPt[a]>ptLo[0] && leadPt[a]<ptHi[2]) ) { continue; }
      // if ( !(leadPt[a]>ptLo[1] && leadPt[a]<ptHi[1]) ) { continue; }
      else {

	
	hRho[a]->Fill(rho[a]);

	hNchg[a]->Fill(nBGpart_chg[a]);
	hNneu[a]->Fill(nBGpart_neu[a]);

	
	// if ( leadPt >= 10.0 && leadPt <= 30.0 ) {
	// for ( int p=0; p<3; ++p ) {
	//   if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
	// }
      }
      
    }
  }

  
  for (int a=0; a<nEAbins; ++a) {

    hChgUEpt[a]->GetZaxis()->SetRangeUser(10.0, 30.0);
    hChgUEpt2D[a] = (TH2D*)hChgUEpt[a]->Project3D("YX");
    hChgUEpt2D[a]->Scale(1./hRho[a]->GetEntries());
    name = "hChgUEpt2D_"; name += EAbinString[a]; name += "EA";
    hChgUEpt2D[a]->SetName( name );
    
    hNeuUEpt[a]->GetZaxis()->SetRangeUser(10.0, 30.0);
    hNeuUEpt2D[a] = (TH2D*)hNeuUEpt[a]->Project3D("YX");
    hNeuUEpt2D[a]->Scale(1./hRho[a]->GetEntries()); 
    name = "hNeuUEpt2D_"; name += EAbinString[a]; name += "EA";
    hNeuUEpt2D[a]->SetName( name );
 }
  
  TH1D *hChgPtDist[ybins][nEAbins];
  
  int ecolor[ybins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  TString ybinEdgeString[ybins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };
  
  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS 0

  c0->SetLogy();


  // HERE
  TH1D *hTreeChgRho[nEAbins];
  TH1D *hChgDist[nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    for (int i=0; i<ybins; ++i) {
      int binno = i+1;
      name = "hChgpTdist_10_15GeV_"; name += EAbinString[a]; name += "EA_etabin"; name += binno;
      hChgUEpt2D[a]->GetYaxis()->SetRange( binno, binno+1 );
      hChgUEpt2D[a]->SetLineColor( ecolor[i] );
      hChgUEpt2D[a]->SetMarkerColor( ecolor[i] );
      hChgPtDist[i][a] = (TH1D*) hChgUEpt2D[a]->ProjectionX( name, binno, binno+1 );
      title = "Charged #rho dist. ("; title += ybinEdgeString[i]; title +="<#eta<"; title += ybinEdgeString[binno]; title += ");#rho (GeV)";
      hChgPtDist[i][a]->SetTitle( title );
      hChgPtDist[i][a]->GetYaxis()->SetRangeUser(0.000001,0.3);
      hChgPtDist[i][a]->Draw("SAME");
    }
  }
  c0->Close();



  for (int j=0; j<nEAbins; ++j) {
    TCanvas *c1 = new TCanvas( "c1" , "" ,700 ,500 );
    c1->SetLogy();

    name = "hTreeChgRho_"; name += EAbinString[j];
    hTreeChgRho[j] = new TH1D(name,"",30,0,15);
    drawString = "((0.35*chgEastRho)+(0.3*chgMidRho)+(0.35*chgWestRho))>>"; drawString += name;
    jetTree[j]->Draw( drawString ,"leadPt>10.0 && leadPt<30.0");
    hTreeChgRho[j]->Scale(1./hTreeChgRho[j]->GetEntries());

    hChgDist[j] = (TH1D*)hChgPtDist[0][j]->DrawCopy();
    name = "hChgDist_"; name += EAbinString[j];
    title = EAbinString[j]; title += "EA: "; title += "Charged #rho (p_{T}) distribution;#rho or p_{T} (GeV)";
    hChgDist[j]->SetNameTitle( name , title );
    for (int i=1; i<ybins; ++i) { hChgDist[j]->Add( hChgPtDist[i][j] ); }
    hChgDist[j]->SetStats(0);


    TLegend *leg = new TLegend(0.4, 0.6, 0.89, 0.89,NULL,"brNDC");    // LEGEND 0
    leg->SetBorderSize(0); leg->SetLineColor(1); leg->SetLineStyle(1); leg->SetLineWidth(1); leg->SetFillColor(0); leg->SetFillStyle(1001);
    leg->SetNColumns(3);
    leg->AddEntry((TObject*)0,"name", "");
    leg->AddEntry((TObject*)0,"#bf{<#rho_{chg}>}", "");
    leg->AddEntry((TObject*)0,"#bf{<n_{chg}>}", "");

    hChgDist[j]->SetAxisRange(0.000005,3.0, "Y");
    hChgDist[j]->Draw();
    hTreeChgRho[j]->SetStats(0);
    hTreeChgRho[j]->SetAxisRange(0.000005,3.0, "Y");
    hTreeChgRho[j]->Draw("SAME");
    //c1->BuildLegend();
	
    title = "Avg. UE p_{T} dist";
    leg->AddEntry( hChgDist[j], title , "lpf" );                            // ADD TO LEGEND
    avg = "";
    avg += hChgDist[j]->GetMean(1)*(hChgDist[j]->Integral()/AREA);                       // 1 denotes x-axis
    avg = avg(0,5);
    leg->AddEntry((TObject*)0,avg, "");
    avg = "";
    avg += hChgDist[j]->Integral();//hRho[j]->Integral();                       // 1 denotes x-axis
    //cout<<hChgDist[j]->Integral()<<endl;
    avg = avg(0,5);
    leg->AddEntry((TObject*)0,avg, "");

    title = "Event-by-event #rho dist.";
    leg->AddEntry( hTreeChgRho[j], title, "lpf" );                            // ADD TO LEGEND
    avg = "";
    avg += hTreeChgRho[j]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,5);
    leg->AddEntry((TObject*)0,avg, "");
    avg = "";
    avg += hNchg[j]->GetMean();                       // 1 denotes x-axis
    avg = avg(0,5);
    leg->AddEntry((TObject*)0,avg, "");
    
    leg->Draw();
    saveName = "plots/UE/rho_debug/UEptMeanCompare_"; saveName += EAbinString[j]; saveName += "EA.pdf";
    c1->SaveAs( saveName ,"PDF");
  
  }
  // c1->Close();


  TCanvas *c2 = new TCanvas( "c2" , "" ,700 ,500 );
  c2->SetLogy();
  
  for (int a=0; a<nEAbins; ++a) {
    for (int i=0; i<ybins; ++i) {
      int binno = i+1;
      name = "hChgpTdist_10_15GeV_"; name += EAbinString[a]; name += "EA__etabin"; name += binno;
      hChgUEpt2D[a]->GetYaxis()->SetRange( binno, binno+1 );
      // hChgPtDist[i] = (TH1D*)
      hChgUEpt2D[a]->SetLineColor( ecolor[i] );
      hChgUEpt2D[a]->SetMarkerColor( ecolor[i] );
      hChgPtDist[i][a] = (TH1D*) hChgUEpt2D[a]->ProjectionX( name, binno, binno+1 );
      if (i==0) { title = "";  title = EAbinString[a]; title += "EA Charged #rho dist.";  hChgPtDist[i][a]->SetTitle( title ); }
      else { title = "";  title = ybinEdgeString[i]; title += "<UE #eta<"; title += ybinEdgeString[binno];  hChgPtDist[i][a]->SetTitle( title ); }
      hChgPtDist[i][a]->SetStats(0);
      hChgPtDist[i][a]->SetAxisRange(0.000001,0.3, "Y");
      if (i==0) { hChgPtDist[i][a]->Draw(); }
      else { hChgPtDist[i][a]->Draw("SAME"); }
    }

    title = "";  title += ybinEdgeString[0]; title += "<UE #eta<"; title += ybinEdgeString[1];  hChgPtDist[0][a]->SetTitle( title );
    hChgPtDist[0][a]->SetTitle( title );
    TLegend *leg0 = (TLegend*) c2->BuildLegend(0.4,0.45,0.75,0.8);
    leg0->SetFillColorAlpha(0,0);  leg0->SetLineWidth(0);
    title = "";  title = EAbinString[a]; title += " EA Charged UE p_{T} dist.";  hChgPtDist[0][a]->SetTitle( title );
    hChgPtDist[0][a]->SetTitle( title );
    name = "plots/UE/rho_debug/UEptDist_"; name+= EAbinString[a]; name += "EA.pdf";
    c2->SaveAs( name ,"PDF");

  }


  // //hEff[0]->Draw();
  // for (int i=0; i<ybins; ++i) {
  //   int binno = i+1;

  //   hChgPtDist[i]->Divide( hEff_fit1[i] );
  //   hChgPtDist[i]->Divide( hEff_fit2[i] );
  //   //hChgPtDist[i]->Draw("SAME");
  //   // hEff_fit[i]->Draw("SAME");
  //   // new TCanvas;
  //   // hChgPtDist[i]->Draw();
  //   hChgPtDist[i]->GetXaxis()->SetRangeUser(0.01,30);
  // }

  // TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 0
  // c1->SetLogy();
  // for (int i=0; i<ybins; ++i) { hChgPtDist[i]->Draw("SAME"); }
  // c1->BuildLegend();
  // c1->SaveAs("plots/UE/rho_debug/UEptDist_efficcorrected_loEA.pdf","PDF");

  // TH1D* hChgUEdist = (TH1D*)hChgPtDist[0]->Clone();
  // hChgUEdist->SetName("hChgUEdist");
  
  // for (int i=1; i<ybins; ++i) {
  //   hChgUEdist->Add( hChgPtDist[i] );
  // }

  // hChgUEdist->Draw();

}
