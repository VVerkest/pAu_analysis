// Veronica Verkest
// November 11, 2019

void UEjetPlot(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };

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
  const double BBCEsumLo[nEA] = { 4107, 11503, 28537 };   // LO: 4107-11503;  HI: 28537+
  const double BBCEsumHi[nEA] = { 11503, 28537, 64000 };
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, rho;
  
  TString fileName = "out/UE/pAuHTjetUE.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH2D *hBGchg = (TH2D*) inFile->Get("hChgBgEtaPhi");
  TH2D *hBGneu = (TH2D*) inFile->Get("hNeuBgEtaPhi");

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, eastRho, midRho, westRho;

  jetTree->SetBranchAddress( "RunID", &RunID );      	       	jetTree->SetBranchAddress( "EventID", &EventID );					jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );       	jetTree->SetBranchAddress( "nGlobal", &nGlobal );					jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );		jetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );	       		jetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );	       	jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );		jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );	jetTree->Branch( "leadArea", &leadArea );

  int nEntries = jetTree->GetEntries();

  TH1D *hRho = new TH1D("hRho","Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hRho_HI = new TH1D("hRho_HI","High EA Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hRho_LO = new TH1D("hRho_LO","Low EA Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hLeadPhi = new TH1D("hLeadPhi","Lead Jet #phi;#phi_{lead}",12,0,2*pi);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 60,0,15, 200,0,200 );
  TH3D *hLeadJetPtRhoEta = new TH3D( "hLeadJetPtRhoEta", "Lead Jet p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );  
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 140,0,70000 );
  TH1D *hBBCsumE = new TH1D("hBBCsumE","BBC ADC East Sum;BBCE Sum", 140,0,70000 );
  TH1D *hBBCsumE_integral = new TH1D("integral","BBC ADC East Sum Integral;Int(BBCE Sum)", 280,0,70000 );
  TH1D *hleadEta_LoEA = new TH1D( "hleadEta_LoEA", "Low Event Activity", 40,-1,1 );
  TH1D *hleadEta_HiEA = new TH1D( "hleadEta_HiEA", "High Event Activity", 40,-1,1 );

  TH2D *hLeadPtVsPtCorrected = new TH2D("hLeadPtVsPtCorrected","Lead Jet p_{T} vs. Background-Subtracted p_{T} ;p_{T} (GeV);p_{T}^{corrected} (GeV)", 80,10,30, 120,0,30 );
  TH1D *hPtCorrectedRatio = new TH1D("hPtCorrectedRatio","Lead Jet:   p_{T}^{corrected} / p_{T};p_{T}^{corrected}/p_{T}",10,0.9,1.0);
  
  TH1D *hLeadEta[nPtBins];
  TH1D *hBBCEastSum[nEtaBins];
  TH1D *hEAdist[nPtBins][nEtaBins];

  TH1D *hRhoDist[nEA][nEtaBins][nPtBins][nEtaBins];

  for ( int ea=0; ea<nEA; ++ea ) {
    for ( int bge=0; bge<nEtaBins; ++bge ) {
      for ( int p=0; p<nPtBins; ++p ) {
	for ( int je=0; je<nEtaBins; ++je ) {
	  name = "hRho_" + EAstring[ea] + jetEtaBinName[je] + ptBinName[p] + etaBinName[bge];
	  title = "Underlying Event" + EAstring[ea] + jetEtaBinName[je] + ptBinName[p] + etaBinName[bge] +";#rho (GeV)";
	  hRhoDist[ea][bge][p][je]= new TH1D(name,title,120,0,30);
	}
      }
    }
  }
    
  for ( int e=0; e<nEtaBins; ++e ) {

    name = "hBBCEastSum" + etaBinName[e];
    title = "BBC ADC East Sum:  " + etaBinString[e] + ";BBC East Sum";
    hBBCEastSum[e] = new TH1D( name, title, 18,0,70000 );
    hBBCEastSum[e]->SetMarkerStyle( etaMarker[e] );
    hBBCEastSum[e]->SetMarkerColor( etaColor[e] );
    hBBCEastSum[e]->SetLineColor( etaColor[e] );

    for ( int p=0; p<nPtBins; ++p ) {
      name = "hEAdist" + ptBinName[p] + etaBinName[e];
      title = "BBC ADC East Sum:  " + ptBinString[p] + "     " + etaBinString[e] + ";BBC East Sum";
      hEAdist[p][e] = new TH1D( name, title, 18,0,70000 );
      hEAdist[p][e]->SetLineColor( etaColor[e] );
      hEAdist[p][e]->SetMarkerColor( etaColor[e] );
      hEAdist[p][e]->SetMarkerStyle( etaMarker[e] );
    }
  }
  
  for ( int p=0; p<nPtBins; ++p ) {
    name = "hLeadEta" + ptBinName[p];
    title = "Lead Jet #eta:  " + ptBinString[p] + ";#eta_{lead}";
    hLeadEta[p] = new TH1D( name, title, 40, -1.0, 1.0 );
    hLeadEta[p]->SetLineColor( ptColor[p] );
    hLeadEta[p]->SetMarkerColor( ptColor[p] );
  }

  double rhoByEta[nEtaBins];
  double BBCEintegral;
  
  for ( int i=0; i<nEntries; ++i ) {

    jetTree->GetEntry(i);

    rhoByEta[0] = chgEastRho + neuEastRho;
    rhoByEta[1] = chgMidRho + neuMidRho;
    rhoByEta[2] = chgWestRho + neuWestRho;
    
    chgRho = ( chgEastRho + chgMidRho + chgWestRho )/3;
    neuRho = ( neuEastRho + neuMidRho + neuWestRho )/3;
    rho = chgRho + neuRho;

    hRho->Fill(rho);
    hLeadPhi->Fill(leadPhi);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadJetPtRhoEta->Fill(leadPt,rho,leadEta);
    hLeadPtEtaPhi->Fill(leadPt,leadEta,leadPhi);
    hPt_UE_BBCsumE->Fill(leadPt,rho,BbcAdcEastSum);
    hBBCsumE->Fill(BbcAdcEastSum);
    BBCEintegral = hBBCsumE->Integral(0,i);
    hBBCsumE_integral->Fill( BBCEintegral );
  
    if ( BbcAdcEastSum>4107 && BbcAdcEastSum<11503 ) { // LO: 4107-11503;  HI: 28537+
      hleadEta_LoEA->Fill(leadEta);
    }

    if ( BbcAdcEastSum>28537 ) {
      hleadEta_HiEA->Fill(leadEta);
    }
	
    pval = 99;    jeval = 99;    eaval = 99;
    
    for ( int ea=0; ea<3; ++ea ) {
      if ( BbcAdcEastSum > BBCEsumLo[ea]  &&  BbcAdcEastSum < BBCEsumHi[ea] ) { eaval = ea; }
    }    
    for ( int p=0; p<3; ++p ) {
      if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
    }
    for ( int je=0; je<3; ++je ) {
      if ( leadEta >= etaLo[je]  &&  leadEta <= etaHi[je] ) { jeval = je; }
    }
    if ( eaval==99 ) { continue; }
    if ( pval==99 || jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<pval<<endl<<jeval<<endl<<bgeval<<endl<<leadEta<<endl<<endl; }
    
    for ( int bge=0; bge<nEtaBins; ++bge ) {
      hRhoDist[eaval][bge][pval][jeval]->Fill( rhoByEta[bge] );
    }
  
    hBBCEastSum[jeval]->Fill( BbcAdcEastSum );
    hLeadEta[pval]->Fill( leadEta );
    hEAdist[pval][jeval]->Fill( BbcAdcEastSum );
    
  }


  double avgRho[nEA][nEtaBins][nPtBins][nEtaBins];
  double stdevRho[nEA][nEtaBins][nPtBins][nEtaBins];
  
  for ( int ea=0; ea<nEA; ++ea ) {
    for ( int bge=0; bge<nEtaBins; ++bge ) {
      for ( int p=0; p<nPtBins; ++p ) {
  	for ( int je=0; je<nEtaBins; ++je ) {
  	  avgRho[ea][bge][p][je] = hRhoDist[ea][bge][p][je]->GetMean(1);
  	  stdevRho[ea][bge][p][je] = hRhoDist[ea][bge][p][je]->GetMeanError(1);
  	}
      }
    }
  }

  for ( int i=0; i<nEntries; ++i ) {
    jetTree->GetEntry(i);
    
    for ( int ea=0; ea<nEA; ++ea ) {
      if ( BbcAdcEastSum >= BBCEsumLo[ea]  &&  BbcAdcEastSum <= BBCEsumHi[ea] ) { eaval = ea; }
    }    
    for ( int p=0; p<nPtBins; ++p ) {
      if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
    }
    for ( int je=0; je<nEtaBins; ++je ) {
      if ( leadEta >= etaLo[je]  &&  leadEta <= etaHi[je] ) { jeval = je; }
    }
    if ( pval==99 || jeval==99 || bgeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl; }
    if ( eaval==99 ) { continue; }

    double bsPt = leadPt - pi*(0.4)*(0.4)*avgRho[eaval][jeval][pval][jeval];
    hLeadPtVsPtCorrected->Fill( leadPt, bsPt );
    hPtCorrectedRatio->Fill( bsPt/leadPt );
  }


  
  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS 0

  hPtCorrectedRatio->Scale(1./hPtCorrectedRatio->GetEntries());
  hPtCorrectedRatio->Draw();
  c0->SaveAs( "plots/UE/correctedPtRatio.pdf" , "PDF" );

  hLeadPhi->Scale(1./hLeadPhi->Integral("WIDTH"));
  hLeadPhi->Draw();
  c0->SaveAs( "plots/UE/leadPhi.pdf" , "PDF" );

  c0->SetLogz();

  hLeadPtVsPtCorrected->Scale(1./hLeadPtVsPtCorrected->Integral("WIDTH"));
  hLeadPtVsPtCorrected->Draw("COLZ");
  c0->SaveAs( "plots/UE/leadPtVsPtCorrected.pdf" , "PDF" );
  
  hBGchg->Scale(1./hBGchg->Integral("WIDTH"));
  hBGchg->Draw("COLZ");
  c0->SaveAs( "plots/UE/chgBgEtaPhi.pdf" , "PDF" );
  
  hBGneu->Scale(1./hBGneu->Integral("WIDTH"));
  hBGneu->Draw("COLZ");
  c0->SaveAs( "plots/UE/neuBgEtaPhi.pdf" , "PDF" );

  hTowersVsRho->Scale(1./hTowersVsRho->Integral("WIDTH"));
  hTowersVsRho->GetZaxis()->SetRangeUser(0.000001,1);
  hTowersVsRho->Draw("COLZ");
  c0->SaveAs( "plots/UE/towersVsRho.pdf" , "PDF" );

  c0->SetLogy();

  hRho->Scale(1./hRho->Integral("WIDTH"));
  hRho->Draw();
  c0->SaveAs( "plots/UE/rho.pdf" , "PDF" );


  TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 1
  
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);
  leg0->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
  leg0->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");
  
  TH2D *sBBCbyEta = new TH2D("sBBCbyEta", "BBC ADC East Sum by Lead Jet #eta;BBC East Sum", 35,0,70000, 10,0,0.15);
  sBBCbyEta->SetStats(0);
  sBBCbyEta->Draw();
  for ( int e=0; e<nEtaBins; ++e ) {
    hBBCEastSum[e]->SetStats(0);
    hBBCEastSum[e]->Scale(1./hBBCEastSum[e]->Integral());
    hBBCEastSum[e]->Draw("SAME");
    avg = "";
    avg += hBBCEastSum[e]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,5);
    name = "hBBCEastSum" + etaBinName[e];
    title = etaBinString[e];
    leg0->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg0->AddEntry((TObject*)0,avg, "");

  }
  leg0->Draw();
  c1->SaveAs( "plots/UE/BBCEastSum_by_eta.pdf" , "PDF" );
  c1->SetLogy(0);


  TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  leg1->SetNColumns(2);
  leg1->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg1->AddEntry((TObject*)0,"#bf{<#eta>}", "");
  
  TH2D *sLeadEtaByPt = new TH2D("sLeadEtaByPt", "Lead Jet #eta by p_{T};#eta_{lead}", 40,-1,1, 20,0.5,1);
  sLeadEtaByPt->SetStats(0);
  sLeadEtaByPt->Draw();

  for ( int p=0; p<nPtBins; ++p ) {
    hLeadEta[p]->SetStats(0);
    hLeadEta[p]->Scale(1./hLeadEta[p]->Integral("WIDTH"));
    hLeadEta[p]->SetMarkerStyle( 20 );
    hLeadEta[p]->Draw("SAME");
    avg = "";
    avg += hLeadEta[p]->GetMean(1);
    avg = avg(0,6);
    name = "hLeadEta" + ptBinName[p];
    title = ptBinString[p];
    leg1->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg1->AddEntry((TObject*)0,avg, "");

  }
  leg1->Draw();
  c1->SaveAs( "plots/UE/LeadEta_by_pt.pdf" , "PDF" );


  jetTree->Draw("leadPt:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRho2d","","COLZ");
  TH2D *hRho2d = (TH2D*)gDirectory->Get("hRho2d");
  c1->SetLogy();
  TLegend *leg2 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(3);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho[nPtBins];
  TH2D *sPtRho = new TH2D( "sPtRho", "Underlying Event by Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sPtRho->SetStats(0);
  sPtRho->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRho2d->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hPtRho[p] = (TH1D*) hRho2d->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hPtRho[p]->SetNameTitle( name, title );
    hPtRho[p]->SetLineColor( ptColor[p] );
    hPtRho[p]->SetMarkerColor( ptColor[p] );
    hPtRho[p]->SetMarkerStyle( 20 );
    hPtRho[p]->Scale(1./hPtRho[p]->GetEntries());
    hPtRho[p]->SetStats(0);
    hPtRho[p]->Draw("SAME");

    avg = "";
    avg+=hPtRho[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hPtRho[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg2->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg2->AddEntry((TObject*)0,avg, "");
    leg2->AddEntry((TObject*)0,sigma, "");
  }
  leg2->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt.pdf","PDF");





  // LO: 4107-11503
  jetTree->Draw("leadPt:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRho2d_LO","BbcAdcEastSum>4107 && BbcAdcEastSum<11503","COLZ");
  TH2D *hRho2d_LO = (TH2D*)gDirectory->Get("hRho2d_LO");
  c1->SetLogy();
  TLegend *leg3 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg3->SetBorderSize(1);   leg3->SetLineColor(1);   leg3->SetLineStyle(1);   leg3->SetLineWidth(1);   leg3->SetFillColor(0);   leg3->SetFillStyle(1001);
  leg3->SetNColumns(3);
  leg3->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg3->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg3->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho_LO[nPtBins];
  TH2D *sPtRho_LO = new TH2D( "sPtRho_LO", "Low EA: Underlying Event by Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sPtRho_LO->SetStats(0);
  sPtRho_LO->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRho2d_LO->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hPtRho_LO[p] = (TH1D*) hRho2d_LO->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hPtRho_LO[p]->SetNameTitle( name, title );
    hPtRho_LO[p]->SetLineColor( ptColor[p] );
    hPtRho_LO[p]->SetMarkerColor( ptColor[p] );
    hPtRho_LO[p]->SetMarkerStyle( 20 );
    hPtRho_LO[p]->Scale(1./hPtRho_LO[p]->GetEntries());
    hPtRho_LO[p]->SetStats(0);
    hPtRho_LO[p]->Draw("SAME");

    avg = "";
    avg+=hPtRho_LO[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hPtRho_LO[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg3->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg3->AddEntry((TObject*)0,avg, "");
    leg3->AddEntry((TObject*)0,sigma, "");
  }
  leg3->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt_LOEA.pdf","PDF");




  // HI: 28537+
  jetTree->Draw("leadPt:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRho2d_HI","BbcAdcEastSum>28537","COLZ");
  TH2D *hRho2d_HI = (TH2D*)gDirectory->Get("hRho2d_HI");
  c1->SetLogy();
  TLegend *leg4 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg4->SetBorderSize(1);   leg4->SetLineColor(1);   leg4->SetLineStyle(1);   leg4->SetLineWidth(1);   leg4->SetFillColor(0);   leg4->SetFillStyle(1001);
  leg4->SetNColumns(3);
  leg4->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg4->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg4->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho_HI[nPtBins];
  TH2D *sPtRho_HI = new TH2D( "sPtRho_HI", "High EA: Underlying Event by Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sPtRho_HI->SetStats(0);
  sPtRho_HI->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRho2d_HI->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hPtRho_HI[p] = (TH1D*) hRho2d_HI->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hPtRho_HI[p]->SetNameTitle( name, title );
    hPtRho_HI[p]->SetLineColor( ptColor[p] );
    hPtRho_HI[p]->SetMarkerColor( ptColor[p] );
    hPtRho_HI[p]->SetMarkerStyle( 20 );
    hPtRho_HI[p]->Scale(1./hPtRho_HI[p]->GetEntries());
    hPtRho_HI[p]->SetStats(0);
    hPtRho_HI[p]->Draw("SAME");

    avg = "";
    avg+=hPtRho_HI[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hPtRho_HI[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg4->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg4->AddEntry((TObject*)0,avg, "");
    leg4->AddEntry((TObject*)0,sigma, "");
  }
  leg4->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt_HIEA.pdf","PDF");







  
  //  corrected lead pt
  
  jetTree->Draw("leadPtCorrected:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRhoCorr2d","","COLZ");
  TH2D *hRhoCorr2d = (TH2D*)gDirectory->Get("hRhoCorr2d");
  c1->SetLogy();
  TLegend *leg12 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg12->SetBorderSize(1);   leg12->SetLineColor(1);   leg12->SetLineStyle(1);   leg12->SetLineWidth(1);   leg12->SetFillColor(0);   leg12->SetFillStyle(1001);
  leg12->SetNColumns(3);
  leg12->AddEntry((TObject*)0,"#bf{p_{T}^{corrected}}", "");
  leg12->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg12->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hCorrPtRho[nPtBins];
  TH2D *sCorrPtRho = new TH2D( "sCorrPtRho", "Underlying Event by Corrected Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sCorrPtRho->SetStats(0);
  sCorrPtRho->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho[p] = (TH1D*) hRhoCorr2d->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hCorrPtRho[p]->SetNameTitle( name, title );
    hCorrPtRho[p]->SetLineColor( ptColor[p] );
    hCorrPtRho[p]->SetMarkerColor( ptColor[p] );
    hCorrPtRho[p]->SetMarkerStyle( 20 );
    hCorrPtRho[p]->Scale(1./hCorrPtRho[p]->GetEntries());
    hCorrPtRho[p]->SetStats(0);
    hCorrPtRho[p]->Draw("SAME");

    avg = "";
    avg+=hCorrPtRho[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hCorrPtRho[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg12->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg12->AddEntry((TObject*)0,avg, "");
    leg12->AddEntry((TObject*)0,sigma, "");
  }
  leg12->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt_correctedPt.pdf","PDF");





  // LO: 4107-11503
  jetTree->Draw("leadPtCorrected:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRhoCorr2d_LO","BbcAdcEastSum>4107 && BbcAdcEastSum<11503","COLZ");
  TH2D *hRhoCorr2d_LO = (TH2D*)gDirectory->Get("hRhoCorr2d_LO");
  c1->SetLogy();
  TLegend *leg13 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg13->SetBorderSize(1);   leg13->SetLineColor(1);   leg13->SetLineStyle(1);   leg13->SetLineWidth(1);   leg13->SetFillColor(0);   leg13->SetFillStyle(1001);
  leg13->SetNColumns(3);
  leg13->AddEntry((TObject*)0,"#bf{p_{T}^{corrected}}", "");
  leg13->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg13->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hCorrPtRho_LO[nPtBins];
  TH2D *sCorrPtRho_LO = new TH2D( "sCorrPtRho_LO", "Low EA: Underlying Event by Corrected Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sCorrPtRho_LO->SetStats(0);
  sCorrPtRho_LO->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d_LO->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho_LO[p] = (TH1D*) hRhoCorr2d_LO->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hCorrPtRho_LO[p]->SetNameTitle( name, title );
    hCorrPtRho_LO[p]->SetLineColor( ptColor[p] );
    hCorrPtRho_LO[p]->SetMarkerColor( ptColor[p] );
    hCorrPtRho_LO[p]->SetMarkerStyle( 20 );
    hCorrPtRho_LO[p]->Scale(1./hCorrPtRho_LO[p]->GetEntries());
    hCorrPtRho_LO[p]->SetStats(0);
    hCorrPtRho_LO[p]->Draw("SAME");

    avg = "";
    avg+=hCorrPtRho_LO[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hCorrPtRho_LO[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg13->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg13->AddEntry((TObject*)0,avg, "");
    leg13->AddEntry((TObject*)0,sigma, "");
  }
  leg13->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt_LOEA_correctedPt.pdf","PDF");




  // HI: 28537+
  jetTree->Draw("leadPtCorrected:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRhoCorr2d_HI","BbcAdcEastSum>28537","COLZ");
  TH2D *hRhoCorr2d_HI = (TH2D*)gDirectory->Get("hRhoCorr2d_HI");
  c1->SetLogy();
  TLegend *leg14 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg14->SetBorderSize(1);   leg14->SetLineColor(1);   leg14->SetLineStyle(1);   leg14->SetLineWidth(1);   leg14->SetFillColor(0);   leg14->SetFillStyle(1001);
  leg14->SetNColumns(3);
  leg14->AddEntry((TObject*)0,"#bf{p_{T}^{corrected}}", "");
  leg14->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg14->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hCorrPtRho_HI[nPtBins];
  TH2D *sCorrPtRho_HI = new TH2D( "sCorrPtRho_HI", "High EA: Underlying Event by Corrected Lead Jet p_{T};#rho (GeV)", 20,0,8, 10,0.000001,1 );
  sCorrPtRho_HI->SetStats(0);
  sCorrPtRhoyyyyyy_HI->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d_HI->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho_HI[p] = (TH1D*) hRhoCorr2d_HI->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinString[p];
    hCorrPtRho_HI[p]->SetNameTitle( name, title );
    hCorrPtRho_HI[p]->SetLineColor( ptColor[p] );
    hCorrPtRho_HI[p]->SetMarkerColor( ptColor[p] );
    hCorrPtRho_HI[p]->SetMarkerStyle( 20 );
    hCorrPtRho_HI[p]->Scale(1./hCorrPtRho_HI[p]->GetEntries());
    hCorrPtRho_HI[p]->SetStats(0);
    hCorrPtRho_HI[p]->Draw("SAME");

    avg = "";
    avg+=hCorrPtRho_HI[p]->GetMean(1);
    avg = avg(0,5);
    sigma="";
    sigma+=hCorrPtRho_HI[p]->GetStdDev(1);
    sigma = sigma(0,5);
    leg14->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg14->AddEntry((TObject*)0,avg, "");
    leg14->AddEntry((TObject*)0,sigma, "");
  }
  leg14->Draw();
  c1->SaveAs("plots/UE/rhoByLeadPt_HIEA_correctedPt.pdf","PDF");





  c1->SetLogy();

  hBBCsumE->SetBinError(4107,100000);
  hBBCsumE->SetBinError(4107,100000000000);
  hBBCsumE->SetBinError(7802,100000000000);
  hBBCsumE->SetBinError(11503,100000000000);
  hBBCsumE->SetBinError(15399,100000000000);
  hBBCsumE->SetBinError(19525,100000000000);;
  hBBCsumE->SetBinError(23781,100000000000);
  hBBCsumE->SetBinError(28537,100000000000);
  hBBCsumE->SetBinError(34082,100000000000);
  hBBCsumE->SetBinError(41108,100000000000);
  hBBCsumE->SetBinError(63974,100000000000);
  hBBCsumE->Draw();

  TCanvas * c2 = new TCanvas( "c2" , "" ,700 ,500 );              // CANVAS 1
  
  TLegend *leg5[nPtBins];
  
  for ( int p=0; p<nPtBins; ++p ) {
    leg5[p] = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
    leg5[p]->SetBorderSize(1);   leg5[p]->SetLineColor(1);   leg5[p]->SetLineStyle(1);   leg5[p]->SetLineWidth(1);   leg5[p]->SetFillColor(0);   leg5[p]->SetFillStyle(1001);
    leg5[p]->SetNColumns(2);
    leg5[p]->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
    leg5[p]->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");

    title = "BBC ADC East Sum by Lead Jet #eta:    " + ptBinString[p] + ";BBC East Sum";
    sBBCbyEta->SetTitle( title );
    //TH2D *sBBCbyEta = new TH2D("sBBCbyEta", title, 35,0,70000, 10,0,0.15);
    sBBCbyEta->SetStats(0);
    sBBCbyEta->Draw();
    for ( int e=0; e<nEtaBins; ++e ) {
      hEAdist[p][e]->SetStats(0);
      hEAdist[p][e]->Scale(1./hEAdist[p][e]->Integral());
      hEAdist[p][e]->Draw("SAME");
      avg = "";
      avg += hEAdist[p][e]->GetMean(1);                                           // 1 denotes x-axis
      avg = avg(0,5);
      name = "hEAdist" + ptBinName[p] + etaBinName[e];
      title = etaBinString[e];
      leg5[p]->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
      leg5[p]->AddEntry((TObject*)0,avg, "");

    }
    leg5[p]->Draw();
    saveName = "plots/UE/BBCEastSum_by_eta" + ptBinName[p] +".pdf";
    c2->SaveAs( saveName , "PDF" );
  }
}
