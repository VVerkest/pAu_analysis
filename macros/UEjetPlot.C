// Veronica Verkest
// November 11, 2019

void UEjetPlot(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };
  const TString ptBinStringRaw[nPtBins] = { "10 < p_{T,lead}^{raw} < 15", "15 < p_{T,lead}^{raw} < 20",  "20 < p_{T,lead}^{raw} < 30" };
  const TString ptCorrectedBinString[nPtBins] = { "10 < p_{T}^{corrected} < 15", "15 < p_{T}^{corrected} < 20",  "20 < p_{T}^{corrected} < 30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };
  const int ptMarker[nPtBins] = { 20, 21, 29 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6 < #eta_{jet}^{lead} < -0.3", "-0.3 < #eta_{jet}^{lead} < 0.3", "0.3 < #eta_{jet}^{lead} < 0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

  const int nEAbins = 2;
  TString EAbinString[nEAbins] = { "Low EA", "High EA" };

  const int nEA = 3;
  const TString EAstring[nEA] = { "low", "mid", "high" };
  const double BBCEsumLo[nEA] = { 3559.12, 10126.1, 26718.1 };   // LO: 3559.12-10126.1;  HI: 26718.1+
  const double BBCEsumHi[nEA] = { 10126.1, 26718.1, 64000.0 };  // 0, 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000
  // { 0, 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000 }

  TString textString[5] = {"p+Au #sqrt{#it{s}_{NN}} = 200 GeV","E^{trig}_{T} > 5.4 GeV","anti-k_{T} R=0.4 jets","#left|#eta^{jets}#right| < 0.6",
			   "10 < p_{T,lead}^{reco} < 30 [GeV/#it{c}]"};
  TLatex *tex[5];

  for (int i=0; i<5; ++i) {
    tex[i] = new TLatex(0.2,0.85-0.05*i,textString[i]);
    tex[i]->SetTextFont(63);
    tex[i]->SetTextSize(16);
    tex[i]->SetTextColor(kBlack);
    tex[i]->SetLineWidth(1);
    tex[i]->SetNDC();
  }

  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, rho;
  
  //TString fileName = "out/UE/pAuHTjetUE.root";
  TString fileName = "out/UE/pAuHTjetUE_10cmVzCut_diffEta.root";
  //TString fileName = "out/UE/pAuHTjetUE_inLead.root";
  TFile* inFile = new TFile( fileName, "READ" );
  string plotdir = "plots/UE/10cmVzCut_diffEta/";

  // TH3D *hBGchg3D = (TH3D*) inFile->Get("hChgBgPtEtaPhi");
  // TH3D *hBGneu3D = (TH3D*) inFile->Get("hNeuBgPtEtaPhi");

  // TH2D *hBGchg = (TH2D*)hBGchg3D->Project3D("ZY");
  // TH2D *hBGneu = (TH2D*)hBGneu3D->Project3D("ZY");
  // hBGchg->SetName("hBGchg");
  // hBGneu->SetName("hBGneu");
  
  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nBGpart_chg, nBGpart_neu;
  double Vz, BbcAdcSumEast, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, eastRho, midRho, westRho, leadPtCorrected, chgEastRho_te, chgMidRho_te, chgWestRho_te, rho_te;

  jetTree->SetBranchAddress( "RunID", &RunID );
  jetTree->SetBranchAddress( "EventID", &EventID );
  jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );
  jetTree->SetBranchAddress( "nGlobal", &nGlobal );
  jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );
  jetTree->SetBranchAddress( "gRefMult", &gRefMult );
  jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );
  jetTree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );
  jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );
  jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );
  jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );
  jetTree->SetBranchAddress( "chgEastRho_te", &chgEastRho_te );
  jetTree->SetBranchAddress( "chgMidRho_te", &chgMidRho_te );
  jetTree->SetBranchAddress( "chgWestRho_te", &chgWestRho_te );
  jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );
  jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );
  jetTree->SetBranchAddress( "leadArea", &leadArea );
  jetTree->SetBranchAddress( "leadPtCorrected", &leadPtCorrected );
  jetTree->SetBranchAddress( "nBGpart_chg", &nBGpart_chg );
  jetTree->SetBranchAddress( "nBGpart_neu", &nBGpart_neu );
  jetTree->SetBranchAddress( "rho", &rho );
  jetTree->SetBranchAddress( "rho_te", &rho_te );

  int nEntries = jetTree->GetEntries();

  TH1D *hRho = new TH1D("hRho","Underlying Event;#rho [GeV]",60,0,15);
  TH1D *hRho_HI = new TH1D("hRho_HI","High EA Underlying Event;#rho [GeV]",120,0,40);
  TH1D *hRho_LO = new TH1D("hRho_LO","Low EA Underlying Event;#rho [GeV]",120,0,30);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho [GeV];# of Towers", 100,0,25, 200,0,200 );
  TH3D *hLeadJetPtRhoEta = new TH3D( "hLeadJetPtRhoEta", "Lead Jet p_{T}, #rho, #eta;Jet p_{T} [GeV];#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );  
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} [GeV];#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} [GeV];Underlying Event [GeV];BBC ADC East Sum", 500,0,125, 50,0,25, 140,0,70000 );
  TH1D *hBBCsumE = new TH1D("hBBCsumE","BBC ADC East Sum;BBCE Sum", 140,0,70000 );
  TH1D *hBBCsumE_integral = new TH1D("integral","BBC ADC East Sum Integral;Int(BBCE Sum)", 280,0,70000 );
  TH1D *hleadEta_LoEA = new TH1D( "hleadEta_LoEA", "Low Event Activity", 40,-1,1 );
  TH1D *hleadEta_HiEA = new TH1D( "hleadEta_HiEA", "High Event Activity", 40,-1,1 );

  TH2D *hLeadPtVsPtCorrected = new TH2D("hLeadPtVsPtCorrected","Lead Jet p_{T} vs. Background-Subtracted p_{T} ;p_{T} [GeV];p_{T}^{corrected} [GeV]", 80,10,30, 120,0,30 );
  
  TH1D *hLeadEta[nPtBins];
  TH1D *hBBCEastSum_byEta[nEtaBins];
  TH1D *hBBCEastSum_byPt[nPtBins];
  TH1D *hBBCEastSum_byPtCorrected[nPtBins];
  TH1D *hEAdist[nPtBins][nEtaBins];
  TH1D *hEAdistCORRECTED[nPtBins][nEtaBins];

  TH1D *hpt = new TH1D("hpt",";p_{T} [GeV];#frac{1}{#it{N}_{jets}}#frac{d#it{N}_{jets}}{dp_{T}}",50,8,30);
  TH1D *hptc = new TH1D("hptc",";p_{T} [GeV];#frac{1}{#it{N}_{jets}}#frac{d#it{N}_{jets}}{dp_{T}}",50,8,30);
  hpt->SetStats(0);
  hptc->SetStats(0);  

  for ( int e=0; e<nEtaBins; ++e ) {

    name = "hBBCEastSum_byEta" + etaBinName[e];
    title = "BBC ADC East Sum:  " + etaBinString[e] + ";BBC East Sum";
    hBBCEastSum_byEta[e] = new TH1D( name, title, 12,3000,70000 );
    hBBCEastSum_byEta[e]->SetMarkerStyle( etaMarker[e] );
    hBBCEastSum_byEta[e]->SetMarkerColor( etaColor[e] );
    hBBCEastSum_byEta[e]->SetLineColor( etaColor[e] );

    name = "hBBCEastSum_byPt" + ptBinName[e];
    title = "BBC ADC East Sum:  " + ptBinString[e] + ";BBC East Sum";
    hBBCEastSum_byPt[e] = new TH1D( name, title, 12,3000,70000 );
    hBBCEastSum_byPt[e]->SetMarkerStyle( ptMarker[e] );
    hBBCEastSum_byPt[e]->SetMarkerColor( ptColor[e] );
    hBBCEastSum_byPt[e]->SetLineColor( ptColor[e] );
    
    name = "hBBCEastSum_byPtCorrected" + ptBinName[e];
    title = "BBC ADC East Sum:  " + ptCorrectedBinString[e] + ";BBC East Sum";
    hBBCEastSum_byPtCorrected[e] = new TH1D( name, title, 12,3000,70000 );
    hBBCEastSum_byPtCorrected[e]->SetMarkerStyle( 20 );
    hBBCEastSum_byPtCorrected[e]->SetMarkerColor( ptColor[e] );
    hBBCEastSum_byPtCorrected[e]->SetLineColor( ptColor[e] );
    
    for ( int p=0; p<nPtBins; ++p ) {
      name = "hEAdist" + ptBinName[p] + etaBinName[e];
      title = "BBC ADC East Sum:  " + ptBinString[p] + "     " + etaBinString[e] + ";BBC East Sum";
      hEAdist[p][e] = new TH1D( name, title, 18,0,70000 );
      hEAdist[p][e]->SetLineColor( etaColor[e] );
      hEAdist[p][e]->SetMarkerColor( etaColor[e] );
      hEAdist[p][e]->SetMarkerStyle( etaMarker[e] );


      name = "hEAdistCORRECTED" + ptBinName[p] + etaBinName[e];
      title = "BBC ADC East Sum:  " + ptCorrectedBinString[p] + "     " + etaBinString[e] + ";BBC East Sum";
      hEAdistCORRECTED[p][e] = new TH1D( name, title, 18,0,70000 );
      hEAdistCORRECTED[p][e]->SetLineColor( etaColor[e] );
      hEAdistCORRECTED[p][e]->SetMarkerColor( etaColor[e] );
      hEAdistCORRECTED[p][e]->SetMarkerStyle( etaMarker[e] );
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
    
    // chgRho = ( chgEastRho + chgMidRho + chgWestRho )/3;
    // neuRho = ( neuEastRho + neuMidRho + neuWestRho )/3;
    // rho = chgRho + neuRho;

    hRho->Fill(rho_te);
    hTowersVsRho->Fill(rho_te,nTowers);
    hLeadJetPtRhoEta->Fill(leadPt,rho_te,leadEta);
    hLeadPtEtaPhi->Fill(leadPt,leadEta,leadPhi);
    hPt_UE_BBCsumE->Fill(leadPt,rho_te,BbcAdcSumEast);
    hBBCsumE->Fill(BbcAdcSumEast);
    BBCEintegral = hBBCsumE->Integral(0,i);
    hBBCsumE_integral->Fill( BBCEintegral );

    hpt->Fill(leadPt);
    hptc->Fill(leadPtCorrected);

    if ( BbcAdcSumEast>3559.12 && BbcAdcSumEast<10126.1 ) { // LO: 3559.12-10126.1;  HI: 26718.1+
      hleadEta_LoEA->Fill(leadEta);
    }

    if ( BbcAdcSumEast>26718.1 ) {
      hleadEta_HiEA->Fill(leadEta);
    }
	
    pval = 99;    jeval = 99;    eaval = 99;
    

    if ( leadPtCorrected >= 10.0 && leadPtCorrected <= 30.0 ) {
      for ( int ea=0; ea<3; ++ea ) {
	if ( BbcAdcSumEast > BBCEsumLo[ea]  &&  BbcAdcSumEast < BBCEsumHi[ea] ) { eaval = ea; }
      }    
      for ( int p=0; p<3; ++p ) {
	if ( leadPtCorrected >= ptLo[p]  &&  leadPtCorrected <= ptHi[p] ) { pval = p; }
      }
      for ( int je=0; je<3; ++je ) {
	if ( leadEta >= etaLo[je]  &&  leadEta <= etaHi[je] ) { jeval = je; }
      }
      if ( eaval==99 ) { continue; }
      if ( pval==99 || jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<leadPt<<endl<<endl; }
    }
    else {continue;}
    
    hBBCEastSum_byEta[jeval]->Fill( BbcAdcSumEast );
    hBBCEastSum_byPt[pval]->Fill( BbcAdcSumEast );
    hLeadEta[pval]->Fill( leadEta );
    hEAdist[pval][jeval]->Fill( BbcAdcSumEast );

    pval = 99;
    for ( int p=0; p<3; ++p ) {      if ( leadPtCorrected >= ptLo[p]  &&  leadPtCorrected <= ptHi[p] ) { pval = p; }    }
    if ( pval==99 ) { continue; }
    hBBCEastSum_byPtCorrected[pval]->Fill( BbcAdcSumEast );
    hEAdistCORRECTED[pval][jeval]->Fill( BbcAdcSumEast );
  }


  for ( int i=0; i<nEntries; ++i ) {
    jetTree->GetEntry(i);
    hLeadPtVsPtCorrected->Fill( leadPt, leadPtCorrected );
  }


  
  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS 0


  c0->SetLogz();

  hLeadPtVsPtCorrected->Scale(1./hLeadPtVsPtCorrected->Integral("WIDTH"));
  hLeadPtVsPtCorrected->Draw("COLZ");
  c0->SaveAs( (plotdir + "leadPtVsPtCorrected.pdf").c_str() , "PDF" );
  
  // hBGchg->Scale(1./hBGchg->Integral("WIDTH"));
  // hBGchg->Draw("COLZ");
  // c0->SaveAs( "plots/UE/10cmVzCut/chgBgEtaPhi.pdf").c_str() , "PDF" );
  
  // hBGneu->Scale(1./hBGneu->Integral("WIDTH"));
  // hBGneu->Draw("COLZ");
  // c0->SaveAs( "plots/UE/10cmVzCut/neuBgEtaPhi.pdf").c_str() , "PDF" );

  hTowersVsRho->Scale(1./hTowersVsRho->Integral("WIDTH"));
  hTowersVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hTowersVsRho->Draw("COLZ");
  c0->SaveAs( (plotdir + "towersVsRho.pdf").c_str() , "PDF" );

  c0->SetLogy();

  hRho->Scale(1./hRho->Integral("WIDTH"));
  hRho->Draw();
  c0->SaveAs( (plotdir + "rho.pdf").c_str() , "PDF" );
  

  TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 1
  c1->SetMargin(0.14,0.05,0.13,0.05);
  TLegend *leg0 = new TLegend(0.48, 0.66, .92, .95,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(0);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColorAlpha(0,0.0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);
  leg0->AddEntry((TObject*)0,"", "");
  leg0->AddEntry((TObject*)0,"<iBBCE sum>", "");

  TString lmean1[nEtaBins] = { "23430", "23080", "22700" };
  TString lerr1[nEtaBins] = { "50", "40", "50" };
  
  TH2D *sBBCbyEta = new TH2D("sBBCbyEta", ";iBBC East Sum;#frac{1}{#it{N}_{jets}} #frac{d#it{N}_{jets}}{diBBCEsum}", 35,0,70000, 10,0,0.25);
  sBBCbyEta->SetStats(0);
  sBBCbyEta->GetXaxis()->SetLabelSize(.04);
  sBBCbyEta->GetYaxis()->SetLabelSize(.045);
  sBBCbyEta->GetXaxis()->SetTitleSize(0.045);
  sBBCbyEta->GetXaxis()->SetTitleOffset(1.2);
  sBBCbyEta->GetYaxis()->SetTitleOffset(1.4);
  sBBCbyEta->GetYaxis()->SetTitleSize(0.04);
  sBBCbyEta->Draw();
  for ( int e=0; e<nEtaBins; ++e ) {
    hBBCEastSum_byEta[e]->SetStats(0);
    hBBCEastSum_byEta[e]->Scale(1./hBBCEastSum_byEta[e]->Integral());
    hBBCEastSum_byEta[e]->Draw("SAME");
    avg = "";
    avg += hBBCEastSum_byEta[e]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    name = "hBBCEastSum_byEta" + etaBinName[e];
    title = etaBinString[e];
    leg0->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND


    //avg = lmean1[e] + " #pm " + lerr1[e];
    leg0->AddEntry((TObject*)0,avg, "");

  }
  leg0->Draw();
  
  TPaveText *sp1 = new TPaveText(0.1,0.8,0.5,1.0,"NDC");
  //sp1->AddText("STAR Preliminary");
  sp1->AddText("");
  sp1->SetTextSize(0.04);
  sp1->SetTextColor(kRed);
  sp1->SetFillColorAlpha(0,0.0);
  sp1->SetLineColorAlpha(0,0.0);
  sp1->Draw("NB");

  auto jetText = new TLatex(0.57,0.52,"jets uncorrected for detector effects");
  jetText->SetTextFont(63);
  jetText->SetTextSize(16);
  jetText->SetTextColor(kBlack);
  jetText->SetLineWidth(1);
  jetText->SetNDC();
  jetText->Draw();
  jetText->DrawLatex(0.65,0.46,"statistical errors only");

    
  for (int i=0; i<5; ++i) {
    tex[i]->DrawLatex(0.2,0.45-0.05*i,textString[i]);
  }
  c1->SaveAs( (plotdir + "BBCEastSum_by_eta.pdf").c_str() , "PDF" );

  c1->SetLogy(0);

  
  TString bmean[nPtBins] = { "22780", "21870", "21200" };
  TString berror[nPtBins] = { "30", "60", "100" };

  TLegend *leg0a = new TLegend(0.48, 0.66, .92, .95,NULL,"brNDC");    // LEGEND 0
  leg0a->SetBorderSize(0);   leg0a->SetLineColor(1);   leg0a->SetLineStyle(1);   leg0a->SetLineWidth(1);   leg0a->SetFillColorAlpha(0,0.0);   leg0a->SetFillStyle(1001);
  leg0a->SetNColumns(2);
  leg0a->AddEntry((TObject*)0,"", "");
  leg0a->AddEntry((TObject*)0,"<BBCE sum>", "");
  
  TH2D *sBBCbyPt = new TH2D("sBBCbyPt", ";iBBC East Sum;#frac{1}{#it{N}_{jets}} #frac{d#it{N}_{jets}}{diBBCEsum}", 35,0,70000, 10,0,0.25);
  sBBCbyPt->GetXaxis()->SetLabelSize(.04);
  sBBCbyPt->GetYaxis()->SetLabelSize(.045);
  sBBCbyPt->GetXaxis()->SetTitleSize(0.045);
  sBBCbyPt->GetXaxis()->SetTitleOffset(1.2);
  sBBCbyPt->GetYaxis()->SetTitleOffset(1.4);
  sBBCbyPt->GetYaxis()->SetTitleSize(0.04);
  sBBCbyPt->SetStats(0);
  sBBCbyPt->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hBBCEastSum_byPtCorrected[p]->SetStats(0);
    hBBCEastSum_byPtCorrected[p]->Scale(1./hBBCEastSum_byPtCorrected[p]->Integral());
    hBBCEastSum_byPtCorrected[p]->Draw("SAME");
    avg = "";
    // avg += hBBCEastSum_byPtCorrected[p]->GetMean(1);                                           // 1 denotes x-axis
    // avg = avg(0,5);
    name = "hBBCEastSum_byPtCorrected" + ptBinName[p];
    title = ptBinString[p];
    leg0a->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND

    TString meanError = "";
    // meanError += hBBCEastSum_byPtCorrected[p]->GetMeanError();
    // meanError = meanError(0,6);
    avg += hBBCEastSum_byPtCorrected[p]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    // avg = bmean[p] + " #pm " + berror[p];
    leg0a->AddEntry((TObject*)0,avg, "");

  }
  leg0a->Draw();

  jetText->DrawLatex(0.57,0.52,"jets uncorrected for detector effects");
  jetText->DrawLatex(0.65,0.46,"statistical errors only");

  sp1->Draw("NB");

  for (int i=0; i<4; ++i) {
    tex[i]->DrawLatex(0.2,0.4-0.05*i,textString[i]);
  }



  
  c1->SaveAs( (plotdir + "BBCEastSum_by_ptCorrected.pdf").c_str() , "PDF" );
  
  c1->SetLogy(0);


  
  // TLegend *leg0b = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  // leg0b->SetBorderSize(1);   leg0b->SetLineColor(1);   leg0b->SetLineStyle(1);   leg0b->SetLineWidth(1);   leg0b->SetFillColor(0);   leg0b->SetFillStyle(1001);
  // leg0b->SetNColumns(2);
  // leg0b->AddEntry((TObject*)0,"#bf{p_{T}^{corrected}}", "");
  // leg0b->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");
  
  // TH2D *sBBCbyPtCorrected = new TH2D("sBBCbyPtCorrected", "BBC ADC East Sum by Corrected Lead Jet p_{T};BBC East Sum", 35,0,70000, 10,0,0.15);
  // sBBCbyPtCorrected->SetStats(0);
  // sBBCbyPtCorrected->Draw();
  // for ( int p=0; p<nPtBins; ++p ) {
  //   hBBCEastSum_byPtCorrected[p]->SetStats(0);
  //   hBBCEastSum_byPtCorrected[p]->Scale(1./hBBCEastSum_byPtCorrected[p]->Integral());
  //   hBBCEastSum_byPtCorrected[p]->Draw("SAME");
  //   avg = "";
  //   avg += hBBCEastSum_byPtCorrected[p]->GetMean(1);                                           // 1 denotes x-axis
  //   avg = avg(0,5);
  //   name = "hBBCEastSum_byPtCorrected" + ptBinName[p];
  //   title = ptCorrectedBinString[p];
  //   leg0b->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
  //   leg0b->AddEntry((TObject*)0,avg, "");

  // }
  // leg0b->Draw();
  // c1->SaveAs( (plotdir + "BBCEastSum_by_ptCorrected.pdf").c_str() , "PDF" );
  
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
  c1->SaveAs( (plotdir + "LeadEta_by_pt.pdf").c_str() , "PDF" );


  jetTree->Draw("leadPt:rho>>hRho2d","","COLZ");
  TH2D *hRho2d = (TH2D*)gDirectory->Get("hRho2d");
  c1->SetLogy();
  TLegend *leg2 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(3);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho[nPtBins];
  TH2D *sPtRho = new TH2D( "sPtRho", "Underlying Event by Lead Jet p_{T};#rho [GeV]", 40,0,12, 10,0.000001,1 );
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
  c1->SaveAs((plotdir + "rhoByLeadPt.pdf").c_str(),"PDF");





  // LO: 3559.12-10126.1
  jetTree->Draw("leadPt:rho_te>>hRho2d_LO","BbcAdcSumEast>3559.12 && BbcAdcSumEast<10126.1","COLZ");
  TH2D *hRho2d_LO = (TH2D*)gDirectory->Get("hRho2d_LO");
  c1->SetLogy();
  TLegend *leg3 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg3->SetBorderSize(1);   leg3->SetLineColor(1);   leg3->SetLineStyle(1);   leg3->SetLineWidth(1);   leg3->SetFillColor(0);   leg3->SetFillStyle(1001);
  leg3->SetNColumns(3);
  leg3->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg3->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg3->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho_LO[nPtBins];
  TH2D *sPtRho_LO = new TH2D( "sPtRho_LO", "Low EA: Underlying Event by Lead Jet p_{T};#rho [GeV]", 40,0,12, 10,0.000001,1 );
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

  // TLatex *EAtext = new TLatex(0.65, 0.85,EAbinString[0]);
  // EAtext->SetTextFont(63);
  // EAtext->SetTextSize(16);
  // EAtext->SetTextColor(kBlack);
  // EAtext->SetLineWidth(1);
  // EAtext->SetNDC();
  // EAtext->Draw();

  c1->SaveAs((plotdir + "rhoByLeadPt_LOEA.pdf").c_str(),"PDF");




  // HI: 26718.1+
  jetTree->Draw("leadPt:rho_te>>hRho2d_HI","BbcAdcSumEast>26718.1","COLZ");
  TH2D *hRho2d_HI = (TH2D*)gDirectory->Get("hRho2d_HI");
  //c1->SetLogy();
  TLegend *leg4 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg4->SetBorderSize(1);   leg4->SetLineColor(1);   leg4->SetLineStyle(1);   leg4->SetLineWidth(1);   leg4->SetFillColor(0);   leg4->SetFillStyle(1001);
  leg4->SetNColumns(3);
  leg4->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg4->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg4->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho_HI[nPtBins];
  TH2D *sPtRho_HI = new TH2D( "sPtRho_HI", "High EA: Underlying Event by Lead Jet p_{T};#rho [GeV]", 40,0,12, 10,0.000001,1 );
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
  c1->SaveAs((plotdir + "rhoByLeadPt_HIEA.pdf").c_str(),"PDF");





  
  //  corrected lead pt
  
  jetTree->Draw("leadPtCorrected:rho_te>>hRhoCorr2d","","COLZ");
  TH2D *hRhoCorr2d = (TH2D*)gDirectory->Get("hRhoCorr2d");
  c1->SetLogy();
  TLegend *leg12 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg12->SetBorderSize(1);   leg12->SetLineColor(1);   leg12->SetLineStyle(1);   leg12->SetLineWidth(1);   leg12->SetFillColor(0);   leg12->SetFillStyle(1001);
  leg12->SetNColumns(3);
  leg12->AddEntry((TObject*)0,"#bf{p_{T}^{corrected}}", "");
  leg12->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg12->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hCorrPtRho[nPtBins];
  TH2D *sCorrPtRho = new TH2D( "sCorrPtRho", "Underlying Event by Corrected Lead Jet p_{T};#rho [GeV]", 40,0,12, 10,0.000001,1 );
  sCorrPtRho->SetStats(0);
  sCorrPtRho->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho[p] = (TH1D*) hRhoCorr2d->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptCorrectedBinString[p];
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
  c1->SaveAs((plotdir + "rhoByLeadPt_correctedPt.pdf").c_str(),"PDF");
  


  TCanvas * c1a = new TCanvas( "c1a" , "" ,700 ,500 );              // CANVAS 1
  c1a->SetMargin(0.15,0.1,0.12,0.1);

  // LO: 3559.12-10126.1
  jetTree->Draw("leadPtCorrected:rho>>hRhoCorr2d_LO","BbcAdcSumEast>3559.12 && BbcAdcSumEast<10126.1","COLZ");
  TH2D *hRhoCorr2d_LO = (TH2D*)gDirectory->Get("hRhoCorr2d_LO");
  //c1a->SetLogy();
  c1a->SetMargin(0.125,0.05,0.1,0.05);

  TLegend *leg13 = new TLegend(0.57, 0.62, 0.90, 0.95,NULL,"brNDC");    // LEGEND 0
  leg13->SetBorderSize(0);   leg13->SetLineColor(1);   leg13->SetLineStyle(1);   leg13->SetLineWidth(1);   leg13->SetFillColor(0);   leg13->SetFillStyle(1001);
  leg13->SetNColumns(2);
  leg13->AddEntry((TObject*)0,"lead jet p_{T}", "");
  leg13->AddEntry((TObject*)0,"uncorr. <#rho>", "");
  leg13->SetFillColorAlpha(0,0.0);

  TH1D *hCorrPtRho_LO[nPtBins];
  TH2D *sCorrPtRho_LO = new TH2D( "sCorrPtRho_LO", ";#rho [GeV];#frac{1}{#it{N}_{jets}} #frac{d#it{N}_{jets}}{d#rho}", 40,0,4.5, 10,0,0.25 );
  sCorrPtRho_LO->SetStats(0);
  sCorrPtRho_LO->GetXaxis()->SetLabelSize(.04);
  sCorrPtRho_LO->GetYaxis()->SetLabelSize(.04);
  sCorrPtRho_LO->GetXaxis()->SetTitleOffset(1);
  sCorrPtRho_LO->GetYaxis()->SetTitleOffset(1.25);
  sCorrPtRho_LO->GetXaxis()->SetTitleSize(.04);
  sCorrPtRho_LO->GetYaxis()->SetTitleSize(.04);
  sCorrPtRho_LO->Draw();

  //0.7099 ± 0.0027  0.6770 ± 0.0055  0.6419 ± 0.0085
  double lmean[nPtBins] = { 0.710, 0.677, 0.642 };
  double lerr[nPtBins] = { 0.003, 0.006, 0.009 };
    
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d_LO->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho_LO[p] = (TH1D*) hRhoCorr2d_LO->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinStringRaw[p];
    hCorrPtRho_LO[p]->SetNameTitle( name, title );
    hCorrPtRho_LO[p]->SetLineColor( ptColor[p] );
    hCorrPtRho_LO[p]->SetMarkerColor( ptColor[p] );
    hCorrPtRho_LO[p]->SetMarkerStyle( 20 );
    hCorrPtRho_LO[p]->Scale(1./hCorrPtRho_LO[p]->GetEntries());
    hCorrPtRho_LO[p]->SetStats(0);
    hCorrPtRho_LO[p]->Draw("SAME");

    leg13->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND

    avg = "";
    avg+=hCorrPtRho_LO[p]->GetMean(1);
    avg = avg(0,6);
    TString meanError;
    meanError+=hCorrPtRho_LO[p]->GetMeanError();
    meanError = meanError(0,6);

    // avg = "";
    // avg += lmean[p];
    // avg = avg(0,5);
    // TString meanError = "";
    // meanError += lerr[p];
    // meanError = meanError(0,5);
    // avg += " #pm " + meanError;
    leg13->AddEntry((TObject*)0,avg, "");  // FOR PRELIMINARIES!
  }

  TLatex *EAtext = new TLatex(0.3, 0.7,EAbinString[0]);
  EAtext->SetTextFont(63);
  EAtext->SetTextSize(32);
  EAtext->SetTextColor(kBlack);
  EAtext->SetLineWidth(1);
  EAtext->SetNDC();
  EAtext->Draw();

  for (int i=0; i<4; ++i) {
    tex[i]->DrawLatex(0.65,0.35-0.05*i,textString[i]);
  }
  sp1->Draw("NB");

  jetText->DrawLatex(0.57,0.52,"jets uncorrected for detector effects");
  jetText->DrawLatex(0.65,0.46,"statistical errors only");

  leg13->Draw();
  c1a->SaveAs((plotdir + "rhoByLeadPt_LOEA_lin.pdf").c_str(),"PDF");



  // HI: 26718.1+
  jetTree->Draw("leadPtCorrected:rho>>hRhoCorr2d_HI","BbcAdcSumEast>26718.1","COLZ");
  TH2D *hRhoCorr2d_HI = (TH2D*)gDirectory->Get("hRhoCorr2d_HI");
  //c1a->SetLogy();
  c1a->SetMargin(0.125,0.05,0.1,0.05);
  TLegend *leg14 = new TLegend(0.57, 0.62, 0.90, 0.95,NULL,"brNDC");    // LEGEND 0
  leg14->SetBorderSize(0);   leg14->SetLineColor(1);   leg14->SetLineStyle(1);   leg14->SetLineWidth(1);   leg14->SetFillColor(0);   leg14->SetFillStyle(1001);
  leg14->SetNColumns(2);
  leg14->AddEntry((TObject*)0,"lead jet p_{T}", "");
  leg14->AddEntry((TObject*)0,"uncorr. <#rho>", "");

  //  1.1900 ± 0.0027   1.1495 ± 0.0046   1.1001 ± 0.0098
  double lmean2[nPtBins] = { 1.190, 1.150, 1.10 };
  double lerr2[nPtBins] = { 0.003, 0.005, 0.01 };
  
  TH1D *hCorrPtRho_HI[nPtBins];
  TH2D *sCorrPtRho_HI = new TH2D( "sCorrPtRho_HI", ";#rho [GeV];#frac{1}{#it{N}_{jets}} #frac{d#it{N}_{jets}}{d#rho}", 40,0,4.5, 10,0,0.25 );
  sCorrPtRho_HI->SetStats(0);
  sCorrPtRho_HI->GetXaxis()->SetLabelSize(.04);
  sCorrPtRho_HI->GetYaxis()->SetLabelSize(.04);
  sCorrPtRho_HI->GetXaxis()->SetTitleOffset(1);
  sCorrPtRho_HI->GetYaxis()->SetTitleOffset(1.25);
  sCorrPtRho_HI->GetXaxis()->SetTitleSize(.04);
  sCorrPtRho_HI->GetYaxis()->SetTitleSize(.04);
  sCorrPtRho_HI->Draw();
  for ( int p=0; p<nPtBins; ++p ) {
    hRhoCorr2d_HI->GetYaxis()->SetRangeUser( ptLo[p], ptHi[p] );
    hCorrPtRho_HI[p] = (TH1D*) hRhoCorr2d_HI->ProjectionX();
    name = "hUE" + ptBinName[p];
    title =  ptBinStringRaw[p];
    hCorrPtRho_HI[p]->SetNameTitle( name, title );
    hCorrPtRho_HI[p]->SetLineColor( ptColor[p] );
    hCorrPtRho_HI[p]->SetMarkerColor( ptColor[p] );
    hCorrPtRho_HI[p]->SetMarkerStyle( 20 );
    hCorrPtRho_HI[p]->Scale(1./hCorrPtRho_HI[p]->GetEntries());
    hCorrPtRho_HI[p]->SetStats(0);
    hCorrPtRho_HI[p]->Draw("SAME");

    leg14->SetFillColorAlpha(0,0.0);
    leg14->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND

    avg = "";
    avg += hCorrPtRho_HI[p]->GetMean(1);
    avg = avg(0,6);
    TString meanError = "";
    meanError += hCorrPtRho_HI[p]->GetMeanError();
    meanError = meanError(0,6);
    avg += " #pm " + meanError;

    // avg = "";
    // avg += lmean2[p];
    // if (p==2) { avg = avg(0,4); }
    // else { avg = avg(0,5); }
    // TString meanError = "";
    // meanError += lerr2[p];
    // if (p==2) { meanError = meanError(0,4); }
    // else { meanError = meanError(0,5); }
    // avg += " #pm " + meanError;
    leg14->AddEntry((TObject*)0,avg, "");   //  FOR PRELIMINARIES ONLY!
  }

  sp1->Draw("NB");
  
  EAtext->DrawText(0.3, 0.7,EAbinString[1]);
  //0.7,0.15
  for (int i=0; i<4; ++i) {
    tex[i]->DrawLatex(0.65,0.35-0.05*i,textString[i]);
  }

  jetText->DrawLatex(0.57,0.52,"jets uncorrected for detector effects");
  jetText->DrawLatex(0.65,0.46,"statistical errors only");
		   
  leg14->Draw();
  c1a->SaveAs((plotdir + "rhoByLeadPt_HIEA_lin.pdf").c_str(),"PDF");


  
 

  c1->SetLogy();

  // hBBCsumE->SetBinError(3559.12,100000000000);
  // hBBCsumE->SetBinError(7802,100000000000);
  // hBBCsumE->SetBinError(10126.1,100000000000);
  // hBBCsumE->SetBinError(15399,100000000000);
  // hBBCsumE->SetBinError(19525,100000000000);;
  // hBBCsumE->SetBinError(23781,100000000000);
  // hBBCsumE->SetBinError(26718.1,100000000000);
  // hBBCsumE->SetBinError(34082,100000000000);
  // hBBCsumE->SetBinError(41108,100000000000);
  // hBBCsumE->SetBinError(63974,100000000000);
  // hBBCsumE->Draw();

  TCanvas * c2 = new TCanvas( "c2" , "" ,700 ,500 );              // CANVAS 1
  c2->SetMargin(0.12,0.1,0.12,0.1);
  
  TLegend *leg5[nPtBins];
  
  for ( int p=0; p<nPtBins; ++p ) {
    leg5[p] = new TLegend(0.45, 0.45, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
    leg5[p]->SetBorderSize(0);   leg5[p]->SetLineColor(1);   leg5[p]->SetLineStyle(1);   leg5[p]->SetLineWidth(1);   leg5[p]->SetFillColor(0);   leg5[p]->SetFillStyle(1001);
    leg5[p]->SetNColumns(2);
    leg5[p]->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
    leg5[p]->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");

    title = "BBC ADC East Sum by Lead Jet #eta:    " + ptBinString[p] + ";BBC East Sum";
    sBBCbyEta->SetTitle( title );
    sBBCbyEta->GetXaxis()->SetLabelSize(.04);
    sBBCbyEta->GetYaxis()->SetLabelSize(.04);
    sBBCbyEta->GetXaxis()->SetTitleOffset(1);
    sBBCbyEta->GetYaxis()->SetTitleOffset(1.2);
    sBBCbyEta->GetXaxis()->SetTitleSize(.04);
    sBBCbyEta->GetYaxis()->SetTitleSize(.04);
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
    saveName = plotdir + "BBCEastSum_by_eta" + ptBinName[p] +".pdf";
    c2->SaveAs( saveName , "PDF" );
  }





  TLegend *leg6[nPtBins];
  
  for ( int p=0; p<nPtBins; ++p ) {
    leg6[p] = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
    leg6[p]->SetBorderSize(1);   leg6[p]->SetLineColor(1);   leg6[p]->SetLineStyle(1);   leg6[p]->SetLineWidth(1);   leg6[p]->SetFillColor(0);   leg6[p]->SetFillStyle(1001);
    leg6[p]->SetNColumns(2);
    leg6[p]->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
    leg6[p]->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");

    title = "BBC ADC East Sum by Lead Jet #eta:    " + ptCorrectedBinString[p] + ";BBC East Sum";
    sBBCbyEta->SetTitle( title );
    //TH2D *sBBCbyEta = new TH2D("sBBCbyEta", title, 35,0,70000, 10,0,0.15);
    sBBCbyEta->SetStats(0);
    sBBCbyEta->Draw();
    for ( int e=0; e<nEtaBins; ++e ) {
      hEAdistCORRECTED[p][e]->SetStats(0);
      hEAdistCORRECTED[p][e]->Scale(1./hEAdistCORRECTED[p][e]->Integral());
      hEAdistCORRECTED[p][e]->Draw("SAME");
      avg = "";
      avg += hEAdistCORRECTED[p][e]->GetMean(1);                                           // 1 denotes x-axis
      avg = avg(0,5);
      name = "hEAdistCORRECTED" + ptBinName[p] + etaBinName[e];
      title = etaBinString[e];
      leg6[p]->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
      leg6[p]->AddEntry((TObject*)0,avg, "");

    }
    leg6[p]->Draw();
    saveName = plotdir + "BBCEastSum_by_eta" + ptBinName[p] +"_corrected.pdf";
    c2->SaveAs( saveName , "PDF" );
  }

  hpt->Scale(1./hpt->Integral());
  hpt->SetLineColor(kBlue);
  hpt->SetMarkerColor(kBlue);
  hpt->SetMarkerStyle(20);
  hptc->Scale(1./hptc->Integral());
  hptc->SetLineColor(kRed);
  hptc->SetMarkerColor(kRed);
  hptc->SetMarkerStyle(20);
  
  TCanvas * c3 = new TCanvas( "c3" , "" ,500 ,700 );
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetLogy();
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.18);
  pad1->Draw();
  pad1->cd();
  hpt->GetYaxis()->SetTitleOffset(2);
  hpt->GetYaxis()->SetTitleSize(0.04);
  hpt->Draw("P");
  hptc->Draw("PSAME");

  TH1D *hPtRatio = (TH1D*) hptc->Clone("hPtRatio");
  hPtRatio->GetYaxis()->SetRangeUser(0.8,1.4);
  hPtRatio->GetYaxis()->SetLabelSize(0.1);
  hPtRatio->GetXaxis()->SetLabelSize(0.1);
  hPtRatio->GetYaxis()->SetTitle("p_{T}^{corrected}/p_{T}");
  hPtRatio->GetXaxis()->SetTitle("p_{T}");
  hPtRatio->GetXaxis()->SetTitleSize(0.13);
  hPtRatio->GetYaxis()->SetTitleSize(0.13);
  hPtRatio->Divide(hpt);
  hPtRatio->SetStats(0);
  hPtRatio->SetLineColor(kBlack);
  hPtRatio->SetMarkerColor(kBlack);

  c3->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetLeftMargin(0.18);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  hPtRatio->GetYaxis()->SetTitleSize(0.13);
  hPtRatio->GetXaxis()->SetTitleSize(0.13);
  hPtRatio->GetYaxis()->SetTitleOffset(0.5);
  hPtRatio->GetXaxis()->SetTitleOffset(0.6);
  hPtRatio->Draw("P");

  pad1->cd();
  //pad1->BuildLegend();
  TLegend *leg7;
  
  leg7 = new TLegend(0.4, 0.65, 0.89, 0.89,NULL,"brNDC");    // LEGEND 0
  leg7->SetBorderSize(0);   leg7->SetLineColor(1);   leg7->SetLineStyle(1);   leg7->SetLineWidth(1);   leg7->SetFillColor(0);   leg7->SetFillStyle(1001);
  leg7->AddEntry( hpt, "p_{T}", "lpf" );
  leg7->AddEntry( hptc,"Corrected p_{T}", "lpf" );
  leg7->Draw();
  
  c3->SaveAs((plotdir + "CorrectedPtRatioPlot.pdf").c_str(),"PDF");


}
