// Veronica Verkest
// November 12, 2019

void UEdijetPlot(){

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
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

  int eval, pval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, midRho, eastRho, westRho, rho;
  
  TString fileName = "out/UE/pAuHTdijetUE_noRecoEtaMatchReq.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH2D *hBGchg = (TH2D*) inFile->Get("hChgBgEtaPhi");
  TH2D *hBGneu = (TH2D*) inFile->Get("hNeuBgEtaPhi");

  TTree *dijetTree = (TTree*) inFile->Get("HTdijetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  dijetTree->SetBranchAddress( "RunID", &RunID );      	       	dijetTree->SetBranchAddress( "EventID", &EventID );		       		dijetTree->SetBranchAddress( "nTowers", &nTowers );
  dijetTree->SetBranchAddress( "nPrimary", &nPrimary );       	dijetTree->SetBranchAddress( "nGlobal", &nGlobal );		       		dijetTree->SetBranchAddress( "nVertices", &nVertices );
  dijetTree->SetBranchAddress( "refMult", &refMult );	       	dijetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		dijetTree->SetBranchAddress( "Vz", &Vz );
  dijetTree->SetBranchAddress( "leadPt", &leadPt );	       	       	dijetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	dijetTree->SetBranchAddress( "leadEta", &leadEta );
  dijetTree->SetBranchAddress( "leadPhi", &leadPhi );	       	dijetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		dijetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  dijetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	dijetTree->SetBranchAddress( "neuEastRho", &neuEastRho );     	dijetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  dijetTree->SetBranchAddress( "neuWestRho", &neuWestRho );

  int nEntries = dijetTree->GetEntries();

  TH1D *hRho = new TH1D("hRho","Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hLeadPhi = new TH1D("hLeadPhi","Lead Jet #phi;#phi_{lead}",12,0,2*pi);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 60,0,15, 200,0,200 );
  TH3D *hLeadJetPtRhoEta = new TH3D( "hLeadJetPtRhoEta", "Lead Jet p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );  
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 80,0,80000 );
  TH1D *hleadEta_LoEA = new TH1D( "hleadEta_LoEA", "Low Event Activity", 40,-1,1 );
  TH1D *hleadEta_HiEA = new TH1D( "hleadEta_HiEA", "High Event Activity", 40,-1,1 );

  TH1D *hLeadEta[nPtBins];
  TH1D *hBBCEastSum[nEtaBins];

  for ( int e=0; e<nEtaBins; ++e ) {

    name = "hBBCEastSum" + etaBinName[e];
    title = "BBC ADC East Sum:  " + etaBinString[e] + ";BBC East Sum";
    hBBCEastSum[e] = new TH1D( name, title, 18,0,70000 );
    hBBCEastSum[e]->SetMarkerStyle( etaMarker[e] );
    hBBCEastSum[e]->SetMarkerColor( etaColor[e] );
    hBBCEastSum[e]->SetLineColor( etaColor[e] );
  }

  for ( int p=0; p<nPtBins; ++p ) {
    name = "hLeadEta" + ptBinName[p];
    title = "Lead Jet #eta:  " + ptBinString[p] + ";#eta_{lead}";
    hLeadEta[p] = new TH1D( name, title, 40, -1.0, 1.0 );
    hLeadEta[p]->SetLineColor( ptColor[p] );
    hLeadEta[p]->SetMarkerColor( ptColor[p] );
  }
    
  for ( int i=0; i<nEntries; ++i ) {

    dijetTree->GetEntry(i);

    chgRho = ( chgEastRho + chgMidRho + chgWestRho )/3;
    neuRho = ( neuEastRho + neuMidRho + neuWestRho )/3;
    rho = chgRho + neuRho;

    hRho->Fill(rho);
    hLeadPhi->Fill(leadPhi);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadJetPtRhoEta->Fill(leadPt,rho,leadEta);
    hLeadPtEtaPhi->Fill(leadPt,leadEta,leadPhi);
    hPt_UE_BBCsumE->Fill(leadPt,rho,BbcAdcEastSum);

    if ( BbcAdcEastSum>8000 && BbcAdcEastSum<16000 ) {
      hleadEta_LoEA->Fill(leadEta);
    }

    if ( BbcAdcEastSum>30000 ) {
      hleadEta_HiEA->Fill(leadEta);
    }
	
    pval = 99;    eval = 99;
    
    for ( int p=0; p<3; ++p ) {
      if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
    }
    for ( int e=0; e<3; ++e ) {
      if ( leadEta >= etaLo[e]  &&  leadEta <= etaHi[e] ) { eval = e; }
    }
    if ( pval==99 || eval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl; }

    hBBCEastSum[eval]->Fill( BbcAdcEastSum );
    hLeadEta[pval]->Fill( leadEta );
    
  }

  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS 0

  hLeadPhi->Scale(1./hLeadPhi->Integral("WIDTH"));
  hLeadPhi->Draw();
  c0->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_leadPhi.pdf" , "PDF" );

  c0->SetLogz();

  hBGchg->Scale(1./hBGchg->Integral("WIDTH"));
  hBGchg->Draw("COLZ");
  c0->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_chgBgEtaPhi.pdf" , "PDF" );
  
  hBGneu->Scale(1./hBGneu->Integral("WIDTH"));
  hBGneu->Draw("COLZ");
  c0->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_neuBgEtaPhi.pdf" , "PDF" );

  hTowersVsRho->Scale(1./hTowersVsRho->Integral("WIDTH"));
  hTowersVsRho->GetZaxis()->SetRangeUser(0.000001,1);
  hTowersVsRho->Draw("COLZ");
  c0->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_towersVsRho.pdf" , "PDF" );

  c0->SetLogy();

  hRho->Scale(1./hRho->Integral("WIDTH"));
  hRho->Draw();
  c0->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_rho.pdf" , "PDF" );


  TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 1
  
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);
  leg0->AddEntry((TObject*)0,"#bf{#eta_{lead}}", "");
  leg0->AddEntry((TObject*)0,"#bf{<BBCE sum>}", "");
  
  TH2D *sBBCbyEta = new TH2D("sBBCbyEta", "BBC ADC East Sum by Lead Jet #eta;BBC East Sum", 10,0,70000, 10,0,0.15);
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
  c1->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_BBCEastSum_by_eta.pdf" , "PDF" );
  c1->SetLogy(0);


  TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  leg1->SetNColumns(2);
  leg1->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg1->AddEntry((TObject*)0,"#bf{<#eta>}", "");
  
  TH2D *sLeadEtaByPt = new TH2D("sLeadEtaByPt", "Lead Jet #eta by p_{T};#eta_{lead}", 40,-1,1, 20,0.5,1.5);
  sLeadEtaByPt->SetStats(0);
  sLeadEtaByPt->Draw();

  for ( int p=0; p<nPtBins; ++p ) {
    hLeadEta[p]->SetStats(0);

    for (int j=1;j<15;++j)   {double val=hLeadEta[p]->GetBinContent(j); hLeadEta[p]->SetBinContent(j,val/0.3);}
    for (int j=15;j<27;++j) {double val=hLeadEta[p]->GetBinContent(j); hLeadEta[p]->SetBinContent(j,val/0.6);}
    for (int j=27;j<40;++j) {double val=hLeadEta[p]->GetBinContent(j); hLeadEta[p]->SetBinContent(j,val/0.3);}    
    
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
  c1->SaveAs( "plots/UE/dijet_noRecoEtaMatchReq_LeadEta_by_pt.pdf" , "PDF" );


  dijetTree->Draw("leadPt:((chgEastRho+neuEastRho)+(chgMidRho+neuMidRho)+(chgWestRho+neuWestRho))/3>>hRho2d","","COLZ");
  TH2D *hRho2d = (TH2D*)gDirectory->Get("hRho2d");
  c1->SetLogy();
  TLegend *leg2 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND 0
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(3);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{lead}}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#rho>}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#sigma>}", "");

  TH1D *hPtRho[nPtBins];
  TH2D *sPtRho = new TH2D( "sPtRho", "Underlying Event by Lead Jet p_{T};#rho (GeV)", 20,0,6, 10,0.0001,1 );
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
  c1->SaveAs("plots/UE/dijet_noRecoEtaMatchReq_rhoByLeadPt.pdf","PDF");

}
