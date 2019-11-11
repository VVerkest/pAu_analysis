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
  const int ptColor[nPtBins] = { 797, 593, 875 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-1.0<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<1.0" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

  int eval, pval;
  TString name, saveName, title;
  
  TString fileName = "out/UE/pAuHTjetUE.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH2D *hBGchg = (TH2D*) inFile->Get("hChgBgEtaPhi");
  TH2D *hBGneu = (TH2D*) inFile->Get("hNeuBgEtaPhi");

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  jetTree->SetBranchAddress( "RunID", &RunID );      	       		jetTree->SetBranchAddress( "EventID", &EventID );					jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );       		jetTree->SetBranchAddress( "nGlobal", &nGlobal );					jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );			jetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );	       			jetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );	       		jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );			jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );

  int nEntries = jetTree->GetEntries();

  TH1D *hRho = new TH1D("hRho","Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hLeadPhi = new TH1D("hLeadPhi","Lead Jet #phi;#phi_{lead}",12,0,2*pi);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 60,0,15, 200,0,200 );
  TH3D *hLeadJetPtRhoEta = new TH3D( "hLeadJetPtRhoEta", "Lead Jet p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );  
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 160,0,80000 );

  TH1D *hLeadEta[nPtBins];
  TH1D *hBBCEastSum[nEtaBins];

  for ( int e=0; e<nEtaBins; ++e ) {

    name = "hBBCEastSum" + etaBinName[e];
    title = "BBC ADC East Sum:  " + etaBinString[e] + ";BBC East Sum";
    hBBCEastSum[e] = new TH1D( name, title, 70,0,70000 );
    hBBCEastSum[e]->SetMarkerStyle( etaMarker[e] );
    //hBBCEastSum[e]->SetMarkerSize(2);
  }

  for ( int p=0; p<nPtBins; ++p ) {
    name = "hLeadEta" + ptBinName[p];
    title = "Lead Jet #eta:  " + ptBinString[p] + ";#eta_{lead}";
    hLeadEta[p] = new TH1D( name, title, 40, -1.0, 1.0 );
    hLeadEta[p]->SetLineColor( ptColor[p] );
    hLeadEta[p]->SetMarkerColor( ptColor[p] );
    hLeadEta[p]->SetMarkerSize(2);
  }
  
  double chgRho, neuRho, midRho, eastRho, westRho, rho, sigma;
  
  for ( int i=0; i<nEntries; ++i ) {

    jetTree->GetEntry(i);

    chgRho = ( chgEastRho + chgMidRho + chgWestRho )/3;
    neuRho = ( neuEastRho + neuMidRho + neuWestRho )/3;
    rho = chgRho + neuRho;

    hRho->Fill(rho);
    hLeadPhi->Fill(leadPhi);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadJetPtRhoEta->Fill(leadPt,rho,leadEta);
    hLeadPtEtaPhi->Fill(leadPt,leadEta,leadPhi);
    hPt_UE_BBCsumE->Fill(leadPt,rho,BbcAdcEastSum);

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
  c0->SaveAs( "plots/UE/leadPhi.pdf" , "PDF" );

  c0->SetLogz();

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
    

  TCanvas * c1 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS 1
  
  for ( int e=0; e<nEtaBins; ++e ) {
    hBBCEastSum[e]->SetStats(0);
    hBBCEastSum[e]->Scale(1./hBBCEastSum[e]->Integral("WIDTH"));
    hBBCEastSum[e]->Draw("SAME");
  }
  c1->SaveAs( "plots/UE/BBCEastSum_by_eta.pdf" , "PDF" );

  TCanvas * c2 = new TCanvas( "c2" , "" ,700 ,500 );              // CANVAS 2
  for ( int e=0; e<nEtaBins; ++e ) {
    for ( int p=0; p<nPtBins; ++p ) {
      hLeadEta[p]->Scale(1./hLeadEta[p]->Integral("WIDTH"));
      hLeadEta[p]->Draw("SAME");
    }
  }
  c1->SaveAs( "plots/UE/LeadEta_by_pteta.pdf" , "PDF" );
    
}
