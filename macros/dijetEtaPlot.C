{
  gStyle->SetErrorX(0.0);
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };
  const TString emw[nEtaBins] = {"east", "mid", "west"};
  const TString leadEtaBinString[nEtaBins] = { "-0.6<#eta_{jet}^{lead}<-0.3", "-0.3<#eta_{jet}^{lead}<0.3", "0.3<#eta_{jet}^{lead}<0.6" };
  const TString recoEtaBinString[nEtaBins] = { "-0.6<#eta_{jet}^{reco}<-0.3", "-0.3<#eta_{jet}^{reco}<0.3", "0.3<#eta_{jet}^{reco}<0.6" };
  TString name, title;
  
  TFile *inFile = new TFile("out/UE/pAuHTdijetUE_10cmVzCut_6cmVzDiff.root", "READ");

  TTree *dijetTree = (TTree*)inFile->Get("HTdijetTree");

  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, eval, pval, nUEpart_chg, nUEpart_neu, nHTtrig;
  double Vz, BbcAdcSumEast, BbcAdcSumEastOuter, BbcAdcSumWest, BbcAdcSumWestOuter, leadPt, recoPt, leadEta, recoEta, leadPhi, recoPhi,
    chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, recoArea, dRTrigLead, dRTrigReco;

  dijetTree->SetBranchAddress( "RunID", &RunID );
  dijetTree->SetBranchAddress( "EventID", &EventID );
  dijetTree->SetBranchAddress( "nTowers", &nTowers );
  dijetTree->SetBranchAddress( "nPrimary", &nPrimary );
  dijetTree->SetBranchAddress( "nGlobal", &nGlobal );
  dijetTree->SetBranchAddress( "nVertices", &nVertices );
  dijetTree->SetBranchAddress( "refMult", &refMult );
  dijetTree->SetBranchAddress( "gRefMult", &gRefMult );
  dijetTree->SetBranchAddress( "Vz", &Vz );
  dijetTree->SetBranchAddress( "leadPt", &leadPt );
  dijetTree->SetBranchAddress( "leadEta", &leadEta );
  dijetTree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );
  dijetTree->SetBranchAddress( "BbcAdcSumEastOuter", &BbcAdcSumEastOuter );
  dijetTree->SetBranchAddress( "BbcAdcSumWest", &BbcAdcSumWest );
  dijetTree->SetBranchAddress( "BbcAdcSumWestOuter", &BbcAdcSumWestOuter );
  dijetTree->SetBranchAddress( "leadPhi", &leadPhi );
  dijetTree->SetBranchAddress( "chgEastRho", &chgEastRho );
  dijetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  dijetTree->SetBranchAddress( "chgWestRho", &chgWestRho );
  dijetTree->SetBranchAddress( "neuEastRho", &neuEastRho );
  dijetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  dijetTree->SetBranchAddress( "neuWestRho", &neuWestRho );
  dijetTree->SetBranchAddress( "leadArea", &leadArea );
  dijetTree->SetBranchAddress( "recoArea", &recoArea );
  dijetTree->SetBranchAddress( "recoPt", &recoPt );
  dijetTree->SetBranchAddress( "recoEta", &recoEta );
  dijetTree->SetBranchAddress( "recoPhi", &recoPhi );
  dijetTree->SetBranchAddress( "nHTtrig", &nHTtrig );

  TH1D *hBBCE[nEtaBins][nEtaBins];
  
  for (int i=0; i<nEtaBins; ++i) {
    for (int j=0; j<nEtaBins; ++j) {
      name = "hBBCE_" + emw[i] + "lead_" + emw[j] + "reco";
      title = leadEtaBinString[i] + "         " + recoEtaBinString[j] +";iBBCE Sum";
      hBBCE[i][j] = new TH1D(name, title, 12,3000,70000 );
    }
  }

  for ( int i=0; i<dijetTree->GetEntries(); ++i ) {

    dijetTree->GetEntry(i);

    int ival = 99; int jval = 99;
    for (int ebin=0; ebin<nEtaBins; ++ebin) {
      if ( leadEta>=etaLo[ebin] && leadEta<=etaHi[ebin]) {ival = ebin;}
      if ( recoEta>=etaLo[ebin] && recoEta<=etaHi[ebin]) {jval = ebin;}
    }
    
    hBBCE[ival][jval]->Fill(BbcAdcSumEast);
    
  }
  
  for (int i=0; i<nEtaBins; ++i) {
    for (int j=0; j<nEtaBins; ++j) { hBBCE[i][j]->Scale(1./hBBCE[i][j]->GetEntries()); }
  }
  
  TCanvas * c1 = new TCanvas( "c1" , "" ,700 ,500 );              // CANVAS 1
  c1->Divide(3,3);

  for (int i=0; i<nEtaBins; ++i) {
    for (int j=0; j<nEtaBins; ++j) {
      int panel = (3*j) + i + 1;
      c1->cd(panel);
      hBBCE[i][j]->Draw();
    }
  }

  c1->SaveAs("plots/UE/BBCEsum_diffEta.pdf","PDF");

  TCanvas * c2 = new TCanvas( "c2" , "" ,700 ,500 );              // CANVAS 2
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.05);
  
  TH1D *hLeadEta[nEtaBins];
  TH1D *hRecoEta[nEtaBins];

  for (int i=0; i<nEtaBins; ++i) {
    name = "hLeadEta_" + emw[i] + "RecoEta";
    hLeadEta[i] = new TH1D(name,";#LT iBBCEsum#GT",3,0,3);
    hLeadEta[i]->SetLineColor( etaColor[i] );
    hLeadEta[i]->SetMarkerColor( etaColor[i] );
    hLeadEta[i]->SetMarkerStyle( etaMarker[i] );

    name = "hRecoEta_" + emw[i] + "LeadEta";
    hRecoEta[i] = new TH1D(name,";#LT iBBCEsum#GT",3,0,3);
    hRecoEta[i]->SetLineColor( etaColor[i] );
    hRecoEta[i]->SetMarkerColor( etaColor[i] );
    hRecoEta[i]->SetMarkerStyle( etaMarker[i] );
    
    for (int j=0; j<nEtaBins; ++j) {
      name = emw[j] + " reco jet";
      hLeadEta[i]->GetXaxis()->SetBinLabel(j+1,name);
      hLeadEta[i]->SetBinContent(j+1,hBBCE[i][j]->GetMean(1));
      hLeadEta[i]->SetBinError(j+1,hBBCE[i][j]->GetMeanError(1));
      
      name = emw[j] + " lead jet";
      hRecoEta[i]->GetXaxis()->SetBinLabel(j+1,name);
      hRecoEta[i]->SetBinContent(j+1,hBBCE[j][i]->GetMean(1));
      hRecoEta[i]->SetBinError(j+1,hBBCE[j][i]->GetMeanError(1));
    }
  }


  TH2D *hScale = new TH2D("hScale",";;#LT iBBCEsum#GT",3,0,3,20,21000,25000);
  hScale->SetStats(0);
  hScale->GetYaxis()->SetTitleOffset(1.5);
  hScale->SetLineColorAlpha(1,0.0);
  hScale->SetMarkerColorAlpha(1,0.0);
    
  for (int j=0; j<nEtaBins; ++j) {
    name = emw[j] + " reco jet";
    hRecoEta[j]->SetName(name);
    hScale->GetXaxis()->SetBinLabel(j+1,name);
    name = emw[j] + " lead jet";
    hLeadEta[j]->SetName(name);
  }

  hScale->Draw();
  for (int i=0; i<nEtaBins; ++i) {
    hLeadEta[i]->Draw("SAME");
  }
  hScale->SetName(" ");
  TLegend *leadleg = (TLegend*) c2->BuildLegend();
  leadleg->SetLineColorAlpha(1,0.0);
  c2->SaveAs("plots/UE/BBCEsum_recoEtaByLeadEta.pdf","PDF");
  
  for (int j=0; j<nEtaBins; ++j) {
    name = emw[j] + " lead jet";
    hScale->GetXaxis()->SetBinLabel(j+1,name);
  }
  hScale->SetTitle(";;#LT iBBCEsum#GT");
  hScale->Draw();
  for (int i=0; i<nEtaBins; ++i) {
    hRecoEta[i]->Draw("SAME");
  }
  hScale->SetName(" ");
  TLegend *recoleg = (TLegend*) c2->BuildLegend();
  recoleg->SetLineColorAlpha(1,0.0);
  c2->SaveAs("plots/UE/BBCEsum_leadEtaByRecoEta.pdf","PDF");
  
}
