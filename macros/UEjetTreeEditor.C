// Veronica Verkest
// December 9, 2019

void UEjetTreeEditor(){

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
  const double BBCEsumLo[nEA] = { 3559.12, 10126.1, 26718.1 };   // LO: 3559.12-10126.1;  HI: 26718.1+
  const double BBCEsumHi[nEA] = { 10126.1, 26718.1, 64000.0 };  // 0, 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho;
  
  // TString fileName = "out/UE/pAuHTjetUE.root";
  TString fileName = "out/UE/pAuHTjetUE_1HTtrig_inLead.root";
  TFile* inFile = new TFile( fileName, "UPDATE" );

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;  //  Tree variables
  double Vz, BbcAdcSumEast, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, chgEastRho_te, chgMidRho_te, chgWestRho_te, neuEastRho, neuMidRho, neuWestRho, leadArea, eastRho, midRho, westRho, leadPtCorrected, rho, rho_te;

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

  auto leadPtCorrectedBranch = jetTree->Branch("leadPtCorrected", &leadPtCorrected, "leadPtCorrected/D");
  auto rhoBranch = jetTree->Branch("rho", &rho, "rho/D");
  auto rho_teBranch = jetTree->Branch("rho_te", &rho_te, "rho_te/D");
  
  int nEntries = jetTree->GetEntries();

  jetTree->Draw("(chgEastRho+neuEastRho):BbcAdcSumEast>>hEastRhoByBBCEsum","BbcAdcSumEast>3559.12  &&  BbcAdcSumEast<64000  &&  leadEta<-0.3","COLZ");
  TH2D* hEastRhoByBBCEsum = (TH2D*)gDirectory->Get("hEastRhoByBBCEsum");
  TH1D* hEastRhoProfile = (TH1D*)hEastRhoByBBCEsum->ProfileX();
  int nEastProfBins = hEastRhoProfile->GetNbinsX();
  jetTree->Draw("(chgMidRho+neuMidRho):BbcAdcSumEast>>hMidRhoByBBCEsum","BbcAdcSumEast>3559.12  &&  BbcAdcSumEast<64000  &&  leadEta>=-0.3  &&  leadEta<=0.3","COLZ");
  TH2D* hMidRhoByBBCEsum = (TH2D*)gDirectory->Get("hMidRhoByBBCEsum");
  TH1D* hMidRhoProfile = (TH1D*)hMidRhoByBBCEsum->ProfileX();
  int nMidProfBins = hMidRhoProfile->GetNbinsX();
  jetTree->Draw("(chgWestRho+neuWestRho):BbcAdcSumEast>>hWestRhoByBBCEsum","BbcAdcSumEast>3559.12  &&  BbcAdcSumEast<64000  &&  leadEta>0.3","COLZ");
  TH2D* hWestRhoByBBCEsum = (TH2D*)gDirectory->Get("hWestRhoByBBCEsum");
  TH1D* hWestRhoProfile = (TH1D*)hWestRhoByBBCEsum->ProfileX();
  int nWestProfBins = hWestRhoProfile->GetNbinsX();

  for ( int i=0; i<nEntries; ++i ) {
    jetTree->GetEntry(i);

    double rhoValue;

    if ( leadEta >= etaLo[0]  &&  leadEta <= etaHi[0] ) {  //EAST
      for ( int j=0; j<nEastProfBins; ++j ) {
	double BBCElo = hEastRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hEastRhoProfile->GetBinWidth(j);
	if ( BbcAdcSumEast >= BBCElo &&  BbcAdcSumEast <= BBCEhi ) { rhoValue = hEastRhoProfile->GetBinContent(j); }
      }
    }
    else if ( leadEta >= etaLo[1]  &&  leadEta <= etaHi[1] ) {  //MID
      for ( int j=0; j<nMidProfBins; ++j ) {
	double BBCElo = hMidRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hMidRhoProfile->GetBinWidth(j);
	if ( BbcAdcSumEast >= BBCElo &&  BbcAdcSumEast <= BBCEhi ) { rhoValue = hMidRhoProfile->GetBinContent(j); }
      }
    }
    else if ( leadEta >= etaLo[2]  &&  leadEta <= etaHi[2] ) {  //WEST
      for ( int j=0; j<nWestProfBins; ++j ) {
	double BBCElo = hWestRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hWestRhoProfile->GetBinWidth(j);
	if ( BbcAdcSumEast >= BBCElo &&  BbcAdcSumEast <= BBCEhi ) { rhoValue = hWestRhoProfile->GetBinContent(j); }
      }
    }
    else { cerr<<"error with finding lead jet eta"<<endl; }
    
    leadPtCorrected = leadPt - leadArea*rhoValue;
    rho = (chgEastRho+neuEastRho+chgMidRho+neuMidRho+chgWestRho+neuWestRho)/3;
    rho_te = (chgEastRho_te+neuEastRho+chgMidRho_te+neuMidRho+chgWestRho_te+neuWestRho)/3;
    
    leadPtCorrectedBranch->Fill();
    rhoBranch->Fill();
    rho_teBranch->Fill();
  }

  jetTree->Write("", TObject::kOverwrite);
}
