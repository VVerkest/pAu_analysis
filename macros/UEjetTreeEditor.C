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
  const double BBCEsumLo[nEA] = { 4107, 11503, 28537 };   // LO: 4107-11503;  HI: 28537+
  const double BBCEsumHi[nEA] = { 11503, 28537, 64000 };
  
  int jeval, bgeval, pval, eaval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, rho;
  
  TString fileName = "out/UE/pAuHTjetUE.root";
  TFile* inFile = new TFile( fileName, "UPDATE" );

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;  //  Tree variables
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, eastRho, midRho, westRho, leadPtCorrected;

  jetTree->SetBranchAddress( "RunID", &RunID );      	       	jetTree->SetBranchAddress( "EventID", &EventID );					jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );       	jetTree->SetBranchAddress( "nGlobal", &nGlobal );					jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );		jetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );	       		jetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );	       	jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );		jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );	//jetTree->Branch( "leadArea", &leadArea );					jetTree->Branch( "leadPtCorrected", &leadPtCorrected );

  auto leadAreaBranch = jetTree->Branch("leadArea", &leadArea, "leadArea/D");
  auto leadPtCorrectedBranch = jetTree->Branch("leadPtCorrected", &leadPtCorrected, "leadPtCorrected/D");
  
  int nEntries = jetTree->GetEntries();

  // double rhoByEta[nEtaBins];

  // TH1D *hRhoDist[nEA][nEtaBins][nPtBins][nEtaBins];
  
  // for ( int ea=0; ea<nEA; ++ea ) {
  //   for ( int bge=0; bge<nEtaBins; ++bge ) {
  //     for ( int p=0; p<nPtBins; ++p ) {
  // 	for ( int je=0; je<nEtaBins; ++je ) {
  // 	  name = "hRho_" + EAstring[ea] + jetEtaBinName[je] + ptBinName[p] + etaBinName[bge];
  // 	  title = "Underlying Event" + EAstring[ea] + jetEtaBinName[je] + ptBinName[p] + etaBinName[bge] +";#rho (GeV)";
  // 	  hRhoDist[ea][bge][p][je]= new TH1D(name,title,120,0,30);
  // 	}
  //     }
  //   }
  // }
  
  // for ( int i=0; i<nEntries; ++i ) {

  //   jetTree->GetEntry(i);

  //   rhoByEta[0] = chgEastRho + neuEastRho;
  //   rhoByEta[1] = chgMidRho + neuMidRho;
  //   rhoByEta[2] = chgWestRho + neuWestRho;
    
  //   pval = 99;    jeval = 99;    eaval = 99;
    
  //   for ( int ea=0; ea<3; ++ea ) {
  //     if ( BbcAdcEastSum > BBCEsumLo[ea]  &&  BbcAdcEastSum < BBCEsumHi[ea] ) { eaval = ea; }
  //   }    
  //   for ( int p=0; p<3; ++p ) {
  //     if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
  //   }
  //   for ( int je=0; je<3; ++je ) {
  //     if ( leadEta >= etaLo[je]  &&  leadEta <= etaHi[je] ) { jeval = je; }
  //   }
  //   if ( eaval==99 ) { continue; }
  //   if ( pval==99 || jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<pval<<endl<<jeval<<endl<<bgeval<<endl<<leadEta<<endl<<endl; }
    
  //   for ( int bge=0; bge<nEtaBins; ++bge ) {
  //     hRhoDist[eaval][bge][pval][jeval]->Fill( rhoByEta[bge] );
  //   }
  // }


  // double avgRho[nEA][nEtaBins][nPtBins][nEtaBins];
  
  // for ( int ea=0; ea<nEA; ++ea ) {
  //   for ( int bge=0; bge<nEtaBins; ++bge ) {
  //     for ( int p=0; p<nPtBins; ++p ) {
  // 	for ( int je=0; je<nEtaBins; ++je ) {
  // 	  avgRho[ea][bge][p][je] = hRhoDist[ea][bge][p][je]->GetMean(1);
  // 	}
  //     }
  //   }
  // }


  jetTree->Draw("(chgEastRho+neuEastRho):BbcAdcEastSum>>hEastRhoByBBCEsum","BbcAdcEastSum>4107","COLZ");
  TH2D* hEastRhoByBBCEsum = (TH2D*)gDirectory->Get("hEastRhoByBBCEsum");
  TH1D* hEastRhoProfile = (TH1D*)hEastRhoByBBCEsum->ProfileX();
  int nEastProfBins = hEastRhoProfile->GetNbinsX();
  jetTree->Draw("(chgMidRho+neuMidRho):BbcAdcEastSum>>hMidRhoByBBCEsum","BbcAdcEastSum>4107","COLZ");
  TH2D* hMidRhoByBBCEsum = (TH2D*)gDirectory->Get("hMidRhoByBBCEsum");
  TH1D* hMidRhoProfile = (TH1D*)hMidRhoByBBCEsum->ProfileX();
  int nMidProfBins = hMidRhoProfile->GetNbinsX();
  jetTree->Draw("(chgWestRho+neuWestRho):BbcAdcEastSum>>hWestRhoByBBCEsum","BbcAdcEastSum>4107","COLZ");
  TH2D* hWestRhoByBBCEsum = (TH2D*)gDirectory->Get("hWestRhoByBBCEsum");
  TH1D* hWestRhoProfile = (TH1D*)hWestRhoByBBCEsum->ProfileX();
  int nWestProfBins = hWestRhoProfile->GetNbinsX();

  for ( int i=0; i<nEntries; ++i ) {
    jetTree->GetEntry(i);

    leadArea = pi*(0.4)*(0.4);
    
    // pval=99; jeval=99;
    // for ( int ea=0; ea<nEA; ++ea ) {
    //   if ( BbcAdcEastSum >= BBCEsumLo[ea]  &&  BbcAdcEastSum <= BBCEsumHi[ea] ) { eaval = ea; }
    // }    
    // for ( int p=0; p<nPtBins; ++p ) {
    //   if ( leadPt >= ptLo[p]  &&  leadPt <= ptHi[p] ) { pval = p; }
    // }
    // for ( int je=0; je<nEtaBins; ++je ) {
    //   if ( leadEta >= etaLo[je]  &&  leadEta <= etaHi[je] ) { jeval = je; }
    // }
    // if ( pval==99 || jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl; }
    // if ( eaval==99 ) { continue; }

    double rhoValue;

    if ( leadEta >= etaLo[0]  &&  leadEta <= etaHi[0] ) {  //EAST
      for ( int j=0; j<nEastProfBins; ++j ) {
	double BBCElo = hEastRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hEastRhoProfile->GetBinWidth(j);
	if ( BbcAdcEastSum >= BBCElo &&  BbcAdcEastSum <= BBCEhi ) { rhoValue = hEastRhoProfile->GetBinContent(j); }
      }
    }
    else if ( leadEta >= etaLo[1]  &&  leadEta <= etaHi[1] ) {  //MID
      for ( int j=0; j<nMidProfBins; ++j ) {
	double BBCElo = hMidRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hMidRhoProfile->GetBinWidth(j);
	if ( BbcAdcEastSum >= BBCElo &&  BbcAdcEastSum <= BBCEhi ) { rhoValue = hMidRhoProfile->GetBinContent(j); }
      }
    }
    else if ( leadEta >= etaLo[2]  &&  leadEta <= etaHi[2] ) {  //WEST
      for ( int j=0; j<nWestProfBins; ++j ) {
	double BBCElo = hWestRhoProfile->GetBinLowEdge(j);
	double BBCEhi = BBCElo + hWestRhoProfile->GetBinWidth(j);
	if ( BbcAdcEastSum >= BBCElo &&  BbcAdcEastSum <= BBCEhi ) { rhoValue = hWestRhoProfile->GetBinContent(j); }
      }
    }
    else { cerr<<"error with finding lead jet eta"<<endl; }

    cout<<rhoValue<<endl;
    
    leadPtCorrected = leadPt - leadArea*rhoValue;
    
    // leadPtCorrected = leadPt - leadArea*avgRho[eaval][jeval][pval][jeval];
    leadPtCorrectedBranch->Fill();
    leadAreaBranch->Fill();
  }

  jetTree->Write("", TObject::kOverwrite);
}
