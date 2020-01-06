// Veronica Verkest
// December 17, 2019

void UEdijetTreeEditor(){

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
  
  TString fileName = "out/UE/pAuHTdijetUE.root";
  TFile* inFile = new TFile( fileName, "UPDATE" );

  TTree *dijetTree = (TTree*) inFile->Get("HTdijetTree");

  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;  //  Tree variables
  double Vz, BbcAdcSumEast, leadPt, leadEta, leadPhi, recoPt, recoEta, recoPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, recoArea,
    eastRho, midRho, westRho, leadPtCorrected, recoPtCorrected;

  dijetTree->SetBranchAddress( "RunID", &RunID );      	       	dijetTree->SetBranchAddress( "EventID", &EventID );		       		dijetTree->SetBranchAddress( "nTowers", &nTowers );
  dijetTree->SetBranchAddress( "nPrimary", &nPrimary );       	dijetTree->SetBranchAddress( "nGlobal", &nGlobal );	       			dijetTree->SetBranchAddress( "nVertices", &nVertices );
  dijetTree->SetBranchAddress( "refMult", &refMult );		dijetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		dijetTree->SetBranchAddress( "Vz", &Vz );
  dijetTree->SetBranchAddress( "leadPt", &leadPt );	       		dijetTree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );	dijetTree->SetBranchAddress( "leadEta", &leadEta );
  dijetTree->SetBranchAddress( "leadPhi", &leadPhi );	       	dijetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		dijetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  dijetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	dijetTree->SetBranchAddress( "neuEastRho", &neuEastRho );		dijetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  dijetTree->SetBranchAddress( "neuWestRho", &neuWestRho );	dijetTree->SetBranchAddress( "leadArea", &leadArea );			dijetTree->SetBranchAddress( "recoArea", &recoArea );
  dijetTree->SetBranchAddress( "recoPt", &recoPt );    			dijetTree->SetBranchAddress( "recoEta", &recoEta );			dijetTree->SetBranchAddress( "recoPhi", &recoPhi );	       	
  
  auto leadPtCorrectedBranch = dijetTree->Branch("leadPtCorrected", &leadPtCorrected, "leadPtCorrected/D");
  auto recoPtCorrectedBranch = dijetTree->Branch("recoPtCorrected", &recoPtCorrected, "recoPtCorrected/D");
  
  int nEntries = dijetTree->GetEntries();
  
  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho+neuEastRho)", "(chgMidRho+neuMidRho)", "(chgWestRho+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString correctedPtSelection[nPtBins] = {"leadPtCorrected>10.0&&leadPtCorrected<15.0", "leadPtCorrected>=15.0&&leadPtCorrected<=20.0", "leadPtCorrected>20.0&&leadPtCorrected<30.0"};
  TString etaSelection[nEtaBins] = { "leadEta<-0.3", "leadEta>=-0.3 && leadEta<=0.3", "leadEta>0.3" };

  TH2D *hRho[nEtaBins][nEtaBins];
  TH1D *hRhoProfile[nEtaBins][nEtaBins];
  TString hname[nEtaBins][nEtaBins];
  int nBins[nEtaBins][nEtaBins];

  for ( int je=0; je<nEtaBins; ++je ) {   //  create rhoVsBBCEsum profiles!
    for ( int bge=0; bge<nEtaBins; ++bge ) {
      hname[je][bge] = "hRho" + eastmidwest[bge] + jetEtaBinName[je];
      TString drawString = rhoVal[bge] + ":BbcAdcSumEast>>" + hname[je][bge];
      TString selectionString = "BbcAdcSumEast>4107  &&  BbcAdcSumEast<64000  &&  " + etaSelection[je];
      dijetTree->Draw( drawString, selectionString, "COLZ" );
      hRho[je][bge] = (TH2D*)gDirectory->Get( hname[je][bge] );
      hRhoProfile[je][bge] = (TH1D*)hRho[je][bge]->ProfileX();
      nBins[je][bge] = hRho[je][bge]->GetNbinsX();
    }
  }


  for ( int i=0; i<nEntries; ++i ) {    // loop over tree and add corrected pt values
    dijetTree->GetEntry(i);

    double leadRhoValue = 99;    int le = 99;   // lead eta value
    double recoRhoValue = 99;    int re = 99;   // reco eta value

    for ( int e=0; e<nEtaBins; ++e ) {
      if ( leadEta >= etaLo[e]  &&  leadEta <= etaHi[e] ) { le = e; }
      if ( recoEta >= etaLo[e]  &&  recoEta <= etaHi[e] ) { re = e; }
    }
    if ( le == 99  ||  re == 99 ) { cerr<<"error with finding lead or reco jet eta"<<endl; }

    for ( int j=0; j<nBins[le][le]; ++j ) { //  lead jet
      double BBCElo = hRhoProfile[le][le]->GetBinLowEdge(j);
      double BBCEhi = BBCElo + hRhoProfile[le][le]->GetBinWidth(j);
      if ( BbcAdcSumEast >= BBCElo &&  BbcAdcSumEast <= BBCEhi ) { leadRhoValue = hRhoProfile[le][le]->GetBinContent(j); }
    }
    
    for ( int j=0; j<nBins[le][re]; ++j ) { //  recoil jet
      double BBCElo = hRhoProfile[le][re]->GetBinLowEdge(j);
      double BBCEhi = BBCElo + hRhoProfile[le][re]->GetBinWidth(j);
      if ( BbcAdcSumEast >= BBCElo &&  BbcAdcSumEast <= BBCEhi ) { recoRhoValue = hRhoProfile[le][re]->GetBinContent(j); }
    }    
  
    leadPtCorrected = leadPt - leadArea*leadRhoValue;
    recoPtCorrected = recoPt - recoArea*recoRhoValue;
    
    leadPtCorrectedBranch->Fill();
    recoPtCorrectedBranch->Fill();
  }
  

  dijetTree->Write("", TObject::kOverwrite);   // only write new version of the tree!
}
