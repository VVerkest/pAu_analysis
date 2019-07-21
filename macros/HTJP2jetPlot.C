
void HTJP2jetPlot(){

  TFile* inFile = new TFile( "out/HTJP2jets/pAu_2015_200_HTJP2jets.root", "READ" );

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-1.0<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<1.0" };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 24 };

  gStyle->SetOptStat(0);
  TH2D *hRhoByEta[nPtBins][nEtaBins][nChgBins];


  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {

	TString name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	hRhoByEta[p][e][c] = (TH2D*)inFile->Get(name);
      }
    }
  }


  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  c0->Divide(3,3);

    for ( int p=0; p<3; ++p ) {
      for ( int e=0; e<3; ++e ) {

	int dir = 1+e+(3*p);
	c0->cd(dir);

	for ( int c=0; c<3; ++c ) {

	  hRhoByEta[p][e][c]->ProfileX()->Draw("SAME");


	}
      }
    }

}
