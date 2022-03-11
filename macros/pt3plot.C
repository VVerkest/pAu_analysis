// Veronica Verkest
// May 21, 2020

void pt3plot(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  const double AREA = 4.*(pi/3.);

  TString name, title;
  
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

  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 23, 22 };

  int nChg[nEtaBins][nEAbins][nPtBins];
  double ptSum[nEtaBins][nEAbins][nPtBins];
  double te_ptSum[nEtaBins][nEAbins][nPtBins];
  //  [nEAbins][nPtBins][nEtaBins]

  TFile *inFile = new TFile( "src/differentialUE.root", "READ" );

  TH2D* hChgPt2D[nEAbins][nPtBins];
  TH2D* hChgPt2D_te[nEAbins][nPtBins];

  TH1D* hChgPt[nEAbins][nPtBins];
  TH1D* hChgPt_te[nEAbins][nPtBins];

  TH1D *hNchg[nEAbins][nPtBins];
  TH1D *hNneu[nEAbins][nPtBins];
  
  for ( int a=0; a<nEAbins; ++a) {
    for ( int p=0; p<nPtBins; ++p) {

      name = "hChgUEpt"; name += ptBinName[p]; name += EAbinString[a];
      hChgPt2D[a][p] = (TH2D*)inFile->Get( name );
      
      name += "_te";
      hChgPt2D_te[a][p] = (TH2D*)inFile->Get( name );

      //name = "hNchg_" + EAbinString[a] + ptBinName[p];

      for ( int e=0; e<nEtaBins; ++e) {

	hChgPt2D[a][p]->GetYaxis()->SetRangeUser( etaLo[e], etaHi[e] );
	ptSum[e][a][p] = hChgPt2D[a][p]->GetMean(1);
	cout<<ptSum[e][a][p]<<"   ";
	
	//nChg[e][a][p] = ;
		
      }      
    }
  }



  // TH2D *hscale = new TH2D("hscale","",3,-1.5,1.5,10,0.0,1.0);
  // hscale->SetStats(0);
  // hscale->GetXaxis()->SetLabelSize(0);
  // hscale->GetYaxis()->SetLabelSize(0.08);

  
  
}
