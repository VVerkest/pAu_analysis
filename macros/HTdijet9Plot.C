
void HTdijet9Plot(){

  TFile* inFile = new TFile( "out/HTdijets/HiEA/pAu_2015_HTdijets.root", "READ" );
  TString saveName = "plots/HTdijets/HiEAdijet9plot.pdf";
  
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
  const int marker[nChgBins] = { 20, 20, 21 };

  gStyle->SetOptStat(0);
  TH2D *hRhoByEta[nPtBins][nEtaBins][nChgBins];
  TH2D *sRhoByEta = new TH2D("sRhoByEta","", 3,-1.5,1.5, 10,0.2,1.2);
  sRhoByEta->GetXaxis()->SetLabelSize(0);
  sRhoByEta->GetYaxis()->SetLabelSize(0.06);
  sRhoByEta->GetYaxis()->SetNdivisions(12);
  
  TH1D *pRhoByEta[nPtBins][nEtaBins][nChgBins];

  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {

	TString name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	hRhoByEta[p][e][c] = (TH2D*)inFile->Get(name);
      }
    }
  }

  double stdev;
  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  c0->SetTopMargin(0.4);
  TPaveText *cTitle = new TPaveText(0.345843,.881306,0.655712,.980712,"NB");
  cTitle->AddText("Underlying Event");
  cTitle->SetFillStyle(0);
  cTitle->SetLineWidth(0);
  cTitle->SetTextAlign(21);
  cTitle->Draw();
    
  c0->Divide(nEtaBins,nPtBins,0,0);

  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {

      int dir = 1+e+(3*p);
      c0->cd(dir);
      sRhoByEta->Draw();
      for ( int c=0; c<3; ++c ) {

	stdev = hRhoByEta[p][e][c]->GetMeanError(2);
	
	TString name = "pRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	hRhoByEta[p][e][c]->Scale(1./hRhoByEta[p][e][c]->GetEntries());
	pRhoByEta[p][e][c] = (TH1D*) hRhoByEta[p][e][c]->ProfileX(name,1,-1,"i");
	//pRhoByEta[p][e][c]->SetError( (const Double_t) stdev );
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridy();
	// pRhoByEta[p][e][c]->SetLineColorAlpha(0,0.000001);
	pRhoByEta[p][e][c]->SetLineWidth(0);
	pRhoByEta[p][e][c]->SetMarkerSize(2);
	pRhoByEta[p][e][c]->SetMarkerStyle( marker[c] );
	pRhoByEta[p][e][c]->GetXaxis()->SetLabelSize(0);
	pRhoByEta[p][e][c]->GetYaxis()->SetLabelSize(0.06);
	pRhoByEta[p][e][c]->GetYaxis()->SetNdivisions(10);
	pRhoByEta[p][e][c]->Draw("SAME");


      }
    }
  }

  c0->SaveAs( saveName,"PDF");
  
}
