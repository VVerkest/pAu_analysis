
void EAratioPlots(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  gStyle->SetErrorX(0.0001);

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

  const int c = 3;
  const TString BackgroundChargeBias[c] = { "_chgBG", "_neuBG", "_allBG" };
  TString name;
  
  const int EA = 2;
  const int DJ = 2;
  
  TFile* File[EA][DJ];
  TString activity[EA] = {"LoEA", "HiEA"};
  TString jetdijet[DJ] = {"jets", "dijets"};
  TString fileName[EA][DJ];
  TString us = "_";

  int marker[EA][DJ] ={ { 20, 24 } , { 21, 25 } };  // LoEA: circle,  HiEA: square
  EColor color[DJ] = { kBlack, kRed };

  gStyle->SetOptStat(0);
  TH2D *sRhoByEta = new TH2D("sRhoByEta","", 3,-1.5,1.5, 10,0,2.0);
  sRhoByEta->GetXaxis()->SetLabelSize(0);
  sRhoByEta->GetYaxis()->SetLabelSize(0.06);
  sRhoByEta->GetYaxis()->SetNdivisions(12);
  TH2D *hRhoByEta[nPtBins][nEtaBins][EA][DJ];
  TH1D *pRhoByEta[nPtBins][nEtaBins][EA][DJ];


  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1000 ,700 );
  c0->SetTopMargin(0.4);
  TPaveText *cTitle = new TPaveText(0.345843,.881306,0.655712,.980712,"NB");
  cTitle->AddText( "Event Activity" );
  cTitle->SetFillStyle(0);
  cTitle->SetLineWidth(0);
  cTitle->SetTextAlign(21);
  cTitle->Draw();
    
  c0->Divide(nEtaBins,nPtBins,0,0);

  
  for (int ea=0; ea<EA; ++ea) {
    for (int dj=0; dj<DJ; ++dj) {
      fileName[ea][dj] = "out/HTdijets/" + activity[ea] + "/pAu_2015_HT" + jetdijet[dj] + ".root";
      cout << fileName[ea][dj] <<endl;
      File[ea][dj] = new TFile( fileName[ea][dj], "READ" );
      cout << File[ea][dj] << endl;
      for ( int p=0; p<3; ++p ) {
	for ( int e=0; e<3; ++e ) {
	
	  name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	  hRhoByEta[p][e][ea][dj] = (TH2D*)File[ea][dj]->Get(name);
	  cout << hRhoByEta[p][e][ea][dj] << endl;
	}
      }
      
      File[ea][dj]->Close();
    }
  }


  for (int ea=0; ea<EA; ++ea) {
    for (int dj=0; dj<DJ; ++dj) {
      
      for ( int p=0; p<3; ++p ) {
	for ( int e=0; e<3; ++e ) {

	  double entries = hRhoByEta[p][e][ea][dj]->GetEntries();
	  hRhoByEta[p][e][ea][dj]->Scale(1./entries);
	  hRhoByEta[p][e][ea][dj]->SetMarkerStyle( marker[ea][dj] );
	  hRhoByEta[p][e][ea][dj]->SetMarkerColor( color[dj] );
	  hRhoByEta[p][e][ea][dj]->SetLineColor( color[dj] );

	  name = "pRho" + ptBinName[p] + etaBinName[e] + us + activity[ea] + us + jetdijet[dj] + ".root";
	  pRhoByEta[p][e][ea][dj] = (TH1D*) hRhoByEta[p][e][ea][dj]->ProfileX(name,1,-1,"i");
	  
	  gPad->SetTickx();
	  gPad->SetTicky();
	  gPad->SetGridy();
	  // pRhoByEta[p][e][ea][dj]->SetLineWidth(0);
	  pRhoByEta[p][e][ea][dj]->SetMarkerSize(2);
	  pRhoByEta[p][e][ea][dj]->GetXaxis()->SetLabelSize(0);
	  pRhoByEta[p][e][ea][dj]->GetYaxis()->SetLabelSize(0.06);
	  pRhoByEta[p][e][ea][dj]->GetYaxis()->SetNdivisions(10);
	  pRhoByEta[p][e][ea][dj]->Draw("SAME");


	}
      }
      
    }
  }
  



  

  

}
