
void HTjet9Plot(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  gStyle->SetErrorX(0.0001);
  
  TFile* inFile = new TFile( "out/HTdijets/LoEA/pAu_2015_HTjets.root", "READ" );
  TString saveName = "plots/HTjets/LoEAjet9plot.pdf";
  TString title = "Low EA Jet Underlying Event";
  
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
  TH2D *sRhoByEta = new TH2D("sRhoByEta","", 3,-1.5,1.5, 10,0,2.0);
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
  cTitle->AddText( title );
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
	//pRhoByEta[p][e][c]->SetLineWidth(0);
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

/*

  void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2, Float_t lMargin = 0.15, Float_t rMargin = 0.05, Float_t bMargin = 0.15, Float_t tMargin = 0.05) {
  if (!C) return;

  // Setup Pad layout:
  Float_t vSpacing = 0.0;
  Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

  Float_t hSpacing = 0.0;
  Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

  Float_t vposd,vposu,vmard,vmaru,vfactor;
  Float_t hposl,hposr,hmarl,hmarr,hfactor;

  for (Int_t i=0;i<Nx;i++) {

  if (i==0) {
  hposl = 0.0;
  hposr = lMargin + hStep;
  hfactor = hposr-hposl;
  hmarl = lMargin / hfactor;
  hmarr = 0.0;
  } else if (i == Nx-1) {
  hposl = hposr + hSpacing;
  hposr = hposl + hStep + rMargin;
  hfactor = hposr-hposl;
  hmarl = 0.0;
  hmarr = rMargin / (hposr-hposl);
  } else {
  hposl = hposr + hSpacing;
  hposr = hposl + hStep;
  hfactor = hposr-hposl;
  hmarl = 0.0;
  hmarr = 0.0;
  }

  for (Int_t j=0;j<Ny;j++) {

  if (j==0) {
  vposd = 0.0;
  vposu = bMargin + vStep;
  vfactor = vposu-vposd;
  vmard = bMargin / vfactor;
  vmaru = 0.0;
  } else if (j == Ny-1) {
  vposd = vposu + vSpacing;
  vposu = vposd + vStep + tMargin;
  vfactor = vposu-vposd;
  vmard = 0.0;
  vmaru = tMargin / (vposu-vposd);
  } else {
  vposd = vposu + vSpacing;
  vposu = vposd + vStep;
  vfactor = vposu-vposd;
  vmard = 0.0;
  vmaru = 0.0;
  }

  C->cd(0);

  char name[16];
  sprintf(name,"pad_%i_%i",i,j);
  TPad *pad = (TPad*) gROOT->FindObject(name);
  if (pad) delete pad;
  pad = new TPad(name,"",hposl,vposd,hposr,vposu);
  pad->SetLeftMargin(hmarl);
  pad->SetRightMargin(hmarr);
  pad->SetBottomMargin(vmard);
  pad->SetTopMargin(vmaru);

  pad->SetFrameBorderMode(0);
  pad->SetBorderMode(0);
  pad->SetBorderSize(0);

  pad->Draw();
  }
  }
  }
*/