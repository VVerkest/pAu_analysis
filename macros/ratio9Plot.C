
void ratio9Plot(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  gStyle->SetErrorX(0.0001);

  // TFile* File1 = new TFile( "out/HTdijets/LoEA/pAu_2015_HTdijets.root", "READ" );
  // TFile* File2 = new TFile( "out/HTdijets/LoEA/pAu_2015_HTjets.root", "READ" );
  // TString ratioTitle = "Low EA #rho_{dijet}/#rho_{jet}";
  // TString saveName = "plots/HTdijets/ratio9plot__LoEA.pdf";

  TFile* File1 = new TFile( "out/HTdijets/LoEA/pAu_2015_HTdijets.root", "READ" );
  TFile* File2 = new TFile( "out/HTdijets/HiEA/pAu_2015_HTdijets.root", "READ" );
  TString ratioTitle = "#rho_{dijet}^{lo EA} / #rho_{dijet}^{hi EA}";
  TString saveName = "plots/HTdijets/ratio9plot__dijetLoHiEA.pdf";
  
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

  TString name;
  
  gStyle->SetOptStat(0);
  TH1D *pRhoRatio[nPtBins][nEtaBins][nChgBins];
  TH2D *sRhoRatio = new TH2D("sRhoRatio","", 3,-1.5,1.5, 10,0.4,1.0);
  sRhoRatio->GetXaxis()->SetLabelSize(0);
  sRhoRatio->GetYaxis()->SetLabelSize(0.06);
  sRhoRatio->GetYaxis()->SetNdivisions(12);
  
  TH2D *file1RhoByEta[nPtBins][nEtaBins][nChgBins];
  TH2D *file2RhoByEta[nPtBins][nEtaBins][nChgBins];
  TH1D *hRhoRatio[nPtBins][nEtaBins][nChgBins];
  
  for ( int p=0; p<3; ++p ) {
    for ( int c=0; c<3; ++c ) {

      for ( int e=0; e<3; ++e ) {

	name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	file1RhoByEta[p][e][c] = (TH2D*)File1->Get(name);
	file2RhoByEta[p][e][c] = (TH2D*)File2->Get(name);
      
	name = "hRhoRatio" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	hRhoRatio[p][e][c] = new TH1D(name,"", 3,-1.5,1.5);

      }
    }
  }

  double ratio, error, var1, var2;
  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  c0->SetTopMargin(0.6);
  TPaveText *cTitle = new TPaveText(0.345843,.881306,0.655712,.980712,"NB");
  // cTitle->AddText("#rho_{dijet}/#rho_{jet}");
  cTitle->AddText( ratioTitle );
  cTitle->SetFillStyle(0);
  cTitle->SetLineWidth(0);
  cTitle->SetTextAlign(21);
  cTitle->Draw();
    
  c0->Divide(nEtaBins,nPtBins,0,0);

  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {

	int dir = 1+e+(3*p);
	c0->cd(dir);
	sRhoRatio->Draw();

      	for ( int i=1; i<4; ++i ) {
	
	  file1RhoByEta[p][e][c]->GetXaxis()->SetRange(i,i);
	  file2RhoByEta[p][e][c]->GetXaxis()->SetRange(i,i);
	  ratio = ( file1RhoByEta[p][e][c]->GetMean(2) )/( file2RhoByEta[p][e][c]->GetMean(2) );
	  hRhoRatio[p][e][c]->SetBinContent( i, ratio );
	  var1 = (file1RhoByEta[p][e][c]->GetMeanError(2))*(file1RhoByEta[p][e][c]->GetMeanError(2));
	  var2 = (file2RhoByEta[p][e][c]->GetMeanError(2))*(file2RhoByEta[p][e][c]->GetMeanError(2));
	  error = sqrt( var1 + var2 );
	  hRhoRatio[p][e][c]->SetBinContent( i, ratio );
	  //cout << i << "        " << ratio << endl;
	}
      }
    }
  }


  for ( int c=0; c<3; ++c ) {
    for ( int p=0; p<3; ++p ) {
      for ( int e=0; e<3; ++e ) {

	TString name = "pRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];

	// int i = e+1;
	
	// file1RhoByEta[p][e][c]->GetXaxis()->SetRange(i,i);
	// file2RhoByEta[p][e][c]->GetXaxis()->SetRange(i,i);
	// ratio = ( file1RhoByEta[p][e][c]->GetMean(2) )/( file2RhoByEta[p][e][c]->GetMean(2) );
	// hRhoRatio[p][e][c]->Fill( i, ratio );
	// cout << i << "        " << ratio << endl;

	//hRhoRatio[p][e][c]->SetError( (const Double_t) stdev );
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGridy();
	// hRhoRatio[p][e][c]->SetLineColorAlpha(0,0.000001);
	//hRhoRatio[p][e][c]->SetLineWidth(0);
	hRhoRatio[p][e][c]->SetMarkerSize(2);
	hRhoRatio[p][e][c]->SetMarkerStyle( marker[c] );
	hRhoRatio[p][e][c]->SetMarkerColor( color[c] );
	hRhoRatio[p][e][c]->GetXaxis()->SetLabelSize(0);
	hRhoRatio[p][e][c]->GetYaxis()->SetLabelSize(0.06);
	//hRhoRatio[p][e][c]->GetYaxis()->SetNdivisions(10);

	int dir = 1+e+(3*p);
	c0->cd(dir);

	hRhoRatio[p][e][c]->Draw("PSAME");


      }
    }
  }

  c0->SaveAs(saveName,"PDF");
  
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
