void UEdiffRho(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  gStyle->SetErrorX(0.0001);
  
  const double pi = 3.14159265;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const int ptLoVal[nPtBins] = { 10, 15, 20 };
  const int ptHiVal[nPtBins] = { 15, 20, 30 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };

  int eval, pval;
  TString name, saveName, title, avg, sigma;
  double chgRho, neuRho, midRho, eastRho, westRho, rho, stdev;
  
  TString fileName = "out/UE/pAuHTjetUE.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  jetTree->SetBranchAddress( "RunID", &RunID );      	       		jetTree->SetBranchAddress( "EventID", &EventID );					jetTree->SetBranchAddress( "nTowers", &nTowers );
  jetTree->SetBranchAddress( "nPrimary", &nPrimary );       		jetTree->SetBranchAddress( "nGlobal", &nGlobal );					jetTree->SetBranchAddress( "nVertices", &nVertices );
  jetTree->SetBranchAddress( "refMult", &refMult );			jetTree->SetBranchAddress( "gRefMult", &gRefMult );		       		jetTree->SetBranchAddress( "Vz", &Vz );
  jetTree->SetBranchAddress( "leadPt", &leadPt );	       			jetTree->SetBranchAddress( "BbcAdcEastSum", &BbcAdcEastSum );	jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );	       		jetTree->SetBranchAddress( "chgEastRho", &chgEastRho );	       		jetTree->SetBranchAddress( "chgMidRho", &chgMidRho );
  jetTree->SetBranchAddress( "chgWestRho", &chgWestRho );	jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );			jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );

  int nEntries = jetTree->GetEntries();

  const int nEAbins = 2;
  TString EAbinString[nEAbins] = { "Lo", "Hi" };
  TString BBCselection[nEAbins] = { "BbcAdcEastSum>4107 && BbcAdcEastSum<11503", "BbcAdcEastSum>28537" };

  TH1D *hRho_pt[nPtBins][nEtaBins][nEAbins];

  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho+neuEastRho)", "(chgMidRho+neuMidRho)", "(chgWestRho+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };

  int EAcolor[nEAbins] = { 810, 884 };
  int EAmarker[nEAbins] = { 24, 25 };

  TH1D *hRhoByEta_pt[nPtBins][nEAbins];

  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {

      name = "hRhoByEta" + EAbinString[a] + ptBinName[p];
      hRhoByEta_pt[p][a] = new TH1D( name, "", 3,-1.5,1.5 );

      for ( int e=0; e<nEtaBins; ++e ) {
	//jetTree->Draw("(chgEastRho+neuEastRho)>>hEastHi_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>28537","");
	name = "h" + eastmidwest[e] + EAbinString[a] + ptBinName[p];
	TString drawString = rhoVal[e] + ">>h" + name;
	TString drawCuts = ptSelection[p] + " && " + BBCselection[a];
	jetTree->Draw( drawString, drawCuts, "" );
	hRho_pt[p][e][a] = (TH1D*)gDirectory->Get( name );
	hRho_pt[p][e][a]->SetStats(0);

	int ebin = e+1;

	rho = hRho_pt[p][e][a]->GetMean(1);
	stdev = hRho_pt[p][e][a]->GetMeanError(1);
      
	hRhoByEta_pt[p][a]->SetBinContent( e, rho );
	hRhoByEta_pt[p][a]->SetBinError( e, stdev );
	hRhoByEta_pt[p][a]->SetMarkerStyle( EAmarker[e] );
	hRhoByEta_pt[p][a]->SetMarkerColor( EAcolor[e] );
	hRhoByEta_pt[p][a]->SetLineColor( EAcolor[e] );	

      }
    }
  }


  TCanvas *cpt = new TCanvas( 0, 23, 700, 1280 );
  cpt->Divide(1,3,0,0);

  for ( int p=0; p<nPtBins; ++p ) {
    cpt->cd(p+1);

    for ( int a=0; a<nEAbins; ++a ) {
      hRhoByEta_pt[p][a]->Draw("PSAME");
    }
  }




  
  // TH1D *hEastHi_eta[nEtaBins];
  // TH1D *hMidHi_eta[nEtaBins];
  // TH1D *hWestHi_eta[nEtaBins];

  // TH1D *hEastLo_eta[nEtaBins];
  // TH1D *hMidLo_eta[nEtaBins];
  // TH1D *hWestLo_eta[nEtaBins];


  

  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastHi_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>28537","");
  // TH1D *hEastHi_10_15 = (TH1D*)gDirectory->Get("hEastHi_10_15");  hEastHi_10_15->Scale(1./hEastHi_10_15->Integral());
  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastHi_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>28537","");
  // TH1D *hEastHi_15_20 = (TH1D*)gDirectory->Get("hEastHi_15_20");  hEastHi_15_20->Scale(1./hEastHi_15_20->Integral());
  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastHi_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>28537","");
  // TH1D *hEastHi_20_30 = (TH1D*)gDirectory->Get("hEastHi_20_30");  hEastHi_20_30->Scale(1./hEastHi_20_30->Integral());

  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidHi_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>28537","");
  // TH1D *hMidHi_10_15 = (TH1D*)gDirectory->Get("hMidHi_10_15");  hMidHi_10_15->Scale(1./hMidHi_10_15->Integral());
  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidHi_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>28537","");
  // TH1D *hMidHi_15_20 = (TH1D*)gDirectory->Get("hMidHi_15_20");  hMidHi_15_20->Scale(1./hMidHi_15_20->Integral());
  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidHi_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>28537","");
  // TH1D *hMidHi_20_30 = (TH1D*)gDirectory->Get("hMidHi_20_30");  hMidHi_20_30->Scale(1./hMidHi_20_30->Integral());

  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestHi_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>28537","");
  // TH1D *hWestHi_10_15 = (TH1D*)gDirectory->Get("hWestHi_10_15");  hWestHi_10_15->Scale(1./hWestHi_10_15->Integral());
  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestHi_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>28537","");
  // TH1D *hWestHi_15_20 = (TH1D*)gDirectory->Get("hWestHi_15_20");  hWestHi_15_20->Scale(1./hWestHi_15_20->Integral());
  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestHi_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>28537","");
  // TH1D *hWestHi_20_30 = (TH1D*)gDirectory->Get("hWestHi_20_30");  hWestHi_20_30->Scale(1./hWestHi_20_30->Integral());


  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastLo_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hEastLo_10_15 = (TH1D*)gDirectory->Get("hEastLo_10_15");  hEastLo_10_15->Scale(1./hEastLo_10_15->Integral());
  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastLo_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hEastLo_15_20 = (TH1D*)gDirectory->Get("hEastLo_15_20");  hEastLo_15_20->Scale(1./hEastLo_15_20->Integral());
  // jetTree->Draw("(chgEastRho+neuEastRho)>>hEastLo_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hEastLo_20_30 = (TH1D*)gDirectory->Get("hEastLo_20_30");  hEastLo_20_30->Scale(1./hEastLo_20_30->Integral());

  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidLo_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hMidLo_10_15 = (TH1D*)gDirectory->Get("hMidLo_10_15");  hMidLo_10_15->Scale(1./hMidLo_10_15->Integral());
  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidLo_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hMidLo_15_20 = (TH1D*)gDirectory->Get("hMidLo_15_20");  hMidLo_15_20->Scale(1./hMidLo_15_20->Integral());
  // jetTree->Draw("(chgMidRho+neuMidRho)>>hMidLo_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hMidLo_20_30 = (TH1D*)gDirectory->Get("hMidLo_20_30");  hMidLo_20_30->Scale(1./hMidLo_20_30->Integral());

  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestLo_10_15","leadPt>10.0 && leadPt<15.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hWestLo_10_15 = (TH1D*)gDirectory->Get("hWestLo_10_15");  hWestLo_10_15->Scale(1./hWestLo_10_15->Integral());
  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestLo_15_20","leadPt>15.0 && leadPt<20.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hWestLo_15_20 = (TH1D*)gDirectory->Get("hWestLo_15_20");  hWestLo_15_20->Scale(1./hWestLo_15_20->Integral());
  // jetTree->Draw("(chgWestRho+neuWestRho)>>hWestLo_20_30","leadPt>20.0 && leadPt<30.0 && BbcAdcEastSum>4107 && BbcAdcEastSum<11503","");
  // TH1D *hWestLo_20_30 = (TH1D*)gDirectory->Get("hWestLo_20_30");  hWestLo_20_30->Scale(1./hWestLo_20_30->Integral());
  
}
