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
  
  TString fileName = "out/UE/pAuHTjetUE_inLead.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TTree *jetTree = (TTree*) inFile->Get("HTjetTree");

  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcSumEast, leadPt, leadPtCorrected, leadEta, leadPhi, chgEastRho_te, chgMidRho_te, chgWestRho_te, neuEastRho, neuMidRho, neuWestRho;

  jetTree->SetBranchAddress( "leadPt", &leadPt );
  jetTree->SetBranchAddress( "leadPtCorrected", &leadPtCorrected );
  jetTree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );
  jetTree->SetBranchAddress( "leadEta", &leadEta );
  jetTree->SetBranchAddress( "leadPhi", &leadPhi );
  jetTree->SetBranchAddress( "chgEastRho_te", &chgEastRho_te );
  jetTree->SetBranchAddress( "chgMidRho_te", &chgMidRho_te );
  jetTree->SetBranchAddress( "chgWestRho_te", &chgWestRho_te );
  jetTree->SetBranchAddress( "neuEastRho", &neuEastRho );
  jetTree->SetBranchAddress( "neuMidRho", &neuMidRho );
  jetTree->SetBranchAddress( "neuWestRho", &neuWestRho );

  int nEntries = jetTree->GetEntries();

  const int nEAbins = 2;
  TString EAbinString[nEAbins] = { "Lo", "Hi" };
  TString BBCselection[nEAbins] = { "BbcAdcSumEast>3559.12 && BbcAdcSumEast<11503", "BbcAdcSumEast>26718.1" };
 
  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString rhoVal[nEtaBins] = { "(chgEastRho_te+neuEastRho)", "(chgMidRho_te+neuMidRho)", "(chgWestRho_te+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPtCorrected>10.0 && leadPtCorrected<15.0", "leadPtCorrected>=15.0 && leadPtCorrected<=20.0", "leadPtCorrected>20.0 && leadPtCorrected<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 23, 22 };

  TH1D *hRho_pt[nPtBins][nEtaBins][nEAbins];
  TH1D *hRhoByEta_pt[nPtBins][nEAbins][nEtaBins];


  TH2D *hscale = new TH2D("hscale","",3,-1.5,1.5,10,0.6,1.7);
  hscale->SetStats(0);
  hscale->GetXaxis()->SetLabelSize(0);
  hscale->GetYaxis()->SetLabelSize(0.08);
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int p=0; p<nPtBins; ++p ) {

      for ( int e=0; e<nEtaBins; ++e ) {

	name = "hRhoByEta" + EAbinString[a] + ptBinName[p];
	hRhoByEta_pt[p][a][e] = new TH1D( name, "", 3,-1.5,1.5 );
	hRhoByEta_pt[p][a][e]->SetStats(0);
	hRhoByEta_pt[p][a][e]->SetMarkerStyle( EAmarker[a] );
	hRhoByEta_pt[p][a][e]->SetMarkerSize( 1.5 );
	hRhoByEta_pt[p][a][e]->SetMarkerColor( etaColor[e] );
	hRhoByEta_pt[p][a][e]->SetLineColor( etaColor[e] );
	hRhoByEta_pt[p][a][e]->GetYaxis()->SetRangeUser(0.6,1.7);
	hRhoByEta_pt[p][a][e]->GetXaxis()->SetLabelSize(0);
	hRhoByEta_pt[p][a][e]->GetYaxis()->SetLabelSize(0.06);
      
	//jetTree->Draw("(chgEastRho_te+neuEastRho)>>hEastHi_10_15","leadPtCorrected>10.0 && leadPtCorrected<15.0 && BbcAdcSumEast>26718.1","");
	hscale->Draw();
	name = "h" + eastmidwest[e] + EAbinString[a] + ptBinName[p];
	TString drawString = rhoVal[e] + ">>" + name;
	TString drawCuts = ptSelection[p] + " && " + BBCselection[a];
	jetTree->Draw( drawString, drawCuts, "" );
	hRho_pt[p][e][a] = (TH1D*)gDirectory->Get( name );
	hRho_pt[p][e][a]->SetStats(0);

	rho = hRho_pt[p][e][a]->GetMean(1);
	stdev = hRho_pt[p][e][a]->GetMeanError(1);
	hRhoByEta_pt[p][a][e]->SetBinContent( e+1, rho );
	hRhoByEta_pt[p][a][e]->SetBinError( e+1, stdev );

      }
    }
  }

  TCanvas *cpt = new TCanvas("cpt","", 200 ,600 );
  cpt->Divide(1,3,0,0);

  for ( int p=0; p<nPtBins; ++p ) {
    cpt->cd(p+1);
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetGridy();
    for ( int a=0; a<nEAbins; ++a ) {
      hscale->Draw("SAME");

      for (int binno=0; binno<3; ++binno) {  
	hRhoByEta_pt[p][a][binno]->GetYaxis()->SetLabelSize(0.08);
	hRhoByEta_pt[p][a][binno]->SetMarkerColor( etaColor[binno] );
	hRhoByEta_pt[p][a][binno]->Draw("SAME");
      }

    }
  }



  new TCanvas;
  
  const TString etaSuffix[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };

  TH1D *hRho_eta[nEtaBins][nPtBins][nEAbins];
  TH1D *hRhoByPt_eta[nPtBins][nEAbins][nEtaBins];
  
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int je=0; je<nEtaBins; ++je ) {
      for (int binno=0; binno<3; ++binno) {  

	name = "hRhoByEta" + EAbinString[a] + etaSuffix[je];
	hRhoByPt_eta[je][a][binno] = new TH1D( name, "", 3,-1.5,1.5 );
	hRhoByPt_eta[je][a][binno]->SetStats(0);
	hRhoByPt_eta[je][a][binno]->SetMarkerStyle( EAmarker[a] );
	hRhoByPt_eta[je][a][binno]->SetMarkerSize( 1.5 );
	hRhoByPt_eta[je][a][binno]->SetMarkerColor( etaColor[binno] );
	hRhoByPt_eta[je][a][binno]->SetLineColor( etaColor[binno] );
	hRhoByPt_eta[je][a][binno]->GetYaxis()->SetRangeUser(0.6,1.7);
	hRhoByPt_eta[je][a][binno]->GetXaxis()->SetLabelSize(0);
	hRhoByPt_eta[je][a][binno]->GetYaxis()->SetLabelSize(0.06);
      }
    }
  }
  for ( int a=0; a<nEAbins; ++a ) {
    for ( int je=0; je<nEtaBins; ++je ) {
      for ( int e=0; e<nEtaBins; ++e ) {
	//jetTree->Draw("(chgEastRho_te+neuEastRho)>>hEastHi_10_15","leadPtCorrected>10.0 && leadPtCorrected<15.0 && BbcAdcSumEast>26718.1","");
	name = "h" + eastmidwest[e] + EAbinString[a] + ptBinName[je];
	TString drawString = rhoVal[e] + ">>" + name;
	TString drawCuts = etaSelection[e] + " && " + BBCselection[a];
	//cout<<name<<endl<<drawString<<endl<<drawCuts<<endl<<endl;
	jetTree->Draw( drawString, drawCuts, "" );
	hRho_eta[je][e][a] = (TH1D*)gDirectory->Get( name );
	hRho_eta[je][e][a]->SetStats(0);

	rho = hRho_eta[je][e][a]->GetMean(1);
	stdev = hRho_eta[je][e][a]->GetMeanError(1);
	hRhoByPt_eta[je][a][e]->SetBinContent( e+1, rho );
	hRhoByPt_eta[je][a][e]->SetBinError( e+1, stdev );

      }
    }
  }


  TCanvas *ceta = new TCanvas("ceta","", 700, 350 );
  ceta->Divide(3,1,0,0);

  for ( int je=0; je<nEtaBins; ++je ) {
    ceta->cd(je+1);
    hscale->Draw("SAME");
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->SetGridy();
    for ( int a=0; a<nEAbins; ++a ) {
      for (int binno=0; binno<3; ++binno) {        
      hRhoByPt_eta[je][a][binno]->GetYaxis()->SetLabelSize(0.08);
      hRhoByPt_eta[je][a][binno]->Draw("SAME");
      }
    }
  }


  cpt->SaveAs("plots/UE/pt3plot.pdf","PDF");
  ceta->SaveAs("plots/UE/eta3plot.pdf","PDF");
  
}
