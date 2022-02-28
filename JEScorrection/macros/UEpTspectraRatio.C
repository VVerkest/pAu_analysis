// Veronica Verkest
// July 21, 2021

void UEpTspectraRatio(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };
  const TString ptCorrectedBinString[nPtBins] = { "10 < p_{T} < 15", "15 < p_{T} < 20",  "20 < p_{T} < 30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };
  const int ptMarker[nPtBins] = { 20, 21, 29 };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta^{lead}_{jet}<-0.3", "-0.3<#eta^{lead}_{jet}<0.3", "0.3<#eta^{lead}_{jet}<0.6" };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };

  const int nEAbins = 2;
  TString EAbinName[nEAbins] = { "Lo", "Hi" };
  TString EAbinString[nEAbins] = { "Low EA", "High EA" };
  TString BBCselection[nEAbins] = { "BbcAdcSumEast>3559.12 && BbcAdcSumEast<11503", "BbcAdcSumEast>26718.1" };
 
  TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  TString emw[nEtaBins] = { "east", "mid", "west" };
  TString rhoVal[nEtaBins] = { "(chgEastRho_te+neuEastRho)", "(chgMidRho_te+neuMidRho)", "(chgWestRho_te+neuWestRho)" };
  TString ptSelection[nPtBins] = { "leadPt>10.0 && leadPt<15.0", "leadPt>=15.0 && leadPt<=20.0", "leadPt>20.0 && leadPt<30.0" };
  TString etaSelection[nEtaBins] = { "leadEta>=-0.6 && leadEta<=-0.3", "leadEta>-0.3 && leadEta<0.3", "leadEta>=0.3 && leadEta<=0.6" };

  const double eastArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  double area[nEtaBins] = { eastArea, midArea, westArea };
  
  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 24, 20 };
  int EAmarker2[nEAbins] = { 25, 21 };
  const TString lohi[nEAbins] = { "lo", "hi" };

  const int xbins = 55;
  const double xbinEdge[xbins+1] = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0};
  const int zbins = 20;
  const double zbinEdge[zbins+1] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const int ybins = 30;
  const double ybinEdge[ybins+1] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
  				     10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0 };
  
  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  // TString directory = "plots/compare/";
  // TString directory = "SYSTEMATICS/plots/";
  TString directory = "SYSTEMATICS/plots/JEScorrection_neg_";

  TFile* file1 = new TFile("out/JEScorrection.root", "READ");
  // TFile *file2 = new TFile("SYSTEMATICS/out/JEScorrection_det_TU.root","READ");
  TFile *file2 = new TFile("SYSTEMATICS/out/JEScorrection_neg.root","READ");
  // TFile* file2 = new TFile("out/noCorrection.root", "READ");
  // TFile* file2 = new TFile("out/noCorrection_halfGeVbins_prelimTrackEffic.root", "READ");

  TH1D *hPt[nPtBins][nEtaBins][nEAbins];
  TH1D *hPt2[nPtBins][nEtaBins][nEAbins];
  TH1D *hPtRATIO[nPtBins][nEtaBins][nEAbins];

  for (int a=0; a<nEAbins; ++a) {   //   for (int a=1; a>=0; --a) {
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {

	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e] ;
    	hPt[p][e][a] = (TH1D*)file1->Get(name);
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
    	hPt[p][e][a]->SetLineColor(etaColor[e]);
    	hPt[p][e][a]->SetMarkerStyle(EAmarker[a]);
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
	
    	hPt2[p][e][a] = (TH1D*)file2->Get(name);
    	hPt2[p][e][a]->SetMarkerColor(etaColor[e]);
    	hPt2[p][e][a]->SetLineColor(etaColor[e]);
    	hPt2[p][e][a]->SetMarkerStyle(EAmarker2[a]);
    	hPt2[p][e][a]->SetMarkerColor(etaColor[e]);
	
	name += "_RATIO";
	hPtRATIO[p][e][a] = (TH1D*)hPt[p][e][a]->Clone(name);
	hPtRATIO[p][e][a]->Divide(hPt2[p][e][a]);
      }

    }
  }

  double canFrac[4] = { 0.0, 0.33333, 0.66667, 1.0 };
  // double canFrac[4] = { 1.0, 0.66667, 0.33333, 0.0 };

  TCanvas * can = new TCanvas( "can" , "" ,1400 ,1000 );
  TPad *pad[nPtBins][nEtaBins];

  for (int p=0; p<nPtBins; ++p) {  // pT along horizontal
    for (int e=0; e<nEtaBins; ++e) {  // eta along vertical
      name = "pad" + ptBinString[p] + etaBinString[e];
      pad[p][e] = new TPad(name,"",canFrac[p],canFrac[e],canFrac[p+1],canFrac[e+1]);
      pad[p][e]->Draw();
    }
  }

  // TH2D *hScale = new TH2D("hScale","", 15,0.2,15., 5, 0.0, 2. );
  // TH2D *hScale = new TH2D("hScale","", 15,0.2,15., 5, 0.5, 1.5 );
  TH2D *hScale = new TH2D("hScale","", 15,0.2,15., 5, 0.99, 1.01 );
  hScale->SetStats(0);
  
  for (int p=0; p<nPtBins; ++p) {  // pT along horizontal
    for (int e=0; e<nEtaBins; ++e) {  // eta along vertical
      pad[p][e]->cd();
      hScale->Draw();
      for (int a=1; a>=0; --a) {
	hPtRATIO[p][e][a]->Draw("psame");
      }
    }
  }

  name = directory + "UEpTspectraRatios.pdf";
  can->SaveAs(name,"PDF");

  
}
