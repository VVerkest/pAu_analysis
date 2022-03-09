// Veronica Verkest
// June 17, 2021
// makes plots using the output file from newJEScorrection

void FakesAndMissesComparison(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  // gStyle->SetErrorX(0.0001);
  
  const double pi = 3.14159265;
  const double AREA = 4.*(pi/3.);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

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

  const double eastArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )
  const double midArea = 2.*(0.6)*(pi/3.);   // (  0.6 in eta  ) X (  2pi/3 in phi  )
  const double westArea = 2.*(0.7)*(pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )
  double area[nEtaBins] = { eastArea, midArea, westArea };
  
  int EAcolor[nEAbins] = { 884, 810 };
  int EAmarker[nEAbins] = { 24, 20 };
  const TString lohi[nEAbins] = { "lo", "hi" };

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  const int nFiles = 4;
  TString directory = "plots/JEScorrection/";
  TFile* inFile[nFiles];
  TString filePrefix[nFiles] = { "_noMissOrFake", "_MissCorrected", "_FakeCorrected", "" };

  for (int f=0; f<nFiles; ++f) {
    name = "out/JEScorrection" + filePrefix[f] + ".root";
    inFile[f] = new TFile(name, "READ");
  }


  TH1D *hPt[nFiles][nPtBins][nEtaBins][nEAbins], *hNch[nFiles][nEtaBins][nEAbins], *hMeanPt[nFiles][nEtaBins][nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    for (int f=0; f<nFiles; ++f) {
      for (int e=0; e<nEtaBins; ++e) {
	for (int p=0; p<nPtBins; ++p) {
	  name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e];
	  hPt[f][p][e][a] = (TH1D*)inFile[f]->Get(name);
	  name += filePrefix[f];
	  hPt[f][p][e][a]->SetName(name);
	}
      }
    }
  }

  hPt[0][0][2][1]->Divide(hPt[3][0][2][1]);
  hPt[0][0][2][1]->Draw();
  
  // for (int a=1; a>=0; --a) {
  //   for (int e=0; e<nEtaBins; ++e) {
  //     nCh_hscale->GetXaxis()->SetBinLabel(e+1,ptBinString[e]);
  //     for (int i=0; i<=n_bins; ++i) { shiftedBins[i] = bin_edge[i] + 0.3*(e-1); }
  //     name = "hNch" + etaBinName[e] + EAbinName[a];
  //     hNch[e][a] = new TH1D(name,";leading jet p_{T} (GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT",n_bins,shiftedBins);
  //     hNch[e][a]->SetLineColor(etaColor[e]);
  //     hNch[e][a]->SetMarkerColor(etaColor[e]);
  //     hNch[e][a]->SetMarkerStyle(EAmarker[a]);
  //     hNch[e][a]->SetMarkerColor(etaColor[e]);

  //     name = "hMeanPt" + etaBinName[e] + EAbinName[a];
  //     hMeanPt[e][a] = new TH1D(name,";;#LT p_{T}^{ch}#GT (GeV)",n_bins,shiftedBins);
  //     hMeanPt[e][a]->SetLineColor(etaColor[e]);
  //     hMeanPt[e][a]->SetMarkerColor(etaColor[e]);
  //     hMeanPt[e][a]->SetMarkerStyle(EAmarker[a]);
  //     hMeanPt[e][a]->SetMarkerColor(etaColor[e]);
  //     for (int p=0; p<nPtBins; ++p) {

  // 	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e] ;
  //   	hPt[p][e][a] = (TH1D*)inFile->Get(name);

  // 	name += "_dbw";
  //   	hPt_dbw[p][e][a] = (TH1D*)inFile->Get(name);
	
  //   	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
  //   	hPt[p][e][a]->SetLineColor(etaColor[e]);
  //   	hPt[p][e][a]->SetMarkerStyle(EAmarker[a]);
  //   	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
	
  //   	hNch[e][a]->SetBinContent(p+1,hPt[p][e][a]->Integral()/area[e]);
  //   	// hNch[e][a]->SetBinError(p+1,hPt[p][e][a]->GetMeanError(1));

  // 	hMeanPt[e][a]->SetBinContent(p+1,hPt[p][e][a]->GetMean(1));
  //     }
  //     hNch[e][a]->GetYaxis()->SetRangeUser(0.5,1.8);
  //     hs_n->Add(hNch[e][a]);
  //     hs_pt->Add(hMeanPt[e][a]);
  //   }
  // }

  
}
