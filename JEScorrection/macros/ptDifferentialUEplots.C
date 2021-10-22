// Veronica Verkest
// November 15, 2020

void ptDifferentialUEplots(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  gStyle->SetErrorX(0.0001);
  
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
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta^{lead}_{jet}<-0.3", "-0.3<#eta^{lead}_{jet}<0.3", "0.3<#eta^{lead}_{jet}<0.6" };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgUE", "_neuUE", "_allUE" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

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
  const TString lohi[nEAbins] = { "lo", "hi" };

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  string directory = "plots/JEScorrection/";

  auto hs_n = new THStack("hs_n","dNch/dEtadPhi");
  auto hs_pt = new THStack("hs_pt","mean pT");

  const int etaColor[nEtaBins] = { 877, 596, 814 };

  const int n_bins = 3;
  double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  double shiftedBins[n_bins+1];
  const int n_ybins = 3;
  double y_bin_edge[n_ybins+1] = { 0.55,0.65,0.75,0.85 };
  TH2D *meanPt_hscale = new TH2D("meanPt_hscale",";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge);

  double bin_edge2[n_bins+1] = { 0.0, 1.0, 2.0, 3.0 };
  double y_bin_edge2[n_ybins+1] = { 0.5,1.2,1.5,1.8 };
  TH2D *nCh_hscale = new TH2D("nCh_hscale",";leading jet p_{T} (GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT",n_bins,bin_edge,n_ybins,y_bin_edge2);

  TH1D *hPt[nPtBins][nEtaBins][nEAbins];

  // TFile *inFile = new TFile("out/pAuUE_pt.root","READ");
  TFile *inFile = new TFile("out/JEScorrection.root","READ");
  
  TH1D *hNch[nEtaBins][nEAbins]; 
  TH1D *hMeanPt[nEtaBins][nEAbins];
  
  for (int a=1; a>=0; --a) {
    for (int e=0; e<nEtaBins; ++e) {
      nCh_hscale->GetXaxis()->SetBinLabel(e+1,ptBinString[e]);
      for (int i=0; i<=n_bins; ++i) { shiftedBins[i] = bin_edge[i] + 0.3*(e-1); }
      name = "hNch" + etaBinName[e] + EAbinName[a];
      hNch[e][a] = new TH1D(name,";leading jet p_{T} (GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT",n_bins,shiftedBins);
      hNch[e][a]->SetLineColor(etaColor[e]);
      hNch[e][a]->SetMarkerColor(etaColor[e]);
      hNch[e][a]->SetMarkerStyle(EAmarker[a]);
      hNch[e][a]->SetMarkerColor(etaColor[e]);

      name = "hMeanPt" + etaBinName[e] + EAbinName[a];
      hMeanPt[e][a] = new TH1D(name,";;#LT p_{T}^{ch}#GT (GeV)",n_bins,shiftedBins);
      hMeanPt[e][a]->SetLineColor(etaColor[e]);
      hMeanPt[e][a]->SetMarkerColor(etaColor[e]);
      hMeanPt[e][a]->SetMarkerStyle(EAmarker[a]);
      hMeanPt[e][a]->SetMarkerColor(etaColor[e]);
      for (int p=0; p<nPtBins; ++p) {

    	// name = "hPartJetUE" + ptBinName[p] + etaBinName[e] + "_" + lohi[a];
    	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e];
    	hPt[p][e][a] = (TH1D*)inFile->Get(name);
	
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
    	hPt[p][e][a]->SetLineColor(etaColor[e]);
    	hPt[p][e][a]->SetMarkerStyle(EAmarker[a]);
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
	
    	hNch[e][a]->SetBinContent(p+1,hPt[p][e][a]->Integral("WIDTH")/area[e]);
    	// hNch[e][a]->SetBinError(p+1,hPt[p][e][a]->GetMeanError(1));

	double mean = 0.;
	double totalContent = 0.;
	for (int i=0; i<hPt[p][e][a]->GetNbinsX(); ++i) {
	  int binno = i+1;
	  double width = hPt[p][e][a]->GetXaxis()->GetBinWidth(binno);
	  mean += hPt[p][e][a]->GetBinContent( binno )*hPt[p][e][a]->GetBinCenter( binno )*width;///hPt[p][e][a]->Integral("width");
	  totalContent += hPt[p][e][a]->GetBinContent( binno )*width;
	}
	mean /= hPt[p][e][a]->Integral("width");
	hMeanPt[e][a]->SetBinContent(p+1,mean);
      }
      hNch[e][a]->GetYaxis()->SetRangeUser(0.5,1.8);
      hs_n->Add(hNch[e][a]);
      hs_pt->Add(hMeanPt[e][a]);
    }
  }
  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
  
  nCh_hscale->SetStats(0);
  nCh_hscale->Draw();
  hs_n->Draw("SAMEnostackEX0pf");   // draw the stack
  name = directory + "nCh.pdf";
  can->SaveAs( name, "PDF");

  meanPt_hscale->SetStats(0);
  meanPt_hscale->Draw();
  hs_pt->Draw("SAMEnostackEX0pf");   // draw the stack
  name = directory + "meanPt.pdf";
  can->SaveAs( name, "PDF");

  
  can->Destructor();

}

    

