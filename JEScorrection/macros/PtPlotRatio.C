// Veronica Verkest
// July 19, 2021

void PtPlotRatio(){

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

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  // TString directory = "plots/compare/";
  // TString directory = "SYSTEMATICS/plots/";
  TString directory = "SYSTEMATICS/plots/JEScorrection_pos_";

  TFile* file1 = new TFile("out/JEScorrection.root", "READ");
  TFile *file2 = new TFile("SYSTEMATICS/out/JEScorrection_pos.root","READ");
  // TFile *file2 = new TFile("SYSTEMATICS/out/JEScorrection_det_nom.root","READ");
  // TFile* file2 = new TFile("out/noCorrection.root", "READ");
  // TFile* file1 = new TFile("out/noCorrection_halfGeVbins_new.root", "READ");
  // TFile* file2 = new TFile("out/noCorrection_halfGeVbins_prelimTrackEffic.root", "READ");
  
  auto hs_n = new THStack("hs_n","dNch/dEtadPhi");
  auto hs_pt = new THStack("hs_pt","mean pT");

  auto hs_nRATIO = new THStack("hs_nRATIO","");
  auto hs_ptRATIO = new THStack("hs_ptRATIO","");

  const int etaColor[nEtaBins] = { 877, 596, 814 };

  const int n_bins = 3;
  double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  double shiftedBins[n_bins+1];
  const int n_ybins = 3;
  // double y_bin_edge[n_ybins+1] = { 0.55,0.65,0.75,0.85 };
  double y_bin_edge[n_ybins+1] = { 0.45,0.65,0.75,0.95 };
  TH2D *meanPt_hscale = new TH2D("meanPt_hscale",";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge);

  double bin_edge2[n_bins+1] = { 0.0, 1.0, 2.0, 3.0 };
  // double y_bin_edge2[n_ybins+1] = { 0.5,1.2,1.5,1.8 };
  double y_bin_edge2[n_ybins+1] = { 0.6,1.2,1.5,1.9 };
  TH2D *nCh_hscale = new TH2D("nCh_hscale",";leading jet p_{T} (GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT",n_bins,bin_edge,n_ybins,y_bin_edge2);

  TH1D *hPt[nPtBins][nEtaBins][nEAbins];
  TH1D *hPt_dbw[nPtBins][nEtaBins][nEAbins];
  TH1D *hPt2[nPtBins][nEtaBins][nEAbins];
  TH1D *hPt_dbw2[nPtBins][nEtaBins][nEAbins];
  
  TH1D *hNch[nEtaBins][nEAbins]; 
  TH1D *hMeanPt[nEtaBins][nEAbins];
  TH1D *hNch2[nEtaBins][nEAbins];
  TH1D *hMeanPt2[nEtaBins][nEAbins];

  TH1D *hNchRATIO[nEtaBins][nEAbins]; 
  TH1D *hPtRATIO[nEtaBins][nEAbins]; 

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

      name = "hNch2" + etaBinName[e] + EAbinName[a];
      hNch2[e][a] = new TH1D(name,";leading jet p_{T} (GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT",n_bins,shiftedBins);
      hNch2[e][a]->SetLineColor(etaColor[e]);
      hNch2[e][a]->SetMarkerColor(etaColor[e]);
      hNch2[e][a]->SetMarkerStyle(EAmarker2[a]);
      hNch2[e][a]->SetMarkerColor(etaColor[e]);

      name = "hMeanPt2" + etaBinName[e] + EAbinName[a];
      hMeanPt2[e][a] = new TH1D(name,";;#LT p_{T}^{ch}#GT (GeV)",n_bins,shiftedBins);
      hMeanPt2[e][a]->SetLineColor(etaColor[e]);
      hMeanPt2[e][a]->SetMarkerColor(etaColor[e]);
      hMeanPt2[e][a]->SetMarkerStyle(EAmarker2[a]);
      hMeanPt2[e][a]->SetMarkerColor(etaColor[e]);
      for (int p=0; p<nPtBins; ++p) {

	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e] ;
    	hPt[p][e][a] = (TH1D*)file1->Get(name);
	name += "_dbw";
    	hPt_dbw[p][e][a] = (TH1D*)file1->Get(name);
	
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);
    	hPt[p][e][a]->SetLineColor(etaColor[e]);
    	hPt[p][e][a]->SetMarkerStyle(EAmarker[a]);
    	hPt[p][e][a]->SetMarkerColor(etaColor[e]);

	double intErr;
	double integral = hPt[p][e][a]->IntegralAndError(1,hPt[p][e][a]->GetNbinsX(),intErr);
    	hNch[e][a]->SetBinContent(p+1,integral/area[e]);
    	hNch[e][a]->SetBinError(p+1,intErr/area[e]);

	hMeanPt[e][a]->SetBinContent(p+1,hPt[p][e][a]->GetMean(1));
	hMeanPt[e][a]->SetBinError(p+1,hPt[p][e][a]->GetMeanError(1));




	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e] ;
    	hPt2[p][e][a] = (TH1D*)file2->Get(name);
	name += "_dbw";
    	hPt_dbw2[p][e][a] = (TH1D*)file2->Get(name);
	
    	hPt2[p][e][a]->SetMarkerColor(etaColor[e]);
    	hPt2[p][e][a]->SetLineColor(etaColor[e]);
    	hPt2[p][e][a]->SetMarkerStyle(EAmarker2[a]);
    	hPt2[p][e][a]->SetMarkerColor(etaColor[e]);

	double intErr2;
	double integral2 = hPt2[p][e][a]->IntegralAndError(1,hPt2[p][e][a]->GetNbinsX(),intErr2);
    	hNch2[e][a]->SetBinContent(p+1,integral2/area[e]);
    	hNch2[e][a]->SetBinError(p+1,intErr2/area[e]);

	hMeanPt2[e][a]->SetBinContent(p+1,hPt2[p][e][a]->GetMean(1));
	hMeanPt2[e][a]->SetBinError(p+1,hPt2[p][e][a]->GetMeanError(1));
	
      }
      hNch[e][a]->GetYaxis()->SetRangeUser(0.5,1.8);
      hNch2[e][a]->GetYaxis()->SetRangeUser(0.5,1.8);
      hs_n->Add(hNch[e][a]);
      hs_pt->Add(hMeanPt[e][a]);
      hs_n->Add(hNch2[e][a]);
      hs_pt->Add(hMeanPt2[e][a]);

      hNchRATIO[e][a] = (TH1D*)hNch[e][a]->Clone();
      hNchRATIO[e][a]->Divide(hNch2[e][a]);
      hs_nRATIO->Add(hNchRATIO[e][a]);
      
      hPtRATIO[e][a] = (TH1D*)hMeanPt[e][a]->Clone();
      hPtRATIO[e][a]->Divide(hMeanPt2[e][a]);
      hs_ptRATIO->Add(hPtRATIO[e][a]);
    }
  }
  TCanvas * can = new TCanvas( "can" , "" ,700 ,700 );
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();

  TH2D *hRatioScale = new TH2D("hRatioScale","", 5,10,30, 5,0.995,1.005 );
  hRatioScale->SetStats(0);
    
  pad1->cd();
  nCh_hscale->SetStats(0);
  nCh_hscale->Draw();
  hs_n->Draw("SAMEnostackEX0pf");   // draw the stack
  pad2->cd();
  hRatioScale->Draw();
  hs_nRATIO->Draw("SAMEnostackEX0pf");
  name = directory + "nCh.pdf";
  can->SaveAs( name, "PDF");


  pad1->cd();
  meanPt_hscale->SetStats(0);
  meanPt_hscale->Draw();
  hs_pt->Draw("SAMEnostackEX0pf");   // draw the stack
  pad2->cd();
  hRatioScale->Draw();
  hs_ptRATIO->Draw("SAMEnostackEX0pf");
  name = directory + "meanPt.pdf";
  can->SaveAs( name, "PDF");

  
  can->Destructor();

}

    

