void makePtDiffPlot(){

  //gStyle->SetErrorX(0.0001);
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString jetEtaBinName[nEtaBins] = { "_eastJet", "_midJet", "_westJet" };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta^{lead}_{jet}<-0.3", "-0.3<#eta^{lead}_{jet}<0.3", "0.3<#eta^{lead}_{jet}<0.6" };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const int etaMarker[nEtaBins] = { 25, 27, 28 };
  int EAptMarker[2][nPtBins] = {{107,108,110},{20,21,33}};

  const int nEAbins = 2;
  TString EAbinName[nEAbins] = { "Lo", "Hi" };
  TString EAbinString[nEAbins] = { "Low EA", "High EA" };
  TString EAstring[2] = { "loEA", "hiEA" };
  int EAmarker[nEAbins] = { 24, 20 };

  TString emw[nEtaBins] = { "east", "mid", "west" };

  const double eastArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  double area[nEtaBins] = { eastArea, midArea, westArea };
    
  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  TString dirName = "plots/prelim";


  
  double nChg_te[nEtaBins][nEAbins][nPtBins];
  double nChg_te_err[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt[nEtaBins][nEAbins][nPtBins];
  double te_meanChgPt_err[nEtaBins][nEAbins][nPtBins];

  TFile *inFile = new TFile("outFile.root","READ");

  TH1D *hUE[nEtaBins][nEAbins][nPtBins];
  
  for ( int p=0; p<nPtBins; ++p) {
    for ( int e=0; e<nEtaBins; ++e) {
      for (int a=nEAbins-1; a>=0; --a) {
	name = EAstring[a] + "_hUE_te_" + emw[e] + "Eta" + ptBinName[p];
	hUE[e][a][p]=(TH1D*)inFile->Get(name);
	nChg_te[e][a][p] = hUE[e][a][p]->Integral()/area[e];
	nChg_te_err[e][a][p] = hUE[e][a][p]->GetMeanError();
	te_meanChgPt[e][a][p] = hUE[e][a][p]->GetMean();
	te_meanChgPt_err[e][a][p] = hUE[e][a][p]->GetMeanError();
	cout<<nChg_te[e][a][p]<<"   "<<te_meanChgPt[e][a][p]<<endl;
      }
    }
  }

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  double binEdge[nEtaBins+1] = { -0.5, 0.5, 1.5, 2.5 };

  //  double ptBinEdge[nPtBins+1] = { 9.5, 14.5, 19.5, 29.5 };
  double ptBinEdge[nPtBins+1] = { 10, 15, 20, 20 };
  double yEdge[2] = {0.5,1.8};
  int startBin[nPtBins] = { 11, 30, 59 };
  

  TH2D *hscale0 = new TH2D("hscale0",";p_{T,lead}^{reco} [GeV/#it{c}];#LT #frac{d#it{N}_{ch}}{d#eta d#phi} #GT", 10,10,30,10,0.5,1.8);
  hscale0->GetYaxis()->SetTitleOffset(1);
  hscale0->GetXaxis()->SetTitleOffset(1.0);
  hscale0->GetYaxis()->SetTitleSize(0.05);
  hscale0->GetXaxis()->SetLabelSize(0.04);
  hscale0->GetYaxis()->SetLabelSize(0.04);
  hscale0->GetXaxis()->SetTitleOffset(1.0);
  hscale0->GetXaxis()->SetTitleSize(0.05);

  hscale0->GetXaxis()->SetNdivisions(6);
  
  hscale0->SetName("");
  TCanvas *c0 = new TCanvas;  
  c0->SetMargin(0.125,0.05,0.125,0.05);
  
  TH1D *hCHGdNdetadphi_te[nEAbins][nEtaBins];
  // TH1D *hCHGdNdetadphi_te_sys[nEAbins][nEtaBins];
  hscale0->Draw();
  hscale0->SetName(" ");
  hscale0->SetLineWidth(0);
  hscale0->SetStats(0);

  TLegend *leg0;
  
  leg0 = new TLegend(0.45, 0.67, 0.95, 0.95,NULL,"brNDC");    // LEGEND 0
  leg0->SetBorderSize(0);  leg0->SetLineColor(1);  leg0->SetLineStyle(1);
  leg0->SetLineWidth(1);  leg0->SetFillColorAlpha(0,0.0);  leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);
  leg0->AddEntry((TObject*)0,"0-30% EA", "");
  leg0->AddEntry((TObject*)0,"70-90% EA", "");

  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hCHGdNdetadphi_te_" + EAbinName[a] + etaBinName[e];
      title = UEetaBinString[e];
      hCHGdNdetadphi_te[a][e] = new TH1D( name, title , 79,9.25,29.75 );

      hCHGdNdetadphi_te[a][e]->SetStats(0);
      hCHGdNdetadphi_te[a][e]->SetLineColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerColor( etaColor[e] );
      hCHGdNdetadphi_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hCHGdNdetadphi_te[a][e]->SetMarkerSize( 2 );

      leg0->AddEntry(hCHGdNdetadphi_te[a][e],title,"lpf");
      
      // name = "hCHGdNdetadphi_te_sys_" + EAbinName[a] + etaBinName[e];
      // hCHGdNdetadphi_te_sys[a][e] = new TH1D( name, "" , 79,9.25,29.75 );
      
      // hCHGdNdetadphi_te_sys[a][e]->SetStats(0);
      // hCHGdNdetadphi_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      // hCHGdNdetadphi_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      // hCHGdNdetadphi_te_sys[a][e]->SetMarkerSize( 0 );

      // if (a==0) {
      // 	hCHGdNdetadphi_te_sys[a][e]->SetFillStyle(3002);
      // 	hCHGdNdetadphi_te_sys[a][e]->SetFillColor( etaColor[e]);
      // }
      // else { hCHGdNdetadphi_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 ); }

      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int binno = startBin[p] + (2*e);
	
	double value = (double) ( nChg_te[e][a][p] / area[e] );
	hCHGdNdetadphi_te[a][e]->SetBinContent( binno, value );
	hCHGdNdetadphi_te[a][e]->SetBinError( binno, nChg_te_err[e][a][p] );

	// for (int i=-3; i<4; ++i) {
	//   hCHGdNdetadphi_te_sys[a][e]->SetBinContent( binno+i, value );
	//   hCHGdNdetadphi_te_sys[a][e]->SetBinError( binno+i, value*nChg_te_sys[e][a][p] );
	// }
	
      }
    }
  }

  // for (int a=nEAbins-1; a>=0; --a) {
  //   for ( int e=0; e<nEtaBins; ++e) {
  //     hCHGdNdetadphi_te_sys[a][e]->Draw("e2SAME");
  //   }
  // }
  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hCHGdNdetadphi_te[a][e]->Draw("lpfX0ESAME");
    }
  }


  
  
  leg0->Draw();
  name = dirName+"/CHGdNdetadphi_te_pt_systematics.pdf";
  c0->SaveAs(name,"PDF");
  //c0->Close();

 
  TCanvas *c1 = new TCanvas;
  c1->SetMargin(0.125,0.05,0.125,0.05);

  double yEdge1[2] = { 0.55, 0.85 };
  TH2D *hscale3 = new TH2D("hscale3",";p_{T,lead}^{reco} [GeV/#it{c}];#LT p_{T}^{ch} #GT [GeV/#it{c}]", 10,10,30,10,0.55,0.85);
  hscale3->GetYaxis()->SetTitleOffset(1.0);
  hscale3->GetYaxis()->SetTitleSize(0.05);
  hscale3->GetXaxis()->SetTitleOffset(1.0);
  hscale3->GetXaxis()->SetTitleSize(0.05);
  hscale3->GetXaxis()->SetLabelSize(0.04);
  hscale3->GetYaxis()->SetLabelSize(0.04);
  hscale3->SetName("");

  hscale3->GetXaxis()->SetNdivisions(6);

  TH1D *hChg_pt_te[nEAbins][nEtaBins];
  // TH1D *hChg_pt_te_sys[nEAbins][nEtaBins];
  hscale3->SetName(" ");
  hscale3->SetLineWidth(0);
  hscale3->SetStats(0);
  hscale3->Draw();


  for ( int e=0; e<nEtaBins; ++e) {
    for (int a=nEAbins-1; a>=0; --a) {

      name = "hChg_MeanPt_te_" + EAbinName[a] + etaBinName[e];
      title = EAbinString[a] + "   "  + emw[e] + "   ";
      hChg_pt_te[a][e] = new TH1D( name, title , 79,9.25,29.75 );

      hChg_pt_te[a][e]->SetStats(0);
      hChg_pt_te[a][e]->SetLineColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerColor( etaColor[e] );
      hChg_pt_te[a][e]->SetMarkerStyle( EAmarker[a] );
      hChg_pt_te[a][e]->SetMarkerSize( 2 );
      
      // name = "hChg_MeanPt_te_sys_" + EAbinName[a] + etaBinName[e];
      // hChg_pt_te_sys[a][e] = new TH1D( name, "" , 79,9.25,29.75 );

      // hChg_pt_te_sys[a][e]->SetStats(0);
      // hChg_pt_te_sys[a][e]->SetMarkerColorAlpha( etaColor[e], 0.0 );
      // hChg_pt_te_sys[a][e]->SetLineColorAlpha( etaColor[e], 0.0 );
      // hChg_pt_te_sys[a][e]->SetMarkerSize( 0.0 );

      // if (a==0) {
      // 	hChg_pt_te_sys[a][e]->SetFillStyle(3002);
      // 	hChg_pt_te_sys[a][e]->SetFillColor( etaColor[e]);
      // }
      // else { hChg_pt_te_sys[a][e]->SetFillColorAlpha( etaColor[e], 0.2 ); }
      //hCHGdNdetadphi_te[a][e];
	
      for ( int p=0; p<nPtBins; ++p) {
	
	int binno = startBin[p] + (2*e);

	double value = te_meanChgPt[e][a][p];
	hChg_pt_te[a][e]->SetBinContent( binno, value );
	hChg_pt_te[a][e]->SetBinError( binno, te_meanChgPt_err[e][a][p] );

	// for (int i=-2; i<3; ++i) {
	//   hChg_pt_te_sys[a][e]->SetBinContent( binno+i, value );
	//   hChg_pt_te_sys[a][e]->SetBinError( binno+i, value*nChg_te_sys[e][a][p] );
	// }
		
      }
    }
  }

  // for (int a=nEAbins-1; a>=0; --a) {
  //   for ( int e=0; e<nEtaBins; ++e) {
  //     hChg_pt_te_sys[a][e]->Draw("E2SAME");
  //   }
  // }
  for (int a=nEAbins-1; a>=0; --a) {
    for ( int e=0; e<nEtaBins; ++e) {
      hChg_pt_te[a][e]->Draw("lpfX0ESAME");
    }
  }  
  leg0->Draw();


  name = dirName+"/Chg_MeanPt_te_pt_systematics.pdf";
  c1->SaveAs(name,"PDF");      

}
