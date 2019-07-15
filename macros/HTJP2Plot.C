//  HTJP2Plot.C
//  Veronica Verkest		July 14, 2019

void HTJP2Plot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  // const int acn = 3;
  // TString BackgroundChargeBias[acn] = { "allBG", "chgBG", "neuBG" };
  // TString JetChargeBias[acn] = { "allJets", "chgJets", "neuJets" };
  TString fileName, name, title, path, Ndj, avg;			double scale;

  TString BackgroundChargeBias = "allBG";     //  [  "allBG",  "chgBG",  "neuBG"  ]
  TString JetChargeBias = "allJets";                  //  [ "allJets", "chgJets", "neuJets" ]

  fileName = "out/HTJP2/pAu_2015_200_" + BackgroundChargeBias + "_" + JetChargeBias + ".root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");			hTowersPerEvent->Scale(1.0/hTowersPerEvent->Integral("WIDTH"));
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");			hPrimaryPerEvent->Scale(1.0/hPrimaryPerEvent->Integral("WIDTH"));
  TH2D *hChgVsNeuBG = (TH2D*) inFile->Get("hChgVsNeuBG");				hChgVsNeuBG->Scale(1.0/hChgVsNeuBG->Integral("WIDTH"));
  TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");	       		hTowersPerRun->Scale(1.0/hTowersPerRun->Integral("WIDTH"));
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");		       	hPrimaryPerRun->Scale(1.0/hPrimaryPerRun->Integral("WIDTH"));
  TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");	hnPrimaryVSnTowers->Scale(1.0/hnPrimaryVSnTowers->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");				hPrimaryVsBBC->Scale(1.0/hPrimaryVsBBC->Integral("WIDTH"));
  TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");		       	hPrimaryVsGlobal->Scale(1.0/hPrimaryVsGlobal->Integral("WIDTH"));
  // TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");				hGlobalVsBBC->Scale(1.0/hGlobalVsBBC->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");			hPrimaryVsBBCE->Scale(1.0/hPrimaryVsBBCE->Integral("WIDTH"));
  // TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");				hGlobalVsBBCE->Scale(1.0/hGlobalVsBBCE->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");		hPrimaryVsBBCsumE->Scale(1.0/hPrimaryVsBBCsumE->Integral("WIDTH"));
  // TH2D *hTowersVsBBCsumE = (TH2D*) inFile->Get("hTowersVsBBCsumE");		hTowersVsBBCsumE->Scale(1.0/hTowersVsBBCsumE->Integral("WIDTH"));
  TH2D *hPrimaryVsRho = (TH2D*) inFile->Get("hPrimaryVsRho");				hPrimaryVsRho->Scale(1.0/hPrimaryVsRho->Integral("WIDTH"));
  TH2D *hGlobalVsRho = (TH2D*) inFile->Get("hGlobalVsRho");				hGlobalVsRho->Scale(1.0/hGlobalVsRho->Integral("WIDTH"));
  TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");					hSubEtaPhi->Scale(1.0/hSubEtaPhi->Integral("WIDTH"));
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");				hTowersVsRho->Scale(1.0/hTowersVsRho->Integral("WIDTH"));
  TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  TH3D *hBG = (TH3D*) inFile->Get("hBG");
  TH3D *hNEUTRAL = (TH3D*) inFile->Get("hBG");
  TH3D *hPt_UE_RefMult = (TH3D*) inFile->Get("hPt_UE_RefMult");
  TH3D *hCHARGED = (TH3D*) inFile->Get("hPt_UE_RefMult");
  TH3D *hAllJetsPtRhoEta = (TH3D*) inFile->Get("hAllJetsPtRhoEta");
  TH3D *hPartPtEtaPhi = (TH3D*) inFile->Get("hPartPtEtaPhi");
  TH3D *hAllPtEtaPhi = (TH3D*) inFile->Get("hAllPtEtaPhi");
  TH3D *hPartPtDEtaDPhi = (TH3D*) inFile->Get("hPartPtDEtaDPhi");
  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  TH3D *hLeadPtEtaPhi = (TH3D*) inFile->Get("hLeadPtEtaPhi");
  TH3D *hAllJetsPtEtaPhi = (TH3D*) inFile->Get("hAllJetsPtEtaPhi");

  hPrimaryVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hPrimaryVsRho->GetYaxis()->SetRangeUser(0.0,125.0);
  hPrimaryVsGlobal->GetXaxis()->SetRangeUser(0.0,2000.0);
  hPrimaryVsGlobal->GetYaxis()->SetRangeUser(0.0,80.0);
  hTowersVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hChgVsNeuBG->GetXaxis()->SetRangeUser(0.0,10.0);
  hChgVsNeuBG->GetYaxis()->SetRangeUser(0.0,10.0);
  hGlobalVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hnPrimaryVSnTowers->GetXaxis()->SetRangeUser(0.0,100.0);
  hnPrimaryVSnTowers->GetYaxis()->SetRangeUser(0.0,100.0);
  hLeadPtVsRho->GetYaxis()->SetRangeUser(0.0,10.0);
  
  const int nPtBins = 4;
  double ptBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0 };
  double ptBinHi[nPtBins] = { 15.0, 20.0, 30.0, 100.0 };
  TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", ">30 GeV" };
  TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30" };
  int color[nPtBins] = { 879, 856, 796, 896 };
  int marker[nPtBins] = { 33, 22, 21, 20 };

  TH1D * hRho[nPtBins];		TH1D * hLeadEta[nPtBins];

  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS

  TH2D *sLeadEta = new TH2D( "sLeadEta", "Lead Jet #eta by p_{T};#eta;", 40,-1.0,1.0, 10,0.0, 1.0 );
  sLeadEta->SetStats(0);
  
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(3);
  leg0->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg0->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg0->AddEntry((TObject*)0,"#bf{<#eta>}", "");

  sLeadEta->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadEta" + ptBinName[i];      title = ptBinString[i];
    hLeadPtEtaPhi->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hLeadEta[i] = (TH1D*) hLeadPtEtaPhi->Project3D( "Z" );       // PROJECT
    hLeadEta[i]->SetNameTitle(name,title);
    hLeadEta[i]->SetStats(0);
    hLeadEta[i]->Scale( 1./hLeadEta[i]->Integral("WIDTH") );                     // NORMALIZE
    hLeadEta[i]->SetLineColor( color[i] );    hLeadEta[i]->SetMarkerStyle( marker[i] );    hLeadEta[i]->SetMarkerColor( color[i] );
    hLeadEta[i]->Draw("SAME");                                                    // DRAW
    Ndj = ""; avg = "";    Ndj += hLeadEta[i]->GetEntries();
    avg += hLeadEta[i]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    leg0->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg0->AddEntry((TObject*)0,Ndj, "");    leg0->AddEntry((TObject*)0,avg, "");
  }

  leg0->Draw();  c0->Modified();  c0->cd();  c0->SetSelected(c0);
  path = "plots/HTJP2/LeadEtaDist_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";
  c0->SaveAs( path , "PDF");
  
  hLeadPtEtaPhi->GetXaxis()->SetRangeUser(0.0,100.0);
  
  // c0->SetLogy();
  
  // TH2D *sLeadPtVsRho = new TH2D( "sLeadPtVsRho", "Underlying Event by Lead Jet p_{T};#rho (GeV);", 50,0,15, 10,0.000001, 1.0 );
  // sLeadPtVsRho->SetStats(0);

  // TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  // leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  // leg1->SetNColumns(3);
  // leg1->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  // leg1->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  // leg1->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");

  // sLeadPtVsRho->Draw();
  // for ( int i=0; i<nPtBins; ++i ) {
  //   name = "LeadPtVsRho" + ptBinName[i];      title = ptBinString[i];
  //   hLeadPtVsRho->GetYaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
  //   hRho[i] = hLeadPtVsRho->ProjectionX( name );       // PROJECT
  //   hRho[i]->SetStats(0);
  //   hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
  //   hRho[i]->SetLineColor( color[i] );    hRho[i]->SetMarkerStyle( marker[i] );    hRho[i]->SetMarkerColor( color[i] );
  //   hRho[i]->Draw("SAME");                                                    // DRAW
  //   Ndj = ""; avg = "";    Ndj += hRho[i]->GetEntries();
  //   avg += hRho[i]->GetMean(1);                                           // 1 denotes x-axis
  //   avg = avg(0,6);
  //   leg1->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
  //   leg1->AddEntry((TObject*)0,Ndj, "");    leg1->AddEntry((TObject*)0,avg, "");
  // }

  // leg1->Draw();  c0->Modified();  c0->cd();  c0->SetSelected(c0);
  // path = "plots/HTJP2/RhoByLeadPt_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";
  // c0->SaveAs( path , "PDF");

  // hLeadPtVsRho->GetXaxis()->SetRangeUser(0.0, 10.0);
  // hLeadPtVsRho->GetYaxis()->SetRangeUser(0.0, 100.0);
  

  
  // hTowersPerEvent->Draw();				name = "plots/HTJP2/TowersPerEvent_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  // hPrimaryPerEvent->Draw();				name = "plots/HTJP2/PrimaryPerEvent_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  // TH1D *hAllJetsPt = (TH1D*) hAllJetsPtEtaPhi->ProjectionX();		hAllJetsPt->Scale(1.0/hAllJetsPt->Integral());			hAllJetsPt->SetTitle("Inclusive Jet p_{T};p_{T} (GeV)");
  // hAllJetsPt->Draw();					name = "plots/HTJP2/AllJetsPt_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";					c0->SaveAs( name, "PDF" );

  // c0->SetLogy(0);		c0->SetLogz();
  
  // TH2D *hAllJetsEtaPhi = (TH2D*) hAllJetsPtEtaPhi->Project3D("ZY");		hAllJetsEtaPhi->Scale(1.0/hAllJetsEtaPhi->Integral());
  // hAllJetsEtaPhi->SetTitle("Inclusive Jet #eta-#phi;#eta;#phi");			hAllJetsEtaPhi->GetZaxis()->SetRangeUser(0.0000001,1);
  // hAllJetsEtaPhi->Draw("COLZ");			name = "plots/HTJP2/AllJetsEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			       	c0->SaveAs( name, "PDF" );
  
  // TH2D *hLeadEtaPhi = (TH2D*) hLeadPtEtaPhi->Project3D("ZY");			hLeadEtaPhi->Scale(1.0/hLeadEtaPhi->Integral());
  // hLeadEtaPhi->SetTitle("Lead Jet #eta-#phi;#eta;#phi");				hLeadEtaPhi->GetZaxis()->SetRangeUser(0.00001,1);
  // hLeadEtaPhi->Draw("COLZ");			name = "plots/HTJP2/LeadEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			       	c0->SaveAs( name, "PDF" );
  
  // TH2D *hAllEtaPhi = (TH2D*) hAllPtEtaPhi->Project3D("ZY");			hAllEtaPhi->Scale(1.0/hAllEtaPhi->Integral());
  // hAllEtaPhi->SetTitle("All Particles #eta-#phi;#eta;#phi");				hAllEtaPhi->GetZaxis()->SetRangeUser(0.000001,1);
  // hAllEtaPhi->Draw("COLZ");			name = "plots/HTJP2/AllEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			  	  	       	c0->SaveAs( name, "PDF" );
    
  // TH2D *hPartEtaPhi = (TH2D*) hPartPtEtaPhi->Project3D("ZY");			hPartEtaPhi->Scale(1.0/hPartEtaPhi->Integral());
  // hPartEtaPhi->SetTitle("Background Particles #eta-#phi;#eta;#phi");		hPartEtaPhi->GetZaxis()->SetRangeUser(0.000001,1);
  // hPartEtaPhi->Draw("COLZ");			name = "plots/HTJP2/PartEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			  	     	    	c0->SaveAs( name, "PDF" );
  
  // hChgVsNeuBG->GetZaxis()->SetRangeUser(0.00001,1);
  // hChgVsNeuBG->Draw("COLZ");			name = "plots/HTJP2/ChgVsNeuBG_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hnPrimaryVSnTowers->GetZaxis()->SetRangeUser(0.00001,1);
  // hnPrimaryVSnTowers->Draw("COLZ");		name = "plots/HTJP2/nPrimaryVSnTowers_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  // // hPrimaryVsBBC->GetZaxis()->SetRangeUser(0.00001,1);
  // // hPrimaryVsBBC->Draw("COLZ");			name = "plots/HTJP2/PrimaryVsBBC_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hPrimaryVsGlobal->GetZaxis()->SetRangeUser(0.00001,1);
  // hPrimaryVsGlobal->Draw("COLZ");		name = "plots/HTJP2/PrimaryVsGlobal_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  // hPrimaryVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // hPrimaryVsRho->Draw("COLZ");			name = "plots/HTJP2/PrimaryVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hGlobalVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // hGlobalVsRho->Draw("COLZ");			name = "plots/HTJP2/GlobalVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hSubEtaPhi->GetZaxis()->SetRangeUser(0.00001,1);
  // hSubEtaPhi->Draw("COLZ");				name = "plots/HTJP2/SubEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";					c0->SaveAs( name, "PDF" );

  // hTowersVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // hTowersVsRho->Draw("COLZ");			name = "plots/HTJP2/TowersVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // scale = hLeadPtVsRho->Integral("WIDTH");			hLeadPtVsRho->Scale(1.0/scale);			hLeadPtVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // hLeadPtVsRho->Draw("COLZ");			name = "plots/HTJP2/LeadPtVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );
  // hLeadPtVsRho->Scale(scale);


  inFile->Close();

}
