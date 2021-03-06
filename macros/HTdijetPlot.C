//  HTJP2Plot.C
//  Veronica Verkest		July 14, 2019

void HTdijetPlot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  // const int acn = 3;
  // TString BackgroundChargeBias[acn] = { "allBG", "chgBG", "neuBG" };
  // TString JetChargeBias[acn] = { "allJets", "chgJets", "neuJets" };
  TString fileName, name, title, Ndj, avg, sigma;			double scale;

  // TString BackgroundChargeBias = "allBG";		TString BackgroundChargeString = "All Background";			//  [  "allBG",  "chgBG",  "neuBG"  ]
  // TString BackgroundChargeBias = "chgBG";		TString BackgroundChargeString = "Charged Background";
  // TString BackgroundChargeBias = "neuBG";		TString BackgroundChargeString = "Neutral Background";
  
  // TString JetChargeBias = "allJets";				TString JetChargeString = "All Jets";						//  [ "allJets", "chgJets", "neuJets" ]
  // TString JetChargeBias = "chgJets";				TString JetChargeString = "Charged Jets";
  // TString JetChargeBias = "neuJets";				TString JetChargeString = "Neutral Jets";

  // TString path = "HTjets/allTowers/";   // good/    allTowers/    towersRemoved/
  TString path = "HTdijets/";
  
  //fileName = "out/" + path + "pAu_2015_200_" + BackgroundChargeBias + "_" + JetChargeBias + ".root";
  fileName = "out/HTdijets/pAu_2015_HTdijets.root";
  //fileName = "out/HTJP2dijets/pAu_dijets_good.root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");
  hTowersPerEvent->Scale(1.0/hTowersPerEvent->Integral("WIDTH"));
  
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");
  hPrimaryPerEvent->Scale(1.0/hPrimaryPerEvent->Integral("WIDTH"));
  
  TH3D *hAllJetsPtEtaPhi = (TH3D*) inFile->Get("hAllJetsPtEtaPhi");
  
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");
  hTowersVsRho->Scale(1.0/hTowersVsRho->Integral("WIDTH"));
  hTowersVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hTowersVsRho->GetYaxis()->SetRangeUser(0.0,200.0);

  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  
  TH3D *hLeadPtEtaPhi = (TH3D*) inFile->Get("hLeadPtEtaPhi");

  TH2D *hLeadPtVsRho = (TH2D*)hPt_UE_BBCE->Project3D("XY");
  hLeadPtVsRho->GetYaxis()->SetRangeUser(0.0,10.0);
  
  // const int nPtBins = 4;
  // double LeadPtBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0 };
  // double LeadPtBinHi[nPtBins] = { 15.0, 20.0, 30.0, 100.0 };
  // TString LeadPtBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", ">30 GeV" };
  // TString LeadPtBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30" };
  // double SubPtBinLo[nPtBins] = { 0.0, 10.0, 15.0, 20.0 };
  // double SubPtBinHi[nPtBins] = { 10.0, 15.0, 20.0, 100.0 };
  // TString SubPtBinString[nPtBins] = { "<10 GeV", "10-15 GeV", "15-20 GeV", ">20 GeV" };
  // TString SubPtBinName[nPtBins] = { "_02_10","_10_15", "_15_20", "_20" };
  // int color[nPtBins] = { 879, 856, 796, 896 };
  // int marker[nPtBins] = { 33, 22, 21, 20 };

  const int nPtBins = 3;
  double LeadPtBinLo[nPtBins] = { 10.0, 15.0, 20.0 };
  double LeadPtBinHi[nPtBins] = { 15.0, 20.0, 30.0 };
  TString LeadPtBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV" };
  TString LeadPtBinName[nPtBins] = { "_10_15", "_15_20", "_20_30" };
  double SubPtBinLo[nPtBins] = { 5, 7.5, 10.0 };
  double SubPtBinHi[nPtBins] = { 7.5, 10.0, 30.0 };
  TString SubPtBinString[nPtBins] = { "5-7.5 GeV", "7.5-10 GeV", "10-30 GeV" };
  TString SubPtBinName[nPtBins] = { "_05_075","_075_10", "_10_30" };
  int color[nPtBins] = { 879, 856, 896 };
  int marker[nPtBins] = { 33, 21, 20 };
  
  TH1D *hRho[nPtBins];		TH1D *hLeadEta[nPtBins];		TH1D *hSubEta[nPtBins];

  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS

  name = "Lead Jet #eta by p_{T};#eta;";
  TH2D *sLeadEta = new TH2D( "sLeadEta", name, 40,-1.0,1.0, 10,0.0, 2.0 );
  sLeadEta->SetStats(0);
  
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(3);
  leg0->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg0->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg0->AddEntry((TObject*)0,"#bf{<#eta>}", "");

  sLeadEta->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadEta" + LeadPtBinName[i];      title = LeadPtBinString[i];
    hLeadPtEtaPhi->GetXaxis()->SetRangeUser(LeadPtBinLo[i], LeadPtBinHi[i]);
    hLeadEta[i] = (TH1D*) hLeadPtEtaPhi->Project3D( "Y" );       // PROJECT
    hLeadEta[i]->SetNameTitle(name,title);
    hLeadEta[i]->SetStats(0);

    for (int j=1;j<15;++j){double val=hLeadEta[i]->GetBinContent(j); hLeadEta[i]->SetBinContent(j,val/0.3);}
    for (int j=15;j<27;++j){double val=hLeadEta[i]->GetBinContent(j); hLeadEta[i]->SetBinContent(j,val/0.6);}
    for (int j=27;j<40;++j){double val=hLeadEta[i]->GetBinContent(j); hLeadEta[i]->SetBinContent(j,val/0.3);}    
    
    hLeadEta[i]->Scale( 1./hLeadEta[i]->Integral("WIDTH") );                     // NORMALIZE
    hLeadEta[i]->SetLineColor( color[i] );    hLeadEta[i]->SetMarkerStyle( marker[i] );    hLeadEta[i]->SetMarkerColor( color[i] );
    hLeadEta[i]->Draw("SAME");                                                    // DRAW
    Ndj = ""; avg = "";    Ndj += (int) hLeadEta[i]->GetEntries();
    avg += hLeadEta[i]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    leg0->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg0->AddEntry((TObject*)0,Ndj, "");    leg0->AddEntry((TObject*)0,avg, "");
  }

  leg0->Draw();  c0->Modified();  c0->cd();  c0->SetSelected(c0);
  name = "plots/" + path + "LeadEtaDist.pdf";
  c0->SaveAs( name , "PDF");
  
  hLeadPtEtaPhi->GetXaxis()->SetRangeUser(0.0,100.0);


  // name = JetChargeString + ": Sub Jet #eta by p_{T};#eta;";
  // TH2D *sSubEta = new TH2D( "sSubEta", name, 40,-1.0,1.0, 10,0.0, 2.0 );
  // sSubEta->SetStats(0);
  
  // TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  // leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  // leg1->SetNColumns(3);
  // leg1->AddEntry((TObject*)0,"#bf{p_{T}^{Sub} (GeV)}", "");
  // leg1->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  // leg1->AddEntry((TObject*)0,"#bf{<#eta>}", "");

  // sSubEta->Draw();
  // for ( int i=0; i<nPtBins; ++i ) {
  //   name = "SubEta" + SubPtBinName[i];      title = SubPtBinString[i];
  //   hSubPtEtaPhi->GetXaxis()->SetRangeUser(SubPtBinLo[i], SubPtBinHi[i]);
  //   hSubEta[i] = (TH1D*) hSubPtEtaPhi->Project3D( "Y" );       // PROJECT
  //   hSubEta[i]->SetNameTitle(name,title);
  //   hSubEta[i]->SetStats(0);
  //   hSubEta[i]->Scale( 1./hSubEta[i]->Integral("WIDTH") );                     // NORMALIZE
  //   hSubEta[i]->SetLineColor( color[i] );    hSubEta[i]->SetMarkerStyle( marker[i] );    hSubEta[i]->SetMarkerColor( color[i] );
  //   hSubEta[i]->Draw("SAME");                                                    // DRAW
  //   Ndj = ""; avg = "";    Ndj += hSubEta[i]->GetEntries();
  //   avg += hSubEta[i]->GetMean(1);                                           // 1 denotes x-axis
  //   avg = avg(0,6);
  //   leg1->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
  //   leg1->AddEntry((TObject*)0,Ndj, "");    leg1->AddEntry((TObject*)0,avg, "");
  // }

  // leg1->Draw();  c0->Modified();  c0->cd();  c0->SetSelected(c0);
  // name = "plots/" + path + "SubEtaDist.pdf";
  // c0->SaveAs( name , "PDF");
  
  
  c0->SetLogy();
  
  name = /*JetChargeString + ": " + BackgroundChargeString + " */ "Underlying Event by Lead Jet p_{T};#rho (GeV);";
  TH2D *sLeadPtVsRho = new TH2D( "sLeadPtVsRho", name, 50,0,15, 10,0.000001, 1.0 );
  sLeadPtVsRho->SetStats(0);

  TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(4);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg2->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#sigma> (GeV)}", "");

  sLeadPtVsRho->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadPtVsRho" + LeadPtBinName[i];      title = LeadPtBinString[i];
    hLeadPtVsRho->GetYaxis()->SetRangeUser(LeadPtBinLo[i], LeadPtBinHi[i]);
    hRho[i] = hLeadPtVsRho->ProjectionX( name );       // PROJECT
    hRho[i]->SetStats(0);
    hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
    hRho[i]->SetLineColor( color[i] );    hRho[i]->SetMarkerStyle( marker[i] );    hRho[i]->SetMarkerColor( color[i] );
    hRho[i]->Draw("SAME");                                                    // DRAW
    Ndj = ""; avg = ""; sigma="";
    Ndj += (int) hRho[i]->GetEntries();
    avg += hRho[i]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    sigma += hRho[i]->GetStdDev(1);
    sigma = sigma(0,5);
    leg2->AddEntry( name, title, "lpf" );                            // ADD TO LEGEND
    leg2->AddEntry((TObject*)0,Ndj, "");
    leg2->AddEntry((TObject*)0,avg, "");
    leg2->AddEntry((TObject*)0,sigma, "");
  }

  leg2->Draw();  c0->Modified();  c0->cd();  c0->SetSelected(c0);
  name = "plots/" + path + "RhoByLeadPt.pdf";
  c0->SaveAs( name , "PDF");

  hLeadPtVsRho->GetXaxis()->SetRangeUser(0.0, 10.0);
  hLeadPtVsRho->GetYaxis()->SetRangeUser(0.0, 100.0);
  

  
  hTowersPerEvent->Draw();		name = "plots/" + path + "TowersPerEvent.pdf";
  // title = JetChargeString + ": " + hTowersPerEvent->GetTitle();		hTowersPerEvent->SetTitle( title );
  c0->SaveAs( name, "PDF" );

  hPrimaryPerEvent->Draw();		name = "plots/" + path + "PrimaryPerEvent.pdf";
  // title = JetChargeString + ": " + hPrimaryPerEvent->GetTitle();		hPrimaryPerEvent->SetTitle( title );
  c0->SaveAs( name, "PDF" );

  TH1D *hAllJetsPt = (TH1D*) hAllJetsPtEtaPhi->ProjectionX();		hAllJetsPt->Scale(1.0/hAllJetsPt->Integral());
  title = "Inclusive Jet p_{T};p_{T} (GeV)";			hAllJetsPt->SetTitle(title);
  hAllJetsPt->Draw();					name = "plots/" + path + "AllJetsPt.pdf";											c0->SaveAs( name, "PDF" );

  c0->SetLogy(0);		c0->SetLogz();

  // hTowersVsBBCsumE->Draw("COLZ");
  // name = "plots/" + path + "TowersVsBBCsumE.pdf";			       	c0->SaveAs( name, "PDF" );

  TH2D *hAllJetsEtaPhi = (TH2D*) hAllJetsPtEtaPhi->Project3D("ZY");		hAllJetsEtaPhi->Scale(1.0/hAllJetsEtaPhi->Integral());
  title = "Inclusive Jet #eta-#phi;#eta;#phi";			hAllJetsEtaPhi->SetTitle( title );			hAllJetsEtaPhi->GetZaxis()->SetRangeUser(0.000001,0.01);
  hAllJetsEtaPhi->Draw("COLZ");			name = "plots/" + path + "AllJetsEtaPhi.pdf";			       	c0->SaveAs( name, "PDF" );
  
  TH2D *hLeadEtaPhi = (TH2D*) hLeadPtEtaPhi->Project3D("ZY");			hLeadEtaPhi->Scale(1.0/hLeadEtaPhi->Integral());
  title = "Lead Jet #eta-#phi;#eta;#phi";				hLeadEtaPhi->SetTitle( title );				hLeadEtaPhi->GetZaxis()->SetRangeUser(0.00001,0.01);
  hLeadEtaPhi->Draw("COLZ");			name = "plots/" + path + "LeadEtaPhi.pdf";				       	c0->SaveAs( name, "PDF" );
  
  // TH2D *hAllEtaPhi = (TH2D*) hAllPtEtaPhi->Project3D("ZY");			hAllEtaPhi->Scale(1.0/hAllEtaPhi->Integral());
  // title = JetChargeString + ": All Particles #eta-#phi;#eta;#phi";			hAllEtaPhi->SetTitle( title );				hAllEtaPhi->GetZaxis()->SetRangeUser(0.000001,0.01);
  // hAllEtaPhi->Draw("COLZ");			name = "plots/" + path + "AllEtaPhi.pdf";			  	  	       	c0->SaveAs( name, "PDF" );
    
  // TH2D *hPartEtaPhi = (TH2D*) hPartPtEtaPhi->Project3D("ZY");			hPartEtaPhi->Scale(1.0/hPartEtaPhi->Integral());
  // title = JetChargeString + ", " + BackgroundChargeString + ": Background Particles #eta-#phi;#eta;#phi";		hPartEtaPhi->SetTitle( title );
  // hPartEtaPhi->GetZaxis()->SetRangeUser(0.000001,0.01);
  // hPartEtaPhi->Draw("COLZ");			name = "plots/" + path + "PartEtaPhi_" + BackgroundChargeBias + ".pdf";			  	     	    	c0->SaveAs( name, "PDF" );

  // if ( BackgroundChargeBias == "allBG" ) {
  //   hChgVsNeuBG->GetZaxis()->SetRangeUser(0.00001,1);		title = JetChargeString + ": " + hChgVsNeuBG->GetTitle();		hChgVsNeuBG->SetTitle( title );
  //   hChgVsNeuBG->Draw("COLZ");			name = "plots/" + path + "ChgVsNeuBG.pdf";				c0->SaveAs( name, "PDF" );
  // }
  
  // hnPrimaryVSnTowers->GetZaxis()->SetRangeUser(0.00001,0.01);
  // title = JetChargeString + ": " + hnPrimaryVSnTowers->GetTitle();			hnPrimaryVSnTowers->SetTitle( title );
  // hnPrimaryVSnTowers->Draw("COLZ");		name = "plots/" + path + "nPrimaryVSnTowers.pdf";			c0->SaveAs( name, "PDF" );

  // hPrimaryVsBBC->GetZaxis()->SetRangeUser(0.00001,1);
  // hPrimaryVsBBC->Draw("COLZ");			name = "plots/" + path + "PrimaryVsBBC_" + BackgroundChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hPrimaryVsGlobal->GetZaxis()->SetRangeUser(0.00001,0.01);
  // title = JetChargeString + ": " + hPrimaryVsGlobal->GetTitle();			hPrimaryVsGlobal->SetTitle( title );
  // hPrimaryVsGlobal->Draw("COLZ");		name = "plots/" + path + "PrimaryVsGlobal.pdf";			c0->SaveAs( name, "PDF" );

  // hPrimaryVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // title = JetChargeString + ", " + BackgroundChargeString + ": " + hPrimaryVsRho->GetTitle();			hPrimaryVsRho->SetTitle( title );
  // hPrimaryVsRho->Draw("COLZ");			name = "plots/" + path + "PrimaryVsRho_" + BackgroundChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hGlobalVsRho->GetZaxis()->SetRangeUser(0.00001,0.01);
  // title = JetChargeString + ", " + BackgroundChargeString + ": " + hGlobalVsRho->GetTitle();			hGlobalVsRho->SetTitle( title );
  // hGlobalVsRho->Draw("COLZ");			name = "plots/" + path + "GlobalVsRho_" + BackgroundChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  // hTowersVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  // // title = JetChargeString + ", " + BackgroundChargeString + ": " + hTowersVsRho->GetTitle();			hTowersVsRho->SetTitle( title );
  // hTowersVsRho->Draw("COLZ");			name = "plots/" + path + "TowersVsRho_" + BackgroundChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  scale = hLeadPtVsRho->Integral("WIDTH");			hLeadPtVsRho->Scale(1.0/scale);			hLeadPtVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  //title = JetChargeString + ", " + BackgroundChargeString + ": " + hLeadPtVsRho->GetTitle();			hLeadPtVsRho->SetTitle( title );
  hLeadPtVsRho->Draw("COLZ");			name = "plots/" + path + "LeadPtVsRho.pdf";				c0->SaveAs( name, "PDF" );
  hLeadPtVsRho->Scale(scale);


  inFile->Close();

}
