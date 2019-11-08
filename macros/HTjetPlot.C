//  HTmonojetPlot.C
//  Veronica Verkest		July 14, 2019

void HTjetPlot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  // const int acn = 3;
  // TString BackgroundChargeBias[acn] = { "allBG", "chgBG", "neuBG" };
  // TString JetChargeBias[acn] = { "allJets", "chgJets", "neuJets" };
  TString fileName, name, title, Ndj, avg, sigma, jdj, ea;			double scale;
  
  jdj = "jets";
  // jdj = "dijets";
  ea = "HiEA";
  // ea = "LoEA";

  TString slash = "/";
  
  TString path = "HTdijets/HiEA/";
  // TString path = "HTdijets/MidEA/";
  // TString path = "HTdijets/LoEA/";
  
  fileName = "out/" + path + "pAu_2015_HT" + jdj + ".root";
  TFile* inFile = new TFile( fileName, "READ" );

  path = ea + slash + jdj + slash;

  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");
  hTowersPerEvent->Scale(1.0/hTowersPerEvent->Integral("WIDTH"));
  
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");
  hPrimaryPerEvent->Scale(1.0/hPrimaryPerEvent->Integral("WIDTH"));
  
  TH3D *hAllJetsPtEtaPhi = (TH3D*) inFile->Get("hAllJetsPtEtaPhi");
  
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");
  hTowersVsRho->Scale(1.0/hTowersVsRho->Integral("WIDTH"));
  hTowersVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hTowersVsRho->GetYaxis()->SetRangeUser(0.0,200.0);

  //TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  
  TH3D *hLeadPtEtaPhi = (TH3D*) inFile->Get("hLeadPtEtaPhi");

  TH2D *hLeadPtVsRho = (TH2D*)hPt_UE_BBCsumE->Project3D("XY");
  hLeadPtVsRho->GetYaxis()->SetRangeUser(0.0,10.0);

  // TH2D *hChgVsNeuBG = (TH2D*) inFile->Get("hChgVsNeuBG");				hChgVsNeuBG->Scale(1.0/hChgVsNeuBG->Integral("WIDTH"));
  // TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");	       		hTowersPerRun->Scale(1.0/hTowersPerRun->Integral("WIDTH"));
  // TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");		       	hPrimaryPerRun->Scale(1.0/hPrimaryPerRun->Integral("WIDTH"));
  // TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");	hnPrimaryVSnTowers->Scale(1.0/hnPrimaryVSnTowers->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");				hPrimaryVsBBC->Scale(1.0/hPrimaryVsBBC->Integral("WIDTH"));
  // TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");		       	hPrimaryVsGlobal->Scale(1.0/hPrimaryVsGlobal->Integral("WIDTH"));
  // TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");				hGlobalVsBBC->Scale(1.0/hGlobalVsBBC->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");			hPrimaryVsBBCE->Scale(1.0/hPrimaryVsBBCE->Integral("WIDTH"));
  // TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");				hGlobalVsBBCE->Scale(1.0/hGlobalVsBBCE->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");		hPrimaryVsBBCsumE->Scale(1.0/hPrimaryVsBBCsumE->Integral("WIDTH"));
  // TH2D *hTowersVsBBCsumE = (TH2D*) inFile->Get("hTowersVsBBCsumE");		hTowersVsBBCsumE->Scale(1.0/hTowersVsBBCsumE->Integral("WIDTH"));
  // TH2D *hPrimaryVsRho = (TH2D*) inFile->Get("hPrimaryVsRho");				hPrimaryVsRho->Scale(1.0/hPrimaryVsRho->Integral("WIDTH"));
  // TH2D *hGlobalVsRho = (TH2D*) inFile->Get("hGlobalVsRho");				hGlobalVsRho->Scale(1.0/hGlobalVsRho->Integral("WIDTH"));
  // TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  // TH3D *hBG = (TH3D*) inFile->Get("hBG");
  // TH3D *hPt_UE_RefMult = (TH3D*) inFile->Get("hPt_UE_RefMult");
  // TH3D *hAllJetsPtRhoEta = (TH3D*) inFile->Get("hAllJetsPtRhoEta");
  // TH3D *hPartPtEtaPhi = (TH3D*) inFile->Get("hPartPtEtaPhi");
  // TH3D *hAllPtEtaPhi = (TH3D*) inFile->Get("hAllPtEtaPhi");
  // TH3D *hPartPtDEtaDPhi = (TH3D*) inFile->Get("hPartPtDEtaDPhi");
  // TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  // TH3D *hSubPtEtaPhi = (TH3D*) inFile->Get("hSubPtEtaPhi");
  
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
  
  // TCanvas * cBG0 = new TCanvas( "cBG0" , "" ,0 ,23 ,1280 ,700 );
  // cBG0->SetTopMargin(0.4);
  // TPaveText *cTitle = new TPaveText(0.345843,.881306,0.655712,.980712,"NB");
  // cTitle->AddText( ea + jdj + ": Charged Background" );
  // cTitle->SetFillStyle(0);
  // cTitle->SetLineWidth(0);
  // cTitle->SetTextAlign(21);
  // cTitle->Draw();
    
  // cBG0->Divide(nEtaBins,2,0,0);

  TH1D *hRho[nPtBins];		TH1D *hLeadEta[nPtBins];		TH1D *hSubEta[nPtBins];

  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS

  name = "Lead Jet #eta by p_{T};#eta;";
  TH2D *sLeadEta = new TH2D( "sLeadEta", name, 40,-1.0,1.0, 10,0.0, 2.0 );
  sLeadEta->SetStats(0);
  
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(3);
  leg0->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg0->AddEntry((TObject*)0,"#bf{# of Jets}", "");
  leg0->AddEntry((TObject*)0,"#bf{<#eta>}", "");

  sLeadEta->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadEta" + LeadPtBinName[i];      title = LeadPtBinString[i];
    hLeadPtEtaPhi->GetXaxis()->SetRangeUser(LeadPtBinLo[i], LeadPtBinHi[i]);
    hLeadEta[i] = (TH1D*) hLeadPtEtaPhi->Project3D( "Y" );       // PROJECT
    hLeadEta[i]->SetNameTitle(name,title);
    hLeadEta[i]->SetStats(0);
    
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
  
  c0->SetLogy();
  
  name = /*JetChargeString + ": " + BackgroundChargeString + " */ "Underlying Event by Lead Jet p_{T};#rho (GeV);";
  TH2D *sLeadPtVsRho = new TH2D( "sLeadPtVsRho", name, 50,0,15, 10,0.000001, 1.0 );
  sLeadPtVsRho->SetStats(0);

  TLegend *leg2 = new TLegend(0.6, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(4);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg2->AddEntry((TObject*)0,"#bf{# of Jets}", "");
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
    
  // hTowersPerEvent->Draw();		name = "plots/" + path + "TowersPerEvent.pdf";
  // // title = JetChargeString + ": " + hTowersPerEvent->GetTitle();		hTowersPerEvent->SetTitle( title );
  // c0->SaveAs( name, "PDF" );

  // hPrimaryPerEvent->Draw();		name = "plots/" + path + "PrimaryPerEvent.pdf";
  // // title = JetChargeString + ": " + hPrimaryPerEvent->GetTitle();		hPrimaryPerEvent->SetTitle( title );
  // c0->SaveAs( name, "PDF" );

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
  
  scale = hLeadPtVsRho->Integral("WIDTH");			hLeadPtVsRho->Scale(1.0/scale);			hLeadPtVsRho->GetZaxis()->SetRangeUser(0.00001,1);
  //title = JetChargeString + ", " + BackgroundChargeString + ": " + hLeadPtVsRho->GetTitle();			hLeadPtVsRho->SetTitle( title );
  hLeadPtVsRho->Draw("COLZ");			name = "plots/" + path + "LeadPtVsRho.pdf";				c0->SaveAs( name, "PDF" );
  hLeadPtVsRho->Scale(scale);


  inFile->Close();

}
