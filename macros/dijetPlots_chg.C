//  dijetPlots_chg.C
//  Veronica Verkest     May 22, 2019

void dijetPlots_chg() {
  
  const float pi = 3.141592;
  const double twopi = 2*3.14159265358979;
  TString Ndj, avg;      TString lpf = "lpf";
  double scale;
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  //TString dir = "HTjets/allTowers/chg/";
  TString dir = "HTjets/good/chg/";
  //TString dir = "HTjets/towersRemoved/chg/";
  TString path = "out/" + dir + "pAu_HT_dijets_CHG.root";
  TFile* inFile = new TFile( path, "READ" );
  
  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");
  TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");
  TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");
  TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");
  TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");
  TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");
  TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");
  TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");
  TH2D *hTowersVsBBCsumE = (TH2D*) inFile->Get("hTowersVsBBCsumE");
  TH2D *hLeadEtaPhi =(TH2D*) inFile->Get("hLeadEtaPhi");
  TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");
  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");  
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");
  TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH3D *hPartPtEtaPhi = (TH3D*) inFile->Get("hPartPtEtaPhi");
  TH3D *hAllPtEtaPhi = (TH3D*) inFile->Get("hAllPtEtaPhi");
  TH3D *hBG = (TH3D*) inFile->Get("hBG");

  hPt_UE_BBCsumE->GetZaxis()->SetRangeUser( 0.0, 80000.0 );
  // hPt_UE_BBCsumE->GetYaxis()->SetRangeUser( 0.0,10 );
  // hPt_UE_BBCE->GetYaxis()->SetRangeUser( 0.0,10 );
  hTowersVsRho->GetXaxis()->SetRangeUser(0.0,10.0);
  hTowersVsRho->GetYaxis()->SetRangeUser(0.0,250.0);
  
  TH2D *hscale0 = new TH2D( "hscale0", "Charged Underlying Event by Lead Jet p_{T};#rho (GeV);", 50,0,25, 10,0.000001, 1.0 );
  TH2D *hscale1 = new TH2D( "hscale1", "Charged Underlying Event vs. BBC East Rate", 140,0,7000000, 20,0,10 );
  TH2D *hscale2 = new TH2D( "hscale2", "Charged Underlying Event vs. BBC ADC East Sum", 150,0,100000, 20,0,10 );

  //  const int nPtBins = 5;
  const int nPtBins = 4;

  hscale0->SetStats(0);  hscale1->SetStats(0);  hscale2->SetStats(0);  
  TH1D * hRho[nPtBins];  TH2D * hUE_BBCE[nPtBins];  TH2D * hUE_BBCsumE[nPtBins];
  TProfile *pUE_BBCE[nPtBins];  TProfile *pUE_BBCsumE[nPtBins];
  
  // double ptBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0, 40.0 };
  // double ptBinHi[nPtBins] = { 15.0, 20.0, 30.0, 40.0, 100.0 };
  // TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", "30-40 GeV", ">40 GeV" };
  // TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30_40", "_40" };
  // int color[nPtBins] = { 879, 856, 835, 796, 896 };
  // int marker[nPtBins] = { 33, 34, 22, 21, 20 };
  // TString name, title;

  double ptBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0 };
  double ptBinHi[nPtBins] = { 15.0, 20.0, 30.0, 100.0 };
  TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", ">30 GeV" };
  TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30" };
  int color[nPtBins] = { 879, 856, 796, 896 };
  int marker[nPtBins] = { 33, 22, 21, 20 };
  TString name, title;
  
  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );              // CANVAS

  TH2D* hBGEtaPhi[nPtBins];
  TH1D* hBGEta[nPtBins];


  c0->SetLogz();
  for ( int i=0; i<nPtBins; ++i ) {
    c0->cd();
    hBG->GetXaxis()->SetRangeUser( ptBinLo[i], ptBinHi[i] );
    hBGEtaPhi[i] = (TH2D*) hBG->Project3D("YZ");
    name = "BGEtaPhi_" + ptBinName[i];
    title = "p_{T}^{lead}: " + ptBinString[i];
    hBGEtaPhi[i]->SetNameTitle( name, title );
    hBGEtaPhi[i]->Scale(1./hBGEtaPhi[i]->Integral("width"));
    hBGEtaPhi[i]->Draw("COLZ");
    name = "plots/" + dir + "BG_eta_phi_" + ptBinName[i] + ".pdf";
    c0->SaveAs( name ,"PDF");
    name = "BGEta_" + ptBinName[i];
    title = "p_{T}^{lead}: " + ptBinString[i];
    hBGEta[i] = (TH1D*) hBGEtaPhi[i]->ProjectionY();
    hBGEta[i]->SetNameTitle( name, title );
  }

  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );              // CANVAS

  TLegend *leg3 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg3->SetBorderSize(1);   leg3->SetLineColor(1);   leg3->SetLineStyle(1);   leg3->SetLineWidth(1);   leg3->SetFillColor(0);   leg3->SetFillStyle(1001);
  leg3->SetNColumns(2);  leg3->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");  leg3->AddEntry((TObject*)0,"#bf{Particle} #bf{#eta}", "");

  c1->cd();
  c1->SetLogy();
  TH2D *hscale5  = new TH2D( "hscale5", "Charged Background Particle #eta by Lead Jet p_{T};#eta;", 40,-1.0,1.0, 10,0.2,1.0 );
  hscale5->SetStats(0);
  hscale5->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    hBGEta[i]->Scale( 1./hBGEta[i]->Integral("width") );
    hBGEta[i]->SetLineColor( color[i] );     hBGEta[i]->SetMarkerColor( color[i] );     hBGEta[i]->SetMarkerStyle( marker[i] );
    hBGEta[i]->SetStats(0);
    hBGEta[i]->Draw("SAME");
    avg = "";     avg += hBGEta[i]->GetMean();
    avg = avg(0,8);
    title = "p_{T}^{lead}: " + ptBinString[i];
    leg3->AddEntry( hBGEta[i], title, "lpf" );    leg3->AddEntry((TObject*)0,avg, "");
  }
  leg3->Draw();        c1->Modified();        c1->cd();        c1->SetSelected(c1);
  
  path = "plots/" + dir + "bg_eta_byPtBin_logY.pdf";
  c1->SaveAs( path ,"PDF");

  c1->SetLogy(0);
  path = "plots/" + dir + "bg_eta_byPtBin.pdf";
  c1->SaveAs( path ,"PDF");


  hBG->GetXaxis()->SetRangeUser( 0, 100.0 );

  c0->SetLogz();
  c1->SetLogz();

  TH2D* hPartEtaPhi = (TH2D*) hPartPtEtaPhi->Project3D("YZ");
  hPartEtaPhi->SetTitle("Charged BG Particles #eta vs. #phi");
  hPartEtaPhi->Scale(1./hPartEtaPhi->Integral("WIDTH"));
  hPartEtaPhi->Draw("COLZ");
  path = "plots/" + dir + "bg_eta_phi.pdf";
  c1->SaveAs( path ,"PDF");

  TH2D* hAllEtaPhi = (TH2D*) hAllPtEtaPhi->Project3D("YZ");
  hAllEtaPhi->SetTitle("All Particles #eta vs. #phi");
  hAllEtaPhi->Scale(1./hAllEtaPhi->Integral("WIDTH"));
  hAllEtaPhi->Draw("COLZ");
  path = "plots/" + dir + "all_eta_phi.pdf";
  c1->SaveAs( path ,"PDF");
 
  for ( int i=0; i<nPtBins; ++i ) {
    c0->cd();
    name = "UEvsBBCE" + ptBinName[i];
    title = "Charged Underlying Event vs. BBC East Rate - p_{T}^{lead}: " + ptBinString[i];
    hPt_UE_BBCE->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hUE_BBCE[i] = (TH2D*)hPt_UE_BBCE->Project3D( "yz" );       // PROJECT
    hUE_BBCE[i]->GetYaxis()->SetRangeUser(0,10);
    scale = hUE_BBCE[i]->Integral("width");
    hUE_BBCE[i]->Scale( 1./scale );                     // NORMALIZE
    hUE_BBCE[i]->SetNameTitle(name,title);
    //hUE_BBCE[i]->Write();
    hUE_BBCE[i]->Draw("COLZ");
    //hUE_BBCE[i]->RebinX();
    hUE_BBCE[i]->SetLineColor(kRed);
    name = "pUE_BBCE" + ptBinName[i];
    pUE_BBCE[i] = hUE_BBCE[i]->ProfileX("",1,-1,"S");
    pUE_BBCE[i]->SetStats(0);
    pUE_BBCE[i]->SetName(name);
    pUE_BBCE[i]->Draw("SAME");
    name = "plots/" + dir + "UEvsBBCE" + ptBinName[i] + ".pdf";
    c0->SaveAs( name,"PDF");
    hUE_BBCE[i]->Scale( scale );                     // UN-NORMALIZE
    hUE_BBCE[i]->GetYaxis()->SetRangeUser(0,25);


    c1->cd();
    name = "UEvsBBCsumE" + ptBinName[i];
    title = "Charged Underlying Event vs. BBC ADC East Sum - p_{T}^{lead}: " + ptBinString[i];
    hPt_UE_BBCsumE->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hUE_BBCsumE[i] = (TH2D*)hPt_UE_BBCsumE->Project3D( "yz" );       // PROJECT
    hUE_BBCsumE[i]->GetYaxis()->SetRangeUser(0,10);
    scale = hUE_BBCsumE[i]->Integral("width");
    hUE_BBCsumE[i]->Scale( 1./scale );                     // NORMALIZE
    hUE_BBCsumE[i]->SetNameTitle(name,title);
    //hUE_BBCsumE[i]->Write();
    hUE_BBCsumE[i]->Draw("COLZ");
    hUE_BBCsumE[i]->SetLineColor(kRed);
    name = "pUE_BBCsumE" + ptBinName[i];
    pUE_BBCsumE[i] = hUE_BBCsumE[i]->ProfileX("",1,-1,"S");
    pUE_BBCsumE[i]->SetStats(0);
    pUE_BBCsumE[i]->SetName(name);
    pUE_BBCsumE[i]->Draw("SAME");
    name = "plots/" + dir + "UEvsBBCsumE" + ptBinName[i] + ".pdf";
    c1->SaveAs( name,"PDF");
    hUE_BBCsumE[i]->Scale( scale );                     // UN-NORMALIZE
    hUE_BBCsumE[i]->GetYaxis()->SetRangeUser(0,25);
  }

  c0->Clear();  c1->Clear();

  TLegend *leg4 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg4->SetBorderSize(1);   leg4->SetLineColor(1);   leg4->SetLineStyle(1);   leg4->SetLineWidth(1);   leg4->SetFillColor(0);   leg4->SetFillStyle(1001);
  leg4->SetNColumns(2);  leg4->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");  leg4->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  
  hscale1->SetStats(0);  c1->cd();    hscale1->GetYaxis()->SetRangeUser(0.0,4.0);     hscale1->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    pUE_BBCE[i]->SetLineColor( color[i] );    pUE_BBCE[i]->SetMarkerStyle( marker[i] );    pUE_BBCE[i]->SetMarkerColor( color[i] );    
    pUE_BBCE[i]->GetYaxis()->SetRangeUser(0.0,4.0);         pUE_BBCE[i]->Draw("SAME");
    title = "p_{T}^{lead}: " + ptBinString[i];
    avg = "";    avg += hUE_BBCE[i]->GetMean(2);    avg = avg(0,6);
    leg4->AddEntry( pUE_BBCE[i], title, "lpf" );    leg4->AddEntry((TObject*)0,avg, "");
  }
  leg4->Draw();        c1->Modified();        c1->cd();        c1->SetSelected(c1);
  path = "plots/" + dir + "UE_BBCE_profile.pdf";
  c1->SaveAs(  path ,"PDF");          c1->Clear();


  TLegend *leg5 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg5->SetBorderSize(1);   leg5->SetLineColor(1);   leg5->SetLineStyle(1);   leg5->SetLineWidth(1);   leg5->SetFillColor(0);   leg5->SetFillStyle(1001);
  leg5->SetNColumns(2);  leg5->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");  leg5->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  
  hscale2->SetStats(0);  c1->cd();    hscale2->GetYaxis()->SetRangeUser(0.0,5.0);     hscale2->Draw();
  for ( int i=0; i<nPtBins; ++i ) {
    pUE_BBCsumE[i]->SetLineColor( color[i] );    pUE_BBCsumE[i]->SetMarkerStyle( marker[i] );    pUE_BBCsumE[i]->SetMarkerColor( color[i] );    
    pUE_BBCsumE[i]->GetYaxis()->SetRangeUser(0.0,5.0);         pUE_BBCsumE[i]->Draw("SAME");
    title = "p_{T}^{lead}: " + ptBinString[i];
    avg = "";    avg += hUE_BBCsumE[i]->GetMean(2);    avg = avg(0,6);
    leg5->AddEntry( pUE_BBCsumE[i], title, "lpf" );    leg5->AddEntry((TObject*)0,avg, "");
  }
  leg5->Draw();        c1->Modified();        c1->cd();        c1->SetSelected(c1);
  path = "plots/" + dir + "UE_BBCsumE_profile.pdf";
  c1->SaveAs( path ,"PDF");          c1->Clear();


  TH2D *hscale_1 = new TH2D( "hscale_1", "Charged Underlying Event vs. BBC East Rate", 140,0,7000000, 50,0,25 );
  TH2D *hscale_2 = new TH2D( "hscale_2", "Charged Underlying Event vs. BBC ADC East Sum", 150,0,100000, 50,0,25 );

  hPt_UE_BBCsumE->GetYaxis()->SetRangeUser( 0.0,25 );
  hPt_UE_BBCE->GetYaxis()->SetRangeUser( 0.0,25 );

  TCanvas * c5 = new TCanvas( "c5" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  TCanvas * c6 = new TCanvas( "c6" , "" ,0 ,23 ,1280 ,700 );              // CANVAS

  TH2D *hscale3 = new TH2D( "hscale3", "Charged Underlying Event vs. BBC East Rate", 50,0,25, 100,0.000001,1 );
  TH2D *hscale4 = new TH2D( "hscale4", "Charged Underlying Event vs. BBC ADC East Sum", 50,0,25, 100,0.000001,1 );
  hscale3->SetStats(0);     hscale4->SetStats(0);

  hPt_UE_BBCsumE->GetYaxis()->SetRangeUser( 0.0,10 );
  hPt_UE_BBCE->GetYaxis()->SetRangeUser( 0.0,10 );
  
  TH1D *hUE_BBCE_py[nPtBins];        TH1D *hUE_BBCsumE_py[nPtBins];
    
  TLegend *leg0 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg0->SetBorderSize(1);   leg0->SetLineColor(1);   leg0->SetLineStyle(1);   leg0->SetLineWidth(1);   leg0->SetFillColor(0);   leg0->SetFillStyle(1001);
  leg0->SetNColumns(2);  leg0->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");  leg0->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  
  c5->SetLogy();  c5->cd();  c5->Clear();  hscale3->GetYaxis()->SetRangeUser(0,1);  // hscale3->GetXaxis()->SetRangeUser(0,5);
 
  TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  leg1->SetNColumns(2);  leg1->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");  leg1->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");
  
  c6->SetLogy();  c6->cd();  c6->Clear();  hscale4->GetYaxis()->SetRangeUser(0,1);  // hscale4->GetXaxis()->SetRangeUser(0.0,5.0);
  
  hLeadPtVsRho->GetXaxis()->SetRange(1,100);

  
  TCanvas * c3 = new TCanvas( "c3" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c3->SetLogy();
  hscale0->SetStats(0);
  hscale0->Draw();

  TLegend *leg2 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg2->SetBorderSize(1);   leg2->SetLineColor(1);   leg2->SetLineStyle(1);   leg2->SetLineWidth(1);   leg2->SetFillColor(0);   leg2->SetFillStyle(1001);
  leg2->SetNColumns(3);
  leg2->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg2->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg2->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");

  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadPtVsRho" + ptBinName[i];      title = ptBinString[i];
    hLeadPtVsRho->GetYaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hRho[i] = hLeadPtVsRho->ProjectionX( name );       // PROJECT
    hRho[i]->SetStats(0);
    hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
    hRho[i]->SetLineColor( color[i] );    hRho[i]->SetMarkerStyle( marker[i] );    hRho[i]->SetMarkerColor( color[i] );
    hRho[i]->Draw("SAME");                                                    // DRAW
    Ndj = ""; avg = "";    Ndj += hRho[i]->GetEntries();
    avg += hRho[i]->GetMean(1);                                           // 1 denotes x-axis
    avg = avg(0,6);
    leg2->AddEntry( name, title, lpf );                            // ADD TO LEGEND
    leg2->AddEntry((TObject*)0,Ndj, "");    leg2->AddEntry((TObject*)0,avg, "");    lpf += "lpf";
  }

  leg2->Draw();  c3->Modified();  c3->cd();  c3->SetSelected(c3);
  path = "plots/" + dir + "RhoByLeadPt.pdf";
  c3->SaveAs( path , "PDF");



  
  gStyle->SetOptStat(1);
  TCanvas * c4 = new TCanvas( "c4" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  
  hLeadEtaPhi->Scale( 1./hLeadEtaPhi->Integral("width") );
  hLeadEtaPhi->Draw("COLZ");
  path = "plots/" + dir + "LeadEtaPhi.pdf";
  c4->SaveAs( path ,"PDF");

  hSubEtaPhi->Scale( 1./hSubEtaPhi->Integral("width") );
  hSubEtaPhi->Draw("COLZ");
  path = "plots/" + dir + "SubEtaPhi.pdf";
  c4->SaveAs( path ,"PDF");

  c4->SetLogz();
  hTowersVsRho->Scale( 1./hTowersVsRho->Integral("width") );
  // hTowersVsRho->GetXaxis()->SetRangeUser( 0,7 );
  // hTowersVsRho->GetYaxis()->SetRangeUser( 0,200 );
  hTowersVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hTowersVsRho->Draw("COLZ");
  hTowersVsRho->SetLineColor(kRed);
  hTowersVsRho->ProfileX("",1,-1,"S")->Draw("SAME");
  path = "plots/" + dir + "TowersVsRho.pdf";
  c4->SaveAs( path ,"PDF");
  c4->Clear();
  
  hPrimaryVsGlobal->Scale( 1./hPrimaryVsGlobal->Integral("width") );
  hPrimaryVsGlobal->Draw("COLZ");
  path = "plots/" + dir + "PrimaryVsGlobal.pdf";
  c4->SaveAs( path ,"PDF");

  hPrimaryVsBBCsumE->Scale( 1./hPrimaryVsBBCsumE->Integral("width") );
  hPrimaryVsBBCsumE->Draw("COLZ");
  hPrimaryVsBBCsumE->SetLineColor(kRed);
  hPrimaryVsBBCsumE->ProfileX("",1,-1,"S")->Draw("SAME");
  path = "plots/" + dir + "PrimaryVsBBCsumE.pdf";
  c4->SaveAs( path ,"PDF");

  hTowersVsBBCsumE->GetXaxis()->SetRangeUser( 0,80000 );  
  hTowersVsBBCsumE->GetYaxis()->SetRangeUser( 0,200 );  
  hTowersVsBBCsumE->Scale( 1./hTowersVsBBCsumE->Integral("width") );
  hTowersVsBBCsumE->Draw("COLZ");
  hTowersVsBBCsumE->SetLineColor(kRed);
  hTowersVsBBCsumE->ProfileX("",0,80000,"S")->Draw("SAME");
  path = "plots/" + dir + "hTowersVsBBCsumE.pdf";
  c4->SaveAs( path ,"PDF");

  
  inFile->Close();

}
