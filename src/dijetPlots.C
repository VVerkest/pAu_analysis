//  dijetPlots.C
//  Veronica Verkest     May 22, 2019

void dijetPlots() {
  
  const float pi = 3.141592;
  const double twopi = 2*3.14159265358979;
  TString Ndj, avg;
  TString lpf = "lpf";
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile* inFile = new TFile( "out/HTjets/pAu_HT_dijets.root", "READ" );

  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");
  TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");
  TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");
  TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");
  TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");
  TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");
  TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");
  TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");
  TH2D *hLeadEtaPhi =(TH2D*) inFile->Get("hLeadEtaPhi");
  TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");
  TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");

  hPt_UE_BBCsumE->GetYaxis()->SetRange( 2,50 );
  hPt_UE_BBCE->GetYaxis()->SetRange( 2,50 );
  hPt_UE_BBCsumE->GetZaxis()->SetRangeUser( 0.0, 80000.0 );

  TH2D *hscale0 = new TH2D( "hscale0", "Underlying Event by Lead Jet p_{T};#rho (GeV);", 60,0,30, 10,0.00000001, 1.0 );

  const int nPtBins = 5;
  TH1D * hRho[nPtBins];
  TH2D * hUE_BBCE[nPtBins];
  TH2D * hUE_BBCsumE[nPtBins];
  
  double ptBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0, 40.0 };
  double ptBinHi[nPtBins] = { 15.0, 20.0, 30.0, 40.0, 100.0 };
  TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", "30-40 GeV", ">40 GeV" };
  TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30_40", "_40" };
  int color[nPtBins] = { 633, 613, 596, 414, 797 };
  int marker[nPtBins] = { 33, 34, 22, 21, 20 };
  TString name, title;




  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c0->SetLogz();
  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c1->SetLogz();
  
  for ( int i=0; i<nPtBins; ++i ) {
    c0->cd();
    name = "plots/UEvsBBCE" + ptBinName[i] + ".pdf";
    title = "Underlying Event vs. BBC East Rate - p_{T}^{lead}: " + ptBinString[i];
    hPt_UE_BBCE->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hUE_BBCE[i] = (TH2D*)hPt_UE_BBCE->Project3D( "yz" );       // PROJECT
    hUE_BBCE[i]->Scale( 1./hUE_BBCE[i]->Integral("width") );                     // NORMALIZE
    hUE_BBCE[i]->SetTitle(title);
    hUE_BBCE[i]->Draw("COLZ");
    c0->SaveAs( name,"PDF");

    c1->cd();
    name = "plots/UEvsBBCsumE" + ptBinName[i] + ".pdf";
    title = "Underlying Event vs. BBC ADC East Sum - p_{T}^{lead}: " + ptBinString[i];
    hPt_UE_BBCsumE->GetXaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hUE_BBCsumE[i] = (TH2D*)hPt_UE_BBCsumE->Project3D( "yz" );       // PROJECT
    hUE_BBCsumE[i]->Scale( 1./hUE_BBCsumE[i]->Integral("width") );                     // NORMALIZE
    hUE_BBCsumE[i]->SetTitle(title);
    hUE_BBCsumE[i]->Draw("COLZ");
    c1->SaveAs( name,"PDF");
  }
  
  TCanvas * c2 = new TCanvas( "c2" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  hLeadPtVsRho->GetXaxis()->SetRange(2,80);
  double LeadVrhoScale = hLeadPtVsRho->Integral("width");
  hLeadPtVsRho->Scale(1./LeadVrhoScale);
  c2->SetLogz();
  hLeadPtVsRho->Draw("COLZ");
  c2->SaveAs("plots/LeadPtVsRho.pdf","PDF");
  hLeadPtVsRho->Scale(LeadVrhoScale);
  
  TCanvas * c3 = new TCanvas( "c3" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  c3->SetLogy();
  hscale0->SetStats(0);
  hscale0->Draw();

  TLegend *leg1 = new TLegend(0.65, 0.65, 0.9, 0.9,NULL,"brNDC");    // LEGEND
  leg1->SetBorderSize(1);   leg1->SetLineColor(1);   leg1->SetLineStyle(1);   leg1->SetLineWidth(1);   leg1->SetFillColor(0);   leg1->SetFillStyle(1001);
  TLegendEntry *entry;
  leg1->SetNColumns(3);
  leg1->AddEntry((TObject*)0,"#bf{p_{T}^{Lead} (GeV)}", "");
  leg1->AddEntry((TObject*)0,"#bf{# of Dijets}", "");
  leg1->AddEntry((TObject*)0,"#bf{<#rho> (GeV)}", "");

  for ( int i=0; i<nPtBins; ++i ) {
    name = "LeadPtVsRho" + ptBinName[i];      title = ptBinString[i];
    hLeadPtVsRho->GetYaxis()->SetRangeUser(ptBinLo[i], ptBinHi[i]);
    hRho[i] = hLeadPtVsRho->ProjectionX( name );       // PROJECT
    hRho[i]->SetStats(0);
    hRho[i]->Scale( 1./hRho[i]->Integral() );                     // NORMALIZE
    hRho[i]->SetLineColor( color[i] );    hRho[i]->SetMarkerStyle( marker[i] );    hRho[i]->SetMarkerColor( color[i] );
    hRho[i]->Draw("SAME");                                                    // DRAW
    /*entry=*/
    Ndj = ""; avg = "";
    Ndj += hRho[i]->GetEntries();
    avg += hRho[i]->GetMean(1);                                           // 1 denotes x-axis
    leg1->AddEntry( name, title, lpf );                            // ADD TO LEGEND
    leg1->AddEntry((TObject*)0,Ndj, "");
    leg1->AddEntry((TObject*)0,avg, "");
    lpf += "lpf";
  }

  leg1->Draw();
  c3->Modified();
  c3->cd();
  c3->SetSelected(c3);
  c3->SaveAs("plots/RhoByLeadPt.pdf","PDF");

  
  gStyle->SetOptStat(1);
  TCanvas * c4 = new TCanvas( "c4" , "" ,0 ,23 ,1280 ,700 );              // CANVAS
  
  hLeadEtaPhi->Scale( 1./hLeadEtaPhi->Integral("width") );
  hLeadEtaPhi->Draw("COLZ");
  c4->SaveAs("plots/LeadEtaPhi.pdf","PDF");

  hSubEtaPhi->Scale( 1./hSubEtaPhi->Integral("width") );
  hSubEtaPhi->Draw("COLZ");
  c4->SaveAs("plots/SubEtaPhi.pdf","PDF");

  c4->SetLogz();
  hTowersVsRho->Scale( 1./hTowersVsRho->Integral("width") );
  hTowersVsRho->GetXaxis()->SetRangeUser( 0,7 );
  hTowersVsRho->GetYaxis()->SetRangeUser( 0,200 );
  hTowersVsRho->Draw("COLZ");
  c4->SaveAs("plots/TowersVsRho.pdf","PDF");
  
  hPrimaryVsGlobal->Scale( 1./hPrimaryVsGlobal->Integral("width") );
  hPrimaryVsGlobal->Draw("COLZ");
  c4->SaveAs("plots/PrimaryVsGlobal.pdf","PDF");

  hPrimaryVsBBCsumE->Scale(1./hPrimaryVsBBCsumE->Integral("width"));
  c4->SaveAs("plots/PrimaryVsBBCsumE.pdf","PDF");
  c4->SetLogy();
  c4->SaveAs("plots/PrimaryVsBBCsumE_logy.pdf","PDF");


  
  inFile->Close();

  // TFile *djFile = new TFile( "plots/pAu_plots.root" ,"RECREATE");
  
  // djFile->Write();
  // djFile->Close();
}
