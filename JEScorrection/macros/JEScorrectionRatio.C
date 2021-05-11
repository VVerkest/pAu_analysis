// Veronica Verkest
// March 25, 2021

void JEScorrectionRatio(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile* inFile = new TFile("out/test0.root", "READ");

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };
  const int nEAbins = 2;
  const TString EAbinName[nEAbins] = { "Lo", "Hi" };
  const TString lohi[nEAbins] = { "lo", "hi" };
  const TString EAbinString[nEAbins] = { "Low EA", "High EA" }; 

  TString name;
  
  double x1[4] = {10,15,20,30};
  double y1[2] = {0.5,1.8};
  TH2D *hScaleMult = new TH2D("hScaleMult","",3,x1,1,y1);

  double y2[2] = {0.55,0.85};
  TH2D *hScalePt = new TH2D("hScalePt","",3,x1,1,y2);

  double y3[2] = {0.8, 1.3};
  TH2D *hScaleRatio = new TH2D("hScaleRatio","",3,x1,1,y2);

  TH1D *hMeanPt[nEAbins], *hMeanMult[nEAbins], *hMeanPt_noJES[nEAbins], *hMeanMult_noJES[nEAbins], *rMeanPt[nEAbins], *rMeanMult[nEAbins];
  
  for (int a=0; a<nEAbins; ++a) {
    name = "hMeanPt_" + lohi[a]; 
    hMeanPt[a] = (TH1D*)inFile->Get(name);
    name = "hMeanMult_" + lohi[a]; 
    hMeanMult[a] = (TH1D*)inFile->Get(name);
    hMeanPt[a]->SetLineColor(kRed);
    hMeanPt[a]->SetMarkerColor(kRed);
    hMeanMult[a]->SetLineColor(kRed);
    hMeanMult[a]->SetMarkerColor(kRed);
    
    name = "hMeanPt_noJES_" + lohi[a]; 
    hMeanPt_noJES[a] = (TH1D*)inFile->Get(name);
    name = "hMeanMult_noJES_" + lohi[a]; 
    hMeanMult_noJES[a] = (TH1D*)inFile->Get(name);
    hMeanPt_noJES[a]->SetLineColor(kBlue);
    hMeanPt_noJES[a]->SetMarkerColor(kBlue);
    hMeanMult_noJES[a]->SetLineColor(kBlue);
    hMeanMult_noJES[a]->SetMarkerColor(kBlue);

    name = "rMeanPt_" + lohi[a]; 
    rMeanPt[a] = (TH1D*)hMeanPt_noJES[a]->Clone();
    rMeanPt[a]->Divide(hMeanPt[a]);
    // rMeanPt[a]->SetMakerStyle(kCircle);
    // rMeanPt[a]->SetMakerColor(kBlack);
    // rMeanPt[a]->SetLineColor(kBlack);
    
    name = "rMeanMult_" + lohi[a]; 
    rMeanMult[a] = (TH1D*)hMeanMult_noJES[a]->Divide(hMeanMult[a]);
    // rMeanMult[a]->SetMakerStyle(kCircle);
    // rMeanMult[a]->SetMakerColor(kBlack);
    // rMeanMult[a]->SetLineColor(kBlack);
  }



  TCanvas *c1 = new TCanvas("c1","",600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();
  


  for (int a=0; a<nEAbins; ++a) {
    
    pad1->cd();
    hScalePt->Draw();
    hMeanPt_noJES[a]->Draw("psame");
    hMeanPt[a]->Draw("psame");
    c1->cd();

    pad2->cd();
    hScaleRatio->GetYaxis()->SetRangeUser(0.8,1.2);
    hScaleRatio->Draw();
    rMeanPt[a]->SetStats(0);
    rMeanPt[a]->GetYaxis()->SetRangeUser(0.8,1.2);
    rMeanPt[a]->SetAxisRange(0.8,1.2,"Y");
    rMeanPt[a]->Draw("epsame");
    c1->cd();
    
  }
  
}
