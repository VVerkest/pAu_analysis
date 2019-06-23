//  particlePlot.C
//  Veronica Verkest     June 12, 2019

void particlePlot() {
  
  const float pi = 3.141592;
  const double twopi = 2*3.14159265358979;
  const int nEtaBins = 4;
  const double etaBinLo[nEtaBins] = { -1.0, -0.5, 0.0, 0.5 };
  const double etaBinHi[nEtaBins] = { -0.5, 0.0, 0.5, 1.0 };
  const TString etaBinName[nEtaBins] = { "_n10_n05", "_n05_00", "_00_05", "_05_10" };
  const TString etaBinString[nEtaBins] = { "-1.0<#eta<-0.5", "-0.5<#eta<0.0", "0.0<#eta<0.5", "0.5<#eta<1.0" };
  const TString etaColor[nEtaBins] = { "kMagenta+3", "kBlue+1", "kTeal-6", "kGreen+3" };
  const int etaMarker[nEtaBins] = { 20, 21, 22, 23 };
  const int etaMarker2[nEtaBins] = { 24, 25, 26, 32 };

  const int nPtBins = 5;
  const double ptBinLo[nPtBins] = { 10.0, 15.0, 20.0, 30.0, 40.0 };
  const double ptBinHi[nPtBins] = { 15.0, 20.0, 30.0, 40.0, 100.0 };
  const TString ptBinString[nPtBins] = { "10-15 GeV", "15-20 GeV",  "20-30 GeV", "30-40 GeV", ">40 GeV" };
  const TString ptBinName[nPtBins] = { "_10_15", "_15_20", "_20_30", "_30_40", "_40" };
  const int ptColor[nPtBins] = { 633, 613, 596, 414, 797 };
  const int ptMarker[nPtBins] = { 33, 34, 22, 21, 20 };

  TString avg, name, title;      TString lpf = "lpf";        double scale;
  double rho, rho_avg;
  vector<double> partPt, partEta, partPhi, partEt, deltaPhi;
  double leadPt, subPt, nTowers;
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile* inFile = new TFile( "out/HTjets/pAu_HT_dijets.root", "READ" );

  TTree *sp[nEtaBins];

  for ( int i=0; i<nEtaBins; ++i ) {
    name = "sp" + etaBinName[i];    sp[nEtaBins] = (TTree*) inFile->Get(name);
    sp[i]->SetBranchAddress("rho",&rho);               sp[i]->SetBranchAddress("leadPt",&leadPt);
    sp[i]->SetBranchAddress("subPt",&subPt);        sp[i]->SetBranchAddress("nTowers",&nTowers);
    sp[i]->SetLineColor( etaColor[i] );                     sp[i]->SetMarkerStyle( etaMarker[i] );               sp[i]->SetMarkerColor( etaColor[i] );
  }

  TH1D *hUE = new TH1D("hUE", "Underlying Event by Eta");

  inFile->Close();

}


/*
  TTree *sp[nEtaBins];
  for ( int i=0; i<nEtaBins; ++i ) {
    name = "sp" + etaBinName[i];          title = "Selected Particles: " + etaBinString[i];          sp[i] = new TTree( name, title );
    sp[i]->Branch("partPt",&partPt);        sp[i]->Branch("partEta",&partEta);        sp[i]->Branch("partPhi",&partPhi);
    sp[i]->Branch("partEt",&partEt);        sp[i]->Branch("rho",&rho);        sp[i]->Branch("sigma",&sigma);
    sp[i]->Branch("deltaPhi",&deltaPhi);    sp[i]->Branch("leadPt",&leadPt);        sp[i]->Branch("subPt",&subPt);
    sp[i]->Branch("nTowers",&nTowers);        sp[i]->Branch("partChg",&partChg);
  }
