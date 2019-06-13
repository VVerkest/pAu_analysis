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

  TString Ndj, avg, name, title;      TString lpf = "lpf";        double scale;

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile* inFile = new TFile( "out/HTjets/pAu_HT_dijets.root", "READ" );

  TTree *sp[nEtaBins];

  for ( int i=0; i<nEtaBins; ++i ) {
    name = "sp" + etaBinName[i];    sp[nEtaBins] = (TTree*) inFile->Get(name);
    sp[i]->SetBranchAddress("partPt",&partPt);        sp[i]->SetBranchAddress("partEta",&partEta);        sp[i]->SetBranchAddress("partPhi",&partPhi);
    sp[i]->SetBranchAddress("partEt",&partEt);        sp[i]->SetBranchAddress("rho",&rho);        sp[i]->SetBranchAddress("sigma",&sigma);
    sp[i]->SetBranchAddress("BMErho",&BMErho);        sp[i]->SetBranchAddress("BMEsigma",&BMEsigma);        sp[i]->SetBranchAddress("deltaPhi",&deltaPhi);
    sp[i]->SetBranchAddress("leadPt",&leadPt);        sp[i]->SetBranchAddress("subPt",&subPt);        sp[i]->SetBranchAddress("nTowers",&nTowers);
  }

   
  // inFile->Close();


}
