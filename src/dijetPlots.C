//  dijetPlots.C
//  Veronica Verkest     May 22, 2019

void dijetPlots() {
  
  const float pi = 3.141592;
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile* inFile = new TFile( "out/MB/pAu_HT_dijets.root", "READ" );


  TH3D *hPrimaryTracks = (TH3D*) inFile->Get("hPrimaryTracks");
  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  TH3D *hTowEtEtaPhi = (TH3D*) inFile->Get("hTowEtEtaPhi");
  TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");
  TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");
  TH2D *hDCAvsPt = (TH2D*) inFile->Get("hDCAvsPt");
  TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");
  TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");
  TH2D *hTowEt = (TH2D*) inFile->Get("hTowEt");
  TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");
  TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");
  TH2D *hLeadEtaPhi =(TH2D*) inFile->Get("hLeadEtaPhi");
  TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");
  TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH1D *hTowersPerEvent = (TH12D*) inFile->Get("hTowersPerEvent");
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");
  TH1D *hTowerMult =(TH1D*) inFile->Get("hTowerMult");

  // hPrimaryTracks
  // hVertex
  // hTowersPerEvent
  // hTowersPerRun
  // hPrimaryPerEvent
  // hPrimaryPerRun
  // hnPrimaryVSnTowers
  // hPrimaryVsBBC
  // hGlobalVsBBC
  // hTowEt
  // hTowEtEtaPhi
  // hDCAvsPt
  // hTowerMult
  // hPrimaryVsBBCE
  // hGlobalVsBBCE
  // hLeadEtaPhi
  // hSubEtaPhi
  // hPt_UE_BBCE
  // hTowersVsRho
  // hLeadPtVsRho

  inFile->Close();

  TFile *djFile = new TFile( "plots/pAu_plots.root" ,"RECREATE");

  // TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  // c0->SetLogz();
  // hPrimaryPerEvent->Draw("COLZ");
  // hVertex->Draw();
  
  djFile->Close();
}
