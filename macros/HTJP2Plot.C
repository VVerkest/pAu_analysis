//  HTJP2Plot.C
//  Veronica Verkest		July 14, 2019

void HTJP2Plot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TString BackgroundChargeBias = "allBG";     // options: chgBG, neuBG, allBG
  TString JetChargeBias = "allJets";     // options: chgJets, neuJets, allJets
  TString fileName;			TString name, title;

  fileName = "out/HTJP2/pAu_2015_200_" + BackgroundChargeBias + "_" + JetChargeBias + ".root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");		TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");
  TH2D *hChgVsNeuBG = (TH2D*) inFile->Get("hChgVsNeuBG");			TH2D *hTowersPerRun = (TH2D*) inFile->Get("hChgVsNeuBG");
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");		TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");
  TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");			TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");
  TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");			TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");
  TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");			TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");
  TH2D *hTowersVsBBCsumE = (TH2D*) inFile->Get("hTowersVsBBCsumE");	TH2D *hPrimaryVsRho = (TH2D*) inFile->Get("hPrimaryVsRho");
  TH2D *hGlobalVsRho = (TH2D*) inFile->Get("hGlobalVsRho");			TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");			TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");			TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  TH3D *hBG = (TH3D*) inFile->Get("hBG");							TH3D *hNEUTRAL = (TH3D*) inFile->Get("hBG");
  TH3D *hPt_UE_RefMult = (TH3D*) inFile->Get("hPt_UE_RefMult");		TH3D *hCHARGED = (TH3D*) inFile->Get("hPt_UE_RefMult");
  TH3D *hAllJetsPtRhoEta = (TH3D*) inFile->Get("hAllJetsPtRhoEta");		TH3D *hPartPtEtaPhi = (TH3D*) inFile->Get("hPartPtEtaPhi");
  TH3D *hAllPtEtaPhi = (TH3D*) inFile->Get("hAllPtEtaPhi");				TH3D *hPartPtDEtaDPhi = (TH3D*) inFile->Get("hPartPtDEtaDPhi");
  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");					TH3D *hLeadPtEtaPhi = (TH3D*) inFile->Get("hLeadPtEtaPhi");
  TH3D *hAllJetsPtEtaPhi = (TH3D*) inFile->Get("hAllJetsPtEtaPhi");


























}
