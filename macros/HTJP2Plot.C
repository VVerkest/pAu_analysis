//  HTJP2Plot.C
//  Veronica Verkest		July 14, 2019

void HTJP2Plot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const int acn = 3;
  TString BackgroundChargeBias[acn] = { "allBG", "chgBG", "neuBG" };
  TString JetChargeBias[acn] = { "allJets", "chgJets", "neuJets" };
  TString fileName;			TString name, title;			TFile* inFile[acn][acn];
  
  TH1D *hTowersPerEvent[bb][jb];
  TH1D *hPrimaryPerEvent[bb][jb];
  TH2D *hChgVsNeuBG[bb][jb];
  TH2D *hTowersPerRun[bb][jb];
  TH2D *hPrimaryPerRun[bb][jb];
  TH2D *hnPrimaryVSnTowers[bb][jb];
  TH2D *hPrimaryVsBBC[bb][jb];
  TH2D *hPrimaryVsGlobal[bb][jb];
  TH2D *hGlobalVsBBC[bb][jb];
  TH2D *hPrimaryVsBBCE[bb][jb];
  TH2D *hGlobalVsBBCE[bb][jb];
  TH2D *hPrimaryVsBBCsumE[bb][jb];
  TH2D *hTowersVsBBCsumE[bb][jb];
  TH2D *hPrimaryVsRho[bb][jb];
  TH2D *hGlobalVsRho[bb][jb];
  TH2D *hSubEtaPhi[bb][jb];
  TH2D *hTowersVsRho[bb][jb];
  TH2D *hLeadPtVsRho[bb][jb];
  TH3D *hPt_UE_BBCE[bb][jb];
  TH3D *hPt_UE_BBCsumE[bb][jb];
  TH3D *hBG[bb][jb];
  TH3D *hNEUTRAL[bb][jb];
  TH3D *hPt_UE_RefMult[bb][jb];
  TH3D *hCHARGED[bb][jb];
  TH3D *hAllJetsPtRhoEta[bb][jb];
  TH3D *hPartPtEtaPhi[bb][jb];
  TH3D *hAllPtEtaPhi[bb][jb];
  TH3D *hPartPtDEtaDPhi[bb][jb];
  TH3D *hVertex[bb][jb];
  TH3D *hLeadPtEtaPhi[bb][jb];
  TH3D *hAllJetsPtEtaPhi[bb][jb];

  for ( int bb=0; bb<acn; ++bb ) {
    for ( int jb=0; jb<acn; ++jb ) {

      if ( bb!=jb && bb!=0 && jb!=0 ) { continue; }  // skip (chgBG,neuJet) and (neuBG,chgJet)
      
      fileName = "out/HTJP2/pAu_2015_200_" + BackgroundChargeBias + "_" + JetChargeBias + ".root";
      inFile[bb][jb] = new TFile( fileName, "READ" );

      hTowersPerEvent[bb][jb] = (TH1D*) inFile[bb][jb]->Get("hTowersPerEvent");
      hPrimaryPerEvent[bb][jb] = (TH1D*) inFile[bb][jb]->Get("hPrimaryPerEvent");
      hChgVsNeuBG[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hChgVsNeuBG");
      hTowersPerRun[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hTowersPerRun");
      hPrimaryPerRun[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryPerRun");
      hnPrimaryVSnTowers[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hnPrimaryVSnTowers");
      hPrimaryVsBBC[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryVsBBC");
      hPrimaryVsGlobal[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryVsGlobal");
      hGlobalVsBBC[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hGlobalVsBBC");
      hPrimaryVsBBCE[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryVsBBCE");
      hGlobalVsBBCE[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hGlobalVsBBCE");
      hPrimaryVsBBCsumE[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryVsBBCsumE");
      hTowersVsBBCsumE[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hTowersVsBBCsumE");
      hPrimaryVsRho[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hPrimaryVsRho");
      hGlobalVsRho[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hGlobalVsRho");
      hSubEtaPhi[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hSubEtaPhi");
      hTowersVsRho[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hTowersVsRho");
      hLeadPtVsRho[bb][jb] = (TH2D*) inFile[bb][jb]->Get("hLeadPtVsRho");
      hPt_UE_BBCE[bb][jb] = (TH3D*) inFile[bb][jb]->Get(hPt_UE_BBCE"");
      hPt_UE_BBCsumE[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hPt_UE_BBCsumE");
      hBG[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hBG");
      hNEUTRAL[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hBG");
      hPt_UE_RefMult[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hPt_UE_RefMult");
      hCHARGED[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hPt_UE_RefMult");
      hAllJetsPtRhoEta[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hAllJetsPtRhoEta");
      hPartPtEtaPhi[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hPartPtEtaPhi");
      hAllPtEtaPhi[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hAllPtEtaPhi");
      hPartPtDEtaDPhi[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hPartPtDEtaDPhi");
      hVertex[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hVertex");
      hLeadPtEtaPhi[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hLeadPtEtaPhi");
      hAllJetsPtEtaPhi[bb][jb] = (TH3D*) inFile[bb][jb]->Get("hAllJetsPtEtaPhi");
    }
  }



}
