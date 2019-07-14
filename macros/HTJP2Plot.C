//  HTJP2Plot.C
//  Veronica Verkest		July 14, 2019

void HTJP2Plot() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  // const int acn = 3;
  // TString BackgroundChargeBias[acn] = { "allBG", "chgBG", "neuBG" };
  // TString JetChargeBias[acn] = { "allJets", "chgJets", "neuJets" };
  TString fileName;			TString name, title;			double scale;

  TString BackgroundChargeBias = "allBG";     //  [  "allBG",  "chgBG",  "neuBG"  ]
  TString JetChargeBias = "allJets";                  //  [ "allJets", "chgJets", "neuJets" ]

  fileName = "out/HTJP2/pAu_2015_200_" + BackgroundChargeBias + "_" + JetChargeBias + ".root";
  TFile* inFile = new TFile( fileName, "READ" );

  TH1D *hTowersPerEvent = (TH1D*) inFile->Get("hTowersPerEvent");			hTowersPerEvent->Scale(1.0/hTowersPerEvent->Integral("WIDTH"));
  TH1D *hPrimaryPerEvent = (TH1D*) inFile->Get("hPrimaryPerEvent");			hPrimaryPerEvent->Scale(1.0/hPrimaryPerEvent->Integral("WIDTH"));
  TH2D *hChgVsNeuBG = (TH2D*) inFile->Get("hChgVsNeuBG");				hChgVsNeuBG->Scale(1.0/hChgVsNeuBG->Integral("WIDTH"));
  TH2D *hTowersPerRun = (TH2D*) inFile->Get("hTowersPerRun");	       		hTowersPerRun->Scale(1.0/hTowersPerRun->Integral("WIDTH"));
  TH2D *hPrimaryPerRun = (TH2D*) inFile->Get("hPrimaryPerRun");		       	hPrimaryPerRun->Scale(1.0/hPrimaryPerRun->Integral("WIDTH"));
  TH2D *hnPrimaryVSnTowers = (TH2D*) inFile->Get("hnPrimaryVSnTowers");	hnPrimaryVSnTowers->Scale(1.0/hnPrimaryVSnTowers->Integral("WIDTH"));
  TH2D *hPrimaryVsBBC = (TH2D*) inFile->Get("hPrimaryVsBBC");				hPrimaryVsBBC->Scale(1.0/hPrimaryVsBBC->Integral("WIDTH"));
  TH2D *hPrimaryVsGlobal = (TH2D*) inFile->Get("hPrimaryVsGlobal");		       	hPrimaryVsGlobal->Scale(1.0/hPrimaryVsGlobal->Integral("WIDTH"));
  // TH2D *hGlobalVsBBC = (TH2D*) inFile->Get("hGlobalVsBBC");				hGlobalVsBBC->Scale(1.0/hGlobalVsBBC->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCE = (TH2D*) inFile->Get("hPrimaryVsBBCE");			hPrimaryVsBBCE->Scale(1.0/hPrimaryVsBBCE->Integral("WIDTH"));
  // TH2D *hGlobalVsBBCE = (TH2D*) inFile->Get("hGlobalVsBBCE");				hGlobalVsBBCE->Scale(1.0/hGlobalVsBBCE->Integral("WIDTH"));
  // TH2D *hPrimaryVsBBCsumE = (TH2D*) inFile->Get("hPrimaryVsBBCsumE");		hPrimaryVsBBCsumE->Scale(1.0/hPrimaryVsBBCsumE->Integral("WIDTH"));
  // TH2D *hTowersVsBBCsumE = (TH2D*) inFile->Get("hTowersVsBBCsumE");		hTowersVsBBCsumE->Scale(1.0/hTowersVsBBCsumE->Integral("WIDTH"));
  TH2D *hPrimaryVsRho = (TH2D*) inFile->Get("hPrimaryVsRho");				hPrimaryVsRho->Scale(1.0/hPrimaryVsRho->Integral("WIDTH"));
  TH2D *hGlobalVsRho = (TH2D*) inFile->Get("hGlobalVsRho");				hGlobalVsRho->Scale(1.0/hGlobalVsRho->Integral("WIDTH"));
  TH2D *hSubEtaPhi = (TH2D*) inFile->Get("hSubEtaPhi");					hSubEtaPhi->Scale(1.0/hSubEtaPhi->Integral("WIDTH"));
  TH2D *hTowersVsRho = (TH2D*) inFile->Get("hTowersVsRho");				hTowersVsRho->Scale(1.0/hTowersVsRho->Integral("WIDTH"));
  TH2D *hLeadPtVsRho = (TH2D*) inFile->Get("hLeadPtVsRho");
  TH3D *hPt_UE_BBCE = (TH3D*) inFile->Get("hPt_UE_BBCE");
  TH3D *hPt_UE_BBCsumE = (TH3D*) inFile->Get("hPt_UE_BBCsumE");
  TH3D *hBG = (TH3D*) inFile->Get("hBG");
  TH3D *hNEUTRAL = (TH3D*) inFile->Get("hBG");
  TH3D *hPt_UE_RefMult = (TH3D*) inFile->Get("hPt_UE_RefMult");
  TH3D *hCHARGED = (TH3D*) inFile->Get("hPt_UE_RefMult");
  TH3D *hAllJetsPtRhoEta = (TH3D*) inFile->Get("hAllJetsPtRhoEta");
  TH3D *hPartPtEtaPhi = (TH3D*) inFile->Get("hPartPtEtaPhi");
  TH3D *hAllPtEtaPhi = (TH3D*) inFile->Get("hAllPtEtaPhi");
  TH3D *hPartPtDEtaDPhi = (TH3D*) inFile->Get("hPartPtDEtaDPhi");
  TH3D *hVertex = (TH3D*) inFile->Get("hVertex");
  TH3D *hLeadPtEtaPhi = (TH3D*) inFile->Get("hLeadPtEtaPhi");
  TH3D *hAllJetsPtEtaPhi = (TH3D*) inFile->Get("hAllJetsPtEtaPhi");




  TCanvas * c0 = new TCanvas( "c0" , "" ,700 ,500 );              // CANVAS
  
  c0->SetLogy();
  
  hTowersPerEvent->Draw();				name = "plots/HTJP2/TowersPerEvent_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  hPrimaryPerEvent->Draw();				name = "plots/HTJP2/PrimaryPerEvent_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  TH1D *hAllJetsPt = (TH1D*) hAllJetsPtEtaPhi->ProjectionX();		hAllJetsPt->Scale(1.0/hAllJetsPt->Integral());			hAllJetsPt->SetTitle("Inclusive Jet p_{T};p_{T} (GeV)");
  hAllJetsPt->Draw();					name = "plots/HTJP2/AllJetsPt_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";					c0->SaveAs( name, "PDF" );

  c0->SetLogy(0);		c0->SetLogz();

  hChgVsNeuBG->GetZaxis()->SetRangeUser(0.0000001,1);
  hChgVsNeuBG->Draw("COLZ");			name = "plots/HTJP2/ChgVsNeuBG_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  hnPrimaryVSnTowers->GetZaxis()->SetRangeUser(0.0000001,1);
  hnPrimaryVSnTowers->Draw("COLZ");		name = "plots/HTJP2/nPrimaryVSnTowers_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  hPrimaryVsBBC->GetZaxis()->SetRangeUser(0.0000001,1);
  hPrimaryVsBBC->Draw("COLZ");			name = "plots/HTJP2/PrimaryVsBBC_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  hPrimaryVsGlobal->GetZaxis()->SetRangeUser(0.0000001,1);
  hPrimaryVsGlobal->Draw("COLZ");		name = "plots/HTJP2/PrimaryVsGlobal_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";			c0->SaveAs( name, "PDF" );

  hPrimaryVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hPrimaryVsRho->Draw("COLZ");			name = "plots/HTJP2/PrimaryVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  hGlobalVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hGlobalVsRho->Draw("COLZ");			name = "plots/HTJP2/GlobalVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  hSubEtaPhi->GetZaxis()->SetRangeUser(0.0000001,1);
  hSubEtaPhi->Draw("COLZ");				name = "plots/HTJP2/SubEtaPhi_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";					c0->SaveAs( name, "PDF" );

  hTowersVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hTowersVsRho->Draw("COLZ");			name = "plots/HTJP2/TowersVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );

  scale = hLeadPtVsRho->Integral("WIDTH");			hLeadPtVsRho->Scale(1.0/scale);			hLeadPtVsRho->GetZaxis()->SetRangeUser(0.0000001,1);
  hLeadPtVsRho->Draw("COLZ");			name = "plots/HTJP2/LeadPtVsRho_" + BackgroundChargeBias + "_" + JetChargeBias + ".pdf";				c0->SaveAs( name, "PDF" );
  hLeadPtVsRho->Scale(scale);


  inFile->Close();

}
