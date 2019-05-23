//  pAuPlot.C
//  Veronica Verkest     May 22, 2019

void pAuPlot() {
  
  const float pi = 3.141592;
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH2D *hDCAvsPt = new TH2D("hDCAvsPt","DCA vs. p_{T};p_{T} GeV;DCA (cm)??", 20,0,10, 60,-3.0,3.0);
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hTowEt = new TH2D("hTowEt","Tower E_{T} by ID;Tower ID;E_{T} (GeV)", 4800,0,4800, 60,-3.0,3.0);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta,Tower #phi", 1000,0,200, 200,-1.0,1.0, 2*pi,-pi,pi );
  //  TH2D * = new TH2D("","",);

  TFile* inFile = new TFile( "out/MB/pAu_analysis.root", "READ" );

  TTree *MBtree = (TTree*) inFile->Get("MBTree");
  TTree *MBtowers = (TTree*) inFile->Get("MBTowers");
  TTree *MBtracks = (TTree*) inFile->Get("MBTracks");
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, towID, Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz, towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz;

  MBtree->SetBranchAddress("EventID", &EventID);            MBtree->SetBranchAddress("RunID", &RunID);          MBtree->SetBranchAddress("Vx", &Vx);
  MBtree->SetBranchAddress("Vy", &Vy);                             MBtree->SetBranchAddress("Vz", &Vz);                     MBtree->SetBranchAddress("nTowers", &nTowers);
  MBtree->SetBranchAddress("nPrimary", &nPrimary);         MBtree->SetBranchAddress("nGlobal", &nGlobal);     MBtree->SetBranchAddress("nVertices", &nVertices);
  MBtree->SetBranchAddress("refMult", &refMult);              MBtree->SetBranchAddress("gRefMult", &gRefMult);
  MBtree->SetBranchAddress("BbcCoincidenceRate", &BbcCoincidenceRate);                    MBtree->SetBranchAddress("BbcEastRate", &BbcEastRate);
  MBtree->SetBranchAddress("BbcWestRate", &BbcWestRate);                                           MBtree->SetBranchAddress("vpdVz", &vpdVz);
  
  MBtowers->SetBranchAddress("EventID", &EventID);	  MBtowers->SetBranchAddress("RunID", &RunID);  	  MBtowers->SetBranchAddress("towEt", &towEt);
  MBtowers->SetBranchAddress("towEta", &towEta);	  MBtowers->SetBranchAddress("towPhi", &towPhi);	  MBtowers->SetBranchAddress("towID", &towID);

  MBtracks->SetBranchAddress("EventID", &EventID);	  MBtracks->SetBranchAddress("RunID", &RunID);	  MBtracks->SetBranchAddress("nHitsPoss", &nHitsPoss);
  MBtracks->SetBranchAddress("nHitsFit", &nHitsFit);	  MBtracks->SetBranchAddress("trPx", &trPx);   	  MBtracks->SetBranchAddress("trPy", &trPy);
  MBtracks->SetBranchAddress("trPz", &trPz);           	  MBtracks->SetBranchAddress("trPt", &trPt);    	  MBtracks->SetBranchAddress("trEta", &trEta);
  MBtracks->SetBranchAddress("trPhi", &trPhi);        	  MBtracks->SetBranchAddress("DCA", &DCA);
  
  int mbEntries = MBtree->GetEntries();     int towEntries = MBtowers->GetEntries();     int trkEntries = MBtracks->GetEntries();
  //double vx, vy, vz, et, eta, phi, px, py, pz, pt, dca, bbcCoinc, bbcEast, bbcWest, vpdvz;
  
  for ( int i=0; i<mbEntries; ++i ) {
    MBtree->GetEntry(i);
    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
  }

  for ( int i=0; i<towEntries; ++i ) {
    MBtowers->GetEntry(i);
    hTowEt->Fill( towID, towEt );
    hTowEtEtaPhi->Fill( towEt, towEta, towPhi );
  }

  for ( int i=0; i<trkEntries; ++i ) {
    MBtracks->GetEntry(i);
    hDCAvsPt->Fill( trPt, DCA );
  }
  
  inFile->Close();

  TFile *pAuPlots = new TFile( "plots/pAu_plots.root" ,"RECREATE");
  hPrimaryVsBBC->Write();
  hGlobalVsBBC->Write();
  hTowEt->Write();
  hTowEtEtaPhi->Write();
  hDCAvsPt->Write();
  
  // TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  // c0->SetLogz();
  // hPrimaryPerEvent->Draw("COLZ");
  // hVertex->Draw();

  // TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );
  // c1->SetLogz();
  // hPrimaryPerRun->Draw("COLZ");
  // hVertex->Project3D("XY")->Draw("colz");

  // TCanvas * c2 = new TCanvas( "c2" , "" ,0 ,23 ,1280 ,700 );
  // c2->SetLogz();
  // hnPrimaryVSnTowers->Draw("COLZ");
  // hVertex->Project3D("YZ")->Draw("colz");

  // TCanvas * c3 = new TCanvas( "c3" , "" ,0 ,23 ,1280 ,700 );
  // c3->SetLogz();
  // hVertex->Project3D("XZ")->Draw("colz");

  // TCanvas * c4 = new TCanvas( "c4" , "" ,0 ,23 ,1280 ,700 );
  // c4->SetLogz();
  // hTowersPerEvent->Draw("COLZ");
  
  // TCanvas * c5 = new TCanvas( "c5" , "" ,0 ,23 ,1280 ,700 );
  // c5->SetLogz();
  // hTowersPerRun->Draw("COLZ");
  
  pAuPlots->Close();
}
