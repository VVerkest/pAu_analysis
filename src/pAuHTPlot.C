//  pAuPlot.C
//  Veronica Verkest     May 22, 2019

void pAuHTPlot() {
  
  const float pi = 3.141592;
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-2,2, 16,-pi,pi );
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicity (per event);# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Multiplicity (per run);Run no.;# of Towers", 60,16120000,16160000, 700,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 60,16120000,16160000, 200,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 700,0,700, 200,0,200);
  TH2D *hDCAvsPt = new TH2D("hDCAvsPt","DCA vs. p_{T};p_{T} GeV;DCA (cm)", 100,0.0,50.0, 70,0,3.5);
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hTowEt = new TH2D("hTowEt","Tower E_{T} by ID;Tower ID;E_{T} (GeV)", 4800,0,4800, 40,0,20.0);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta;Tower #phi", 100,0,20, 200,-1.0,1.0, 2*pi,-pi,pi );
  TH1D *hTowerMult = new TH1D("hTowerMult","Tower Multiplicity by ID;Tower ID", 4800,0,4800);

  TFile* inFile = new TFile( "out/HT/pAu_analysis_HT.root", "READ" );

  TTree *HTtree = (TTree*) inFile->Get("HTTree");
  TTree *HTtowers = (TTree*) inFile->Get("HTTowers");
  TTree *HTtracks = (TTree*) inFile->Get("HTTracks");
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, towID, Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz, towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz;

  HTtree->SetBranchAddress("EventID", &EventID);            HTtree->SetBranchAddress("RunID", &RunID);          HTtree->SetBranchAddress("Vx", &Vx);
  HTtree->SetBranchAddress("Vy", &Vy);                             HTtree->SetBranchAddress("Vz", &Vz);                     HTtree->SetBranchAddress("nTowers", &nTowers);
  HTtree->SetBranchAddress("nPrimary", &nPrimary);         HTtree->SetBranchAddress("nGlobal", &nGlobal);     HTtree->SetBranchAddress("nVertices", &nVertices);
  HTtree->SetBranchAddress("refMult", &refMult);              HTtree->SetBranchAddress("gRefMult", &gRefMult);
  HTtree->SetBranchAddress("BbcCoincidenceRate", &BbcCoincidenceRate);                    HTtree->SetBranchAddress("BbcEastRate", &BbcEastRate);
  HTtree->SetBranchAddress("BbcWestRate", &BbcWestRate);                                           HTtree->SetBranchAddress("vpdVz", &vpdVz);
  
  HTtowers->SetBranchAddress("EventID", &EventID);	  HTtowers->SetBranchAddress("RunID", &RunID);  	  HTtowers->SetBranchAddress("towEt", &towEt);
  HTtowers->SetBranchAddress("towEta", &towEta);	  HTtowers->SetBranchAddress("towPhi", &towPhi);	  HTtowers->SetBranchAddress("towID", &towID);

  HTtracks->SetBranchAddress("EventID", &EventID);	  HTtracks->SetBranchAddress("RunID", &RunID);	  HTtracks->SetBranchAddress("nHitsPoss", &nHitsPoss);
  HTtracks->SetBranchAddress("nHitsFit", &nHitsFit);	  HTtracks->SetBranchAddress("trPx", &trPx);   	  HTtracks->SetBranchAddress("trPy", &trPy);
  HTtracks->SetBranchAddress("trPz", &trPz);           	  HTtracks->SetBranchAddress("trPt", &trPt);    	  HTtracks->SetBranchAddress("trEta", &trEta);
  HTtracks->SetBranchAddress("trPhi", &trPhi);        	  HTtracks->SetBranchAddress("DCA", &DCA);
  
  int mbEntries = HTtree->GetEntries();     int towEntries = HTtowers->GetEntries();     int trkEntries = HTtracks->GetEntries();
  //double vx, vy, vz, et, eta, phi, px, py, pz, pt, dca, bbcCoinc, bbcEast, bbcWest, vpdvz;
  
  for ( int i=0; i<mbEntries; ++i ) {
    HTtree->GetEntry(i);
    hVertex->Fill( Vx, Vy, Vz );                                        //  FILL HISTOGRAMS
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
  }

  for ( int i=0; i<towEntries; ++i ) {
    HTtowers->GetEntry(i);
    hTowEt->Fill( towID, towEt );
    hTowEtEtaPhi->Fill( towEt, towEta, towPhi );
    hTowerMult->Fill( towID );
  }

  for ( int i=0; i<trkEntries; ++i ) {
    HTtracks->GetEntry(i);
    hDCAvsPt->Fill( trPt, DCA );
    hPrimaryTracks->Fill( trPt, trEta, trPhi );
  }
  
  inFile->Close();

  hPrimaryTracks->Write();
  hVertex->Write();
  hTowersPerEvent->Write();
  hTowersPerRun->Write();
  hPrimaryPerEvent->Write();
  hPrimaryPerRun->Write();
  hnPrimaryVSnTowers->Write();
  hDCAvsPt->Write();
  hPrimaryVsBBC->Write();
  hGlobalVsBBC->Write();
  hTowEt->Write();
  hTowEtEtaPhi->Write();
  hTowerMult->Write();
  
  pAuPlots->Close();
}
