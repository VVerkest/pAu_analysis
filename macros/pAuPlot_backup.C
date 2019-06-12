//  pAuPlot.C
//  Veronica Verkest     May 22, 2019

void pAuPlot() {
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicity (per event);# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Multiplicity (per run);Run no.;# of Towers", 60,16120000,16160000, 700,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 60,16120000,16160000, 200,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 700,0,700, 200,0,200);
  //TH2D *hDCAvsPt = new TH2D("hDCAvsPt","DCA vs. p_{T};DCA;p_{T} GeV");
  const float pi = 3.141592;

  TFile* inFile = new TFile( "out/MB/pAu_analysis.root", "READ" );
  
  TTree *MBtree = (TTree*) inFile->Get("MBTree");
  TTree *MBtowers = (TTree*) inFile->Get("MBTowers");
  TTree *MBtracks = (TTree*) inFile->Get("MBTracks");

  int RunID, EventID, nTowers, nPrimary;
  int Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz;
  double towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA;

  MBtree->SetBranchAddress("EventID", &EventID);            MBtree->SetBranchAddress("RunID", &RunID);          MBtree->SetBranchAddress("Vx", &Vx);
  MBtree->SetBranchAddress("Vy", &Vy);                             MBtree->SetBranchAddress("Vz", &Vz);                     MBtree->SetBranchAddress("nTowers", &nTowers);
  MBtree->SetBranchAddress("nPrimary", &nPrimary);

  // MBtowers->SetBranchAddress("EventID", &EventID);	  MBtowers->SetBranchAddress("RunID", &RunID);  	  MBtowers->SetBranchAddress("towEt", &towEt);
  // MBtowers->SetBranchAddress("towEta", &towEta);          MBtowers->SetBranchAddress("towPhi", &towPhi);

  // MBtracks->SetBranchAddress("EventID", &EventID);	  MBtracks->SetBranchAddress("RunID", &RunID);	  MBtracks->SetBranchAddress("nHitsPoss", &nHitsPoss);
  // MBtracks->SetBranchAddress("nHitsFit", &nHitsFit);	  MBtracks->SetBranchAddress("trPx", &trPx);   	  MBtracks->SetBranchAddress("trPy", &trPy);
  // MBtracks->SetBranchAddress("trPz", &trPz);           	  MBtracks->SetBranchAddress("trPt", &trPt);    	  MBtracks->SetBranchAddress("trEta", &trEta);
  // MBtracks->SetBranchAddress("trPhi", &trPhi);        	  MBtracks->SetBranchAddress("DCA", &DCA);

  int mbEntries = MBtree->GetEntries();     int towEntries = MBtowers->GetEntries();     int trkEntries = MBtracks->GetEntries();

  for ( int i=0; i<mbEntries; ++i ) {
    MBtree->GetEntry(i);
    hVertex->Fill( Vx, Vy, Vz );
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
  }

  inFile->Close();
  
  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  c0->SetLogz();
  hPrimaryPerEvent->Draw("COLZ");
  // hVertex->Draw();

  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );
  c1->SetLogz();
  hPrimaryPerRun->Draw("COLZ");
  // hVertex->Project3D("XY")->Draw("colz");

  TCanvas * c2 = new TCanvas( "c2" , "" ,0 ,23 ,1280 ,700 );
  c2->SetLogz();
  hnPrimaryVSnTowers->Draw("COLZ");
  // hVertex->Project3D("YZ")->Draw("colz");

  // TCanvas * c3 = new TCanvas( "c3" , "" ,0 ,23 ,1280 ,700 );
  // c3->SetLogz();
  // hVertex->Project3D("XZ")->Draw("colz");

  TCanvas * c4 = new TCanvas( "c4" , "" ,0 ,23 ,1280 ,700 );
  c4->SetLogz();
  hTowersPerEvent->Draw("COLZ");
  
  TCanvas * c5 = new TCanvas( "c5" , "" ,0 ,23 ,1280 ,700 );
  c5->SetLogz();
  hTowersPerRun->Draw("COLZ");

}
