//  pAuPlot.C
//  Veronica Verkest     May 22, 2019

void pAuPlot() {
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 120,-30,30 );
  TH3D *hTowersPerEvent = new TH2D("hTowersPerEvent","Tower Multiplicity (per event);Event no.;# of Towers", 15000,0,1500000, 700,0,700 );
  TH3D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Multiplicity (per run);Run no.;# of Towers", 15000,2000000,17000000, 700,0,700 );
  // 2119440-16156005      2000000,17000000
  const float pi = 3.141592;

  TFile* inFile = new TFile( "out/MB/pAu_2015_200_MB_125_130_0_MB.root", "READ" );
  
  TTree *MBtree = (TTree*) inFile->Get("MBTree");
  TTree *MBtowers = (TTree*) inFile->Get("MBTowers");
  TTree *MBtracks = (TTree*) inFile->Get("MBTracks");

  int RunID, EventID, nTowers, nPrimary, nTracks;
  int Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz;
  double towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA;

  MBtree->SetBranchAddress("EventID", &EventID);            MBtree->SetBranchAddress("RunID", &RunID);          MBtree->SetBranchAddress("Vx", &Vx);
  MBtree->SetBranchAddress("Vy", &Vy);                             MBtree->SetBranchAddress("Vz", &Vz);                     MBtree->SetBranchAddress("nTowers", &nTowers);
  MBtree->SetBranchAddress("nTracks", &nTracks);           MBtree->SetBranchAddress("nPrimary", &nPrimary);

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
    hTowersPerEvent->Fill( EventID, nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
  }

  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,700 );
  hVertex->Draw();

  TCanvas * c1 = new TCanvas( "c1" , "" ,0 ,23 ,1280 ,700 );
  hVertex->Project3D("XY")->Draw("colz");

  TCanvas * c2 = new TCanvas( "c2" , "" ,0 ,23 ,1280 ,700 );
  hVertex->Project3D("YZ")->Draw("colz");

  TCanvas * c3 = new TCanvas( "c3" , "" ,0 ,23 ,1280 ,700 );
  hVertex->Project3D("XZ")->Draw("colz");
}
