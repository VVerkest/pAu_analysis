//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019

//  (from Isaac's code:)
//for pAu, need to ask if the event has the trigger. The trigger IDs are:
//HT2*BBCMB : 500205, 500215
//JP2 : 500401, 500411
//BBCMB : 500008, 500018
//VPDMB :  500904


#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;


int main ( int argc, const char** argv ) {
  
  string inFile, outFile;
  int nEvents;
  if ( argc ==  4 ) {
    vector<string> arguments( argv+1, argv+argc );
    inFile = arguments[0];
    outFile = arguments[1];
    nEvents = atoi(arguments[2].c_str());
  }
  else if ( argc ==  1 ) {
    vector<string> arguments( argv+1, argv+argc );
    inFile = "production_pAu200_2015/MB/pAu_2015_200_MB*.root";
    outFile = "out/pAu.root";
    nEvents = 10000;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in MB: 59388132
  InitReader( Reader, Chain, numEvents );
  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-2,2, 16,-pi,pi );
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicity (per event);# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Multiplicity (per run);Run no.;# of Towers", 60,16120000,16160000, 700,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 60,16120000,16160000, 200,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 700,0,700, 200,0,200);
  TTree *MBtree = new TTree( "MBTree", "MBtree" );
  TTree *MBtowers = new TTree( "MBTowers", "MBtowers" );
  TTree *MBtracks = new TTree( "MBTracks", "MBtracks" );
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, towID, Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz, towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz;

  MBtree->Branch("EventID", &EventID);            MBtree->Branch("RunID", &RunID);          MBtree->Branch("Vx", &Vx);
  MBtree->Branch("Vy", &Vy);                             MBtree->Branch("Vz", &Vz);                     MBtree->Branch("nTowers", &nTowers);
  MBtree->Branch("nPrimary", &nPrimary);         MBtree->Branch("nGlobal", &nGlobal);     MBtree->Branch("nVertices", &nVertices);
  MBtree->Branch("refMult", &refMult);              MBtree->Branch("gRefMult", &gRefMult);
  MBtree->Branch("BbcCoincidenceRate", &BbcCoincidenceRate);                    MBtree->Branch("BbcEastRate", &BbcEastRate);
  MBtree->Branch("BbcWestRate", &BbcWestRate);                                           MBtree->Branch("vpdVz", &vpdVz);
  
  MBtowers->Branch("EventID", &EventID);	  MBtowers->Branch("RunID", &RunID);  	  MBtowers->Branch("towEt", &towEt);
  MBtowers->Branch("towEta", &towEta);	          MBtowers->Branch("towPhi", &towPhi);	  MBtowers->Branch("towID", &towID);

  MBtracks->Branch("EventID", &EventID);	  MBtracks->Branch("RunID", &RunID);	  MBtracks->Branch("nHitsPoss", &nHitsPoss);
  MBtracks->Branch("nHitsFit", &nHitsFit);	  MBtracks->Branch("trPx", &trPx);   	  MBtracks->Branch("trPy", &trPy);
  MBtracks->Branch("trPz", &trPz);           	  MBtracks->Branch("trPt", &trPt);    	  MBtracks->Branch("trEta", &trEta);
  MBtracks->Branch("trPhi", &trPhi);        	  MBtracks->Branch("DCA", &DCA);

  vector<PseudoJet> rawParticles, rawJets;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;

  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();  rawJets.clear();  //  clear containers
    
    Reader.PrintStatus(5); 

    eID = Reader.GetNOfCurrentEvent();
    rID = header->GetRunId();

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    EventID = eID;        RunID = rID;
    int npt = header->GetNOfPrimaryTracks();      nPrimary = npt;
    int ntow = header->GetNOfTowers();               nTowers = ntow;
    Vx = header->GetPrimaryVertexX();    Vy = header->GetPrimaryVertexY();    Vz = header->GetPrimaryVertexZ();
    nGlobal = header->GetNGlobalTracks();                                 nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();                        gRefMult = header->GetGReferenceMultiplicity();
    BbcCoincidenceRate = header->GetBbcCoincidenceRate();     vpdVz = header->GetVpdVz();
    BbcEastRate = header->GetBbcEastRate();                               BbcWestRate = header->GetBbcEastRate();

    hVertex->Fill( Vx, Vy, Vz );
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    
    for ( int i=0; i<npt; ++i ) {
      double primTrackPt = (double) event->GetPrimaryTrack(i)->GetPt();          trPt = primTrackPt;
      double primTrackPx = (double) event->GetPrimaryTrack(i)->GetPx();         trPx = primTrackPx;
      double primTrackPy = (double) event->GetPrimaryTrack(i)->GetPy();         trPy = primTrackPy;
      double primTrackPz = (double) event->GetPrimaryTrack(i)->GetPz();         trPz = primTrackPz;
      double primTrackEta = (double) event->GetPrimaryTrack(i)->GetEta();       trEta = primTrackEta;
      double primTrackPhi = (double) event->GetPrimaryTrack(i)->GetPhi();       trPhi = primTrackPhi;
      
      nHitsFit = event->GetPrimaryTrack(i)->GetNOfFittedHits();
      nHitsPoss = event->GetPrimaryTrack(i)->GetNOfPossHits();
      DCA = event->GetPrimaryTrack(i)->GetDCA();
      
      hPrimaryTracks->Fill( primTrackPt, primTrackEta, primTrackPhi );
      MBtracks->Fill();
    }



    for ( int i=0; i<ntow; ++i ) {
      towEt = event->GetTower(i)->GetEt();
      towEta = event->GetTower(i)->GetEta();
      towPhi = event->GetTower(i)->GetPhi();
      towID = event->GetTower(i)->GetID();
      MBtowers->Fill();
    }
    
    // GatherParticles( container, rawParticles);        //cout<<rawParticles.size()<<endl;


    MBtree->Fill();
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  hPrimaryTracks->Write();
  hVertex->Write();
  hTowersPerEvent->Write();
  hTowersPerRun->Write();
  hPrimaryPerEvent->Write();
  hPrimaryPerRun->Write();
  hnPrimaryVSnTowers->Write();
  
  MBtree->Write();
  MBtowers->Write();
  MBtracks->Write();
  
  pAuFile->Close();
  
  return 0;
}
