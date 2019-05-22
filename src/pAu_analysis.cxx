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
    nEvents = 1000000;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-2,2, 16,-pi,pi );

  TTree *MBtree = new TTree( "MBTree", "MBtree" );
  TTree *MBtowers = new TTree( "MBTowers", "MBtowers" );
  TTree *MBtracks = new TTree( "MBTracks", "MBtracks" );

  int RunID, EventID, nTowers, nPrimary, nTracks;
  vector<int> Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz;
  vector<double> towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA;

  MBtree->Branch("EventID", &EventID);      MBtree->Branch("RunID", &RunID);                MBtree->Branch("Vx", &Vx);
  MBtree->Branch("Vy", &Vy);                       MBtree->Branch("Vz", &Vz);                           MBtree->Branch("nTowers", &nTowers);
  MBtree->Branch("nTracks", &nTracks);      MBtree->Branch("nPrimary", &nPrimary);

  MBtowers->Branch("EventID", &EventID);	  MBtowers->Branch("RunID", &RunID);  	  MBtowers->Branch("towEt", &towEt);
  MBtowers->Branch("towEta", &towEta);	          MBtowers->Branch("towPhi", &towPhi);

  MBtracks->Branch("EventID", &EventID);	  MBtracks->Branch("RunID", &RunID);	  MBtracks->Branch("nHitsPoss", &nHitsPoss);
  MBtracks->Branch("nHitsFit", &nHitsFit);	  MBtracks->Branch("trPx", &trPx);   	  MBtracks->Branch("trPy", &trPy);
  MBtracks->Branch("trPz", &trPz);           	  MBtracks->Branch("trPt", &trPt);    	  MBtracks->Branch("trEta", &trEta);
  MBtracks->Branch("trPhi", &trPhi);        	  MBtracks->Branch("DCA", &DCA);

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in MB: 59388132
  InitReader( Reader, Chain, numEvents );

  vector<PseudoJet> rawParticles, rawJets;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;

  

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();  rawJets.clear();  //  clear containers
    Reader.PrintStatus(20); 

    eID = Reader.GetNOfCurrentEvent();
    rID = header->GetRunId();

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    EventID = eID;
    RunID = rID;
    int npt = header->GetNOfPrimaryTracks();      nPrimary = npt;
    int ntow = header->GetNOfTowers();               nTowers = ntow;

    
    for ( int i=0; i<npt; ++i ) {
      double primTrackPt = (double) event->GetPrimaryTrack(i)->GetPt();         trPt.push_back(primTrackPt);
      double primTrackPx = (double) event->GetPrimaryTrack(i)->GetPx();         trPx.push_back(primTrackPx);
      double primTrackPy = (double) event->GetPrimaryTrack(i)->GetPy();         trPy.push_back(primTrackPy);
      double primTrackPz = (double) event->GetPrimaryTrack(i)->GetPz();         trPz.push_back(primTrackPz);
      double primTrackEta = (double) event->GetPrimaryTrack(i)->GetEta();      trEta.push_back(primTrackEta);
      double primTrackPhi = (double) event->GetPrimaryTrack(i)->GetPhi();      trPhi.push_back(primTrackPhi);

      nHitsFit.push_back(event->GetPrimaryTrack(i)->GetNOfFittedHits());
      nHitsPoss.push_back(event->GetPrimaryTrack(i)->GetNOfPossHits());
      DCA.push_back(event->GetPrimaryTrack(i)->GetDCA());
      
      hPrimaryTracks->Fill( primTrackPt, primTrackEta, primTrackPhi );
    }



    for ( int i=0; i<ntow; ++i ) {
      
    }
    
    GatherParticles( container, rawParticles);        //cout<<rawParticles.size()<<endl;


    MBtree->Fill();
    MBtowers->Fill();
    MBtracks->Fill();
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  hPrimaryTracks->Write();
  MBtree->Write();
  MBtowers->Write();
  MBtracks->Write();
  pAuFile->Close();
  
  return 0;
}
