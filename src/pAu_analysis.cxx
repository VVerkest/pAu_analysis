//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019


//for pAu, need to ask if the event has the trigger. The trigger IDs are:
//HT2*BBCMB : 500205, 500215
//JP2 : 500401, 500411
//BBCMB : 500008, 500018
//VPDMB :  500904


#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

const double pi = 3.14159;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-2,2, 16,-pi,pi );


  TTree *MBtree, *MBtowers, *MBtracks;
  int RunID, EventID, nTowers, nPrimary;
  double Vx, Vy, Vz;
  vector<int> nTracks, Charge, nHits;
  vector<double> towE, towEta, towPhi, trPx, trPy, trPz, DCA;
  
  MBtree->Branch("EventID", &EventID_MB);
  MBtree->Branch("RunID", &RunID_MB);
  MBtree->Branch("Vx", &Vx_MB);
  MBtree->Branch("Vy", &Vy_MB);
  MBtree->Branch("Vz", &Vz_MB);
  MBtree->Branch("nTowers", &nTowers_MB);
  MBtree->Branch("nPrimary", &nPrimary_MB);

  MBtowers->Branch("EventID", &EventID_MB);
  MBtowers->Branch("RunID", &RunID_MB);
  MBtowers->Branch("nTracks", &nTracks_MB);
  MBtowers->Branch("towE", &towE_MB);
  MBtowers->Branch("towEta", &towEta_MB);
  MBtowers->Branch("towPhi", &towPhi_MB);

  MBtracks->Branch("EventID", &EventID_MB);
  MBtracks->Branch("RunID", &RunID_MB);
  MBtracks->Branch("nHits", &nHits_MB);
  MBtracks->Branch("trPx", &trPx_MB);
  MBtracks->Branch("trPy", &trPy_MB);
  MBtracks->Branch("trPz", &trPz_MB);
  MBtracks->Branch("DCA", &DCA_MB);
  
  TChain* Chain = new TChain( "JetTree" );
  // Chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  Chain->Add( "production_pAu200_2015/MB/pAu_2015_200_MB*.root" );
  TStarJetPicoReader Reader;
  int numEvents = 5000000;        // total events in MB: 59388132
  InitReader( Reader, Chain, numEvents );

  vector<PseudoJet> rawParticles, rawJets;
  int ID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;

  

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();  rawJets.clear();  //  clear containers

    ID = Reader.GetNOfCurrentEvent();
    Reader.PrintStatus(20); 

    event = Reader.GetEvent();    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    int npt = header->GetNOfPrimaryTracks();
    for ( int i=0; i<npt; ++i ) {
      double primTrackPt = (double) event->GetPrimaryTrack(i)->GetPt();
      double primTrackEta = (double) event->GetPrimaryTrack(i)->GetEta();
      double primTrackPhi = (double) event->GetPrimaryTrack(i)->GetPhi();

      hPrimaryTracks->Fill( primTrackPt, primTrackEta, primTrackPhi );
    }


    
    GatherParticles ( container, rawParticles);        //cout<<rawParticles.size()<<endl;
    
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile("out/pAu.root","RECREATE");
  hPrimaryTracks->Write();
  MBtree->Write();
  MBtowers->Write();
  MBtracks->Write();
  
  return 0;
}
