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

// const int pi = 3.14159;

vector<PseudoJet> rawParticles, rawJets;
int ID;

TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #{eta}, and #{phi};p_{T} (GeV);#{eta};#{phi}", 40,0,20, 40,-2,2, 16r,-3.141,3.141 );


int main () {

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  TChain* Chain = new TChain( "JetTree" );
  // Chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  Chain->Add( "production_pAu200_2015/MB/pAu_2015_200_MB*.root" );
  TStarJetPicoReader Reader;                               int numEvents = -1;
  InitReader( Reader, Chain, numEvents );


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
  
  return 0;
}
