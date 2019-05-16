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

const int pi = 3.14159;

vector<PseudoJet> rawParticles, rawJets;
int ID;

TH3D *PrimaryTracks = new TH3D( hPrimaryTracks, "Primary Tracks: p_{T}, #{eta}, and #{phi}", 50,0,20, 20,-2,2, 20,0,2*pi );

int main () {

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();
  
  TChain* Chain = new TChain( "JetTree" );         Chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  TStarJetPicoReader Reader;                               int numEvents = 1000;
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
      double trackPt = (double*) event->GetPrimaryTrack(i)->GetPt();
    }


    
    GatherParticles ( container, rawParticles);        //cout<<rawParticles.size()<<endl;
    
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  return 0;
}
