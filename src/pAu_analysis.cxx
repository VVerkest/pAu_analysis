//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019

#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

vector<PseudoJet> rawParticles, rawJets;
int ID;


int main () {

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

    GatherParticles ( container, rawParticles);

    cout<<rawParticles.size()<<endl;
    
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  return 0;
}
