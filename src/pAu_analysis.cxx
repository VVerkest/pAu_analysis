//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019

#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

vector<PseudoJet> rawParticles, rawJets;
int ID;


int main () {


  
  TChain* chain = new TChain( "JetTree" );         chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  TStarJetPicoReader reader;

  // set the chain
  reader.SetInputChain( chain );
  // apply hadronic correction - subtract 100% of charged track energy from towers
  reader.SetApplyFractionHadronicCorrection( false );
  reader.SetFractionHadronicCorrection( 0.9999 );
  reader.SetRejectTowerElectrons( kFALSE );
    
  // Event and track selection
  // -------------------------
    
  TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
  evCuts->SetVertexZCut ( 1000 );
  evCuts->SetRefMultCut( 0 );
  evCuts->SetMaxEventPtCut( 1000 );
  evCuts->SetMaxEventEtCut( 1000 );

  // Tracks cuts
  TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
  trackCuts->SetDCACut( 100 );
  trackCuts->SetMinNFitPointsCut( -1 );
  trackCuts->SetFitOverMaxPointsCut( -1 );    
    
  std::cout << "Using these track cuts:" << std::endl;
  std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
  std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
  std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
  // Towers
  TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
  towerCuts->SetMaxEtCut( 9999.0 );
  towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );
    
  // V0s: Turn off
  reader.SetProcessV0s(false);
    
  // Initialize the reader
  reader.Init( 1000 ); //runs through all events with -1


  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;





  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( reader.NextEvent() ) {

    rawParticles.clear();  rawJets.clear();  //  clear containers

    ID = reader.GetNOfCurrentEvent();
    reader.PrintStatus(20); 

    event = reader.GetEvent();    header = event->GetHeader();

    container = reader.GetOutputContainer();

    GatherParticles ( container, rawParticles);

    cout<<rawParticles.size()<<endl;
    
  }
    
  return 0;
}
