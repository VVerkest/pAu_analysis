//  pAuFunctions.cxx
//  Veronica Verkest       May 15, 2019

#include "pAuFunctions.hh"

namespace pAuAnalysis {

  
//  INITIATE READER
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {

    std::string collisionType = "AuAu";
      
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( true );
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

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }

}
