//  pAuFunctions.cxx
//  Veronica Verkest       May 15, 2019

#include "pAuFunctions.hh"

namespace pAuAnalysis {

  double ConvertPhi( double &phi ) {
    double pi = 3.14159265;
    if ( phi<0 ) { phi += 2*pi; return phi; }
    else if ( phi>2*pi ) { phi -= 2*pi; return phi; }
    else if ( phi>0 && phi<2*pi ) { return phi; }
    else { std::cerr<<"FUNTION 'ConvertPhi' IS BROKEN"<<std::endl; return -1; }
  }
  

  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & rawParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );

      std::cout<<sv->GetCharge()<<std::cout<<"      "<<current.get_user_index()<<std::endl;
      
      if ( std::abs(current.eta()) > etaCut )      { continue; }  // removes particles with |eta|>1
      if ( current.pt() < partMinPt )      { continue; }  // removes particles with pt<0.2GeV

      rawParticles.push_back(current);
    }
    return rawParticles;
  }

  
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( false );
    //  reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( 30.0 );
    evCuts->SetRefMultCut( 0 );
    evCuts->SetMaxEventPtCut( 30.0 );
    evCuts->SetMaxEventEtCut( 30.0 );
    evCuts->SetMinEventEtCut( 0.2 );
    evCuts->SetVertexZDiffCut( 3.0 );
    
    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( 3.0 );
    trackCuts->SetMinNFitPointsCut( 20 );
    trackCuts->SetFitOverMaxPointsCut( 0.52 );
    trackCuts->SetMaxPtCut ( 30.0 );

    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( 30.0 );
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
