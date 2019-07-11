//  pAuFunctions.cxx
//  Veronica Verkest       May 15, 2019

#include "pAuFunctions.hh"

namespace pAuAnalysis {
  

  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & rawParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );
      if ( std::abs(current.eta()) > etaCut )      { continue; }  // removes particles with |eta|>1
      if ( current.pt() < partMinPt )      { continue; }  // removes particles with pt<0.2GeV

      rawParticles.push_back(current);
    }
    return rawParticles;
  }

  
  std::vector<fastjet::PseudoJet> GatherChargedBG ( fastjet::PseudoJet trigJet, TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & chgParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      if ( sv->GetCharge() == 0 ) { continue; }

      double absDeltaPhi = fabs( current.delta_phi_to( trigJet ) );
      if ( absDeltaPhi < 1.0  ||  absDeltaPhi > 2.14159265 ) { continue; }       //  1 < delta phi < ( pi - 1 )
      
      if ( std::abs(current.eta()) > etaCut )      { continue; }  // removes particles with |eta|>1
      if ( current.pt() < partMinPt )      { continue; }  // removes particles with pt<0.2GeV

      current.set_user_index( sv->GetCharge() );

      chgParticles.push_back(current);
    }
    return chgParticles;
  }

  
  std::vector<fastjet::PseudoJet> GatherNeutralBG ( fastjet::PseudoJet trigJet, TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & neuParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      if ( sv->GetCharge() != 0 ) { continue; }

      double absDeltaPhi = fabs( current.delta_phi_to( trigJet ) );
      if ( absDeltaPhi < 1.0  ||  absDeltaPhi > 2.14159265 ) { continue; }       //  1 < delta phi < ( pi - 1 )
      
      if ( std::abs(current.eta()) > etaCut )      { continue; }  // removes particles with |eta|>1
      if ( current.pt() < partMinPt )      { continue; }  // removes particles with pt<0.2GeV

      current.set_user_index( sv->GetCharge() );

      neuParticles.push_back(current);
    }
    return neuParticles;
  }


  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    // reader.SetApplyFractionHadronicCorrection( false );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( VertexZCut );
    evCuts->SetRefMultCut( RefMultCut );
    evCuts->SetMaxEventPtCut( MaxEventPtCut );
    evCuts->SetMaxEventEtCut( MaxEventEtCut );
    evCuts->SetMinEventEtCut( MinEventEtCut );
    evCuts->SetVertexZDiffCut( VertexZDiffCut );
    
    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( DCACut );
    trackCuts->SetMinNFitPointsCut( MinNFitPointsCut );
    trackCuts->SetFitOverMaxPointsCut( FitOverMaxPointsCut );
    trackCuts->SetMaxPtCut ( MaxPtCut );

    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( MaxEtCut );
    //towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );
    towerCuts->AddBadTowers( "/nfs/rhi/STAR/Data/P16id/resources/bad_towers_pAu2015.list" );

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }


  bool UseEvent( TStarJetPicoEventHeader* Header, double vz_cut, double vz ) {
    if (Header->GetRunId() >= 16142059 && Header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (Header->GetRunId() == 16135031 || Header->GetRunId() == 16135032) { continue; }
    if (!(Header->HasTriggerId(500401) || Header->HasTriggerId(500411))) {continue;}   //  ONLY SELECT JP2 TRIGGER EVENTS
    if ( abs(vz) > vz_cut ) { continue; }
    if ( Header->GetBbcAdcSumEast() > 64000 ) { continue; } 
  }

}
