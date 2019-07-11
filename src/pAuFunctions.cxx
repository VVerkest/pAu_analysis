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


  void InitHistograms( ) {
    TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
    TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Frequency;# of Towers", 700,0,700 );
    TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Frequency (per run);Run no.;# of Towers", 400,16124000,16164000, 140,0,700 );
    TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
    TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 4000,16124000,16164000, 40,0,200 );
    TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;# Primary Tracks", 140,0,700, 40,0,200);
    TH2D *hPrimaryVsGlobal = new TH2D("hPrimaryVsGlobal","# Primary Tracks vs. # Global Tracks;# Global Tracks;# Primary Tracks", 150,0,3000, 150,0,150 );
    TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 3500,0,3500000, 40,0,200 );
    TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 3500,0,3500000, 150,0,3000 );
    TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 140,0,7000000, 40,0,200 );
    TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 80,0,80000, 150,0,3000 );
    TH2D *hPrimaryVsBBCsumE = new TH2D("hPrimaryVsBBCsumE","# Primary Tracks vs. BBC ADC East Sum;BBC ADC East Sum;# Primary Tracks", 160,0,80000, 40,0,200 );
    TH2D *hTowersVsBBCsumE = new TH2D("hTowersVsBBCsumE","# Towers vs. BBC ADC East Sum;BBC ADC East Sum;# Towers", 160,0,80000, 140,0,700 );
    TH2D *hChgVsNeuBG = new TH2D("hChgVsNeuBG","Charged vs. Neutral Background;Neutral #rho;Charged #rho",100,0,25,100,0,25);
    TH3D *hAllJetsPtEtaPhi = new TH3D( "hAllJetsPtEtaPhi", "Inclusive Jets p_{T}, #eta, #phi;Jet p_{T} (GeV);Jet #eta;Jet #phi", 400,0.0,100.0, 40,-1.0,1.0, 120,0.0,2*pi  );
    TH3D *hAllJetsPtRhoEta = new TH3D( "hAllJetsPtRhoEta", "Inclusive Jets p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );
    TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#phi;#eta", 280,0,70, 120,0,6.3, 40,-1.0,1.0);
    TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 120,0,2*pi, 40,-1.0,1.0);
    TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 120,0,2*pi, 40,-1.0,1.0);
    TH2D *hPrimaryVsRho = new TH2D("hPrimaryVsRho","# Primary Tracks vs. Underlying Event;#rho (GeV);# Primary Tracks", 100,0,25, 40,0,200);
    TH2D *hGlobalVsRho = new TH2D("hGlobalVsRho","# Global Tracks vs. Underlying Event;#rho (GeV);# Global Tracks", 100,0,25, 150,0,3000 );
    TH3D *hPt_UE_RefMult = new TH3D("hPt_UE_RefMult","UE vs. Ref. Mult;Lead Jet p_{T} (GeV);Underlying Event (GeV);Ref. Mult.", 500,0,125, 50,0,25, 100,0,1000 );  
    TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 500,0,125, 50,0,25, 140,0,7000000 );
    TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 160,0,80000 );
    TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 100,0,25, 140,0,700 );
    TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 140,0,35, 280,0,70);
    TH3D *hBG = new TH3D("hBG","Background Particle #eta-#phi vs. p_{T};Lead Jet p_{T}(GeV);Particle #eta;Particle #phi", 280,0,70, 40,-1,1, 120,0,2*pi );
    TH3D *hCHARGED = new TH3D("hCHARGED","CHARGED: Background Particle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
    TH3D *hNEUTRAL = new TH3D("hNEUTRAL","NEUTRAL: Background Particle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
    TH3D *hPartPtDEtaDPhi = new TH3D("hPartPtDEtaDPhi","Background Particle p_{T} vs. #Delta#eta vs. #Delta#phi;Particle p_{T} (GeV);#Delta#eta;#Delta#phi", 120,0,30, 80,-2.0,2.0, 120,-pi,pi );
    TH3D *hPartPtEtaPhi = new TH3D("hPartPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;Lead Jet p_{T} (GeV);Particle #eta;Particle #phi", 120,0,30, 40,-1.0,1.0, 120,0,2*pi );
    TH3D *hAllPtEtaPhi = new TH3D("hAllPtEtaPhi","All Particles p_{T} vs. #eta vs. #phi;Particle p_{T} (GeV);Particle #eta;Particle #phi", 120,0,30, 40,-1.0,1.0, 120,0,2*pi );
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

}
