//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace fastjet;
using namespace std;




int main () {

  // TString inFileName = "pAu_2015_200_MB_156_160_2.root";
  // TFile* inFile = new TFile( inFileName, "READ" );

  cout << "uuuuuuggghhhhhh"<<endl;
  
  // TChain* Chain = new TChain( "JetTree" );         Chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  // TStarJetPicoReader Reader;
  // InitReader( Reader, Chain, -1 );
  // TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  // TStarJetVectorContainer<TStarJetVector> * container;



  

  return 0;
}


  // //  INITIATE READER
  // void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents ) {

  //   std::string collisionType = "pAu";
      
  //   // First tolower() on the analysisType
  //   // shouldnt be necessary....
  //   std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
  //   // set the chain
  //   reader.SetInputChain( chain );
  //   // apply hadronic correction - subtract 100% of charged track energy from towers
  //   reader.SetApplyFractionHadronicCorrection( true );
  //   reader.SetFractionHadronicCorrection( 0.9999 );
  //   reader.SetRejectTowerElectrons( kFALSE );
    
  //   // Event and track selection
  //   // -------------------------
    
  //   TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
  //   evCuts->SetVertexZCut ( 1000 );
  //   evCuts->SetRefMultCut( 0 );
  //   evCuts->SetMaxEventPtCut( 1000 );
  //   evCuts->SetMaxEventEtCut( 1000 );

  //   // Tracks cuts
  //   TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
  //   trackCuts->SetDCACut( 100 );
  //   trackCuts->SetMinNFitPointsCut( -1 );
  //   trackCuts->SetFitOverMaxPointsCut( -1 );    
    
  //   std::cout << "Using these track cuts:" << std::endl;
  //   std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
  //   std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
  //   std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
  //   // Towers
  //   TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
  //   towerCuts->SetMaxEtCut( 9999.0 );
  //   towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );

  //   std::cout << "Using these tower cuts:" << std::endl;
  //   std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
  //   std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
  //   // V0s: Turn off
  //   reader.SetProcessV0s(false);
    
  //   // Initialize the reader
  //   reader.Init( nEvents ); //runs through all events with -1
  // }
