
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

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

#ifndef pAuFunctions_hh
#define pAuFunctions_hh

namespace pAuAnalysis {

  const int RefMultCut = 0;
  const int MinNFitPointsCut = 20;

  const double VertexZCut = 30.0;
  const double VertexZDiffCut = 3.0;
  const double MaxEventPtCut = 30.0;
  const double MaxEventEtCut = 30.0;
  const double MinEventEtCut = 0.2;
  const double DCACut = 1.0;
  const double FitOverMaxPointsCut = 0.52;
  const double MaxPtCut = 30.0;
  const double MaxEtCut = 30.0;
  
  const double R = 0.4;
  const double vzCut = 30.0;   // |Vz|<=30 cm
  const double dcaCut = 3.0;   // DCA<=3.0 cm
  const double etaCut = 1.0;       // |eta|<=1
  const double partMinPt = 0.2;      //  particle Pt >= 0.2 GeV
  const double jetMinPt = 2.0;      //  Jet Pt >= 2.0 GeV
  const double qpi = 3.141592653589793238462643383279502884197/4;

  const double AREA = 4*(fastjet::pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  
  const int nEtaBins = 4;
  const double etaBinLo[nEtaBins] = { -1.0, -0.5, 0.0, 0.5 };
  const double etaBinHi[nEtaBins] = { -0.5, 0.0, 0.5, 1.0 };
  const TString etaBinName[nEtaBins] = { "_n10_n05", "_n05_00", "_00_05", "_05_10" };
  const TString etaBinString[nEtaBins] = { "-1.0<#eta<-0.5", "-0.5<#eta<0.0", "0.0<#eta<0.5", "0.5<#eta<1.0" };


  void BackGroundEstimation( std::vector<fastjet::PseudoJet> *chgPart, std::vector<fastjet::PseudoJet> *neuPart, fastjet::PseudoJet *leadJet,
			     TH3D *PartPtDEtaDPhi, TH3D *PartPtEtaPhi, TH3D *CHARGED, TH3D *NEUTRAL, TH3D *BG, double chgSum, double neuSum );

  int CountTowers( TList *selectedtowers );
    
  std::vector<fastjet::PseudoJet> GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & rawParticles );

  std::vector<fastjet::PseudoJet> GatherChargedBG (  fastjet::PseudoJet trigJet, TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & chgParticles );

  std::vector<fastjet::PseudoJet> GatherNeutralBG (  fastjet::PseudoJet trigJet, TStarJetVectorContainer<TStarJetVector> * container , std::vector<fastjet::PseudoJet> & newParticles );

  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents );

  bool UseEvent( TStarJetPicoEventHeader* Header, double vz_cut, double vz );

}

#endif
