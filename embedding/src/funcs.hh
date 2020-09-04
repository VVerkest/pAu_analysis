// Isaac Mooney, WSU - June 2019
// functions.hh

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include "fastjet/contrib/SoftDrop.hh"

// #include "RooUnfold.h"
// #include "RooUnfoldTUnfold.h"

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"
#include "TCanvas.h"

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

//#include "ktTrackEff.hh"
#include <string>
#include <utility>      // std::pair
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <memory>
#include <set>
#include <fstream>

#ifndef funcs_hh
#define funcs_hh

using std::unordered_map; using std::make_shared; using std::shared_ptr; 

namespace Analysis {

  const int nEffEtaBins = 20;
                                        
  int CountTowers( TList* );
  
  bool DiscardpAuEmbedEvent(const TString, const std::vector<fastjet::PseudoJet>, const std::vector<fastjet::PseudoJet> );
  
  void DrawText(const char*, float, float, int);

  std::vector<fastjet::PseudoJet> GatherParticles( TStarJetVectorContainer<TStarJetVector> *, std::vector<fastjet::PseudoJet> & );

  void GenerateWeightedPtResponse( TH1D *[3], TH1D *[21], TH1D *, TString );

  void GetEmbeddingHistograms( TFile *, TH2D *, TH1D *, TH1D *, TString );

  bool HasEnding (std::string const &, std::string const &);  // Helper to build the TChain, used to decide which input format

  void InitReader( TStarJetPicoReader &, TChain*, int, double, double, double, double, double, double, int, double, double, double, TString );
  
  double LookupRun6Xsec(TString);
  double LookupRun12Xsec(TString);
  double LookupRun15Xsec(TString);

  void MissesFakesAndMatches(std::vector<fastjet::PseudoJet> &, std::vector<fastjet::PseudoJet> &,
			     std::vector<fastjet::PseudoJet> &, std::vector<fastjet::PseudoJet> &);

  void ProjectAndSaveFinalUEPlots( TH1D *, TString, TString );
  
  void ProjectAndScaleUEHistogramForAllPt( TH2D *, TH1D *, TH1D *[55], TH1D *[nPtBins], TString );
  
  void ProjectPartLevelJetPt( TH2D*, TH1D*[21], TString );

  TH2D *ProjectUEHistograms( TH3D *, TString );

  void TrackingEfficiency2DCorrection( TH2D*, TH2D*, TH1D * );

  bool UseTriggerTower( int );

  void WeightUEPtByLeadPtAndFakes( TH1D *[nPtBins], TH1D *[55], TH1D *[nPtBins], TH1D *, TH1D *, TString );
  
}

#endif
