// Veronica Verkest
// October 13, 2020

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
                                        
  void DrawText(const char*, float, float, int);

  double CalculateGeometricMean( TH1D* );

  TH1D* GenerateFractionalContribution( TH2D*, double, double, TString, TString );
  
  void GenerateWeightedPtResponse( TH1D *[3], TH1D *[21], TH1D *, TString );

  void GetEmbeddingHistograms( TFile *, TH2D *, TH1D *, TH1D *, TString );
  
  void ProjectAndPlotByEta( TH2D *, TH1D *, int, int, int, TString, TFile* );

  void ProjectAndSaveFinalUEPlots( TH1D *, TString, double, TString );

  void ProjectPartLevelJetPt( TH2D*, TH1D *[21], TString );
  
  // void StackAndSaveNchPlots( TH1D *[3][3], TH2D *, TString, TString );

  void StackAndSavePtPlots( TH1D *[3][3], TH2D *, TString, TString );

  void ProjectScaleAndSaveUE1D( TH2D*, TH1D &, TH1D*, TH1D &, double, double, TString, TString );

  TString RoundDecimal(double, int );
  
  void TrackingEfficiencyByPtAndEta( TH2D* [nPtBins], TH2D*[nPtBins], TFile *, TString, TString );

  void WeightAndAddCorrected2Ds( TH2D *, TH1D *, TH2D *[55], TString );

  TH1D* WeightAndSumByFC1D( TH1D*, TH1D*[55] );
  
  void WeightUEPtByLeadPtAndFakes( TH1D *[nPtBins], TH1D *[55], TH1D *[nPtBins], TH1D *, TH1D *, TString, TString );
  
}

#endif
