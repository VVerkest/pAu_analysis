#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"
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
#include <vector>
#include <algorithm>
#include <iterator>

#ifndef pAuParameters_hh
#define pAuParameters_hh

namespace pAuAnalysis {

const int RefMultCut = 0;
const int MinNFitPointsCut = 20;

const double vzCut = 30.0;   // |Vz|<=30 cm
const double VertexZDiffCut = 6.0;
const double MaxEventPtCut = 30.0;
const double MaxEventEtCut = 30.0;
const double MinEventEtCut = 0.2;
const double DCACut = 1.0;
const double FitOverMaxPointsCut = 0.52;
const double MaxPtCut = 30.0;
const double MaxEtCut = 30.0;

const double R = 0.4;
const double etaCut = 1.0;       // |eta|<=1
const double partMinPt = 0.2;      //  particle Pt >= 0.2 GeV
const double jetMinPt = 2.0;      //  Jet Pt >= 2.0 GeV
const double qpi = 3.141592653589793238462643383279502884197/4;
const double AREA = 4.*(fastjet::pi/3.);

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };

  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  
  const double eastArea = 2.*(0.7)*(fastjet::pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )
  const double midArea = 2.*(0.6)*(fastjet::pi/3.);   // (  0.6 in eta  ) X (  2pi/3 in phi  )
  const double westArea = 2.*(0.7)*(fastjet::pi/3.);   // (  0.7 in eta  ) X (  2pi/3 in phi  )

  const int nChgBins = 3;
  const TString BackgroundChargeBias[nChgBins] = { "_chgBG", "_neuBG", "_allBG" };
  const TString BackgroundChargeString[nChgBins] = { "Charged", "Neutral", "Chg+Neu" };
  const int color[nChgBins] = { 807, 823, 874 };
  const int marker[nChgBins] = { 22, 23, 20 };

  const double leadJetMinPt = 10.0;
  //const double leadJetMaxPt = 30.0;
  
}

#endif
