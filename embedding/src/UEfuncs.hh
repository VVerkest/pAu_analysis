#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
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

#ifndef UEfuncs_hh
#define UEfuncs_hh


namespace Analysis {

  void DrawText(const char*, float, float, int);

  TString GetEfficHistoName( TString );
  
}

#endif
