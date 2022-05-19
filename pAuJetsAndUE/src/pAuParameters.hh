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




//  THnSparse binning

const int NspJetDims = 6; // leadPt;leadEta;leadPhi;leadNcons;iBBCEsum;zdcx
const int NspJetBins[NspJetDims] = { 50, 3, 30, 25, 10, 20 };

const double bin_leadPt[NspJetBins[0]+1] = {4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54.};
const double bin_leadEta[NspJetBins[1]+1] = {-1., -0.3, 0.3, 1.};
//{-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
const double bin_leadPhi[NspJetBins[2]+1] = {0., 0.209439, 0.418879, 0.628318, 0.837757, 1.0472, 1.25664, 1.46608, 1.67551, 1.88495, 2.09439, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159, 3.35103, 3.56047, 3.76991, 3.97935, 4.18879, 4.39823, 4.60767, 4.8171, 5.02654, 5.23598, 5.44542, 5.65486, 5.8643, 6.07374, 6.28318};
const double bin_leadNcons[NspJetBins[3]+1] = {0., 4., 8., 12., 16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 64., 68., 72., 76., 80., 84., 88., 92., 96., 100.};
const double bin_iBBCEsum[NspJetBins[4]+1] = {0., 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000.};
const double bin_ZDCx[NspJetBins[5]+1] = {0., 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 9000., 10000., 11000., 12000., 13000., 14000., 15000., 16000., 17000., 18000., 19000., 20000.};

const int NspUEdims = 6;
const int NspUEbins[NspUEdims] = { 50, 15, 20, 30, 20, 10 };

const double bin_UEpt[NspUEbins[1]+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };
const double bin_UEeta[NspUEbins[2]+1] = {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
const double bin_UEphi[NspUEbins[3]+1] = {0., 0.209439, 0.418879, 0.628318, 0.837757, 1.0472, 1.25664, 1.46608, 1.67551, 1.88495, 2.09439, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159, 3.35103, 3.56047, 3.76991, 3.97935, 4.18879, 4.39823, 4.60767, 4.8171, 5.02654, 5.23598, 5.44542, 5.65486, 5.8643, 6.07374, 6.28318};
  
}

#endif
