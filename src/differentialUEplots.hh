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

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"

using namespace pAuAnalysis;

#ifndef differentialUEplots_hh
#define differentialUEplots_hh

namespace differentialUEplots {

  const int xbins = 55;
  const double xbinEdge[xbins+1] = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0};
  const int zbins = 20;
  const double zbinEdge[zbins+1] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const int ybins = 14;
  const double ybinEdge[ybins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0 };

  
  const int nEAbins = 2;
  const TString EAbinName[nEAbins] = { "Lo", "Hi" };
  const double area[nEtaBins] = { eastArea, midArea, westArea };
  const int EAcolor[nEAbins] = { 884, 810 };
  const int EAmarker[nEAbins] = { 24, 20 };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };
  const TString EAbinString[nEAbins] = { "Low EA", "High EA" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const TString emw[nEtaBins] = { "east", "mid", "west" };
  const int ptMarker[nPtBins] = { 20, 21, 29 };
  const TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };
  
  const int eta_bins = 10;
  const double etabinEdge[eta_bins+1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0 }; // {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

  const int ecolor[eta_bins] = { 618, 633, 807, 800, 819, 419, 433, 862, 884, 619 };
  const TString etabinEdgeString[eta_bins+1] = { "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8", "1.0" };



}

#endif
