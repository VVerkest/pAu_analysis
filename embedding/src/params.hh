// Isaac Mooney, WSU - June 2019
// Parameters for jet mass analysis

//#include "TStopwatch.h"
//#include "ktTrackEff.hh"
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <utility>      // std::pair
#include <string>
#include <iostream>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TRandom.h>
#include "TLegend.h"
#include "THStack.h"


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

// Define a namespace for the variables

#ifndef PARAMS_HH
#define PARAMS_HH

#define __ERR(message) {std::cerr << "[" << __FILE__ << "::" << __func__ << "()] -- ERR: " << message << std::endl;}
#define __OUT(message) {std::cout << "[" << __FILE__ << "::" << __func__ << "()] -- OUT: " << message << std::endl;}

namespace Analysis {

  const double Pi = 3.141592653;

  const int xbins = 55;
  const double xbinEdge[xbins+1] = {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0};
  const int zbins = 20;
  const double zbinEdge[zbins+1] = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  const int ybins = 14;
  const double ybinEdge[ybins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0 };
  
  const std::string det_badTowers = "lists/bad_towers_pAu2015.list";

  const double detVertexZCut = 10.0;
  const int RefMultCut = 0;
  const double detMaxEventPtCut = 30.0;
  const double detMaxEventEtCut = 30.0;
  const double detMinEventEtCut = 0.2;
  const double detVertexZDiffCut = 6.0;
  const double detDCACut = 1.0;
  const int detMinNFitPointsCut = 20;
  const double detFitOverMaxPointsCut = 0.52;
  const double detMaxPtCut = 30.0;
  const double detMaxEtCut = 30.0;

  const double R = 0.4;
  
  const double maxTrackEta = 1.0;
  const double maxJetEta = 1.0 - R;
  const double minJetPt = 5.0;
  const double partMinPt = 0.2;


  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const int ptMarker[nPtBins] = { 20, 34, 33 };

  const double AREA = 4*(1.14159265);
  const double eastArea = 2*(0.7)*(fastjet::pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(fastjet::pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(fastjet::pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  
  const int marker[55] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
  
  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString emw[nEtaBins] = { "east", "mid", "west" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const double area[nEtaBins] = { eastArea, midArea, westArea };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };
  const int nEAbins = 2;
  const TString EAbinName[nEAbins] = { "Lo", "Hi" };
  const TString EAbinString[nEAbins] = { "Low EA", "High EA" }; 
  const TString eastmidwest[nEtaBins] = { "East", "Mid", "West" };

  const double NEF_max = 0.9;             //neutral energy fraction of jet must be < 90% (not used for PYTHIA) !!!!!

  
  /*
  const std::string outFileName = "out/analysis.root";
  const double nEvents = -1;       //  NUMBER OF EVENTS  (-1 runs all)

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~consts~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                                                                                                                                                                        
  const double chPionMass = 0.13957018;
  const double Pi0Mass = 0.1349766;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//                 

  //~~~~~~~~~~~~~~~~~~~~~~~~~~triggerIDs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  const int tppJP2 = 370621; //ppY12 JP2
  const int tppHT2a = 370521, tppHT2b = 370522, tppHT2c = 370531; //ppY12 BHT2*BBCMB, BHT2*BBCMB, BHT2
  const int tppVPDMB_nobsmd = 370011; //ppY12 VPDMB-nobsmd
  const int tpAuJP2a = 500401, tpAuJP2b = 500411; //pAuY15 JP2 //CAUTION USING THIS - ABNORMAL JP2 THRESHOLD FOR NEGATIVE RAPIDITY
  const int tpAuHT2a = 500205, tpAuHT2b = 500215; //pAuY15 BHT2*BBCMB
  const int tpAuBBCMBa = 500008, tpAuBBCMBb = 500018; //pAuY15 BBCMB
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //NOTE: "det" means detector-level (i.e. both Geant and data), while "dat" means only data.
  //Similarly, "truth" means particle-level (i.e. Pythia or Herwig), while "sim" means simulation (i.e. Pythia, Herwig, Geant)
  
  //~~~~~~~~~~~~~~~~~~~~~~~pp quality cuts~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  //bad runs and towers
  
  const int refMultCut = 0;

  //detector: event, track, tower cuts
  const std::string det_triggerString = "All"; //these trigger strings are not really used!! It's done by hand with the triggerIDs
  const double det_absMaxVz = 30.0;//cm //|Vz|<=30 cm
  const double det_vZDiff = 6.0;//cm //max diff btwn selected TPC vertex & most probable VPD vertex (in ppRun6 VPD vz = 0, so vZDiff should be > absMaxVz)
  const double det_evEtMin = -1;//GeV
  const double det_evEtMax = 30.0;//GeV
  const double det_evPtMax = 30.0;//GeV
  const double dat_maxEtTow = 9999;//GeV
  const double det_DCA = 1.0;//cm //The 2012 pp Picos are already limited to 2 and the 2015 pA to 3
  const double det_NFitPts = 20;
  const double det_FitOverMaxPts = 0.52;

  //particle cuts
  const double max_track_rap = 1.0;       //detector acceptance - mimicked in Pythia as well
  const double partMinPt = 0.2;           //30.0 GeV >= particle pT >= 0.2 GeV 
  const double partMaxPt = 30.0;          

  //jet cuts
  const double R = 0.4;                   //jet resolution parameter                 
  const double jet_ptmin = 5.0;//GeV      //gen-jet pT >= 5.0 GeV                 
  const double det_jet_ptmin = 5.0;//GeV //detector-level jet pT >= 15 GeV
  const double jet_ptmax = 1000.0;//GeV   //DEBUG
  const double max_rap = max_track_rap-R; //|eta_jet| < 1-R
  const double mass_min = 1.0;//GeV       //det-jet M >= 1.0 GeV !!!!!
  
  //ghosts - not used for pp (jets aren't bkground subtracted), so once I have AuAu running as well, this can be moved there
  const int ghost_repeat = 1;
  const double ghost_area = 0.01;
  const double ghost_maxrap = max_rap + 2.0 * R;

  //Soft Drop params
  const double z_cut = 0.10;
  const double Beta = 0.0;
  const double R0 = 1.0;
  
  //~~~~~~~~~~~~~~~~~~~~~~~~pAu quality cuts~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  const std::string pAu_triggerString = "All";
  const std::string pAu_badTowers = "lists/bad_towers_pAu2015.list";//preliminary; once calibrations and status tables are redone, list will likely be updated
  const std::string pAu_bad_run_list = "lists/dummy_badrun.list";
  
  //event, track, tower cuts
  const double pAu_vZDiff = 3.0; //this may be too tight of a cut - to be revisited
  const double pAu_BBCE_ADC_sum_max = 64000.0; //tower multiplicity in event explodes above this number. Detectors probably saturated.
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  */
}

#endif
