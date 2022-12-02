//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 13 14:00:31 2022 by ROOT version 6.26/00
// from TTree events/embedding events
// found on file: hadd_piKp_track_emb.root
//////////////////////////////////////////////////////////

#ifndef Reader_h
#define Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class Reader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxmc_Cjet = 5;
   static constexpr Int_t kMaxCjet = 22;
   static constexpr Int_t kMaxFjet = 25;
   static constexpr Int_t kMaxmc_track = 5;
   static constexpr Int_t kMaxtrack = 110;
   static constexpr Int_t kMaxtower = 120;

   // Declaration of leaf types
   Int_t           mc_Cjet_;
   UInt_t          mc_Cjet_fUniqueID[kMaxmc_Cjet];   //[mc_Cjet_]
   UInt_t          mc_Cjet_fBits[kMaxmc_Cjet];   //[mc_Cjet_]
   Float_t         mc_Cjet_pt[kMaxmc_Cjet];   //[mc_Cjet_]
   Float_t         mc_Cjet_eta[kMaxmc_Cjet];   //[mc_Cjet_]
   Float_t         mc_Cjet_phi[kMaxmc_Cjet];   //[mc_Cjet_]
   Float_t         mc_Cjet_area[kMaxmc_Cjet];   //[mc_Cjet_]
   vector<short>   mc_Cjet_index_track[kMaxmc_Cjet];
   vector<short>   mc_Cjet_index_tower[kMaxmc_Cjet];
   Float_t         mc_Cjet_rho;
   Float_t         mc_Cjet_rho_sigma;
   Int_t           Cjet_;
   UInt_t          Cjet_fUniqueID[kMaxCjet];   //[Cjet_]
   UInt_t          Cjet_fBits[kMaxCjet];   //[Cjet_]
   Float_t         Cjet_pt[kMaxCjet];   //[Cjet_]
   Float_t         Cjet_eta[kMaxCjet];   //[Cjet_]
   Float_t         Cjet_phi[kMaxCjet];   //[Cjet_]
   Float_t         Cjet_area[kMaxCjet];   //[Cjet_]
   vector<short>   Cjet_index_track[kMaxCjet];
   vector<short>   Cjet_index_tower[kMaxCjet];
   Float_t         Cjet_rho;
   Float_t         Cjet_rho_sigma;
   Int_t           Fjet_;
   UInt_t          Fjet_fUniqueID[kMaxFjet];   //[Fjet_]
   UInt_t          Fjet_fBits[kMaxFjet];   //[Fjet_]
   Float_t         Fjet_pt[kMaxFjet];   //[Fjet_]
   Float_t         Fjet_eta[kMaxFjet];   //[Fjet_]
   Float_t         Fjet_phi[kMaxFjet];   //[Fjet_]
   Float_t         Fjet_area[kMaxFjet];   //[Fjet_]
   vector<short>   Fjet_index_track[kMaxFjet];
   vector<short>   Fjet_index_tower[kMaxFjet];
   Float_t         Fjet_rho;
   Float_t         Fjet_rho_sigma;
 //mupicoEventHeader_v2 *mu_event;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           runId;
   Int_t           eventId;
   Float_t         ZDCx;
   Double_t        BBCx;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           BBC_Ein;
   Int_t           BBC_Eout;
   Int_t           BBC_Win;
   Int_t           BBC_Wout;
   Float_t         vzVpd;
   Float_t         ranking;
   Int_t           ZdcSumAdcEast;
   Int_t           ZdcSumAdcWest;
   Short_t         ZdcSmdEastHorizontal[8];
   Short_t         ZdcSmdEastVertical[8];
   Short_t         ZdcSmdWestHorizontal[8];
   Short_t         ZdcSmdWestVertical[8];
   Short_t         EastBBC[24];
   Short_t         WestBBC[24];
   Short_t         refMult;
   Int_t           nGlobalTracks;
   Short_t         pthat_min;
   Short_t         pthat_max;
   Short_t         pthat_bin;
   Int_t           mc_track_;
   UInt_t          mc_track_fUniqueID[kMaxmc_track];   //[mc_track_]
   UInt_t          mc_track_fBits[kMaxmc_track];   //[mc_track_]
   Short_t         mc_track_geantId[kMaxmc_track];   //[mc_track_]
   Short_t         mc_track_id[kMaxmc_track];   //[mc_track_]
   Float_t         mc_track_pt[kMaxmc_track];   //[mc_track_]
   Float_t         mc_track_eta[kMaxmc_track];   //[mc_track_]
   Float_t         mc_track_phi[kMaxmc_track];   //[mc_track_]
   Int_t           track_;
   UInt_t          track_fUniqueID[kMaxtrack];   //[track_]
   UInt_t          track_fBits[kMaxtrack];   //[track_]
   Float_t         track_pt[kMaxtrack];   //[track_]
   Float_t         track_eta[kMaxtrack];   //[track_]
   Float_t         track_phi[kMaxtrack];   //[track_]
   Float_t         track_dcaXY[kMaxtrack];   //[track_]
   Float_t         track_dcaXYZ[kMaxtrack];   //[track_]
   Bool_t          track_TOF_match[kMaxtrack];   //[track_]
   Bool_t          track_BEMC_match[kMaxtrack];   //[track_]
   Short_t         track_towerID[kMaxtrack];   //[track_]
   Float_t         track_towerEt[kMaxtrack];   //[track_]
   Short_t         track_nHitsFit[kMaxtrack];   //[track_]
   Short_t         track_nHitsPoss[kMaxtrack];   //[track_]
   Short_t         track_nHitsDedx[kMaxtrack];   //[track_]
   Bool_t          track_pass_cuts[kMaxtrack];   //[track_]
   Int_t           tower_;
   UInt_t          tower_fUniqueID[kMaxtower];   //[tower_]
   UInt_t          tower_fBits[kMaxtower];   //[tower_]
   Float_t         tower_Et[kMaxtower];   //[tower_]
   Float_t         tower_eta[kMaxtower];   //[tower_]
   Float_t         tower_phi[kMaxtower];   //[tower_]
   Float_t         tower_Et_hadroncorr[kMaxtower];   //[tower_]
   Short_t         tower_towerID[kMaxtower];   //[tower_]

   // List of branches
   TBranch        *b_mc_Cjet_;   //!
   TBranch        *b_mc_Cjet_fUniqueID;   //!
   TBranch        *b_mc_Cjet_fBits;   //!
   TBranch        *b_mc_Cjet_pt;   //!
   TBranch        *b_mc_Cjet_eta;   //!
   TBranch        *b_mc_Cjet_phi;   //!
   TBranch        *b_mc_Cjet_area;   //!
   TBranch        *b_mc_Cjet_index_track;   //!
   TBranch        *b_mc_Cjet_index_tower;   //!
   TBranch        *b_mc_Cjet_rho;   //!
   TBranch        *b_mc_Cjet_rho_sigma;   //!
   TBranch        *b_Cjet_;   //!
   TBranch        *b_Cjet_fUniqueID;   //!
   TBranch        *b_Cjet_fBits;   //!
   TBranch        *b_Cjet_pt;   //!
   TBranch        *b_Cjet_eta;   //!
   TBranch        *b_Cjet_phi;   //!
   TBranch        *b_Cjet_area;   //!
   TBranch        *b_Cjet_index_track;   //!
   TBranch        *b_Cjet_index_tower;   //!
   TBranch        *b_Cjet_rho;   //!
   TBranch        *b_Cjet_rho_sigma;   //!
   TBranch        *b_Fjet_;   //!
   TBranch        *b_Fjet_fUniqueID;   //!
   TBranch        *b_Fjet_fBits;   //!
   TBranch        *b_Fjet_pt;   //!
   TBranch        *b_Fjet_eta;   //!
   TBranch        *b_Fjet_phi;   //!
   TBranch        *b_Fjet_area;   //!
   TBranch        *b_Fjet_index_track;   //!
   TBranch        *b_Fjet_index_tower;   //!
   TBranch        *b_Fjet_rho;   //!
   TBranch        *b_Fjet_rho_sigma;   //!
   TBranch        *b_mu_event_fUniqueID;   //!
   TBranch        *b_mu_event_fBits;   //!
   TBranch        *b_mu_event_runId;   //!
   TBranch        *b_mu_event_eventId;   //!
   TBranch        *b_mu_event_ZDCx;   //!
   TBranch        *b_mu_event_BBCx;   //!
   TBranch        *b_mu_event_vx;   //!
   TBranch        *b_mu_event_vy;   //!
   TBranch        *b_mu_event_vz;   //!
   TBranch        *b_mu_event_BBC_Ein;   //!
   TBranch        *b_mu_event_BBC_Eout;   //!
   TBranch        *b_mu_event_BBC_Win;   //!
   TBranch        *b_mu_event_BBC_Wout;   //!
   TBranch        *b_mu_event_vzVpd;   //!
   TBranch        *b_mu_event_ranking;   //!
   TBranch        *b_mu_event_ZdcSumAdcEast;   //!
   TBranch        *b_mu_event_ZdcSumAdcWest;   //!
   TBranch        *b_mu_event_ZdcSmdEastHorizontal;   //!
   TBranch        *b_mu_event_ZdcSmdEastVertical;   //!
   TBranch        *b_mu_event_ZdcSmdWestHorizontal;   //!
   TBranch        *b_mu_event_ZdcSmdWestVertical;   //!
   TBranch        *b_mu_event_EastBBC;   //!
   TBranch        *b_mu_event_WestBBC;   //!
   TBranch        *b_mu_event_refMult;   //!
   TBranch        *b_mu_event_nGlobalTracks;   //!
   TBranch        *b_pthat_min;   //!
   TBranch        *b_pthat_max;   //!
   TBranch        *b_pthat_bin;   //!
   TBranch        *b_mc_track_;   //!
   TBranch        *b_mc_track_fUniqueID;   //!
   TBranch        *b_mc_track_fBits;   //!
   TBranch        *b_mc_track_geantId;   //!
   TBranch        *b_mc_track_id;   //!
   TBranch        *b_mc_track_pt;   //!
   TBranch        *b_mc_track_eta;   //!
   TBranch        *b_mc_track_phi;   //!
   TBranch        *b_track_;   //!
   TBranch        *b_track_fUniqueID;   //!
   TBranch        *b_track_fBits;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_dcaXY;   //!
   TBranch        *b_track_dcaXYZ;   //!
   TBranch        *b_track_TOF_match;   //!
   TBranch        *b_track_BEMC_match;   //!
   TBranch        *b_track_towerID;   //!
   TBranch        *b_track_towerEt;   //!
   TBranch        *b_track_nHitsFit;   //!
   TBranch        *b_track_nHitsPoss;   //!
   TBranch        *b_track_nHitsDedx;   //!
   TBranch        *b_track_pass_cuts;   //!
   TBranch        *b_tower_;   //!
   TBranch        *b_tower_fUniqueID;   //!
   TBranch        *b_tower_fBits;   //!
   TBranch        *b_tower_Et;   //!
   TBranch        *b_tower_eta;   //!
   TBranch        *b_tower_phi;   //!
   TBranch        *b_tower_Et_hadroncorr;   //!
   TBranch        *b_tower_towerID;   //!

   Reader(TTree *tree=0);
   virtual ~Reader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Reader_cxx
Reader::Reader(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("hadd_piKp_track_emb.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("hadd_piKp_track_emb.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

Reader::~Reader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Reader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Reader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Reader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mc_Cjet", &mc_Cjet_, &b_mc_Cjet_);
   fChain->SetBranchAddress("mc_Cjet.fUniqueID", mc_Cjet_fUniqueID, &b_mc_Cjet_fUniqueID);
   fChain->SetBranchAddress("mc_Cjet.fBits", mc_Cjet_fBits, &b_mc_Cjet_fBits);
   fChain->SetBranchAddress("mc_Cjet.pt", mc_Cjet_pt, &b_mc_Cjet_pt);
   fChain->SetBranchAddress("mc_Cjet.eta", mc_Cjet_eta, &b_mc_Cjet_eta);
   fChain->SetBranchAddress("mc_Cjet.phi", mc_Cjet_phi, &b_mc_Cjet_phi);
   fChain->SetBranchAddress("mc_Cjet.area", mc_Cjet_area, &b_mc_Cjet_area);
   fChain->SetBranchAddress("mc_Cjet.index_track", mc_Cjet_index_track, &b_mc_Cjet_index_track);
   fChain->SetBranchAddress("mc_Cjet.index_tower", mc_Cjet_index_tower, &b_mc_Cjet_index_tower);
   fChain->SetBranchAddress("mc_Cjet_rho", &mc_Cjet_rho, &b_mc_Cjet_rho);
   fChain->SetBranchAddress("mc_Cjet_rho_sigma", &mc_Cjet_rho_sigma, &b_mc_Cjet_rho_sigma);
   fChain->SetBranchAddress("Cjet", &Cjet_, &b_Cjet_);
   fChain->SetBranchAddress("Cjet.fUniqueID", Cjet_fUniqueID, &b_Cjet_fUniqueID);
   fChain->SetBranchAddress("Cjet.fBits", Cjet_fBits, &b_Cjet_fBits);
   fChain->SetBranchAddress("Cjet.pt", Cjet_pt, &b_Cjet_pt);
   fChain->SetBranchAddress("Cjet.eta", Cjet_eta, &b_Cjet_eta);
   fChain->SetBranchAddress("Cjet.phi", Cjet_phi, &b_Cjet_phi);
   fChain->SetBranchAddress("Cjet.area", Cjet_area, &b_Cjet_area);
   fChain->SetBranchAddress("Cjet.index_track", Cjet_index_track, &b_Cjet_index_track);
   fChain->SetBranchAddress("Cjet.index_tower", Cjet_index_tower, &b_Cjet_index_tower);
   fChain->SetBranchAddress("Cjet_rho", &Cjet_rho, &b_Cjet_rho);
   fChain->SetBranchAddress("Cjet_rho_sigma", &Cjet_rho_sigma, &b_Cjet_rho_sigma);
   fChain->SetBranchAddress("Fjet", &Fjet_, &b_Fjet_);
   fChain->SetBranchAddress("Fjet.fUniqueID", Fjet_fUniqueID, &b_Fjet_fUniqueID);
   fChain->SetBranchAddress("Fjet.fBits", Fjet_fBits, &b_Fjet_fBits);
   fChain->SetBranchAddress("Fjet.pt", Fjet_pt, &b_Fjet_pt);
   fChain->SetBranchAddress("Fjet.eta", Fjet_eta, &b_Fjet_eta);
   fChain->SetBranchAddress("Fjet.phi", Fjet_phi, &b_Fjet_phi);
   fChain->SetBranchAddress("Fjet.area", Fjet_area, &b_Fjet_area);
   fChain->SetBranchAddress("Fjet.index_track", Fjet_index_track, &b_Fjet_index_track);
   fChain->SetBranchAddress("Fjet.index_tower", Fjet_index_tower, &b_Fjet_index_tower);
   fChain->SetBranchAddress("Fjet_rho", &Fjet_rho, &b_Fjet_rho);
   fChain->SetBranchAddress("Fjet_rho_sigma", &Fjet_rho_sigma, &b_Fjet_rho_sigma);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_mu_event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_mu_event_fBits);
   fChain->SetBranchAddress("runId", &runId, &b_mu_event_runId);
   fChain->SetBranchAddress("eventId", &eventId, &b_mu_event_eventId);
   fChain->SetBranchAddress("ZDCx", &ZDCx, &b_mu_event_ZDCx);
   fChain->SetBranchAddress("BBCx", &BBCx, &b_mu_event_BBCx);
   fChain->SetBranchAddress("vx", &vx, &b_mu_event_vx);
   fChain->SetBranchAddress("vy", &vy, &b_mu_event_vy);
   fChain->SetBranchAddress("vz", &vz, &b_mu_event_vz);
   fChain->SetBranchAddress("BBC_Ein", &BBC_Ein, &b_mu_event_BBC_Ein);
   fChain->SetBranchAddress("BBC_Eout", &BBC_Eout, &b_mu_event_BBC_Eout);
   fChain->SetBranchAddress("BBC_Win", &BBC_Win, &b_mu_event_BBC_Win);
   fChain->SetBranchAddress("BBC_Wout", &BBC_Wout, &b_mu_event_BBC_Wout);
   fChain->SetBranchAddress("vzVpd", &vzVpd, &b_mu_event_vzVpd);
   fChain->SetBranchAddress("ranking", &ranking, &b_mu_event_ranking);
   fChain->SetBranchAddress("ZdcSumAdcEast", &ZdcSumAdcEast, &b_mu_event_ZdcSumAdcEast);
   fChain->SetBranchAddress("ZdcSumAdcWest", &ZdcSumAdcWest, &b_mu_event_ZdcSumAdcWest);
   fChain->SetBranchAddress("ZdcSmdEastHorizontal[8]", ZdcSmdEastHorizontal, &b_mu_event_ZdcSmdEastHorizontal);
   fChain->SetBranchAddress("ZdcSmdEastVertical[8]", ZdcSmdEastVertical, &b_mu_event_ZdcSmdEastVertical);
   fChain->SetBranchAddress("ZdcSmdWestHorizontal[8]", ZdcSmdWestHorizontal, &b_mu_event_ZdcSmdWestHorizontal);
   fChain->SetBranchAddress("ZdcSmdWestVertical[8]", ZdcSmdWestVertical, &b_mu_event_ZdcSmdWestVertical);
   fChain->SetBranchAddress("EastBBC[24]", EastBBC, &b_mu_event_EastBBC);
   fChain->SetBranchAddress("WestBBC[24]", WestBBC, &b_mu_event_WestBBC);
   fChain->SetBranchAddress("refMult", &refMult, &b_mu_event_refMult);
   fChain->SetBranchAddress("nGlobalTracks", &nGlobalTracks, &b_mu_event_nGlobalTracks);
   fChain->SetBranchAddress("pthat_min", &pthat_min, &b_pthat_min);
   fChain->SetBranchAddress("pthat_max", &pthat_max, &b_pthat_max);
   fChain->SetBranchAddress("pthat_bin", &pthat_bin, &b_pthat_bin);
   fChain->SetBranchAddress("mc_track", &mc_track_, &b_mc_track_);
   fChain->SetBranchAddress("mc_track.fUniqueID", mc_track_fUniqueID, &b_mc_track_fUniqueID);
   fChain->SetBranchAddress("mc_track.fBits", mc_track_fBits, &b_mc_track_fBits);
   fChain->SetBranchAddress("mc_track.geantId", mc_track_geantId, &b_mc_track_geantId);
   fChain->SetBranchAddress("mc_track.id", mc_track_id, &b_mc_track_id);
   fChain->SetBranchAddress("mc_track.pt", mc_track_pt, &b_mc_track_pt);
   fChain->SetBranchAddress("mc_track.eta", mc_track_eta, &b_mc_track_eta);
   fChain->SetBranchAddress("mc_track.phi", mc_track_phi, &b_mc_track_phi);
   fChain->SetBranchAddress("track", &track_, &b_track_);
   fChain->SetBranchAddress("track.fUniqueID", track_fUniqueID, &b_track_fUniqueID);
   fChain->SetBranchAddress("track.fBits", track_fBits, &b_track_fBits);
   fChain->SetBranchAddress("track.pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track.eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track.phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track.dcaXY", track_dcaXY, &b_track_dcaXY);
   fChain->SetBranchAddress("track.dcaXYZ", track_dcaXYZ, &b_track_dcaXYZ);
   fChain->SetBranchAddress("track.TOF_match", track_TOF_match, &b_track_TOF_match);
   fChain->SetBranchAddress("track.BEMC_match", track_BEMC_match, &b_track_BEMC_match);
   fChain->SetBranchAddress("track.towerID", track_towerID, &b_track_towerID);
   fChain->SetBranchAddress("track.towerEt", track_towerEt, &b_track_towerEt);
   fChain->SetBranchAddress("track.nHitsFit", track_nHitsFit, &b_track_nHitsFit);
   fChain->SetBranchAddress("track.nHitsPoss", track_nHitsPoss, &b_track_nHitsPoss);
   fChain->SetBranchAddress("track.nHitsDedx", track_nHitsDedx, &b_track_nHitsDedx);
   fChain->SetBranchAddress("track.pass_cuts", track_pass_cuts, &b_track_pass_cuts);
   fChain->SetBranchAddress("tower", &tower_, &b_tower_);
   fChain->SetBranchAddress("tower.fUniqueID", tower_fUniqueID, &b_tower_fUniqueID);
   fChain->SetBranchAddress("tower.fBits", tower_fBits, &b_tower_fBits);
   fChain->SetBranchAddress("tower.Et", tower_Et, &b_tower_Et);
   fChain->SetBranchAddress("tower.eta", tower_eta, &b_tower_eta);
   fChain->SetBranchAddress("tower.phi", tower_phi, &b_tower_phi);
   fChain->SetBranchAddress("tower.Et_hadroncorr", tower_Et_hadroncorr, &b_tower_Et_hadroncorr);
   fChain->SetBranchAddress("tower.towerID", tower_towerID, &b_tower_towerID);
   Notify();
}

Bool_t Reader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Reader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Reader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Reader_cxx
