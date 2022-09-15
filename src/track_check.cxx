// tower_count.cxx
// David Stewart  2022.03.04 -> modified from V. Verkest code

#include "pAuFunctions.hh"
#include "pAu_params.h"
/* #include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!! */
#include <iostream>
#include <fstream>
#include <map>
#include "pAu_consts.h"
#include "THnSparse.h"


using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

// functions and cuts specified in pAuFunctions.hh
// count the trigger's available on a per-run basis
int main ( int argc, const char** argv ) {

    int    number_of_events { (argc>1) ? atoi(argv[1]) : 100                                                  };
    int    triggerId        { (argc>2) ? atoi(argv[2]) : 500206                                               };
    string inFile           { (argc>3) ? argv[3]       : "production_pAu200_2015_newTriggers/pAu_2015_0.root" };
    string outFile          { (argc>4) ? argv[4]       : "track_check"                                         };
    string runPeriod        { (argc>5) ? argv[5]       : "ALL"                                                };// 1, 2, or 3

    outFile += ".root";
    cout << " outFile: " << outFile << endl;
    TFile fout { outFile.c_str(),"recreate" };

    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;

    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    TChain* Chain = new TChain( "JetTree" );
    Chain->Add( inFile.c_str() );
    TStarJetPicoReader Reader;
    InitReader( Reader, Chain, number_of_events, 1. );

    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/dummy_tower.list");
    /* towerCuts->SetMaxEtCut( 1.0 ); */

    TStarJetPicoTrackCuts* trackCuts = Reader.GetTrackCuts();
    trackCuts->SetDCACut(6.);
    trackCuts->SetDCACut(3.);
    trackCuts->SetMinNFitPointsCut(1.);
    trackCuts->SetFitOverMaxPointsCut(0.);
    trackCuts->SetMaxPtCut(30.);


    /* cout << " a0 " << endl; */
    /* int numEvents = number_of_events;        // total events in HT: 152,007,032 */


    // Wanted:
    // 1. TH3D:
    //  0: ZDCx
    //  1: nGlobal
    //  2: BBCES
    // 2. THnSparse : overall, transverse
    //  0: ZDCx
    //  1: nGlobal
    //  2: BBCES
    //  3: track pT
    //  4: nhitsposs
    //  5: nhitsact
    //  6: dca
    //  7: pass_cuts

    ioBinVec bin_zdcx      {{ 4000., 4000., 8, 28000., 33000. }};
    ioBinVec bin_nGlobal   {{  300., 300., 10, 2500. }};
  // bins_EAbbc3  :: const ioBinVec bins_EAbbc3 {{     0.0000,  8315.4606, 24292.8207, 64000.0000 }};
    ioBinVec bin_trackpT   {{ 0.2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1, 1.2, 1.3, 1.4, 1.6, 1.8, 2., 5., 10 }};
    ioBinVec bin_nHitsPoss {{  0., 0., 46, 46. }};
    ioBinVec bin_nHitsFit  {{  0., 0., 46, 46. }};
    ioBinVec bin_dca       {{  0., 0.2, 0.4, 0.7, 1., 1.5, 2., 3. }};
    ioBinVec bin_passcuts  {{ -0.5, 0.5, 1.5 }}; // cut for nhitsfit >= 15 and ratio >=52
    ioBinVec bin_hastrig   {{ -0.5, 0.5, 1.5 }}; // cut for nhitsfit >= 15 and ratio >=52

    int nbins[8];
    nbins[0] = bin_zdcx;
    nbins[1] = bin_nGlobal;
    nbins[2] = bins_EAbbc3;
    nbins[3] = bin_trackpT;
    nbins[4] = bin_nHitsPoss;
    nbins[5] = bin_nHitsFit;
    nbins[6] = bin_dca;
    nbins[7] = bin_passcuts;

    THnSparseD* nevents = new THnSparseD("nevents","nevents;ZDCx;nGlobal;EA_{BBC-East};has 8GeV trig;", 3, nbins, NULL, NULL);
    nevents->SetBinEdges(0,bin_zdcx);
    nevents->SetBinEdges(1,bin_nGlobal);
    nevents->SetBinEdges(2,bins_EAbbc3);
    nevents->SetBinEdges(3,bin_hastrig);

    double hopper[8];

    THnSparseD* trans_track = new THnSparseD("trans_tracks", "transverse tracks;ZDCx;nGlobal;EA_{BBC-East};p_{T};"
            "nHitsPoss;nHitsFit;DCA;pass_cuts:(0,1)=(N,Y);",
        8, nbins, NULL, NULL);
    trans_track->SetBinEdges(0,bin_zdcx);
    trans_track->SetBinEdges(1,bin_nGlobal);
    trans_track->SetBinEdges(2,bins_EAbbc3);
    trans_track->SetBinEdges(3,bin_trackpT);
    trans_track->SetBinEdges(4,bin_nHitsPoss);
    trans_track->SetBinEdges(5,bin_nHitsFit);
    trans_track->SetBinEdges(6,bin_dca);
    trans_track->SetBinEdges(7,bin_passcuts);

    THnSparseD* all_track = new THnSparseD("all_tracks", "transverse tracks;ZDCx;nGlobal;EA_{BBC-East};p_{T};"
            "nHitsPoss;nHitsFit;DCA;pass_cuts:(0,1)=(N,Y);",
        8, nbins, NULL, NULL);
    all_track->SetBinEdges(0,bin_zdcx);
    all_track->SetBinEdges(1,bin_nGlobal);
    all_track->SetBinEdges(2,bins_EAbbc3);
    all_track->SetBinEdges(3,bin_trackpT);
    all_track->SetBinEdges(4,bin_nHitsPoss);
    all_track->SetBinEdges(5,bin_nHitsFit);
    all_track->SetBinEdges(6,bin_dca);
    all_track->SetBinEdges(7,bin_passcuts);

    std::set<int> which_trigs;

    while ( Reader.NextEvent() ) {
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();

        if (!GOOD_RUNS_v0.has(header->GetRunId())
        ||  !header->HasTriggerId(triggerId)) continue;
        if (TMath::Abs(header->GetPrimaryVertexZ())>10.) continue;

        // get the lead tower
        double trig_Et{0};
        double trig_phi{0};
        if (triggerId=500206) {
            TList *l_tows = Reader.GetListOfSelectedTowers();//GetTowers();
            TIter next_tow(l_tows);
            while (TStarJetPicoTower* tower = (TStarJetPicoTower *) next_tow()) {
                double Et = tower->GetEt();
                if (Et > trig_Et && Et > 4.) {
                    trig_Et  = Et;
                    trig_phi = tower->GetPhi();
                }
            }
            /* if (trig_Et < 8.) continue; // failed to find a 8 GeV trigger */
        }
        bool has_trig = trig_Et > 8.;

        float Nch=0;
        TList *p_prim_trks = Reader.GetListOfSelectedTracks();
        TIter nextptrk(p_prim_trks);
        double zdcx  = header->GetZdcCoincidenceRate();
        double nglob = header->GetNGlobalTracks();
        double bbc   = header->GetBbcAdcSumEastInner();


        hopper[0] = zdcx;
        hopper[1] = nglob;
        hopper[2] = bbc;
        hopper[3] = has_trig ? 1. : 0.;

        nevents->Fill(hopper);
        
        while (TStarJetPicoPrimaryTrack* p_track = (TStarJetPicoPrimaryTrack *) nextptrk()) {
            hopper[3] = p_track->GetPt();
            hopper[4] = p_track->GetNOfPossHits();
            hopper[5] = p_track->GetNOfFittedHits();
            hopper[6] = p_track->GetDCA();
            hopper[7] = (hopper[5] >= 14 && hopper[5]/hopper[4] >= 0.52) ? 1. : 0.;

            all_track->Fill(hopper);
            if (has_trig && is_trans(p_track->GetPhi(),trig_phi)) trans_track->Fill(hopper);
        }
    }
    nevents->Write();
    all_track->Write();
    trans_track->Write();
    fout.Save();

    cout << "finished. Recommed: " << endl;
    cout << " scl wsu:/tier2/home/groups/rhi/vverkest/pAu_analysis/" << outFile << " ." << endl;

    return 0;
}
