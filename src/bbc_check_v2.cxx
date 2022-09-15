// tower_count.cxx
// David Stewart  2022.03.04 -> modified from V. Verkest code

#include "pAuFunctions.hh"
#include "pAu_params.h"
/* #include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!! */
#include <iostream>
#include <fstream>
#include <map>
#include "pAu_consts.h"


using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

// functions and cuts specified in pAuFunctions.hh
// count the trigger's available on a per-run basis
int main ( int argc, const char** argv ) {

    int    number_of_events { (argc>1) ? atoi(argv[1]) : 100                                                  };
    int    triggerId        { (argc>2) ? atoi(argv[2]) : 500004                                               };
    string inFile           { (argc>3) ? argv[3]       : "production_pAu200_2015_newTriggers/pAu_2015_0.root" };
    string outFile          { (argc>4) ? argv[4]       : "bbc_check_v2"                                       };
    string runPeriod        { (argc>5) ? argv[5]       : "ALL"                                                };// 1, 2, or 3

    outFile += ".root";
    cout << " outFile: " << outFile << endl;

    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;

    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    TChain* Chain = new TChain( "JetTree" );
    Chain->Add( inFile.c_str() );
    TStarJetPicoReader Reader;
    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/dummy_tower.list");
    towerCuts->SetMaxEtCut( 1.0 );
    /* cout << " a0 " << endl; */
    /* int numEvents = number_of_events;        // total events in HT: 152,007,032 */
    InitReader( Reader, Chain, number_of_events, 1. );

    ioBinVec bin_zdcx     {{ 0.,   0.,   28,  28000 }};
    ioBinVec bin_bbcEin   {{ 0.,   0.,   80., 80000 }};
    /* ioBinVec bin_gntracks {{ 0.,   0.,   300, 3000  }}; */
    ioBinVec bin_gntracks {{  0., 0., 26, 2600 }};

    ioBinVec bin_Nch      {{ -0.5, -0.5, 81, 81.5  }};
    ioBinVec bin_dca      {{   0.,   0., 30, 3.    }};

    ioBinVec bin3_pt   {{ 0., 0., 300, 30. }};
    ioBinVec bin3_glob {{ 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 
                          1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600 }};
    /* ioBinVec bin3_zdcx {{ 0., 4000., 10000., 16000., 22000., 50000 }}; */
    ioBinVec bin3_zdcx {{ 0., 0., 25, 25000. }};


    int nbins_bbcEout  = 25;  double lobin_bbcEout  = 0.; double hibin_bbcEout  = 49000;
    int nbins_bbcWin   = 25;  double lobin_bbcWin   = 0.; double hibin_bbcWin   = 45000;
    int nbins_bbcWout  = 25;  double lobin_bbcWout  = 0.; double hibin_bbcWout  = 12000;
    int nbins_refmult  = 100; double lobin_refmult  = 0.; double hibin_refmult  = nbins_refmult;
    int nbins_grefmult = 100; double lobin_grefmult = 0.; double hibin_grefmult = nbins_grefmult;
    int nbins_gntracks = 300; double lobin_gntracks = 0.; double hibin_gntracks = 3000;

    TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

    TH1D     hg_zdcx          { "zdcx",           ";zdcx;N_{events}",                bin_zdcx,       bin_zdcx        };
    TProfile prNch_zdcx       { "pr_Nch_zdcx",    ";zdcx;<Nch>",                     bin_zdcx,       bin_zdcx        };

    TH3D  hg3_tracks { "hg_tracks", "pT;glob;zdcx", bin3_glob, bin3_glob, bin3_zdcx, bin3_zdcx, bin3_pt, bin3_pt };
    TH2D  hg2_events { "hg_events", "n_events;glob;zdcx", bin3_glob, bin3_glob, bin3_zdcx, bin3_zdcx };

    TH1D     hg_bbcEin        { "bbcEin",         ";bbcIn;N_{events}",               bin_bbcEin,     bin_bbcEin      };
    TH1D     hg_bbcEout       { "bbcEout",        ";bbcOut;N_{events}",              nbins_bbcEout,  lobin_bbcEout,  hibin_bbcEout   };
    TH1D     hg_bbcWin        { "bbcWin",         ";bbcIn;N_{events}",               nbins_bbcWin,   lobin_bbcWin,   hibin_bbcWin    };
    TH1D     hg_bbcWout       { "bbcWout",        ";bbcOut;N_{events}",              nbins_bbcWout,  lobin_bbcWout,  hibin_bbcWout   };
    TH1D     hg_refmult       { "refmult",        ";refmult;N_{events}",             nbins_refmult,  lobin_refmult,  hibin_refmult   };
    TH1D     hg_grefmult      { "grefmult",       ";grefmult;N_{events}",            nbins_grefmult, lobin_grefmult, hibin_grefmult  };

    TH1D     hg_gntracks      { "gntracks",       ";N Global Tracks;N_{events}",     bin_gntracks,   bin_gntracks    };
    TProfile prNch_gntracks   { "prNch_gntracks", ";N Global Tracks;<Nch>",          bin_gntracks,   bin_gntracks    };

    TH2D hg_bbcEin_zdcx   { "bbcEin_zdcx",   "N_{events};bbcIn;zdcx",           bin_bbcEin,     bin_bbcEin,     bin_zdcx,       bin_zdcx  };
    TH2D hg_bbcEout_zdcx  { "bbcEout_zdcx",  "N_{events};bbcOut;zdcx",          nbins_bbcEout,  lobin_bbcEout,  hibin_bbcEout,  bin_zdcx, bin_zdcx };
    TH2D hg_bbcWin_zdcx   { "bbcWin_zdcx",   "N_{events};bbcIn;zdcx",           nbins_bbcWin,   lobin_bbcWin,   hibin_bbcWin,   bin_zdcx, bin_zdcx };
    TH2D hg_bbcWout_zdcx  { "bbcWout_zdcx",  "N_{events};bbcOut;zdcx",          nbins_bbcWout,  lobin_bbcWout,  hibin_bbcWout,  bin_zdcx, bin_zdcx };
    TH2D hg_refmult_zdcx  { "refmult_zdcx",  "N_{events};refmult;zdcx",         nbins_refmult,  lobin_refmult,  hibin_refmult,  bin_zdcx, bin_zdcx };
    TH2D hg_grefmult_zdcx { "grefmult_zdcx", "N_{events};grefmult;zdcx",        nbins_grefmult, lobin_grefmult, hibin_grefmult, bin_zdcx, bin_zdcx };
    TH2D hg_gntracks_zdcx { "gntracks_zdcx", "N_{events};N Global Tracks;zdcx", nbins_gntracks, lobin_gntracks, hibin_gntracks, bin_zdcx, bin_zdcx };

    TH2D hg_Nch_zdcx     { "Nch_zdcx",     "N_{events};Nch;zdcx",            bin_Nch, bin_Nch, bin_zdcx,     bin_zdcx     };
    TH2D hg_Nch_gntracks { "Nch_gntracks", "N_{events};Nch;N Global Tracks", bin_Nch, bin_Nch, bin_gntracks, bin_gntracks };

    TH2D dca_zdcx     { "dca_zdcx",     ";zdcX;DCA_{3D}",            bin_zdcx,     bin_zdcx,     bin_dca, bin_dca };
    TH2D dca_gntracks { "dca_gntracks", ";N global tracks;DCA_{3D}", bin_gntracks, bin_gntracks, bin_dca, bin_dca };

    TH2D dca_zdcx_TOF     { "dca_zdcx_TOF",     "TOF matched;zdcX;DCA_{3D}",            bin_zdcx,     bin_zdcx,     bin_dca, bin_dca };
    TH2D dca_gntracks_TOF { "dca_gntracks_TOF", "TOF matched;N global tracks;DCA_{3D}", bin_gntracks, bin_gntracks, bin_dca, bin_dca };

    

    bool is_jet_trig = triggerId > 500004;

    while ( Reader.NextEvent() ) {
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();

        if (!GOOD_RUNS_v0.has(header->GetRunId())
        ||  !header->HasTriggerId(triggerId)) continue;
        if (TMath::Abs(header->GetPrimaryVertexZ())>10.) continue;

        // get the trigger tower

        // get the lead tower
        double trig_Et{0};
        double trig_phi{0};
        if (is_jet_trig) {
            TList *l_tows = Reader.GetListOfSelectedTowers();//GetTowers();
            TIter next_tow(l_tows);
            while (TStarJetPicoTower* tower = (TStarJetPicoTower *) next_tow()) {
                double Et = tower->GetEt();
                if (Et > trig_Et && Et > 4.) {
                    trig_Et  = Et;
                    trig_phi = tower->GetPhi();
                }
            }
            if (trig_Et < 4.) continue; // failed to find a 4 GeV trigger
        }

        float Nch=0;
        TList *p_prim_trks = Reader.GetListOfSelectedTracks();
        TIter nextptrk(p_prim_trks);
        double zdcx  = header->GetZdcCoincidenceRate();
        double nglob = header->GetNGlobalTracks();
        while (TStarJetPicoPrimaryTrack* p_track = (TStarJetPicoPrimaryTrack *) nextptrk()) {
            if (is_jet_trig && !is_trans(p_track->GetPhi(),trig_phi)) continue;
            Nch += 1.;
            if (p_track->GetDCA() <= 1.) hg3_tracks.Fill(nglob, zdcx, p_track->GetPt());
            dca_zdcx.Fill(header->GetZdcCoincidenceRate(), p_track->GetDCA());
            dca_gntracks.Fill(header->GetNGlobalTracks(), p_track->GetDCA());

            cout << " " << p_track->GetTofMatchFlag() << " " << p_track->GetBemcMatchFlag() 
                 << " " << p_track->GetTofTime() << endl;

            if (p_track->GetTofMatchFlag()) {
                dca_zdcx_TOF    .Fill(header->GetZdcCoincidenceRate(), p_track->GetDCA());
                dca_gntracks_TOF.Fill(header->GetNGlobalTracks(),      p_track->GetDCA());
            }


        }
        hg2_events.Fill(nglob, zdcx);

        
    

        hg_zdcx       .Fill(header->GetZdcCoincidenceRate());
        prNch_zdcx    .Fill(header->GetZdcCoincidenceRate(), Nch);

        hg_bbcEin     .Fill(header->GetBbcAdcSumEastInner());
        hg_bbcEout    .Fill(header->GetBbcAdcSumEastOuter());
        hg_bbcWin     .Fill(header->GetBbcAdcSumWestInner());
        hg_bbcWout    .Fill(header->GetBbcAdcSumWestOuter());

        hg_refmult    .Fill(header->GetReferenceMultiplicity());
        hg_grefmult   .Fill(header->GetGReferenceMultiplicity());
        hg_gntracks   .Fill(header->GetNGlobalTracks());
        prNch_gntracks.Fill(header->GetNGlobalTracks(),Nch);

        hg_Nch_gntracks .Fill(Nch, header->GetNGlobalTracks());
        hg_Nch_zdcx     .Fill(Nch, header->GetZdcCoincidenceRate());

        hg_bbcEin_zdcx  .Fill(header->GetBbcAdcSumEastInner(), header->GetZdcCoincidenceRate());
        hg_bbcEout_zdcx .Fill(header->GetBbcAdcSumEastOuter(), header->GetZdcCoincidenceRate());
        hg_bbcWin_zdcx  .Fill(header->GetBbcAdcSumWestInner(), header->GetZdcCoincidenceRate());
        hg_bbcWout_zdcx .Fill(header->GetBbcAdcSumWestOuter(), header->GetZdcCoincidenceRate());

        hg_refmult_zdcx .Fill(header->GetReferenceMultiplicity(), header->GetZdcCoincidenceRate());
        hg_grefmult_zdcx.Fill(header->GetGReferenceMultiplicity(), header->GetZdcCoincidenceRate());
        hg_gntracks_zdcx.Fill(header->GetNGlobalTracks(), header->GetZdcCoincidenceRate());
    }

    double scale =  is_jet_trig ? 1/(4*M_PI/3.) : 1/(4*M_PI) ;
    prNch_zdcx       .Scale(scale);
    prNch_gntracks   .Scale(scale);
    hg_Nch_zdcx      .Scale(scale);
    hg_Nch_gntracks  .Scale(scale);

    hg_zdcx          .Write();
    prNch_zdcx       .Write();

    hg_bbcEin        .Write();
    hg_bbcEout       .Write();
    hg_bbcWin        .Write();
    hg_bbcWout       .Write();
    hg_refmult       .Write();
    hg_grefmult      .Write();

    hg_gntracks      .Write();
    prNch_gntracks   .Write();

    hg_bbcEin_zdcx   .Write();
    hg_bbcEout_zdcx  .Write();
    hg_bbcWin_zdcx   .Write();
    hg_bbcWout_zdcx  .Write();
    hg_refmult_zdcx  .Write();
    hg_grefmult_zdcx .Write();
    hg_gntracks_zdcx .Write();

    hg_Nch_zdcx      .Write();
    hg_Nch_gntracks  .Write();

    dca_zdcx         .Write();
    dca_gntracks     .Write();

    dca_zdcx_TOF     .Write();
    dca_gntracks_TOF .Write();

    hg3_tracks .Write();
    hg2_events .Write();

    pAuFile->Save();
    pAuFile->Close();

    cout << "finished. Recommed: " << endl;
    cout << " scl wsu:/tier2/home/groups/rhi/vverkest/pAu_analysis/" << outFile << " ." << endl;

    return 0;
}
