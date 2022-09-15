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

TH1D makeTH1D(const char* tag_meas, const char* tag_against, ioBinVec& bins) {
    return { Form("%s_%s",tag_meas, tag_against), Form("%s;%s;N_{events}",tag_against, tag_meas), bins, bins };
};

static ioBinVec bin_nHitsPoss {{  0., 0., 46, 46. }};
static ioBinVec bin_nHitsFit  {{  0., 0., 46, 46. }};
static ioBinVec bin_nch       {{ -0.5, -0.5, 51, 50.0 }};

static ioBinVec bin_zdcx      {{ 1000.0000, 6017.6527, 8225.4080, 9393.6762,
        10569.276, 11631.295, 12812.422, 14551.313, 17155.391, 19668.749,
        28000.000}};
static ioBinVec bin_nGlobal   {{  10.000000, 515.73632, 656.28995, 757.77145,
        847.70074, 935.17289, 1027.9222, 1132.6310, 1258.1754, 1428.6614,
        2990.0000 }};

// bins_EAbbc3  :: const ioBinVec bins_EAbbc3 {{     0.0000,  8315.4606, 24292.8207, 64000.0000 }};
static ioBinVec bin_dca       {{  0., 0.2, 0.4, 0.7, 1., 1.5, 2., 3. }};
static ioBinVec run_period    {{  16124017 , 16133086 , 16148019 , 16159024 }};

struct TH1D_TProf_Set {
    double &bbc, &zdcx, &nglob, &runid;
    TH1D     hg_all, hg_run0, hg_run1, hg_run2, hg_EAhi, hg_EAlo, hg_ZDCxhi, hg_ZDCxlo,
             hg_nGlobhi, hg_nGloblo;
    TProfile pr_E, pr_S; // 1,2,3 RunPeriod, 4,5 EA_BBC HiLow, 6,7 ZDCx HiLo, 8,9 nGlobal HiLo
    array<TH1D*,10> arr_hg;
    TH1D_TProf_Set(const char* tag, ioBinVec& bins, double& _bbc, double& _zdcx, double& _nglob, double& _runid) 
        : bbc{_bbc}
        , zdcx{_zdcx}
        , nglob{_nglob}
        , runid{_runid}
        , hg_all     { makeTH1D(tag, "all",     bins) }
        , hg_run0    { makeTH1D(tag, "run0",    bins) }
        , hg_run1    { makeTH1D(tag, "run1",    bins) }
        , hg_run2    { makeTH1D(tag, "run2",    bins) }
        , hg_EAhi    { makeTH1D(tag, "EAhi",    bins) }
        , hg_EAlo    { makeTH1D(tag, "EAlo",    bins) }
        , hg_ZDCxhi  { makeTH1D(tag, "ZDCxhi",  bins) }
        , hg_ZDCxlo  { makeTH1D(tag, "ZDCxlo",  bins) }
        , hg_nGlobhi { makeTH1D(tag, "nGlobhi", bins) }
        , hg_nGloblo { makeTH1D(tag, "nGloblo", bins) }
        , pr_E{ Form("pr_E_%s",tag), Form(
            ";all[0],run(0,1,2)[1,2,3] EAlo/hi[4,5] ZDCxlo/hi[6,7] nGloblo/hi[8,9];mean(%s),mean(%s)",tag),
                10, -0.5, 9.5, "E" }
        , pr_S{ Form("pr_S_%s",tag), Form(
            ";all[0],run(0,1,2)[1,2,3] EAlo/hi[4,5] ZDCxlo/hi[6,7] nGloblo/hi[8,9];mean(%s),mean(%s)",tag),
                10, -0.5, 9.5, "S" }
        , arr_hg { &hg_all, &hg_run0, &hg_run1, &hg_run2, &hg_EAlo, &hg_EAhi, &hg_ZDCxlo, &hg_ZDCxhi, &hg_nGloblo, &hg_nGlobhi }
    {}
    void fill (const double val) {
        bool is_run0 { runid >= run_period[0] && runid <= run_period[1] };
        bool is_run1 { runid >= run_period[1] && runid <= run_period[2] };
        bool is_run2 { runid >= run_period[2] && runid <= run_period[3] };
        bool isEAlo  { bbc < bins_EAbbc10[3] };
        bool isEAhi  { bbc > bins_EAbbc10[7] };
        bool isZDCxlo { zdcx < bin_zdcx[3] };
        bool isZDCxhi { zdcx > bin_zdcx[7] };
        bool isGloblo { nglob < bin_nGlobal[3] };
        bool isGlobhi { nglob > bin_nGlobal[7] };

        if (true)     hg_all.Fill(val);
        if (is_run0)  hg_run0.Fill(val);
        if (is_run1)  hg_run1.Fill(val);
        if (is_run2)  hg_run2.Fill(val);
        if (isEAlo)   hg_EAlo.Fill(val);
        if (isEAhi)   hg_EAhi.Fill(val);
        if (isZDCxlo) hg_ZDCxlo.Fill(val);
        if (isZDCxhi) hg_ZDCxhi.Fill(val);
        if (isGloblo) hg_nGloblo.Fill(val);
        if (isGlobhi) hg_nGlobhi.Fill(val);

        if (true)     pr_E.Fill(0., val); pr_S.Fill(0., val); 
        if (is_run0)  { pr_E.Fill(1., val); pr_S.Fill(1., val); }
        if (is_run1)  { pr_E.Fill(2., val); pr_S.Fill(2., val); }
        if (is_run2)  { pr_E.Fill(3., val); pr_S.Fill(3., val); } 
        if (isEAlo)   { pr_E.Fill(4., val); pr_S.Fill(4., val); } 
        if (isEAhi)   { pr_E.Fill(5., val); pr_S.Fill(5., val); } 
        if (isZDCxlo) { pr_E.Fill(6., val); pr_S.Fill(6., val); } 
        if (isZDCxhi) { pr_E.Fill(7., val); pr_S.Fill(7., val); } 
        if (isGloblo) { pr_E.Fill(8., val); pr_S.Fill(8., val); } 
        if (isGlobhi) { pr_E.Fill(9., val); pr_S.Fill(9., val); } 

    }
    void write() { 
        for (auto& hg : arr_hg) hg->Write();
        hg_EAlo.Write();
        pr_E.Write();
        pr_S.Write();
    }
};

// functions and cuts specified in pAuFunctions.hh
// count the trigger's available on a per-run basis
int main ( int argc, const char** argv ) {
    int    number_of_events { (argc>1) ? atoi(argv[1]) : 100                                                  };
    int    triggerId        { (argc>2) ? atoi(argv[2]) : 500206                                               };
    string inFile           { (argc>3) ? argv[3]       : "production_pAu200_2015_newTriggers/pAu_2015_0.root" };
    string outFile          { (argc>4) ? argv[4]       : "runQA"                                         };

    outFile += ".root";
    cout << " outFile: " << outFile << endl;
    TFile fout { outFile.c_str(),"recreate" };

    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;

    ioBinVec bin_nHitsPoss {{ 0.,   0.,   46,  46.  }};
    ioBinVec bin_nHitsFit  {{ 0.,   0.,   46,  46.  }};
    ioBinVec bin_EtpT      {{ 0.,   0.,   300, 30.  }};
    ioBinVec bin_nch       {{ -0.5, -0.5, 101, 100. }};
    ioBinVec bin_nTow      {{ -0.5, -0.5, 801, 800. }};

  // bins_EAbbc3  :: const ioBinVec bins_EAbbc3 {{     0.0000,  8315.4606, 24292.8207, 64000.0000 }};
    ioBinVec bin_dca       {{ 0., 0., 30, 3. }};
    ioBinVec bin_nHitsRat  {{ 0., 0., 100., 1. }};

    double bbc;
    double zdcx;
    double nglob;
    double runid;

    TH1D_TProf_Set data_Et        { "Et",        bin_EtpT     , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_nETp2     { "nETp2",     bin_nTow     , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_nET4GeV   { "nET4GeV",   bin_nTow     , bbc, zdcx, nglob, runid };

    TH1D_TProf_Set data_nhitsfit  { "nHitsFit",  bin_nHitsFit , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_nhitsposs { "nHitsPoss", bin_nHitsFit , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_nhitsrat  { "nHitsRat",  bin_nHitsRat , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_dca       { "DCA",       bin_dca      , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_pt        { "pt",        bin_EtpT     , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_nch       { "nch",       bin_nch      , bbc, zdcx, nglob, runid };

    TH1D_TProf_Set data_nhitsfitT  { "nHitsFit_trans",  bin_nHitsFit, bbc, zdcx, nglob, runid  };
    TH1D_TProf_Set data_nhitspossT { "nHitsPoss_trans", bin_nHitsFit, bbc, zdcx, nglob, runid  };
    TH1D_TProf_Set data_nhitsratT  { "nHitsRat_trans",  bin_nHitsRat , bbc, zdcx, nglob, runid };
    TH1D_TProf_Set data_dcaT       { "DCA_trans",       bin_dca     , bbc, zdcx, nglob, runid  };
    TH1D_TProf_Set data_ptT        { "pt_trans",        bin_EtpT    , bbc, zdcx, nglob, runid  };
    TH1D_TProf_Set data_nchT       { "nch_trans",       bin_nch     , bbc, zdcx, nglob, runid  };

    ioBinVec _bin_zdcx   {{ 0., 0., 28, 28000 }};
    ioBinVec _bin_bbcEin {{ 0., 0., 80., 80000 }};
    ioBinVec bin_nGlob  {{ 0., 0., 300, 3000. }};

    TH1D hg_zdcx   { "zdcx",   ";zdcx;N_{events}",            _bin_zdcx,   _bin_zdcx   };
    TH1D hg_bbcEin { "bbcEin", ";bbcIn;N_{events}",           _bin_bbcEin, _bin_bbcEin };
    TH1D hg_nGlob  { "nGlob",  ";N Global Tracks;N_{events}", bin_nGlob, bin_nGlob  };
    TH1D hg_runPeriod { "nEvent", ";runId;N_{events}",        run_period, run_period };


    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    TChain* Chain = new TChain( "JetTree" );
    Chain->Add( inFile.c_str() );
    TStarJetPicoReader Reader;
    InitReader( Reader, Chain, number_of_events, 1. );

    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/BAD_TOWERS_v0.list");
    /* towerCuts->SetMaxEtCut( 1.0 ); */

    TStarJetPicoTrackCuts* trackCuts = Reader.GetTrackCuts();
    trackCuts->SetDCACut(3.);
    trackCuts->SetMinNFitPointsCut(1.);
    trackCuts->SetFitOverMaxPointsCut(0.);
    trackCuts->SetMaxPtCut(30.);

    /* std::set<int> which_trigs; */

    while ( Reader.NextEvent() ) {
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();

        if (!GOOD_RUNS_v0.has(header->GetRunId())
        ||  !header->HasTriggerId(triggerId)) continue;
        if (TMath::Abs(header->GetPrimaryVertexZ())>10.) continue;

        runid = header->GetRunId();
        zdcx  = header->GetZdcCoincidenceRate();
        nglob = header->GetNGlobalTracks();
        bbc   = header->GetBbcAdcSumEastInner();

        hg_zdcx.Fill(zdcx);
        hg_bbcEin.Fill(bbc);
        hg_nGlob.Fill(nglob);
        hg_runPeriod.Fill(runid);

        double nTowp2   {0.};
        double nTow4GeV {0.};

        double maxTow_Et {0};
        double maxTow_phi{0};
        /* if (triggerId=500206) { */
        TList *l_tows = Reader.GetListOfSelectedTowers();//GetTowers();
        TIter next_tow(l_tows);
        while (TStarJetPicoTower* tower = (TStarJetPicoTower *) next_tow()) {
            double Et = tower->GetEt();
            if (Et > maxTow_Et) {
                maxTow_Et  = Et;
                maxTow_phi = tower->GetPhi();
            }
            if (Et > 0.2) nTowp2   += 1.;
            if (Et > 4.0) nTow4GeV += 1.;
            data_Et.fill(Et);
        }
        data_nETp2.fill(nTowp2);
        data_nET4GeV.fill(nTow4GeV);

        double nch       {0.};
        double nch_trans {0.};

        TList *p_prim_trks = Reader.GetListOfSelectedTracks();
        TIter nextptrk(p_prim_trks);
        while (TStarJetPicoPrimaryTrack* p_track = (TStarJetPicoPrimaryTrack *) nextptrk()) {
            data_dca       .fill(p_track->GetDCA());
            if (p_track->GetDCA()<1.) {
                nch += 1.;
                data_nhitsfit  .fill(p_track->GetNOfFittedHits());
                data_nhitsposs .fill(p_track->GetNOfPossHits());
                data_nhitsrat  .fill( (float)p_track->GetNOfFittedHits()/p_track->GetNOfPossHits());
                data_pt        .fill(p_track->GetPt());
            }

            if (maxTow_Et > 4. && is_trans(p_track->GetPhi(), maxTow_phi)) {
                data_dcaT       .fill(p_track->GetDCA());
                if (p_track->GetDCA()<1.) {
                    nch_trans += 1.;
                    data_nhitsfitT  .fill(p_track->GetNOfFittedHits());
                    data_nhitspossT .fill(p_track->GetNOfPossHits());
                    data_nhitsratT  .fill( (float)p_track->GetNOfFittedHits()/p_track->GetNOfPossHits());
                    data_ptT        .fill(p_track->GetPt());
                }
            }
        }
        data_nch.fill(nch);
        if (maxTow_Et > 4.) data_nchT.fill(nch_trans);
    }

    hg_zdcx      .Write();
    hg_bbcEin    .Write();
    hg_nGlob     .Write();
    hg_runPeriod .Write();

    data_Et        .write();
    data_nETp2     .write();
    data_nET4GeV   .write();

    data_nhitsfit  .write();
    data_nhitsposs .write();
    data_nhitsrat  .write();
    data_dca       .write();
    data_pt        .write();
    data_nch       .write();

    data_nhitsfitT  .write();
    data_nhitspossT .write();
    data_nhitsratT  .write();
    data_dcaT       .write();
    data_ptT        .write();
    data_nchT       .write();


    fout.Save();

    cout << "finished. Recommed: " << endl;
    cout << " scp wsu:/tier2/home/groups/rhi/vverkest/pAu_analysis/" << outFile << " ." << endl;

    return 0;
}
