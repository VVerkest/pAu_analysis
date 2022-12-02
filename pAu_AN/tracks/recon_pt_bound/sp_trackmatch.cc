#include "events.h"
#include <sstream>
#include <fstream>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "tuClass.h"
#include "myDiscreteLinFnc.h"
#include "TRandom3.h"
#include "TFile.h"
#include "THnSparse.h"

// make a track matching histogram

using namespace std;

int main(int nargs, char** argv) { // eff_study_glob
    cout << " Running fn \"sp_trackmatch\"" << endl;

    string o_name_tag  { (nargs>1) ? argv[1] : "sp_trackmatch" };
    int    n_events    { (nargs>2) ? atoi(argv[2]) : -1              };
    string inp_list    { (nargs>3) ? argv[3] : "in-lists/list_global.list" };

    ofstream log;
    log.open((o_name_tag + ".log").c_str());
    log << "Starting output."  << endl
        << "Command line input:" << endl;
    for (int i{0}; i<nargs; ++i) log << "arg("<<i<<")  " << argv[i] << endl;
    log << endl << endl;

    events dat{log, n_events, inp_list};
    TFile fout  { (o_name_tag+".root").c_str(), "recreate" };

    /* /1* istringstream options ( _options ); *1/ */
    /* /1* int n_options = 0; *1/ */
    /* /1* string arg; *1/ */

    ioMsgTree msg_tree{false};
    msg_tree.slurp_file("src/sp_trackmatch.cc");
    msg_tree.write();

    /* /1* ioBinVec bins_pt    {{ 0.0, 0., 10, 1., 2., 5, 30 }}; *1/ */
    /* ioBinVec bins_pt    {{ 0.0, 0., 1200, 30.}}; */

    /* ioBinVec bins_zdcx2 {{ 4000., 4000., 3, 22000.  }}; */
    /* TProfile2D p2_match   { "p2_match","Ratio matching;ZDCx;#it{p}_{T}", 6, 4000., 22000., bins_pt, bins_pt }; */
    /* array<TProfile2D*, 6> p2_arr; */
    /* /1* TProfile2D p2_match_s { "p2_match_s","Ratio matching;ZDCx;#it{p}_{T}", 6, 4000., 22000., bins_pt, bins_pt, "s" }; *1/ */

    ioBinVec bin_zdcx     {{ 2000.,  2000., 22, 30000.  }};
    ioBinVec bin_bbc      {{ 0., 0., 100., 64000 }};
    ioBinVec bin_vz       {{ -10., -10., 20., 10. }};
    ioBinVec bin_vz_big   {{ -50., -50., 100., 50. }};
    ioBinVec bin_nglob    {{  0., 0., 26, 2600 }};
    ioBinVec bin_dca      {{ 0., 0., 30, 3. }};
    ioBinVec bin_pt       {{ 0., 0., 300, 30. }};
    ioBinVec bin_id       {{ -0.5, -0.5, 6, 5.5 }};

    ioBinVec bin_zdcx_corse     {{ 0.,  4000., 10000., 16000., 22000., 50000.  }};

    ioBinVec bin_dca_sparse {{ 0, 0, 3., 3 }};

    // Histogram declarations here:
    // event
    TH1D h_vz      { "vz" , ";Vz;N_{events}", bin_vz, bin_vz };
    TH1D h_vz_big  { "vz_big" , ";Vz;N_{events}", bin_vz_big, bin_vz_big };
    TH1D h_bbc     { "bbc" , ";EA_{BBC};N_{events}", bin_bbc, bin_bbc };
    TH1D h_zdcx    { "zdcx", ";ZDCx;N_{events}",     bin_zdcx, bin_zdcx };
    TH1D h_nglob   { "nglob", ";ZDCx;N_{events}",    bin_nglob, bin_nglob };

    TH2D h_nglob_zdcx  { "nglob_zdcx", ";n_glob;zdcx",    bin_nglob, bin_nglob, bin_zdcx, bin_zdcx };
    TProfile pr_glob_zdcx { "pr_zdcx_nglob", ";ZDCx; nglobal tracks", bin_zdcx, bin_zdcx };

    // axes:
    // 0. particle-id
    // 1. n_Global_tracks
    // 2. pT_Truth
    // 3. pT_Matched
    // 4. DCA
    int nbins[5];
    nbins[0] = bin_id;
    nbins[1] = bin_nglob;
    nbins[2] = bin_zdcx_corse;
    nbins[3] = bin_pt;
    nbins[4] = bin_pt;
    nbins[5] = bin_dca;

    double hopper[6];
    
    auto truth = new THnSparseD ("truth_tracks", "truth_tracks;particle-id;n_Global_tracks;ZDCx;pT_{Truth};",
        4, nbins, NULL, NULL);
    truth->SetBinEdges(0,bin_id);
    truth->SetBinEdges(1,bin_nglob);
    truth->SetBinEdges(2,bin_zdcx_corse);
    truth->SetBinEdges(3,bin_pt);

    auto match = new THnSparseD ("match_tracks", "match_tracks;particle-id;n_Global_tracks;ZDCx;pT_{Truth};pT_{Measured};Measured DCA;",
        6, nbins, NULL, NULL);
    match->SetBinEdges(0,bin_id);
    match->SetBinEdges(1,bin_nglob);
    match->SetBinEdges(2,bin_zdcx_corse);
    match->SetBinEdges(3,bin_pt);
    match->SetBinEdges(4,bin_pt);
    match->SetBinEdges(5,bin_dca);

    // OK, now read in the 5 sigma outlier terms
    TFile *s_current = gDirectory->GetFile();
    TFile *fin = new TFile("in_data/bin_fits.root","read");
    vector<double> *_bins, *m_lo, *b_lo, *m_hi, *b_hi;
    fin->GetObject("bins",_bins);
    array<myDiscreteLinFnc*,6> fn_lo, fn_hi;
    for (int i=0; i<6; ++i) {
        fin->GetObject(Form("m_lo_%s",io_geant05_ascii(i)), m_lo);
        fin->GetObject(Form("m_hi_%s",io_geant05_ascii(i)), m_hi);
        fin->GetObject(Form("b_lo_%s",io_geant05_ascii(i)), b_lo);
        fin->GetObject(Form("b_hi_%s",io_geant05_ascii(i)), b_hi);
        fn_lo[i] = new myDiscreteLinFnc(*_bins, *b_lo, *m_lo);
        fn_hi[i] = new myDiscreteLinFnc(*_bins, *b_hi, *m_hi);
    };
    s_current->cd();    

    while (dat.next()) {
        double val_vz   = dat.mu_event->vz;
        double val_zdcx = dat.mu_event->ZDCx;
        double val_bbc  = dat.mu_event->BBC_Ein;
        double val_nglob = dat.mu_event->nGlobalTracks;

        h_vz    .Fill(val_vz);
        h_vz_big.Fill(val_vz);
        h_zdcx  .Fill(val_zdcx);
        h_bbc   .Fill(val_bbc);
        h_nglob .Fill(val_nglob);

        h_nglob_zdcx.Fill(val_nglob, val_zdcx);
        pr_glob_zdcx.Fill( val_zdcx, val_nglob);

        for (auto& mc : dat.iter_mc_track) {
            int k = io_geant05(mc.geantId);
            double T = mc.pt;

            hopper[0] = (double) k;
            hopper[1] = val_nglob;
            hopper[2] = val_zdcx;
            hopper[3] = T;
            if (mc.id >= 0) {
                double M = dat.get_track(mc.id)->pt;
                if (M < (*fn_lo[k])(T) || M > (*fn_hi[k])(T))  continue; // just cut the event overall
                truth->Fill(hopper,1.);
                if (dat.get_track(mc.id)->dcaXYZ>3 || 
                    dat.get_track(mc.id)->nHitsFit<15 ||
                    (((float) dat.get_track(mc.id)->nHitsFit /
                     (float) dat.get_track(mc.id)->nHitsPoss) < 0.52)
                ) continue;
                double dca = dat.get_track(mc.id)->dcaXYZ;
                hopper[4] = M;
                hopper[5] = dca;
                match->Fill(hopper,1.);
            } else {
                truth->Fill(hopper,1.);
            }
        }
    }
    
     cout << " match " << match->Projection(1)->Integral() << endl;

    h_vz      .Write();
    h_vz_big  .Write();
    h_bbc     .Write();
    h_zdcx    .Write();
    h_nglob   .Write();

    h_nglob_zdcx.Write();
    pr_glob_zdcx.Write();

    truth ->Write();
    match ->Write();
    
    /* // Wrap-up work here: */
    dat.log  << " Done running events" << endl;
    cout     << " Done running events" << endl;

    fout.Close();
    log.close();

    cout << " Done running, suggest: " << endl;
    cout << " scp wsu:/wsu/home/hg/hg80/hg8093/gs_track_emb/"<<o_name_tag<<".root ."<< endl;
    return 0;
}
