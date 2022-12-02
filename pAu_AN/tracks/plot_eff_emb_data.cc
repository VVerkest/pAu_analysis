root -l <<EOF
    .L ../loc_lib/loc_libs.h

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kCMYK);

    /* vector<int> colors { 1179, 1207, 1235, 1263, 1291, 1319, 1347, 1375, 1403, 1431, kRed }; */
    /* vector<int> shapes { kOpenCircle, kOpenSquare, kOpenDiamond, */
    /*                      kOpenCircle, kOpenSquare, kOpenDiamond, */
    /*                      kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCircle, kOpenSquare }; */


    noiGetter got{};
    auto tracks = (THnSparseD*) got("hadd_sparse_wsu.root","sp_tracks");
    auto events = (THnSparseD*) got("hadd_sparse_wsu.root","sp_events");

    auto reco  = (THnSparseD*) got("make_sparse_dAu.root","sp_reco");
    auto truth = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");

    set_track_cuts({tracks,reco});

    array<TH1D*, 7> hg_eff;
    array<TH1D*, 7> hg_data;

    array<int,9> i0_zdcx { 5, 8,  11, 14, 17, 20, 23, 26,  5 }; // binning by 3 for zdcx
    array<int,9> i1_zdcx { 7, 10, 13, 16, 19, 22, 25, 28,  28};
    for (int i=0;i<9;++i) {
        /* int i0 = i0bin_ZDCx[i];//3*i+4; */
        /* int i1 = i1bin_ZDCx[i];//3*i+4; */
        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {reco, truth}, false);
        auto h_reco  = (TH1D*) reco->Projection(_ptreco);
        auto h_truth = (TH1D*) truth->Projection(_pt);  
        div_by_W (h_reco, h_truth);
        h_reco->SetName(Form("reco_%i",i));
        hg_eff[i] = h_reco;
        delete h_truth;

        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {tracks, events}, false);
        auto h_raw  = (TH1D*)  tracks->Projection(_pt);
        auto h_event = (TH1D*) events->Projection(_zdcx);
        /* auto nevents = h_event->Integral(); */
        h_raw->Scale(1./h_event->Integral(3,-1));
        delete h_event;
        h_raw->SetName(Form("data_%i",i));
        hg_data[i] = h_raw;
    }

    TH2D recon { "recon", "eff. change;ZDCx-bin emb.; ZDCx-bin data", 9, 0.5, 9.5, 9, 0.5, 9.5 };
    noiPads pads{1, {{kRight, 0.1}} };
    pads(0);

    array<array<float,9>,9> aval;
    for (int i=0;i<9;++i) { // i for data
    for (int j=0;j<9;++j) { // j for emb
        auto hdata = (TH1D*) hg_data[i]->Clone("hg_eff");
        double val0 = hdata->Integral(3,-1);
        div_by_W(hdata,hg_eff[j]);
        double val1 = hdata->Integral(3,-1);
        delete hdata ;
        aval[j][i] = val0/val1;
        recon.Fill( 1.+j, 1.+i, aval[j][i] );
    }
    }
    recon.Draw("COLZ");
    for (int i=0;i<9;++i) { // i for data
    for (int j=0;j<9;++j) { // j for emb
        noiDrawTLatex( Form("%8.6f", aval[j][i]), 0.65+ j, 1.+i, {{"TextColorAlpha", kWhite,1., "TextSize", 19}});
    }}
    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");
EOF

