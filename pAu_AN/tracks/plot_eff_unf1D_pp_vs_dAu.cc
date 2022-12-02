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
    /* auto events = (THnSparseD*) got("hadd_sparse_wsu.root","sp_events"); */


    auto reco_pp  = (THnSparseD*) got("make_sparse_pp.root","sp_reco");
    auto truth_pp = (THnSparseD*) got("make_sparse_pp.root","sp_truth");

    auto reco_dAu  = (THnSparseD*) got("make_sparse_dAu.root","sp_reco");
    auto truth_dAu = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");

    set_track_cuts({tracks, reco_pp, reco_dAu});

    // for each ZDCx bin, unfold with 1d and 2d and see what the results are


    TH1D hg_pp {"hg_pp",   ";ZDCx bin; eff-1D", 9, 0.5, 9.5};
    TH1D hg_dAu {"hg_dAu", ";ZDCx bin; eff-1D", 9, 0.5, 9.5};
    // reco_pp eff.
    array<int,9> i0_zdcx { 5, 8,  11, 14, 17, 20, 23, 26,  5 }; // binning by 3 for zdcx
    array<int,9> i1_zdcx { 7, 10, 13, 16, 19, 22, 25, 28,  28};

    const int iRepUnfold { 4 };
    for (int i=0;i<9;++i) {
        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {tracks});

        // data
        auto h_fpp   = (TH1D*)  tracks->Projection(_pt);
        auto h_fdAu  = (TH1D*)  h_fpp->Clone("h_fdAu");
        double raw_val = h_fpp->Integral(3,-1);

        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {reco_pp, truth_pp});
        auto h_reco  = (TH1D*) reco_pp->Projection(_ptreco);
        auto h_truth = (TH1D*) truth_pp->Projection(_pt);  

        // div_by_W divides values and errors in h_reco by the tracks truth values
        // without any additional uncertainty from the truth embedding statistics
        div_by_W (h_reco, h_truth);
        div_by_W (h_fpp, h_reco);
        double v1, e1;
        v1 = h_fpp->IntegralAndError(3, h_fpp->GetNbinsX(), e1);
        hg_pp.SetBinContent(i+1, raw_val / v1);
        hg_pp.SetBinError  (i+1, raw_val / v1 * e1 / v1);
        delete h_fpp;
        delete h_reco;
        delete h_truth;

        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {reco_dAu, truth_dAu});
        h_reco  = (TH1D*) reco_dAu->Projection(_ptreco);
        h_truth = (TH1D*) truth_dAu->Projection(_pt);  

        div_by_W (h_reco, h_truth);
        div_by_W (h_fdAu, h_reco);
        v1 = h_fdAu->IntegralAndError(3, h_fdAu->GetNbinsX(), e1);
        hg_dAu.SetBinContent(i+1, raw_val / v1);
        hg_dAu.SetBinError  (i+1, raw_val / v1 * e1 / v1);
        
        delete h_reco;
        delete h_truth;
        delete h_fdAu;
    }

    noi_fmt(&hg_pp, {{"MarkerStyle",kOpenCircle, "MarkerColor", kBlue, "Title", "Eff using pp weighting"}});
    noi_fmt(&hg_dAu, {{"MarkerStyle",kOpenSquare, "MarkerColor", kBlack, "Title", "Eff using dAu weighting"}});

    noiPads pads{2};
    pads(0);
    hg_pp.Draw("PE");
    hg_dAu.Draw("PE same");
    pads(0)->BuildLegend();

    pads(1);
    auto _rat = (TH1D*) noiDivide(&hg_pp, &hg_dAu);
    _rat->Draw("PE");
    _rat->SetTitle("Ratio eff. pp:dAu");
    _rat->GetYaxis()->SetTitle("Ratio pp:dAu eff values");
    pads(1)->BuildLegend();

    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");
EOF

