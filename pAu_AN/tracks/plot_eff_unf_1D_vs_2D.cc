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

    // for each ZDCx bin, unfold with 1d and 2d and see what the results are


    TH1D hg_1D {"hg_1D", ";ZDCx bin; eff-1D", 9, 0.5, 9.5};
    TH1D hg_2D {"hg_2D", ";ZDCx bin; eff-1D", 9, 0.5, 9.5};
    // reco eff.
    const int iRepUnfold { 4 };
    /* array<int, 7> i0bin_ZDCx { 4, 7,  10, 13, 16, 19 }; */
    /* array<int, 7> i1bin_ZDCx { 7, 10, 13, 16, 19, 22 }; */
    array<int,9> i0_zdcx { 5, 8,  11, 14, 17, 20, 23, 26,  5 }; // binning by 3 for zdcx
    array<int,9> i1_zdcx { 7, 10, 13, 16, 19, 22, 25, 28,  28};

    for (int i=0;i<9;++i) {
        // data
        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {tracks, events}, false);
        auto h_raw1d  = (TH1D*)  tracks->Projection(_pt);
        auto h_raw2d  = (TH1D*)  h_raw1d->Clone(noiUniqueName());

        const double raw_val = h_raw1d->Integral(3,-1);

        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {reco, truth}, false);
        auto h_reco  = (TH1D*) reco->Projection(_ptreco);
        auto h_truth = (TH1D*) truth->Projection(_pt);  

        // div_by_W divides values and errors in h_reco by the tracks truth values
        // without any additional uncertainty from the truth embedding statistics
        div_by_W (h_reco, h_truth);
        div_by_W (h_raw1d, h_reco);
        /* h_raw1d->Scale(1./raw_val); */
        double v1, e1;
        v1 = h_raw1d->IntegralAndError(3, h_raw1d->GetNbinsX(), e1);
        hg_1D.SetBinContent(i+1, raw_val / v1);
        hg_1D.SetBinError  (i+1, raw_val / v1 * e1 / v1);
        delete h_reco;
        delete h_truth;
        delete h_raw1d;

        auto R = (TH2D*) reco->Projection(_pt, _ptreco);
        auto T = (TH1D*) truth->Projection(_pt);
        auto M = (TH1D*) R->ProjectionX();
        auto rooUnfRes = new RooUnfoldResponse(M,T,R,noiUniqueName());
        auto bayes = new RooUnfoldBayes(rooUnfRes, h_raw2d, iRepUnfold);
        auto unf = (TH1D*) bayes->Hreco();
        v1 = unf->IntegralAndError(3, h_raw2d->GetNbinsX(), e1);
        cout << " v1: " << v1 << endl;
        hg_2D.SetBinContent(i+1, raw_val / v1 );
        hg_2D.SetBinError  (i+1, raw_val / v1 * e1 / v1);

        delete R, T, M, rooUnfRes, bayes, unf;
    }

    noi_fmt(&hg_1D, {{"MarkerStyle",kOpenCircle, "MarkerColor", kBlue, "Title", "Eff using single track eff correction"}});
    noi_fmt(&hg_2D, {{"MarkerStyle",kOpenSquare, "MarkerColor", kBlack, "Title", "Eff using RooUnfoldBayes"}});

    noiPads pads{2};
    pads(0);
    hg_1D.Draw("PE");
    hg_2D.Draw("PE same");
    pads(0)->BuildLegend();

    pads(1);
    auto _rat = (TH1D*) noiDivide(&hg_1D, &hg_2D);
    _rat->Draw("PE");
    _rat->SetTitle("Ratio eff. single:bayes");
    _rat->GetYaxis()->SetTitle("Ratio 1D:2D eff values");
    pads(1)->BuildLegend();

    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");
EOF

