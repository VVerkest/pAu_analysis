root -l <<EOF
    .L ../loc_lib/loc_libs.h

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kCMYK);

    vector<int> colors { 1179, 1207, 1235, 1263, 1291, 1319, 1347, 1375, 1403, 1431, kRed };
    vector<int> shapes { kOpenCircle, kOpenSquare, kOpenDiamond,
                         kOpenCircle, kOpenSquare, kOpenDiamond,
                         kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCircle, kOpenSquare };

    noiGetter got{};


    auto reco  = (THnSparseD*) got("make_sparse_dAu.root","sp_reco");
    auto truth = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");

    set_track_cuts({reco});

    /* set_range(_nhitrat, 6,  10, reco); */
    /* set_range(_nhitfit, 16, 47, reco); */

    array<TH1D*,11> eff_arr;
    array<TF1*,11> fit;
    array<double[2],11> pvals;
    array<double[2],11> evals;
    for (int i=0;i<11;++i) {
        int i0, i1;
        if (i<10) {
            i0 = i+1;
            i1 = i+1;
        } else {
            i0 = 1;
            i1 = 10;
        }
        set_range(_bbc,i0,i1, {reco, truth});

        eff_arr[i] = emb_eff_per_zdcx_binnedby3 (reco, truth, Form("EA_%i",i+1));
        noi_fmt(eff_arr[i], {{"MarkerStyle",shapes[i],"MarkerColor", colors[i],  "Title",
                Form("EA bin (%i-%i)",i0,i1), "yAxisTitle","Tracking efficiency", "yAxisRangeLo",0., "yAxisRangeHi",0.9}});
        fit[i] = new TF1(Form("line%i",i),Form("pol1"));
        fit[i]->SetLineColorAlpha(colors[i],0.8);
        fit[i]->SetLineStyle(1);
        eff_arr[i]->Fit(Form("line%i",i));
        auto params = fit[i]->GetParameters();
        auto errors = fit[i]->GetParErrors();
        pvals[i][0] = params[0];
        pvals[i][1] = params[1];

        evals[i][0] = errors[0];
        evals[i][1] = errors[1];
    }

    noiPads pads{3, {800,1200, kBottom, 0.10, 0.02, kTop, 0.02, 0.02}};
    pads(0);

    TH1D hg_m {"hg_m",";EA bin; eff. fit \"m\"", 11, 0.5, 11.5 };
    TH1D hg_b {"hg_b",";EA bin; eff. fit \"b\"", 11, 0.5, 11.5 };

    for (int i=0;i<11;++i) {
        int i0, i1;
        if (i<10) {
            i0 = i+1;
            i1 = i+1;
        } else {
            i0 = 1;
            i1 = 10;
        }
        eff_arr[i]->Draw(i==0?"PE":"PE same");
        noiDrawTLatex(Form("Fit EA bin %i-%i: Eff = %6.3f^{#pm%6.4f} + %13.9g^{#pm%10.9g} * zdcX", i0, i1,
                pvals[i][0], evals[i][0], pvals[i][1], evals[i][1]), 5000., 0.05+(0.05*i), 
                {{"TextColor",colors[i],"LineStyle",2,"TextSize",12}});
        hg_m.SetBinContent(i+1, pvals[i][1]);
        hg_m.SetBinError  (i+1, evals[i][1]);

        hg_b.SetBinContent(i+1, pvals[i][0]);
        hg_b.SetBinError  (i+1, evals[i][0]);
    }
    pads(1);
    hg_m.Draw("PE");
    pads(2);
    hg_b.Draw("PE");

    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");

EOF
