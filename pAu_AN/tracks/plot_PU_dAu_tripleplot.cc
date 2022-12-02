root -l <<EOF
    .L ../loc_lib/loc_libs.h

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kCMYK);

    noiGetter got{};
    auto tracks = (THnSparseD*) got("hadd_sparse_wsu.root","sp_tracks");
    auto events = (THnSparseD*) got("hadd_sparse_wsu.root","sp_events");

    auto reco   = (THnSparseD*) got("make_sparse_dAu.root","sp_reco");
    auto truth = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");

    set_track_cuts({tracks,reco});

    array<TH1D*,4> hg_raw;
    array<TH1D*,4> hg_cor;

    array<TF1*,4> fit;
    array<double[2],4> pvals;
    array<double[2],4> evals;

    tuBinVec bin_zdcx_d3 {{ 5000., 5000., 5,    20000. }};

    vector<int> shapes_raw { kOpenCircle, kOpenSquare, kOpenDiamond, kOpenStar };
    vector<int> shapes_cor { kFullCircle, kFullSquare, kFullDiamond, kFullStar };
    vector<int> colors { kGreen+2, kBlue, kBlack, kRed };

    array<double,4> arr_eta  {{ -1., -0.3, 0.3, 1. }};
    for (int eta_bin=0; eta_bin<4; ++eta_bin) {

        double eta0, eta1, ieta0, ieta1; 
        if (eta_bin<3) {
            eta0  = arr_eta[eta_bin];
            eta1  = arr_eta[eta_bin+1];
            ieta0 = eta_bin+1;
            ieta1 = eta_bin+1;
        } else {
            eta0 = -1.;
            eta1 = 1.;
            ieta0 = 1;
            ieta1 = 3;
        }

        hg_raw[eta_bin] = new TH1D(Form("RawData_eta_%.1f,%.1f",eta0,eta1),
                Form(";ZDCx; dNch/dEta/d#phi, #eta#in[%.1f,%.1f]", eta0, eta1),
                bin_zdcx_d3, bin_zdcx_d3);
        hg_cor[eta_bin] = new TH1D(Form("CorrData_eta_%.1f,%.1f",eta0,eta1),
                Form(";ZDCx; dNch/dEta/d#phi, #eta#in[%.1f,%.1f]", eta0, eta1),
                bin_zdcx_d3, bin_zdcx_d3);

        set_range(_eta, ieta0, ieta1, {tracks, reco, truth});

        for (int i=0;i<5;++i) {

            set_range(_zdcx, _i0_zdcx[i], _i1_zdcx[i], {tracks, reco, truth, events});
            auto n_events = sparse_integral(events);

            auto h_raw = (TH1D*) tracks->Projection(_pt);
            h_raw->Scale(1./(eta1-eta0)/n_events.first/2/M_PI);

            double v_raw, e_raw;
            v_raw = h_raw->IntegralAndError(3, h_raw->GetNbinsX(), e_raw);
            hg_raw[eta_bin]->SetBinContent (i+1, v_raw);
            hg_raw[eta_bin]->SetBinError   (i+1, e_raw);

            // do the correction
            auto h_true = (TH1D*) truth->Projection(_pt);
            auto h_reco = (TH1D*) reco->Projection(_ptreco);
            div_by_W (h_reco, h_true);
            auto h_cor = (TH1D*) noiDivide(h_raw, h_reco);

            double v_cor, e_cor;
            v_cor = h_cor->IntegralAndError(3, -1, e_cor);

            hg_cor[eta_bin]->SetBinContent (i+1, v_cor);
            hg_cor[eta_bin]->SetBinError   (i+1, e_cor);

            delete h_raw;
            delete h_true;
            delete h_reco;
            delete h_cor;
        }

        fit[eta_bin] = new TF1(Form("line_%i",eta_bin),"pol1");
        fit[eta_bin]->SetLineColorAlpha(colors[eta_bin],0.8);
        fit[eta_bin]->SetLineStyle(2);
        fit[eta_bin]->SetLineWidth(2);
        hg_cor[eta_bin]->Fit(Form("line_%i",eta_bin));

        auto params = fit[eta_bin]->GetParameters();
        auto errors = fit[eta_bin]->GetParErrors();
        pvals[eta_bin][0] = params[0];
        pvals[eta_bin][1] = params[1];

        evals[eta_bin][0] = errors[0];
        evals[eta_bin][1] = errors[1];

        noi_fmt( hg_raw[eta_bin], {{"MarkerColor", colors[eta_bin], "MarkerStyle", shapes_raw[eta_bin]}});
        noi_fmt( hg_cor[eta_bin], {{"MarkerColor", colors[eta_bin], "MarkerStyle", shapes_cor[eta_bin]}});
    }

    noiPads pads{3, {800,1200, kBottom, 0.10, 0.02, kTop, 0.02, 0.02}};
    pads(0);

    TH1D hg_m {"hg_m",";Eta bins [-1,-.3] [-.3,.3] [.3,1] [-1,1]; PU fit \"m\"", 4, 0.5, 4.5 };
    TH1D hg_b {"hg_b",";Eta bins [-1,-.3] [-.3,.3] [.3,1] [-1,1]; PU fit \"b\"", 4, 0.5, 4.5 };

    for (int i=0; i<4; ++i) {
        hg_cor[i]->GetYaxis()->SetRangeUser(0.9, 1.2);
        hg_cor[i]->Draw( i == 0 ? "PE" : "PE same" );

        hg_m.SetBinContent(i+1, pvals[i][1]);
        hg_m.SetBinError  (i+1, evals[i][1]);

        hg_b.SetBinContent(i+1, pvals[i][0]);
        hg_b.SetBinError  (i+1, evals[i][0]);
    }
    pads(0)->BuildLegend();


    pads(1);
    hg_m.Draw("PE");
    pads(2);
    hg_b.Draw("PE");

    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");
EOF

