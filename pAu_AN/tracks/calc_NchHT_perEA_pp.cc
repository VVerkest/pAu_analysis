root -l <<EOF
    .L ../loc_lib/loc_libs.h


    noiGetter got{};
    auto tracks  = (THnSparseD*) got("hadd_sparse_wsu_HT.root", "sp_tracks");
    auto events  = (THnSparseD*) got("hadd_sparse_wsu_HT.root", "sp_events");
    auto zdcxsum = (THnSparseD*) got("hadd_sparse_wsu_HT.root", "sp_zdcx");

    auto reco   = (THnSparseD*) got("make_sparse_pp.root","sp_reco");
    auto truth  = (THnSparseD*) got("make_sparse_pp.root","sp_truth");

    string inp1  = "$1";
    if (inp1 == "") inp1 = "2"; // triget are {{ 0.,    4.,     8.,    12.,  30. }};
    const int iEt = atoi(inp1.c_str());
    const int iEtGeV = static_cast<int>(tracks->GetAxis(_trigEt)->GetBinLowEdge(iEt));

    set_range(_trigEt, iEt, iEt, {tracks, events, zdcxsum});

    set_track_cuts({tracks,reco});

    // output:
    /* array<int, 7> i0bin_ZDCx { 4, 7,  10, 13, 16, 19 }; */
    /* array<int, 7> i1bin_ZDCx { 7, 10, 13, 16, 19, 22 }; */
    /* tuBinVec bin_zdcx_d3 {{ 4000., 4000., 6, 28000. }}; */
    array<double,4> arr_eta  {{ -1., -0.3, 0.3, 1. }};
    
    TFile f_save {Form("%s_%i.root",noiStripExtension("$0").c_str(),iEtGeV),"recreate"};
    noiMsgTree msg_tree{false};
    msg_tree.slurp_file("$0");
    msg_tree.write();

    auto blank_pt = tracks->Projection(_pt);
    blank_pt->Reset();
    blank_pt->SetTitle("");
    blank_pt->SetName("blank");

    array<TH1D*, 10> arrpt_raw;
    array<TH1D*, 10> arrpt_cor;
    array<double,10> arr_nevents;
    array<double,10> arr_sumzdcx;

    set_range(_zdcx, 1, 23, {events, zdcxsum} );
    for (int i=0;i<10;++i) {
        arrpt_raw[i] = (TH1D*) blank_pt->Clone(Form("pt_raw_bbc%i",i));
        arrpt_cor[i] = (TH1D*) blank_pt->Clone(Form("pt_cor_bbc%i",i));

        set_range(_bbc, i+1, i+1, {events, zdcxsum});
        arr_nevents [i] = sparse_integral(events).first;
        arr_sumzdcx [i] = sparse_integral(zdcxsum).first;
    }

    for (int izdcx=0;izdcx<5;++izdcx) {
        set_range(_zdcx, _i0_zdcx[izdcx], _i1_zdcx[izdcx], {tracks, reco, truth});

        for (int ieta=0; ieta<3; ++ieta) {
            set_range(_eta, ieta+1, ieta+1, {tracks, reco, truth});
            double eta0 = arr_eta[ieta];
            double eta1 = arr_eta[ieta+1];

            auto h_true = (TH1D*) truth->Projection(_pt);
            auto h_eff = (TH1D*) reco->Projection(_ptreco);
            div_by_W (h_eff, h_true);
            delete h_true;

            for (int iEA=0;iEA<10;++iEA) {
                set_range(_bbc, iEA+1, iEA+1, {tracks});
                auto h_raw = (TH1D*) tracks->Projection(_pt);
                auto h_cor = (TH1D*) noiDivide(h_raw, h_eff);

                arrpt_raw[iEA]->Add(h_raw);
                arrpt_cor[iEA]->Add(h_cor);

                delete h_raw;
                delete h_cor;
                cout << " izdcx:iEA  "<<izdcx<<":"<<iEA << endl;
            }
            delete h_eff;
        }
    }

    // calculate the final values
    TH1D* hg_raw = new TH1D("hg_raw",";EA_{BBC};d^{2}N_{ch}/d#eta/d#phi", bin_bbc, bin_bbc);    
    TH1D* hg_cor = new TH1D("hg_cor",";EA_{BBC};d^{2}N_{ch}/d#eta/d#phi", bin_bbc, bin_bbc);    
    for (int iEA=0; iEA<10; ++iEA) {
        double sum_val, sum_err; 

        const double delta_phi = (2./3.)*M_PI;
        const double delta_eta = 2.;
        arrpt_raw[iEA]->Scale(1./arr_nevents[iEA]/delta_phi/delta_eta);
        sum_val = arrpt_raw[iEA]->IntegralAndError(3, -1, sum_err);
        hg_raw->SetBinContent(iEA+1,sum_val);
        hg_raw->SetBinError  (iEA+1,sum_err);

        arrpt_cor[iEA]->Scale(1./arr_nevents[iEA]/delta_phi/delta_eta);
        sum_val = arrpt_cor[iEA]->IntegralAndError(3, -1, sum_err);
        hg_cor->SetBinContent(iEA+1,sum_val);
        hg_cor->SetBinError  (iEA+1,sum_err);
    }
    hg_raw->Write();
    hg_cor->Write();
    f_save.WriteObject(&arr_sumzdcx, "sumZDCx");
    f_save.WriteObject(&arr_nevents, "sumNevents");
EOF

