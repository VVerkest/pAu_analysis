root -l <<EOF
    .x tu_loadlibs.C
    .L loc_lib.h

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kCMYK);

    tuGetter got{};
    auto tracks  = (THnSparseD*) got("hadd_sparse_wsu.root", "sp_tracks");
    auto events  = (THnSparseD*) got("hadd_sparse_wsu.root", "sp_events");
    auto zdcxsum = (THnSparseD*) got("hadd_sparse_wsu.root", "sp_zdcx");

    auto reco   = (THnSparseD*) got("make_sparse_dAu.root","sp_reco");
    auto truth  = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");

    set_range(_nhitrat, 6,  10, {tracks, reco} );
    set_range(_nhitfit, 16, 47, {tracks, reco} );
    set_range(_dca,     1,  5,  {tracks, reco} );



    // save: sumZDCx[iEA] -> dNch_PU[iEA]
    //       nEvents[iEA][iZDCx]
                   // TH1D* pt_{raw,corr}[ieta][iEA][iZDCx]
    // TH1D* Nch_vs_ZDCx_etaI_bbcI{raw,corr}[ieta][iEA] -> 30 of these
    // double sum_ZDCx_bbcI -> 10 of these
    // double nevents_bbcI  -> 10 of these
    
    TFile f_save ((tuStripExtension("$0")+".root").c_str(),"recreate");
    tuMsgTree msg_tree{false};
    msg_tree.slurp_file("$0");
    msg_tree.write();

    //data out to fill
    array< array<TH1D*, 10>, 3 > arr_NchZDCx_raw; // arr_NchZDCx[ieta][ibbc]
    array< array<TH1D*, 10>, 3 > arr_NchZDCx_cor; // arr_NchZDCx[ieta][ibbc]
    array< double, 10 > arr_sumZDCx;   // [ibbc]
    array< array<double, 6>, 10 > arr_sumevents; // [ibbc][izdcx]
                                    
    // fill in arr_sum_ZDCx
    set_range(_zdcx, 4, 22, {zdcxsum, events});
    for (int i=0;i<10;++i) {
        set_range(_bbc,i+1,i+1, {zdcxsum, events});
        arr_sumZDCx[i] = sparse_integral(zdcxsum).first;
        /* arr_sumevents[i] = sparse_integral(events).first; */
    }
    f_save.WriteObject(&arr_sumZDCx, "sumZDCx_perEA"); 
    f_save.WriteObject(&arr_sumZDCx, "sumEvents_perEA"); 
    
    array<int, 7> i0bin_ZDCx { 4, 7,  10, 13, 16, 19 };
    array<int, 7> i1bin_ZDCx { 7, 10, 13, 16, 19, 22 };
    tuBinVec bin_zdcx_d3 {{ 4000., 4000., 6,    22000. }};
    array<double,4> arr_eta  {{ -1., -0.3, 0.3, 1. }};

    TH1D Nch_cor {"Nch_cor",";EA_{bbc};d^2N_{ch}/d#eta/d#phi corrected", bins_EAbbc10, bins_EAbbc10 };
    TH1D Nch_raw {"Nch_raw",";EA_{bbc};d^2N_{ch}/d#eta/d#phi raw",       bins_EAbbc10, bins_EAbbc10 };

    for (int iEA=0;iEA<10;++iEA) {
        set_range(_bbc, iEA+1, iEA+1, {tracks,events});

        TH1D* sumpt_raw = tracks->Projection(_pt);
        sumpt_raw->SetName(Form("sumpt_raw_bbc%i",iEA));
        sumpt_raw->Reset();
        auto sumpt_cor = (TH1D*) sumpt_raw->Clone(Form("sumpt_cor_bbc%i",iEA));

        for (int ieta=0; ieta<3; ++ieta) {
            set_range(_eta, ieta+1, ieta+1, {tracks, reco, truth});
            double eta0 = arr_eta[ieta];
            double eta1 = arr_eta[ieta+1];

            arr_NchZDCx_raw[ieta][iEA] = new TH1D(Form("NchZDCx_raw_eta%i_EA%i",ieta,iEA),
                Form("EA_{decile}(%i);ZDCx; dNch(raw)/dEta/d#phi, #eta#in[%.1f,%.1f]", iEA, eta0, eta1), bin_zdcx_d3, bin_zdcx_d3);
            arr_NchZDCx_cor[ieta][iEA] = new TH1D(Form("NchZDCx_cor_eta%i_EA%i",ieta,iEA),
                Form("EA_{decile}(%i);ZDCx; dNch(cor)/dEta/d#phi, #eta#in[%.1f,%.1f]", iEA, eta0, eta1), bin_zdcx_d3, bin_zdcx_d3);

            auto hg_raw = arr_NchZDCx_raw[ieta][iEA];
            auto hg_cor = arr_NchZDCx_cor[ieta][iEA];

            for (int i=0;i<6;++i) {
                int i0 = i0bin_ZDCx[i];//3*i+4;
                int i1 = i1bin_ZDCx[i];//3*i+4;

                set_range(_zdcx, i0, i1, {tracks, reco, truth, events});
                arr_sumevents[iEA][i] = sparse_integral(events).first;

                auto n_events = sparse_integral(events).first;

                auto h_raw = (TH1D*) tracks->Projection(_pt);
                sumpt_raw->Add(h_raw);
                /* h_raw->Scale(1./(eta1-eta0)/n_events/2/M_PI); */

                double v_raw, e_raw;
                v_raw = h_raw->IntegralAndError(3, h_raw->GetNbinsX(), e_raw);
                hg_raw->SetBinContent (i+1, v_raw);
                hg_raw->SetBinError   (i+1, e_raw);



                // do the correction
                auto h_true = (TH1D*) truth->Projection(_pt);
                auto h_reco = (TH1D*) reco->Projection(_ptreco);
                div_by_W (h_reco, h_true);
                auto h_cor = (TH1D*) tuDivide(h_raw, h_reco);

                double v_cor, e_cor;
                v_cor = h_cor->IntegralAndError(3, -1, e_cor);

                sumpt_cor->Add(h_cor);

                hg_cor->SetBinContent (i+1, v_cor);
                hg_cor->SetBinError   (i+1, e_cor);

                delete h_raw;
                delete h_true;
                delete h_reco;
                delete h_cor;
            }
            hg_raw->Write();
            hg_cor->Write();
        }
    }
EOF

