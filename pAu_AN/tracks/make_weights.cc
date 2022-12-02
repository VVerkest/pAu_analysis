root -l<<EOF
    .L ../loc_lib/loc_libs.h
    .L TF1_Tsallis.C
    //get the embedding pT spectra
    noiGetter got;
    auto truth = (THnSparseD*) got("make_sparse.root","sp_truth");
    auto hg_pt = (TH1D*) truth->Projection(_pt);

    TAxis* axis = hg_pt->GetXaxis();
    
    TsallisAll pp {"pp"};
    TsallisAll dAu{"dAu"};

    auto fout = new TFile("make_weights.root","recreate");
    array<TH1D*,6> hg_pp, hg_dAu;
    TH1D* h;
    for (int p=0;p<6;++p) {
        hg_pp[p] = (TH1D*) hg_pt->Clone(Form("w_%i_pp",p));
        h = hg_pp[p];  
        for (int k = 1; k <= h->GetNbinsX();++k) {
            /* if (p==0 && k==1) cout << " first: " << pp.integral(p, axis->GetBinLowEdge(k), axis->GetBinUpEdge(k)) << " and " << h->GetBinContent(k) << endl; */
            h->SetBinContent(k, pp.integral(p, axis->GetBinLowEdge(k), axis->GetBinUpEdge(k))/h->GetBinContent(k));
            h->SetBinError(k,0);
        }
        h->Write();

        hg_dAu[p] = (TH1D*) hg_pt->Clone(Form("w_%i_dAu",p));
        h = hg_dAu[p];  
        for (int k = 1; k <= h->GetNbinsX();++k) {
            h->SetBinContent(k, dAu.integral(p, axis->GetBinLowEdge(k), axis->GetBinUpEdge(k))/h->GetBinContent(k));
            h->SetBinError(k,0);
        }
        h->Write();
    }
    fout->Save();
    fout->Close();

EOF
