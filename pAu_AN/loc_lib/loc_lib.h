#ifndef LOC_LIB_H
#define LOC_LIB_H

// local library for jet unfolding

#include "noiBinVec.h"
#include "pAu_bins.h"

#include "noi_fnc.h"

void div_by_W (TH1D* hg, TH1D* w) {
    for (int i=1; i<=hg->GetNbinsX(); ++i) {
        double W = w->GetBinContent(i);
        if (W==0) {
            hg->SetBinContent(i, 0.);
            hg->SetBinError  (i, 0.);
        } else {
            hg->SetBinContent(i, hg->GetBinContent(i)/W);
            hg->SetBinError  (i, hg->GetBinError  (i)/W);
        }
    }
}

void set_range(int axis, int bin0, int bin1, vector<THnSparseD*> vec_sparse, bool print=false) {
    bool first = true;
    for (auto& sparse : vec_sparse) {
        TAxis* ax = sparse->GetAxis(axis);
        sparse->GetAxis(axis)->SetRange(bin0, bin1);
        if (first && print) {
            first = false;
            cout << Form("Set range(%2i) to %7.2f-%7.2f",axis,ax->GetBinLowEdge(bin0),ax->GetBinUpEdge(bin1)) << endl;
        }
    }
}

pair<double,double> sparse_integral(THnSparseD* sp) {
    auto hg = sp->Projection(0);
    double _val, _err;
    _val = hg->IntegralAndError(0, hg->GetNbinsX(), _err);
    delete hg;
    return {_val,_err};
}

void set_track_cuts(vector<THnSparseD*> data) {
    set_range(_dca,     1,  5,  data);
    set_range(_nhitfit, 20, 48, data);
    set_range(_nhitrat, 6,  10, data);
}

TH1D* emb_eff_per_zdcx_binnedby3(THnSparseD* reco, THnSparseD* truth, string name) {
    array<int,5> i0_zdcx { 5, 8,  11, 14, 17 }; // binning by 3 for zdcx
    array<int,5> i1_zdcx { 7, 10, 13, 16, 19 };
   
    tuBinVec bins_5by3to20 {{ 5000., 5000., 5, 20000 }};
    //
    /* array<int,7> i0_zdcx { 4, 7, 10, 13, 16, 19, 22 }; */
    /* array<int,7> i1_zdcx { 4, 7, 10, 13, 16, 19, 22 }; */
    
    /* tuBinVec bins_5by3to20 {{ 5000., 5000., 5, 20000 }}; */
    TH1D* hg_out = new TH1D(name.c_str(), ";ZDCx;Tracking Efficiency", bins_5by3to20, bins_5by3to20 );

    for (int i=0; i<5; ++i) {
        set_range(_zdcx, i0_zdcx[i], i1_zdcx[i], {reco, truth});
        auto h_reco  = (TH1D*) reco->Projection(_pt);  h_reco ->SetName(noiUniqueName());
        auto h_truth = (TH1D*) truth->Projection(_pt); h_truth->SetName(noiUniqueName());
        double _val, _err;
        _val = h_reco->IntegralAndError(3,h_reco->GetNbinsX(),_err);
        cout << " i("<<i<<") val("<<_val<<" +/- "<<_err<<")   ----> " << _err/_val << endl;
        double norm = h_truth->Integral(3,h_truth->GetNbinsX());
        delete h_reco;
        delete h_truth;
        if (norm == 0) continue;
        _val /= norm;
        _err /= norm;
        hg_out->SetBinContent(i+1,_val);
        hg_out->SetBinError  (i+1,_err);
    }
    return hg_out;
}

// make tuBinVec avialable locally

#endif
