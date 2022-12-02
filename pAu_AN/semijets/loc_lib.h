#ifndef LOC_LIB_H
#define LOC_LIB_H

// local library for jet unfolding

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

#endif
