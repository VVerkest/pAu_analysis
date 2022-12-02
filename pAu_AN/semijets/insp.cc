root -l <<EOF
    .x tu_loadlibs.C
    .L RooUnfold_RFM.h

    auto _all = RooUnfold_RFM_Array("hadd_jetemb.root", "response0to60");
    
    tuBinVec bin_M = {{ 5.,  10, 55. }};
    tuBinVec bin_T = {{ 15., 45. }};

    cout << _all[0].miss->GetNbinsX()<< endl;

    auto new_all = _all.rebin(bin_M, bin_T);

    cout << new_all[0].miss->GetNbinsX()<< endl;

EOF
