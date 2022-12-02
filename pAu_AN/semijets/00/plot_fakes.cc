root -l <<EOF
    .x tu_loadlibs.C
    .L RooUnfold_RFM.h

    auto _all     = RooUnfold_RFM_Array("hadd_jetemb_$1_$2.root", "response120to180");
    auto _A_arr   = RooUnfold_RFM_Array("hadd_jetemb_$1_$2.root", "response120to180_A");
    auto _B_arr   = RooUnfold_RFM_Array("hadd_jetemb_$1_$2.root", "response120to180_B");

    _A_arr.scrub($3);
    _B_arr.scrub($3);
    
    tuBinVec bin_M = {{ 10., 15., 20., 25., 30., 35., 40., 45.}};
    /* tuBinVec bin_M = {{ 10., 15., 20., 25., 30., 35., 40., 45.}}; */

    /* tuBinVec bin_T = {{ 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55. }}; */
    tuBinVec bin_T = {{ 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55. }};

    auto A_arr = _A_arr.rebin(bin_M, bin_T);
    auto B_arr = _B_arr.rebin(bin_M, bin_T);

    auto canv = A_arr.plot_fakes();
    canv->SaveAs("pdf/A_$1$2_$3_fakes.pdf");

    canv = B_arr.plot_fakes();
    canv->SaveAs("pdf/B_$1$2_$3_fakes.pdf");

    canv = A_arr.plot_misses();
    canv->SaveAs("pdf/A_$1$2_$3_misses.pdf");

    canv = B_arr.plot_misses();
    canv->SaveAs("pdf/B_$1$2_$3_misses.pdf");

EOF
