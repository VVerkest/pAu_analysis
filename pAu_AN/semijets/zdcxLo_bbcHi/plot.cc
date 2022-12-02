root -l <<EOF
    .x tu_loadlibs.C
    .L RooUnfold_RFM.h

    auto _all = RooUnfold_RFM_Array("hadd_jetemb.root", "response120to180");
    auto _A_arr   = RooUnfold_RFM_Array("hadd_jetemb.root", "response120to180_A");
    auto _B_arr   = RooUnfold_RFM_Array("hadd_jetemb.root", "response120to180_B");

    _A_arr.scrub($1);
    _B_arr.scrub($1);
    
    tuBinVec bin_M = {{ 10., 15., 20., 25., 30., 35., 40., 45.}};
    /* tuBinVec bin_M = {{ 10., 15., 20., 25., 30., 35., 40., 45.}}; */

    /* tuBinVec bin_T = {{ 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55. }}; */
    tuBinVec bin_T = {{ 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55. }};

    auto A_arr = _A_arr.rebin(bin_M, bin_T);
    auto B_arr = _B_arr.rebin(bin_M, bin_T);

    auto pa = A_arr.plot();
    pa.pdf("pdf/A_$1");
    tuPause();
    pb = B_arr.plot();
    pb.pdf("pdf/B_$1");
    tuPause();

    /* A_arr.set_sqrt(); */
    /* B_arr.set_sqrt(); */

    auto _A = A_arr.hadd();
    auto _B = B_arr.hadd();


    // use B to unfold A and check
    auto hadd_resp_B = _B();
    auto measured = _A.measured();
    auto truth    = _A.truth();

    /* truth->Draw("PE"); */
    /* tuPause(); */

    tuPads pads{1};
    pads(0)->SetGridy();
    TLegend *leg = new TLegend(0.0008687949,0.3057933,0.1507785,0.738973,NULL,"brNDC");
    TAxis* x = measured->GetXaxis();
    double xlo = x->GetBinLowEdge(1);
    double xhi = x->GetBinUpEdge(x->GetNbins());
    TAxis* y = truth->GetXaxis();
    double ylo = y->GetBinLowEdge(1);
    double yhi = y->GetBinUpEdge(y->GetNbins());
    double yrlo = 0.501;
    double yrhi = 1.299;

    for (int rep = 1;rep<6;++rep) {
        RooUnfoldBayes*    bayes  = new RooUnfoldBayes(hadd_resp_B, measured, rep);
        auto unf = (TH1D*) bayes->Hreco();
        unf->SetName(tuUniqueName());
        auto _rat = (TH1D*) tuDivide(unf,truth);
        _rat->SetName(tuUniqueName());
        tu_fmt(_rat,{{"MakerStyle",kOpenCircle,"yAxisTitle","Ratio unfolded/truth",
            "xAxisRangeLo",xlo,
            "xAxisRangeHi",xhi,
            "Title",Form("Recoil Side: Range_{truth}#inc[%.0f,%.0f] Range_{meas}#inc[%.0f,%.0f]",ylo,yhi,xlo,xhi),
            "yAxisRangeLo",yrlo, "yAxisRangeHi",yrhi}});
        const char* cmd = (rep == 1 ? "PE1 PLC PMC " : "PE1 PLC PMC same");
        /* _rat->GetXaxis()->SetRangeUser(0.,60.); */
        /* _rat->Scale(1+(double)rep*0.1); */
        _rat->Draw(cmd);
        leg->AddEntry(_rat,Form("rep: %i",rep));
    }
    leg->Draw();
    tuPause();
    pads.pdf("./pdf/$0");

    /* hadd.response->Draw("colz"); */
EOF
