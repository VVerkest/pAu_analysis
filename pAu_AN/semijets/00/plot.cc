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

    auto pa = A_arr.plot();
    pa.pdf("pdf/A$1$2_$3");
    /* tuPause(); */
    pb = B_arr.plot();
    pb.pdf("pdf/B$1$2_$3");
    /* tuPause(); */


    /* A_arr.set_sqrt(); */
    /* B_arr.set_sqrt(); */

    auto _A = A_arr.hadd();
    auto _B = B_arr.hadd();
    //
    // plot the ratios of T, M, miss, fakes
    auto T_A = _A.truth();
    auto T_B = _B.truth();
    auto T_R = tuDivide(T_A, T_B);
    delete T_A, T_B;

    auto M_A = _A.measured();
    auto M_B = _B.measured();
    auto M_R = tuDivide(M_A, M_B);
    delete M_A, M_B;

    auto F_A = _A.fakes;
    auto F_B = _B.fakes;
    auto F_R = tuDivide(F_A, F_B);

    auto I_A = _A.miss;
    auto I_B = _B.miss;
    auto I_R = tuDivide(I_A, I_B);


    tu_fmt(T_R, {{"Title", "Ratio: Truth (A:B)",     "MarkerStyle", kFullCircle,       "MarkerSize", 1.4, "MarkerAlpha", 0.6, "MarkerColor", kRed, 
            "yAxisRangeLo", 0.0, "yAxisRangeHi", 2.0}});
    tu_fmt(M_R, {{"Title", "Ratio: Meausured (A:B)", "MarkerStyle", kFullSquare,       "MarkerSize", 1.4, "MarkerAlpha", 0.6, "MarkerColor", kBlue}});
    tu_fmt(F_R, {{"Title", "Ratio: Fakes (A:B)",     "MarkerStyle", kFullTriangleUp,   "MarkerSize", 1.4, "MarkerAlpha", 0.6, "MarkerColor", kGreen+2}});
    tu_fmt(I_R, {{"Title", "Ratio: Misses (A:B)",    "MarkerStyle", kFullTriangleDown, "MarkerSize", 1.4, "MarkerAlpha", 0.6, "MarkerColor", kBlack}});

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas* canv = new TCanvas("canv","canv", 800, 500 );
    canv->cd();
    canv->SetBottomMargin(0.2);
    T_R->Draw("PE");
    M_R->Draw("PE same");
    F_R->Draw("PE same");
    I_R->Draw("PE same");
    canv->BuildLegend();
    canv->SaveAs("pdf/ratios_$1$2_$3.pdf");
    string test = "$4";
    if (test!="") {
        tuPause();
    }


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
    /* tuPause(); */
    pads.pdf("./pdf/plot_$1$2_$3");

    /* hadd.response->Draw("colz"); */
EOF
