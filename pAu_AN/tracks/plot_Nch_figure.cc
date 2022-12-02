root -l<<EOF
    .x tu_loadlibs.C

    tuPads pads {1,{900, 450}};    
    gStyle->SetPalette(kCMYK)

    tuGetter got{};

    auto mb  = (TH1D*) got("plot_NchMB_fin.root",    "rhoNch");
    auto ht0 = (TH1D*) got("plot_NchHT_fin_4.root",  "rhoNch");
    auto ht1 = (TH1D*) got("plot_NchHT_fin_8.root",  "rhoNch");
    auto ht2 = (TH1D*) got("plot_NchHT_fin_12.root", "rhoNch");

    auto mb_sys  = (TH1D*) got("plot_NchMB_fin.root",    "rhoNch_syserr");
    auto ht0_sys = (TH1D*) got("plot_NchHT_fin_4.root",  "rhoNch_syserr");
    auto ht1_sys = (TH1D*) got("plot_NchHT_fin_8.root",  "rhoNch_syserr");
    auto ht2_sys = (TH1D*) got("plot_NchHT_fin_12.root", "rhoNch_syserr");


    /* mb->Scale(4.*M_PI); */
    /* ht0->Scale(4.*M_PI); */
    /* ht1->Scale(4.*M_PI); */
    /* ht2->Scale(4.*M_PI); */

    vector<int> shapes { kOpenCircle, kOpenTriangleDown, kOpenSquare, kOpenTriangleUp };
    /* vector<int> colors { 1179, 1230, 1281, 1332, 1383, 1433 }; */
    vector<int> colors { 1179, 1383, 1332, 1281, 1332, 1383, 1433 };

    vector<TH1D*> v_hg   { mb, ht0, ht1, ht2 };
    vector<TH1D*> v_err   { mb_sys, ht0_sys, ht1_sys, ht2_sys };
    vector<const char*> tags { "MB", "HT E_{T}^{trig} #in[ 4, 8] GeV/#it{c}",
                                     "HT E_{T}^{trig} #in[ 8,12] GeV/#it{c}",
                                     "HT E_{T}^{trig} #in[12,30] GeV/#it{c}"};

    const double yRangeHi = 1.79;
    auto hg_blank = (TH1D*) mb->Clone("hg_blank");
    hg_blank->Reset();
    tu_fmt( hg_blank, {{ 
        "yAxisTitle",        "#LT#frac{d^{2}#it{N}_{ch}}{d#etad#phi}#GT",
        "yAxisTitleSize",      18,
        "yAxisTitleOffset",  2.40,
        "xAxisTitleOffset",  1.31,

        "MarkerSize",         1.0,
        "Title",               "",
        "xAxisTitle",        "EA_{BBC}",
        "xAxisTitleSize",      18,
        "xAxisTitleOffset",  1.45,

        "MarkerStyle",kOpenCircle,
        "yAxisRangeLo",      0.01,
        "yAxisRangeHi",      yRangeHi}});
    pads(0);

    hg_blank->Draw("PE");
    TLegend *leg = new TLegend(0.6406578,0.3019208,0.9526623,0.6572629,NULL,"brNDC");
    tu_fmt(leg);

    for (int i=0; i<4; ++i) {
        auto hg = v_hg[i];

        const double x_cen = (i+1.5)/6.;
        tuSysErrors pts { hg, {0.,1., x_cen } };
        pts.tgase->SetMarkerStyle(shapes[i]);
        pts.tgase->SetMarkerColorAlpha(colors[i],1.);
        pts.tgase->SetLineColorAlpha(colors[i],1.);
        pts.tgase->Draw("PE1");

        leg->AddEntry(pts.tgase, tags[i]);

        auto h_err = v_err[i];
        tuSysErrors p_err { h_err, {x_cen-1./12., x_cen+1./12., x_cen} };
        /* tuSysErrors p_err { h_err, {0.,1.,x_cen}}; */
        p_err.tgase->SetMarkerStyle(kDot);
        p_err.tgase->SetMarkerColorAlpha(kWhite,0.);
        p_err.tgase->SetFillColorAlpha(colors[i],0.2);
        p_err.tgase->Draw("PE2");
    };
    
    for (int i=1; i<10; ++i) {
        double x = hg_blank->GetXaxis()->GetBinUpEdge(i);
        tuDrawTLine( x, 0.01, x, yRangeHi, {{"LineColor",kGray,"LineStyle",2}});
    }
    leg->Draw();
    pads.stamp("`date` `pwd`/$0 $*");
    tuPause();
    pads.pdf("./pdf/$0");
    



EOF
