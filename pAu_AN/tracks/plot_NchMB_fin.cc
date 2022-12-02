root -l <<EOF
    .L ../loc_lib/loc_libs.h

    noiGetter got{};

    auto hg_dAu = (TH1D*) got("calc_NchMB_perEA_dAu.root", "hg_MB_cor");
    auto hg_pp  = (TH1D*) got("calc_NchMB_perEA_pp.root",  "hg_MB_cor");


    TFile PU_file { "plot_PU_dAu.root","read" };
    array<double[2],4> *PU_fit_values, *PU_fit_errors;
    PU_file.GetObject( "PU_fit_values", PU_fit_values );
    PU_file.GetObject( "PU_fit_errors", PU_fit_errors );

    TFile dat_file { "calc_NchMB_perEA_dAu.root", "read" };
    array<double,10> *sumZDCx, *sumNevents;
    dat_file.GetObject( "sumZDCx", sumZDCx );
    dat_file.GetObject( "sumNevents", sumNevents );


    // where to save output files
    TFile f_save ((noiStripExtension("$0")+".root").c_str(),"recreate");

    // calculate the PU errors
    auto hg_PU_nom    = (TH1D*) hg_dAu->Clone("PU_nom");
    auto hg_PU_zero   = (TH1D*) hg_dAu->Clone("PU_zero");
    auto hg_PU_m_max  = (TH1D*) hg_dAu->Clone("PU_m_max");

    double m_etaAu  = (*PU_fit_values)[_etaAu][_PUfit_m];
    double m_etaPP  = (*PU_fit_values)[_etaPP][_PUfit_m];
    double m_etaMid = (*PU_fit_values)[_etaPP][_PUfit_m];

    double merr_etaAu  = (*PU_fit_errors)[_etaAu][_PUfit_m];
    double merr_etaPP  = (*PU_fit_errors)[_etaPP][_PUfit_m];
    double merr_etaMid = (*PU_fit_errors)[_etaPP][_PUfit_m];

    for (int iEA{0}; iEA<10; ++iEA) {
        double nEvents = (*sumNevents)[iEA];

        double sum_zdcx   = (*sumZDCx)[iEA]-2000*nEvents;
        double rho_PU_Au  = sum_zdcx * m_etaAu  / nEvents;
        double rho_PU_PP  = sum_zdcx * m_etaPP  / nEvents;
        double rho_PU_Mid = sum_zdcx * m_etaMid / nEvents;
        double rho_PU = (rho_PU_Au*0.7 + rho_PU_PP*0.6 + rho_PU_Mid*0.6)/2.;
        hg_PU_nom->SetBinContent(iEA+1, hg_PU_nom->GetBinContent(iEA+1) - rho_PU);

        // calculate PU from 0 intead of 2 kHz
        double sum_zdcx_zero = (*sumZDCx)[iEA];
        rho_PU_Au  = sum_zdcx_zero * m_etaAu  / nEvents;
        rho_PU_PP  = sum_zdcx_zero * m_etaPP  / nEvents;
        rho_PU_Mid = sum_zdcx_zero * m_etaMid / nEvents;
        double rho_PU_zero = (rho_PU_Au*0.7 + rho_PU_PP*0.6 + rho_PU_Mid*0.6)/2.;
        hg_PU_zero->SetBinContent(iEA+1, hg_PU_zero->GetBinContent(iEA+1) - rho_PU_zero);

        // calculate PU with m + 1 sigma
        rho_PU_Au  = sum_zdcx * (m_etaAu  + merr_etaAu)  / nEvents;
        rho_PU_PP  = sum_zdcx * (m_etaPP  + merr_etaPP)  / nEvents;
        rho_PU_Mid = sum_zdcx * (m_etaMid + merr_etaMid) / nEvents;
        double rho_PU_errmax = (rho_PU_Au*0.7 + rho_PU_PP*0.6 + rho_PU_Mid*0.6)/2.;
        hg_PU_m_max->SetBinContent(iEA+1, hg_PU_zero->GetBinContent(iEA+1) - rho_PU_errmax);
    }

    TH1D* syserr_ppdAu_prior = hg_abs_deltarat (hg_dAu, hg_pp); // about 0.7% error for pp vs dAu priors
    TH1D* syserr_unfmethod   = hg_const        (hg_dAu,0.005,"UnfMethod");         // 0.5% for using single bin vs bayesian unfolding error
    TH1D* syserr_m_fit       = hg_abs_deltarat (hg_PU_nom, hg_PU_m_max);
    TH1D* syserr_zero        = hg_abs_deltarat (hg_PU_nom, hg_PU_zero );
    TH1D* syserr_const_tr    = hg_const        (hg_dAu, 0.045,"other");

    auto syserr_RSS        = hg_RSS( { syserr_ppdAu_prior, syserr_unfmethod, syserr_m_fit, syserr_zero, syserr_const_tr}, "syserr_RSS");
    
    // draw the data
    noiPads pads{2};
    pads(0);
    noi_fmt(hg_PU_nom, {{"yAxisTitle", "#frac{d^{2}N_{ch}}{d#etad#phi}", "MarkerStyle", kFullCircle, "MarkerColor", kBlue+2, "yAxisRangeLo",0.6, "yAxisRangeHi",2.2}});
    hg_PU_nom->Draw("PE");
    auto hg_boxerrs = hg_syserr(hg_PU_nom, syserr_RSS, "SystematicErrors");
    noi_fmt(hg_boxerrs, {{"MarkerStyle", kDot, "MarkerColor", kWhite, "MarkerAlpha", 0., "FillColor", kBlue+2, "FillAlpha", 0.2}});
    hg_boxerrs->Draw("PE2 same");

    gStyle->SetPalette(kCMYK);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    vector<int> colors { 1179, 1242, 1305, 1368, 1431 };
    noi_fmt(syserr_ppdAu_prior, {{"Title","Prior: pp vs dAu","MarkerColor",    colors[0], "MarkerStyle", kOpenCircle,"no_errors","",
            "xAxisTitle","EA_{BBC}","yAxisTitle","Systematic Error Ratio", "yAxisRangeLo",0.001, "yAxisRangeHi", 0.1205}});
    noi_fmt(syserr_unfmethod,   {{"Title","Unfold Method","MarkerColor",       colors[1], "MarkerStyle", kOpenTriangleUp, "no_errors",""}});
    noi_fmt(syserr_m_fit,       {{"Title","PU to ZDCx fit","MarkerColor",      colors[2], "MarkerStyle", kOpenTriangleDown, "no_errors",""}});
    noi_fmt(syserr_zero,        {{"Title","PU to zero","MarkerColor",          colors[3], "MarkerStyle", kOpenSquare, "no_errors",""}});
    noi_fmt(syserr_const_tr,    {{"Title","STAR tracking unc.", "MarkerColor", colors[4], "MarkerStyle", kOpenDiamond,"no_errors", ""}});
    noi_fmt(syserr_RSS,    {{"Title","Overall Systematic Error", "MarkerColor", kRed+1, "MarkerStyle", kFullCircle,"no_errors", ""}});


    // plot results and error
    pads(1);
    syserr_ppdAu_prior->Draw("PE");
    syserr_m_fit      ->Draw("PE same");
    syserr_const_tr   ->Draw("PE same");
    syserr_unfmethod  ->Draw("PE same");
    syserr_zero       ->Draw("PE same");
    syserr_RSS        ->Draw("PE same");
    pads(1)->BuildLegend();


    pads.stamp("`date` `pwd`/$0 $*");
    noiPause();
    pads.pdf("./pdf/$0");

    /* hg_PU_nom         ->SetMarkerColorAlpha(kBlack,1.); */
    /* hg_boxerrs        ->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_ppdAu_prior->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_m_fit      ->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_const_tr   ->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_unfmethod  ->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_zero       ->SetMarkerColorAlpha(kBlack,1.); */
    /* syserr_RSS        ->SetMarkerColorAlpha(kBlack,1.); */

    hg_PU_nom->SetMarkerColor(kBlack);
    hg_boxerrs->SetMarkerColor(kBlack);
    syserr_ppdAu_prior->SetMarkerColor(kBlack);
    syserr_m_fit      ->SetMarkerColor(kBlack);
    syserr_const_tr   ->SetMarkerColor(kBlack);
    syserr_unfmethod  ->SetMarkerColor(kBlack);
    syserr_zero       ->SetMarkerColor(kBlack);
    syserr_RSS        ->SetMarkerColor(kBlack);

    f_save.WriteObject(hg_PU_nom,         "rhoNch");
    f_save.WriteObject(hg_boxerrs,        "rhoNch_syserr");
    f_save.WriteObject(syserr_ppdAu_prior,"syserr_ppdAu_prior");
    f_save.WriteObject(syserr_m_fit      ,"syserr_m_fit");
    f_save.WriteObject(syserr_const_tr   ,"syserr_const_tr");
    f_save.WriteObject(syserr_unfmethod  ,"syserr_unfmethod");
    f_save.WriteObject(syserr_zero       ,"syserr_zero");
    f_save.WriteObject(syserr_RSS        ,"syserr_RSS");

EOF
