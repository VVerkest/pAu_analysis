root -l<<EOF
    .x tu_loadlibs.C
    .L ../common/pAu_consts.h

    /* ioPads pads { {{0,0.21,0.99}},{{0,0.10,0.4,.5},{0.5,0.6,0.90,1.}},905,474}; */
    tuPads pads { 1, {{kLeft, 0.15, kRight, 0.2, 800, 700}} };
    /* gStyle->SetPalette(kCMYK); */

    pads(0)->SetLogz();

    tuGetter got{};
    
    /* const char* bin_dat = "/home/dsjohnny/root_macros/io_lib/pAu2015_common/bin_edges.txt"; */
    // get statistics

    auto hg = (TH2D*) got("hadd_UE_luminosity.root","Trans_h2UE");

    // rotate the data

    tu_fmt(hg, {{"yAxisTitle","EA_{UE}",
            "yAxisTitleOffset",1.4,
            "xAxisTitleOffset",1.3,
            "yAxisTitleSize",27,
            "xAxisTitleSize",27,
            "yAxisLabelSize",25,
            "xAxisLabelSize",25,
            "Title","","xAxisTitle","EA_{BBC}"}});
    auto px = (TProfile*) got("hadd_UE_luminosity.root","EA_tpc_vs_bbc");
    auto py = (TProfile*) got("hadd_UE_luminosity.root","EA_bbc_vs_tpc");

    auto points_y = tuMakeTGraphErrors(py->ProjectionX(),true);

    tuOptMap fmt_bbc {{"MarkerStyle",kFullSquare,"MarkerColor",kRed,
                       "MarkerSize",0.9,"MarkerAlpha",0.8}};
    tuOptMap fmt_tpc {{"MarkerStyle",kFullCircle,"MarkerColor",kBlack,
                       "MarkerSize",0.9,"MarkerAlpha",0.8}};

    tu_fmt(points_y,fmt_bbc);

    hg->Draw("colz");
    points_y->Draw("PE");
    auto px = hg->ProfileX("px");
    tu_fmt(px,fmt_tpc);
    px->Draw("PE same");

   TLegend *leg = new TLegend(0.5613779,0.9028723,0.7800368,0.9950869,NULL,"brNDC");
    leg->SetFillColorAlpha(kWhite,0.8);
    leg->SetLineColorAlpha(kWhite,0.);
    leg->AddEntry(points_y,"#LTEA_{BBC}#GT(EA_{UE})");
    leg->AddEntry(px,      "#LTEA_{UE}#GT(EA_{BBC})");
    leg->Draw();
    /* tuDrawTLatex("MB Events",4.15,62000.,{{"TextSize",13}}); */

    tuBinVec bin_bbc = bins_EAbbc3;
    tuBinVec bin_tpc = bins_UEtpc3;


    /* bin_bbc.vec[3] = 70000; */
    /* bin_tpc.vec[3] = 1.5.; */

    /* array<double,3> loc_x { (bin_tpc[0]+bin_tpc[1])*0.5, */ 
                             /* (bin_tpc[1]+bin_tpc[2])*0.5, */ 
                             /* 2.0 }; */
    double loc_y { (bin_tpc[0]+bin_tpc[1])*0.5 };
    array<double,3> loc_x { (bin_bbc[0]+bin_bbc[1])*0.5, 
                             (bin_bbc[1]+bin_bbc[2])*0.5, 
                             35000. };

    for (int i{1};i<3;++i){
        tuDrawTLine(bin_bbc[i],0.,bin_bbc[i],bin_tpc[3],{{"LineColor",kBlack,"LineStyle",2}});
    }
    /* for (int i{1};i<3;++i){ */
    /*     tuDrawTLine(bin_tpc[i],0.,bin_tpc[i],bin_bbc[3],{{"LineColor",kBlack,"LineStyle",2}}); */
    /* } */
    /* for (int i_tpc{0};i_tpc<3;++i_tpc) */
    /*     tuDrawTLineHorizontal(bin_tpc[i_tpc+1],{{"LineColor",kBlack,"LineStyle",2}}); */

    auto cnt = (TH2D*) got("hadd_UE_luminosity.root","h2UE_9grid");
    cnt->Scale(1./cnt->Integral());
    for (int i_bbc{0};i_bbc<3;++i_bbc){
        double val = 0;
        for (int i_tpc{0};i_tpc<3;++i_tpc) val += cnt->GetBinContent(i_tpc+1,i_bbc+1);
        tuDrawTLatex(Form("%2.1f%%",100*val),
                    loc_x[i_bbc],loc_y, {{"TextColor",(kGray+2+0),"ColorAlpha",0.7,
                    "TextSize",22,"TextAngle",0,"TextAlign",22}});
    }

    /* hg_old = (TH2D*) got("prelim","mb_rho_bbc_All"); */
    /* hg_old->Rebin2D(2,2); */
    /* auto py_old = hg_old->ProfileY(tuUniqueName()); */
    /* auto points_old = tuMakeTGraphErrors(py_old->ProjectionX(),true); */
    /* tu_fmt(points_old,{{"MarkerStyle",kOpenSquare,"MarkerColor",kCyan}}); */
    /* points_old->Draw("PE"); */
    /* auto px_old = hg_old->ProfileX(tuUniqueName()); */
    /* tu_fmt(px_old,{{"MarkerStyle",kOpenCircle,"MarkerColor",kMagenta}}); */
    /* px_old->Draw("PE same"); */

    pads.stamp("`date` `pwd`/$0 $*");
    pads.save("$0");

    tuPause();

EOF
