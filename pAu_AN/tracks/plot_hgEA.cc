root -l <<EOF
    .x tu_loadlibs.C
    .L loc_libs.h


    tuGetter got{};
    auto hg = (TH2D*) got("hadd_sparse_wsu.root","h2UE");

    tu_fmt(hg, {{"yAxisTitle","EA_{UE}",
            "yAxisTitleOffset",1.4,
            "xAxisTitleOffset",1.3,
            "yAxisTitleSize",27,
            "xAxisTitleSize",27,
            "yAxisLabelSize",25,
            "xAxisLabelSize",25,
            "Title","","xAxisTitle","EA_{BBC}"}});
    auto px = (TProfile*) got("hadd_sparse_wsu.root","EA_tpc_vs_bbc");
    auto py = (TProfile*) got("hadd_sparse_wsu.root","EA_bbc_vs_tpc");

    auto points_y = tuMakeTGraphErrors(py->ProjectionX(),true);

    tuOptMap fmt_bbc {{"MarkerStyle",kFullSquare,"MarkerColor",kRed,
                       "MarkerSize",0.9,"MarkerAlpha",0.8}};
    tuOptMap fmt_tpc {{"MarkerStyle",kFullCircle,"MarkerColor",kBlack,
                       "MarkerSize",0.9,"MarkerAlpha",0.8}};

    tu_fmt(points_y,fmt_bbc);

    tuPads pads { 1, {{kLeft, 0.15, kRight, 0.2, 800, 700}} };
    pads(0)->SetLogz();

    hg->Draw("colz");
    points_y->Draw("PE");
    /* auto px = hg->ProfileX("px"); */
    tu_fmt(px,fmt_tpc);
    px->Draw("PE same");

    TLegend *leg = new TLegend(0.5613779,0.9028723,0.7800368,0.9950869,NULL,"brNDC");
    leg->SetFillColorAlpha(kWhite,0.8);
    leg->SetLineColorAlpha(kWhite,0.);
    leg->AddEntry(points_y,"#LTEA_{BBC}#GT(EA_{UE})");
    leg->AddEntry(px,      "#LTEA_{UE}#GT(EA_{BBC})");
    leg->Draw();

    /* tuBinVec bin_bbc = bins_EAbbc3070; */
    auto cnt = (TH2D*) got("hadd_sparse_wsu.root","hg_bbc3070");
    auto axis = cnt->GetXaxis();
    array<double,3> loc_x { axis->GetBinCenter(1),
                            axis->GetBinCenter(2),
                             35000. };

    const float y_hi = (78+0.5)/(4.0*M_PI)
    for (int i{1};i<3;++i){
        tuDrawTLine(axis->GetBinUpEdge(i),0.,axis->GetBinUpEdge(i),y_hi,{{"LineColor",kBlack,"LineStyle",2}});
    }

    cnt->Scale(1./cnt->Integral());
    const double loc_y = 1.25;
    for (int i_bbc{0};i_bbc<3;++i_bbc){
        tuDrawTLatex(Form("%2.1f%%",100.*cnt->GetBinContent(i_bbc+1)),
                    loc_x[i_bbc],loc_y, {{"TextColor",(kGray+2+0),"ColorAlpha",0.7,
                    "TextSize",22,"TextAngle",0,"TextAlign",22}});
    }

    pads.stamp("`date` `pwd`/$0 $*");
    tuPause();
    pads.pdf("pdf/$0");
EOF
