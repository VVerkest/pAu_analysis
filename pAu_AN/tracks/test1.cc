root -l << EOF
    .L ../loc_lib/noiPads.h
    .L ../loc_lib/noiDict.h
    .L ../loc_lib/noi_fmt.h
    .L ../loc_lib/noi_fnc.h

    noiPads pads{3};

    TH1D hg {"hg","a;b;c",10,0.,10.};
    TRandom3 _rand;
    for (int i=0;i<100;++i) {
        hg.Fill(_rand.Uniform(0.,10.));
    }
    noi_fmt(&hg,{{"MarkerStyle",kFullSquare,"LineColorAlpha",kBlue,1.,"FillColorAlpha",kBlue,0.9,"FillColor",kGreen+2}});
    /* hg.SetMarkerColorAlpha(kRed,0.3); */
    pads(0);
    hg.Draw("PE1");
    pads(1);
    hg.Draw("PE2");
    pads(2);
    hg.Draw("PE3");

    pads.stamp("`date` `pwd`/$0 $*");
    pads.pdf("./pdf/$0");

    noiPause();

EOF
