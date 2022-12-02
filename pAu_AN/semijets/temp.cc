root -l<<EOF

    .x tu_loadlibs.C
    .L myEnv.h // shouldn't need
    .L RooUnfold_RFM.h
    .L unf_arr.h
    

    .x tu_loadlibs.C
    tuGetter got {};
    gStyle->SetPalette(kCMYK);

    string which = "$2";
    if (which=="") which="0to60";

    auto resp   = (RooUnfoldResponse*) got("hadd_xsec$1.root",Form("response%s_B",which.c_str()));
    auto resp_B = (RooUnfoldResponse*) got("hadd_xsec$1.root",Form("response%s_A",which.c_str()));
    auto truth  = (TH1D*) resp_B->Htruth();
    auto meas   = (TH1D*) resp_B->Hmeasured();
    tuPads pads{1};
    pads(0);
   TLegend *leg = new TLegend(0.0008687949,0.3057933,0.1507785,0.738973,NULL,"brNDC");
    TAxis* x = meas->GetXaxis();
    double xlo = x->GetBinLowEdge(1);
    double xhi = x->GetBinUpEdge(x->GetNbins());
    TAxis* y = truth->GetXaxis();
    double ylo = y->GetBinLowEdge(1);
    double yhi = y->GetBinUpEdge(y->GetNbins());

    double yrlo = 0.801;
    double yrhi = 1.299;
    
    for (int rep = 1;rep<6;++rep) {
        RooUnfoldBayes*    bayes  = new RooUnfoldBayes(resp, meas, rep);
        auto unf = (TH1D*) bayes->Hreco();
        unf->SetName(tuUniqueName());
        auto _rat = (TH1D*) tuDivide(unf,truth);
        _rat->SetName(tuUniqueName());
        /* tu_fmt(_rat,{{"MarkerStyle",kOpenCircle}}); */
    tu_fmt(_rat,{{"MakerStyle",kOpenCircle,"yAxisTitle","Ratio unfolded/truth",
            "xAxisTitle","Jet pt truth",
            "xAxisRangeLo",0.,
            "xAxisRangeHi",60.,
            "Title",Form("Range_{truth}#inc[%.0f,%.0f] Range_{meas}#inc[%.0f,%.0f]",ylo,yhi,xlo,xhi),
            "yAxisRangeLo",yrlo, "yAxisRangeHi",yrhi}});
        const char* cmd = (rep == 1 ? "PE1 PLC PMC " : "PE1 PLC PMC same");
        _rat->GetXaxis()->SetRangeUser(0.,60.);
        _rat->Draw(cmd);
        leg->AddEntry(_rat,Form("rep: %i",rep));
    }
    double x= 15;
    tuDrawTLine ( x,yrlo,x,yrhi, {{"LineColor",kRed,"LineStyle",2}});
    x = 40.;
    tuDrawTLine ( x,yrlo,x,yrhi, {{"LineColor",kRed,"LineStyle",2}});
    tuDrawTLine ( 0.,1.,60.,1., {{"LineColor",kRed,"LineStyle",2}});
    leg->Draw();

    pads.stamp("`date` `pwd`/$0 $*");
    pads.print("$0 $* .eps")

EOF
