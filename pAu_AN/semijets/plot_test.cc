root -l <<EOF
    .x tu_loadlibs.C
    .L myEnv.h // "myEnvelope"
    .L RooUnfold_RFM.h
    .L unf_arr.h

    // input:
    int    bins        = (strcmp("$1","")==0) ? 0:       atoi("$1");
    int    cuts        = (strcmp("$2","")==0) ? 0:       atoi("$2");
    string zdcx_or_bbc = (strcmp("$3","")==0) ? "zdcx":  "$3";
    string set_lim     = (strcmp("$4","")==0) ? "hi":    "$4";
    string set0        = (strcmp("$5","")==0) ? "hiA":   "$5";
    string set1        = (strcmp("$6","")==0) ? "hiB":   "$6";
    string side        = (strcmp("$7","")==0) ? "0to60": "$7";
    string method      = "$8";

    if (side=="t") {
        side = "0to60";
    } else if (side == "r") { 
        side = "120to180";
    }
    
    gStyle->SetPalette(kCMYK);

    /* tuBinVec bin_M {{ 15., 20., 25., 30., 35., 40., 50 }}; */
    /* tuBinVec bin_T {{ 12, 15., 20., 25., 30., 35., 40., 50., 60. }}; */
    /* string which_bins = "$1"; */
    /* if (which_bins=="") which_bins="0"; */
    /* string which_cuts = "$2"; */
    /* if (which_cuts=="") which_cuts="0"; */

    ostringstream which_file, desc_tag, stamp_tag, print_tag;//, minitag;
    desc_tag << "b"<<bins<<"c"<<cuts<<"_"<<side<<zdcx_or_bbc<<"_"<<set_lim<<set0<<set1<<method;
    which_file << "in-root/hadd" << side << "_" << zdcx_or_bbc << method << ".root";
    /* desc_tag   <<minitag.str()<<"_"<< side<<"_lim"<<set_lim<<"_"<<set0<<set1<<method; */
    stamp_tag << "`date` `pwd`/$0 $* || "<<desc_tag.str();
    print_tag << tuFileName("pdf/$0",{desc_tag.str()});

    array<vector<myEnv>,9> env;
    auto resp_A = hadd_draw_rebin(env, which_file.str(), set0,    desc_tag.str(), stamp_tag.str(), print_tag.str(), cuts, bins);
    auto resp   = hadd_draw_rebin(env, which_file.str(), set_lim, desc_tag.str(), stamp_tag.str(), print_tag.str(), cuts, bins);
    auto resp_B = hadd_draw_rebin(env, which_file.str(), set1,    desc_tag.str(), stamp_tag.str(), print_tag.str(), cuts, bins);

    ostringstream fout_loc; fout_loc<<"out-root/"<<desc_tag.str()<<".root";
    TFile* fout = new TFile( fout_loc.str().c_str(),"recreate");

                       draw_save_resp(resp.first,   resp.second,   set_lim, stamp_tag.str(), print_tag.str());
    auto hadd_resp_A = draw_save_resp(resp_A.first, resp_A.second, set0, stamp_tag.str(), print_tag.str());
    auto hadd_resp_B = draw_save_resp(resp_B.first, resp_B.second, set1, stamp_tag.str(), print_tag.str());

    fout->Save();
    fout->Close();
   
    /* // draw the closure test */
    auto truth    = (TH1D*) hadd_resp_A->Htruth();
    auto measured = (TH1D*) hadd_resp_A->Hmeasured();
    tuPads pads{1};
    pads.prefix="closure_test";
    
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
        /* tu_fmt(_rat,{{"MarkerStyle",kOpenCircle}}); */
    tu_fmt(_rat,{{"MakerStyle",kOpenCircle,"yAxisTitle","Ratio unfolded/truth",
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
    /* tuDrawTLine ( x,yrlo,x,yrhi, {{"LineColor",kRed,"LineStyle",2}}); */
    /* tuDrawTLine ( xlo,1.,xhi,1., {{"LineColor",kRed,"LineStyle",2}}); */
    leg->Draw();
    /* tuDrawTLatex(resp.second->Hresponse()->GetXaxis()->GetTitle(),20.,.75,{{"TextSize",25,"TextColor",kBlue}}); */

    pads.stamp(stamp_tag.str());
    pads.print("pdf/",{desc_tag.str()+"_closure .pdf"});

    /* Form("%s %s %s .pdf","$0 $*", which.c_str(),f_in.c_str())); */
    /* pads.print(Form("%s %s %s .pdf","$0 $*", which.c_str(),f_in.c_str())); */
    /* tuPause(); */
EOF
