root -l <<EOF
    // Draw the profile pT per ZDCx
    .L ../loc_lib/loc_libs.h
    /* .L ../loc_lib/pAu_bins.h */
    /* .L ../loc_lib/noiGetter.h */
    /* .L ../loc_lib/noiPads.h */
    /* .L ../loc_lib/noi_fmt.h */
    /* .L ../loc_lib/noi_fnc.h */

    noiGetter got{};
    gStyle->SetPalette(kCMYK);

    vector<int> colors { 1179, 1264, 1349, 1433 }

    noiPads pads {1,{kLeft,0.15,kRight,0.15}};
    pads(0);

    TLegend *leg = new TLegend(0.5894348,0.7323897,0.780229,0.9483213,NULL,"brNDC");
    auto events = (THnSparseD*) got("hadd_sparse_wsu.root","sp_events");
    auto mb = events->Projection(_zdcx);
    mb->SetName("mb");
    noi_fmt(mb,{{"MarkerStyle",kFullCircle,"MarkerColor",colors[0],"yAxisTitle","N_{events}, arb. scaling"}});
    mb->SetFillColorAlpha(colors[0],0.6);
    mb->SetFillStyle(3325);
            /* "yAxisTitle","#frac{dN_{ch}}{dN_{event}} #frac{1}{4#pi}","xAxisTitle","ZDCx", */
            /* "yAxisRangeLo",-0.09,"yAxisRangeHi",1.99,"Title",""}}); */
    mb->Scale(1./mb->GetMaximum());
    mb->SetTitle("");
    mb->Draw("hist F");
    leg->AddEntry(mb,"Min Bias.");

    auto events_hi = (THnSparseD*) got("hadd_sparse_wsu_HT.root","sp_events");
    set_range(_trigEt, 2, 4, {events_hi});
    auto ht = events_hi->Projection(_zdcx);
    ht->SetName("ht");
    noi_fmt(ht,{{"MarkerStyle",kFullSquare,"MarkerColor",colors[1]}});
            /* "yAxisTitle","#frac{dN_{ch}}{dN_{event}} #frac{1}{4#pi}","xAxisTitle","ZDCx", */
            /* "yAxisRangeLo",-0.09,"yAxisRangeHi",1.99,"Title",""}}); */
    ht->Scale(1./ht->GetMaximum());
    ht->SetFillColorAlpha(colors[1],0.6);
    ht->SetFillStyle(3354);
    ht->Draw("hist same");
    leg->AddEntry(ht,"High Tower (w/ tower #GE 4 GeV)");

    auto emb_sp = (THnSparseD*) got("make_sparse_dAu.root","sp_truth");
    auto emb = emb_sp->Projection(_zdcx);
    emb->SetName("emb");
    noi_fmt(emb,{{"MarkerStyle",kFullTriangleUp,"MarkerColor",colors[2]}});
            /* "yAxisTitle","#frac{dN_{ch}}{dN_{event}} #frac{1}{4#pi}","xAxisTitle","ZDCx", */
            /* "yAxisRangeLo",-0.09,"yAxisRangeHi",1.99,"Title",""}}); */
    emb->Scale(1./emb->GetMaximum());
    emb->SetFillColorAlpha(colors[2],0.2);
    emb->SetFillStyle(3001);
    emb->Draw("hist same");
    leg->AddEntry(emb,"Embedded Tracks (st_physics)");

    auto emb_jets = (TH1D*) got("jet_emb_ssd_zdcx.root","zdcx");
    noi_fmt(emb_jets, {{"MarkerColor",colors[3],"FillColorAlpha",colors[3],0.2,"FillStyle",3001}});
    emb_jets->Scale(1./emb_jets->GetMaximum());
    emb_jets->Draw("hist same");
    leg->AddEntry(emb_jets,"Embedded jets (st_ssd)");
    leg->Draw();

    pads.stamp("`date` `pwd`/$0 $*");
    pads.pdf("./pdf/$0");
    noiPause();

EOF
