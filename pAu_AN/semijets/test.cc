root -l <<EOF
    .x tu_loadlibs.C
    .L RooUnfold_RFM.h

    // ok let's check closure
    auto _data = RooUnfold_RFM_Array("hadd_jetemb.root", "response0to60");
    _data.scrub(10);

    auto hadd = _data.hadd();
    tuPads pads{1};
    pads(0)->SetLogz();
    hadd.response->Draw("colz");
    tuPause();
EOF
