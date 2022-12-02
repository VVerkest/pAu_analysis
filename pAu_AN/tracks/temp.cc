
opt << "xAxisTitle" << "This thing is obsurd!" 
    << "MarkerStyle" << kOpenSquare 
    << "MarkerAlpha" << 0.6<<"xAxisRangeLo"<<1.2<<"xAxisRangeHi"<<4.5;

TH1* tu_fmt (TH1* hg, tuOptMap _override, tuOptMap dict) {
    // if don't want defaults, add {} to override value
    dict += _override;

    if (dict["normalize"])    hg->Scale(1./hg->Integral());
    if (dict["Rebin"])        hg->Rebin(dict("Rebin"));

    if (dict["noTitle"])      hg->SetTitle("");
    if (dict["SetStats"])     hg->SetStats((int)dict("SetStats"));

    if (dict["no_errors"]) for (int i=1; i<=hg->GetNbinsX(); ++i) hg->SetBinError(i,0.);

    dict["normalize"] -- check if it has it;
    dict("normalize"); -- defaults to dict("normalize",'d'); // double, string, int
    dict("noramlize",'s');
    dict("normalize",'i');
    dict("normalize",'d');
    dict("axes",'d',0), dict("axes",'d',1) // additional offset
    dict("

