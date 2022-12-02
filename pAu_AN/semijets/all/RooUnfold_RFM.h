// note requires that .x tu_loadlibs.C has already been run
/* #ifndef RooUnfold_RFM */
/* #define RooUnfold_RFM */
struct RooUnfold_RFM { // RFM for response fakes miss
    TH2D* response;
    TH1D* miss;
    TH1D* fakes;
    TH1D* match_T() { 
        tuInflate(response);
        return (TH1D*) response->ProjectionY(tuUniqueName());
    };
    TH1D* match_M() { 
        tuInflate(response);
        return (TH1D*) response->ProjectionX(tuUniqueName());  
    };
    TH1D* truth() {
        auto _truth = match_T();
        /* cout << " _truth " << _truth->GetXaxis()->GetNbins() << " miss: " << miss->GetXaxis()->GetNbins() << endl; // lion */
        _truth->Add(miss,1.);
        return _truth;
    };
    TH1D* measured() {
        auto _measured = match_M();
        _measured->Add(fakes,1.);
        /* cout << " _measured " << _measured->GetXaxis()->GetNbins() << " fakes: " << fakes->GetXaxis()->GetNbins() << endl; // lion */
        return _measured;
    };
    RooUnfold_RFM(RooUnfoldResponse* ruu) {
        response = (TH2D*) ruu->Hresponse();
        TH1D* _;

        miss = (TH1D*) ruu->Htruth();
        _ = match_T();
        miss->Add(_,-1.);
        delete _;

        fakes = (TH1D*) ruu->Hmeasured();
        _ = match_M();
        fakes->Add(_,-1.);
        delete _;
    };
    RooUnfold_RFM() : response{nullptr}, miss{nullptr}, fakes{nullptr} {};

    RooUnfoldResponse* operator()(string name="") {
    /* operator RooUnfoldResponse*(string name) { */
        auto _truth    = truth();
        auto _measured = measured();
        auto ruu = new RooUnfoldResponse(_measured, _truth, response, name=="" ? tuUniqueName() : name.c_str());
        delete _truth;
        delete _measured;
        return ruu;
    }; // caste to a RooUnfoldResponse& from response, miss and fakes
    RooUnfold_RFM& scale(double val) {
        miss->Scale(val);
        fakes->Scale(val);
        response->Scale(val);
        return *this;
    };
    RooUnfold_RFM clone(double weight=1.) {
        RooUnfold_RFM _clone{};
        _clone.response = (TH2D*) response->Clone(tuUniqueName());
        _clone.fakes    = (TH1D*) fakes->Clone(tuUniqueName());
        _clone.miss     = (TH1D*) miss->Clone(tuUniqueName());
        if (weight != 1.) _clone.scale(weight);
        return _clone;
    };
    void add(RooUnfold_RFM& ruu, double weight=1.) {
        miss->Add(ruu.miss,weight);
        response->Add(ruu.response,weight);
        fakes->Add(ruu.fakes,weight);
    };
    void set_sqrt() {
        tuSqrtErr(miss,true);
        tuSqrtErr(fakes,true);
        tuSqrtErr(response,true);
    };
    RooUnfold_RFM rebin(vector<double>& bin_M, vector<double> bin_T, bool set_sqrt=true) {
        RooUnfold_RFM ruu {};
        ruu.fakes = tuNaiveRebin1D(fakes, bin_M);
        /* cout << " rebin fakes: " << ruu.fakes->GetXaxis()->GetNbins() << endl; // lion */
        ruu.miss  = tuNaiveRebin1D(miss,  bin_T);
        /* cout << " rebin miss: : " << ruu.miss->GetXaxis()->GetNbins() << endl; // lion */
        ruu.response =  tuNaiveRebin2D(response, bin_M, bin_T);
        /* cout << " res Y: " << ruu.response->GetYaxis()->GetNbins() << " X: " << ruu.response->GetXaxis()->GetNbins() << endl; // lion */
        if (set_sqrt) ruu.set_sqrt();
        return ruu;
    }
};

struct RooUnfold_RFM_Array {

    array<RooUnfold_RFM,9> data{};
    array<vector<int>,9> vec_kleft{}, vec_kright{};
    array<int,9> fake_left{}, fake_right{}, miss_left{}, miss_right{};
    RooUnfold_RFM_Array(){};
    RooUnfold_RFM_Array (string file, string tag) {
        tuGetter got{};
        for (int i=0;i<9;++i) {
            auto _got = (RooUnfoldResponse*) got(file.c_str(),Form("%s_%i",tag.c_str(),i));
            data[i] = RooUnfold_RFM(_got);
        }
    };
    void print() {
        for (int i=0;i<9;++i) {
            cout << " data["<<i<<"].response -> " << data[i].response->GetEntries() << endl;
            cout << " data["<<i<<"].miss     -> " << data[i].miss->GetEntries() << endl;
            cout << " data["<<i<<"].fakes    -> " << data[i].fakes->GetEntries() << endl;
        }
    };
    RooUnfold_RFM hadd(bool app_weight=true) {
        tuXsec_pAu2015 Xsec{};
        RooUnfold_RFM _hadd;
        for (int i=0;i<9;++i) {
            double weight = app_weight ? Xsec.XsecPyth6(i) : 1.;
            if (i==0) {
                _hadd = data[i].clone(weight);
            } else {
                _hadd.add(data[i],weight);
            }
        }
        return _hadd;
    };
    void scrub( int nbins ) {
        for (int i=0;i<9;++i) {
            tuScrubBins(data[i].response, nbins);
            tuScrubBins(data[i].fakes,    nbins);
            tuScrubBins(data[i].miss ,    nbins);
        }
    };
    void scrub_left_quant( double quant ) {
        for (int i=0;i<9;++i) {
            vec_kleft[i] = tuVecScrubQuant(data[i].response, quant, kLeft, true);
            fake_left[i] = tuScrubQuant(data[i].fakes,    quant, kLeft, true);
            miss_left[i] = tuScrubQuant(data[i].miss ,    quant, kLeft, true);
        }
    };
    void scrub_right_quant( double quant ) {
        for (int i=0;i<9;++i) {
            vec_kright[i] = tuVecScrubQuant(data[i].response, quant, kRight, true);
            fake_right[i] = tuScrubQuant(data[i].fakes,    quant, kRight, true);
            miss_right[i] = tuScrubQuant(data[i].miss ,    quant, kRight, true);
        }
    };
    void scrub_right_nsig( double nsig, double q0, double q1 ) {
        for (int i=0;i<9;++i) {
            vec_kright[i] = tuVecScrubNsig(data[i].response, nsig, q0, q1, kRight, true);
            fake_right[i] = tuVecScrubNsig(data[i].fakes,    nsig, q0, q1, kRight, true);
            miss_right[i] = tuVecScrubNsig(data[i].miss ,    nsig, q0, q1, kRight, true);
        }
    };
    void copy_scrub(RooUnfold_RFM_Array* _in) {
        for (int i=0;i<9;++i) {
            if (_in->vec_kleft[i].size()==0) return;
            tuVecScrub( data[i].response, _in->vec_kleft[i],  kLeft);
            tuVecScrub( data[i].response, _in->vec_kright[i], kRight);

            tuScrubBlock( data[i].fakes, _in->fake_left[i],  _in->fake_right[i]);
            tuScrubBlock( data[i].miss,  _in->miss_left[i],  _in->miss_right[i]);
        }
    };
    RooUnfold_RFM_Array rebin(vector<double> bin_M, vector<double> bin_T, bool set_sqrt=true) {
        RooUnfold_RFM_Array r_val{};
        for (int i=0;i<9;++i) {
            r_val.data[i] = data[i].rebin(bin_M, bin_T, set_sqrt);
        }
        return r_val;
    };
    RooUnfold_RFM& operator[](int i) { return data[i]; };
    tuPads plot(tuOptMap opts={}) {
        tuOptMap plot_opt {{ "Title","",
                    "xAxisTitle","#it{p}_{T,full-jet}^{rec.} [GeV/#it{c}]",
                    "yAxisTitle","#it{p}_{T,full-jet}^{emb.} [GeV/#it{c}]",
                    "xAxisTitleSize",15,"yAxisTitleSize",15,
                    "yAxisLabelSize",15, "yAxisTitleOffset",1.8,
                    "xAxisLabelSize",15, "xAxisTitleOffset",1.2,
                    /* "zAxisNdivisions",6, */
                    /* "xAxisNdivisions",6, */
                    "yAxisNdivisions",6,
                    "xAxisRangeLo",1.,"xAxisRangeHi",59.,
                    "yAxisRangeLo",1.,"yAxisRangeHi",59.,
        }};
        plot_opt += opts;
        tuPads pads  {3,{600,400,kRight,0.2}, 3 };
        auto pthat_bins {vector<int>{5,7,9,11,15,25,35,45,55,65}};
        for (int i=0;i<9;++i) { 
            pads(i)->SetLogz();//->SetLogz();
            auto& h2 = data[i].response;
            tu_fmt(h2,plot_opt);
            h2->Draw("colz");
            tuDrawTLatex(Form("#hat{#it{p}}_{T} #in [%i,%i]",pthat_bins[i],pthat_bins[i+1]),10.,48.,{{"TextSize",12}});
        }
        return pads;
    }
};

/* #endif */

