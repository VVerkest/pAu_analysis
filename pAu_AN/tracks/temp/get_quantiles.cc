root -l hadd_sparse_eta.root<<EOF
    .x tu_loadlibs.C
    vector<double> rats {{ 0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.}};
    auto loc = tuQuantiles(hg_bbc, rats);
    for (int i=0;i<11;++i) {
        cout << " i(" << i << ",  " << rats[i] << ", " << loc[i] << endl;
    }
    cout << " static ioBinVec bin_bbc {{";
    for (auto val : loc) cout << ", " << val;
    cout << "}};" << endl;
    

EOF
