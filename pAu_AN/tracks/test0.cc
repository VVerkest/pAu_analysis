root -l << EOF
    .x tu_loadlibs.C
    .L loc_lib.h

    TFile fout {"plot_PU_dAu.root","read"};
    array<double[2],4>* pt, *pt_err;
    fout.GetObject("PU_fit_values", pt);
    fout.GetObject("PU_fit_errors", pt_err);
    cout << " this m: " << (*pt)[_etaAu][_PUfit_m]  << "+/- " << (*pt_err)[_etaAu][_PUfit_m]  <<endl;
    cout << " this m: " << (*pt)[_etaMid][_PUfit_m] << "+/- " << (*pt_err)[_etaMid][_PUfit_m]  <<endl;
    cout << " this m: " << (*pt)[_etaPP][_PUfit_m]  << "+/- " << (*pt_err)[_etaPP][_PUfit_m]  <<endl;

EOF
