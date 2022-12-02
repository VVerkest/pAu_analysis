
root -l <<EOF
    .L ../loc_lib/loc_libs.h
    .L recon_pt_bound/myDiscreteLinFnc.cc
    .L make_sparse_pp.cc
    Reader R;
    R.Loop()
EOF
