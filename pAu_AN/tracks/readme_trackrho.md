# To generate supporting (and final figures) for track density and EA



------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------
### Figure Nch for the paper
------------------------------------------------------------------------------------------------------------------
Input files: common input files are located in `../loc_lib/`
    They are:
 * `loc_libs.h`  : a header file that includes all the other header files
 * `pAu_bins.h`  : binning for the analysis
 * `noiBinVec.h` : conveince struct for binning
 * `noiDict.h`   : conveince struct to pass options around (used in `noi{Pads,_fnc}`)
 * `noiPads.h`   : conveince struct for making TPads on TCanvas
 * `noi_fmt.h`   : conveince struct for formatting options in drawing things
 * `noi_fnc.h`   : conveince struct for little functions
 * `noiGetter.h` : conveince struct for getting objects from TFiles

 Files to pull from WSU and rcas:
`hadd_sparse_wsu.root`      : `scp wsu:/tier2/home/groups/rhi/vverkest/pAu_GraceAnalysis/out-data/sparse_eta/hadd.root hadd_sparse_wsu.root`
`hadd_sparse_wsu_HT.root`   : `scp wsu:/tier2/home/groups/rhi/vverkest/pAu_GraceAnalysis/out-data/sparse_eta_HT.root hadd_sparse_wsu_HT.root`
`hadd_piKp_track_emb.root`  : rcas at `/gpfs/mnt/gpfs01/star/pwg/dstewart/outroot/AN/track-embedding-SL16d_embed/v0/hadd_piKp_track_emb.root`
Also a zdcx TH1D from the jet embedding

user runs:               comments(//), inputs(<-), output (->)
---------------------    -----------------------------------------------------------------------------------------
                        - the input file `pAu_bins.h` is used in almost all of the following

//--------------------------------------------------------------
// On WSU `/tier2/home/groups/rhi/vverkest/pAu_GraceAnalysis`:
//--------------------------------------------------------------
`./make`                            -> `bin/main`
                                    <- `src/{main.cc, events.h, sparse_eta.cc, ...}`

`./sbatch sub-slurm/sparse_eta{,_HT}.sh`  <- `bin/main`
                                    <- `./in-lists/` // point to data at `/tier2/home/groups/rhi/STAR/Data/pAu2015_Grace/`
                                    <- `./in-data/{bad_tower.list, good_run.list}`
                                    -> `out-data/sparse_eta/out-files/sparse_eta_{0, 2, .., 39}.root`

`./grace_loops.py --hadd sparse_eta{,_HT}` <- `out-data/sparse_eta/out-files/sparse_eta_{0, 2, .., 39}.root`
                                     -> `out-data/sparse_eta/hadd.root`

                        
//------------
// run locally
//------------

`plot_zdcx_data_emb.cc`  <- `../loc_lib/loc_libs.h`
                         -> `pdf/$0.pdf`
                         // result : our jet embedded is limiting the analysis to zdcx<=20kHz
                   
`run_make_sparse.cc`     // generate the `make_sparse.root` file
                         <- `../loc_lib/loc_libs.h`
                         <- `hadd_piKp_track_emb.root` // From RCAS
                         <- `recon_pt_bound/{bin_fits.root,  myDiscreteLinFnc.cc, myDiscreteLinFnc.h, sp_trackmatch.cc}` // matches outliers
                         <- `{Reader.h, make_sparse.cc}`  // code to read `hadd_piKp_track_emb.root`
                         <- `make_sparse.cc` // main function to generate `make_sparse.root`
                         -> `make_sparse.root`
                         // skipped: 26525 vs 1998900  rat: 0.0132698 --> dropped 1.3% of events as outliers

`make_weights.cc`        // generate the proper weights for dAu and pp for the priors, using the `make_sparse.root`
                         <- `../loc_lib/loc_libs.h`
                         <- `make_sparse.root` (generated above)
                         <- `TF1_Tsallis.C`
                         <- `Tsallis_params.txt`
                         -> `$0.root`

`run_make_sparse_dAu.cc` // generate the `make_sparse_dAu.root` file
                         <- (all the same inputs as `run_make_sparse.cc` above)
                         <- `make_weights.root`
                         -> `make_sparse_dAu.root`

`run_make_sparse_pp.cc` // generate the `make_sparse_pp.root` file
                         <- (all the same inputs as `run_make_sparse.cc` above)
                         <- `make_weights.root`
                         -> `make_sparse_pp.root`

// using embedding prior as efficiency spectra, calculate total tracking efficiency, with error, w.r.t. ZDCx bins

// using the data pT spectra, unfold for tracking efficiency:
//  (1) w.r.t. ZDCx in bins of EA
//      -> result shows approximate linear dependence, and parameters `m` and `b` are not well determined w.r.t. EA

`./plot_eff10plusall_dAu.cc`         <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root}`
                                     -> `pdf/plot_eff10plusall_dAu.pdf`
                                     // shows that efficiency is not significantly dependent on EA, so use one common efficiency for all EA bins

`./plot_eff10plusall_dAu_eta.cc`     <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root}`
                                     -> `pdf/plot_eff10plusall_dAu_eta.pdf` 
                                     // shows that efficiency is eta dependent and should be used as such

`./plot_eff_unf_1D_vs_2D.cc`         <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_eff_unf_1D_vs_2D.cc` 
                                     // 2D unfolding of data results are within about 0.5% of the 1D unfolding
                                     // therefore wish to use 1D unfolding

`./plot_eff_emb_data.cc`             <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `pdf/$0.pdf`
                                     // shows that efficiency is only embedded on the embedding, not on the pT distributino of the data itself
                                     // (for data at DCA<1)

// show equivalency of PU at 1cm and 3cm in unfolding
`./plot_eff_emb_vs_data_1and3cm.cc`  <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `pdf/plot_eff_emb_vs_data_1and3cm.pdf`
                                     // compares the efficiency correction applied to tracks DCA 0-1 cm vs DCA 1-3 cm
                                     // the tracks at 1-3 cm are much more contaminated by PU, but have within 3% the same
                                     // efficiency. 
                                     // Therefore, we are justified in saying that ``PU is essentially extra MinBias events with
                                     // about the same pT distribution'', and therefore we are ok to ``unfold'' it with the regular
                                     // tracks.

// compare pp vs dAu values -> unfold priors and find efficiency at ZDCx bins with both pp and dAu
`./plot_eff_unf1D_pp_vs_dAu.cc`      <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `./plot_pp_vs_dAu_ZDCx_Eff.cc` : The result is that the pp weightings are about 0.75% when used in the weighting.

// Gives PU functions.
`./plot_PU_dAu.cc`                   <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `pdf/plot_PU_dAu.pdf`
                                     -> `plot_PU_dAu.root` (data points for PU pileup fit)

`./plot_PU_dAu_tripleplot.cc`        <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `pdf/plot_PU_dAu_tripleplot.pdf` (and hand zoomed)
                                     // Same information as `plot_PU_dAu.pdf` but alternate presentation


`calc_NchMB_perEA_{pp,dAu}.cc`       <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root}` :++:
                                     -> `calc_NchMB_perEA_{pp,dAu}.root`

`plot_NchMB_fin.cc`                  <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root, calc_NchMB_perEA_{pp,dAu}.root}`
                                     -> `pdf/$0.root`  // figure of the final values
                                     -> `$0.root` the final values and systemmatic errors

`calc_NchHT_perEA_{pp,dAu}.cc` {2,3,4} <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root}`
                                     -> `$0_{4,8,12}.root`
                                     // the command line input parameter is the trigEt bin, with bin 2\in[4,8]GeV, 3\in[8,12]GeV, 4\in[12,30]GeV

`plot_NchHT_fin.cc` {2,3,4}          <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu_HT.root, calc_NchMB_perEA_{pp,dAu}_{2,3,4}.root}`
                                     -> `pdf/$0.root`  // figure of the final values
                                     -> `$0.root` the final values and systemmatic errors

`plot_Nch_figure.cc`                 <- `{../loc_lib/loc_libs.h, plot_NchMB_fin.root, plot_NchHT_fin_{4,8,12}.root}
                                     -> `pdf/$0.pdf` // this is the figure that goes in the paper

------------------------------------------------------------------------------------------------------------------
// end figure Nch for the paper
------------------------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------------------------
### Figure Nch for the paper
------------------------------------------------------------------------------------------------------------------

`
`plot_hg`
user runs:               comments(//), inputs(<-), output (->)
---------------------    -----------------------------------------------------------------------------------------

                        - `pAu_bins.h` : some comming binning for histograms and THnSparse axes
                        - `../loc_lib/loc_libs.h` : set of local convenience functions
                        The above three files are used almost universally, and lists as `../loc_lib/loc_libs.h` below

//------------
// run locally
//------------
`plot_zdcx_data_emb.cc`  <- `hadd_sparse_wsu.root`
                         -> `pdf/$0.pdf`
                   
`run_make_sparse.cc`     // generate the `make_sparse.root` file
                         <- `../loc_lib/loc_libs.h`
                         <- `hadd_piKp_track_emb.root` // From RCAS (I think)
                         <- `recon_pt_bound/{bin_fits.root,  myDiscreteLinFnc.cc, myDiscreteLinFnc.h, sp_trackmatch.cc}` // matches outliers
                         <- `{Reader.h, make_sparse.cc}`  // code to read `hadd_piKp_track_emb.root`
                         <- `make_sparse.cc` // main function to generate `make_sparse.root`
                         -> `make_sparse.root`
                         // skipped: 26525 vs 1998900  rat: 0.0132698 --> dropped 1.3% of events as outliers

`make_weights.cc`        // generate the proper weights for dAu and pp for the priors, using the `make_sparse.root`
                         <- `../loc_lib/loc_libs.h`
                         <- `make_sparse.root` (generated above)
                         <- `TF1_Tsallis.C`
                         <- `Tsallis_params.txt`
                         -> `$0.root`

`run_make_sparse_dAu.cc` // generate the `make_sparse_dAu.root` file
                         <- (all the same inputs as `run_make_sparse.cc` above)
                         <- `make_weights.root`
                         -> `make_sparse_dAu.root`

`run_make_sparse_pp.cc` // generate the `make_sparse_pp.root` file
                         <- (all the same inputs as `run_make_sparse.cc` above)
                         <- `make_weights.root`
                         -> `make_sparse_pp.root`


//--------------------------------------------------------------
// On WSU `/tier2/home/groups/rhi/vverkest/pAu_GraceAnalysis`:
//--------------------------------------------------------------
`./make`                            -> `bin/main`
                                    <- `src/{main.cc, events.h, sparse_eta.cc, ...}`

`./sbatch sub-slurm/sparse_eta.sh`  <- `bin/main`
                                    <- `./in-lists/` // point to data at `/tier2/home/groups/rhi/STAR/Data/pAu2015_Grace/`
                                    <- `./in-data/{bad_tower.list, good_run.list}`
                                    -> `out-data/sparse_eta/out-files/sparse_eta_{0, 2, .., 39}.root`

`./grace_loops.py --hadd sparse_eta` <- `out-data/sparse_eta/out-files/sparse_eta_{0, 2, .., 39}.root`
                                     -> `out-data/sparse_eta/hadd.root`


//------------
// run locally
//------------
`scp {userid}@grid.wayne.edu:/tier2/home/groups/rhi/vverkest/pAu_GraceAnalysis/out-data/sparse_eta/hadd.root hadd_sparse_wsu.root` // moves the data THnSparse locally

// using embedding prior as efficiency spectra, calculate total tracking efficiency, with error, w.r.t. ZDCx bins

// using the data pT spectra, unfold for tracking efficiency:
//  (1) w.r.t. ZDCx in bins of EA
//      -> result shows approximate linear dependence, and parameters `m` and `b` are not well determined w.r.t. EA

`./plot_eff10plusall_dAu.cc`         <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root}`
                                     -> `pdf/plot_eff10plusall_dAu.pdf`
                                     // shows that efficiency is not significantly dependent on EA, so use one common efficiency for all EA bins

`./plot_eff10plusall_dAu_eta.cc`     <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root}`
                                     -> `pdf/plot_eff10plusall_dAu_eta.pdf` 
                                     // shows that efficiency is eta dependent and should be used as such

`./plot_eff_unf_1D_vs_2D.cc`         <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_eff_unf_1D_vs_2D.cc` 
                                     // 2D unfolding of data results are within about 0.5% of the 1D unfolding
                                     // therefore wish to use 1D unfolding

`./plot_eff_emb_data.cc`             <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_Eff_Emb_data.cc`
                                     // shows that efficiency is only embedded on the embedding, not on the pT distributino of the data itself
                                     // (for data at DCA<1)

// show equivalency of PU at 1cm and 3cm in unfolding
`./plot_eff_emb_vs_data_1and3cm.cc`  <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_eff_emb_vs_data_1and3cm.pdf`
                                     // compares the efficiency correction applied to tracks DCA 0-1 cm vs DCA 1-3 cm
                                     // the tracks at 1-3 cm are much more contaminated by PU, but have within 3% the same
                                     // efficiency. 
                                     // Therefore, we are justified in saying that ``PU is essentially extra MinBias events with
                                     // about the same pT distribution'', and therefore we are ok to ``unfold'' it with the regular
                                     // tracks.

// compare pp vs dAu values -> unfold priors and find efficiency at ZDCx bins with both pp and dAu
`./plot_eff_unf1D_pp_vs_dAu.cc`      <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}` +::+
                                     -> `./plot_pp_vs_dAu_ZDCx_Eff.cc` : The result is that the pp weightings are about 0.75% when used in the weighting.

// Gives PU functions.
`./plot_PU_dAu.cc`                   <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_PU_dAu.pdf`
                                     -> `plot_PU_dAu.root` (data points for PU pileup fit)

`./plot_PU_dAu_tripleplot.cc`        <- `{../loc_lib/loc_libs.h, make_sparse_dAu.root, hadd_sparse_wsu.root}`
                                     -> `pdf/plot_PU_dAu_tripleplot.pdf` (and hand zoomed)
                                     // Same information as `plot_PU_dAu.pdf` but alternate presentation


`calc_NchMB_perEA_{pp,dAu}.cc`       <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root}`
                                     -> `calc_NchMB_perEA_{pp,dAu}.root`

`plot_NchMB_fin.cc`                  <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root, calc_NchMB_perEA_{pp,dAu}.root}`
                                     -> `pdf/$0.root`  // figure of the final values
                                     -> `$0.root` the final values and systemmatic errors

`calc_NchHT_perEA_{pp,dAu}.cc` {2,3,4} <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu.root}`
                                     -> `$0_{4,8,12}.root`
                                     // the command line input parameter is the trigEt bin, with bin 2\in[4,8]GeV, 3\in[8,12]GeV, 4\in[12,30]GeV

`plot_NchHT_fin.cc` {2,3,4}          <- `{../loc_lib/loc_libs.h, make_sparse_{pp,dAu}.root, hadd_sparse_wsu_HT.root, calc_NchMB_perEA_{pp,dAu}_{2,3,4}.root}`
                                     -> `pdf/$0.root`  // figure of the final values
                                     -> `$0.root` the final values and systemmatic errors

`plot_Nch_figure.cc`                 <- `{../loc_lib/loc_libs.h, plot_NchMB_fin.root, plot_NchHT_fin_{4,8,12}.root}`
                                     -> `pdf/$0.pdf` // this is the figure that goes in the paper

------------------------------------------------------------------------------------------------------------------
// end figure Nch for the paper
------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------
`plot_hgEA{,_HT}.cc`                <- `{../loc_lib/loc_libs.h, loc_libs.h}`
                                    -> `pdf/$0.pdf`
------------------------------------------------------------------------------------------------------------------



## [ ] uncertainty correlation: for `pdf/plot_Nch_figure.pdf`, need a statement for the correlation of the uncertainties 
between high and low EA. (They are not totally correlated, need a general value). 

##  [ ] Fix `_nhitfix` from 16 to 15 on the minimum

# Generate THnSparse* objects for embedded tracks

