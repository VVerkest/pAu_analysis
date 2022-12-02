#define Reader_cxx
#include "Reader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Reader::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Reader.C
//      root> Reader t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

// --Get out of bounds for reconstructed tracks--------------------------
    TFile *s_current = gDirectory->GetFile();
    TFile *fin = new TFile("recon_pt_bound/bin_fits.root","read");
    vector<double> *_bins, *m_lo, *b_lo, *m_hi, *b_hi;
    fin->GetObject("bins",_bins);
    array<myDiscreteLinFnc*,6> fn_lo, fn_hi;
    for (int i=0; i<6; ++i) {
        fin->GetObject(Form("m_lo_%s",noi_geant05_ascii(i)), m_lo);
        fin->GetObject(Form("m_hi_%s",noi_geant05_ascii(i)), m_hi);
        fin->GetObject(Form("b_lo_%s",noi_geant05_ascii(i)), b_lo);
        fin->GetObject(Form("b_hi_%s",noi_geant05_ascii(i)), b_hi);
        fn_lo[i] = new myDiscreteLinFnc(*_bins, *b_lo, *m_lo);
        fn_hi[i] = new myDiscreteLinFnc(*_bins, *b_hi, *m_hi);
    };
    s_current->cd();    
// ----------------------------------------------------------------------

// --Get weighting for the particle types -------------------------------
array<TH1D*, 6> weight;
noiGetter got;
for (int i=0;i<6;++i) {
    weight[i] = (TH1D*) got("make_weights.root", Form("w_%i_pp",i));
}
// ----------------------------------------------------------------------

    TFile fout { "make_sparse_pp.root", "recreate" };
    fout.cd();
    noiMsgTree msg_tree{false};
    msg_tree.slurp_file("make_sparse_pp.cc");

    int nbins[9];
    nbins[_geantid] = bin_id;      // reco + emb.
    nbins[_bbc]     = bin_bbc;     // reco + emb.
    nbins[_zdcx]    = bin_zdcx;    // reco + emb.
    nbins[_pt]      = bin_pt;      // reco + emb.
    nbins[_eta]     = bin_eta;     // reco + emb.
    nbins[_dca]     = bin_dca;     // reco
    nbins[_nhitfit] = bin_hitsfit; // reco
    nbins[_nhitrat] = bin_fitrat;  // reco
    nbins[_ptreco]  = bin_pt;      // reco

    auto sp_truth = new THnSparseD ("sp_truth", "emb_tracks;particle-id;bbcES;ZDCx;pT_{truth};eta_{truth};",
         5, nbins, NULL, NULL);
    sp_truth->SetBinEdges(_geantid, bin_id);
    sp_truth->SetBinEdges(_bbc,     bin_bbc);
    sp_truth->SetBinEdges(_zdcx,    bin_zdcx);
    sp_truth->SetBinEdges(_pt,      bin_pt);
    sp_truth->SetBinEdges(_eta,     bin_eta);
    sp_truth->Sumw2();

    auto sp_reco = new THnSparseD ("sp_reco", 
            "emb_reco;particle-id;bbcES;ZDCx;pT_{truth};eta;dca;NHitsFit;NHitsFit/NHitsPoss;pT_{reco};",
         9, nbins, NULL, NULL);
    sp_reco->SetBinEdges(_geantid, bin_id);
    sp_reco->SetBinEdges(_bbc,     bin_bbc);
    sp_reco->SetBinEdges(_zdcx,    bin_zdcx);
    sp_reco->SetBinEdges(_pt,      bin_pt);
    sp_reco->SetBinEdges(_eta,     bin_eta);
    sp_reco->SetBinEdges(_dca,     bin_dca);
    sp_reco->SetBinEdges(_nhitfit, bin_hitsfit);
    sp_reco->SetBinEdges(_nhitrat, bin_fitrat);
    sp_reco->SetBinEdges(_ptreco,  bin_pt);
    sp_reco->Sumw2();

    double hopper [9];

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
  long long int skip = 0;
  long long int keep = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      hopper[_bbc] = BBC_Ein;
      hopper[_zdcx] = ZDCx;

      for (int i=0; i<mc_track_; ++i) {
          int k = noi_geant05(mc_track_geantId[i]);
          if (k==-1) continue;
          hopper[_geantid] = k;
          float T = mc_track_pt[i];
          double W = weight[k]->GetBinContent( bin_pt.bin_from_1(T) );
          hopper[_pt] = T;
          hopper[_eta] = mc_track_eta[i];

          int id = mc_track_id[i];
          float M = (id >= 0) ? track_pt[id] : 0.; 
          if (id >=0 && (M < (*fn_lo[k])(T) || M > (*fn_hi[k])(T)))  { // cut pT outlier events
              ++skip;
              continue;
          }
          ++keep;
          sp_truth->Fill(hopper,W);
          if (id >= 0 && track_dcaXYZ[id]<3.) {
              hopper[_dca]     = track_dcaXYZ[id];
              hopper[_nhitfit] = track_nHitsFit[id];
              hopper[_nhitrat] = (float) track_nHitsFit[id] / (float) track_nHitsPoss[id];
              hopper[_ptreco]  = M;
              sp_reco->Fill(hopper,W);
          }
       }
   }
   cout << "skipped: " << skip << " vs " << keep << "  rat: " << (double) skip / keep << endl;
   sp_truth->Write();
   sp_reco->Write();
   msg_tree.write();
   fout.Write();
   fout.Save();
}
