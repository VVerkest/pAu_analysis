// Veronica Verkest
// August 31, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString name, title;
  TString directory = "plots/test/";

  TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);
  TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.0,59.0);
  TH1D *hDet[21];
  TH1D *hDetWt[nPtBins];
  TH1D *hLeadJetPt[nEtaBins];
  TH3D *hUE3D[nEtaBins];
  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.0,59.0);
  TH2D *hChgUE2D[20];
  TH2D *hChgUE2D_corr[20];
  TH2D *hAddedChgUE2D_corr[nEtaBins];
  TH1D *hUEpt[55];
  TH1D *hWtUEpt[nEtaBins][nPtBins];
  // TH1D *hWtUEpt[nEAbins][nEtaBins][nPtBins];

  const int n_bins = 3;
  double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  const int n_ybins = 3;
  double y_bin_edge[n_ybins+1] = { 0.55,0.65,0.75,0.85 };
  TH2D *meanPt_hscale = new TH2D("meanPt_hscale",";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge);
  double y_bin_edge2[n_ybins+1] = { 0.5,1.2,1.5,1.8 };
  TH2D *nCh_hscale = new TH2D("nCh_hscale",";(GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge2);

  for (int p=0; p<nPtBins; ++p) {
    for (int e=0; e<nEtaBins; ++e) {
      name = "hWtUEpt"; name += ptBinName[p] + etaBinName[e];
      TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
      hWtUEpt[e][p] = new TH1D(name,title,ybins,ybinEdge);
    }
  }

  for (int e=0; e<nEtaBins; ++e) {
    name = "hChgUE2D_corr_" + emw[e] + "UE";
    hAddedChgUE2D_corr[e] = new TH2D( name,";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)", xbins,xbinEdge,ybins,ybinEdge);
  }

  

  // for (int a=1; a>=0; --a) {
      
    ////////////////////////////////// EMBEDDING FILE FOR LEADING JET PT CORRECTION //////////////////////////////////
    TFile *inFile = new TFile("out/sim/pAu2015embedding_hiEA.root","READ");  // OPEN RESPONSE FILE AND COLLECT HISTOS

    TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

    GetEmbeddingHistograms( inFile, hPtResponse, hFakeJets, hMissedJets, directory );

    ProjectPartLevelJetPt( hPtResponse, hDet, directory );

    GenerateWeightedPtResponse( hDetWt, hDet, hMissedJets, directory );

    ////////////////////////////////////////// pAu DATA FILE FOR UE HISTOGRAMS //////////////////////////////////////////
    TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_hiEA_uncorrected.root","READ");  // OPEN pAu FILE AND GATHER HISTOGRAMS OF UE AND LEAD PT
  
    for (int e=0; e<nEtaBins; ++e) {
      name = "hLeadPt" + etaBinName[e] + "Jet";
      hLeadJetPt[e] = (TH1D*)UEfile->Get(name);
      name = "hChgUE" + etaBinName[e] + "Jet";
      hUE3D[e] = (TH3D*)UEfile->Get(name);
      hUE3D[e]->Scale(hUE3D[e]->GetEntries()/hUE3D[e]->Integral());
    }

    for (int e=0; e<nEtaBins; ++e) {
      hChgUE3D->Add(hUE3D[e]);
      hleadPt->Add(hLeadJetPt[e]);
      hUE3D[e]->Delete();
      hLeadJetPt[e]->Delete();
    }
  
    for (int z=0; z<zbins; ++z) {
      hChgUE2D[z] = (TH2D*)ProjectUEHistograms( hChgUE3D, z, "plots/test/" );
      // cout<<hChgUE2D[z]->GetName()<<"  "<<hChgUE2D[z]->GetMean(1)<<endl;
      // cout<<hChgUE2D[z]<<endl;
      name = "hChgUE2D_corr_etaBin"; name += z;
      hChgUE2D_corr[z] = new TH2D( name,";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)", xbins,xbinEdge,ybins,ybinEdge);    
    }

    hChgUE3D->Delete();
    ////////////////////////////////// TRACKING EFFICIENCY FILE FOR UE pT CORRECTION //////////////////////////////////  
    TFile *ef = new TFile( "src/trackeffic.root", "READ" );

    for (int e=0; e<10; ++e) {
      name = "eff_s_bin_8_10_bbc__"; name += e+1; name += "_"; name += e+1; name += "_eta";
      TH1D *hEffic = (TH1D*)ef->Get( name );//_3_3, _8_10
      TrackingEfficiency2DCorrection( hChgUE2D_corr[2*e], hChgUE2D[e], hEffic, directory );
      TrackingEfficiency2DCorrection( hChgUE2D_corr[2*e+1], hChgUE2D[2*e+1], hEffic, directory );
    }

    for (int i=0; i<=7; ++i) { hAddedChgUE2D_corr[0]->Add(hChgUE2D_corr[i]); }
    for (int i=7; i<=13; ++i) { hAddedChgUE2D_corr[1]->Add(hChgUE2D_corr[i]); }
    for (int i=13; i<=20; ++i) { hAddedChgUE2D_corr[2]->Add(hChgUE2D_corr[i]); }
  
    for (int e=0; e<nEtaBins; ++e) {
      ProjectAndScaleUEHistogramForAllPt( hAddedChgUE2D_corr[e], hleadPt, hUEpt, hWtUEpt[e], "plots/test/", emw[e] );
      WeightUEPtByLeadPtAndFakes( hWtUEpt[e], hUEpt, hDetWt, hleadPt, hFakeJets, "plots/test/", emw[e] );
    }
  
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	TString suf = ptBinName[p]+etaBinName[e];
	ProjectAndSaveFinalUEPlots( hWtUEpt[e][p], suf, area[e], "plots/test/" );
      }
    }

    StackAndSavePtPlots( hWtUEpt, meanPt_hscale, "meanPt.pdf", directory );
    StackAndSaveNchPlots( hWtUEpt, nCh_hscale, "nCh.pdf", directory );

    TFile *outFile = new TFile("out/ptResponse.root","RECREATE");
  
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	hWtUEpt[e][p]->Write();
      }
    }
  // }

  return 0;

}
