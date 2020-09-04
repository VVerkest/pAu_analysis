// Veronica Verkest
// August 31, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString name, title;

  // ########################################################################################################
  //                           EMBEDDING FILE FOR LEADING JET PT CORRECTION
  // ########################################################################################################
  
  // OPEN RESPONSE FILE AND COLLECT HISTOS
  TFile *inFile = new TFile("out/sim/pAu2015embedding.root","READ");

  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.5,59.5, 55,4.5,59.5);
  TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.5,59.5);

  TString directory = "plots/test/";
  TH1D *hDet[21];
  TH1D *hDetWt[nPtBins];

  GetEmbeddingHistograms( inFile, hPtResponse, hFakeJets, hMissedJets, directory );

  ProjectPartLevelJetPt( hPtResponse, hDet, directory );

  GenerateWeightedPtResponse( hDetWt, hDet, hMissedJets, directory );


  // ########################################################################################################
  //                                 pAu DATA FILE FOR UE HISTOGRAMS
  // ########################################################################################################

  // OPEN pAu FILE AND GATHER HISTOGRAMS OF UE AND LEAD PT
  TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_allEA.root","READ");

  TH1D *hLeadJetPt[nEtaBins];
  TH3D *hUE3D[nEtaBins];
  
  for (int e=0; e<nEtaBins; ++e) {
  
    name = "hLeadPt" + etaBinName[e] + "Jet";
    hLeadJetPt[e] = (TH1D*)UEfile->Get(name);

    name = "hChgUE" + etaBinName[e] + "Jet";
    hUE3D[e] = (TH3D*)UEfile->Get(name);    
  }

  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", 55,4.5,59.5, 30,0.0,30.0, 40,-1.0,1.0);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.5,59.5);

  for (int e=0; e<nEtaBins; ++e) {
    hChgUE3D->Add(hUE3D[e]);
    hleadPt->Add(hLeadJetPt[e]);
  }
  
  TH2D *hChgUE2D = (TH2D*)ProjectUEHistograms( hChgUE3D, "plots/test/" );
  TH2D *hChgUE2D_corr = new TH2D("hChgUE2D_corr",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV)", 55,4.5,59.5, 30,0.0,30.0);


  // ########################################################################################################
  //                              TRACKING EFFICIENCY FILE FOR UE pT CORRECTION
  // ########################################################################################################
  
  TFile *ef = new TFile( "src/trackeffic_allEta.root", "READ" );

  TH1D *hEffic = (TH1D*)ef->Get( "eff_s_bin_1_10_bbc__1_10_eta" );

  TrackingEfficiency2DCorrection( hChgUE2D_corr, hChgUE2D, hEffic );
  
  TH1D *hUEpt[55];  TH1D *hWtUEpt[nPtBins];
  
  for (int p=0; p<nPtBins; ++p) {
    name = "hWtUEpt"; name += ptBinName[p];
    TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
    hWtUEpt[p] = new TH1D(name,title,30,0.0,30.0);
  }

  ProjectAndScaleUEHistogramForAllPt( hChgUE2D_corr, hleadPt, hUEpt, hWtUEpt, "plots/test/" );

  WeightUEPtByLeadPtAndFakes( hWtUEpt, hUEpt, hDetWt, hleadPt, hFakeJets, "plots/test/" );


  for (int p=0; p<nPtBins; ++p) { ProjectAndSaveFinalUEPlots( hWtUEpt[p], ptBinName[p], "plots/test/" ); }

  return 0;

}
