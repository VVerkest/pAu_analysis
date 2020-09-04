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

  TH3D *hChgUE3D = new TH3D("hChgUE3D",";chg. UE part. p_{T} (GeV);chg. UE part. #eta", 55,4.5,59.5, 30,0.0,30.0, 40,-1.0,1.0);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.5,59.5);

  for (int e=0; e<nEtaBins; ++e) {
    hChgUE3D->Add(hUE3D[e]);
    hleadPt->Add(hLeadJetPt[e]);
  }
  
  TH2D *hChgUE2D = (TH2D*)ProjectUEHistograms( hChgUE3D, "plots/test/" );  
  TH1D *hUEpt[55];  TH1D *hWtUEpt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    name = "hWtUEpt"; name += ptBinName[p];
    TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
    hWtUEpt[p] = new TH1D(name,title,30,0.0,30.0);
  }

  ProjectAndScaleUEHistogramForAllPt( hChgUE2D, hleadPt, hUEpt, hWtUEpt, "plots/test/" );

  WeightUEPtByLeadPtAndFakes( hWtUEpt, hUEpt, hDetWt, hleadPt, hFakeJets, "plots/test/" );


  // ########################################################################################################
  //                              TRACKING EFFICIENCY FILE FOR UE pT CORRECTION
  // ########################################################################################################
  
  TFile *ef = new TFile( "src/trackeffic_allEta.root", "READ" );

  TH1D *hEffic = (TH1D*)ef->Get( "eff_s_bin_1_10_bbc__1_10_eta" );

  TCanvas * can8 = new TCanvas( "can8" , "" ,700 ,500 );              // CANVAS 8
  can8->SetLogy();
  
  double pt, effic, corrPt, corrErr;
  int ptBin;
  for (int p=0; p<nPtBins; ++p) {

    for ( int i=1; i<=hWtUEpt[p]->GetNbinsX(); ++i ) {

      pt = hWtUEpt[p]->GetBinCenter(i);  // loop over chg UE pT bins
      if ( pt > 3.0 ) { pt = 3.0; }
      ptBin = hEffic->FindBin( pt );    // find histogram bin corresponding to track pt
      effic = hEffic->GetBinContent( ptBin );
      corrPt = hWtUEpt[p]->GetBinContent(i)/effic;                     // calculate corrected bin content and error
      corrErr = hWtUEpt[p]->GetBinError(i)/effic;                      // (divide bin content and error by efficiency)

      hWtUEpt[p]->SetBinContent( i, corrPt );
      hWtUEpt[p]->SetBinError( i, corrErr );
    }

    hWtUEpt[p]->SetStats(1);
    hWtUEpt[p]->SetAxisRange( 0.000001,10,"Y");
    hWtUEpt[p]->SetLineColor(kBlack);
    hWtUEpt[p]->SetMarkerColor(kBlack);
    hWtUEpt[p]->Draw();

    TString text = "#LT p_{T}^{ch}#GT = "; text += hWtUEpt[p]->GetMean(1);
    text = text(0,26);
    DrawText(text, 0.6, 0.7, 20);

    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += hWtUEpt[p]->Integral()/AREA;
    text = text(0,42);
    DrawText( text, 0.6, 0.55, 20 );
      
    name = "plots/test/CorrectedWtUEpt"; name += ptBinName[p]; name += ".pdf";
    can8->SaveAs(name,"PDF");
  }

  
  return 0;

}
