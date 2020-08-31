// Veronica Verkest
// August 31, 2020

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TAttMarker.h"

#include <string>
#include <iostream>
#include <sstream>

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

  TH2D *hResponse[nEtaBins];
  TH1D *hFakes[nEtaBins];
  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  int nAccepted = 0; int nFakes = 0;
  
  for (int e=0; e<nEtaBins; ++e) {    
    name = "hPtResponse" + etaBinName[e] + "Jet";     hResponse[e] = (TH2D*)inFile->Get(name);

    name = "hFakes" + etaBinName[e] + "Jet";          hFakes[e] = (TH1D*)inFile->Get(name);

    nAccepted += hResponse[e]->GetEntries();
    nFakes += hFakes[e]->GetEntries();
  }

  TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.5,59.5, 55,4.5,59.5);
  TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.5,59.5);
  
  for (int e=0; e<nEtaBins; ++e) {     hPtResponse->Add(hResponse[e]);     hFakeJets->Add(hFakes[e]);   }

  TCanvas * can1 = new TCanvas( "can1" , "" ,700 ,500 );              // CANVAS 1
  can1->SetLogz();
  TCanvas * can2 = new TCanvas( "can2" , "" ,700 ,500 );              // CANVAS 2
  can2->SetLogy();

  int nMissed = hMissedJets->GetEntries();
  int nEvents = nAccepted + nMissed + nFakes;
  
  double scale = (double)nMissed/nEvents;
  hMissedJets->Scale(scale/hMissedJets->Integral());
  hMissedJets->Draw();
  can2->SaveAs("plots/test/MissedJets.pdf","PDF");

  TCanvas * can3 = new TCanvas( "can3" , "" ,700 ,500 );              // CANVAS 3
  can3->SetLogz();
  hPtResponse->Draw("COLZ");
  can3->SaveAs("plots/test/pTresponse.pdf","PDF");  // SAVE 2D PT RESPONSE
  
  TH1D *hDet[55]; TH1D *hDetWt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    name = "hPtResponse"; name += ptBinName[p];
    hDetWt[p] = new TH1D( name, ";det-level leading jet p_{T} (GeV)",55,4.5,59.5);
  }

  TCanvas * can4 = new TCanvas( "can4" , "" ,700 ,500 );              // CANVAS 4
  
  for (int i=5; i<26; ++i) {  // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
    int ptVal = i + 5;
    int binno = i + 1;
    
    name = "hPtResponse_"; name += ptVal; name += "GeV";
    hDet[i] = (TH1D*) hPtResponse->ProjectionY(name,binno,binno);
    hDet[i]->Scale(1./hDet[i]->Integral());
    hDet[i]->SetMarkerStyle(marker[i]);
    hDet[i]->Draw("SAME PLC PMC");
    name = "plots/test/PtResponse_"; name += ptVal; name += "GeV"; name += ".pdf";
    hDet[i]->SetMarkerStyle(marker[i]);
    hDet[i]->SetMarkerSize(1);
    hDet[i]->SetStats(0);
    name = ""; name += ptVal; name += " GeV part. jet";
    hDet[i]->SetNameTitle(name,";det-level leading jet p_{T} (GeV)");
  }
  
  can4->BuildLegend(0.68,0.1,0.9,0.9);
  can4->SaveAs("plots/test/detPtResponses.pdf","PDF");

  TCanvas * can5 = new TCanvas( "can5" , "" ,700 ,500 );              // CANVAS 5

  for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
    for (int i=5; i<26; ++i) {
      int ptVal = i + 5;
      int binno = i + 1;
      if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	double wt = 1./hMissedJets->GetBinContent(binno);  // weight according to misses at part-level
	hDetWt[p]->Add(hDet[i],wt);
      }
    }    
    hDetWt[p]->Scale(1./hDetWt[p]->Integral());
    hDetWt[p]->GetYaxis()->SetRangeUser(0.0,0.2);
    hDetWt[p]->SetMarkerStyle(ptMarker[p]);
    hDetWt[p]->Draw("SAME PLC PMC");
  }
  can5->BuildLegend(0.4,0.68,0.9,0.9);
  can5->SaveAs("plots/test/WeightedPtResponse.pdf","PDF");



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

  for (int e=0; e<nEtaBins; ++e) {     hChgUE3D->Add(hUE3D[e]);     hleadPt->Add(hLeadJetPt[e]);   }

  can5->SetLogy(0);
  can5->SetLogz();
  
  // SAVE 2D HISTO OF UE PT vs. LEADING JET PT
  TH2D *hChgUE2D = (TH2D*)hChgUE3D->Project3D("YX");
  hChgUE2D->GetXaxis()->SetRangeUser(0.0,60.0);  hChgUE2D->GetYaxis()->SetRangeUser(0.0,15.0);
  hChgUE2D->GetYaxis()->SetTitleOffset(1.25);
  hChgUE2D->Draw("COLZ");
  can5->SaveAs("plots/test/hChgUE2D.pdf","PDF");
  hChgUE2D->GetXaxis()->SetRangeUser(4.5,59.5);  hChgUE2D->GetYaxis()->SetRangeUser(0.0,30.0);
 

  TH1D *hUEpt[55];  TH1D *hWtUEpt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    name = "hWtUEpt"; name += ptBinName[p];
    TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
    hWtUEpt[p] = new TH1D(name,title,30,0.0,30.0);
  }

  TCanvas * can6 = new TCanvas( "can6" , "" ,700 ,500 );              // CANVAS 6

  can6->SetLogy();
  for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
    int ptVal = i + 5;
    int binno = i + 1;
    name = "hUEpt_"; name += ptVal; name += "GeV";
    hUEpt[i] = (TH1D*)hChgUE2D->ProjectionY(name,binno,binno);
    if ( (int)hleadPt->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }
    int nJets = hleadPt->GetBinContent( binno );
    if (nJets>0){
      hUEpt[i]->Scale(1./nJets);
      hUEpt[i]->SetAxisRange( 0.00001,10, "Y");
      hUEpt[i]->Draw("SAME PLC PMC");
      
      for (int p=0; p<nPtBins; ++p) {  // ADD UE DISTRIBUTIONS WEIGHTED BY DET-LEVEL FRACTIONAL CONTRIBUTION
	double wt = hDetWt[p]->GetBinContent(binno)*( 1.0 - hFakeJets->GetBinContent(binno) );
	hWtUEpt[p]->Add( hUEpt[i], wt );
      }
      hUEpt[i]->SetMarkerStyle(marker[i]);
      hUEpt[i]->SetMarkerSize(1);
      hUEpt[i]->SetStats(0);
      name = ""; name += ptVal; name += " GeV det. jet";
      hUEpt[i]->SetNameTitle(name,";chg. UE particle p_{T} (GeV)");
    }
  }
  can6->BuildLegend(0.68,0.1,0.9,0.9);
  can6->SaveAs("plots/test/UEptByLeadPt.pdf","PDF");

  TCanvas * can7 = new TCanvas( "can7" , "" ,700 ,500 );              // CANVAS 7
  for (int p=0; p<nPtBins; ++p) {
    hWtUEpt[p]->SetMarkerStyle(ptMarker[p]);
    hWtUEpt[p]->Draw("PLC PMC SAME");
    // hWtUEpt[p]->;
  }
  can7->SetLogy();
  can7->BuildLegend(0.4,0.68,0.9,0.9);
  can7->SaveAs("plots/test/weightedUEptByLeadPt.pdf","PDF");


  
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

    hWtUEpt[p]->GetXaxis()->SetRangeUser(0.0,15.0);
    hWtUEpt[p]->SetLineColor(kBlack);
    hWtUEpt[p]->SetMarkerColor(kBlack);
    hWtUEpt[p]->Draw();

    TString text = "#LT p_{T}^{ch}#GT = "; text += hWtUEpt[p]->GetMean(1);
    text = text(0,26);
    drawText(text, 0.6, 0.7, 20);

    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += hWtUEpt[p]->Integral()/AREA;
    text = text(0,42);
    drawText( text, 0.6, 0.55, 20 );
      
    name = "plots/test/CorrectedWtUEpt"; name += ptBinName[p]; name += ".pdf";
    can8->SaveAs(name,"PDF");
  }

  hMissedJets->Draw();
  can8->SaveAs("plots/test/missedJets.pdf","PDF");
  hFakeJets->Draw();
  can8->SaveAs("plots/test/fakeJets.pdf","PDF");

  
  return 0;

}
