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
//#include "THistPainter.h."

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TString name;
  
  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10<p_{T}^{lead}<15", "15<p_{T}^{lead}<20",  "20<p_{T}^{lead}<30" };
  const double AREA = 4*(1.14159265);
  const int marker[55] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };

  const int nEtaBins = 3;
  // const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  // const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };

  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };

  // OPEN RESPONSE FILE AND COLLECT HISTOS
  TFile *inFile = new TFile("out/sim/pAu2015embedding.root","READ");

  TH2D *hResponse[nEtaBins];
  TH1D *hFakes[nEtaBins];
  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  int nAccepted = 0; int nFakes = 0;
  
  for (int e=0; e<nEtaBins; ++e) {    
    name = "hPtResponse" + etaBinName[e] + "Jet";
    hResponse[e] = (TH2D*)inFile->Get(name);

    name = "hFakes" + etaBinName[e] + "Jet";
    hFakes[e] = (TH1D*)inFile->Get(name);

    nAccepted += hResponse[e]->GetEntries();
    nFakes += hFakes[e]->GetEntries();

  }

  TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.5,59.5, 55,4.5,59.5);
  TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.5,59.5);
  
  for (int e=0; e<nEtaBins; ++e) {    
    hPtResponse->Add(hResponse[e]);
    hFakeJets->Add(hFakes[e]);
  }
  
  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );              // CANVAS 0
  can->SetLogz();
  hPtResponse->Draw("COLZ");
  can->SaveAs("plots/pTresponse.pdf","PDF");  // SAVE 2D PT RESPONSE
  
  TH1D *hDet[55]; TH1D *hDetWt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    TString name = "hPtResponse"; name += ptBinName[p];
    hDetWt[p] = new TH1D( name, ";det-level leading jet p_{T} (GeV)",55,4.5,59.5);
  }

  
  for (int i=5; i<26; ++i) {  // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
    int ptVal = i + 5;
    int binno = i + 1;
    
    TString name = "hPtResponse_"; name += ptVal; name += "GeV";
    hDet[i] = (TH1D*) hPtResponse->ProjectionY(name,binno,binno);
    hDet[i]->Scale(1./hDet[i]->Integral());
    hDet[i]->Draw();
    name = "plots/PtResponse_"; name += ptVal; name += "GeV"; name += ".pdf";
    can->SaveAs( name, "PDF" );
    hDet[i]->SetMarkerStyle(marker[i]);
    hDet[i]->SetMarkerSize(1);
    hDet[i]->SetStats(0);
    name = ""; name += ptVal; name += " GeV part. jet";
    hDet[i]->SetNameTitle(name,";det-level leading jet p_{T} (GeV)");
  }
  
  hDet[5]->Draw("PLC PMC");
  for (int i=6; i<26; ++i) { hDet[i]->Draw("SAME PLC PMC"); }
  can->BuildLegend(0.68,0.1,0.9,0.9);
  can->SaveAs("plots/detPtResponses.pdf","PDF");

  for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
    for (int i=5; i<26; ++i) {
      int ptVal = i + 5;
      int binno = i + 1;
      if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	double wt = 1./hMissedJets->GetBinContent(binno);  // weight according to misses at part-level
	hDetWt[p]->Add(hDet[i],wt);
      }
    }    
    TString name = "plots/WeightedPtResponse"; name += ptBinName[p]; name += ".pdf";
    hDetWt[p]->Scale(1./hDetWt[p]->Integral());
    hDetWt[p]->Draw();
    can->SaveAs(name,"PDF");
  }


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

  // TH3D *hUE3D = (TH3D*)UEfile->Get("hChgUE");  // leadPt:UEpt:UEeta
  // TH1D *hleadPt = (TH1D*)UEfile->Get("hLeadPt");  // leadPt

  
  // SAVE 2D HISTO OF UE PT vs. LEADING JET PT
  TH2D *hChgUE2D = (TH2D*)hChgUE3D->Project3D("YX");
  hChgUE2D->GetXaxis()->SetRangeUser(0.0,60.0);  hChgUE2D->GetYaxis()->SetRangeUser(0.0,15.0);
  hChgUE2D->GetYaxis()->SetTitleOffset(1.25);
  hChgUE2D->Draw("COLZ");
  can->SaveAs("plots/hChgUE2D.pdf","PDF");
  hChgUE2D->GetXaxis()->SetRangeUser(4.5,59.5);  hChgUE2D->GetYaxis()->SetRangeUser(0.0,30.0);
 

  TH1D *hUEpt[55];  TH1D *hWtUEpt[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    TString name = "hWtUEpt"; name += ptBinName[p];
    TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
    hWtUEpt[p] = new TH1D(name,title,30,0.0,30.0);
  }

  
  can->SetLogy();
  for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
    int ptVal = i + 5;
    int binno = i + 1;
    TString name = "hUEpt_"; name += ptVal; name += "GeV";
    hUEpt[i] = (TH1D*)hChgUE2D->ProjectionY(name,binno,binno);
    if ( (int)hleadPt->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }
    int nJets = hleadPt->GetBinContent( binno );
    if (nJets>0){
      hUEpt[i]->Scale(1./nJets);
      hUEpt[i]->Draw();
      name = "plots/UEpt_"; name += ptVal; name += "GeV"; name += ".pdf";
      can->SaveAs(name,"PDF");
      
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

  hUEpt[0]->GetYaxis()->SetRangeUser(0.000005,10);
  hUEpt[0]->Draw("PLC PMC");

  for (int i=1; i<55; ++i) { if (hUEpt[i]->GetEntries()!=0) { hUEpt[i]->Draw("SAME PLC PMC"); } }
  can->BuildLegend(0.75,0.1,0.9,0.9);
  can->SaveAs("plots/detUEpt.pdf","PDF");


  TFile *ef = new TFile( "src/trackeffic_allEta.root", "READ" );

  TH1D *hEffic = (TH1D*)ef->Get( "eff_s_bin_1_10_bbc__1_10_eta" );

  double pt, effic, corrPt, corrErr; // eta, phi, e, px, py, pz, 
  int ptBin;//, etaBin;

  
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
    
    hWtUEpt[p]->Draw();

    TString text = "#LT p_{T}^{ch}#GT = "; text += hWtUEpt[p]->GetMean(1);
    text = text(0,26);

    TLatex *tex = new TLatex( 0.6, 0.55,text);
    tex->SetTextFont(63);
    tex->SetTextSize(20);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();

    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += hWtUEpt[p]->Integral()/AREA;
    text = text(0,42);
    tex->DrawLatex( 0.6, 0.4, text );
      
    TString name = "plots/wtUEpt"; name += ptBinName[p]; name += ".pdf";
    can->SaveAs(name,"PDF");
  }



  hMissedJets->Draw();
  can->SaveAs("plots/missedJets.pdf","PDF");
  hFakeJets->Draw();
  can->SaveAs("plots/fakeJets.pdf","PDF");

  
  return 0;
}
