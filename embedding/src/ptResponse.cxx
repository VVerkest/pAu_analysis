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
  
  for (int e=0; e<nEtaBins; ++e) {    
    hPtResponse->Add(hResponse[e]);
    hFakeJets->Add(hFakes[e]);
  }

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

  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );              // CANVAS 0
  can->SetLogz();
  hPtResponse->Draw("COLZ");
  can->SaveAs("plots/test/pTresponse.pdf","PDF");  // SAVE 2D PT RESPONSE
  
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
    name = "plots/test/PtResponse_"; name += ptVal; name += "GeV"; name += ".pdf";
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
  can->SaveAs("plots/test/detPtResponses.pdf","PDF");

  for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
    for (int i=5; i<26; ++i) {
      int ptVal = i + 5;
      int binno = i + 1;
      if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	double wt = 1./hMissedJets->GetBinContent(binno);  // weight according to misses at part-level
	hDetWt[p]->Add(hDet[i],wt);
      }
    }    
    TString name = "plots/test/WeightedPtResponse"; name += ptBinName[p]; name += ".pdf";
    hDetWt[p]->Scale(1./hDetWt[p]->Integral());
    hDetWt[p]->Draw();
    can->SaveAs(name,"PDF");
  }


}
