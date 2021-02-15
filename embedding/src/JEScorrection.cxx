// Veronica Verkest
// October 5, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString name, title;
  TString directory = "plots/JES/";

  TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);
  TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.0,59.0);
  TH1D *hDet[21], *hDetWt[nPtBins], *hLeadJetPt[nEtaBins], *hUEpt[55], *hWtUEpt[nPtBins][nEtaBins];
  TH3D *hUE3D[nEtaBins];
  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.0,59.0);
  TH2D *hChgUE2D[55], *hChgUE2D_corr[55], *hAddedChgUE2D_corr[nEtaBins];

  for (int p=0; p<nPtBins; ++p) {
    for (int e=0; e<nEtaBins; ++e) {
      name = "hWtUEpt"; name += ptBinName[p] + etaBinName[e];
      TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
      hWtUEpt[e][p] = new TH1D(name,title,ybins,ybinEdge);
    }
  }

      
  ////////////////////////////////// EMBEDDING FILE FOR LEADING JET PT CORRECTION //////////////////////////////////
  TFile *inFile = new TFile("out/sim/pAu2015embedding_hiEA.root","READ");  // OPEN RESPONSE FILE AND COLLECT HISTOS

  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  GetEmbeddingHistograms( inFile, hPtResponse, hFakeJets, hMissedJets, directory );

  ProjectPartLevelJetPt( hPtResponse, hDet, directory );

  GenerateWeightedPtResponse( hDetWt, hDet, hMissedJets, directory );

  ////////////////////////////////////////// pAu DATA FILE FOR UE HISTOGRAMS //////////////////////////////////////////
  TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_hiEA_uncorrected.root","READ");  // OPEN pAu FILE & GET HISTOGRAMS OF UE & LEAD PT
  
  for (int e=0; e<nEtaBins; ++e) {
    name = "hLeadPt" + etaBinName[e] + "Jet";
    hLeadJetPt[e] = (TH1D*)UEfile->Get(name);
    hleadPt->Add(hLeadJetPt[e]);
    
    name = "hChgUE" + etaBinName[e] + "Jet";
    hUE3D[e] = (TH3D*)UEfile->Get(name);    
    hChgUE3D->Add(hUE3D[e]);
  }

  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
  can->SetLogz();
  for (int i=0; i<55; ++i) {
    int ptLo = i + 4;
    int binno = i + 1;

    hChgUE3D->GetXaxis()->SetRange(binno,binno);
    hChgUE2D[i] = (TH2D*)hChgUE3D->Project3D("ZY");

    if (hChgUE2D[i]->GetEntries()==0){ continue;}
    
    name = "hChgUE2D_"; name += ptLo; name += "to"; name += ptLo+1; name += "GeVdetJets";
    hChgUE2D[i]->SetName(name);
    hChgUE2D[i]->Draw("COLZ");

    name = directory + name + ".pdf";
    can->SaveAs(name,"PDF");
  }



  return 0; 
}
