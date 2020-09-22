// Veronica Verkest
// August 31, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TString EAstring = "hiEA";
  TString inFileName = "out/sim/pAu2015embedding_" + EAstring + ".root";
  TString UEfileName = "../out/UE/pAuHTjetUE_" + EAstring + ".root";
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString name, title;
  TString directory = "plots/test/" + EAstring + "/";
  TH2D *hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.5,59.5, 55,4.5,59.5);
  TH1D *hFakeJets = new TH1D("hFakeJets",";fake det-level leading jet p_{T} (GeV)",55,4.5,59.5);
  TH1D *hDet[21];
  TH1D *hDetWt[nPtBins];
  TH1D *hLeadJetPt[nEtaBins];
  TH3D *hUE3D[nEtaBins];
  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", 55,4.5,59.5, 30,0.0,30.0, 40,-1.0,1.0);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.5,59.5);
  TH2D *hChgUE2D_corr[nEffEtaBins];
  TH1D *hUEpt[nEtaBins][55]; 
  TH1D *hWtUEpt[nEtaBins][nPtBins];
  TH2D *hChgUE2D[nEffEtaBins];
  TH2D *hAddedChgUE2D_corr[nEtaBins];
  
  for ( int f=0; f<nEffEtaBins; ++f) {
    name = "hChgUE2D_eta"; name += f;
    hChgUE2D[f] = new TH2D( name ,";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV)", 55,4.5,59.5, 30,0.0,30.0);
    name = "hChgUE2D_corr_eta"; name += f;
    hChgUE2D_corr[f] = new TH2D( name ,";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)", 55,4.5,59.5, 30,0.0,30.0);
  }
  for (int e=0; e<nEtaBins; ++e) {
    for (int p=0; p<nPtBins; ++p) {
      name = "hWtUEpt"; name += ptBinName[p]; name += "_"; name += emw[e];
      title = ptBinString[p]; title += "    "; title += UEetaBinString[e]; title += ";ch UE p_{T} (GeV)";
      hWtUEpt[e][p] = new TH1D(name,title,30,0.0,30.0);
    }
    name = "hAddedChgUE2D_corr" + etaBinName[e] ;
    hAddedChgUE2D_corr[e] = new TH2D( name ,";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)", 55,4.5,59.5, 30,0.0,30.0);
  }
 
  ////////////////////////////////// EMBEDDING FILE FOR LEADING JET PT CORRECTION //////////////////////////////////
  TFile *inFile = new TFile( inFileName ,"READ");  // OPEN RESPONSE FILE AND COLLECT HISTOS

  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  GetEmbeddingHistograms( inFile, hPtResponse, hFakeJets, hMissedJets, directory );

  ProjectPartLevelJetPt( hPtResponse, hDet, directory );

  GenerateWeightedPtResponse( hDetWt, hDet, hMissedJets, directory );

  ////////////////////////////////////////// pAu DATA FILE FOR UE HISTOGRAMS //////////////////////////////////////////
  TFile *UEfile = new TFile(UEfileName,"READ");  // OPEN pAu FILE

  for (int e=0; e<nEtaBins; ++e) {
    name = "hLeadPt" + etaBinName[e] + "Jet";    //  GATHER HISTOGRAMS OF UE AND LEAD PT
    hLeadJetPt[e] = (TH1D*)UEfile->Get(name);
    name = "hChgUE" + etaBinName[e] + "Jet";
    hUE3D[e] = (TH3D*)UEfile->Get(name);    
  }

  for (int e=0; e<nEtaBins; ++e) {
    hChgUE3D->Add(hUE3D[e]);
    hleadPt->Add(hLeadJetPt[e]);
  }
  
  for (int f=0; f<nEffEtaBins; ++f ) { hChgUE2D[f] = (TH2D*)ProjectUEHistograms( hChgUE3D, f, directory ); }
  ////////////////////////////////// TRACKING EFFICIENCY FILE FOR UE pT CORRECTION //////////////////////////////////
  TFile *ef = new TFile( "src/trackeffic.root", "READ" );

  for (int f=1; f<=nEffEtaBins; ++f ) {
    name = "eff_s_bin_" + GetEfficHistoName(EAstring) + "_bbc__"; name += f; name += "_"; name += f; name += "_eta";
    TrackingEfficiency2DCorrection( hChgUE2D_corr[f-1], hChgUE2D[f-1], ef, name );
  }

  AddCorrected2DHistograms( hChgUE2D_corr, hAddedChgUE2D_corr );
  
  TCanvas *ctemp = new TCanvas( "ctemp" , "" ,700 ,500 );
  ctemp->SetLogz();
  
  for (int e=0; e<nEtaBins; ++e ) {
    ctemp->cd();
    hAddedChgUE2D_corr[e]->SetAxisRange( 0.0,15.0,"Y");
    hAddedChgUE2D_corr[e]->Draw("COLZ");
    name = directory; name += hAddedChgUE2D_corr[e]->GetName(); name += ".pdf";
    ctemp->SaveAs(name,"PDF");
    hAddedChgUE2D_corr[e]->GetYaxis()->SetRangeUser(0.0,30.0);
    cout<<hAddedChgUE2D_corr[e]->ProjectionX()->Integral()/(hleadPt->Integral()*area[e])<<endl;
    
    ProjectAndScaleUEHistogramForAllPt( hAddedChgUE2D_corr[e], hleadPt, hUEpt[e], e, directory );

    WeightUEPtByLeadPtAndFakes( hWtUEpt[e], hUEpt[e], hDetWt, hleadPt, hFakeJets, e, directory );
  }
  
  for (int p=0; p<nPtBins; ++p) {
    for (int e=0; e<nEtaBins; ++e ) {
      TString suf = ptBinName[p] + etaBinName[e];      
      ProjectAndSaveFinalUEPlots( hWtUEpt[e][p], suf, e, directory );
    }
  }

  return 0;

}
