// Veronica Verkest
// October 5, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString name, title;
  TFile *outFile = new TFile("out/pAuUE_pt_halfGeVbins.root","RECREATE");

  const int n_bins = 3;
  double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  const int n_ybins = 3;
  double y_bin_edge[n_ybins+1] = { 0.55,0.65,0.75,0.85 };
  TH2D *meanPt_hscale = new TH2D("meanPt_hscale",";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge);
  double y_bin_edge2[n_ybins+1] = { 0.5,1.2,1.5,1.8 };
  TH2D *nCh_hscale = new TH2D("nCh_hscale",";(GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge2);
  TH1D *hUE[nPtBins][nEtaBins][nEAbins];


  
  for (int a=0; a<nEAbins; ++a) {

    TH2D* hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);
    TH1D* hFakeJets = new TH1D("hFakeJets",";missing part-level leading jet p_{T} (GeV)",55,4.0,59.0);
    TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);
    TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.0,59.0);
    TH1D *hDet[21], *hDetWt[nPtBins], *hLeadJetPt[nEtaBins], *hUEpt[55], *hWtUEpt[nPtBins][nEtaBins], *hPartJetUE[nPtBins][nEtaBins];
    TH2D *hChgUE2D[55], *hChgUE2D_corr[nPtBins], *hAddedChgUE2D[nPtBins];
    TH3D *hUE3D[nEtaBins];
    
    for (int p=0; p<nPtBins; ++p) {
      name = "hChgUE2D" + ptBinName[p];
      hAddedChgUE2D[p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);

      name = "hChgUE2D"+ ptBinName[p] + "_corr";
      hChgUE2D_corr[p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);

      for (int e=0; e<nEtaBins; ++e) {
	name = "hWtUEpt"; name += ptBinName[p] + etaBinName[e];
	TString title = ptBinString[p]; title += ";ch UE p_{T} (GeV)";
	hWtUEpt[e][p] = new TH1D(name,title,ybins,ybinEdge);
	name = "hPartJetUE"; name += ptBinName[p] + etaBinName[e];
	title = ";ch UE p_{T} (GeV)";
	hPartJetUE[e][p] = new TH1D(name,title,ybins,ybinEdge);
      }
    }

    
    // TString directory = "plots/JES/" + lohi[a] + "/";
    TString directory = "plots/JES/halfGeVbins/" + lohi[a] + "/";

    ////////////////////////////////// EMBEDDING FILE FOR LEADING JET PT CORRECTION //////////////////////////////////
    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    TFile *inFile = new TFile(name,"READ");  // OPEN RESPONSE FILE AND COLLECT HISTOS
    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

    GetEmbeddingHistograms( inFile, hPtResponse, hFakeJets, hMissedJets, directory );
    ProjectPartLevelJetPt( hPtResponse, hDet, directory );
    GenerateWeightedPtResponse( hDetWt, hDet, hMissedJets, directory );

    ////////////////////////////////////////// pAu DATA FILE FOR UE HISTOGRAMS //////////////////////////////////////////
    // TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root","READ");  // OPEN pAu FILE & GET HISTOGRAMS OF UE & LEAD PT
    TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_halfGeVbins_" + lohi[a] + "EA.root","READ");  // OPEN pAu FILE & GET HISTOGRAMS OF UE & LEAD PT
    for (int e=0; e<nEtaBins; ++e) {
      name = "hLeadPt" + etaBinName[e] + "Jet";
      hLeadJetPt[e] = (TH1D*)UEfile->Get(name);
      hleadPt->Add(hLeadJetPt[e]);
                                                     // ";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta"
      name = "hChgUE" + etaBinName[e] + "Jet";  // X: lead pT
      hUE3D[e] = (TH3D*)UEfile->Get(name);      // Y: UE pT
      hChgUE3D->Add(hUE3D[e]);                  // Z: UE eta
    }


    // for (int p=0; p<nPtBins; ++p) {
    //   hChgUE3D->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]);
    //   hChgUE2D_corr[p] = (TH2D*) hChgUE3D->Project3D("ZY");
    //   // int ptlobin = hleadPt->FindBin( ptLo[p] + 0.5 );
    //   // int pthibin = hleadPt->FindBin( ptHi[p] - 0.5 );
    //   // cout<<ptBinName[p]<<":  "<<ptlobin<<"  "<<pthibin<<endl;
    //   // cout<<hleadPt->GetBinLowEdge(ptlobin)<<"  "<<hleadPt->GetBinLowEdge(pthibin)+hleadPt->GetBinWidth(pthibin)<<endl;
    //   // hChgUE2D_corr[p]->Scale(1./hleadPt->Integral(ptlobin,pthibin));
    //   // for (int e=0; e<nEtaBins; ++e) { ProjectAndPlotByEta( hChgUE2D_corr[p], hPartJetUE[e][p], e, p, a, directory, outFile ); }
    // }

    TCanvas * can1 = new TCanvas( "can1" , "" ,700 ,500 );
    can1->SetLogy();
    title = directory + "LeadPt_" + lohi[a] + "EA.pdf";
    hleadPt->Draw();
    can1->SaveAs(title,"PDF");
    
    can1->SetLogy(0);
    can1->SetLogz();
    for (int i=0; i<55; ++i) {
      int ptLo = i + 4;
      int binno = i + 1;

      hChgUE3D->GetXaxis()->SetRange(binno,binno);
      hChgUE2D[i] = (TH2D*)hChgUE3D->Project3D("ZY");

      if (hChgUE2D[i]->GetEntries()==0){ continue;}

      double nJets = hleadPt->Integral(binno,binno);
      double realRate = 1.-hFakeJets->GetBinContent(binno);
      hChgUE2D[i]->Scale( realRate/nJets );

      name = "hChgUE2D_"; name += ptLo; name += "to"; name += ptLo+1; name += "GeVdetJets";
      hChgUE2D[i]->SetName(name);
      hChgUE2D[i]->Draw("COLZ");
    
      name = directory + name + ".pdf";
      can1->SaveAs(name,"PDF");
    }

    ////////////////////////////////////////// EFFICIENCY FILES FOR UE TRACKING CORRECTION //////////////////////////////////////////
    TFile *efficFile = new TFile("../embedding/src/trackeffic_20etaBins.root","READ");
    
    for (int p=0; p<nPtBins; ++p) { WeightAndAddCorrected2Ds( hAddedChgUE2D[p], hDetWt[p], hChgUE2D, directory ); }

    TrackingEfficiencyByPtAndEta( hAddedChgUE2D, hChgUE2D_corr, efficFile, lohi[a], directory );

    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) { ProjectAndPlotByEta( hChgUE2D_corr[p], hPartJetUE[e][p], e, p, a, directory, outFile ); }
    }

    hFakeJets->Delete();
    hMissedJets->Delete();
  }
  
  outFile->Close();
  return 0;
}
