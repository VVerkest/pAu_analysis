// Veronica Verkest
// February 24, 2020
// plotting macro: DifferentialUEplots.C

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  
  TFile *inFile[nEAbins], *embFile[nEAbins];

  double efficShift = 0.;  
  TString dirName = "plots/noCorrection";
  TFile* outFile = new TFile("out/noCorrection.root", "RECREATE");

  // double efficShift = 0.04;  
  // TString dirName = "plots/noCorrection/teSys1";
  // TFile* outFile = new TFile("out/noCorrection_teSys1.root", "RECREATE");

  // double efficShift = -0.04;  
  // TString dirName = "plots/noCorrection/teSys2";
  // TFile* outFile = new TFile("out/noCorrection_teSys2.root", "RECREATE");

  TH2D *hUE2D[nEAbins][55];
  TH3D *hUE3D[nEAbins][nEtaBins];
  TH3D *hUE3Dsum[nEAbins];
  TH1D *hLead[nEAbins][nEtaBins];
  TH1D *hLeadSum[nEAbins];
  TH2D *hResponse[nEAbins][nEtaBins];
  TH2D *hResponseSum[nEAbins];
  TH1D *hFakes[nEAbins][nEtaBins];
  TH1D *hFakesSum[nEAbins];
  TH1D *hMisses[nEAbins];

  
  for (int a=0; a<nEAbins; ++a) {

    // name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
    // name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_leadPtCorrected.root"; 
    name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_leadPtUncorrected.root"; 
    
    inFile[a] = new TFile(name, "READ");
    name = "hUE3Dsum_" + lohi[a] + "EA";
    hUE3Dsum[a] = new TH3D(name,";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta",xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);
    name = "hLead_" + lohi[a];
    hLeadSum[a] = new TH1D(name,";leading jet p_{T} (GeV)", 55,4.0,59.0);
    name = "hResponseSum_" + lohi[a] + "EA";
    hResponseSum[a] = new TH2D(name,";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);
    name = "hFakesSum_" + lohi[a] + "EA";
    hFakesSum[a] = new TH1D(name,";det-level leading jet p_{T} (GeV)", 55,4.0,59.0);

    for (int e=0; e<nEtaBins; ++e) {
      name = "hChgUE_" + emw[e] + "EtaJet";
      hUE3D[a][e] = (TH3D*)inFile[a]->Get(name);
      name += "_" + lohi[a];
      hUE3D[a][e]->SetName(name);
      hUE3Dsum[a]->Add(hUE3D[a][e]);

      name = "hLeadPt_" + emw[e] + "EtaJet";
      hLead[a][e] = (TH1D*)inFile[a]->Get(name);
      name += "_" + lohi[a];
      hLead[a][e]->SetName(name);
      hLeadSum[a]->Add(hLead[a][e]);
    }


    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    embFile[a] = new TFile(name, "READ");

    hMisses[a] = (TH1D*)embFile[a]->Get("hMisses");

    for (int e=0; e<nEtaBins; ++e) {
      name = "hPtResponse_" + emw[e] + "EtaJet";
      hResponse[a][e] = (TH2D*)embFile[a]->Get(name);
      name = "hResponse_" + lohi[a] + "EA_" + emw[e] + "Jet";
      hResponse[a][e]->SetName(name);
      hResponseSum[a]->Add(hResponse[a][e]);

      name = "hFakes_" + emw[e] + "EtaJet";
      hFakes[a][e] = (TH1D*)embFile[a]->Get(name);
      name += "_" + lohi[a];
      hFakes[a][e]->SetName(name);
      hFakesSum[a]->Add(hFakes[a][e]);
    }

  }

  TH1D *hMatched_part[nEAbins], *hMatched_det[nEAbins];//
  for (int a=0; a<nEAbins; ++a) {
    hMatched_part[a] = (TH1D*)hResponseSum[a]->ProjectionX();
    name = "hMatched_part_" + lohi[a] + "EA";
    hMatched_part[a]->SetName(name);
    hMatched_part[a]->Scale(hMatched_part[a]->GetEntries()/hMatched_part[a]->Integral()); // normalize to Nmatched

    hMatched_det[a] = (TH1D*)hResponseSum[a]->ProjectionY();
    name = "hMatched_det_" + lohi[a] + "EA";
    hMatched_det[a]->SetName(name);
    hMatched_det[a]->Scale(hMatched_det[a]->GetEntries()/hMatched_det[a]->Integral()); // normalize to Nmatched
  }



  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 1GeV BINS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  //  PERFORM 2D TRACKING EFFICINECY CORRECTION
  TH2D* hUE2D_detCorr[nEAbins][55];
  TFile *efficFile = new TFile("src/trackeffic_oct20.root","READ");

  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) {
      int plo=jp+4;  int phi=jp+5;
      int binno = hLeadSum[a]->FindBin(plo);
      hUE3Dsum[a]->GetXaxis()->SetRange(binno,binno);
      hUE2D[a][jp] = (TH2D*)hUE3Dsum[a]->Project3D("ZY");  // UE PT IS ON X-AXIS
      name = "hUE2Dsum_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D[a][jp]->SetName(name);
      if (hLeadSum[a]->GetBinContent(binno)!=0) { hUE2D[a][jp]->Scale(1./hLeadSum[a]->GetBinContent(binno)); }// NORMALIZE TO NJETS
      hUE3Dsum[a]->GetXaxis()->SetRange(1,-1);

      name = "hUE2Dcorr_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D_detCorr[a][jp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
      name = dirName + "/";
    }
    TrackingEfficiencyByPtAndEta55( hUE2D[a], hUE2D_detCorr[a], efficFile, lohi[a], name, efficShift );
  }


  TH2D *hUE2D_partCorr[nEAbins][nPtBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      name = "hUE2Ddet_" + lohi[a] + ptBinName[p];
      hUE2D_partCorr[a][p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
    }
  }


  for (int a=0; a<nEAbins; ++a) {  // SUM 2D HISTOGRAMS (AND ACCOUNT FOR MISSED JETS)
    for (int pp=0; pp<20; ++pp) {
      int plo = pp + 10;
      int phi = pp + 11;
      int binno = hLeadSum[a]->FindBin(plo);

      int pval = 99;
      for (int p=0; p<20; ++p) {
	if (plo>=ptLo[p] && phi<=ptHi[p]) { pval = p; }
      }

      int lo_int_bin = hLeadSum[a]->GetBin(ptLo[pval]);
      int hi_int_bin = hLeadSum[a]->GetBin(ptHi[pval]) - 1;
      
      double weight = hLeadSum[a]->GetBinContent(plo)/hLeadSum[a]->Integral( lo_int_bin, hi_int_bin);
      hUE2D_partCorr[a][pval]->Add(hUE2D_detCorr[a][pp],weight);
    }
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 1GeV BINS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~





  
  //  PROJECT HISTOGRAMS BY UE ETA
  TH1D *hUE1D_part[nEAbins][nPtBins][nEtaBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) {
	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e];
	hUE2D_partCorr[a][p]->GetYaxis()->SetRangeUser(etaLo[e],etaHi[e]);
	hUE1D_part[a][p][e] = (TH1D*)hUE2D_partCorr[a][p]->ProjectionX(name);
	cout<<hUE1D_part[a][p][e]->GetMean(1)<<" \t"<<hUE1D_part[a][p][e]->Integral()/area[e]<<endl;
	hUE2D_partCorr[a][p]->GetYaxis()->SetRangeUser(etaLo[0],etaHi[2]);
      }
    }
  }

  //  DIVIDE HISTOGRAMS BY BIN WIDTHS
  TH1D *hUE1D_part_dbw[nEAbins][nPtBins][nEtaBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) {

      }
    }
  }


  

  outFile->cd();

  for (int a=0; a<nEAbins; ++a) {
    hMatched_part[a]->Write();
    hMatched_det[a]->Write();
    hUE3Dsum[a]->Write();  // int over nLead per area is good
    hLeadSum[a]->Write();
    hResponseSum[a]->Write();
    hFakesSum[a]->Write();
  }

  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) { hUE2D[a][jp]->Write(); }  // integral is good
  }

  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) { hUE2D_detCorr[a][jp]->Write(); }  // integral is good
  }
    
  for (int a=0; a<nEAbins; ++a) {
    for (int pp=0; pp<nPtBins; ++pp) { hUE2D_partCorr[a][pp]->Write();}  // integral is good
  }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int p=0; p<nPtBins; ++p) { hUE2D_partSum[a][p]->Write(); }
  // }

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) { hUE1D_part[a][p][e]->Write(); }
    }
  }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int p=0; p<nPtBins; ++p) {
  //     for (int e=0; e<nEtaBins; ++e) { hUE1D_part_dbw[a][p][e]->Write(); }
  //   }
  // }

  for (int a=0; a<nEAbins; ++a) {
    embFile[a]->Close();
    inFile[a]->Close();
  }
  outFile->Close();

  return 0;
}
