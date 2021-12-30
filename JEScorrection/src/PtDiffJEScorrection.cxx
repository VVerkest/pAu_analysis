// Veronica Verkest
// December 6, 2021
// plotting macro: DifferentialUEplots.C

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

void ReweightPrior( TH1D *Part, TString FileName, TString suffix ){

  TFile *priorFile = new TFile(FileName,"READ");
  TString name = "ratio" + suffix;
  TH1D *priorWeight = (TH1D*)priorFile->Get(name);
  
  TH1D *temp = (TH1D*)Part->Clone();
  Part->Reset();

  for (int i=1; i<Part->GetNbinsX(); ++i) {
    double wt = priorWeight->GetBinContent( priorWeight->FindBin( Part->GetBinCenter(i) ) );
    cout<<Part->GetBinCenter(i)<<" \t"<<priorWeight->GetBinLowEdge( priorWeight->FindBin( Part->GetBinCenter(i) ) )<<" \t"<<priorWeight->GetBinLowEdge( 1+priorWeight->FindBin( Part->GetBinCenter(i) ) )<<" \t"<<wt<<endl;
    Part->SetBinContent( i, temp->GetBinContent(i)*wt );
  }
  
}


int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();

  const int nFiles = 3;
  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;
  TString detSuffix[3] = { "_det_nom", "_det_TS", "_det_TU" };
  bool PRIOR;
  TString priorFileName = "SYSTEMATICS/pt_ratios_from_Isaac.root";

  // double effic_Shift[nFiles] = { 0., 0.05, -0.05};    //  TRACKING EFFICIENCY SYSTEMATIC UNCERTAINTY
  // TString dir_Name[nFiles] = { "plots/JEScorrection", "plots/JEScorrection/teSys1", "plots/JEScorrection/teSys2" };
  // TString outFileName[nFiles] = { "out/JEScorrection.root", "out/JEScorrection_teSys1.root", "out/JEScorrection_teSys2.root" };
  // PRIOR = false;

  PRIOR = true;
  double effic_Shift[nFiles] = { 0., 0., 0.};  //  JES CORRECTION SYSTEMATIC UNCERTAINTIES
  TString dir_Name[nFiles], outFileName[nFiles], priorHistName[nFiles];
  for ( int f=0; f<nFiles; ++f ) {
    dir_Name[f] =  "SYSTEMATICS/plots/JEScorrection" + detSuffix[f];
    outFileName[f] = "SYSTEMATICS/out/JEScorrection" + detSuffix[f] + ".root";
    priorHistName[f] = "ratio" + detSuffix[f];
  }
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

  for ( int f=0; f<nFiles; ++f ) {
    double efficShift = effic_Shift[f];  
    TString dirName = dir_Name[f];
    TFile* outFile = new TFile(outFileName[f], "RECREATE");
    
    TFile *inFile[nEAbins], *embFile[nEAbins];
  
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

    TH1D *hMatched_part[nEAbins], *hMatched_det[nEAbins], *hPart[nEAbins], *hDet[nEAbins], *hFakeProb[nEAbins];//
    for (int a=0; a<nEAbins; ++a) {
      hMatched_part[a] = (TH1D*)hResponseSum[a]->ProjectionX();
      name = "hMatched_part_" + lohi[a] + "EA";
      hMatched_part[a]->SetName(name);
      // hMatched_part[a]->Scale(hMatched_part[a]->GetEntries()/hMatched_part[a]->Integral()); // normalize to Nmatched

      hMatched_det[a] = (TH1D*)hResponseSum[a]->ProjectionY();
      name = "hMatched_det_" + lohi[a] + "EA";
      hMatched_det[a]->SetName(name);
      // hMatched_det[a]->Scale(hMatched_det[a]->GetEntries()/hMatched_det[a]->Integral()); // normalize to Nmatched

      name = "hPart_" + lohi[a] + "EA";
      hPart[a] = new TH1D( name, ";part-level leading jet p_{T} (GeV)", 55,4.0,59.0 );
      hPart[a]->Add( hMatched_part[a] );
      hPart[a]->Add( hMisses[a] );  // ADD IN MISSED JETS!

      name = "hDet_" + lohi[a] + "EA";
      hDet[a] = new TH1D( name, ";det-level leading jet p_{T} (GeV)", 55,4.0,59.0 );
      hDet[a]->Add( hMatched_det[a] );
      hDet[a]->Add( hFakesSum[a] );  // ADD IN FAKE JETS!

      name = "hFakeProb_" + lohi[a] + "EA";
      hFakeProb[a] = new TH1D( name, "Prob_{fake};det-level leading jet p_{T} (GeV)", 55,4.0,59.0 );
      hFakeProb[a]->Add( hFakesSum[a] );
      hFakeProb[a]->Divide( hDet[a] );
      hFakeProb[a]->Scale(1./hFakeProb[a]->Integral());
    }


    //  APPLY PRIOR WEIGHTING
    if ( PRIOR==true ) {
      for (int a=0; a<nEAbins; ++a) {
	ReweightPrior( hPart[a], priorFileName, detSuffix[f] );
      }
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


    TH1D *FC_part[nEAbins][20];
    TH2D *hUE2D_part[nEAbins][20];

    for (int a=0; a<nEAbins; ++a) {
      for (int pp=0; pp<20; ++pp) {
	int plo = pp+10;  int phi = pp+11;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);
      
	FC_part[a][pp] = GenerateFractionalContribution( hResponseSum[a], p_lo, p_hi, dirName, lohi[a] );
	name = "FC_part_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_part";
	FC_part[a][pp]->SetName(name);
      
	name = "hUE2D_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_part";
	hUE2D_part[a][pp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);

	// WeightAndSumByFC2D( FC_part[a][pp], hUE2D_detCorr[a], hUE2D_part[a][pp] );
	WeightAndSumByFC2D_fakesCorrection( FC_part[a][pp], hFakeProb[a], hUE2D_detCorr[a], hUE2D_part[a][pp] );
      }
    }


    TH2D *hUE2D_partSum[nEAbins][nPtBins];
    for (int a=0; a<nEAbins; ++a) {
      for (int p=0; p<nPtBins; ++p) {
	name = "hUE2Ddet_" + lohi[a] + ptBinName[p];
	hUE2D_partSum[a][p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
      }
    }


  
    for (int a=0; a<nEAbins; ++a) {  // SUM 2D HISTOGRAMS (AND ACCOUNT FOR MISSED JETS)
      for (int pp=0; pp<20; ++pp) {
	int plo = pp + 10;
	int phi = pp + 11;
	int binno = hPart[a]->FindBin(plo);

	int pval = 99;
	for (int p=0; p<20; ++p) {
	  if (plo>=ptLo[p] && phi<=ptHi[p]) { pval = p; }
	}

	int lo_int_bin = hPart[a]->GetBin(ptLo[pval]);
	int hi_int_bin = hPart[a]->GetBin(ptHi[pval]) - 1;
      
	double weight = hPart[a]->GetBinContent(plo)/hPart[a]->Integral( lo_int_bin, hi_int_bin);
	hUE2D_partSum[a][pval]->Add(hUE2D_part[a][pp],weight);
      }
    }
    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 1GeV BINS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~





  
    //  PROJECT HISTOGRAMS BY UE ETA
    TH1D *hUE1D_part[nEAbins][nPtBins][nEtaBins];
    for (int a=0; a<nEAbins; ++a) {
      for (int p=0; p<nPtBins; ++p) {
	for (int e=0; e<nEtaBins; ++e) {
	  name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e];
	  hUE2D_partSum[a][p]->GetYaxis()->SetRangeUser(etaLo[e],etaHi[e]);
	  hUE1D_part[a][p][e] = (TH1D*)hUE2D_partSum[a][p]->ProjectionX(name);
	  cout<<hUE1D_part[a][p][e]->GetMean(1)<<" \t"<<hUE1D_part[a][p][e]->Integral()/area[e]<<endl;
	  hUE2D_partSum[a][p]->GetYaxis()->SetRangeUser(etaLo[0],etaHi[2]);
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
      hPart[a]->Write();
    }

    for (int a=0; a<nEAbins; ++a) {
      for (int jp=0; jp<55; ++jp) { hUE2D[a][jp]->Write(); }  // integral is good
    }

    for (int a=0; a<nEAbins; ++a) {
      for (int jp=0; jp<55; ++jp) { hUE2D_detCorr[a][jp]->Write(); }  // integral is good
    }
    
    for (int a=0; a<nEAbins; ++a) {
      for (int pp=0; pp<nPtBins; ++pp) { hUE2D_partSum[a][pp]->Write();}  // integral is good
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
    outFile->Write();
    outFile->Close();

  } //end file loop
  
  return 0;
}
