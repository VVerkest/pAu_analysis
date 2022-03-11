// Veronica Verkest
// February 24, 2020
// plotting macro: EtaDifferentialUEplots.C

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();

  const double pi = 3.14159265;
  const double AREA = 4.*(pi/3.);

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead}^{reco} < 15", "15 < p_{T,lead}^{reco} < 20",  "20 < p_{T,lead}^{reco} < 30" };

  const int ptMarker[nPtBins] = {20,21,33}; //int EAptMarker[2][nPtBins] = {{107,108,110},{20,21,33}};
  const int marker[55] = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
			   33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  
  TFile *inFile[nEAbins], *embFile[nEAbins];

  double efficShift = 0.0;
  TFile* outFile = new TFile("out/EtaDiffJEScorrection.root", "RECREATE");
  TString dirName = "plots/EtaDiffJEScorrection";

  TH3D *hUE3D[nEAbins][nEtaBins];

  TH2D *hUE2D[nEAbins][nEtaBins][55];
  TH1D *hLead[nEAbins][nEtaBins];
  TH2D *hResponse[nEAbins][nEtaBins];
  TH1D *hFakes[nEAbins][nEtaBins];
  TH1D *hMisses[nEAbins];
  TH1D *hMissPart[nEAbins][nEtaBins];

  
  for (int a=0; a<nEAbins; ++a) {
    //dirName += lohi[a];

    name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
    // name = "../out/UE/pAuHTjetUE_halfGeVbins_" + lohi[a] + "EA_leadPtUncorrected.root";
    inFile[a] = new TFile(name, "READ");
    name = "hLead_" + lohi[a];

    for (int e=0; e<nEtaBins; ++e) {
      name = "hChgUE_" + emw[e] + "EtaJet";
      hUE3D[a][e] = (TH3D*)inFile[a]->Get(name);
      name += "_" + lohi[a];
      hUE3D[a][e]->SetName(name);

      name = "hLeadPt_" + emw[e] + "EtaJet";
      hLead[a][e] = (TH1D*)inFile[a]->Get(name);
      name += "_" + lohi[a];
      hLead[a][e]->SetName(name);
    }

    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    embFile[a] = new TFile(name, "READ");

    hMisses[a] = (TH1D*)embFile[a]->Get("hMisses");
    name = "hMisses_" + lohi[a] + "EA";
    hMisses[a]->SetName(name);
    
    for (int e=0; e<nEtaBins; ++e) {
      name = "hPtResponse_" + emw[e] + "EtaJet";
      hResponse[a][e] = (TH2D*)embFile[a]->Get(name);
      name = "hResponse_" + lohi[a] + "EA_" + emw[e] + "Jet";
      hResponse[a][e]->SetName(name);

      name = "hFakes_" + emw[e] + "EtaJet";
      hFakes[a][e] = (TH1D*)embFile[a]->Get(name);
      name += "_" + lohi[a] + "EA";
      hFakes[a][e]->SetName(name);

      name = "hMissPart" + etaBinName[e] + "Jet";
      hMissPart[a][e] = (TH1D*)embFile[a]->Get(name);
      name += "_" + lohi[a] + "EA";
      hMissPart[a][e]->SetName(name);
    }

  }

  TH1D *hMatched_part[nEAbins][nEtaBins], *hMatched_det[nEAbins][nEtaBins];//
  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      hMatched_part[a][e] = (TH1D*)hResponse[a][e]->ProjectionX();
      name = "hMatched_part_" + lohi[a] + "EA";
      hMatched_part[a][e]->SetName(name);
    
      hMatched_det[a][e] = (TH1D*)hResponse[a][e]->ProjectionY();
      name = "hMatched_det_" + lohi[a] + "EA";
      hMatched_det[a][e]->SetName(name);
    }
  }

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           FAKE AND MISSED JETS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TH1D *hFakeProb[nEAbins][nEtaBins], *hMissProb[nEAbins][nEtaBins], *hPart[nEAbins][nEtaBins];//, *hMatchPlusFake[nEAbins], *hMatchPlusMiss[nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {

      name = "hFakeProb_" + lohi[a] + "EA" + etaBinName[e] + "Jet";
      hFakeProb[a][e] = new TH1D( name, "Fake Jet Probability; det-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0 );
      name = "hMissProb_" + lohi[a] + "EA" + etaBinName[e] + "Jet";
      hMissProb[a][e] = new TH1D( name, "Missed Jet Probability; part-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0);
      
      for (int jp=0; jp<55; ++jp) {
	int binno = jp+1;
	double FakeProb = hFakes[a][e]->GetBinContent(binno) / ( hMatched_det[a][e]->GetBinContent(binno) + hFakes[a][e]->GetBinContent(binno) );
	if ( isnan(FakeProb) ) { FakeProb = 0.; }
	hFakeProb[a][e]->SetBinContent( binno, FakeProb );
	double MissProb = hMissPart[a][e]->GetBinContent(binno) / ( hMatched_part[a][e]->GetBinContent(binno) + hMissPart[a][e]->GetBinContent(binno) );
	if ( isnan(MissProb) ) { MissProb = 0.; }
	hMissProb[a][e]->SetBinContent( binno, MissProb );

	name = "hPart_" + lohi[a] + "EA" + etaBinName[e] + "Jet";
	hPart[a][e] = (TH1D*)hMatched_part[a][e]->Clone(name);
	hPart[a][e]->SetName(name);
	hPart[a][e]->SetTitle("Particle-level p_T");
	hPart[a][e]->Add(hMissPart[a][e]);
      }
     

    }
  }


  
 
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                             DETECTOR-TO-PARTICLE-LEVEL
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  //  PERFORM 2D TRACKING EFFICINECY CORRECTION   hUE2D[nEAbins][nEtaBins][55]
  TH2D* hUE2D_detCorr[nEAbins][nEtaBins][55];
  TFile *efficFile = new TFile("src/trackeffic_oct20.root","READ");

  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      for (int jp=0; jp<55; ++jp) {
	int binno = jp+1;  int plo=jp+4;  int phi=jp+5;
	hUE3D[a][e]->GetXaxis()->SetRange(binno,binno);
	hUE2D[a][e][jp] = (TH2D*)hUE3D[a][e]->Project3D("ZY");  // UE PT IS ON X-AXIS
	name = "hUE2D_" + lohi[a] + etaBinName[e] + "Jet_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
	hUE2D[a][e][jp]->SetName(name);
	hUE2D[a][e][jp]->Scale(1./hUE2D[a][e][jp]->Integral());
	hUE2D[a][e][jp]->Scale(hUE2D[a][e][jp]->GetEntries()/hLead[a][e]->Integral(binno,binno));// NORMALIZE TO NJETS
	// cout<<hUE2D[a][e][jp]->Integral()<<"\t"<<hUE2D[a][e][jp]->GetEntries()<<endl;
	hUE3D[a][e]->GetXaxis()->SetRange(1,-1);
	
	name = "hUE2Dcorr_" + lohi[a] + etaBinName[e] + "Jet_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
	hUE2D_detCorr[a][e][jp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
	name = dirName + "/";
      }
      TrackingEfficiencyByPtAndEta55( hUE2D[a][e], hUE2D_detCorr[a][e], efficFile, lohi[a], name, efficShift );
    }
  }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int e=0; e<nEtaBins; ++e) {
  //     for (int jp=0; jp<55; ++jp) {
  // 	if (hUE2D_detCorr[a][e][jp]->GetEntries()>0) { cout<<hUE2D_detCorr[a][e][jp]->GetName()<<"\t"<<hUE2D_detCorr[a][e][jp]->Integral()/AREA<<endl; }
  //     }
  //   }
  // }

  
  TH1D *FC_part[nEAbins][nEtaBins][20];
  TH2D *hUE2D_part[nEAbins][nEtaBins][20];

  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      for (int pp=0; pp<20; ++pp) {
	int plo = pp+10;  int phi = pp+11;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);
      
	FC_part[a][e][pp] = GenerateFractionalContribution( hResponse[a][e], p_lo, p_hi, dirName, lohi[a] );
	name = "FC_part_" + lohi[a] + etaBinName[e] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_part";
	FC_part[a][e][pp]->SetName(name);
      
	name = "hUE2D_" + lohi[a] + etaBinName[e] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_part";
	hUE2D_part[a][e][pp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);

	// // WeightAndSumByFC2D( FC_part[a][e][pp], hUE2D_detCorr[a][e], hUE2D_part[a][e][pp] );
	WeightAndSumByFC2D_fakesCorrection( FC_part[a][e][pp], hFakeProb[a][e], hUE2D_detCorr[a][e], hUE2D_part[a][e][pp] );
      }
    }
  }



  TH2D *hUE2D_partSum[nEAbins][nEtaBins][nPtBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	name = "hUE2Dpart_" + lohi[a] + etaBinName[e] + ptBinName[p];
	hUE2D_partSum[a][e][p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
      }
    }
  }


    /*


  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 1GeV BINS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  int binRange[nPtBins+1] = {7,12,17,27};  // SUM 2D HISTOGRAMS (AND ACCOUNT FOR MISSED JETS)
  for (int a=0; a<nEAbins; ++a) {
    for (int pp=0; pp<20; ++pp) {
      int plo = pp+10;  int phi = pp+11;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);    int pval = 99;
      int binno = pp+7;
      //  cout<<binno<<endl; cout<<plo<<"-"<<phi<<endl;  cout<<hMatched_part[a]->GetXaxis()->GetBinLowEdge(binno)<<"-"<<hMatched_part[a]->GetXaxis()->GetBinLowEdge(binno+1)<<endl<<endl;    
      for (int p=0; p<nPtBins; ++p) {
  	if ( p_lo>=ptLo[p] && p_hi<=ptHi[p] ) { pval = p; }    // cout<< p_lo <<"\t"<< p_hi <<"\t \t "<< ptBinName[pval] <<endl<<endl;  // THIS HAS BEEN TESTED :)
      }



      // double weight = hMatched_part[a]->GetBinContent(hMatched_part[a]->FindBin(plo))/hMatched_part[a]->Integral( binRange[pval], binRange[pval+1]-1 ); // THIS WEIGHT COMES FROM THE CROSS SECTION

      double weight = hPart[a]->GetBinContent(hPart[a]->FindBin(plo))/hPart[a]->Integral( binRange[pval], binRange[pval+1]-1 ); // THIS IS WITH MISSED JET CORRECTION
      // // double weight = hPart[a]->GetBinContent(hPart[a]->FindBin(plo))/hPart[a]->Integral( hPart[a]->FindBin(ptLo[pval]), hPart[a]->FindBin(ptHi[pval])-1 ); // THIS IS WITH MISSED JET CORRECTION

      // cout<<hPart[a]->FindBin(ptLo[pval])<<"  \t"<<hPart[a]->FindBin(ptHi[pval])<<endl;
      //cout<<hMatched_part[a]->GetBinLowEdge(binRange[pval])<<"\t"<<hMatched_part[a]->GetBinLowEdge(binRange[pval+1])<<"\t"<<ptBinName[pval]<<"\t"<<weight<<endl;
      //cout<<plo<<"-"<<phi<<":  "<<weight<<endl;   cout<<weight<<endl;
      hUE2D_partSum[a][pval]->Add( hUE2D_part[a][pp], weight );
    }
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 1GeV BINS ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~







  
  
  //  PROJECT HISTOGRAMS BY UE ETA
  TH1D *hUE1D_part[nEAbins][nPtBins][nEtaBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      // hMatched_part[a]->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]);
      // hUE2D_partSum[a][p]->Scale(1./hMatched_part[a]->Integral());  cout<<hMatched_part[a]->Integral()<<endl;
      // hMatched_part[a]->GetXaxis()->SetRange(1,-1);
      for (int e=0; e<nEtaBins; ++e) {
	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e];
	hUE2D_partSum[a][p]->GetYaxis()->SetRangeUser(etaLo[e],etaHi[e]);
	hUE1D_part[a][p][e] = (TH1D*)hUE2D_partSum[a][p]->ProjectionX(name);
      }
    }
  }


  //  DIVIDE HISTOGRAMS BY BIN WIDTHS
  TH1D *hUE1D_part_dbw[nEAbins][nPtBins][nEtaBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) {
	name = "hUE1D_part_" + lohi[a] + "EA" + ptBinName[p] + etaBinName[e] + "_dbw";
	hUE1D_part_dbw[a][p][e] = (TH1D*) hUE1D_part[a][p][e]->Clone(name);
	hUE1D_part_dbw[a][p][e]->SetName(name);
	
	for (int i=0; i<hUE1D_part_dbw[a][p][e]->GetNbinsX(); ++i) {  // DIVIDE BIN CONTENT BY WIDTH
	  int binno = i+1;
	  double newBinContent = hUE1D_part_dbw[a][p][e]->GetBinContent(binno)/hUE1D_part_dbw[a][p][e]->GetBinWidth(binno);
	  hUE1D_part_dbw[a][p][e]->SetBinContent(binno,newBinContent);
	}
	hUE1D_part_dbw[a][p][e]->Scale(hUE1D_part[a][p][e]->Integral()/hUE1D_part_dbw[a][p][e]->Integral());  // PRESERVE INTEGRAL VALUE
      }
    }
  }
*/
  

  outFile->cd();

  for (int a=0; a<nEAbins; ++a) {
    hMisses[nEAbins]->Write();
    for (int e=0; e<nEtaBins; ++e) {
      hLead[a][e]->Write();
      hResponse[a][e]->Write();
      hFakes[a][e]->Write();
      hMissPart[a][e]->Write();
      hMatched_part[a][e]->Write();
      hMatched_det[a][e]->Write();
      hFakeProb[a][e]->Write();
      hMissProb[a][e]->Write();
      hPart[a][e]->Write();
    }
  }
  
  for (int a=0; a<nEAbins; ++a) {
    for (int e=0; e<nEtaBins; ++e) {
      for (int jp=0; jp<55; ++jp) {
	hUE2D[a][e][jp]->Write();
	hUE2D_detCorr[a][e][jp]->Write();
      }
    }
  }
  
  for (int a=0; a<nEAbins; ++a) {
    embFile[a]->Close();
    inFile[a]->Close();
  }
  outFile->Close();

  return 0;
}
