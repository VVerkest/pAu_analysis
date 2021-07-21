// Veronica Verkest
// February 24, 2020
// plotting macro: DifferentialUEplots.C

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;


void WeightAndSumByFC2D( TH1D* FC, TH2D* UE2D[55], TH2D *UE2Dpart ){

  double sum = 0.;

  for (int i=0; i<FC->GetNbinsX(); ++i) {
    
    int binno = i+1;
    
    double wt = FC->GetBinContent(binno);
    // std::cout<<FC->GetBinCenter(binno)<<std::endl;
    if (wt==0) { continue; }
    UE2Dpart->Add(UE2D[i], wt);

    sum += wt;

  }

  // std::cout<<sum<<std::endl;
  
}


void WeightAndSumByFC2D_fakesCorrection( TH1D* FC, TH1D* FakeProb, TH2D* UE2D[55], TH2D *UE2Dpart ){

  double sum = 0.;
  
  for (int i=0; i<FC->GetNbinsX(); ++i) {
    int binno = i+1;
    FC->SetBinContent( binno, FC->GetBinContent(binno)*( 1. - FakeProb->GetBinContent(1+i) ) );
  }

  FC->Scale(1./FC->Integral());

  for (int i=0; i<FC->GetNbinsX(); ++i) {
    int binno = i+1;
    double wt = FC->GetBinContent(binno);
    if (wt==0) { continue; }
    UE2Dpart->Add(UE2D[i], wt);
    sum += wt;
  }
  std::cout<<sum<<std::endl;
}



int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();

  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

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
  // TString dirName = "plots/new/5gevBins";

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  
  TFile *inFile[nEAbins], *embFile[nEAbins];
  // TFile* outFile = new TFile("out/5gevBins.root", "RECREATE");
  TFile* outFile = new TFile("out/JEScorrection.root", "RECREATE");
  TString dirName = "plots/JEScorrection";
 
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
    //dirName += lohi[a];

    name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
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
    name = "hMisses_" + lohi[a] + "EA";
    hMisses[a]->SetName(name);
    
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

    hMatched_det[a] = (TH1D*)hResponseSum[a]->ProjectionY();
    name = "hMatched_det_" + lohi[a] + "EA";
    hMatched_det[a]->SetName(name);
  }



  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           FAKE AND MISSED JETS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TH1D *hFakeProb[nEAbins], *hMissProb[nEAbins], *hPart[nEAbins];//, *hMatchPlusFake[nEAbins], *hMatchPlusMiss[nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    name = "hFakeProb_" + lohi[a] + "EA";
    hFakeProb[a] = new TH1D( name, "Fake Jet Probability; det-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0 );
    name = "hMissProb_" + lohi[a] + "EA";
    hMissProb[a] = new TH1D( name, "Missed Jet Probability; part-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0);
    
    for (int jp=0; jp<55; ++jp) {
      int binno = jp+1;
      double FakeProb = hFakesSum[a]->GetBinContent(binno) / ( hMatched_det[a]->GetBinContent(binno) + hFakesSum[a]->GetBinContent(binno) );
      if ( isnan(FakeProb) ) { FakeProb = 0.; }
      hFakeProb[a]->SetBinContent( binno, FakeProb );
      double MissProb = hMisses[a]->GetBinContent(binno) / ( hMatched_part[a]->GetBinContent(binno) + hMisses[a]->GetBinContent(binno) );
      if ( isnan(MissProb) ) { MissProb = 0.; }
      hMissProb[a]->SetBinContent( binno, MissProb );
    }

    name = "hPart_" + lohi[a] + "EA";
    title = "Particle-level p_{#mathrm{T}}; part-level leading jet p_{#mathrm{T}} [GeV]";
    hPart[a] = (TH1D*)hMatched_part[a]->Clone(name);
    hPart[a]->SetNameTitle(name,title);
    hPart[a]->Add(hMisses[a]);
    // hPart[a]->Multiply(hMissProb[a]);
    // hPart[a]->Add(hMatched_part[a]);
  }


  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                             DETECTOR-TO-PARTICLE-LEVEL
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  //  PERFORM 2D TRACKING EFFICINECY CORRECTION
  TH2D* hUE2D_detCorr[nEAbins][55];
  TFile *efficFile = new TFile("src/trackeffic.root","READ");

  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) {
      int binno = jp+1;  int plo=jp+4;  int phi=jp+5;
      hUE3Dsum[a]->GetXaxis()->SetRange(binno,binno);
      hUE2D[a][jp] = (TH2D*)hUE3Dsum[a]->Project3D("ZY");  // UE PT IS ON X-AXIS
      name = "hUE2Dsum_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D[a][jp]->SetName(name);
      hUE2D[a][jp]->Scale(1./hUE2D[a][jp]->Integral());
      hUE2D[a][jp]->Scale(hUE2D[a][jp]->GetEntries()/hLeadSum[a]->Integral(binno,binno));// NORMALIZE TO NJETS
      // cout<<hUE2D[a][jp]->Integral()<<"\t"<<hUE2D[a][jp]->GetEntries()<<endl;
      hUE3Dsum[a]->GetXaxis()->SetRange(1,-1);

      name = "hUE2Dcorr_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D_detCorr[a][jp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
      name = dirName + "/";
    }
    TrackingEfficiencyByPtAndEta55( hUE2D[a], hUE2D_detCorr[a], efficFile, lohi[a], name );
  }
  
  // for (int a=0; a<nEAbins; ++a) {
  //   for (int jp=0; jp<55; ++jp) {
  //     if (hUE2D_detCorr[a][jp]->GetEntries()>0) { cout<<hUE2D_detCorr[a][jp]->GetName()<<"\t"<<hUE2D_detCorr[a][jp]->Integral()/AREA<<endl; }
  //   }
  // }


  
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
      name = "hUE2Dpart_" + lohi[a] + ptBinName[p];
      hUE2D_partSum[a][p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
    }
  }





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
      // double weight = hPart[a]->GetBinContent(hPart[a]->FindBin(plo))/hPart[a]->Integral( hPart[a]->FindBin(ptLo[pval]), hPart[a]->FindBin(ptHi[pval])-1 ); // THIS IS WITH MISSED JET CORRECTION

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








  

  outFile->cd();

  for (int a=0; a<nEAbins; ++a) {
    hMatched_part[a]->Write();
    hMatched_det[a]->Write();
    hMisses[a]->Write();
    hFakesSum[a]->Write();
    hUE3Dsum[a]->Write();  // int over nLead per area is good
    hLeadSum[a]->Write();
    hResponseSum[a]->Write();
    hFakeProb[a]->Write();
    hMissProb[a]->Write();
  }

  c0->SetLogy(0);
  c0->SetLogz();
  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) {
      hUE2D[a][jp]->Write();
      name = dirName + "/" + hUE2D[a][jp]->GetName() + ".pdf";
      hUE2D[a][jp]->Draw("COLZ");
      c0->SaveAs(name,"PDF");
    }  // integral is good
  }

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<55; ++p) { hUE2D_detCorr[a][p]->Write(); }
  }
  
  for (int a=0; a<nEAbins; ++a) {
    for (int pp=0; pp<20; ++pp) { FC_part[a][pp]->Write(); }
  }
  
  for (int a=0; a<nEAbins; ++a) {
    for (int pp=0; pp<20; ++pp) {
      hUE2D_part[a][pp]->Write();
      name = dirName + "/" + hUE2D_part[a][pp]->GetName() + ".pdf";
      hUE2D_part[a][pp]->Draw("COLZ");
      c0->SaveAs(name,"PDF");
    }  // 1GeV part pT bins
  }

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      hUE2D_partSum[a][p]->Write();
      name = dirName + "/" + hUE2D_partSum[a][p]->GetName() + ".pdf";
      hUE2D_partSum[a][p]->Draw("COLZ");
      c0->SaveAs(name,"PDF");
    }
  }

  c0->SetLogy();
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) {
	hUE1D_part[a][p][e]->Write();
	name = dirName + "/" + hUE1D_part[a][p][e]->GetName() + ".pdf";
	hUE1D_part[a][p][e]->Draw("COLZ");
	c0->SaveAs(name,"PDF");
	// cout<<hUE1D_part[a][p][e]->Integral()<<endl;
      }
    }
  }

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) { hUE1D_part_dbw[a][p][e]->Write(); }
    }
  }

  for (int a=0; a<nEAbins; ++a) {
    embFile[a]->Close();
    inFile[a]->Close();
  }
  outFile->Close();

  return 0;
}
