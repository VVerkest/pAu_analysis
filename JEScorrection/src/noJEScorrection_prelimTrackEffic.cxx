// Veronica Verkest
// February 24, 2020
// plotting macro: DifferentialUEplots.C

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;


void WeightAndSumByFC2D( TH1D* FC, TH2D* UE2D[55], TH2D *UE2Dpart ){

  for (int i=0; i<FC->GetNbinsX(); ++i) {
    
    int binno = i+1;
    
    double wt = FC->GetBinContent(binno);
    if (wt==0) { continue; }
    UE2Dpart->Add(UE2D[i], wt);

  }
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
  // TString dirName = "plots/noCorrection";
  TString dirName = "plots/noCorrection_halfGeVbins_prelimTrackEffic";

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  
  TFile *inFile[nEAbins], *embFile[nEAbins];
  // TFile* outFile = new TFile("out/noCorrection_1GeVbins.root", "RECREATE");
  TFile* outFile = new TFile("out/noCorrection_halfGeVbins_prelimTrackEffic.root", "RECREATE");

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

    // name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
    name = "../out/UE/pAuHTjetUE_halfGeVbins_" + lohi[a] + "EA_uncorrected.root";
    // name = "../out/UE/pAuHTjetUE_halfGeVbins_" + lohi[a] + "EA_leadPtUncorrected.root";
    
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
  TFile* efficFile[nEAbins];
  TString efficFileName[nEAbins] = { "../src/trackeffic_loEA.root", "../src/trackeffic_hiEA.root" };

  const int eta_bins = 10;
  double etabinEdge[eta_bins+1] = { -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0 };

  for (int a=0; a<nEAbins; ++a) {
    efficFile[a] = new TFile( efficFileName[a], "READ" );
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
    // TrackingEfficiencyByPtAndEta55( hUE2D[a], hUE2D_detCorr[a], efficFile, lohi[a], name );
  }


  TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  TH1D *hEff[eta_bins][nEAbins];
  TF1* efficFit[eta_bins][nEAbins];

  TCanvas *can = new TCanvas;  

  for (int i=0; i<eta_bins; ++i) {
    for (int j=0; j<nEAbins; ++j){    
      int binno = i+1;
    
      if (j==0 ) { name = "eff_s_bin_1_3_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
      else if (j==1) {name = "eff_s_bin_7_10_bbc__"; name += binno; name += "_"; name += binno; name += "_eta"; }
    
      hEff[i][j] = (TH1D*)efficFile[j]->Get(name);

      name = EAbinName[j] + "EA_effic_etaBin"; name += i;

      hEff[i][j]->SetAxisRange(0.0,1.3,"Y");
      hEff[i][j]->SetAxisRange(0.0,15.0,"X");
      hEff[i][j]->Fit( "eff", "EM" );

      efficFit[i][j] = hEff[i][j]->GetFunction("eff");
      efficFit[i][j]->SetName(name);
      
      hEff[i][j]->Draw();
      
      saveName = dirName + "/" + name + ".pdf";
      can->SaveAs(saveName,"PDF");

    }
  }

  



  //  TRACKING EFFICIENCY CORRECTIONS!
  for (int a=nEAbins-1; a>=0; --a) {
    for (int p=0; p<55; ++p) {

      for (int iy=0; iy<hUE2D_detCorr[a][p]->GetNbinsY(); ++iy) {  // loop over 20 eta bins

  	int fileNo = 99;
  	int ybin = iy+1;
	
  	double binCenter = hUE2D_detCorr[a][p]->GetYaxis()->GetBinCenter( ybin );

	for (int j=0; j<eta_bins; ++j) {
	  if ( binCenter>etabinEdge[j] && binCenter<etabinEdge[j+1] ) { fileNo = j; };
	}
  	if ( fileNo==99 ) { cerr<<"CANNOT FIND EFFICIENCY HISTOGRAM"<<endl; }

	// cout<<fileNo<<"  "<<binCenter<<"  "<<hEff[fileNo][a]->GetTitle()<<endl;

  	for(int ix = 1; ix <= hUE2D[a][p]->GetNbinsX(); ++ix){
  	  int xbin = ix;
	  
  	  double pt = hUE2D[a][p]->GetXaxis()->GetBinCenter(xbin);
  	  // if ( pt > 3.0 ) { pt = 3.0; }
  	  double ptbin = hEff[a][fileNo]->FindBin(pt);
	  
  	  // double eff = hEff[a][fileNo]->GetBinContent(ptbin);
	  double eff = efficFit[a][fileNo]->Eval( hUE2D[a][p]->GetXaxis()->GetBinCenter(xbin) );
	  // double eff_err = efficFit[a][fileNo]->Eval( hUE2D[a][p]->GetXaxis()->GetBinCenter(xbin) );
	  
	  double old_value = hUE2D[a][p]->GetBinContent(xbin,ybin);
  	  double corr_value = (double) old_value/eff;
  	  double old_err = hUE2D[a][p]->GetBinError(xbin,ybin);
	  double relativeError = sqrt( (old_err*old_err) + ( hEff[a][fileNo]->GetBinError(ptbin)*hEff[a][fileNo]->GetBinError(ptbin)) );
  	  double corr_err = corr_value*relativeError;
  	  hUE2D_detCorr[a][p]->SetBinContent( xbin, ybin, corr_value);
  	  hUE2D_detCorr[a][p]->SetBinError( xbin, ybin, corr_err );
  	}

	
      }
	
    }
  }



  




  



  TH2D *hUE2D_partCorr[nEAbins][nPtBins];
  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      name = "hUE2Ddet_" + lohi[a] + ptBinName[p];
      hUE2D_partCorr[a][p] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
    }
  }
  
  int binRange[nPtBins+1] = {7,12,17,27};  // SUM 2D HISTOGRAMS (AND ACCOUNT FOR MISSED JETS)
  for (int a=0; a<nEAbins; ++a) {
    for (int pp=0; pp<20; ++pp) {
      int plo = pp+10;  int phi = pp+11;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);    int pval = 99;
      int binno = pp+7;
      //  cout<<binno<<endl; cout<<plo<<"-"<<phi<<endl;  cout<<hLeadSum[a]->GetXaxis()->GetBinLowEdge(binno)<<"-"<<hLeadSum[a]->GetXaxis()->GetBinLowEdge(binno+1)<<endl<<endl;    
      for (int p=0; p<nPtBins; ++p) {
  	if ( p_lo>=ptLo[p] && p_hi<=ptHi[p] ) { pval = p; }    // cout<< p_lo <<"\t"<< p_hi <<"\t \t "<< ptBinName[pval] <<endl<<endl;  // THIS HAS BEEN TESTED :)
      }
      double weight = hLeadSum[a]->Integral(hLeadSum[a]->FindBin(plo),hLeadSum[a]->FindBin(plo))/hLeadSum[a]->Integral( binRange[pval], binRange[pval+1]-1 );
      // cout<<hLeadSum[a]->GetBinLowEdge(binRange[pval])<<"\t"<<hLeadSum[a]->GetBinLowEdge(binRange[pval+1])<<"\t"<<ptBinName[pval]<<"\t"<<weight<<endl;
      // cout<< binRange[pval] <<"\t"<< binRange[pval+1]-1 <<endl;
      // cout<<plo<<endl;//"-"<<phi<<":  "<<weight<<endl;   cout<<weight<<endl;
      hUE2D_partCorr[a][pval]->Add( hUE2D_detCorr[a][pp], weight );
      // hUE2D_partCorr[a][pval]->Add( hUE2D[a][pp], weight );  // WITHOUT TRACKING EFFICIENCY CORRECTION
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
	hUE2D_partCorr[a][p]->GetYaxis()->SetRangeUser(etaLo[e],etaHi[e]);
	hUE1D_part[a][p][e] = (TH1D*)hUE2D_partCorr[a][p]->ProjectionX(name);
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
    hUE3Dsum[a]->Write();  // int over nLead per area is good
    hLeadSum[a]->Write();
    hResponseSum[a]->Write();
    hFakesSum[a]->Write();
  }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int jp=0; jp<55; ++jp) { hUE2D[a][jp]->Write(); }  // integral is good
  // }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int pp=0; pp<20; ++pp) { hUE2D_part[a][pp]->Write();}  // integral is good
  // }

  // for (int a=0; a<nEAbins; ++a) {
  //   for (int p=0; p<nPtBins; ++p) { hUE2D_partSum[a][p]->Write(); }
  // }

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {
      for (int e=0; e<nEtaBins; ++e) { hUE1D_part[a][p][e]->Write(); }
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
