// #include "RooUnfold.h"
// #include "RooUnfoldTUnfold.h"

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

#include "params.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TString name, title;
  
  TLatex *tex = new TLatex( 0.4, 0.95," ");     tex->SetTextFont(63);     tex->SetTextSize(30);
  tex->SetTextColor(kBlack);     tex->SetLineWidth(1);     tex->SetNDC();

  
  // OPEN RESPONSE FILE AND COLLECT HISTOS
  TFile *inFile = new TFile("out/sim/pAu2015embedding.root","READ");

  TH2D *hResponse[nEtaBins];
  TH1D *hFakeJets[nEtaBins];
  TH1D *hMissedJets = (TH1D*)inFile->Get("hMisses");

  int nAccepted = 0; int nFakes = 0;
  
  for (int e=0; e<nEtaBins; ++e) {    
    name = "hPtResponse" + etaBinName[e] + "Jet";
    hResponse[e] = (TH2D*)inFile->Get(name);

    name = "hFakes" + etaBinName[e] + "Jet";
    hFakeJets[e] = (TH1D*)inFile->Get(name);

    nAccepted += hResponse[e]->GetEntries();
    nFakes += hFakeJets[e]->GetEntries();

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
  can2->SaveAs("plots/loEA/MissedJets.pdf","PDF");
  
  for (int e=0; e<nEtaBins; ++e) {  
    scale = (double)hFakeJets[e]->GetEntries()/nEvents;
    hFakeJets[e]->Scale(scale/hFakeJets[e]->Integral());
    can2->cd();
    hFakeJets[e]->SetMarkerColor(etaColor[e]);
    hFakeJets[e]->SetLineColor(etaColor[e]);
    hFakeJets[e]->Draw();
    name = "plots/loEA/FakeJets" + etaBinName[e] + "Jet" + ".pdf";
    can2->SaveAs(name,"PDF");

    can1->cd();
    hResponse[e]->Draw("COLZ");
    name = "plots/loEA/pTresponse" + etaBinName[e] + "Jet" + ".pdf";
    can1->SaveAs(name,"PDF");  // SAVE 2D PT RESPONSE
  }

  can1->Close();
  can2->Close();

  TH1D *hDet[55][nEtaBins]; TH1D *hDetWt[nPtBins][nEtaBins];

  TCanvas * can3 = new TCanvas( "can3" , "" ,700 ,500 );              // CANVAS 3
  
  for (int e=0; e<nEtaBins; ++e) {
    
    for (int p=0; p<nPtBins; ++p) {
      name = "hResponse" + etaBinName[e] + "Jet" + ptBinName[p];
      hDetWt[p][e] = new TH1D( name, ";det-level leading jet p_{T} (GeV)",55,4.5,59.5);
      hDetWt[p][e]->SetMarkerColor(etaColor[e]);
      hDetWt[p][e]->SetLineColor(etaColor[e]);
      hDetWt[p][e]->SetMarkerStyle(ptMarker[p]);
    }
    
    for (int i=5; i<26; ++i) {  // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
      int ptVal = i + 5;
      int binno = i + 1;
      
      TString name = "hResponse_" + etaBinName[e] + ptVal + "GeV";
      hDet[i][e] = (TH1D*) hResponse[e]->ProjectionY(name,binno,binno);
      hDet[i][e]->Scale(1./hDet[i][e]->Integral());
      hDet[i][e]->SetMarkerStyle(marker[i]);
      hDet[i][e]->SetMarkerSize(1);
      hDet[i][e]->SetStats(0);
      name = ""; name += ptVal; name += " GeV part. jet";
      hDet[i][e]->SetNameTitle(name,";det-level leading jet p_{T} (GeV)");
    }

    hDet[5][e]->GetYaxis()->SetRangeUser(0.0,0.6);
    hDet[5][e]->Draw("PLC PMC");
    for (int i=6; i<26; ++i) { hDet[i][e]->Draw("SAME PLC PMC"); }
    can3->BuildLegend(0.68,0.1,0.9,0.9);

    tex->DrawLatex(0.4,0.95,etaBinString[e]);
    
    name = "plots/loEA/detPtResponses" + etaBinName[e] + ".pdf";
    can3->SaveAs(name,"PDF");

  }



  for (int e=0; e<nEtaBins; ++e) {
    for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
      for (int i=5; i<26; ++i) {
	int ptVal = i + 5;
	int binno = i + 1;
	if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	  double wt = 1./hMissedJets->GetBinContent(binno);  // weight according to misses at part-level
	  hDetWt[p][e]->Add(hDet[i][e],wt);
	}
      }
      hDetWt[p][e]->Scale(1./hDetWt[p][e]->Integral());
      hDetWt[p][e]->GetYaxis()->SetRangeUser(0.0,0.225);
      hDetWt[p][e]->SetStats(0);
      if (p==0) { hDetWt[p][e]->Draw(); }
      else { hDetWt[p][e]->Draw("SAME"); }
    }
    name = "plots/loEA/WeightedPtResponse"; name += etaBinName[e] + ".pdf";
    can3->BuildLegend(0.65,0.7,0.9,0.9);
    tex->DrawLatex( 0.4, 0.95, etaBinString[e] );
    can3->SaveAs(name,"PDF");
  }


  for (int p=0; p<nPtBins; ++p) {
    for (int e=0; e<nEtaBins; ++e) {
      if (e==0) { hDetWt[p][e]->Draw(); }
      else { hDetWt[p][e]->Draw("SAME"); }
    }
    name = "plots/loEA/WeightedPtResponse"; name += ptBinName[p] + ".pdf";
    can3->BuildLegend(0.65,0.7,0.9,0.9);
    tex->DrawLatex( 0.4, 0.92, ptBinString[p] );
    can3->SaveAs(name,"PDF");
  }

  
  // OPEN pAu FILE AND GATHER HISTOGRAMS OF UE AND LEAD PT
  TFile *UEfile = new TFile("../out/UE/pAuHTjetUE_loEA.root","READ");

  const int nEffEtaBins = 20;
  
  TH3D *hChgUE3D[nEtaBins];
  TH2D *hChgUE2D_uncorr[nEtaBins][nEffEtaBins];
  TH2D *hChgUE2D_te[nEtaBins][nEffEtaBins];

  for (int e=0; e<nEtaBins; ++e) {
    name = "hChgUE" + etaBinName[e] + "Jet";
    hChgUE3D[e] = (TH3D*)UEfile->Get(name);    
    for (int f=0; f<nEffEtaBins; ++f) {
      double binLo = (double)((0.1*f)-1.0);
      double binHi = (double)((0.1*f)-0.9);
      hChgUE3D[e]->GetZaxis()->SetRangeUser(binLo,binHi);
      hChgUE2D_uncorr[e][f] = (TH2D*)hChgUE3D[e]->Project3D("YX");
      name = "hChgUE2D_te" + etaBinName[e] + "_" + f;
      title = ";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)";
      hChgUE2D_te[e][f] = new TH2D( name, title, 55,4.5,59.5, 15,0.0,15.0 );
    }
    hChgUE3D[e]->GetZaxis()->SetRangeUser(-1.0,1.0);
  }

  TH1D *hLeadJetPt[nEtaBins];
  TH2D *hChgUE2D[nEtaBins];

  TCanvas *can4[nEtaBins];
  
  for (int e=0; e<nEtaBins; ++e) {

    name = "can4_" + e;
    can4[e] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 4
    can4[e]->cd();   can4[e]->SetLogz();   can4[e]->SetLogy(0);
  
    name = "hLeadPt" + etaBinName[e] + "Jet";
    hLeadJetPt[e] = (TH1D*)UEfile->Get(name);

    hChgUE2D[e] = (TH2D*)hChgUE3D[e]->Project3D("YX"); // SAVE 2D HISTO OF UE PT vs. LEADING JET PT
    hChgUE2D[e]->GetXaxis()->SetRangeUser(4.5,59.5);
    hChgUE2D[e]->GetYaxis()->SetRangeUser(0.0,15.0);
    hChgUE2D[e]->GetYaxis()->SetTitleOffset(1.25);
    can4[e]->SetLogy(0);
    hChgUE2D[e]->Draw("COLZ");
    name = "plots/loEA/hChgUE2D" + etaBinName[e] + "Jet.pdf";
    can4[e]->SaveAs(name,"PDF");
    hChgUE2D[e]->GetXaxis()->SetRangeUser(4.5,59.5);
    hChgUE2D[e]->GetYaxis()->SetRangeUser(0.0,15.0);

  }

  TH1D *hUEpt[55][nEtaBins];  TH1D *hWtUEpt[nPtBins][nEtaBins];
  TCanvas *can5[nEtaBins];

  for (int e=0; e<nEtaBins; ++e) {

    name = "can5_" + e;
    can5[e] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 4
    can5[e]->cd();   can5[e]->SetLogz();   can5[e]->SetLogy(0);

    for (int p=0; p<nPtBins; ++p) {
      name = "hWtUEpt" + etaBinName[e] + "Jet" + ptBinName[p];
      title = ptBinString[p] + " ~ " + etaBinString[e] + ";ch UE p_{T} (GeV)";
      hWtUEpt[p][e] = new TH1D(name,title,15,0.0,15.0);
    }
  
    can5[e]->cd();
    can5[e]->SetLogy();
    for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
      int ptVal = i + 5;
      int binno = i + 1;
      name = emw[e] + " " + ptVal + " GeV det. jet ";
      hUEpt[i][e] = (TH1D*)hChgUE2D[e]->ProjectionY(name,binno,binno);
      hUEpt[i][e]->GetXaxis()->SetRangeUser(4.5,59.5);
      if ( (int)hLeadJetPt[e]->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }
      int n_Jets = hLeadJetPt[e]->GetBinContent( binno );
      if (n_Jets>0){
	hUEpt[i][e]->Scale(1./n_Jets);
	// name = "plots/loEA/UEpt" + etaBinName[e] + "_" + ptVal + "GeV" + ".pdf";
	// hUEpt[i][e]->Draw();          can5[e]->SaveAs(name,"PDF");

	for (int p=0; p<nPtBins; ++p) {  // ADD UE DISTRIBUTIONS WEIGHTED BY DET-LEVEL FRACTIONAL CONTRIBUTION
	  if (hDetWt[p][e]->GetBinContent(binno)==0) { continue; }
	  double wt = hDetWt[p][e]->GetBinContent(binno)*( 1.0 - hFakeJets[e]->GetBinContent(binno) );
	  hWtUEpt[p][e]->Add( hUEpt[i][e], wt );
	}
	hUEpt[i][e]->SetMarkerStyle(marker[i]);
	hUEpt[i][e]->SetMarkerSize(1);
	hUEpt[i][e]->SetStats(0);
	// TString hName = ""; hName += ptVal + " GeV " + emw[e] + " det. jet ";
	// //title = hName + ";chg. UE particle p_{T} (GeV)";
	hUEpt[i][e]->SetNameTitle(name,";ch UE p_{T} (GeV)");
      }
    }

    hUEpt[0][e]->GetYaxis()->SetRangeUser(0.0000005,2);
    hUEpt[0][e]->Draw("PLC PMC");

    for (int i=1; i<55; ++i) { if (hUEpt[i][e]->GetEntries()!=0) { hUEpt[i][e]->Draw("SAME PLC PMC"); } }
    can5[e]->BuildLegend(0.65,0.1,0.9,0.9);
    name = "plots/loEA/detUEpt" + etaBinName[e] + "Jet.pdf";
    tex->DrawLatex( 0.4, 0.95, etaBinString[e] );
    can5[e]->SaveAs(name,"PDF");

  }

  //hChgUE2D_uncorr[e][f]
  TH1D *hEffic[nEffEtaBins];
  TFile *ef = new TFile( "src/trackeffic.root", "READ" );

  for (int f=1; f<=nEffEtaBins; ++f) {
    
    TString histoName = "eff_s_bin_1_10_bbc__";
    histoName += f; histoName +=  "_"; histoName += f; histoName += "_eta";
    
    hEffic[f] = (TH1D*)ef->Get( histoName );
  
    double pt, effic, ptVal, corrVal, corrErr, jetPt; // eta, phi, e, px, py, pz, 
    int ptBin;

    for (int e=0; e<nEtaBins; ++e) { // ~ ~ ~ ~ ~ ~ ~ ~ BEGIN ETA LOOP ~ ~ ~ ~ ~ ~ ~ ~
      
      for ( int iy=1; iy<=hChgUE2D_uncorr[e][f-1]->GetNbinsY(); ++iy ) {  // loop over chg UE pT bins (y)
    	pt = hChgUE2D_uncorr[e][f-1]->GetYaxis()->GetBinCenter(iy);
    	if ( pt > 3.0 ) { pt = 3.0; }
    	ptBin = hEffic[f]->FindBin( pt );    // find histogram bin corresponding to track pt
    	effic = hEffic[f]->GetBinContent( ptBin );
	
    	for ( int ix=1; ix<=hChgUE2D_uncorr[e][f-1]->GetNbinsX(); ++ix ) { // loop over all lead jet pt bins (x)

	  ptVal = hChgUE2D_uncorr[e][f-1]->GetYaxis()->GetBinCenter(iy);
    	  corrVal = hChgUE2D_uncorr[e][f-1]->GetBinContent(ix,iy)/effic;      // calculate corrected bin content and error
    	  corrErr = hChgUE2D_uncorr[e][f-1]->GetBinError(ix,iy)/effic;       // (divide bin content and error by efficiency)
	  jetPt = hChgUE2D_uncorr[e][f-1]->GetXaxis()->GetBinCenter(ix);
	  ptBin = hChgUE2D_uncorr[e][f-1]->GetYaxis()->FindBin( ptVal );
	  
    	  // hChgUE2D_te[e][f-1]->Fill( jetPt, ptVal, corrVal );
    	  hChgUE2D_te[e][f-1]->SetBinContent( ix, ptBin, corrVal );
    	  hChgUE2D_te[e][f-1]->SetBinError( ix, ptBin, corrErr );

    	}
      }
	
    }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ END ETA LOOP ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 

    
  }
  

  TH2D *hChgUE2D_corr[nEffEtaBins];

  for (int e=0; e<nEtaBins; ++e) {
    name = "hChgUE2D_corr" + etaBinName[e];
    title = ";leading jet p_{T} (GeV);corrected chg. UE part. p_{T} (GeV)";
    hChgUE2D_corr[e] = new TH2D( name, title, 55,4.5,59.5, 15,0.0,15.0 );
  }

  for (int e=0; e<nEtaBins; ++e) {
    int nHistos = 0;
    for (int f=0; f<nEffEtaBins; ++f) {
      
      double eLo = -1.00 + (f*2.00/nEffEtaBins);
      double eHi = eLo + 0.10;
      
      int eval = 99;
      for (int e=0; e<nEtaBins; ++e) {    if (eLo>=etaLo[e]) {eval = e;}    }
      if (eval==99) { cerr<<"error finding eta bin of histogram!   "<<eLo<<endl; continue; }

      double wt2D = hChgUE2D_te[e][f]->Integral()/hChgUE2D[e]->Integral();

      hChgUE2D_corr[e]->Add( hChgUE2D_te[e][f], wt2D );
      nHistos +=1;
      // cout<<etaBinString[e]<<"\t"<<eLo<<"\t"<<eHi<<endl;
    }

    hChgUE2D_corr[e]->Scale(1./nHistos);
    // cout<<nHistos<<endl<<endl;
  }
  
  TCanvas * can6 = new TCanvas( "can6" , "" ,700 ,500 );              // CANVAS 6
  can6->SetLogz();
  for (int e=0; e<nEtaBins; ++e) {
    hChgUE2D_corr[e]->Draw("COLZ");
    name = "plots/loEA/ChgUE2D_corr" + etaBinName[e] + "Jet.pdf";
    can6->SaveAs(name,"PDF");
  }


  TCanvas * can7[nEtaBins];            // CANVAS 7
  TH1D *hUEpt_corr[55][nEtaBins];  TH1D *hWtUEpt_corr[nPtBins][nEtaBins];

  int nHistos[nPtBins][nEtaBins] = { 0,0,0,0,0,0,0,0,0 };

  double meanSum = 0;
  
  for (int e=0; e<nEtaBins; ++e) {

    name = "can7_" + e;
    can7[e] = new TCanvas( name , "" ,700 ,500 );              // CANVAS 4
    can7[e]->cd();   can7[e]->SetLogz();   can7[e]->SetLogy(0);

    for (int p=0; p<nPtBins; ++p) {
      name = "hWtUEpt_corr" + etaBinName[e] + "Jet" + ptBinName[p];
      title = ptBinString[p] + " ~ " + etaBinString[e] + ";ch UE p_{T} (GeV)";
      hWtUEpt_corr[p][e] = new TH1D(name,title,15,0.0,15.0);
    }
  
    can7[e]->cd();
    can7[e]->SetLogy();
    for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
      int ptVal = i + 5;
      int binno = i + 1;
      name = emw[e] + " " + ptVal + " GeV det. jet ";
      hUEpt_corr[i][e] = (TH1D*)hChgUE2D_corr[e]->ProjectionY(name,binno,binno);
      if ( (int)hLeadJetPt[e]->GetBinCenter(binno) != ptVal ) { cerr<<"failed pT matching"<<endl; }
      int nJets = hLeadJetPt[e]->GetBinContent( binno );
      if (nJets>0){
	hUEpt_corr[i][e]->Scale(1./nJets);
	// name = "plots/loEA/UEpt" + etaBinName[e] + "_" + ptVal + "GeV" + ".pdf";
	// hUEpt_corr[i][e]->Draw();          can7[e]->SaveAs(name,"PDF");

	for (int p=0; p<nPtBins; ++p) {  // ADD UE DISTRIBUTIONS WEIGHTED BY DET-LEVEL FRACTIONAL CONTRIBUTION
	  if (hDetWt[p][e]->GetBinContent(binno)==0) { continue; }
	  double wt = (hDetWt[p][e]->GetBinContent(binno))*( 1.0 - hFakeJets[e]->GetBinContent(binno) );
	  hWtUEpt_corr[p][e]->Add( hUEpt_corr[i][e], wt );
	 
	  nHistos[p][e] += 1;
	}
	// cout<<nHistos<<endl;
	  
	hUEpt_corr[i][e]->SetMarkerStyle(marker[i]);
	hUEpt_corr[i][e]->SetMarkerSize(1);
	hUEpt_corr[i][e]->SetStats(0);
	hUEpt_corr[i][e]->SetNameTitle(name,";ch UE p_{T} (GeV)");
      }
    }

    hUEpt_corr[0][e]->GetYaxis()->SetRangeUser(0.000005,10);
    hUEpt_corr[0][e]->Draw("PLC PMC");

    for (int i=1; i<55; ++i) { if (hUEpt_corr[i][e]->GetEntries()!=0) { hUEpt_corr[i][e]->Draw("SAME PLC PMC"); } }
    can7[e]->BuildLegend(0.65,0.1,0.9,0.9);
    name = "plots/loEA/detUEpt_corr" + etaBinName[e] + "Jet.pdf";
    tex->DrawLatex( 0.4, 0.95, etaBinString[e] );
    can7[e]->SaveAs(name,"PDF");

  }

  for (int p=0; p<nPtBins; ++p) {
    for (int e=0; e<nEtaBins; ++e) {
      //hWtUEpt_corr[p][e]->Scale(1./nHistos[p][e]);
      cout<<hWtUEpt_corr[p][e]->GetMean(1)<<"   "<<hWtUEpt_corr[p][e]->Integral()/AREA<<endl;
    }
  }
  
  TCanvas * can8 = new TCanvas( "can8" , "" ,700 ,500 );              // CANVAS 8
  can8->SetLogy();
  for (int e=0; e<nEtaBins; ++e) {
    for (int p=0; p<nPtBins; ++p) {

      hWtUEpt_corr[p][e]->Draw();

      TString text = "#LT p_{T}^{ch}#GT = "; text += hWtUEpt_corr[p][e]->GetMean(1);
      text = text(0,26);
      
      TLatex *tex = new TLatex( 0.6, 0.55,text);
      tex->SetTextFont(63);
      tex->SetTextSize(20);
      tex->SetTextColor(kBlack);
      tex->SetLineWidth(1);
      tex->SetNDC();
      tex->Draw();
      
      meanSum += hWtUEpt_corr[p][e]->Integral()/AREA;
      text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += hWtUEpt_corr[p][e]->Integral()/AREA;
      text = text(0,42);
      tex->DrawLatex( 0.6, 0.4, text );
      
      TString name = "plots/loEA/wtUEpt_corr" + ptBinName[p] + etaBinName[e] + ".pdf";
      can8->SaveAs(name,"PDF");
    }
  }

  cout<<meanSum/9.0<<endl;

  hMissedJets->Draw();
  can8->SaveAs("plots/loEA/missedJets.pdf","PDF");

  
  return 0;
}
