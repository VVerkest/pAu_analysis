// Veronica Verkest
// September 10, 2020

#include "params.hh"
#include "UEfuncs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString directory = "plots/ptDiff/";
  TString EAstring = "hiEA";
  TString name, saveName, title, avg, sigma, drawString;
  int jeval, ueeval, pval, eaval;
  const int nEffEtaBins = 10;
  
  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", 55,4.5,59.5, 30,0.0,30.0, 20,-1.0,1.0);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.5,59.5);
  TH1D *hLeadJetPt[nEtaBins];
  TH3D *hUE3D[nEtaBins];
  TH1D *hUE[55][nEffEtaBins];
  TH1D *hEffic[nEffEtaBins];
  TH1D *hUE_te[55][nEffEtaBins];
  TH1D *hAddedUE_te[nPtBins][nEtaBins];

  for (int e=0; e<nEtaBins; ++e) {
    for (int p=0; p<nPtBins; ++p) {
      name = "Added_" + emw[e] + "UE_te" + ptBinName[p];
      hAddedUE_te[p][e] = new TH1D(name,";te corrected chg. UE part. p_{T} (GeV)",15,0.0,15.0);
    }
  }
  
  TString inFileName = "../out/UE/pAuHTjetUE_" + EAstring + ".root";

  // eventually add loop to go over loEA and hiEA (use canvas destructors)
  
  ////////////////////////////////////////// pAu DATA FILE FOR UE HISTOGRAMS //////////////////////////////////////////
  TFile *UEfile = new TFile(inFileName,"READ");  // OPEN pAu FILE AND GATHER HISTOGRAMS OF UE AND LEAD PT
  
  for (int e=0; e<nEtaBins; ++e) {
    name = "hLeadPt" + etaBinName[e] + "Jet";  // X=leadPt, Y=UEpt, Z=UEeta
    hLeadJetPt[e] = (TH1D*)UEfile->Get(name);
    hleadPt->Add(hLeadJetPt[e]);
    name = "hChgUE" + etaBinName[e] + "Jet";
    hUE3D[e] = (TH3D*)UEfile->Get(name);    
    hChgUE3D->Add(hUE3D[e]);
  }

  
  for ( int e=0; e<nEffEtaBins; ++e ) { // Z=UEeta

    TCanvas *can0 = new TCanvas( "can0" , "" ,700 ,500 );              // CANVAS 0
    can0->SetLogy();

    hChgUE3D->GetZaxis()->SetRange(e+1,e+1);
    hChgUE3D->GetYaxis()->SetRangeUser(0.0,15.0);
    for ( int jp=0; jp<55; ++jp ) {
      int ptVal = jp + 5;
      int binno = jp + 6;
      
      name = EAstring + "_eta"; name += e; name += "UE_"; name += ptVal; name += "GeVdetJet";
      hUE[jp][e] = (TH1D*)hChgUE3D->ProjectionY(name,binno,binno,e+1,e+1);
      hUE[jp][e]->Draw("SAME PLC PMC");
      saveName = directory + EAstring + "_eta"; saveName += e; saveName += "UE.pdf";

      name = EAstring + "_eta"; name += e; name += "UE_TRACKEFFIC_"; name += ptVal; name += "GeVdetJet";      
      hUE_te[jp][e] = new TH1D(name,";te corrected chg. UE part. p_{T} (GeV)",15,0.0,15.0);
    }
    // can0->SaveAs(saveName);
    can0->Destructor();
  }

  ////////////////////////////////// TRACKING EFFICIENCY FILE FOR UE pT CORRECTION //////////////////////////////////
  TFile *ef = new TFile( "src/trackeffic.root", "READ" );

  
  for (int e=0; e<nEffEtaBins; ++e) {  // loop over 10 eta bins

    name = "eff_s_bin_" + GetEfficHistoName(EAstring) + "_bbc__"; name += e+1; name += "_"; name += e+1; name += "_eta";
    hEffic[e] = (TH1D*)ef->Get( name );

    for ( int jp=0; jp<55; ++jp ) {

      for(int ix = 0; ix <15; ++ix){
	int xbin = ix+1;
      
	double pt = hUE[jp][e]->GetBinCenter(xbin);
	if ( pt > 3.0 ) { pt = 3.0; }
	double ptbin = hEffic[e]->FindBin(pt);
	double eff = hEffic[e]->GetBinContent(ptbin);
	double old_value = hUE[jp][e]->GetBinContent(xbin);
	double old_err = hUE[jp][e]->GetBinError(xbin);
	double corr_value = (double) old_value/eff;
	double corr_err = (double) old_err/eff;
	hUE_te[jp][e]->SetBinContent( xbin, corr_value);
	hUE_te[jp][e]->SetBinError( xbin, corr_err );
      }
      hUE_te[jp][e]->SetEntries(hUE[jp][e]->GetEntries());    
    }
  }


  for (int e=0; e<=nEffEtaBins; ++e) {  // loop over 10 eta bins
    double etaValue = -10 + (2*e);
    etaValue /= 10;
    // cout<<etaValue<<endl;
    int eval = 99;      int pval = 99;
    for (int ie=0; ie<nEtaBins; ++ie) {    if (etaValue>=etaLo[ie] && etaValue<=etaHi[ie]) {eval = ie;}    }
    // cout<<eval<<endl<<endl;
    for ( int jp=0; jp<55; ++jp ) {
      double ptValue = jp + 5.0;
      for (int p=0; p<nPtBins; ++p) {    if ( ptValue>=ptLo[p] && ptValue<=ptHi[p] ) {pval = p;}    }
      if (eval==99 ) { std::cerr<<"error finding eta bin of histogram!   "<<etaValue<<std::endl; continue; }
      if (pval==99 ) { continue; }

      // cout<<eval<<"  "<<pval<<endl;
      
      hAddedUE_te[pval][eval]->Add(hUE_te[jp][e]);
    }
  }

  for (int e=0; e<nEtaBins; ++e) {
    TCanvas *can1 = new TCanvas( "can1" , "" ,700 ,500 );              // CANVAS 0
    can1->SetLogy();
    
    for (int p=0; p<nPtBins; ++p) {
      hAddedUE_te[p][e]->SetLineColor(etaColor[e]);
      hAddedUE_te[p][e]->SetMarkerColor(etaColor[e]);
      hAddedUE_te[p][e]->SetMarkerStyle(ptMarker[p]);
      int binlo = hleadPt->FindBin( ptLo[p] );
      int binhi = hleadPt->FindBin( ptHi[p] );
      hAddedUE_te[p][e]->SetStats(0);
      hAddedUE_te[p][e]->Scale(1./hleadPt->Integral( binlo, binhi ));
      cout<<hAddedUE_te[p][e]->GetName()<<"  \t"<<binlo<<"  "<<binhi<<endl;
      hAddedUE_te[p][e]->GetYaxis()->SetRangeUser(0.000001,100);
      hAddedUE_te[p][e]->Draw("SAME");
    }
    name = directory + "AddedUE_te" + etaBinName[e] + ".pdf";
    can1->BuildLegend();
    can1->SaveAs(name,"PDF");
    can1->Destructor();
  }

  TFile *outFile = new TFile("outFile.root","RECREATE");
  for (int e=0; e<nEtaBins; ++e) {
    for (int p=0; p<nPtBins; ++p) {
      hAddedUE_te[p][e]->Write();
    }
  }
  
  return 0;
}
