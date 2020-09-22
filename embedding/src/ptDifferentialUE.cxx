// Veronica Verkest
// September 10, 2020

#include "params.hh"
#include "UEfuncs.hh"

using namespace std;
using namespace Analysis;

int main () {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TString directory = "plots/ptDiff/";
  TString EAstring[2] = { "loEA", "hiEA" };
  TString name, saveName, title, avg, sigma, drawString;
  int jeval, ueeval, pval, eaval;
  const int nEffEtaBins = 10;
  const int ybins = 10;
  double ybinEdge[ybins+1] = { -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1 };
  int nJets[nPtBins];
  int etaStartBin[nEtaBins] = {1,8,14};
  int EAptMarker[2][nPtBins] = {{107,108,110},{20,21,33}};
  
  TH3D *hChgUE3D = new TH3D("hChgUE3D",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta", 55,4.0,59.0, 30,0.0,30.0, 20,-1.0,1.0);
  TH1D *hleadPt = new TH1D("hleadPt",";leading jet p_{T} (GeV)", 55,4.0,59.0);
  TH1D *hLeadJetPt[nEtaBins];
  TH3D *hUE3D[nEtaBins];
  TH2D *hUE2D[55][nEtaBins];
  TH1D *hEffic[nEffEtaBins];
  TH2D *hUE2D_te[55][nEtaBins];
  TH2D *hAdded2DUE_te[nPtBins][nEtaBins];
  TH1D *hAddedUE_te[nPtBins][nEtaBins];
  TFile *outFile;
  
  for (int a=0; a<2; ++a) {
  
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	name = EAstring[a] + "_Added_" + emw[e] + "UE_te" + ptBinName[p];
	hAdded2DUE_te[p][e] = new TH2D(name,";te corrected chg. UE part. p_{T} (GeV);chg. UE part. #eta",15,0.0,15.0,20,-1.0,1.0);
      }
    }
  
    TString inFileName = "../out/UE/pAuHTjetUE_30cmVzCut_3cmVzDiff_" + EAstring[a] + "__binShift.root";
  
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


    hChgUE3D->GetYaxis()->SetRangeUser(0.0,15.0);  
    for (int e=0; e<nEtaBins; ++e) { // Z=UEeta

      TCanvas *can0 = new TCanvas( "can0" , "" ,700 ,500 );              // CANVAS 0
      can0->SetLogz();

      hChgUE3D->GetZaxis()->SetRangeUser(etaLo[e],etaHi[e]);

      for ( int jp=0; jp<55; ++jp ) {
	int binno = jp + 1;
	// int ml = jp + 5;
	int ptVal = jp + 4;  // gnome
	
	hChgUE3D->GetXaxis()->SetRange(binno,binno);
      
	name = EAstring[a]; name += "_UE2D"; name += etaBinName[e] + "_"; name += ptVal; name += "GeVdetJet";
	hUE2D[jp][e] = (TH2D*)hChgUE3D->Project3D("ZY");
	hUE2D[jp][e]->SetName(name);
	if (hleadPt->Integral(binno,binno)==0) { hUE2D[jp][e]->Scale(0); }
	else { hUE2D[jp][e]->Scale(1./hleadPt->Integral(binno,binno)); }
	if (hUE2D[jp][e]->GetEntries()>0) {
	  hUE2D[jp][e]->GetYaxis()->SetRangeUser(etaLo[0],etaHi[2]);
	  hUE2D[jp][e]->SetAxisRange(etaLo[0],etaHi[2],"Y");
	  hUE2D[jp][e]->Draw("COLZ");
	  name = directory + EAstring[a]; name += "_UE2D"; name += etaBinName[e] + "_"; name += ptVal; name += "GeVdetJet.pdf";
	  can0->SaveAs(name,"PDF");
	  // cout<<hUE2D[jp][e]->ProjectionX()->Integral()/area[e]<<endl;
	}
	name = EAstring[a]; name += "_UE2D_TRACKEFFIC"; name += etaBinName[e] + "_"; name += ptVal; name += "GeVdetJet";      
	hUE2D_te[jp][e] = new TH2D(name,";te corrected chg. UE part. p_{T} (GeV);chg. UE part. #eta",15,0.0,15.0,20,-1.0,1.0);
	hUE2D_te[jp][e]->SetName(name);
      }
      // can0->SaveAs(saveName);
      can0->Destructor();
    }

    ////////////////////////////////// TRACKING EFFICIENCY FILE FOR UE pT CORRECTION //////////////////////////////////
    TFile *ef = new TFile( "src/trackeffic.root", "READ" );

    for (int e=0; e<nEffEtaBins; ++e) {  // loop over 10 eta bins
      name = "eff_s_bin_" + GetEfficHistoName(EAstring[a]) + "_bbc__"; name += e+1; name += "_"; name += e+1; name += "_eta";
      hEffic[e] = (TH1D*)ef->Get( name );
    }

    for (int e=0; e<nEtaBins; ++e) {
      for ( int jp=0; jp<55; ++jp ) {

	for(int ix = 0; ix <hUE2D[jp][e]->GetNbinsX(); ++ix){
	  int xbin = ix+1;
      
	  double pt = hUE2D[jp][e]->GetXaxis()->GetBinCenter(xbin);
	  if ( pt > 3.0 ) { pt = 3.0; }

	  for (int iy=0; iy<hUE2D[jp][e]->GetNbinsY(); ++iy) {  // loop over 20 eta bins

	    int fileNo = 99;
	    int ybin = iy+1;
	  
	    double binCenter = hUE2D[jp][e]->GetYaxis()->GetBinCenter( ybin );
	    for (int j=0; j<ybins; ++j) {
	      if ( binCenter>ybinEdge[j] && binCenter<ybinEdge[j+1] ) { fileNo = j; };
	    }
	    if ( fileNo==99 ) { cerr<<"CANNOT FIND EFFICIENCY HISTOGRAM"<<endl; return -1; }

	    double ptbin = hEffic[fileNo]->FindBin(pt);
	    double eff = hEffic[fileNo]->GetBinContent(ptbin);

	    double old_value = hUE2D[jp][e]->GetBinContent(xbin,ybin);
	    double old_err = hUE2D[jp][e]->GetBinError(xbin,ybin);
	    double corr_value = (double) old_value/eff;
	    double corr_err = (double) old_err/eff;
	    hUE2D_te[jp][e]->SetBinContent( xbin, iy+etaStartBin[e], corr_value);
	    hUE2D_te[jp][e]->SetBinError( xbin,iy+etaStartBin[e],corr_err );
	  }
	}
	hUE2D_te[jp][e]->SetEntries( hUE2D[jp][e]->GetEntries() );
      }
    }

    int nHistos[nPtBins][nEtaBins] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int e=0; e<nEtaBins; ++e) {
      for ( int jp=0; jp<55; ++jp ) {
	
	// double ptValue = jp + 5.0;
	double ptValue = jp + 4.0; // gnome
	
	if ( ptValue<ptLo[0] || ptValue>ptHi[2] ) {continue;}
	if (hUE2D_te[jp][e]->GetEntries()<1) { continue; }
      
	for (int p=0; p<nPtBins; ++p) {
	  if ( ptValue>=ptLo[p] && ptValue<=ptHi[p] ) {
	    pval = p;    
	    
	    // cout<<ptValue<<"  "<<pval<<endl;
	    
	    hUE2D_te[pval][e]->GetXaxis()->SetRangeUser(0.0,15.0);
	    hUE2D_te[pval][e]->GetYaxis()->SetRangeUser(-1.0,1.0);
	    hAdded2DUE_te[pval][e]->GetXaxis()->SetRangeUser(0.0,15.0);
	    hAdded2DUE_te[pval][e]->GetYaxis()->SetRangeUser(-1.0,1.0);
      
	    if ( (hUE2D_te[pval][e]->GetNbinsY()!=hAdded2DUE_te[pval][e]->GetNbinsY()) || (hUE2D_te[pval][e]->GetNbinsX()!=hAdded2DUE_te[pval][e]->GetNbinsX()) ) { cout<<hUE2D_te[pval][e]->GetName()<<endl<<ptValue<<endl; break; }

	    // cout<<hUE2D_te[jp][e]->GetEntries()<<endl;
      
	    hAdded2DUE_te[pval][e]->Add(hUE2D_te[jp][e]);
	    nHistos[pval][e] += 1;
	  }
	}
	if (pval==99 ) { continue; }
      }
    }
    

    for (int e=0; e<nEtaBins; ++e) {
      TCanvas *can1 = new TCanvas( "can1" , "" ,700 ,500 );              // CANVAS 0
      can1->SetLogy();
    
      for (int p=0; p<nPtBins; ++p) {
	name = "hUE_te" + etaBinName[e] + ptBinName[p];
	// hAdded2DUE_te[p][e]->Scale(1./nHistos[p][e]);
	// cout<<nHistos[p][e]<<endl;
	int binlo = hleadPt->GetBin(ptLo[p]);
	int binhi = hleadPt->GetBin(ptHi[p]);	
	// hAdded2DUE_te[p][e]->Scale(1./hleadPt->Integral(binlo,binhi));
	cout<<EAstring[a]<<"  "<<ptBinString[p]<<"  "<<hleadPt->Integral(binlo,binhi)<<endl;
	hAddedUE_te[p][e] = (TH1D*)hAdded2DUE_te[p][e]->ProjectionX(name);
	hAddedUE_te[p][e]->SetLineColor(etaColor[e]);
	hAddedUE_te[p][e]->SetMarkerColor(etaColor[e]);
	hAddedUE_te[p][e]->SetMarkerStyle(EAptMarker[a][p]);
	hAddedUE_te[p][e]->SetStats(0);
	// cout<<hAddedUE_te[p][e]->GetName()<<"  \t"<<binlo<<"  "<<binhi<<endl;
	hAddedUE_te[p][e]->GetYaxis()->SetRangeUser(0.000001,100);
	hAddedUE_te[p][e]->Draw("SAME");
      }
      name = directory + "AddedUE_te" + etaBinName[e] + ".pdf";
      can1->BuildLegend();
      can1->SaveAs(name,"PDF");
      can1->Destructor();
    }



  
    TCanvas *can1 = new TCanvas( "can1" , "" ,700 ,500 );              // CANVAS 0
    can1->SetLogz();

    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	hAdded2DUE_te[p][e]->Draw("COLZ");
	name = directory + hAdded2DUE_te[p][e]->GetName() + ".pdf";
	can1->SaveAs(name,"PDF");
      }
    }

    can1->Destructor();

  
    if (a==0) { outFile = new TFile("outFile.root","RECREATE"); }
    else { outFile = new TFile("outFile.root","UPDATE"); }
    for (int e=0; e<nEtaBins; ++e) {
      for (int p=0; p<nPtBins; ++p) {
	hAdded2DUE_te[p][e]->Write();
	name = EAstring[a] + "_" + hAddedUE_te[p][e]->GetName();
	hAddedUE_te[p][e]->SetName(name);
	hAddedUE_te[p][e]->Write();
      }
    }
    
    outFile->Close();
  }
  
  return 0;
}
