// Veronica Verkest
// February 24, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;


TH1D* WeightAndSumByFC( TH1D* FC, TH1D* UE1D[55] ){

  // const int bins = 15;
  // double binEdge[bins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };
  const int bins = 14;
  const double binEdge[bins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0 };  TH1D* UE_part = new TH1D("UE_part","",bins,binEdge);

  for (int i=0; i<FC->GetNbinsX(); ++i) {
    int binno = i+4;
    double wt = FC->GetBinContent(binno);
    if (wt==0) { continue; }
    UE_part->Add(UE1D[i], wt);
  }

  UE_part->Write();
  return UE_part;
}




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
  TString dirName = "plots/test";

  TFile *inFile[nEAbins], *embFile[nEAbins];
  TFile* outFile = new TFile("out/test0.root", "RECREATE");

  TH2D *hUE2D[nEAbins];
  TH3D *hUE3D[nEAbins][nEtaBins];
  TH3D *hUE3Dsum[nEAbins];
  TH1D *hLead[nEAbins][nEtaBins];
  TH1D *hLeadSum[nEAbins];
  TH2D *hResponse[nEAbins][nEtaBins];
  TH2D *hResponseSum[nEAbins];
  
  for (int a=0; a<nEAbins; ++a) {
    //dirName += lohi[a];

    name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
    inFile[a] = new TFile(name, "READ");

    name = "hUE3Dsum_" + lohi[a] + "EA";
    hUE3Dsum[a] = new TH3D(name,";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta",xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);

    name = "hLead_" + lohi[a];
    hLeadSum[a] = new TH1D(name,";leading jet p_{T} (GeV)", 55,4.0,59.0);

    name = "hResponseSum_" + lohi[a] + "EA";
    hResponseSum[a] = new TH2D(name,"part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);

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

      hUE2D[a] = (TH2D*)hUE3Dsum[a]->Project3D("YX");
      name = "hUE2Dsum_" + lohi[a];
      hUE2D[a]->SetName(name);
    }

    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    embFile[a] = new TFile(name, "READ");
    
    for (int e=0; e<nEtaBins; ++e) {
      name = "hPtResponse_" + emw[e] + "EtaJet";
      hResponse[a][e] = (TH2D*)embFile[a]->Get(name);
      name = "hResponse_" + lohi[a] + "EA_" + emw[e] + "Jet";
      hResponse[a][e]->SetName(name);
      hResponseSum[a]->Add(hResponse[a][e]);
    }

  }




  
  TH1D *hUE1D[nEAbins][55], *hUE1D_dbw[nEAbins][55], *hUE1D_noJES_dbw[nEAbins][nPtBins], *hUE1D_noJES[nEAbins][nPtBins];

  for (int a=0; a<nEAbins; ++a) {

    for (int p=0; p<nPtBins; ++p) {
      name = "hUE1D_noJES" + lohi[a] + "EA"; name +=ptBinName[p];
      hUE1D_noJES[a][p] = new TH1D(name,"",ybins,ybinEdge);
      name += "_dbw";
      hUE1D_noJES_dbw[a][p] = new TH1D(name,"",ybins,ybinEdge);
      
      ProjectScaleAndSaveUE1D( hUE2D[a], *hUE1D_noJES[a][p], hLeadSum[a], *hUE1D_noJES_dbw[a][p], ptLo[p], ptHi[p], dirName, lohi[a] );
    }
    
    for (int i=0; i<55; ++i) {
        name = "hUE1D_" + lohi[a] + "EA_"; name +=(i+4); name+="_"; name+=(i+5);
        hUE1D[a][i] = new TH1D(name,"",ybins,ybinEdge);
        name += "_dbw";
        hUE1D_dbw[a][i] = new TH1D(name,"",ybins,ybinEdge);

        ProjectScaleAndSaveUE1D( hUE2D[a], *hUE1D[a][i], hLeadSum[a], *hUE1D_dbw[a][i], i+4., i+5., dirName, lohi[a] );
    }
    
  }


  
  outFile->cd();

  TH1D *hResponseProj[nEAbins][nPtBins]; TH1D* hUE_partProj[nEAbins][nPtBins];

  for (int a=0; a<nEAbins; ++a) {
    for (int p=0; p<nPtBins; ++p) {

        name = "hResponse" + lohi[a] + "EA_" + ptBinName[p];
        hResponseProj[a][p] = GenerateFractionalContribution( hResponseSum[a], ptLo[p], ptHi[p], dirName, lohi[a] );
        hResponseProj[a][p]->SetName(name);

        name = "hUE_part" + lohi[a] + "EA_" + ptBinName[p];
        hUE_partProj[a][p] = (TH1D*) WeightAndSumByFC( hResponseProj[a][p], hUE1D_dbw[a] );
        hUE_partProj[a][p]->SetName(name);
        hUE_partProj[a][p]->SetAxisRange(0.0000001,1,"Y");
        hUE_partProj[a][p]->Draw();
    
        TString text = "#LT #frac{dN^{ch}}{d#eta d#phi} #GT = "; text+=RoundDecimal(hUE_partProj[a][p]->Integral()/AREA,3);
        DrawText(text,0.6,0.7,20);
        text = "#LT p_{T}^{ch} #GT = "; text+=RoundDecimal(CalculateGeometricMean(hUE_partProj[a][p]),3);
        DrawText(text,0.6,0.6,20);
        hUE_partProj[a][p]->Write();
        name = dirName + "/UE_part_" + lohi[a] + "EA_" + ptBinName[p] + ".pdf";
        c0->SaveAs(name,"PDF");
        cout<<hUE_partProj[a][p]->Integral()/AREA<<endl<<CalculateGeometricMean(hUE_partProj[a][p])<<endl;
    
    }
  }

  c0->SetLogy(0);

  int lohiMarker[nEAbins] = {kOpenCircle, 20};
  
  double x1[4] = {10,15,20,30};
  double y1[4] = {0.5,0.93,1.12,1.8};
  TH2D *hScaleMult = new TH2D("hScaleMult","#LT #frac{dN^{ch}}{d#eta d#phi} #GT",3,x1,4,y1);

  double y2[2] = {0.55,0.85};
  TH2D *hScalePt = new TH2D("hScalePt","#LT p_{T}^{ch} #GT",3,x1,1,y2);

  TH1D *hMeanPt[nEAbins], *hMeanMult[nEAbins], *hMeanPt_noJES[nEAbins], *hMeanMult_noJES[nEAbins];






  
  for (int a=0; a<nEAbins; ++a) {
    name = "hMeanPt_" + lohi[a]; 
    hMeanPt[a] = new TH1D(name,"",3,x1);
    name = "hMeanMult_" + lohi[a]; 
    hMeanMult[a] = new TH1D(name,"",3,x1);

    name = "hMeanPt_noJES_" + lohi[a]; 
    hMeanPt_noJES[a] = new TH1D(name,"",3,x1);
    name = "hMeanMult_noJES_" + lohi[a]; 
    hMeanMult_noJES[a] = new TH1D(name,"",3,x1);

    for (int p=0; p<nPtBins; ++p) {
      hMeanPt_noJES[a]->SetBinContent(p+1,CalculateGeometricMean(hUE1D_noJES_dbw[a][p]));
      hMeanMult_noJES[a]->SetBinContent(p+1,hUE1D_noJES_dbw[a][p]->Integral()/AREA);      
      hMeanPt[a]->SetBinContent(p+1,CalculateGeometricMean(hUE1D_dbw[a][p]));
      hMeanMult[a]->SetBinContent(p+1,hUE1D_dbw[a][p]->Integral()/AREA);      
    }
  }


  hScaleMult->SetStats(0);
  hScaleMult->Draw();
  for (int a=0; a<nEAbins; ++a) {
    hMeanMult_noJES[a]->SetMarkerStyle(lohiMarker[a]);
    hMeanMult_noJES[a]->Draw("same p");
  }
  name = "plots/test/meanMult_noJES.pdf";
  c0->SaveAs(name,"PDF");
  
  hScalePt->SetStats(0);
  hScalePt->Draw();
  for (int a=0; a<nEAbins; ++a) {
    hMeanPt_noJES[a]->SetMarkerStyle(lohiMarker[a]);
    hMeanPt_noJES[a]->Draw("same p");
  }
  name = "plots/test/meanPt_noJES.pdf";
  c0->SaveAs(name,"PDF");









  

  double y3[2] = {0.93,1.12};
  TH2D *hScaleRatio = new TH2D("hScaleRatio","",3,x1,1,y3);
  hScaleRatio->SetStats(0);
  
  TH1D *rMeanPt[nEAbins], *rMeanMult[nEAbins];
  for (int a=0; a<nEAbins; ++a) {
    name = "rMeanPt_" + lohi[a]; 
    rMeanPt[a] = (TH1D*)hMeanPt[a]->Clone();
    rMeanPt[a]->Divide(hMeanPt_noJES[a]);
    rMeanPt[a]->SetMarkerStyle(kCircle);
    rMeanPt[a]->SetMarkerColor(kBlack);
    rMeanPt[a]->SetLineColor(kBlack);
    
    name = "rMeanMult_" + lohi[a]; 
    rMeanMult[a] = (TH1D*)hMeanMult[a]->Clone();
    rMeanMult[a]->Divide(hMeanMult_noJES[a]);
    rMeanMult[a]->SetMarkerStyle(kCircle);
    rMeanMult[a]->SetMarkerColor(kBlack);
    rMeanMult[a]->SetLineColor(kBlack);
  }



  TCanvas *c1 = new TCanvas("c1","",600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();
  
  pad1->cd();
  hScalePt->Draw();
  pad2->cd();
  hScaleRatio->GetYaxis()->SetRangeUser(0.93,1.12);
  hScaleRatio->Draw();

  TF1 *unityLine = new TF1("unityLine","1.0",10,30);
  unityLine->SetLineColor(kBlack);    unityLine->SetLineWidth(1);    unityLine->SetLineStyle(kDashed);
  unityLine->Draw("same");


  for (int a=0; a<nEAbins; ++a) { 
    pad1->cd();
    hMeanPt[a]->SetLineColor(kRed);
    hMeanPt[a]->SetMarkerColor(kRed);
    hMeanPt[a]->SetMarkerStyle(lohiMarker[a]);
    hMeanPt_noJES[a]->SetLineColor(kBlue);
    hMeanPt_noJES[a]->SetMarkerColor(kBlue);
    hMeanPt_noJES[a]->Draw("psame");
    hMeanPt[a]->Draw("psame");
    c1->cd();
    pad2->cd();
    rMeanPt[a]->SetStats(0);
    rMeanPt[a]->GetYaxis()->SetRangeUser(0.93,1.12);
    rMeanPt[a]->SetMarkerStyle(lohiMarker[a]);
    rMeanPt[a]->SetAxisRange(0.93,1.12,"Y");
    rMeanPt[a]->Draw("epsame");
    c1->cd();
  }

  saveName = dirName + "/meanPtRatio.pdf";
  c1->SaveAs(saveName,"PDF");



  pad1->cd();
  hScaleMult->GetYaxis()->SetRangeUser(0.5,1.8);
  hScaleMult->SetAxisRange(0.5,1.8,"Y");
  hScaleMult->Draw();
  pad2->cd();
  hScaleRatio->GetYaxis()->SetRangeUser(0.93,1.12);
  hScaleRatio->Draw();
  unityLine->Draw("same");

  pad1->cd();
  hScaleMult->Draw();
  
  for (int a=0; a<nEAbins; ++a) { 
    pad1->cd();
    hMeanMult[a]->SetLineColor(kRed);
    hMeanMult[a]->SetMarkerColor(kRed);
    hMeanMult[a]->SetMarkerStyle(lohiMarker[a]);
    hMeanMult_noJES[a]->SetLineColor(kBlue);
    hMeanMult_noJES[a]->SetMarkerColor(kBlue);
    hMeanMult_noJES[a]->Draw("psame");
    hMeanMult[a]->Draw("psame");
    c1->cd();
    pad2->cd();
    rMeanMult[a]->SetStats(0);
    rMeanMult[a]->GetYaxis()->SetRangeUser(0.93,1.12);
    rMeanMult[a]->SetMarkerStyle(lohiMarker[a]);
    rMeanMult[a]->SetAxisRange(0.93,1.12,"Y");
    rMeanMult[a]->Draw("epsame");
    c1->cd();
  }

  saveName = dirName + "/meanMultRatio.pdf";
  c1->SaveAs(saveName,"PDF");






  
  
  // outFile->cd();
  // for (int a=0; a<nEAbins; ++a) {
  //   hLeadSum[a]->Write();
  //   hUE3Dsum[a]->Write();
  //   hResponseSum[a]->Write();
  //   // for (int i=0; i<55; ++i) {
  //   //   hUE1D[a][i]->Write();
  //   //   hUE1D_dbw[a][i]->Write();
  //   // }
  //   for (int p=0; p<nPtBins; ++p) { hUE_partProj[a][p]->Write(); }

  //   hMeanPt[a]->Write();
  //   hMeanMult[a]->Write();
  //   hMeanPt_noJES[a]->Write();
  //   hMeanMult_noJES[a]->Write();
  // }
  
  for (int a=0; a<nEAbins; ++a) {
    embFile[a]->Close();
    inFile[a]->Close();
  }
  outFile->Close();

  return 0;
}
