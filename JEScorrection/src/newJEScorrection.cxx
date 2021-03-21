// Veronica Verkest
// February 24, 2020

#include "params.hh"
#include "funcs.hh"

using namespace std;
using namespace Analysis;

double CalculateGeometricMean( TH1D* h ) {
  double geoMean = 0;
  double totalContent = 0;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    geoMean += h->GetBinCenter(i)*h->GetBinContent(i)*h->GetBinWidth(i);
    totalContent += h->GetBinContent(i)*h->GetBinWidth(i);
  }
  geoMean /= totalContent;
  return geoMean;
}


TH1D* GenerateFractionalContribution( TH2D* response, double pTlo, double pThi, TString dir ) {
  int lo = (int)pTlo;     int hi = (int)pThi;

  TCanvas *can = new TCanvas();
  can->SetLogy();

  response->GetXaxis()->SetRangeUser(pTlo,pThi);
  TString name = "hResponse_"; name+=lo; name+="_"; name+=hi;
  TH1D *response1D = (TH1D*)response->ProjectionY(name);
  response1D->Scale(1./response1D->Integral());
  TString title = "frac. contribution to a "; title+=lo; title+="-"; title+=hi; title+=" GeV part. jet";
  response1D->GetYaxis()->SetTitle(title);
  response1D->Draw();
  // response1D->Write();
  TString saveName = dir + "/response_"; saveName+=lo; saveName+="_"; saveName+=hi; saveName+=".pdf";
  can->SaveAs(saveName,"PDF");
  response->GetXaxis()->SetRange(1, -1);

  can->Destructor();
  
  return response1D;
}


void ProjectScaleAndSaveUE1D( TH2D* UE2D, TH1D &UE1D, TH1D* lead, TH1D &UE1D_dbw, double pTlo, double pThi, TString dir ){

  TString name, saveName;
  
  TCanvas *can = new TCanvas();
  can->SetLogy();
  
  int lo = (int)pTlo;     int hi = (int)pThi;
    
  //  UEpT_varBinning_pTlo_pThi
  //  UE pT distribution (scaled by Njets) for det-level jets pTlo-pThi GeV
  UE2D->GetXaxis()->SetRangeUser( pTlo, pThi );
  name = "hUE1D_"; name+=pTlo; name+="_"; name+=pThi; name +="GeV";
  TH1D* htemp = (TH1D*)UE2D->ProjectionY(name);
  UE1D = *htemp;
  int binlo = lo-3; int binhi = hi-4;
  UE1D.Scale(1./lead->Integral(lo-3,hi-4));
  UE1D.GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
  if (UE1D.GetEntries()>0) {
    UE1D.SetAxisRange(0.0000001,1.,"Y");
    UE1D.Draw();
    saveName = dir + "/UEpT_varBinning_"; saveName+=pTlo; saveName+="_"; saveName+=pThi; saveName +=".pdf";
    can->SaveAs(saveName,"PDF");
    cout<<"mean: "<<UE1D.GetMean()<<endl;
    cout<<"integral: "<<UE1D.Integral()<<endl<<endl;
  }

  //  UEpT_pTlo_pThi
  //  UE pT distribution (scaled by Njets, adjusted for bin width, and re-scaled to preserve integral) for det-level jets pTlo-pThi GeV
  name = "hUE1D_"; name+=pTlo; name+="_"; name+=pThi; name +="_dbw";

  UE1D_dbw = *htemp;
  binlo = lo-3; binhi = hi-4;
  UE1D_dbw.Scale(1./lead->Integral(lo-3,hi-4));
  UE1D_dbw.GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
  
  for (int i=1; i<=UE1D_dbw.GetNbinsX(); ++i) { UE1D_dbw.SetBinContent(i,UE1D_dbw.GetBinContent(i)/UE1D_dbw.GetBinWidth(i)); }
  UE1D_dbw.Scale(UE1D.Integral()/UE1D_dbw.Integral());
  if (UE1D_dbw.GetEntries()>0) {
    cout<<"geometric mean: "<< CalculateGeometricMean(&UE1D_dbw) <<endl;
    cout<<"integral: "<<UE1D_dbw.Integral()<<endl<<endl;
    UE1D_dbw.SetAxisRange(0.0000001,1.,"Y");
    UE1D_dbw.Draw();
    saveName = dir + "/UEpT_"; saveName+=pTlo; saveName+="_"; saveName+=pThi; saveName +="_dbw.pdf";
    can->SaveAs(saveName,"PDF");
  }

  // UE1D.Write();
  // UE1D_dbw.Write();
  can->Destructor();
}


TString RoundDecimal(double val, int nDigitsAfterDecimal ) {
  double n = val;
  n += 5/pow(10,nDigitsAfterDecimal+1);
  TString rounded = "";
  rounded += n;
  int nDigitsBeforeDecimal = 1 + rounded.First(".");  
  int input = nDigitsBeforeDecimal + nDigitsAfterDecimal;
  rounded = rounded(0,input);
  return rounded;
}


TH1D* WeightAndSumByFC( TH1D* FC, TH1D* UE1D[55] ){

  const int bins = 15;
  double binEdge[bins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };
  TH1D* UE_part = new TH1D("UE_part","",bins,binEdge);

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
  TString dirName = "plots/test/";
  // TString dirName = "plots/";

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();


  // for (int a=0; a<nEAbins; ++a) {
  //   dirname += lohi[a]; dirname += "/";


  // }

  
  
  TFile* inFile = new TFile("../out/UE/pAuHTjetUE_inclusive.root", "READ");
  TH2D *hUE2D = (TH2D*)inFile->Get("hLeadPtVsUEpT");
  TH1D *hLead = (TH1D*)inFile->Get("hLead");

  TFile* embFile = new TFile("../embedding/out/sim/pAu2015embedding_inclusive.root", "READ");
  TH2D *hResponse = (TH2D*)embFile->Get("hResponse");
  
  TFile* outFile = new TFile("out/test.root", "RECREATE");
  outFile->cd();
  hUE2D->Write();
  hLead->Write();
  hResponse->Write();

  TH1D *hUE;//[nPtBins][nEtaBins][nEAbins];

  
  TH1D *hUE1D[55], *hUE1D_dbw[55];
  for (int i=0; i<55; ++i) {
    name = "hUE1D_"; name +=(i+4); name+="_"; name+=(i+5);
    hUE1D[i] = new TH1D(name,"",ybins,ybinEdge);
    name += "_dbw";
    hUE1D_dbw[i] = new TH1D(name,"",ybins,ybinEdge);

    ProjectScaleAndSaveUE1D( hUE2D, *hUE1D[i], hLead, *hUE1D_dbw[i], i+4., i+5., dirName );
    hUE1D[i]->Write();
    hUE1D_dbw[i]->Write();
  }


  TH1D *hResponseProj[nPtBins]; TH1D* hUE_partProj[nPtBins];

  for (int p=0; p<nPtBins; ++p) {

    name = "hResponse" + ptBinName[p];
    hResponseProj[p] = GenerateFractionalContribution( hResponse, ptLo[p], ptHi[p], dirName );
    hResponseProj[p]->SetName(name);

    name = "hUE_part" + ptBinName[p];
    hUE_partProj[p] = (TH1D*) WeightAndSumByFC( hResponseProj[p], hUE1D_dbw );
    hUE_partProj[p]->SetName(name);
    hUE_partProj[p]->SetAxisRange(0.0000001,1,"Y");
    hUE_partProj[p]->Draw();
    
    TString text = "#LT #frac{dN^{ch}}{d#eta d#phi} #GT = "; text+=RoundDecimal(hUE_partProj[p]->Integral()/AREA,3);
    DrawText(text,0.6,0.7,20);
    text = "#LT p_{T}^{ch} #GT = "; text+=RoundDecimal(CalculateGeometricMean(hUE_partProj[p]),3);
    DrawText(text,0.6,0.6,20);
    hUE_partProj[p]->Write();
    name = dirName + "/UE_part_" + ptBinName[p] + ".pdf";
    c0->SaveAs(name,"PDF");
    cout<<hUE_partProj[p]->Integral()/AREA<<endl<<CalculateGeometricMean(hUE_partProj[p])<<endl;
    
  }

   
  embFile->Close();
  inFile->Close();
  outFile->Write();
  outFile->Close();

  return 0;
}
