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


TH1D GenerateFractionalContribution( TH2D* response, double pTlo, double pThi, TString dir ) {
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
  
  return *response1D;
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

  // const int ybins = 15;
  // double ybinEdge[ybins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;

  TString dirName = "plots/test";

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

  // const int n_bins = 3;
  // double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  // const int n_ybins = 3;
  // double y_bin_edge[n_ybins+1] = { 0.55,0.65,0.75,0.85 };
  // TH2D *meanPt_hscale = new TH2D("meanPt_hscale",";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge);
  // double y_bin_edge2[n_ybins+1] = { 0.5,1.2,1.5,1.8 };
  // TH2D *nCh_hscale = new TH2D("nCh_hscale",";(GeV);#LT#frac{dN_{ch}}{d#eta d#phi}#GT (GeV)",n_bins,bin_edge,n_ybins,y_bin_edge2);
  TH1D *hUE;//[nPtBins][nEtaBins][nEAbins];

  
  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();


  TH1D *hUE1D_stack[55];

  //  UEpT_varBinning_allLeadPt
  //  projects and saves the UE pT distribution (scaled by Njets) for a part-level leading jet of 4<=pT<=59 GeV
  for (int i=0; i<55; ++i) {
    int ptLo = i+4, ptHi = i+5, binno = i+1;
    name = "hUE1D_stack_"; name += ptLo; name += "to"; name += ptHi; name += "GeV";
    hUE1D_stack[i] = (TH1D*)hUE2D->ProjectionY(name,binno,binno);
    hUE1D_stack[i]->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
    hUE1D_stack[i]->SetMarkerColor(925 + (5*i));
    hUE1D_stack[i]->SetLineColor(925 + (5*i));
    // hUE1D_stack[i]->SetMarkerColor(1170 - (5*i));
    // hUE1D_stack[i]->SetLineColor(1170 - (5*i));
    hUE1D_stack[i]->SetMarkerStyle(marker[i]);
    hUE1D_stack[i]->SetMarkerSize(1.);
    hUE1D_stack[i]->SetStats(0);
    hUE1D_stack[i]->GetXaxis()->CenterTitle();
    hUE1D_stack[i]->GetXaxis()->SetTitleOffset(1.);
    if ( hUE1D_stack[i]->Integral()!=0 ) {
      hUE1D_stack[i]->Scale(1./hLead->Integral(binno,binno));
      hUE1D_stack[i]->SetAxisRange(0.2,17.,"X");
      hUE1D_stack[i]->SetAxisRange(0.000001,5.,"Y");
      hUE1D_stack[i]->Draw("SAME");
      // hUE1D_stack[i]->Write();
    }
  }

  c0->BuildLegend(0.83,0.,1.,1.);
  saveName = dirName; saveName += "/UEpT_varBinning_allLeadPt.pdf";
  c0->SaveAs(saveName,"PDF");


  //  Response_10_30
  //  for a part-level leading jet 10-30 GeV, this is the (normalized to unity) fractional contibution of det-level leading jets by pT
  hResponse->GetXaxis()->SetRangeUser(10.,30.);
  TH1D *hResponse_10_30 = (TH1D*)hResponse->ProjectionY("hResponse_10_30");
  hResponse_10_30->Scale(1./hResponse_10_30->Integral());
  // c0->SetLogy(0);
  hResponse_10_30->GetYaxis()->SetTitle("frac. contribution to a 10-30 GeV part. jet");
  hResponse_10_30->Draw();
  hResponse_10_30->Write();
  saveName = dirName + "/Response_10_30.pdf";
  c0->SaveAs(saveName,"PDF");
  hResponse->GetXaxis()->SetRange(1, -1);


  //  UEpT_varBinning_10_30
  //  UE pT distribution (scaled by Njets) for det-level jets 10-30 GeV
  hUE2D->GetXaxis()->SetRangeUser(10.,30.);
  name = "hUE1D_10_30GeV";
  TH1D *hUE1D_10_30 = (TH1D*)hUE2D->ProjectionY(name);
  int binlo = 10-3; int binhi = 30-4;
  hUE1D_10_30->Scale(1./hLead->Integral(7,26));
  hUE1D_10_30->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
  hUE1D_10_30->Draw();
  saveName = dirName + "/UEpT_varBinning_10_30.pdf";
  c0->SaveAs(saveName,"PDF");
  // cout<<"mean: "<<hUE1D_10_30->GetMean()<<endl;
  // cout<<"integral: "<<hUE1D_10_30->Integral()<<endl<<endl;


  //  UEpT_10_30
  //  UE pT distribution (scaled by Njets, adjusted for bin width, and re-scaled to preserve integral) for det-level jets 10-30 GeV
  TH1D *hUE1D_10_30_dbw = (TH1D*)hUE1D_10_30->Clone("hUE1D_10_30_dbw");
  for (int i=1; i<=hUE1D_10_30_dbw->GetNbinsX(); ++i) { hUE1D_10_30_dbw->SetBinContent(i,hUE1D_10_30_dbw->GetBinContent(i)/hUE1D_10_30_dbw->GetBinWidth(i)); }
  hUE1D_10_30_dbw->Scale(hUE1D_10_30->Integral()/hUE1D_10_30_dbw->Integral());
  // cout<<"geometric mean: "<< CalculateGeometricMean(hUE1D_10_30_dbw) <<endl;
  // cout<<"integral: "<<hUE1D_10_30_dbw->Integral()<<endl<<endl;
  hUE1D_10_30_dbw->Draw();
  saveName = dirName + "/UEpT_10_30.pdf";
  c0->SaveAs(saveName,"PDF");






  
  TH1D *hUE1D[55], *hUE1D_dbw[55];
  for (int i=0; i<55; ++i) {
    name = "hUE1D_"; name +=(i+4); name+="_"; name+=(i+5);
    hUE1D[i] = new TH1D(name,"",ybins,ybinEdge);
    name += "_dbw";
    hUE1D_dbw[i] = new TH1D(name,"",ybins,ybinEdge);

    ProjectScaleAndSaveUE1D( hUE2D, *hUE1D[i], hLead, *hUE1D_dbw[i], i+4., i+5., dirName );
    name = "hUE1D_"; name +=(i+4); name+="_"; name+=(i+5);
    // hUE1D[i]->SetName(name);
    hUE1D[i]->Write();
    name += "_dbw";
    // hUE1D_dbw[i]->SetName(name);
    hUE1D_dbw[i]->Write();
  }



  TH1D hResponse_10_15 = GenerateFractionalContribution( hResponse, 10., 15., dirName );
  TH1D* hUE_part_10_15 = (TH1D*) WeightAndSumByFC( &hResponse_10_15, hUE1D_dbw );
  hUE_part_10_15->SetName("hUE_part_10_15");
  hUE_part_10_15->SetAxisRange(0.0000001,1,"Y");
  hUE_part_10_15->Draw();
  TString text = "#LT #frac{dN^{ch}}{d#eta d#phi} #GT = "; text+=RoundDecimal(hUE_part_10_15->Integral()/AREA,3);
  DrawText(text,0.6,0.7,20);
  text = "#LT p_{T}^{ch} #GT = "; text+=RoundDecimal(CalculateGeometricMean(hUE_part_10_15),3);
  DrawText(text,0.6,0.6,20);
  hUE_part_10_15->Write();
  name = dirName + "/UE_part_10_15.pdf";
  c0->SaveAs(name,"PDF");
  cout<<hUE_part_10_15->Integral()/AREA<<endl<<CalculateGeometricMean(hUE_part_10_15)<<endl;

  

  TH1D hResponse_15_20 = GenerateFractionalContribution( hResponse, 15., 20., dirName );
  TH1D* hUE_part_15_20 = (TH1D*) WeightAndSumByFC( &hResponse_15_20, hUE1D_dbw );
  hUE_part_15_20->SetName("hUE_part_15_20");
  hUE_part_15_20->SetAxisRange(0.0000001,1,"Y");
  hUE_part_15_20->Write();
  hUE_part_15_20->Draw();
  text = "#LT #frac{dN^{ch}}{d#eta d#phi} #GT = "; text+=RoundDecimal(hUE_part_15_20->Integral()/AREA,3);
  DrawText(text,0.6,0.7,20);
  text = "#LT p_{T}^{ch} #GT = "; text+=RoundDecimal(CalculateGeometricMean(hUE_part_15_20),3);
  DrawText(text,0.6,0.6,20);  name = dirName + "/UE_part_15_20.pdf";
  c0->SaveAs(name,"PDF");
  cout<<hUE_part_15_20->Integral()/AREA<<endl<<CalculateGeometricMean(hUE_part_15_20)<<endl;

  

  TH1D hResponse_20_30 = GenerateFractionalContribution( hResponse, 20., 30., dirName );
  TH1D* hUE_part_20_30 = (TH1D*) WeightAndSumByFC( &hResponse_20_30, hUE1D_dbw );
  hUE_part_20_30->SetName("hUE_part_20_30");
  hUE_part_20_30->SetAxisRange(0.0000001,1,"Y");
  hUE_part_20_30->Write();
  hUE_part_20_30->Draw();
  text = "#LT #frac{dN^{ch}}{d#eta d#phi} #GT = "; text+=RoundDecimal(hUE_part_20_30->Integral()/AREA,3);
  DrawText(text,0.6,0.7,20);
  text = "#LT p_{T}^{ch} #GT = "; text+=RoundDecimal(CalculateGeometricMean(hUE_part_20_30),3);
  DrawText(text,0.6,0.6,20);  name = dirName + "/UE_part_20_30.pdf";
  c0->SaveAs(name,"PDF");
  cout<<hUE_part_20_30->Integral()/AREA<<endl<<CalculateGeometricMean(hUE_part_20_30)<<endl;
  
  
  embFile->Close();
  inFile->Close();
  outFile->Write();
  outFile->Close();

  return 0;
}
