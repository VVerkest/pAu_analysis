// Veronica Verkest
// February 15, 2021

void InclusiveJEScorrection(){

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

  TH1D *hUE1D[55];

  TCanvas *c0 = new TCanvas("c0");
  c0->SetLogy();


  //  UEpT_varBinning_allLeadPt
  //  projects and saves the UE pT distribution (scaled by Njets) for a part-level leading jet of 4<=pT<=59 GeV
  for (int i=0; i<55; ++i) {
    int ptLo = i+4, ptHi = i+5, binno = i+1;
    name = "hUE1D_"; name += ptLo; name += "to"; name += ptHi; name += "GeV";
    hUE1D[i] = (TH1D*)hUE2D->ProjectionY(name,binno,binno);
    hUE1D[i]->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
    hUE1D[i]->SetMarkerColor(925 + (5*i));
    hUE1D[i]->SetLineColor(925 + (5*i));
    // hUE1D[i]->SetMarkerColor(1170 - (5*i));
    // hUE1D[i]->SetLineColor(1170 - (5*i));
    hUE1D[i]->SetMarkerStyle(marker[i]);
    hUE1D[i]->SetMarkerSize(1.);
    hUE1D[i]->SetStats(0);
    hUE1D[i]->GetXaxis()->CenterTitle();
      hUE1D[i]->GetXaxis()->SetTitleOffset(1.);
    if ( hUE1D[i]->Integral()!=0 ) {
      hUE1D[i]->Scale(1./hLead->Integral(binno,binno));
      hUE1D[i]->SetAxisRange(0.2,17.,"X");
      hUE1D[i]->SetAxisRange(0.000001,5.,"Y");
      hUE1D[i]->Draw("SAME");
      hUE1D[i]->Write();
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


  //  UEpT_varBinning_10_30
  //  UE pT distribution (scaled by Njets) for det-level jets 10-30 GeV
  hUE2D->GetXaxis()->SetRangeUser(10.,30.);
  name = "hUE1D_10_30GeV";
  TH1D *hUE1D_10_30 = (TH1D*)hUE2D->ProjectionY(name);
  hUE1D_10_30->Scale(1./hLead->Integral(7,26));
  hUE1D_10_30->GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
  hUE1D_10_30->Draw();
  saveName = dirName + "/UEpT_varBinning_10_30.pdf";
  c0->SaveAs(saveName,"PDF");
  cout<<"mean: "<<hUE1D_10_30->GetMean()<<endl;
  cout<<"integral: "<<hUE1D_10_30->Integral()<<endl<<endl;


  //  UEpT_10_30.pdf
  //  UE pT distribution (scaled by Njets, adjusted for bin width, and re-scaled to preserve integral) for det-level jets 10-30 GeV
  TH1D *hUE1D_10_30_dbw = (TH1D*)hUE1D_10_30->Clone("hUE1D_10_30_dbw");
  double geoMean = 0;
  double totalContent = 0;
  for (int i=1; i<=hUE1D_10_30_dbw->GetNbinsX(); ++i) {
    geoMean += hUE1D_10_30_dbw->GetBinCenter(i)*hUE1D_10_30_dbw->GetBinContent(i);///hUE1D_10_30_dbw->GetBinWidth(i);
    totalContent += hUE1D_10_30_dbw->GetBinContent(i);
    hUE1D_10_30_dbw->SetBinContent(i,hUE1D_10_30_dbw->GetBinContent(i)/hUE1D_10_30_dbw->GetBinWidth(i));
  }
  geoMean /= totalContent;
  hUE1D_10_30_dbw->Scale(hUE1D_10_30->Integral()/hUE1D_10_30_dbw->Integral());
  cout<<"geometric mean: "<<geoMean<<endl;
  cout<<"integral: "<<hUE1D_10_30_dbw->Integral()<<endl<<endl;
  hUE1D_10_30_dbw->Draw();
  saveName = dirName + "/UEpT_10_30.pdf";
  c0->SaveAs(saveName,"PDF");
  

  embFile->Close();
  inFile->Close();

}
