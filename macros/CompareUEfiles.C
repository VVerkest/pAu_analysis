// Veronica Verkest
// January 11, 2021

void CompareUEfiles(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile *fVar = new TFile("out/UE/pAuHTjetUE_hiEA_uncorrected.root","READ");
  TFile *fHalf = new TFile("out/UE/pAuHTjetUE_halfGeVbins_hiEA.root","READ");

  TH3D *hEastVar = (TH3D*)fVar->Get("hChgUE_eastEtaJet");
  TH3D *hEastHalf = (TH3D*)fHalf->Get("hChgUE_eastEtaJet");
  TH3D *hLeadVar = (TH3D*)fVar->Get("hLeadPt_eastEtaJet");
  TH3D *hLeadHalf = (TH3D*)fHalf->Get("hLeadPt_eastEtaJet");

  // hEastHalf->GetYaxis()->SetRangeUser(0.,0.5);
  // hEastVar->GetYaxis()->SetRangeUser(0.2,0.5);
  
  TH1D *hUEVar = (TH1D*)hEastVar->ProjectionY("hUEVar");
  double varInt = hUEVar->Integral();
  TH1D *hUEHalf = (TH1D*)hEastHalf->ProjectionY("hUEHalf");
  double halfInt = hUEHalf->Integral();

  double varMean = 0; double halfMean = 0;

  hUEVar->Scale(1./hUEVar->Integral());
  for (int i=1; i<=hUEVar->GetNbinsX(); ++i) {
    varMean += hUEVar->GetBinContent(i)*hUEVar->GetBinCenter(i);
    cout<<hUEVar->GetBinContent(i)<<"  \t"<<hUEVar->GetBinCenter(i)<<"  \t"<<hUEVar->GetBinContent(i)*hUEVar->GetBinCenter(i)<<endl;
    double width = hUEVar->GetXaxis()->GetBinWidth(i);
    hUEVar->SetBinContent( i, hUEVar->GetBinContent(i)/width );
  }
  hUEVar->Scale(hUEVar->GetEntries()/hLeadVar->GetEntries());
  // hUEVar->Scale(varInt/hUEVar->Integral());
  // hUEVar->Scale(1./hUEVar->GetEntries());

  cout<<endl<<endl;
  
  hUEHalf->Scale(1./hUEHalf->Integral());
  for (int i=1; i<=hUEHalf->GetNbinsX(); ++i) {
    halfMean += hUEHalf->GetBinContent(i)*hUEHalf->GetBinCenter(i);
    cout<<hUEHalf->GetBinContent(i)<<" \t"<<hUEHalf->GetBinCenter(i)<<" \t"<<hUEHalf->GetBinContent(i)*hUEHalf->GetBinCenter(i)<<endl;
    double width = hUEHalf->GetXaxis()->GetBinWidth(i);
    hUEHalf->SetBinContent( i, hUEHalf->GetBinContent(i)/width );
  }
  hUEHalf->Scale(hUEHalf->GetEntries()/hLeadHalf->GetEntries());
  // hUEHalf->Scale(halfInt/hUEHalf->Integral());
  // hUEHalf->Scale(1./hUEHalf->GetEntries());

  cout<<varMean<<endl<<halfMean<<endl;
  
  hUEVar->Draw();
  new TCanvas;
  hUEHalf->Draw();
}
