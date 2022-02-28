// Veronica Verkest
// December 10, 2021
// to explore the file "forVeronica.root" received from Isaac Mooney on December 8

void macro() {

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TString detSuffix[3] = { "_det_nom", "_det_TS", "_det_TU" };
  TString name, title, saveName;
  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );

  TFile *inFile = new TFile("forVeronica.root","READ");
  // these histograms have the binning {4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 74, 79}
  TH1D *pt_gen_nom = (TH1D*)inFile->Get("pt_gen_nom");  // matched generator-level
  TH1D *pt_det_nom = (TH1D*)inFile->Get("pt_det_nom");  // matched detector-level
  TH1D *pt_det_TS = (TH1D*)inFile->Get("pt_det_TS");  // detector-level with tower scale
  TH1D *pt_det_TU = (TH1D*)inFile->Get("pt_det_TU");  // detector-level with tracking effects
  TH1D *pt_det[3] = { pt_det_nom, pt_det_TS, pt_det_TU };  // put the 3 det-level spectra into an array for easy looping
  TH1D *ratio[3];

  pt_gen_nom->Scale(1./pt_gen_nom->Integral()); // NORMALIZE GEN SPECTRUM TO UNITY
  pt_gen_nom->SetMarkerStyle(20);
  pt_gen_nom->SetStats(0);
  
  for (int i=0; i<3; ++i){
    pt_det[i]->Scale(1./pt_det[i]->Integral()); // NORMALIZE DET SPECTRA TO UNITY
    pt_det[i]->SetMarkerStyle(20);
    pt_det[i]->SetStats(0);
    name = "ratio" + detSuffix[i];
    
    ratio[i] = (TH1D*)pt_det[i]->Clone(name);
    ratio[i]->Divide( pt_gen_nom );
    ratio[i]->GetYaxis()->SetTitle("det / part");
    
    ratio[i]->GetYaxis()->SetRangeUser(-0.2,2.5);
    ratio[i]->Draw();
    name = "plots/" + name + ".pdf";
    can->SaveAs(name,"PDF");
  }

  ratio[0]->GetYaxis()->SetRangeUser(0.,1.5);
  ratio[0]->Draw("PLCPMC");
  for (int i=1; i<3; ++i){
    ratio[i]->GetYaxis()->SetRangeUser(0.,1.5);
    ratio[i]->Draw("PLCPMCSAME");
  }
  can->BuildLegend();
  can->SaveAs("plots/ptSpectraRATIO.pdf","PDF");
  
  can->SetLogy();
  pt_gen_nom->Draw("PLCPMC");
  for (int i=0; i<3; ++i){ pt_det[i]->Draw("PLCPMCSAME"); }
  can->BuildLegend();
  can->SaveAs("plots/ptSpectra.pdf","PDF");

  TFile *outFile = new TFile("JES_ratios_from_Isaac.root","RECREATE");
  pt_gen_nom->Write();
  for (int i=0; i<3; ++i){
    pt_det[i]->Write();
    ratio[i]->Write();
  }
  
}
