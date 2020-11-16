// Veronica Verkest
// September 23, 2020

void FitTrackingEfficiencyHistograms() {

  gStyle->SetOptFit(1111);
  TString name;
  
  TFile *efficFile = new TFile( "../src/trackeffic_new.root", "READ" );
  TH1D *hEff[10][3];
  TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);

  TCanvas *can = new TCanvas;

  for (int i=1;i<11;++i){
    name = "eff_s_bin_2_3_bbc__"; name += i; name += "_"; name += i; name += "_eta";
    hEff[i-1][1] = (TH1D*)efficFile->Get(name);
    hEff[i-1][1]->Fit("eff","","",0.,15.);
    hEff[i-1][1]->SetAxisRange(0.0,25.0,"X");
    hEff[i-1][1]->SetAxisRange(0.0,1.3,"Y");
    hEff[i-1][1]->Draw();
    name = "plots/efficFit/loEAefficiency_etaBin"; name += i; name += ".pdf";
    can->SaveAs(name,"PDF");

    // TF1 *fit = hEff[i-1][1]->GetFunction("eff");
    // cout<<"lo EA, eta bin "<<i<<": chi^2 = "<<fit->GetChisquare()<<"\t Ndof = "<<fit->GetNDF()<<"\t chi^2/Ndof = "<<fit->GetChisquare()/fit->GetNDF()<<endl;
  }

  for (int i=1;i<11;++i){
    name = "eff_s_bin_8_10_bbc__"; name += i; name += "_"; name += i; name += "_eta";
    hEff[i-1][2] = (TH1D*)efficFile->Get(name);
    hEff[i-1][2]->Fit("eff","","",0.,15.);
    hEff[i-1][2]->SetAxisRange(0.0,25.0,"X");
    hEff[i-1][2]->SetAxisRange(0.0,1.3,"Y");
    hEff[i-1][2]->Draw();
    name = "plots/efficFit/hiEAefficiency_etaBin"; name += i; name += ".pdf";
    can->SaveAs(name,"PDF");
    // TF1 *fit = hEff[i-1][2]->GetFunction("eff");
    // cout<<"hi EA, eta bin "<<i<<": chi^2 = "<<fit->GetChisquare()<<"\t Ndof = "<<fit->GetNDF()<<"\t chi^2/Ndof = "<<fit->GetChisquare()/fit->GetNDF()<<endl;
  }

  cout<<endl;
  
  TString EAstring[3] = {"all","hi","lo"};
  
  for (int j=1;j<3;++j){
    for (int i=1;i<11;++i){
      TF1 *fit = hEff[i-1][j]->GetFunction("eff");
      cout<<EAstring[j]<<" EA, eta bin "<<i<<": chi^2 =  "<<fit->GetChisquare()<<"\t Ndof = "<<fit->GetNDF()<<"\t chi^2/Ndof = "<<fit->GetChisquare()/fit->GetNDF()<<endl;
    }
      cout<<endl;
  }
  
  // for (int i=1;i<11;++i){
  //   name = "eff_s_bin_1_10_bbc__"; name += i; name += "_"; name += i; name += "_eta";
  //   hEff[i-1][0] = (TH1D*)efficFile->Get(name);
  //   hEff[i-1][0]->Fit("eff");
  //   hEff[i-1][0]->SetAxisRange(0.0,25.0,"X");
  //   hEff[i-1][0]->SetAxisRange(0.0,1.3,"Y");
  //   hEff[i-1][0]->Draw();
  //   name = "plots/efficFit/allEAefficiency_etaBin"; name += i; name += ".pdf";
  //   can->SaveAs(name,"PDF");
  // }

  // for (int i=1;i<11;++i){
  //   name = "eff_s_bin_2_3_bbc__"; name += i; name += "_"; name += i; name += "_eta";
  //   hEff[i-1][1] = (TH1D*)efficFile->Get(name);
  //   hEff[i-1][1]->Fit("eff");
  //   hEff[i-1][1]->SetAxisRange(0.0,25.0,"X");
  //   hEff[i-1][1]->SetAxisRange(0.0,1.3,"Y");
  //   hEff[i-1][1]->Draw();
  //   name = "plots/efficFit/loEAefficiency_etaBin"; name += i; name += ".pdf";
  //   can->SaveAs(name,"PDF");
  // }

  // for (int i=1;i<11;++i){
  //   name = "eff_s_bin_8_10_bbc__"; name += i; name += "_"; name += i; name += "_eta";
  //   hEff[i-1][2] = (TH1D*)efficFile->Get(name);
  //   hEff[i-1][2]->Fit("eff");
  //   hEff[i-1][2]->SetAxisRange(0.0,25.0,"X");
  //   hEff[i-1][2]->SetAxisRange(0.0,1.3,"Y");
  //   hEff[i-1][2]->Draw();
  //   name = "plots/efficFit/hiEAefficiency_etaBin"; name += i; name += ".pdf";
  //   can->SaveAs(name,"PDF");
  // }
  
}
