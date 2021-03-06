//  Veronica Verkest		February 14, 2020
//  ~macro for normalizing BBCEsum distribution & calculating EA deciles
//  ~uses output file from EAdistribution.cxx

void EAdist(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TFile *EAfile = new TFile("out/EAdist/pAu_2015_EAdist.root","UPDATE");

  TH1D *hBBCdist = (TH1D*) EAfile->Get("hEAdist");
  hBBCdist->Scale(1./hBBCdist->Integral());
  
  TTree *eatree = (TTree*) EAfile->Get("EAtree");

  double ps, temp;
  int RunID, EventID;		double Vz, BbcAdcSumEast, ZDCx;
  eatree->SetBranchAddress( "RunID", &RunID );		eatree->SetBranchAddress( "EventID", &EventID );
  eatree->SetBranchAddress( "Vz", &Vz );		eatree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );
  eatree->SetBranchAddress( "ZDCx", &ZDCx );

  int nBins = hBBCdist->GetNbinsX();
  TH1D *hBBCint = new TH1D("hBBCint","BBCEsum Distribution Integral",280000,0,70000);

  temp = 0;
  ps = 0;
  for (int i=0;i<nBins;++i) {
    temp = hBBCdist->Integral(i,i);  //partial sum
    ps += temp;
    hBBCint->SetBinContent(i,ps);
  }

  int bin = 0;
  double tenth = (double) 1.0/10.0;
  double deciles[10];
  for (int i=0;i<10;++i) {
    double min = (i+1)*tenth;
    bin = hBBCint->FindFirstBinAbove(min);
    deciles[i] = hBBCint->GetBinCenter(bin);
    hBBCdist->SetBinError(bin,9999);
  }

  bin = 64000*4;
  hBBCdist->SetBinContent(bin,0.000001);
  hBBCdist->SetBinError(bin,9999);
  
  fstream file("src/EAdeciles.txt", fstream::in | fstream::out | fstream::trunc);
  if (!file) {    cerr << "Error in creating file!!!" << endl; exit(1);  }
  else {    cout << "src/EAdeciles.txt created successfully." << endl;  }
  
  file << 0.00 << ", ";
  for (int i=0;i<9;++i) {file << deciles[i] << ", ";}
  file << 64000;
  
  file.close();
  cout << "Closed bad_towers_pAu2015.list" << endl;

  TCanvas * c0 = new TCanvas( "c0" , "" ,0 ,23 ,1280 ,350 );              // CANVAS
  c0->SetLogy();
  hBBCdist->GetYaxis()->SetRangeUser(0.000001, 0.01);
  hBBCdist->SetTitle("");
  hBBCdist->SetStats(0);
  hBBCdist->Rebin(200);
  hBBCdist->Draw();
  c0->SaveAs("plots/EAdist.pdf","PDF");
  
  //EAfile->Close();
}
