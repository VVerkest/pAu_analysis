//  Veronica Verkest		February 14, 2020
//  ~macro for normalizing BBCEsum distribution & calculating EA deciles
//  ~uses output file from EAdistribution.cxx

void EAdist(){

  TFile *EAfile = new TFile("out/EAdist/pAu_2015_EAdist.root","READ");

  TH1D *hBBCdist = (TH1D*) EAfile->Get("hEAdist");
  hBBCdist->Scale(1./hBBCdist->Integral("width"));
  
  TTree *eatree = (TTree*) EAfile->Get("EAtree");

  int RunID, EventID;		double Vz, BbcAdcSumEast, ps;
  eatree->SetBranchAddress( "RunID", &RunID );	eatree->SetBranchAddress( "EventID", &EventID );
  eatree->SetBranchAddress( "Vz", &Vz );		eatree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );

  int nBins = hBBCdist->GetNbinsX();
  TH1D *hBBCint = new TH1D("hBBCint","BBCEsum Distribution Integral",280000,0,70000);

  for (int i=0;i<nBins;++i) {
    ps = hBBCdist->Integral(0,i);  //partial sum
    hBBCint->SetBinContent(i,ps);
  }

  double tenth = 0.10000000000000;
  double deciles[10];
  for (int i=0;i<10;++i) {
    cout<<"Finding Decile: "<<i<<"0-"<<(i+1)<<"0%"<<endl;
    double min = (i+1)*tenth;
    int bin = hBBCint->FindFirstBinAbove(min);
    deciles[i] = hBBCint->GetBinCenter(bin);
  }

  fstream file("src/EAdeciles.txt", fstream::in | fstream::out | fstream::trunc);
  if (!file) {    cerr << "Error in creating file!!!" << endl; exit(1);  }
  else {    cout << "src/EAdeciles.txt created successfully." << endl;  }
  
  file << 0.00 << ", ";
  for (int i=0;i<9;++i) {file << deciles[i] << ", ";}
  file << 64000; //writes that element
  
  file.close();
  cout << "Closed bad_towers_pAu2015.list" << endl;

  
}
