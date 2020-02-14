//  Veronica Verkest		February 14, 2020
//  macro for normalizing BBCEsum distribution & calculating EA deciles

void EAdist(){

  TFile *EAfile = new TFile("out/EAdist/pAu_2015_EAdist.root","READ");

  TH1D *hBBCdist = (TH1D*) EAfile->Get("hEAdist");
  hBBCdist->Scale(1./hBBCdist->Integral());
  
  TTree *eatree = (TTree*) EAfile->Get("EAtree");

  int RunID, EventID;		double Vz, BbcAdcSumEast, ps;
  eatree->SetBranchAddress( "RunID", &RunID );	eatree->SetBranchAddress( "EventID", &EventID );
  eatree->SetBranchAddress( "Vz", &Vz );		eatree->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast );

  int nBins = hEAdist->GetNbinsX();
  TH1D *hBBCint = new TH1D("hBBCint","BBCEsum Dis. Integral",70000,0,70000);

  for (int i=0;i<nBins;++i) {
    ps = hBBCdist->Integral(0,i);  //partial sum
    hBBCint->SetBinContent(i,ps);
  }

}
