// Veronica Verkest
// October 20, 2021

#include "params.hh"
#include "funcs.hh"
#include <ROOT/TIOFeatures.hxx>
using namespace std;
using namespace Analysis;



void TrackingEfficiencyTEST( TH2D* h_ChgUE2D[55], TH2D* h_ChgUE2D_corr[55], TFile *effic_File, TString ea_string, TString dir_name, double EfficShift = 0. ){
  //  X: UE pT,   Y: UE eta
  TH1D *hEffic;
  TString name, saveName, bbcBins;
  double ptVal, oldVal, oldErr, effic, efficErr, relErr, newVal, newErr;

  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
  can->SetLogz();
  gStyle->SetOptFit(1);

  // TF1 *eff = new TF1("eff","((200.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  // TF1 *eff = new TF1("eff","(([2]+log(x/100.0)))*exp([0]+[1]*x)",0.2,15.0);
  // TF1 *eff = new TF1("eff","(1+log(x/200.))*(exp([0]+[1]*x))+[2]",0.2,15.0);
  TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  
  if ( ea_string=="lo" ) { bbcBins = "1_3"; }
  else if ( ea_string=="hi" ) { bbcBins = "8_10"; }
  else { std::cerr<< "invalid EA string provided!" <<std::endl; }
    
  for (int iy=1; iy<zbins+1; ++iy){ // loop over UE eta bins

    name = "eff_s_bin_" + bbcBins + "_bbc__"; name += iy; name += "_"; name += iy; name += "_eta";
    TH1D *hEffic = (TH1D*)effic_File->Get(name);

    // hEffic->Fit( "eff", "EMR" );
    hEffic->Fit( "eff", "","",0.2,15.0 );

    TF1* efficFit = (TF1*)hEffic->GetFunction("eff");

    hEffic->Draw();

    saveName = dir_name + name + ".pdf";
    can->SaveAs(saveName, "PDF");

    can->Destructor();
  }
}


void WeightAndAddCorrected2Ds( TH2D *h_AddedChgUE2D_corr, TH1D *h_DetWt, TH2D *h_ChgUE2D_corr[55], TString dir_name ){

  TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
  can->SetLogz();

  for (int i=0; i<55; ++i) {
    int ptLo = i + 4;
    int binno = i + 1;

    double weight = h_DetWt->GetBinContent(binno);
    h_AddedChgUE2D_corr->Add( h_ChgUE2D_corr[i], weight );
  }

  h_AddedChgUE2D_corr->Draw("COLZ");
  TString name = dir_name + "added" + h_AddedChgUE2D_corr->GetName() + ".pdf";
  // std::cout<<h_AddedChgUE2D_corr->Integral()/AREA<<std::endl;
  can->SaveAs(name, "PDF");

  can->Destructor();

}



void FitEfficHistos( TFile *effic_File, TF1* efficFit[zbins] ){

  TString name;
  TH1D *hEffic[zbins];

  for (int i=1; i<=zbins; ++i) {
    
    TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);

    name = "eff_s_bin_8_10_bbc__"; name += i; name += "_"; name += i; name += "_eta";
    hEffic[i] = (TH1D*)effic_File->Get(name);
    std::cout<<"doot"<<std::endl;
    
    // hEffic->Fit( "eff", "EMR" );
    hEffic[i]->Fit( "eff", "EMR","",0.2,15.0 );

  }
  
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

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           READ IN FILES & HISTOGRAMS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  
  TFile *inFile[nEAbins], *embFile[nEAbins];
  // TFile* outFile = new TFile("out/5gevBins.root", "RECREATE");

  double efficShift = 0.0;
  TFile* outFile = new TFile("out/StopBreaking.root", "RECREATE");
  TString dirName = "plots/StopBreaking";

  TH2D *hUE2D[nEAbins][55];
  TH3D *hUE3D[nEAbins][nEtaBins];
  TH3D *hUE3Dsum[nEAbins];
  TH1D *hLead[nEAbins][nEtaBins];
  TH1D *hLeadSum[nEAbins];
  TH2D *hResponse[nEAbins][nEtaBins];
  TH2D *hResponseSum[nEAbins];
  TH1D *hFakes[nEAbins][nEtaBins];
  TH1D *hFakesSum[nEAbins];
  TH1D *hMisses[nEAbins];

  
  for (int a=0; a<nEAbins; ++a) {
    //dirName += lohi[a];

    name = "../out/UE/pAuHTjetUE_" + lohi[a] + "EA_uncorrected.root";
    // name = "../out/UE/pAuHTjetUE_halfGeVbins_" + lohi[a] + "EA_leadPtUncorrected.root";
    inFile[a] = new TFile(name, "READ");
    name = "hUE3Dsum_" + lohi[a] + "EA";
    hUE3Dsum[a] = new TH3D(name,";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta",xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);
    name = "hLead_" + lohi[a];
    hLeadSum[a] = new TH1D(name,";leading jet p_{T} (GeV)", 55,4.0,59.0);
    name = "hResponseSum_" + lohi[a] + "EA";
    hResponseSum[a] = new TH2D(name,";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.0,59.0, 55,4.0,59.0);
    name = "hFakesSum_" + lohi[a] + "EA";
    hFakesSum[a] = new TH1D(name,";det-level leading jet p_{T} (GeV)", 55,4.0,59.0);

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
    }


    name = "../embedding/out/sim/pAu2015embedding_" + lohi[a] + "EA.root";
    embFile[a] = new TFile(name, "READ");

    hMisses[a] = (TH1D*)embFile[a]->Get("hMisses");
    name = "hMisses_" + lohi[a] + "EA";
    hMisses[a]->SetName(name);
    
    for (int e=0; e<nEtaBins; ++e) {
      name = "hPtResponse_" + emw[e] + "EtaJet";
      hResponse[a][e] = (TH2D*)embFile[a]->Get(name);
      name = "hResponse_" + lohi[a] + "EA_" + emw[e] + "Jet";
      hResponse[a][e]->SetName(name);
      hResponseSum[a]->Add(hResponse[a][e]);

      name = "hFakes_" + emw[e] + "EtaJet";
      hFakes[a][e] = (TH1D*)embFile[a]->Get(name);
      name += "_" + lohi[a];
      hFakes[a][e]->SetName(name);
      hFakesSum[a]->Add(hFakes[a][e]);
    }

  }

  TH1D *hMatched_part[nEAbins], *hMatched_det[nEAbins];//
  for (int a=0; a<nEAbins; ++a) {
    hMatched_part[a] = (TH1D*)hResponseSum[a]->ProjectionX();
    name = "hMatched_part_" + lohi[a] + "EA";
    hMatched_part[a]->SetName(name);

    hMatched_det[a] = (TH1D*)hResponseSum[a]->ProjectionY();
    name = "hMatched_det_" + lohi[a] + "EA";
    hMatched_det[a]->SetName(name);
  }



  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                           FAKE AND MISSED JETS
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TH1D *hFakeProb[nEAbins], *hMissProb[nEAbins], *hPart[nEAbins];//, *hMatchPlusFake[nEAbins], *hMatchPlusMiss[nEAbins];

  for (int a=0; a<nEAbins; ++a) {
    name = "hFakeProb_" + lohi[a] + "EA";
    hFakeProb[a] = new TH1D( name, "Fake Jet Probability; det-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0 );
    name = "hMissProb_" + lohi[a] + "EA";
    hMissProb[a] = new TH1D( name, "Missed Jet Probability; part-level leading jet p_{#mathrm{T}} [GeV]", 55,4.0,59.0);
    
    for (int jp=0; jp<55; ++jp) {
      int binno = jp+1;
      double FakeProb = hFakesSum[a]->GetBinContent(binno) / ( hMatched_det[a]->GetBinContent(binno) + hFakesSum[a]->GetBinContent(binno) );
      if ( isnan(FakeProb) ) { FakeProb = 0.; }
      hFakeProb[a]->SetBinContent( binno, FakeProb );
      double MissProb = hMisses[a]->GetBinContent(binno) / ( hMatched_part[a]->GetBinContent(binno) + hMisses[a]->GetBinContent(binno) );
      if ( isnan(MissProb) ) { MissProb = 0.; }
      hMissProb[a]->SetBinContent( binno, MissProb );
    }

    name = "hPart_" + lohi[a] + "EA";
    title = "Particle-level p_{#mathrm{T}}; part-level leading jet p_{#mathrm{T}} [GeV]";
    hPart[a] = (TH1D*)hMatched_part[a]->Clone(name);
    hPart[a]->SetNameTitle(name,title);
    hPart[a]->Add(hMisses[a]);
    // hPart[a]->Multiply(hMissProb[a]);
    // hPart[a]->Add(hMatched_part[a]);
  }


  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
  //                                             DETECTOR-TO-PARTICLE-LEVEL
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  //  PERFORM 2D TRACKING EFFICINECY CORRECTION
  TH2D* hUE2D_detCorr[nEAbins][55];
  TFile *efficFile = new TFile("src/trackeffic_oct20.root","READ");

  TH1D *hEffic[zbins];
  TF1 *hEfficFit[zbins];  

  
  
  for (int i=1; i<=zbins; ++i) {
    
    TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);

    name = "eff_s_bin_8_10_bbc__"; name += i; name += "_"; name += i; name += "_eta";
    hEffic[i] = (TH1D*)efficFile->Get(name);
    std::cout<<name<<" \t"<<hEffic[i]<<" \t"<<(TH1D*)efficFile->Get(name)<<std::endl;
    
    // hEffic->Fit( "eff", "EMR" );
    hEffic[i]->Fit( "eff", "EMR","",0.2,15.0 );

  }
  
 
  for (int a=0; a<nEAbins; ++a) {
    for (int jp=0; jp<55; ++jp) {
      int binno = jp+1;  int plo=jp+4;  int phi=jp+5;
      hUE3Dsum[a]->GetXaxis()->SetRange(binno,binno);
      hUE2D[a][jp] = (TH2D*)hUE3Dsum[a]->Project3D("ZY");  // UE PT IS ON X-AXIS
      name = "hUE2Dsum_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D[a][jp]->SetName(name);
      hUE2D[a][jp]->Scale(1./hUE2D[a][jp]->Integral());
      hUE2D[a][jp]->Scale(hUE2D[a][jp]->GetEntries()/hLeadSum[a]->Integral(binno,binno));// NORMALIZE TO NJETS
      // cout<<hUE2D[a][jp]->Integral()<<"\t"<<hUE2D[a][jp]->GetEntries()<<endl;
      hUE3Dsum[a]->GetXaxis()->SetRange(1,-1);

      name = "hUE2Dcorr_" + lohi[a] + "_"; name+=plo; name+="_"; name+=phi; name+="GeV_det";
      hUE2D_detCorr[a][jp] = new TH2D(name,";chg. UE part. p_{T} (GeV);chg. UE part. #eta",ybins,ybinEdge,zbins,zbinEdge);
      name = dirName + "/";
    }
    // TrackingEfficiencyTEST( hUE2D[a], hUE2D_detCorr[a], efficFile, lohi[a], name, efficShift );
  }

  return 0;
}
