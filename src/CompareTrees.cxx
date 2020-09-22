// Veronica Verkest
// September 14, 2020

#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TClonesArray.h"
#include "TLatex.h"
#include "TMathText.h"
#include "TProfile.h"

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

void PrintInfo(int event, int max){
  double percent = (double) event*100; percent /= (double) max;
  cout<<"event number "<<event<<"/"<<max<<" \t"<<setprecision(3)<<percent<<"%"<<endl;
}

int main( int argc, char* argv[] ) { // specify starting event number

  int inVal;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc == 2 ) {
    //inVal = atoi(arguments[1].c_str());
    inVal = strtol(argv[1], NULL, 10);
  }
  else if ( argc == 1 ) { inVal = 0; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  int startEvent = inVal*50000;
  
  TString outFileName = "CompareTrees/CompareTrees_"; outFileName += startEvent; outFileName += ".root";

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  const int nFiles = 2;

  TString inFileName[nFiles] = { "out/UE/pAuHTjetUE_HILO_diffPt_sept.root", "out/UE/pAuHTjetUE_30cmVzCut_3cmVzDiff_HILO__binShift.root" };

  TFile *inFile[nFiles];
  TTree *jetTree[nFiles+1];
  double nEntries[nFiles+1]; 

  int RunID[nFiles+1], EventID[nFiles+1], nTowers[nFiles+1], nPrimary[nFiles+1], nGlobal[nFiles+1], nVertices[nFiles+1], refMult[nFiles+1],
    gRefMult[nFiles+1], nUEpart_chg[nFiles+1], nUEpart_neu[nFiles+1], nHTtrig[nFiles+1], nJetsAbove5[nFiles+1];
  double Vz[nFiles+1], BbcAdcSumEast[nFiles+1], leadPt[nFiles+1], leadEta[nFiles+1], leadPhi[nFiles+1], chgEastRho[nFiles+1], chgMidRho[nFiles+1],
    chgWestRho[nFiles+1], neuEastRho[nFiles+1], neuMidRho[nFiles+1], neuWestRho[nFiles+1], leadArea[nFiles+1], eastRho[nFiles+1], midRho[nFiles+1],
    westRho[nFiles+1], leadPtCorrected[nFiles+1], chgEastRho_te[nFiles+1], chgMidRho_te[nFiles+1], chgWestRho_te[nFiles+1], rho_te[nFiles+1],
    rho[nFiles+1], dPhiTrigLead[nFiles+1], dRTrigLead[nFiles+1];

  TH1D *hPt = new TH1D("hPt",";leading jet p_{T} GeV",60,0,60);
  
  for (int i=0; i<nFiles; ++i) {
    
    inFile[i] = new TFile( inFileName[i], "READ" );
    jetTree[i] = (TTree*)inFile[i]->Get( "HTjetTree" );

    jetTree[i]->SetBranchAddress( "RunID", &RunID[i] );
    jetTree[i]->SetBranchAddress( "EventID", &EventID[i] );
    jetTree[i]->SetBranchAddress( "nTowers", &nTowers[i] );
    jetTree[i]->SetBranchAddress( "nPrimary", &nPrimary[i] );
    jetTree[i]->SetBranchAddress( "nGlobal", &nGlobal[i] );
    jetTree[i]->SetBranchAddress( "nVertices", &nVertices[i] );
    jetTree[i]->SetBranchAddress( "refMult", &refMult[i] );
    jetTree[i]->SetBranchAddress( "gRefMult", &gRefMult[i] );
    jetTree[i]->SetBranchAddress( "Vz", &Vz[i] );
    jetTree[i]->SetBranchAddress( "leadPt", &leadPt[i] );
    jetTree[i]->SetBranchAddress( "BbcAdcSumEast", &BbcAdcSumEast[i] );
    jetTree[i]->SetBranchAddress( "leadEta", &leadEta[i] );
    jetTree[i]->SetBranchAddress( "leadPhi", &leadPhi[i] );
    jetTree[i]->SetBranchAddress( "chgEastRho", &chgEastRho[i] );
    jetTree[i]->SetBranchAddress( "chgMidRho", &chgMidRho[i] );
    jetTree[i]->SetBranchAddress( "chgWestRho", &chgWestRho[i] );
    jetTree[i]->SetBranchAddress( "chgEastRho_te", &chgEastRho_te[i] );
    jetTree[i]->SetBranchAddress( "chgMidRho_te", &chgMidRho_te[i] );
    jetTree[i]->SetBranchAddress( "chgWestRho_te", &chgWestRho_te[i] );
    jetTree[i]->SetBranchAddress( "neuEastRho", &neuEastRho[i] );
    jetTree[i]->SetBranchAddress( "neuMidRho", &neuMidRho[i] );
    jetTree[i]->SetBranchAddress( "neuWestRho", &neuWestRho[i] );
    jetTree[i]->SetBranchAddress( "leadArea", &leadArea[i] );
    jetTree[i]->SetBranchAddress( "leadPtCorrected", &leadPtCorrected[i] );
    jetTree[i]->SetBranchAddress( "nHTtrig", &nHTtrig[i] );
    jetTree[i]->SetBranchAddress( "nJetsAbove5", &nJetsAbove5[i] );
    jetTree[i]->SetBranchAddress( "dPhiTrigLead", &dPhiTrigLead[i] );
    jetTree[i]->SetBranchAddress( "dRTrigLead", &dRTrigLead[i] );

    nEntries[i] = jetTree[i]->GetEntries();

  }

  TFile *outFile = new TFile(outFileName,"RECREATE");

  jetTree[2] = new TTree("newTree","new_tree");
  jetTree[2]->Branch( "RunID", &RunID[2] );
  jetTree[2]->Branch( "EventID", &EventID[2] );
  jetTree[2]->Branch( "nTowers", &nTowers[2] );
  jetTree[2]->Branch( "nPrimary", &nPrimary[2] );
  jetTree[2]->Branch( "nGlobal", &nGlobal[2] );
  jetTree[2]->Branch( "nVertices", &nVertices[2] );
  jetTree[2]->Branch( "refMult", &refMult[2] );
  jetTree[2]->Branch( "gRefMult", &gRefMult[2] );
  jetTree[2]->Branch( "Vz", &Vz[2] );
  jetTree[2]->Branch( "leadPt", &leadPt[2] );
  jetTree[2]->Branch( "BbcAdcSumEast", &BbcAdcSumEast[2] );
  jetTree[2]->Branch( "leadEta", &leadEta[2] );
  jetTree[2]->Branch( "leadPhi", &leadPhi[2] );
  jetTree[2]->Branch( "leadArea", &leadArea[2] );
  jetTree[2]->Branch( "leadPtCorrected", &leadPtCorrected[2] );
  jetTree[2]->Branch( "nHTtrig", &nHTtrig[2] );
  jetTree[2]->Branch( "nJetsAbove5", &nJetsAbove5[2] );
  jetTree[2]->Branch( "dPhiTrigLead", &dPhiTrigLead[2] );
  jetTree[2]->Branch( "dRTrigLead", &dRTrigLead[2] );

  int nMissed = 0;
  int endEvent = startEvent + 50000;
  int max_entry = jetTree[1]->GetEntries();
  // int max_entry = 50000;
  for ( int j=startEvent; j<endEvent; ++j ) {

    if (j>=max_entry) { continue; }
    
    jetTree[1]->GetEvent(j);
    if (j%1000==0) { PrintInfo(j,max_entry); }

    while ( leadPtCorrected[1]<10.0 || leadPtCorrected[1]>30.0 ) {
      if (j>=max_entry) { continue; }
      j+=1;
      jetTree[1]->GetEvent(j);
      if (j%1000==0) { PrintInfo(j,max_entry); }
    } // find event in new file with pt in range
    if (j>=max_entry) { continue; }

    TString selection = "RunID=="; selection += RunID[1]; selection += " && EventID=="; selection += EventID[1];
    // cout<<drawString<<" \t"<<selection<<endl;
    int hasEvent = jetTree[0]->Scan("",selection);

    if (hasEvent==0) {
      jetTree[1]->GetEvent(j);
      nMissed += 1;
      RunID[2] = RunID[1];
      EventID[2] = EventID[1];
      nTowers[2] = nTowers[1];
      nPrimary[2] = nPrimary[1];
      nGlobal[2] = nGlobal[1];
      nVertices[2] = nVertices[1];
      refMult[2] = refMult[1];
      gRefMult[2] = gRefMult[1];
      Vz[2] = Vz[1];
      leadPt[2] = leadPt[1];
      BbcAdcSumEast[2] = BbcAdcSumEast[1];
      leadEta[2] = leadEta[1];
      leadPhi[2] = leadPhi[1];
      leadArea[2] = leadArea[1];
      leadPtCorrected[2] = leadPtCorrected[1];
      nHTtrig[2] = nHTtrig[1];
      nJetsAbove5[2] = nJetsAbove5[1];
      dPhiTrigLead[2] = dPhiTrigLead[1];
      dRTrigLead[2] = dRTrigLead[1];
      jetTree[2]->Fill();
    }
    
  }

  cout<<nMissed<<endl<<endl;

  jetTree[2]->Write();
  outFile->Write();
  outFile->Close();
  
  return 0;
}
