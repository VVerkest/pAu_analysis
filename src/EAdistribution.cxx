//  EAdistribution.cxx
//  Veronica Verkest		February 11, 2020

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/MB/pAu_2015_200_MB*.root"; outFile="src/EAdistribution.root"; number_of_events=4; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH1D *hEAdist = new TH1D( "hEAdist", "BBC Inner East Sum;BBCEsum", 1000,0,100000 );
  
  TTree *EAtree = new TTree( "EAtree", "EA_Tree" );
  
  //  Tree variables
  int RunID, EventID;
  double Vz, BbcAdcSumEast;

  EAtree->Branch( "RunID", &RunID );				EAtree->Branch( "EventID", &EventID );				EAtree->Branch( "Vz", &Vz );
  EAtree->Branch( "BbcAdcSumEast", &BbcAdcSumEast );

  string trigOpt = "MB";
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in MB: 56,855,250
  InitReader( Reader, Chain, numEvents );
  
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {
    
    Reader.PrintStatus(10);

    event = Reader.GetEvent();
    header = event->GetHeader();
    
    Vz = header->GetPrimaryVertexZ();
    //if ( UseEvent( header, event, vzCut, Vz, trigOpt ) == false ) { cout<<"skip"<<endl; continue; }   //  Skip events based on: Run#, vz cut, BBCSumEast;    only accept MB events
    
    cout<<"0"<<endl;
    RunID = header->GetRunId();
    EventID = Reader.GetNOfCurrentEvent();
    BbcAdcSumEast = header->GetBbcAdcSumEast();

    EAtree->Fill();
    hEAdist->Fill( BbcAdcSumEast );

  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  TFile *EAdistFile = new TFile( outFile.c_str() ,"RECREATE");
  hEAdist->Write();
  EAtree->Write();
  
  return 0;
}
