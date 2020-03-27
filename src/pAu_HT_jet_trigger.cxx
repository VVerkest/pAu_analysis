//  pAu_HT_jet_trigger.cxx  (previously pAu_HT_jets.cxx)
//  Veronica Verkest		July 20, 2019

//HT2*BBCMB : 500205, 500215		JP2 : 500401, 500411
//BBCMB : 500008, 500018			VPDMB :  500904

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/HT/pAuJets_trigger.root"; number_of_events=-1; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  TStarJetPicoTriggerInfo *trig;
  
  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  int nTowers, nEvents, nPrimary;
  double Vz, leadPt, leadEta, leadPhi, eastRho, midRho, westRho, rho, ptSum;
  const double pi = 3.14159;
  
  TH3D *hTrigEtEtaPhi = new TH3D("hTrigEtEtaPhi","Trigger E_{T} vs. #eta vs. #phi;E_{T} (GeV);#eta;#phi", 120,0,30, 40,-1.0,1.0, 120,-pi,pi);
  TH3D *hTrackPtEtaPhi = new TH3D("hTrackPtEtaPhi","Track p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 120,0,30, 40,-1.0,1.0, 120,-pi,pi);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Trigger E_{T} vs. #eta vs. #phi;E_{T} (GeV);#eta;#phi", 120,0,30, 40,-1.0,1.0, 120,-pi,pi);
  TH2D *hTrigEtId = new TH2D("hTrigEtId","Trigger E_{T} vs. ID", 4800,0,4800, 120,0,30 );

  TStarJetPicoTriggerInfo *trig;
  TStarJetPicoTower *tow;
  TStarJetPicoPrimaryTrack *track;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    Vz = header->GetPrimaryVertexZ();
    
    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }
    if ( abs(Vz) > vzCut ) { continue; }
    if (!(header->HasTriggerId(500205) || header->HasTriggerId(500215))) {continue;}   //  ONLY SELECT HT TRIGGER EVENTS
    if ( header->GetBbcAdcSumEast() > 64000 ) { continue; }
    if ( header->GetBbcAdcSumEast() < 3559.12 ) { continue; }     //  neglect 0-10% event activity
    else {

      TList *SelectedTowers = Reader.GetListOfSelectedTowers();
      nTowers = CountTowers( SelectedTowers );
      for (int i=0; i<nTowers; ++i){
	tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	hTowEtEtaPhi->Fill( tow->GetEt(), tow->GetEta(), tow->GetPhi() );
      }
      
      for ( int i=0; i<event->GetTrigObjs()->GetEntries(); ++i ) {
	trig = (TStarJetPicoTriggerInfo *)event->GetTrigObj(i);
	int trigTowId = trig->GetId();
	
	for (int i=0; i<nTowers; ++i){
	  tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	  if ( tow->GetId()==trigTowId ){
	    hTrigEtEtaPhi->Fill( tow->GetEt(), trig->GetEta, trig->GetPhi() );
	  }
	}
      }

      nPrimary = header->GetNOfPrimaryTracks();
      for (int i=0; i<nPrimary; ++i) {
	track = (TStarJetPicoPrimaryTrack*) event->GetPrimaryTrack(i);
	hTrackPtEtaPhi->Fill( track->GetPt(), track->GetEta(), track->GetPhi() );
      }
      
    }
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hTrigEtEtaPhi->Write();
  hTrigEtId->Write();
  
  pAuFile->Close();

  return 0;
}
