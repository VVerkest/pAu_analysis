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

  const double pi = 3.14159;
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/HT/pAuJets_trigger.root"; number_of_events=-1; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetVectorContainer<TStarJetVector> * container;
  TStarJetPicoTriggerInfo *trig;
  TStarJetPicoTower *tow;
  TStarJetPicoPrimaryTrack *track;
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  int nTowers, nEvents, nPrimary;
  double Vz, leadPt, leadEta, leadPhi, eastRho, midRho, westRho, rho, ptSum, trigTowEt, trigTowEta, trigTowPhi;

  TH3D *hTrigEtEtaPhi = new TH3D("hTrigEtEtaPhi","Trigger E_{T} vs. #eta vs. #phi;E_{T} (GeV);#eta;#phi", 120,0,60, 40,-1.0,1.0, 120,-pi,pi);
  TH2D *hTrigEtId = new TH2D("hTrigEtId","Trigger E_{T} vs. ID", 4800,0,4800, 120,0,60 );
  TH3D *hTrackPtEtaPhi = new TH3D("hTrackPtEtaPhi","Track p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 120,0,60, 40,-1.0,1.0, 120,-pi,pi);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;E_{T} (GeV);#eta;#phi", 120,0,60, 40,-1.0,1.0, 120,-pi,pi);
  TH1D *hDeltaPhiToJet = new TH1D("hDeltaPhiToJet","|#phi_{lead} - #phi_{trig}|;#Delta#phi", 120,0,2*pi );
  TH1D *hDeltaRToJet = new TH1D("hDeltaRToJet","lead jet-trigger #Delta R;#Delta R", 120,0,2*pi);
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(15.0);
  Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;

  PseudoJet leadJet;
  vector<PseudoJet> rawParticles, rawJets;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    rawParticles.clear();    rawJets.clear();    //  CLEAR VECTORS
    
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

      GatherParticles( container, rawParticles );
      GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
      AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
      ClusterSequenceArea jetCluster( rawParticles, jet_def, area_def); 
      
      leadJetSelector = leadPtMinSelector && jetEtaSelector;
      rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >=15 GeV


      if (rawJets.size()>0) {

	trigTowEt = 0.0;
	
	
	TList *SelectedTowers = Reader.GetListOfSelectedTowers(); //  TRACK, TOWER, TRIGGER
	nTowers = CountTowers( SelectedTowers );
	for (int i=0; i<nTowers; ++i){
	  tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	  hTowEtEtaPhi->Fill( tow->GetEt(), tow->GetEta(), tow->GetPhi() );
	}
	for ( int i=0; i<event->GetTrigObjs()->GetEntries(); ++i ) {
	  
	  trig = (TStarJetPicoTriggerInfo *)event->GetTrigObj(i);
	  int tempID = trig->GetId();
	  for (int i=0; i<nTowers; ++i){
	    tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	    if ( tow->GetId()==tempID ){
	      hTrigEtEtaPhi->Fill( tow->GetEt(), tow->GetEta(), tow->GetPhi() );
	      hTrigEtId->Fill( tempID, tow->GetEt() );
	      if ( tow->GetEt() > trigTowEt ) {   //
		trigTowEt = tow->GetEt();
		trigTowEta = tow->GetEta();
		trigTowPhi = tow->GetPhi() + pi; // get trigger tower phi in range [0,2pi]		
	      }
	    }
	  }
	}
	nPrimary = header->GetNOfPrimaryTracks();
	for (int i=0; i<nPrimary; ++i) {
	  track = (TStarJetPicoPrimaryTrack*) event->GetPrimaryTrack(i);
	  hTrackPtEtaPhi->Fill( track->GetPt(), track->GetEta(), track->GetPhi() );
	}    //  TRACK, TOWER, TRIGGER


	leadJet = rawJets[0];
	double dphi, deta, dR;
	dphi = fabs( trigTowPhi - leadJet.phi() );
	deta = fabs( trigTowEta - leadJet.eta() );
	dR = sqrt( (dphi*dphi) + (deta*deta) );
	hDeltaPhiToJet->Fill( dphi );
	hDeltaRToJet->Fill( dR );

	
      }
      
    }
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hTrigEtEtaPhi->Write();
  hTowEtEtaPhi->Write();
  hTrackPtEtaPhi->Write();
  hTrigEtId->Write();
  hDeltaPhiToJet->Write();
  hDeltaRToJet->Write();
  
  pAuFile->Close();

  return 0;
}
