// Veronica Verkest - July 15, 2020

#include "params.hh"
#include "funcs.hh"
#include "TStarJetPicoDefinitions.h"

using namespace fastjet;
using namespace std;
using namespace Analysis;

int main (int argc, const char ** argv) {
  TH1::SetDefaultSumw2( );    TH2::SetDefaultSumw2( );    TH3::SetDefaultSumw2( );

  //! Command line arguments:
  //! [0]: output filename
  //! [1]: flag: require match? Options: "nomatch" (match = 0), or "match" (match = 1).
  //! [2]: input data

  // Defaults
  string outFileName = "out/test.root"; //output file
  string chainList = "lists/simlist.txt"; // input file: can be .root, .txt, .list
  bool match = 1; //match = 0 => no match between Pythia & Pythia+Geant events.
  
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
  case 1: // Default case
    __OUT("Using Default Settings");
    break;
  case 4: { // Custom case
    __OUT("Using Custom Settings");
    std::vector<std::string> arguments( argv+1, argv+argc );

    // Set non-default values    
    outFileName       = arguments[0];
            
    if (arguments[1] == "matched" || arguments[1] == "match") {match = 1;}
    else if (arguments[4] == "unmatched") {match = 0;}
    else {cerr << "Not a valid flag!" << endl; exit(1);}
    chainList = arguments[2];
      
    //Printing settings:
    cout << "Outputting to: " << (outFileName).c_str() << "\nSettings:\n " << arguments[2] << " jets;\n input file: " << chainList << "\n";
    break;
  }
  default: { // Error: invalid custom settings
    __ERR("Invalid number of command line arguments");
    return -1;
    break;
  }
  }


    //  Initialize readers and provide chains
  TStarJetPicoReader partReader;
  TChain* partChain = new TChain( "JetTreeMc" ); // PURE PYTHIA (particle)
  TStarJetPicoReader detReader;
  TChain* detChain = new TChain( "JetTree" ); // PYTHIA + GEANT + 2015pAuZB
  
  // If its a recognized file type, build the chain, if it's not recognized, exit
  if ( Analysis::HasEnding( chainList.c_str(), ".root" ) ) { partChain->Add( chainList.c_str()); detChain->Add( chainList.c_str());}
  else if ( Analysis::HasEnding( chainList.c_str(), ".txt"  ) || Analysis::HasEnding( chainList.c_str(), ".list" ) )  { partChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str()); detChain = TStarJetPicoUtils::BuildChainFromFileList(chainList.c_str());}
  else { __ERR("data file is not recognized type: .root or .txt only.") return -1; }
  
  TFile *fout = new TFile( ( outFileName ).c_str() ,"RECREATE");
  fout->cd();

  TTree *eventTree = new TTree( "eventTree", "event_Tree" );

  int RunID, EventID, refMult;
  double Vz, BbcAdcSumEast, p_leadPt, p_leadEta, p_leadPhi, d_leadPt, d_leadEta, d_leadPhi, mc_weight;

  eventTree->Branch( "RunID", &RunID );
  eventTree->Branch( "EventID", &EventID );
  eventTree->Branch( "refMult", &refMult );
  eventTree->Branch( "Vz", &Vz );
  eventTree->Branch( "BbcAdcSumEast", &BbcAdcSumEast );
  eventTree->Branch( "p_leadPt", &p_leadPt );
  eventTree->Branch( "p_leadEta", &p_leadEta );
  eventTree->Branch( "p_leadPhi", &p_leadPhi );
  eventTree->Branch( "d_leadPt", &d_leadPt );
  eventTree->Branch( "d_leadEta", &d_leadEta );
  eventTree->Branch( "d_leadPhi", &d_leadPhi );
  eventTree->Branch( "mc_weight", &mc_weight );
  
  TString partFilename, detFilename;  //define relevant data structures
  TStarJetPicoEventHeader* p_header; TStarJetPicoEventHeader* d_header;
  TStarJetPicoEvent* p_event; TStarJetPicoEvent* d_event;
  TStarJetVectorContainer<TStarJetVector>* p_container; TStarJetVectorContainer<TStarJetVector>* d_container;
  TStarJetVector *p_sv = 0; TStarJetVector *d_sv = 0;
  
  JetDefinition jet_def( antikt_algorithm, R );
  Selector jetSelector = SelectorAbsRapMax( maxJetEta ) && SelectorPtMin( minJetPt );

  vector<PseudoJet> p_Particles, d_Particles, p_Jets, d_Jets, p_leadMatches, p_matches, d_matches, fakeJets, missedJets;
  PseudoJet p_leadJet, d_leadJet;

  TH2D *hPtResponse = new TH2D("hPtResponse",";part-level leading jet p_{T} (GeV);det-level leading jet p_{T} (GeV)",55,4.5,59.5, 55,4.5,59.5);
  TH1D *hMisses = new TH1D( "hMisses","Missses;missing part-level leading jet p_{T} (GeV)",55,4.5,59.5);
  TH1D *hFakes = new TH1D( "hFakes","Fakes;fake det-level leading jet p_{T} (GeV)",55,4.5,59.5);
  TH3D *hTrigEtEtaPhi = new TH3D( "hTrigEtEtaPhi",";Trigger E_{T} (GeV);Trigger #eta;Trigger #phi",30,0,30, 40,-1.0,1.0, 120,0,2*pi);

  InitReader( partReader, partChain, -1, 9999, 1000, 1000, -1, 9999, 100, -1, -1, 1000, 9999, "lists/dummy_tower_list.txt" );
  InitReader( detReader, detChain, -1, detVertexZCut, detMaxEventPtCut, detMaxEventEtCut, detMinEventEtCut, detVertexZDiffCut,
	      detDCACut, detMinNFitPointsCut, detFitOverMaxPointsCut, detMaxPtCut, detMaxEtCut, det_badTowers  );

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( partReader.NextEvent() ) {  // loop through all part-level events

    partReader.PrintStatus(10); detReader.PrintStatus(10);     // Print out reader status every 10 seconds

    EventID = partReader.GetNOfCurrentEvent();
    if ( !(detReader.ReadEvent(EventID)) ) {continue;}  // see if corresponding det-level event exists; if not, skip event
    detReader.ReadEvent( EventID );

    p_Particles.clear();    d_Particles.clear();    p_Jets.clear();    d_Jets.clear();    p_leadMatches.clear();  // clear vectors
    p_matches.clear();    d_matches.clear();    fakeJets.clear();    missedJets.clear();    
    
    d_event = detReader.GetEvent();    d_header = d_event->GetHeader();    d_container = detReader.GetOutputContainer();
    p_event = partReader.GetEvent();   p_header = p_event->GetHeader();    p_container = partReader.GetOutputContainer();

    if ( d_header->GetRunId() >= 16142059 && d_header->GetRunId() <= 16149001) { continue; } // TEMPORARILY SKIPPING THESE RUNS
    if ( d_header->GetRunId() == 16135031 || d_header->GetRunId() == 16135032) { continue; }
    if (!(d_header->HasTriggerId(500205) || d_header->HasTriggerId(500215))) {continue;}   //  ONLY SELECT HT TRIGGER EVENTS
    if ( d_header->GetBbcAdcSumEast()>64000 || d_header->GetBbcAdcSumEast()<3559.12 ) { continue; } //  neglect 0-10% and 90-100% event activity
    if ( !(d_header->GetRunId()==p_header->GetRunId() && detReader.GetNOfCurrentEvent()==partReader.GetNOfCurrentEvent()) ) {
      cerr<<"UNMATCHED EVENT OR RUN";        // PART AND DET EVENTS MUST MATCH
      continue;
    }
    
    GatherParticles( d_container, d_Particles );     GatherParticles( p_container, p_Particles );

    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
    
    ClusterSequenceArea d_jetCluster( d_Particles, jet_def, area_def);
    d_Jets = sorted_by_pt( jetSelector( d_jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS WITH pT >= 5 GeV
    ClusterSequenceArea p_jetCluster( p_Particles, jet_def, area_def); 
    p_Jets = sorted_by_pt( jetSelector( p_jetCluster.inclusive_jets() ) );

    if ( d_Jets.size()==0 && p_Jets.size()==0 ) { continue; }

    mc_weight = LookupRun15Xsec( partFilename );

    if ( d_Jets.size()==0 ) {
      hMisses->Fill( p_Jets[0].pt(), mc_weight );
      continue;
    } // missed jet
    if ( p_Jets.size()==0 ) {
      hFakes->Fill( d_Jets[0].pt(), mc_weight );
      continue;
    } // fake jet
    // only events with 1+ det and 1+ part jet pass this line

    partFilename =  partReader.GetInputChain()->GetCurrentFile()->GetName();
    detFilename =  detReader.GetInputChain()->GetCurrentFile()->GetName();
    if ( p_Jets.size()!=0 && DiscardpAuEmbedEvent( partFilename, p_Jets, d_Jets ) ) { continue; }
    
    d_leadJet = d_Jets[0];

    TList *SelectedTowers = detReader.GetListOfSelectedTowers();
    TStarJetPicoTriggerInfo *trig;
    TStarJetPicoTower *tow;
    int nTowers = CountTowers( SelectedTowers );

    vector<int> trigTowers;
    for ( int i=0; i<d_event->GetTrigObjs()->GetEntries(); ++i ) {
      trig = (TStarJetPicoTriggerInfo *)d_event->GetTrigObj(i);
      if ( trig->isBHT2() && UseTriggerTower( trig->GetId()) ) { trigTowers.push_back( trig->GetId() ); }
    }
    sort(trigTowers.begin(), trigTowers.end());

    PseudoJet trigTowerPJ, towPJ;
    double deltaPhi, deltaR;
    int nTowersMatched = 0;
    for (int i=0; i<nTowers; ++i){				// loop throught selected towers in event
      tow = (TStarJetPicoTower*)SelectedTowers->At(i);
      if ( tow->GetEt()>=5.4 && count(trigTowers.begin(), trigTowers.end(), tow->GetId())) { // min 5.4 GeV tower and must be in list of HT towers
	
    	towPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
    	deltaPhi = fabs( d_leadJet.delta_phi_to( towPJ ) );
    	deltaR = d_leadJet.delta_R( towPJ );
    	double trigTowEt = 0;
    	  if ( deltaR<=R || fabs(deltaPhi)>=(pi-R) ) {  // require trigger
    	    if ( nTowersMatched>0 && (tow->GetEt()<trigTowEt) ) { continue; }   // more than 1 trigger tower
    	    else {							// first trigger tower
    	      trigTowEt = tow->GetEt();
    	      trigTowerPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
    	      nTowersMatched += 1;
    	    }
    	  }
	  
      }
    }
    if (nTowersMatched==0) { continue; } // require HT trigger to be within lead det jet

    Selector LeadJetMatcher = SelectorCircle( R );  // Match triggered lead det-level jet to part-level jet
    LeadJetMatcher.set_reference( d_leadJet );
    p_leadMatches = sorted_by_pt( LeadJetMatcher( p_Jets ));
    if ( p_leadMatches.size()==0 ) {
      hFakes->Fill( d_leadJet.pt(), mc_weight );
      continue;
    } // fake jet
    // only events with a det lead jet containing a trigger and a matched part jet pass this line
    
    p_leadJet = p_leadMatches[0];
    
    RunID = d_header->GetRunId();
    EventID = detReader.GetNOfCurrentEvent();
    refMult = d_header->GetReferenceMultiplicity();
    Vz = d_header->GetPrimaryVertexZ();
    BbcAdcSumEast = d_header->GetBbcAdcSumEast();
    p_leadPt = p_leadJet.pt();
    p_leadEta = p_leadJet.eta();
    p_leadPhi = p_leadJet.phi();
    d_leadPt = d_leadJet.pt();
    d_leadEta = d_leadJet.eta();
    d_leadPhi = d_leadJet.phi();
    
    hPtResponse->Fill( p_leadPt, d_leadPt, mc_weight );
    hTrigEtEtaPhi->Fill( trigTowerPJ.e(), trigTowerPJ.eta(), trigTowerPJ.phi(), mc_weight );
    
    eventTree->Fill();
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
  hPtResponse->Write();
  hTrigEtEtaPhi->Write();
  hMisses->Write();
  hFakes->Write();

  eventTree->Write();
  
  fout->Close();
  
  return 0;
}//main
