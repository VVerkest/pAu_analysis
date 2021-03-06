//  pAuHTdijetUE.cxx
//  Veronica Verkest		November 12, 2019

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTdijetUE.root"; number_of_events=100; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTdijetTree = new TTree( "HTdijetTree", "HT_DijetTree" );
  string efficFile = "src/trackeffic.root";
  
  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, eval, pval, nBGpart_chg, nBGpart_neu, nHTtrig;
  double Vz, BbcAdcSumEast, BbcAdcSumEastOuter, BbcAdcSumWest, BbcAdcSumWestOuter, leadPt, recoPt, leadEta, recoEta, leadPhi, recoPhi,
    chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea, recoArea, dRTrigLead, dRTrigReco;

  HTdijetTree->Branch( "RunID", &RunID );
  HTdijetTree->Branch( "EventID", &EventID );
  HTdijetTree->Branch( "nTowers", &nTowers );
  HTdijetTree->Branch( "nPrimary", &nPrimary );
  HTdijetTree->Branch( "nGlobal", &nGlobal );
  HTdijetTree->Branch( "nVertices", &nVertices );
  HTdijetTree->Branch( "refMult", &refMult );
  HTdijetTree->Branch( "gRefMult", &gRefMult );
  HTdijetTree->Branch( "Vz", &Vz );
  HTdijetTree->Branch( "leadPt", &leadPt );
  HTdijetTree->Branch( "leadEta", &leadEta );
  HTdijetTree->Branch( "BbcAdcSumEast", &BbcAdcSumEast );
  HTdijetTree->Branch( "BbcAdcSumEastOuter", &BbcAdcSumEastOuter );
  HTdijetTree->Branch( "BbcAdcSumWest", &BbcAdcSumWest );
  HTdijetTree->Branch( "BbcAdcSumWestOuter", &BbcAdcSumWestOuter );
  HTdijetTree->Branch( "leadPhi", &leadPhi );
  HTdijetTree->Branch( "chgEastRho", &chgEastRho );
  HTdijetTree->Branch( "chgMidRho", &chgMidRho );
  HTdijetTree->Branch( "chgWestRho", &chgWestRho );
  HTdijetTree->Branch( "neuEastRho", &neuEastRho );
  HTdijetTree->Branch( "neuMidRho", &neuMidRho );
  HTdijetTree->Branch( "neuWestRho", &neuWestRho );
  HTdijetTree->Branch( "leadArea", &leadArea );
  HTdijetTree->Branch( "recoArea", &recoArea );
  HTdijetTree->Branch( "recoPt", &recoPt );
  HTdijetTree->Branch( "recoEta", &recoEta );
  HTdijetTree->Branch( "recoPhi", &recoPhi );
  HTdijetTree->Branch( "nBGpart_chg", &nBGpart_chg );
  HTdijetTree->Branch( "nBGpart_neu", &nBGpart_neu );
  HTdijetTree->Branch( "nHTtrig", &nHTtrig );
  
  TH3D *hChgBgPtEtaPhi = new TH3D( "hChgBgPtEtaPhi", "Charged Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-1.0,1.0, 120,0.0,2*pi );
  TH3D *hNeuBgPtEtaPhi = new TH3D( "hNeuBgPtEtaPhi", "Neutral Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-1.0,1.0, 120,0.0,2*pi );

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
  Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
  Selector ptMaxSelector = SelectorPtMax( 30.0 );

  PseudoJet leadJet, recoJet, trigTowerPJ, towPJ;
  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles, BGparticles, recoCandidates;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  double deltaPhi, deltaR;       double trigTowEta, trigTowPhi;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);
    rawParticles.clear();            chgParticles.clear();            neuParticles.clear();
    rawJets.clear();            recoCandidates.clear();             //  CLEAR VECTORS

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }
    if ( abs(Vz) > vzCut ) { continue; }
    if (!(header->HasTriggerId(500205) || header->HasTriggerId(500215))) {continue;}   //  ONLY SELECT HT TRIGGER EVENTS
    if ( header->GetBbcAdcSumEast() > 64000 ) { continue; }
    if ( header->GetBbcAdcSumEast() < 3559.12 ) { continue; }     //  neglect 0-10% event activity

    TList *SelectedTowers = Reader.GetListOfSelectedTowers();
    nTowers = CountTowers( SelectedTowers );

    Vz = header->GetPrimaryVertexZ();
    //if ( UseHTevent( header, event, vzCut, Vz ) == false ) { continue; } // Skip events based on: Run#, vz cut, BBCSumE; only accept HT events

    GatherParticles( container, rawParticles );

    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
    ClusterSequenceArea jetCluster( rawParticles, jet_def, area_def); 

    leadJetSelector = leadPtMinSelector && ptMaxSelector && jetEtaSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN 10-30 GeV RANGE
    
    if ( rawJets.size()>0 ) {
      leadJet = rawJets[0]; 

            
      //   REQUIRE RECOIL JET TO HAVE AT LEAST HALF OF LEAD PT!
      Selector recoPtRangeSelector = SelectorPtRange( leadJet.pt()/2, ptHi[pval] );          //  JET pT RANGE    { 10-15, 15-20, 20-30 }
      Selector etaRangeSelector = SelectorEtaRange( etaLo[eval], etaHi[eval] );          //  JET eta RANGE

      // Selector recoJetSelector = recoPtRangeSelector && etaRangeSelector && jetEtaSelector;
      Selector recoJetSelector = recoPtRangeSelector && jetEtaSelector;
    
      recoCandidates = sorted_by_pt( recoJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    
      bool hasReco = false;
    
      if ( recoCandidates.size()>0 ) {

	for ( int i=0; i<recoCandidates.size(); ++i ) {
	  double deltaPhi_abs = fabs( recoCandidates[i].phi() - leadJet.phi() );
	  if ( fabs( deltaPhi_abs - pi ) <= R ) {
	    recoJet = recoCandidates[i];
	    hasReco = true;
	  }
	  if ( hasReco == true ) { continue; }  // exit loop with highest-pt jet in recoil range
	}

      }
      else { continue; }

      if ( hasReco == true ) {
 
	int trigTowId;
	TStarJetPicoTriggerInfo *trig;
	TStarJetPicoTower *tow, *triggerTower;
	double trigTowEt = 0.0;
	std::vector<int> trigTowers;
	for ( int i=0; i<event->GetTrigObjs()->GetEntries(); ++i ) {
	  trig = (TStarJetPicoTriggerInfo *)event->GetTrigObj(i);
	  if ( trig->isBHT2() && UseTriggerTower( trig->GetId()) ) { trigTowers.push_back( trig->GetId() ); }
	}
	sort(trigTowers.begin(), trigTowers.end());

	int nmatched = 0;
	for (int i=0; i<nTowers; ++i){				// loop throught selected towers in event
	  tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	  if ( tow->GetEt()>=5.4 && count(trigTowers.begin(), trigTowers.end(), tow->GetId())) { // min 5.4 GeV tower and must be in list of HT towers

	    towPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
	    dRTrigLead = leadJet.delta_R( towPJ );
	    dRTrigReco = recoJet.delta_R( towPJ );


	    if ( dRTrigLead<=R || dRTrigReco<=R ) { // require trigger
	      if ( nmatched>0 && (tow->GetEt()<trigTowEt) ) { continue; }   // more than 1 trigger tower
	      else {							// first trigger tower
		trigTowEt = tow->GetEt();
		trigTowerPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
		nmatched += 1;
	      }
	    }
	  
	  }
	}
	if (nmatched>0) {
	  nHTtrig = nmatched;

	  if ( dRTrigLead<=R || dRTrigReco<=R ) {
	    
	    RunID = header->GetRunId();
	    EventID = Reader.GetNOfCurrentEvent();
	    nPrimary = header->GetNOfPrimaryTracks();
	    nGlobal = header->GetNGlobalTracks();
	    nVertices = header->GetNumberOfVertices();
	    refMult = header->GetReferenceMultiplicity();
	    gRefMult = header->GetGReferenceMultiplicity();
	    BbcAdcSumEast = header->GetBbcAdcSumEast();
	    BbcAdcSumEastOuter = header->GetBbcAdcSumEastOuter();
	    BbcAdcSumWest = header->GetBbcAdcSumWest();
	    BbcAdcSumWestOuter = header->GetBbcAdcSumWestOuter();
	    leadPt = leadJet.pt();
	    leadEta = leadJet.eta();
	    leadPhi = leadJet.phi();
	    leadArea = leadJet.area();
	    recoPt = recoJet.pt();
	    recoEta = recoJet.eta();
	    recoPhi = recoJet.phi();
	    recoArea = recoJet.area();

	    //  BACKGROUND ESTIMATION
	    GatherChargedBG( leadJet, container, chgParticles );   // gather BG
	    GatherNeutralBG( leadJet, container, neuParticles );
    
	    nBGpart_chg = chgParticles.size();
	    nBGpart_neu = neuParticles.size();

	    chgEastSum = 0;            chgMidSum = 0;            chgWestSum = 0;            neuEastSum = 0;            neuMidSum = 0;            neuWestSum = 0;
	    CalculateRhoByChargeAndEta(chgParticles,neuParticles,chgEastSum,chgMidSum,chgWestSum,neuEastSum,neuMidSum,neuWestSum);
	    chgEastRho = chgEastSum/eastArea;
	    chgMidRho = chgMidSum/midArea;
	    chgWestRho = chgWestSum/westArea;
	    neuEastRho = neuEastSum/eastArea;
	    neuMidRho = neuMidSum/midArea;
	    neuWestRho = neuWestSum/westArea;
    
	    HTdijetTree->Fill();
	  }
	}
      }
    }
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hChgBgPtEtaPhi->Write();  //  WRITE HISTOGRAMS & TREE
  hNeuBgPtEtaPhi->Write();
  HTdijetTree->Write();  

  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
