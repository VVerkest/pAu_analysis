//  pAuHTjetUE.cxx
//  Veronica Verkest		November 8, 2019

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!!

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTjetUE.root"; number_of_events=-1; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );
  string efficFile = "src/trackeffic.root";

  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nBGpart_chg, nBGpart_neu, nJetsAbove5, nHTtrig;
  double Vz, BbcAdcSumEast, BbcAdcSumEastOuter, BbcAdcSumWest, BbcAdcSumWestOuter, leadPt, leadEta, leadPhi,
    chgEastRho, chgMidRho, chgWestRho, chgEastRho_te, chgMidRho_te, chgWestRho_te, neuEastRho, neuMidRho, neuWestRho, leadArea, dPhiTrigLead, dRTrigLead;

  HTjetTree->Branch( "RunID", &RunID );
  HTjetTree->Branch( "EventID", &EventID );
  HTjetTree->Branch( "nTowers", &nTowers );
  HTjetTree->Branch( "nPrimary", &nPrimary );
  HTjetTree->Branch( "nGlobal", &nGlobal );
  HTjetTree->Branch( "nVertices", &nVertices );
  HTjetTree->Branch( "refMult", &refMult );
  HTjetTree->Branch( "gRefMult", &gRefMult );
  HTjetTree->Branch( "Vz", &Vz );
  HTjetTree->Branch( "leadPt", &leadPt );
  HTjetTree->Branch( "leadEta", &leadEta );
  HTjetTree->Branch( "BbcAdcSumEast", &BbcAdcSumEast );
  HTjetTree->Branch( "BbcAdcSumEastOuter", &BbcAdcSumEastOuter );
  HTjetTree->Branch( "BbcAdcSumWest", &BbcAdcSumWest );
  HTjetTree->Branch( "BbcAdcSumWestOuter", &BbcAdcSumWestOuter );
  HTjetTree->Branch( "leadPhi", &leadPhi );
  HTjetTree->Branch( "chgEastRho", &chgEastRho );
  HTjetTree->Branch( "chgMidRho", &chgMidRho );
  HTjetTree->Branch( "chgWestRho_te", &chgWestRho_te );
  HTjetTree->Branch( "chgEastRho_te", &chgEastRho_te );
  HTjetTree->Branch( "chgMidRho_te", &chgMidRho_te );
  HTjetTree->Branch( "chgWestRho", &chgWestRho );
  HTjetTree->Branch( "neuEastRho", &neuEastRho );
  HTjetTree->Branch( "neuMidRho", &neuMidRho );
  HTjetTree->Branch( "neuWestRho", &neuWestRho );
  HTjetTree->Branch( "leadArea", &leadArea );
  HTjetTree->Branch( "nBGpart_chg", &nBGpart_chg );
  HTjetTree->Branch( "nBGpart_neu", &nBGpart_neu );
  HTjetTree->Branch( "nJetsAbove5", &nJetsAbove5 );
  HTjetTree->Branch( "nHTtrig", &nHTtrig );
  HTjetTree->Branch( "dPhiTrigLead", &dPhiTrigLead );
  HTjetTree->Branch( "dRTrigLead", &dRTrigLead );
       
  TH3D *hChgBgPtEtaPhi = new TH3D( "hChgBgPtEtaPhi", "Charged Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 30,0,15, 20,-1.0,1.0, 60,0.0,2*pi );
  TH3D *hNeuBgPtEtaPhi = new TH3D( "hNeuBgPtEtaPhi", "Neutral Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 30,0,15, 20,-1.0,1.0, 60,0.0,2*pi );
  TH3D *hAllJets = new TH3D( "hAllJets", "All jets p_{T}>=5.0 GeV;p_{T} (GeV);#eta;#phi", 160,0,80, 20,-1.0,1.0, 60,0.0,2*pi );

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
  Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
  Selector ptMaxSelector = SelectorPtMax( 30.0 );

  PseudoJet leadJet, trigTowerPJ, towPJ;
  vector<PseudoJet> rawParticles, rawJets, allJets, chgParticles, neuParticles, BGparticles;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  double deltaPhi, deltaR;       double trigTowEta, trigTowPhi;

  Reader.NextEvent();
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);
    rawParticles.clear();    chgParticles.clear();    neuParticles.clear();    rawJets.clear();    allJets.clear();    //  CLEAR VECTORS

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

    // int trigTowId;
    // TStarJetPicoTriggerInfo *trig;
    // TStarJetPicoTower *tow, *triggerTower;
    // double trigTowEt = 0.0;
    // std::vector<int> trigTowers;
    // for ( int i=0; i<event->GetTrigObjs()->GetEntries(); ++i ) {
    //   trig = (TStarJetPicoTriggerInfo *)event->GetTrigObj(i);
    //   if ( trig->isBHT2() && UseTriggerTower( trig->GetId()) ) { trigTowers.push_back( trig->GetId() ); }
    // }
    // std::sort(trigTowers.begin(), trigTowers.end());

    // int nmatched = 0;
    // for (int i=0; i<nTowers; ++i){
    //   tow = (TStarJetPicoTower*)SelectedTowers->At(i);
    //   if ( tow->GetEt()>=5.4 && std::count(trigTowers.begin(), trigTowers.end(), tow->GetId())) {
    // 	if ( nmatched>0 && (tow->GetEt()<trigTowEt) ) { continue; }
    // 	else {
    // 	  trigTowEt = tow->GetEt();
    // 	  trigTowerPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
    // 	  nmatched += 1;
    // 	}
    //   }
    // }
    // if (nmatched==1) { // only accept events with 1 HT triggers
    //   nHTtrig = nmatched;
	
    Vz = header->GetPrimaryVertexZ();
    //if ( UseHTevent( header, event, vzCut, Vz ) == false ) { continue; } // Skip events based on: Run#, vz cut, BBCSumE; only accept HT events
	
    GatherParticles( container, rawParticles );

    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
    ClusterSequenceArea jetCluster( rawParticles, jet_def, area_def); 
      
    leadJetSelector = leadPtMinSelector && jetEtaSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN >=10 GeV
    
    if ( rawJets.size()>0 ) {
      leadJet = rawJets[0]; 

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
	  deltaPhi = fabs( leadJet.delta_phi_to( towPJ ) );
	  deltaR = leadJet.delta_R( towPJ );

	  if ( deltaR<=R || fabs(deltaPhi)>=(pi-R) ) {  // require trigger
	    if ( nmatched>0 && (tow->GetEt()<trigTowEt) ) { continue; }   // more than 1 trigger tower
	    else {							// first trigger tower
	      trigTowEt = tow->GetEt();
	      trigTowerPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
	      nmatched += 1;
	    }
	  }
	  
	}
      }
      //if (nmatched==1) { // only accept events with 1 HT triggers
      if (nmatched>0) {
	nHTtrig = nmatched;
	
	dPhiTrigLead = fabs( leadJet.delta_phi_to( trigTowerPJ ) );
	dRTrigLead = leadJet.delta_R( trigTowerPJ );

	if ( dRTrigLead<=R || fabs(dPhiTrigLead)>=(pi-R) ) {
	  
	  Selector allJetSelector = SelectorPtMin(2.0) && ptMaxSelector && jetEtaSelector;
	  allJets = sorted_by_pt( allJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN >=5 GeV
	  
	  nJetsAbove5 = allJets.size();
	  for ( int i=0; i<allJets.size(); ++i ) {
	    hAllJets->Fill( allJets[i].pt(), allJets[i].eta(), allJets[i].phi() );
	  }
	  
	  
	  RunID = header->GetRunId();
	  EventID = Reader.GetNOfCurrentEvent();
	  // TList *SelectedTowers = Reader.GetListOfSelectedTowers();
	  // nTowers = CountTowers( SelectedTowers );
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

	  //  BACKGROUND ESTIMATION
	  GatherChargedBG( leadJet, container, chgParticles );
	  GatherNeutralBG( leadJet, container, neuParticles );

	  nBGpart_chg = chgParticles.size();
	  nBGpart_neu = neuParticles.size();
    
	  chgEastSum = 0;  chgMidSum = 0;  chgWestSum = 0;  neuEastSum = 0;  neuMidSum = 0;  neuWestSum = 0;
	  CalculateRhoByChargeAndEta(chgParticles,neuParticles,chgEastSum,chgMidSum,chgWestSum,neuEastSum,neuMidSum,neuWestSum,hChgBgPtEtaPhi,hNeuBgPtEtaPhi);
	  chgEastRho = chgEastSum/eastArea;
	  chgMidRho = chgMidSum/midArea;
	  chgWestRho = chgWestSum/westArea;
	  neuEastRho = neuEastSum/eastArea;
	  neuMidRho = neuMidSum/midArea;
	  neuWestRho = neuWestSum/westArea;

	  chgParticles.clear(); // clear vector!
	  GatherChargedBGwithEfficiency( leadJet, container, chgParticles, efficFile );   // gather BG
	  CalculateBGsubtractedChargedRho(chgParticles,chgEastSum,chgMidSum,chgWestSum);
	  chgEastRho_te = chgEastSum/eastArea;
	  chgMidRho_te = chgMidSum/midArea;
	  chgWestRho_te = chgWestSum/westArea;
    
	  HTjetTree->Fill();
	}
      }
    }
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hChgBgPtEtaPhi->Write();  //  WRITE HISTOGRAMS & TREE
  hNeuBgPtEtaPhi->Write();
  hAllJets->Write();
  HTjetTree->Write();  

  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
