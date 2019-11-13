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
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTdijetUE.root"; number_of_events=100000; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTdijetTree = new TTree( "HTdijetTree", "HT_DijetTree" );
  
  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, eval, pval;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  HTdijetTree->Branch( "RunID", &RunID );			HTdijetTree->Branch( "EventID", &EventID );				HTdijetTree->Branch( "nTowers", &nTowers );
  HTdijetTree->Branch( "nPrimary", &nPrimary );	       	HTdijetTree->Branch( "nGlobal", &nGlobal );				HTdijetTree->Branch( "nVertices", &nVertices );
  HTdijetTree->Branch( "refMult", &refMult );			HTdijetTree->Branch( "gRefMult", &gRefMult );			HTdijetTree->Branch( "Vz", &Vz );
  HTdijetTree->Branch( "leadPt", &leadPt );			HTdijetTree->Branch( "BbcAdcEastSum", &BbcAdcEastSum );	HTdijetTree->Branch( "leadEta", &leadEta );
  HTdijetTree->Branch( "leadPhi", &leadPhi );			HTdijetTree->Branch( "chgEastRho", &chgEastRho );		HTdijetTree->Branch( "chgMidRho", &chgMidRho );
  HTdijetTree->Branch( "chgWestRho", &chgWestRho );	HTdijetTree->Branch( "neuEastRho", &neuEastRho );		HTdijetTree->Branch( "neuMidRho", &neuMidRho );
  HTdijetTree->Branch( "neuWestRho", &neuWestRho );
	
  TH2D *hChgBgEtaPhi = new TH2D( "hChgBgEtaPhi", "Charged Background #phi vs. #eta;#eta;#phi", 40,-1.0,1.0, 120,0.0,2*pi );
  TH2D *hNeuBgEtaPhi = new TH2D( "hNeuBgEtaPhi", "Neutral Background #phi vs. #eta;#eta;#phi", 40,-1.0,1.0, 120,0.0,2*pi );

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
  Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
  Selector ptMaxSelector = SelectorPtMax( 30.0 );

  PseudoJet leadJet, recoJet;  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles, BGparticles, recoCandidates;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );


  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);
    rawParticles.clear();            chgParticles.clear();            neuParticles.clear();            rawJets.clear();            recoCandidates.clear();             //  CLEAR VECTORS

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, event, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept HT Trigger events

    GatherParticles( container, rawParticles );
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS

    leadJetSelector = leadPtMinSelector && ptMaxSelector && jetEtaSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN 10-30 GeV RANGE
    if ( rawJets.size()>0 ) { leadJet = rawJets[0]; }
    else { continue; }

    for ( int p=0; p<3; ++p ) {
      if ( leadJet.pt() >= ptLo[p]  &&  leadJet.pt() <= ptHi[p] ) { pval = p; }
    }
    for ( int e=0; e<3; ++e ) {
      if ( leadJet.eta() >= etaLo[e]  &&  leadJet.eta() <= etaHi[e] ) { eval = e; }
    }
    if ( pval==99 || eval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl; }
    

    //   REQUIRE RECOIL JET TO HAVE AT LEAST HALF OF LEAD PT AND BE IN THE SAME ETA RANGE!
    Selector recoPtRangeSelector = SelectorPtRange( leadJet.pt()/2, ptHi[pval] );          //  JET pT RANGE    { 10-15, 15-20, 20-30 }
    //Selector etaRangeSelector = SelectorEtaRange( etaLo[eval], etaHi[eval] );          //  JET eta RANGE

    Selector recoJetSelector = recoPtRangeSelector && etaRangeSelector && jetEtaSelector;
    
    recoCandidates = sorted_by_pt( recoJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    
    bool hasReco = false;
    
    if ( recoCandidates.size()>0 ) {

      for ( int i=0; i<recoCandidates.size(); ++i ) {
	double deltaPhi_abs = fabs( recoCandidates[i].phi() - leadJet.phi() );
	if ( fabs( deltaPhi_abs - pi ) <= R ) {
	  recoJet = recoCandidates[i];
	  // hRecoVsLeadPt->Fill( leadJet.pt(), recoJet.pt() );
	  // hRecoAbsDeltaPhi->Fill( deltaPhi_abs );
	  // hRecoVsLeadEta->Fill( leadJet.eta(), recoJet.eta() );
	  // hRecoPhi->Fill( recoJet.phi() );

	  hasReco = true;
	}
	if ( hasReco == true ) { continue; }  // exit loop with highest-pt jet in recoil range
      }

    }
    else { continue; }
    if ( hasReco == false ) { continue; }





    
    RunID = header->GetRunId();
    EventID = Reader.GetNOfCurrentEvent();
    TList *SelectedTowers = Reader.GetListOfSelectedTowers();
    nTowers = CountTowers( SelectedTowers );
    nPrimary = header->GetNOfPrimaryTracks();
    nGlobal = header->GetNGlobalTracks();
    nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();
    gRefMult = header->GetGReferenceMultiplicity();
    BbcAdcEastSum = header->GetBbcAdcSumEast();
    leadPt = leadJet.pt();
    leadEta = leadJet.eta();
    leadPhi = leadJet.phi();

    //  BACKGROUND ESTIMATION
    GatherChargedBG( leadJet, container, chgParticles );
    GatherNeutralBG( leadJet, container, neuParticles );
    chgEastSum = 0;            chgMidSum = 0;            chgWestSum = 0;            neuEastSum = 0;            neuMidSum = 0;            neuWestSum = 0;
    CalculateRhoByChargeAndEta( chgParticles, neuParticles, chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum, hChgBgEtaPhi, hNeuBgEtaPhi ); 
    chgEastRho = chgEastSum/eastArea;
    chgMidRho = chgMidSum/midArea;
    chgWestRho = chgWestSum/westArea;
    neuEastRho = neuEastSum/eastArea;
    neuMidRho = neuMidSum/midArea;
    neuWestRho = neuWestSum/westArea;
    
    HTdijetTree->Fill();
    
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hChgBgEtaPhi->Write();  //  WRITE HISTOGRAMS & TREE
  hNeuBgEtaPhi->Write();
  HTdijetTree->Write();  

  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
