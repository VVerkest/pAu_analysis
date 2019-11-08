//  pAuHTjetUE.cxx
//  Veronica Verkest		November 8, 2019

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTjetUE.root"; number_of_events=100000; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );

  double BGeta, chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult;
  double Vz, BbcAdcEastSum, leadPt, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho;

  HTjetTree->Branch( "RunID", &RunID );
  HTjetTree->Branch( "EventID", &EventID );
  HTjetTree->Branch( "nTowers", &nTowers );
  HTjetTree->Branch( "nPrimary", &nPrimary );
  HTjetTree->Branch( "nGlobal", &nGlobal );
  HTjetTree->Branch( "nVertices", &nVertices );
  HTjetTree->Branch( "refMult", &refMult );
  HTjetTree->Branch( "gRefMult", &gRefMult );
  HTjetTree->Branch( "Vz", &Vz );
  HTjetTree->Branch( "BbcAdcEastSum", &BbcAdcEastSum );
  HTjetTree->Branch( "leadPt", &leadPt );
  HTjetTree->Branch( "leadEta", &leadEta );
  HTjetTree->Branch( "leadPhi", &leadPhi );
  HTjetTree->Branch( "chgEastRho", &chgEastRho );
  HTjetTree->Branch( "chgMidRho", &chgMidRho );
  HTjetTree->Branch( "chgWestRho", &chgWestRho );
  HTjetTree->Branch( "neuEastRho", &neuEastRho );
  HTjetTree->Branch( "neuMidRho", &neuMidRho );
  HTjetTree->Branch( "neuWestRho", &neuWestRho );
  
  TH2D *hChgBgEtaPhi = new TH2D( "hChgBgEtaPhi", "Charged Background #phi vs. #eta", 40,-1.0,1.0, 120,0.0,2*pi );
  TH2D *hNeuBgEtaPhi = new TH2D( "hNeuBgEtaPhi", "Nuetral Background #phi vs. #eta", 40,-1.0,1.0, 120,0.0,2*pi );

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
  Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
  Selector ptMaxSelector = SelectorPtMax( 30.0 );

  PseudoJet leadJet;  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles, BGparticles;
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
    rawParticles.clear();            chgParticles.clear();            neuParticles.clear();            rawJets.clear();            //  CLEAR VECTORS

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
    chgEastSum = 0;
    chgMidSum = 0;
    chgWestSum = 0;
    neuEastSum = 0;
    neuMidSum = 0;
    neuWestSum = 0;

    GatherChargedBG( leadJet, container, chgParticles );
    for ( int i=0; i<chgParticles.size(); ++i ) {
      BGeta = chgParticles[i].eta();
      hChgBgEtaPhi->Fill( BGeta, chgParticles[i].phi() );
      if ( BGeta>=etaLoMid && BGeta<etaHiMid ) { chgMidSum+= chgParticles[i].pt(); }
      else if ( BGeta>=etaLoEast && BGeta<=etaHiEast ) { chgEastSum+= chgParticles[i].pt(); }
      else if ( BGeta>etaLoWest && BGeta<=etaHiWest ) { chgWestSum+= chgParticles[i].pt(); }
      else { cerr<<"error with chg BG particle eta"<<endl; }
    }
    
    GatherNeutralBG( leadJet, container, neuParticles );
    for ( int i=0; i<neuParticles.size(); ++i ) {
      BGeta = neuParticles[i].eta();
      hNeuBgEtaPhi->Fill( BGeta, neuParticles[i].phi() );
      if ( BGeta>=etaLoMid && BGeta<etaHiMid ) { neuMidSum+= neuParticles[i].pt(); }
      else if ( BGeta>=etaLoEast && BGeta<=etaHiEast ) { neuEastSum+= neuParticles[i].pt(); }
      else if ( BGeta>etaLoWest && BGeta<=etaHiWest ) { neuWestSum+= neuParticles[i].pt(); }
      else { cerr<<"error with neu BG particle eta"<<endl; }
    }

    chgEastRho = chgEastSum/eastArea;
    chgMidRho = chgMidSum/midArea;
    chgWestRho = chgWestSum/westArea;
    neuEastRho = neuEastSum/eastArea;
    neuMidRho = neuMidSum/midArea;
    neuWestRho = neuWestSum/westArea;

    HTjetTree->Fill();
    
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  //  ~ ~ ~ ~ ~ ~ ~ ~ WRITE HISTOGRAMS & TREE ~ ~ ~ ~ ~ ~ ~ ~
  hChgBgEtaPhi->Write();
  hNeuBgEtaPhi->Write();
  HTjetTree->Write();  
  
  pAuFile->Close();

  return 0;
}
