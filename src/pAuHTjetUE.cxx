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
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTjetUE.root"; number_of_events=100000; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );
  string efficFile = "src/trackeffic.root";

  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nBGpart_chg, nBGpart_neu;
  double Vz, BbcAdcSumEast, BbcAdcSumEastOuter, BbcAdcSumWest, BbcAdcSumWestOuter, leadPt, leadEta, leadPhi,
    chgEastRho, chgMidRho, chgWestRho, neuEastRho, neuMidRho, neuWestRho, leadArea;

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
  HTjetTree->Branch( "chgWestRho", &chgWestRho );
  HTjetTree->Branch( "neuEastRho", &neuEastRho );
  HTjetTree->Branch( "neuMidRho", &neuMidRho );
  HTjetTree->Branch( "neuWestRho", &neuWestRho );
  HTjetTree->Branch( "leadArea", &leadArea );
  HTjetTree->Branch( "nBGpart_chg", &nBGpart_chg );
  HTjetTree->Branch( "nBGpart_neu", &nBGpart_neu );
       
  TH3D *hChgBgPtEtaPhi = new TH3D( "hChgBgPtEtaPhi", "Charged Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-1.0,1.0, 120,0.0,2*pi );
  TH3D *hNeuBgPtEtaPhi = new TH3D( "hNeuBgPtEtaPhi", "Neutral Background #phi vs. #eta;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-1.0,1.0, 120,0.0,2*pi );

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
    if ( UseHTevent( header, event, vzCut, Vz ) == false ) { continue; } // Skip events based on: Run#, vz cut, BBCSumE; only accept HT events
    if ( header->GetBbcAdcSumEast() < 3559.12 ) { continue; }     //  neglect 0-10% event activity
    
    GatherParticles( container, rawParticles );

    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
    ClusterSequenceArea jetCluster( rawParticles, jet_def, area_def); 

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
    BbcAdcSumEast = header->GetBbcAdcSumEast();
    BbcAdcSumEastOuter = header->GetBbcAdcSumEastOuter();
    BbcAdcSumWest = header->GetBbcAdcSumWest();
    BbcAdcSumWestOuter = header->GetBbcAdcSumWestOuter();
    leadPt = leadJet.pt();
    leadEta = leadJet.eta();
    leadPhi = leadJet.phi();
    leadArea = leadJet.area();

    //  BACKGROUND ESTIMATION
    //GatherChargedBGwithEfficiency( leadJet, container, chgParticles, efficFile );   // gather BG
    GatherChargedBG( leadJet, container, chgParticles );
    GatherNeutralBG( leadJet, container, neuParticles );

    nBGpart_chg = chgParticles.size();
    nBGpart_neu = neuParticles.size();
    
    chgEastSum = 0;            chgMidSum = 0;            chgWestSum = 0;            neuEastSum = 0;            neuMidSum = 0;            neuWestSum = 0;
    CalculateRhoByChargeAndEta(chgParticles,neuParticles,chgEastSum,chgMidSum,chgWestSum,neuEastSum,neuMidSum,neuWestSum,hChgBgPtEtaPhi,hNeuBgPtEtaPhi);
    chgEastRho = chgEastSum/eastArea;
    chgMidRho = chgMidSum/midArea;
    chgWestRho = chgWestSum/westArea;
    neuEastRho = neuEastSum/eastArea;
    neuMidRho = neuMidSum/midArea;
    neuWestRho = neuWestSum/westArea;
    
    HTjetTree->Fill();
    
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hChgBgPtEtaPhi->Write();  //  WRITE HISTOGRAMS & TREE
  hNeuBgPtEtaPhi->Write();
  HTjetTree->Write();  

  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
