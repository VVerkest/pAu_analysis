//  pAu_QA.cxx
//  Veronica Verkest       August 23, 2019

//  COMMAND-LINE ARGUMENTS:
//  [0]: inFile (string)
//  [1]: outFile (string)
//  [2]: number_of_events (int ) {-1 for all}
//  [3]: bad_tower_option (string) {"allTowers","noBadTowers"}
//  [4]: trigger_option (string) {"HT","none"}

#include "pAuQAFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuQA;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile, bad_tower_option, trigger_option;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  6 ) {
    inFile = arguments[0];
    outFile = arguments[1];
    number_of_events = atoi(arguments[2].c_str());
    bad_tower_option = arguments[3];
    trigger_option = arguments[4];
  }
  else if ( argc==1 ) {
    inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root";
    outFile="out/pAuQA.root";
    number_of_events=1000000;
    bad_tower_option="noBadTowers";
    trigger_option="HT";
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH1D *hTowerFreq = new TH1D("hTowerFreq","Tower Frequency by ID;Tower ID", 4801,0,4801);
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicity;# of Towers", 700,0,700 );

  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 140,0,700, 40,0,200);
  TH2D *hTowEt = new TH2D("hTowEt","Tower E_{T} by ID;Tower ID;E_{T} (GeV)", 4801,0,4801, 120,-10,50.0);

  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta;Tower #phi", 60,0,30.0, 100,-1.0,1.0, 128,-pi,pi );
  TH3D *hAllJetsPtEtaPhi = new TH3D( "hAllJetsPtEtaPhi", "Inclusive Jets p_{T}, #eta, #phi;Jet p_{T} (GeV);Jet #eta;Jet #phi", 400,0.0,100.0, 40,-1.0,1.0, 120,0.0,2*pi  );
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

  double Vz;  double leadJetMinPt = 10.0;
  
  PseudoJet leadJet;  vector<PseudoJet> rawParticles, rawJets;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  TStarJetPicoTriggerInfo *trig;
  
  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, bad_tower_option, numEvents );

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    rawParticles.clear();            rawJets.clear();            //  CLEAR VECTORS
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, trigger_option, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept HT Trigger events

    GatherParticles( container, rawParticles );
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS

    Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
    Selector ptMinSelector = SelectorPtMin( jetMinPt );
    Selector allJetSelector = jetEtaSelector && ptMinSelector;

    rawJets = sorted_by_pt( allJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >2GeV
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtEtaPhi->Fill( rawJets[i].pt(), rawJets[i].eta(), rawJets[i].phi() ); }

    Selector leadPtMinSelector = SelectorPtMin( leadJetMinPt );
    Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >10GeV

    if ( rawJets.size()>0 ) { leadJet = rawJets[0]; }
    else { continue; }


    TList *SelectedTowers = Reader.GetListOfSelectedTowers();

    TStarJetPicoTower *tow;
    int nTowers = 0;
    for (int i=0; i<SelectedTowers->GetEntries(); ++i) {
      tow = (TStarJetPicoTower *) SelectedTowers->At(i);
      if ( fabs(tow->GetEta())>etaCut ) { continue; }
      nTowers+=1;
      hTowerFreq->Fill( tow->GetId() );
      hTowEt->Fill( tow->GetId(), tow->GetEt() );
      hTowEtEtaPhi->Fill( tow->GetEt(), tow->GetEta(), tow->GetPhi() );
    }
    
    hTowersPerEvent->Fill(nTowers);
    hPrimaryPerEvent->Fill( header->GetNOfPrimaryTracks() );
    hnPrimaryVSnTowers->Fill( nTowers, header->GetNOfPrimaryTracks() );
    hLeadPtEtaPhi->Fill( leadJet.pt(), leadJet.eta(), leadJet.phi() );
    
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  hAllJetsPtEtaPhi->Write();
  hLeadPtEtaPhi->Write();
  hnPrimaryVSnTowers->Write();
  hPrimaryPerEvent->Write();
  hTowerFreq->Write();
  hTowersPerEvent->Write();
  hTowEt->Write();
  hTowEtEtaPhi->Write();

  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
