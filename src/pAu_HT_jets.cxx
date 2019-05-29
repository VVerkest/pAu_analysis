//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019

//  (from Isaac's code:)
//for pAu, need to ask if the event has the trigger. The trigger IDs are:
//HT2*BBCMB : 500205, 500215
//JP2 : 500401, 500411
//BBCMB : 500008, 500018
//VPDMB :  500904


#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

// tracks and towers have eta and DCA cuts specified in pAuFunctions.hh

int main ( int argc, const char** argv ) {
  
  string inFile, outFile;
  int nEvents;
  if ( argc ==  4 ) {
    vector<string> arguments( argv+1, argv+argc );
    inFile = arguments[0];
    outFile = arguments[1];
    nEvents = atoi(arguments[2].c_str());
  }
  else if ( argc ==  1 ) {
    vector<string> arguments( argv+1, argv+argc );
    inFile = "production_pAu200_2015/HT/pAu_2015_200_HT*.root";
    outFile = "out/pAu_HT.root";
    nEvents = 100000;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in HT: ????
  InitReader( Reader, Chain, numEvents );
  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  
  double dca, pt, pz, eta, phi;
  
  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nHitsPoss, nHitsFit, nJets, leadNcons, subNcons;
  double Vx, Vy, Vz, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt;
  vector<int> towID, Charge, nCons;
  vector<double> towEt, towEta, towPhi, trEta, trPhi, trPz, trPt, jetPt, jetEta, jetPhi, jetEt;

  HTjetTree->Branch("EventID", &EventID);            HTjetTree->Branch("RunID", &RunID);                  HTjetTree->Branch("Vz", &Vz);
  HTjetTree->Branch("nTowers", &nTowers);          HTjetTree->Branch("nPrimary", &nPrimary);         HTjetTree->Branch("nGlobal", &nGlobal);
  HTjetTree->Branch("nVertices", &nVertices);        HTjetTree->Branch("refMult", &refMult);              HTjetTree->Branch("gRefMult", &gRefMult);
  HTjetTree->Branch("BbcCoincidenceRate", &BbcCoincidenceRate);                    HTjetTree->Branch("BbcEastRate", &BbcEastRate);
  HTjetTree->Branch("BbcWestRate", &BbcWestRate);                                           HTjetTree->Branch("vpdVz", &vpdVz);
  HTjetTree->Branch("towEt", &towEt);                   HTjetTree->Branch("towEta", &towEta);	          HTjetTree->Branch("towPhi", &towPhi);
  HTjetTree->Branch("towID", &towID);                  HTjetTree->Branch("nHitsPoss", &nHitsPoss);     HTjetTree->Branch("nHitsFit", &nHitsFit);
  HTjetTree->Branch("trPz", &trPz);
  HTjetTree->Branch("trPt", &trPt);                          HTjetTree->Branch("trEta", &trEta);                     HTjetTree->Branch("trPhi", &trPhi);
  HTjetTree->Branch("DCA", &DCA);                       HTjetTree->Branch("nJets", &nJets);                    HTjetTree->Branch("jetPt", &jetPt);
  HTjetTree->Branch("jetEta", &jetEta);                    HTjetTree->Branch("jetPhi", &jetPhi);                       HTjetTree->Branch("jetEt", &jetEt);
  HTjetTree->Branch("nCons", &nCons);                  HTjetTree->Branch("leadPt", &leadPt);
  HTjetTree->Branch("leadEta", &leadEta);                HTjetTree->Branch("leadPhi", &leadPhi);         HTjetTree->Branch("leadEt", &leadEt);
  HTjetTree->Branch("leadNcons", &leadNcons);     HTjetTree->Branch("subPt", &subPt);             HTjetTree->Branch("subEta", &subEta);
  HTjetTree->Branch("subPhi", &subPhi);                 HTjetTree->Branch("subEt", &subEt);             HTjetTree->Branch("subNcons", &subNcons);
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);
  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

  vector<PseudoJet> rawParticles, rawJets;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();     rawJets.clear();     nCons.clear();     jetPt.clear();     jetEta.clear();     jetPhi.clear();     jetEt.clear();     towID.clear();     Charge.clear();     towEt.clear();
    towEta.clear();     towPhi.clear();     trEta.clear();     trPhi.clear();     trPz.clear();     trPt.clear();     jetPt.clear();     jetEta.clear();     jetPhi.clear();     jetEt.clear();       //  clear vectors

    Reader.PrintStatus(5);

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    if (!(header->HasTriggerId(500401) || header->HasTriggerId(500411))) {continue;}   //  ONLY SELECT JP2 TRIGGER EVENTS
    
    Vz = header->GetPrimaryVertexZ();
    if ( abs(Vz) > vzCut ) { continue; }

    eID = Reader.GetNOfCurrentEvent();          EventID = eID;
    rID = header->GetRunId();                        RunID = rID;

    int npt = header->GetNOfPrimaryTracks();      nPrimary = npt;
    int ntow = header->GetNOfTowers();               nTowers = ntow;
    Vx = header->GetPrimaryVertexX();                                           Vy = header->GetPrimaryVertexY();
    nGlobal = header->GetNGlobalTracks();                                    nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();                           gRefMult = header->GetGReferenceMultiplicity();
    BbcCoincidenceRate = header->GetBbcCoincidenceRate();        vpdVz = header->GetVpdVz();
    BbcEastRate = header->GetBbcEastRate();                                  BbcWestRate = header->GetBbcWestRate();
    
    for ( int i=0; i<npt; ++i ) {                                         //  FILL EVENT DATA

      dca = event->GetPrimaryTrack(i)->GetDCA();
      if ( dca > dcaCut ) { nPrimary -= 1;   continue; }                   // track DCA cut
      DCA.push_back( dca );
      eta = (double) event->GetPrimaryTrack(i)->GetEta();
      if ( abs(eta) > etaCut ) { nPrimary -= 1;   continue; }            // track eta cut
      trEta.push_back(eta);
      trPt.push_back((double) event->GetPrimaryTrack(i)->GetPt());
      trPz.push_back((double) event->GetPrimaryTrack(i)->GetPz());
      trPhi.push_back((double) event->GetPrimaryTrack(i)->GetPhi());

      nHitsFit.push_back(event->GetPrimaryTrack(i)->GetNOfFittedHits());
      nHitsPoss.push_back(event->GetPrimaryTrack(i)->GetNOfPossHits());
    }

    
    for ( int i=0; i<ntow; ++i ) {                                         //  FILL TOWER DATA
      towEt.push_back(event->GetTower(i)->GetEt());
      eta = event->GetTower(i)->GetEta());
    if ( abs(eta) > etaCut ) { nTowers -= 1;   continue; }            // tower eta cut
    towEta.push_back(eta);
    towPhi.push_back(event->GetTower(i)->GetPhi());
    towID.push_back(event->GetTower(i)->GetId());
  }



  //   JET-FINDING
  GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT>=2.0GeV    and    |eta|<1.0

  ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
  vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

  nJets=0;
    
  for ( int i=0; i<rawJets.size(); ++i ) {                              //  FILL JET INFO
    if ( rawJets[i].pt()<jetMinPt )  { cerr<<"bad jet pt selector?"<<endl;   continue; }
    jetPt.push_back( rawJets[i].pt() );
    if ( abs(rawJets[i].eta())>etaCut ) { cerr<<"bad jet eta selector?"<<endl;   continue; }
    jetEta.push_back( rawJets[i].eta() );
    jetPhi.push_back( rawJets[i].phi() );
    jetEt.push_back( rawJets[i].Et() );
    vector<PseudoJet> Cons= rawJets[i].constituents();
    nCons.push_back( Cons.size() );
    nJets+=1;
  }


  if ( rawJets.size()>=2 ) {                        //  CREATE DIJET PAIR
    double dphi = fabs( fabs( rawJets.at(0).delta_phi_to( rawJets.at(1) ) ) - pi );
    if ( dphi< R && rawJets[0].pt()>2.0 && rawJets[1].pt()>2.0 ) {  
      leadPt = rawJets[0].pt();
      leadEta = rawJets[0].eta();
      leadPhi = rawJets[0].phi();
      leadEt = rawJets[0].Et();
      vector<PseudoJet> LeadCons= rawJets[0].constituents();
      leadNcons = LeadCons.size();
      subPt = rawJets[0].pt();
      subEta = rawJets[0].eta();
      subPhi = rawJets[0].phi();
      subEt = rawJets[0].Et();
      vector<PseudoJet> SubCons= rawJets[1].constituents();
      subNcons = SubCons.size();
    }
  }

  // double djAxisPhi = ( rawJets[0].phi() + rawJets[1].phi() )/2;
  // Selector bgCircle = SelectorCircle(R);
  // Selector bgRapRange = SelectorRapRange( -0.6, 0.6 );
  // Selector bgSelector = bgCircle && bgRapRange;
  // double ghost_maxrap = 1.0;
  // AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
  // JetMedianBackgroundEstimator UE( bgSelector, jet_def, area_def);
    
  HTjetTree->Fill();                                                           //  FILL TREE
}
// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    
HTjetTree->Write();                            //  WRITE TREES


pAuFile->Close();
  
  
return 0;
}
