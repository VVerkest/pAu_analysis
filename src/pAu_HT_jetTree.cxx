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
    outFile = "out/HTjets/pAu_HT_BBCEsumOver64000.root";
    nEvents = -1;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  
  double dca, pt, pz, eta, phi, pmin1, pmax1, pmin2, pmax2;
  
  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, sigma;
  vector<int> towID, nHitsPoss, nHitsFit, Charge, nCons;
  vector<double> DCA, towEt, towEta, towPhi, trEta, trPhi, trPz, trPt, jetPt, jetEta, jetPhi, jetEt;

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
  HTjetTree->Branch("nCons", &nCons);                 HTjetTree->Branch("leadPt", &leadPt);
  HTjetTree->Branch("leadEta", &leadEta);               HTjetTree->Branch("leadPhi", &leadPhi);         HTjetTree->Branch("leadEt", &leadEt);
  HTjetTree->Branch("leadNcons", &leadNcons);    HTjetTree->Branch("subPt", &subPt);             HTjetTree->Branch("subEta", &subEta);
  HTjetTree->Branch("subPhi", &subPhi);                HTjetTree->Branch("subEt", &subEt);             HTjetTree->Branch("subNcons", &subNcons);
  HTjetTree->Branch("rho", &rho);                           HTjetTree->Branch("sigma", &sigma);
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);
  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  JetDefinition bg_jet_def(kt_algorithm, R);     //  BACKGROUND ESTIMATION JET DEFINITION
  //  cambridge_algorithm??
  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();     rawJets.clear();     nCons.clear();     towID.clear();     Charge.clear();     towEt.clear();     towEta.clear();     towPhi.clear();
    trEta.clear();     trPhi.clear();     trPz.clear();     trPt.clear();     jetPt.clear();     jetEta.clear();     jetPhi.clear();     jetEt.clear();     DCA.clear();
    chgParticles.clear();      neuParticles.clear();     nHitsPoss.clear();     nHitsFit.clear();        //  clear vectors

    Reader.PrintStatus(5);

    event = Reader.GetEvent();    header = event->GetHeader();    container = Reader.GetOutputContainer();
    
    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }
    if (!(header->HasTriggerId(500401) || header->HasTriggerId(500411))) {continue;}   //  ONLY SELECT JP2 TRIGGER EVENTS
    Vz = header->GetPrimaryVertexZ();           if ( abs(Vz) > vzCut ) { continue; }
    if ( header->GetBbcAdcSumEast() < 64000 ) { continue; }    


    //   JET-FINDING
    GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT>=2.0GeV    and    |eta|<1.0
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS
    const double AREA = 4.*(pi/3.);

    nJets=0;

    if ( rawJets.size()<2) { continue; }                                                       //  REQUIRE DIJET
    double phi1 = rawJets[0].phi();     double phi2 = rawJets[1].phi();
    double dphi = fabs( fabs( phi1 - phi2 ) - pi );
    if ( dphi> R || rawJets[0].pt()<10.0 || rawJets[1].pt()<2.0 ) { continue; }
    else {                        //  CREATE DIJET PAIR  
      leadPt = rawJets[0].pt();      leadEta = rawJets[0].eta();      leadPhi = rawJets[0].phi();      leadEt = rawJets[0].Et();
      vector<PseudoJet> LeadCons= rawJets[0].constituents();      leadNcons = LeadCons.size();
      subPt = rawJets[1].pt();      subEta = rawJets[1].eta();      subPhi = rawJets[1].phi();      subEt = rawJets[1].Et();
      vector<PseudoJet> SubCons= rawJets[1].constituents();      subNcons = SubCons.size();
    }

    
    for ( int i=0; i<rawJets.size(); ++i ) {                              //  FILL JET INFO
      jetPt.push_back( rawJets[i].pt() );
      jetEta.push_back( rawJets[i].eta() );
      jetPhi.push_back( rawJets[i].phi() );
      jetEt.push_back( rawJets[i].Et() );
      vector<PseudoJet> Cons= rawJets[i].constituents();
      nCons.push_back( Cons.size() );
      nJets+=1;
    }

    
    
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
      eta = event->GetTower(i)->GetEta();
      if ( abs(eta) > etaCut ) { nTowers -= 1;   continue; }            // tower eta cut
      towEta.push_back(eta);
      towPhi.push_back(event->GetTower(i)->GetPhi());
      towID.push_back(event->GetTower(i)->GetId());
    }


    GatherChargedBG( rawJets[0], container, chgParticles);     GatherNeutralBG( rawJets[0], container, neuParticles);


    //  BACKGROUND ESTIMATION 
    double ptSum = 0;
    for (int i=0; i<chgParticles.size(); ++i) {      ptSum+=chgParticles[i].pt();    }
    for (int i=0; i<neuParticles.size(); ++i) {      ptSum+=neuParticles[i].pt();    }
		      

    rho = ptSum / AREA;
    
    
    HTjetTree->Fill();                                                           //  FILL TREE
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  HTjetTree->Write();                            //  WRITE TREES

  pAuFile->Close();
  
  return 0;
}
