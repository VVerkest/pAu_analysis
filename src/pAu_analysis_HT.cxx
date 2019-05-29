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
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 20,0,10, 24,-1.2,1.2, 128,-pi,pi );
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Frequency;# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Frequency (per run);Run no.;# of Towers", 40000,16124000,16164000, 140,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 40000,16124000,16164000, 40,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 140,0,700, 40,0,200);
  TH2D *hDCAvsPt = new TH2D("hDCAvsPt","DCA vs. p_{T};p_{T} GeV;DCA (cm)", 100,0.0,50.0, 70,0,3.5);
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 3500,0,3500000, 40,0,200 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 3500,0,3500000, 150,0,3000 );
  TH2D *hTowEt = new TH2D("hTowEt","Tower E_{T} by ID;Tower ID;E_{T} (GeV)", 4800,0,4800, 60,0,30.0);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta;Tower #phi", 60,0,30.0,-1.0,1.0, 2*pi,-pi,pi );
  TH1D *hTowerMult = new TH1D("hTowerMult","Tower Multiplicity by ID;Tower ID", 2400,0,4800);
  
  // TTree *HTtree = new TTree( "HTTree", "HTtree" );
  // TTree *HTtowers = new TTree( "HTTowers", "HTtowers" );
  // TTree *HTtracks = new TTree( "HTTracks", "HTtracks" );
  // TTree *HTjets = new TTree( "HTJets", "HTjets" );
  // TTree *HTdijets = new TTree( "HTDijets", "HTdijets" );
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, towID, Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz, towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz;
  double leadPt, leadEta, leadPhi, leadEt, leadNcons, subPt, subEta, subPhi, subEt, subNcons;
  int nJets;            vector<int> nCons;
  vector<double> jetPt, jetEta, jetPhi, jetEt;

  // HTtree->Branch("EventID", &EventID);            HTtree->Branch("RunID", &RunID);          HTtree->Branch("Vx", &Vx);
  // HTtree->Branch("Vy", &Vy);                             HTtree->Branch("Vz", &Vz);                     HTtree->Branch("nTowers", &nTowers);
  // HTtree->Branch("nPrimary", &nPrimary);         HTtree->Branch("nGlobal", &nGlobal);     HTtree->Branch("nVertices", &nVertices);
  // HTtree->Branch("refMult", &refMult);              HTtree->Branch("gRefMult", &gRefMult);
  // HTtree->Branch("BbcCoincidenceRate", &BbcCoincidenceRate);                    HTtree->Branch("BbcEastRate", &BbcEastRate);
  // HTtree->Branch("BbcWestRate", &BbcWestRate);                                           HTtree->Branch("vpdVz", &vpdVz);
  
  // HTtowers->Branch("EventID", &EventID);	  HTtowers->Branch("RunID", &RunID);  	  HTtowers->Branch("towEt", &towEt);
  // HTtowers->Branch("towEta", &towEta);	          HTtowers->Branch("towPhi", &towPhi);	  HTtowers->Branch("towID", &towID);

  // HTtracks->Branch("EventID", &EventID);	  HTtracks->Branch("RunID", &RunID);	  HTtracks->Branch("nHitsPoss", &nHitsPoss);
  // HTtracks->Branch("nHitsFit", &nHitsFit);	  HTtracks->Branch("trPx", &trPx);   	  HTtracks->Branch("trPy", &trPy);
  // HTtracks->Branch("trPz", &trPz);           	  HTtracks->Branch("trPt", &trPt);    	  HTtracks->Branch("trEta", &trEta);
  // HTtracks->Branch("trPhi", &trPhi);        	  HTtracks->Branch("DCA", &DCA);

  // HTjets->Branch("EventID", &EventID);            HTjets->Branch("RunID", &RunID);            HTjets->Branch("nJets", &nJets);
  // HTjets->Branch("jetPt", &jetPt);                      HTjets->Branch("jetEta", &jetEta);              HTjets->Branch("jetPhi", &jetPhi);
  // HTjets->Branch("jetEt", &jetEt);                      HTjets->Branch("nCons", &nCons);            // HTjets->Branch("", &);
  
  // HTdijets->Branch("EventID", &EventID);           HTdijets->Branch("RunID", &RunID);            HTdijets->Branch("leadPt", &leadPt);
  // HTdijets->Branch("leadEta", &leadEta);            HTdijets->Branch("leadPhi", &leadPhi);         HTdijets->Branch("leadEt", &leadEt);
  // HTdijets->Branch("leadNcons", &leadNcons);  HTdijets->Branch("subPt", &subPt);             HTdijets->Branch("subEta", &subEta);
  // HTdijets->Branch("subPhi", &subPhi);              HTdijets->Branch("subEt", &subEt);             HTdijets->Branch("subNcons", &subNcons);
  
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

    rawParticles.clear();  rawJets.clear();  nCons.clear();    jetPt.clear();    jetEta.clear();    jetPhi.clear();    jetEt.clear();  //  clear vectors

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

    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
    
    for ( int i=0; i<npt; ++i ) {                                         //  FILL EVENT DATA

      DCA = event->GetPrimaryTrack(i)->GetDCA();
      if ( DCA > dcaCut ) { nPrimary -= 1;   continue; }                   // track DCA cut
      trEta = (double) event->GetPrimaryTrack(i)->GetEta();
      if ( abs(trEta) > etaCut ) { nPrimary -= 1;   continue; }            // track eta cut

      trPt = (double) event->GetPrimaryTrack(i)->GetPt();
      trPx = (double) event->GetPrimaryTrack(i)->GetPx();
      trPy = (double) event->GetPrimaryTrack(i)->GetPy();
      trPz = (double) event->GetPrimaryTrack(i)->GetPz();
      trPhi = (double) event->GetPrimaryTrack(i)->GetPhi();

      nHitsFit = event->GetPrimaryTrack(i)->GetNOfFittedHits();
      nHitsPoss = event->GetPrimaryTrack(i)->GetNOfPossHits();

      hPrimaryTracks->Fill( trPt, trEta, trPhi );
      hDCAvsPt->Fill( trPt, DCA );
      // HTtracks->Fill();
    }

    
    for ( int i=0; i<ntow; ++i ) {                                         //  FILL TOWER DATA
      towEt = event->GetTower(i)->GetEt();
      towEta = event->GetTower(i)->GetEta();
      if ( abs(towEta) > etaCut ) { nTowers -= 1;   continue; }            // tower eta cut
      towPhi = event->GetTower(i)->GetPhi();
      towID = event->GetTower(i)->GetId();

      hTowEt->Fill( towID, towEt );
      hTowEtEtaPhi->Fill( towEt, towEta, towPhi );
      hTowerMult->Fill( towID );
      // HTtowers->Fill();
    }



    //   JET-FINDING
    GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT>=2.0GeV    and    |eta|<1.0

    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    // nJets=0;
    
    // for ( int i=0; i<rawJets.size(); ++i ) {                              //  FILL JET INFO
    //   if ( rawJets[i].pt()<jetMinPt )  { cerr<<"bad jet pt selector?"<<endl;   continue; }
    //   jetPt.push_back( rawJets[i].pt() );
    //   if ( abs(rawJets[i].eta())>etaCut ) { cerr<<"bad jet eta selector?"<<endl;   continue; }
    //   jetEta.push_back( rawJets[i].eta() );
    //   jetPhi.push_back( rawJets[i].phi() );
    //   jetEt.push_back( rawJets[i].Et() );
    //   vector<PseudoJet> Cons= rawJets[i].constituents();
    //   nCons.push_back( Cons.size() );
    //   nJets+=1;
    // }


    // if ( rawJets.size()>=2 ) {                        //  CREATE DIJET PAIR
    //   double dphi = fabs( fabs( rawJets.at(0).delta_phi_to( rawJets.at(1) ) ) - pi );
    //   if ( dphi< R && rawJets[0].pt()>2.0 && rawJets[1].pt()>2.0 ) {  
    // 	leadPt = rawJets[0].pt();
    // 	leadEta = rawJets[0].eta();
    // 	leadPhi = rawJets[0].phi();
    // 	leadEt = rawJets[0].Et();
    // 	vector<PseudoJet> LeadCons= rawJets[0].constituents();
    // 	leadNcons = LeadCons.size();
    // 	subPt = rawJets[0].pt();
    // 	subEta = rawJets[0].eta();
    // 	subPhi = rawJets[0].phi();
    // 	subEt = rawJets[0].Et();
    // 	vector<PseudoJet> SubCons= rawJets[1].constituents();
    // 	subNcons = SubCons.size();
    // 	// HTdijets->Fill();
    //   }
    // }

    // double djAxisPhi = ( rawJets[0].phi() + rawJets[1].phi() )/2;
    // Selector bgCircle = SelectorCircle(R);
    // Selector bgRapRange = SelectorRapRange( -0.6, 0.6 );
    // Selector bgSelector = bgCircle && bgRapRange;
    // double ghost_maxrap = 1.0;
    // AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
    // JetMedianBackgroundEstimator UE( bgSelector, jet_def, area_def);
    
    
    hVertex->Fill( Vx, Vy, Vz );                                        //  FILL HISTOGRAMS
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    
    // HTtree->Fill();                                                           //  FILL TREES
    // HTjets->Fill();
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  hPrimaryTracks->Write();                //  WRITE HISTOGRAMS
  hVertex->Write();
  hTowersPerEvent->Write();
  hTowersPerRun->Write();
  hPrimaryPerEvent->Write();
  hPrimaryPerRun->Write();
  hnPrimaryVSnTowers->Write();
  hPrimaryVsBBC->Write();
  hGlobalVsBBC->Write();
  hTowEt->Write();
  hTowEtEtaPhi->Write();
  hDCAvsPt->Write();
  hTowerMult->Write();

    
  // HTtree->Write();                            //  WRITE TREES
  // HTtowers->Write();
  // HTtracks->Write();
  // HTjets->Write();
  // HTdijets->Write();

  pAuFile->Close();
  
  
  return 0;
}
