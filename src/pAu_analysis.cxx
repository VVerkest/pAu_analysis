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
    inFile = "production_pAu200_2015/MB/pAu_2015_200_MB*.root";
    outFile = "out/pAu.root";
    nEvents = -1;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in MB: 59388132
  InitReader( Reader, Chain, numEvents );
  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH3D *hPrimaryTracks = new TH3D( "hPrimaryTracks", "Primary Tracks: p_{T}, #eta, and #phi;p_{T} (GeV);#eta;#phi", 40,0,20, 40,-2,2, 16,-pi,pi );
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicity (per event);# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Multiplicity (per run);Run no.;# of Towers", 60,16120000,16160000, 700,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 60,16120000,16160000, 200,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;#Primary Tracks", 700,0,700, 200,0,200);
  TH2D *hDCAvsPt = new TH2D("hDCAvsPt","DCA vs. p_{T};p_{T} GeV;DCA (cm)", 70,0,3.5, 100,0.0,50.0);
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hTowEt = new TH2D("hTowEt","Tower E_{T} by ID;Tower ID;E_{T} (GeV)", 4800,0,4800, 40,0,20.0);
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta;Tower #phi", 100,0,20, 200,-1.0,1.0, 2*pi,-pi,pi );
  
  TTree *MBtree = new TTree( "MBTree", "MBtree" );
  TTree *MBtowers = new TTree( "MBTowers", "MBtowers" );
  TTree *MBtracks = new TTree( "MBTracks", "MBtracks" );
  TTree *MBjets = new TTree( "MBJets", "MBjets" );
  TTree *MBdijets = new TTree( "MBDijets", "MBdijets" );
  
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, towID, Charge, nHitsPoss, nHitsFit;
  double Vx, Vy, Vz, towEt, towEta, towPhi, trEta, trPhi, trPx, trPy, trPz, trPt, DCA, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz;

  MBtree->Branch("EventID", &EventID);            MBtree->Branch("RunID", &RunID);          MBtree->Branch("Vx", &Vx);
  MBtree->Branch("Vy", &Vy);                             MBtree->Branch("Vz", &Vz);                     MBtree->Branch("nTowers", &nTowers);
  MBtree->Branch("nPrimary", &nPrimary);         MBtree->Branch("nGlobal", &nGlobal);     MBtree->Branch("nVertices", &nVertices);
  MBtree->Branch("refMult", &refMult);              MBtree->Branch("gRefMult", &gRefMult);
  MBtree->Branch("BbcCoincidenceRate", &BbcCoincidenceRate);                    MBtree->Branch("BbcEastRate", &BbcEastRate);
  MBtree->Branch("BbcWestRate", &BbcWestRate);                                           MBtree->Branch("vpdVz", &vpdVz);
  
  MBtowers->Branch("EventID", &EventID);	  MBtowers->Branch("RunID", &RunID);  	  MBtowers->Branch("towEt", &towEt);
  MBtowers->Branch("towEta", &towEta);	          MBtowers->Branch("towPhi", &towPhi);	  MBtowers->Branch("towID", &towID);

  MBtracks->Branch("EventID", &EventID);	  MBtracks->Branch("RunID", &RunID);	  MBtracks->Branch("nHitsPoss", &nHitsPoss);
  MBtracks->Branch("nHitsFit", &nHitsFit);	  MBtracks->Branch("trPx", &trPx);   	  MBtracks->Branch("trPy", &trPy);
  MBtracks->Branch("trPz", &trPz);           	  MBtracks->Branch("trPt", &trPt);    	  MBtracks->Branch("trEta", &trEta);
  MBtracks->Branch("trPhi", &trPhi);        	  MBtracks->Branch("DCA", &DCA);

  int nJets;            vector<int> nCons;
  vector<double> jetPt, jetEta, jetPhi, jetEt;

  MBjets->Branch("EventID", &EventID);            MBjets->Branch("RunID", &RunID);            MBjets->Branch("nJets", &nJets);
  MBjets->Branch("jetPt", &jetPt);                      MBjets->Branch("jetEta", &jetEta);              MBjets->Branch("jetPhi", &jetPhi);
  MBjets->Branch("jetEt", &jetEt);                      MBjets->Branch("nCons", &nCons);            // MBjets->Branch("", &);

  double leadPt, leadEta, leadPhi, leadEt, leadNcons, subPt, subEta, subPhi, subEt, subNcons;
  
  MBdijets->Branch("EventID", &EventID);           MBdijets->Branch("RunID", &RunID);            MBdijets->Branch("leadPt", &leadPt);
  MBdijets->Branch("leadEta", &leadEta);            MBdijets->Branch("leadPhi", &leadPhi);         MBdijets->Branch("leadEt", &leadEt);
  MBdijets->Branch("leadNcons", &leadNcons);  MBdijets->Branch("subPt", &subPt);             MBdijets->Branch("subEta", &subEta);
  MBdijets->Branch("subPhi", &subPhi);              MBdijets->Branch("subEt", &subEt);             MBdijets->Branch("subNcons", &subNcons);
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);
  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  //JetMedianBackgroundEstimator UE( selector, jet_def, area_def);

  vector<PseudoJet> rawParticles, rawJets;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;

  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();  rawJets.clear();  nCons.clear();    jetPt.clear();    jetEta.clear();    jetPhi.clear();    jetEt.clear();  //  clear vectors

    Reader.PrintStatus(5);
    
    Vz = header->GetPrimaryVertexZ();
    if ( abs(Vz) > vzCut ) { continue; }

    eID = Reader.GetNOfCurrentEvent();          EventID = eID;
    rID = header->GetRunId();                        RunID = rID;

    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();

    int npt = header->GetNOfPrimaryTracks();      nPrimary = npt;
    int ntow = header->GetNOfTowers();               nTowers = ntow;
    Vx = header->GetPrimaryVertexX();    Vy = header->GetPrimaryVertexY();
    nGlobal = header->GetNGlobalTracks();                                    nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();                           gRefMult = header->GetGReferenceMultiplicity();
    BbcCoincidenceRate = header->GetBbcCoincidenceRate();        vpdVz = header->GetVpdVz();
    BbcEastRate = header->GetBbcEastRate();                                  BbcWestRate = header->GetBbcWestRate();

    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
    
    // From Nick's code:
    // TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    // evCuts->SetTriggerSelection( triggerString.c_str() );
    // evCuts->SetVertexZCut ( vertexZCut );
    // evCuts->SetMaxEventPtCut( eventPtCut );
    // evCuts->SetMaxEventEtCut( eventEtCut );
    // evCuts->SetVertexZDiffCut( vertexZDiffCut );
    // TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    // trackCuts->SetDCACut( DCACut );
    // trackCuts->SetMinNFitPointsCut( minFitPoints );
    // trackCuts->SetFitOverMaxPointsCut( minFitFrac );
    // trackCuts->SetMaxPtCut ( trackPtCut );
    // TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    // towerCuts->SetMaxEtCut( towerEtCut );
   
    
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

      nPrimary +=1;
      hPrimaryTracks->Fill( trPt, trEta, trPhi );
      hDCAvsPt->Fill( trPt, DCA );
      MBtracks->Fill();
    }

    
    for ( int i=0; i<ntow; ++i ) {                                         //  FILL TOWER DATA
      towEt = event->GetTower(i)->GetEt();
      towEta = event->GetTower(i)->GetEta();
      if ( abs(towEta) > etaCut ) { nTowers -= 1;   continue; }            // tower eta cut
      towPhi = event->GetTower(i)->GetPhi();
      towID = event->GetTower(i)->GetId();

      hTowEt->Fill( towID, towEt );
      hTowEtEtaPhi->Fill( towEt, towEta, towPhi );
      MBtowers->Fill();
    }



    //   JET-FINDING
    GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT>=2.0GeV    and    |eta|<1.0

    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

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
      if ( dphi< R ) {  
	leadPt = rawJets[0].pt();
	leadEta = rawJets[0].eta();
	leadPhi = rawJets[0].phi();
	leadEt = rawJets[0].Et();
	vector<PseudoJet> LeadCons= rawJets[0].constituents();
	leadNcons = SubCons.size();
	subPt = rawJets[0].pt();
	subEta = rawJets[0].eta();
	subPhi = rawJets[0].phi();
	subEt = rawJets[0].Et();
	vector<PseudoJet> SubCons= rawJets[1].constituents();
	subNcons = SubCons.size();
      }
    }


    
    hVertex->Fill( Vx, Vy, Vz );                                        //  FILL HISTOGRAMS
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    
    MBtree->Fill();                                                           //  FILL TREES
    MBjets->Fill();
    MBdijets->Fill();
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
  
  MBtree->Write();                            //  WRITE TREES
  MBtowers->Write();
  MBtracks->Write();
  MBjets->Write();
  MBdijets->Write();

  pAuFile->Close();
  
  return 0;
}
