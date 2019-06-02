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
    outFile = "out/pAuJets_HT.root";
    nEvents = 1000;
  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
  
  double dca, pt, pz, eta, phi, pmin1, pmax1, pmin2, pmax2;
    
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, sigma;
  int towID, nHitsPoss, nHitsFit, Charge, nCons;
  double DCA, towEt, towEta, towPhi, trEta, trPhi, trPz, trPt, jetPt, jetEta, jetPhi, jetEt;
  
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
  TH3D *hTowEtEtaPhi = new TH3D("hTowEtEtaPhi","Tower E_{T} vs. #eta vs. #phi;Tower E_{T} (GeV);Tower #eta;Tower #phi", 60,0,30.0, 100,-1.0,1.0, 128,-pi,pi );
  TH1D *hTowerMult = new TH1D("hTowerMult","Tower Multiplicity by ID;Tower ID", 2400,0,4800);
  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 50000,0,5000000, 150,0,150 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 50000,0,5000000, 300,0,3000 );
  TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 180,0,6.3, 100,-1.0,1.0);
  TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 180,0,6.3, 100,-1.0,1.0);
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 300,0,300, 40,0,200, 3500,0,3500000);
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 80,0,35, 100,0,1000);
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 80,0,35, 140,0,70);
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);
  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  JetDefinition bg_jet_def(kt_algorithm, R);     //  BACKGROUND ESTIMATION JET DEFINITION
  //  cambridge_algorithm??
  vector<PseudoJet> rawParticles, rawJets;
  int eID, rID;
  
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    rawParticles.clear();     rawJets.clear();     nCons.clear();     towID.clear();     Charge.clear();     towEt.clear();     towEta.clear();     towPhi.clear();
    trEta.clear();     trPhi.clear();     trPz.clear();     trPt.clear();     jetPt.clear();     jetEta.clear();     jetPhi.clear();     jetEt.clear();     DCA.clear();
    nHitsPoss.clear();     nHitsFit.clear();        //  clear vectors

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
    
    for ( int i=0; i<npt; ++i ) {                                         //  FILL PRIMARY TRACK DATA

      DCA = event->GetPrimaryTrack(i)->GetDCA();
      if ( DCA > dcaCut ) { nPrimary -= 1;   continue; }                   // track DCA cut
      trEta = (double) event->GetPrimaryTrack(i)->GetEta();
      if ( abs(trEta) > etaCut ) { nPrimary -= 1;   continue; }            // track eta cut
      trPt = (double) event->GetPrimaryTrack(i)->GetPt();
      trPz = (double) event->GetPrimaryTrack(i)->GetPz();
      trPhi = (double) event->GetPrimaryTrack(i)->GetPhi();

      nHitsFit = event->GetPrimaryTrack(i)->GetNOfFittedHits();
      nHitsPoss = event->GetPrimaryTrack(i)->GetNOfPossHits();

      hPrimaryTracks->Fill( trPt, trEta, trPhi );
      hDCAvsPt->Fill( trPt, DCA );
    }

    
    for ( int i=0; i<ntow; ++i ) {                                         //  FILL TOWER DATA
      towEt = event->GetTower(i)->GetEt();
      eta = event->GetTower(i)->GetEta();
      if ( abs(eta) > etaCut ) { nTowers -= 1;   continue; }            // tower eta cut
      towEta = eta;
      towPhi = event->GetTower(i)->GetPhi();
      towID = event->GetTower(i)->GetId();

      hTowEt->Fill( towID, towEt );
      hTowEtEtaPhi->Fill( towEt, towEta, towPhi );
      hTowerMult->Fill( towID );
    }

    hVertex->Fill( Vx, Vy, Vz );                                        //  FILL HISTOGRAMS
    hTowersPerEvent->Fill( nTowers );
    hTowersPerRun->Fill( RunID, nTowers );
    hPrimaryPerEvent->Fill( nPrimary );
    hPrimaryPerRun->Fill( RunID, nPrimary );
    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
    hPrimaryVsBBCE->Fill(nPrimary,BbcEastRate);
    hGlobalVsBBCE->Fill(nGlobal,BbcEastRate);

    //   JET-FINDING
    GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT>=2.0GeV    and    |eta|<1.0
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    nJets=0;
    
    for ( int i=0; i<rawJets.size(); ++i ) {                              //  FILL JET INFO
      jetPt = rawJets[i].pt();
      jetEta = rawJets[i].eta();
      jetPhi = rawJets[i].phi();
      jetEt = rawJets[i].Et();
      vector<PseudoJet> Cons= rawJets[i].constituents();
      nCons = Cons.size();
      nJets+=1;
    }


    if ( rawJets.size()<2) { continue; }
    
    double phi1 = rawJets[0].phi();     double phi2 = rawJets[1].phi();
    double dphi = fabs( phi1 - phi2 - pi );

    if ( dphi> R || rawJets[0].pt()<2.0 || rawJets[1].pt()<2.0 ) { continue; }
    else {                        //  CREATE DIJET PAIR  
      leadPt = rawJets[0].pt();      leadEta = rawJets[0].eta();      leadPhi = rawJets[0].phi();      leadEt = rawJets[0].Et();
      vector<PseudoJet> LeadCons= rawJets[0].constituents();      leadNcons = LeadCons.size();
      subPt = rawJets[1].pt();      subEta = rawJets[1].eta();      subPhi = rawJets[1].phi();      subEt = rawJets[1].Et();
      vector<PseudoJet> SubCons= rawJets[1].constituents();      subNcons = SubCons.size();
    }
    
    double tanPhi = tan(phi1) + tan(phi2);          tanPhi /= 2;
    double djAxisPhi = atan( tanPhi );       // qpi = quarter of pi
    pmin1 = djAxisPhi + qpi;           pmax1 = djAxisPhi + (3*qpi);
    pmin2 = djAxisPhi - (3*qpi);      pmax2 = djAxisPhi - qpi;
    Selector bgPhiRange = SelectorPhiRange( pmin1, pmax1 ) || SelectorPhiRange( pmin2, pmax2 );
    Selector bgSelector = bgPhiRange && SelectorAbsEtaMax( 1.0 );

    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition bg_area_def(active_area_explicit_ghosts, gAreaSpec);
    
    ClusterSequenceArea bgCluster( rawParticles, bg_jet_def, bg_area_def);
    JetMedianBackgroundEstimator UE( bgSelector, bgCluster );
    rho = UE.rho();
    sigma = UE.sigma();

    // cout << "rho = "<<rho<<"\t  sigma = "<<sigma<<endl<<endl;


    hLeadEtaPhi->Fill(leadPhi,leadEta);
    hSubEtaPhi->Fill(subPhi,subEta);
    hPt_UE_BBCE->Fill(leadPt,rho,BbcEastRate);
    hTowersVsRho->Fill(rho,nTowers);
    hLeadPtVsRho->Fill(rho,leadPt);


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
  hPrimaryVsBBCE->Write();
  hGlobalVsBBCE->Write();
  hLeadEtaPhi->Write();
  hSubEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hTowersVsRho->Write();
  hLeadPtVsRho->Write();

  pAuFile->Close();
  
  return 0;
}
