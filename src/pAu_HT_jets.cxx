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

  int eID, rID, nEvents;                 string inFile, outFile;                 TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    nEvents = atoi(arguments[2].c_str());  }
  else if ( argc ==  1 ) {    inFile = "production_pAu200_2015/HT/pAu_2015_200_HT*.root";    outFile = "out/pAuJets_HT.root";    nEvents = 10000;  }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  
  double RunID, EventID, pmin1, pmax1, pmin2, pmax2, ptSum, dPhi, dEta, partPt, partEta, partPhi, partEt, deltaPhi;
  int nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons, towID, nHitsPoss, nHitsFit, Charge, nCons, partChg;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, BbcAdcSumEast, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, chgRho, neuRho;
  double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  
  vector<PseudoJet> rawParticles, chgParticles, neuParticles, rawJets, chgBG, neuBG;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;    TStarJetVectorContainer<TStarJetVector> * container;
  
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);        nJets=0;

    chgParticles.clear();    neuParticles.clear();    rawParticles.clear();    rawJets.clear();       //  CLEAR VECTORS
    
    event = Reader.GetEvent();            header = event->GetHeader();            container = Reader.GetOutputContainer();

    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }
    if (!(header->HasTriggerId(500401) || header->HasTriggerId(500411))) {continue;}   //  ONLY SELECT JP2 TRIGGER EVENTS
    Vz = header->GetPrimaryVertexZ();           if ( abs(Vz) > vzCut ) { continue; }
    if ( header->GetBbcAdcSumEast() > 64000 ) { continue; }

      
    //   JET-FINDING
    GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT >= 0.2 GeV    and    |eta|<1.0
    
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS > 2.0 GEV
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtEtaPhi->Fill( rawJets[i].pt(), rawJets[i].eta(), rawJets[i].phi() ); }
    
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

    for ( int i=0; i<rawParticles.size(); ++i ) { hAllPtEtaPhi->Fill( rawParticles[i].pt(), rawParticles[i].eta(), rawParticles[i].phi() ); }
    
    GatherChargedBG( rawJets[0], container, chgParticles);     GatherNeutralBG( rawJets[0], container, neuParticles);

    //  BACKGROUND ESTIMATION 
    double chgPtSum = 0;    double neuPtSum = 0;
    
    for (int i=0; i<chgParticles.size(); ++i) {
      partPt = chgParticles[i].pt();
      partEta = chgParticles[i].eta();
      partPhi = chgParticles[i].phi();
      partEt = chgParticles[i].Et();
      deltaPhi = chgParticles[i].delta_phi_to( rawJets[0] );
      partChg = chgParticles[i].user_index();
      dPhi = chgParticles[i].delta_phi_to( rawJets[0] );
      dEta = rawJets[0].eta() - chgParticles[i].eta();
      hPartPtDEtaDPhi->Fill( chgParticles[i].pt(), dEta, dPhi );
      hPartPtEtaPhi->Fill( leadPt, chgParticles[i].eta(), chgParticles[i].phi() );
      Charge = chgParticles[i].user_index();
      hCHARGED->Fill( leadPt, chgParticles[i].pt(), chgParticles[i].eta() );
      hBG->Fill( leadPt, chgParticles[i].eta(), chgParticles[i].phi() );
      chgPtSum+=chgParticles[i].pt();
      ptSum+=chgParticles[i].pt();
    }

    for (int i=0; i<neuParticles.size(); ++i) {
      partPt = neuParticles[i].pt();
      partEta = neuParticles[i].eta();
      partPhi = neuParticles[i].phi();
      partEt = neuParticles[i].Et();
      deltaPhi = neuParticles[i].delta_phi_to( rawJets[0] );
      partChg = neuParticles[i].user_index();
      dPhi = neuParticles[i].delta_phi_to( rawJets[0] );
      dEta = rawJets[0].eta() - neuParticles[i].eta();
      hPartPtDEtaDPhi->Fill( neuParticles[i].pt(), dEta, dPhi );
      hPartPtEtaPhi->Fill( neuParticles[i].pt(), neuParticles[i].eta(), neuParticles[i].phi() );
      Charge = neuParticles[i].user_index();
      hNEUTRAL->Fill( leadPt, neuParticles[i].pt(), neuParticles[i].eta() );
      hBG->Fill( leadPt, neuParticles[i].eta(), neuParticles[i].phi() );
      neuPtSum+=neuParticles[i].pt();
      ptSum+=neuParticles[i].pt();
    }
		      
    chgRho = chgPtSum / AREA;
    neuRho = neuPtSum / AREA;
    rho = (chgPtSum+neuPtSum) / AREA;


    
    eID = Reader.GetNOfCurrentEvent();          EventID = eID;
    rID = header->GetRunId();                        RunID = rID;
    nPrimary =  header->GetNOfPrimaryTracks();

    TStarJetPicoTower *tow;
    TList *SelectedTowers = reader.GetListOfSelectedTowers();
    nTowers = 0;
    for (int i=0; i<SelectedTowers->GetEntries(); ++i) {
      tow = (TStarJetPicoTower *) SelectedTowers->At(i);
      cout<<tow->GetEta()<<endl;
      if ( fabs(tow->GetEta())<=etaCut ) {nTowers+=1;}
    }
    
    Vx = header->GetPrimaryVertexX();                                           Vy = header->GetPrimaryVertexY();
    nGlobal = header->GetNGlobalTracks();                                    nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();                           gRefMult = header->GetGReferenceMultiplicity();
    BbcCoincidenceRate = header->GetBbcCoincidenceRate();        vpdVz = header->GetVpdVz();
    BbcEastRate = header->GetBbcEastRate();                                  BbcWestRate = header->GetBbcWestRate();
    BbcAdcSumEast = header->GetBbcAdcSumEast();



                                        //  FILL HISTOGRAMS
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtRhoEta->Fill( rawJets[i].pt(), rho,rawJets[i].eta() ); }

    hLeadPtEtaPhi->Fill(leadPt,leadPhi,leadEta);
    hLeadEtaPhi->Fill(leadPhi,leadEta);                                             hSubEtaPhi->Fill(subPhi,subEta);

    hPrimaryVsRho->Fill( rho, nPrimary );                                     hGlobalVsRho->Fill(rho, nGlobal );
    hPt_UE_BBCE->Fill(leadPt,rho,BbcEastRate);                            hPt_UE_BBCsumE->Fill(leadPt,rho,BbcAdcSumEast);
    hTowersVsRho->Fill(rho,nTowers);                                         hLeadPtVsRho->Fill(rho,leadPt);

    hVertex->Fill( Vx, Vy, Vz );
    hTowersPerEvent->Fill( nTowers );    hTowersPerRun->Fill( RunID, nTowers );    hPrimaryPerEvent->Fill( nPrimary );    hPt_UE_RefMult->Fill( leadPt, rho, refMult);
    hPrimaryPerRun->Fill( RunID, nPrimary );    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hPrimaryVsGlobal->Fill( nGlobal, nPrimary );    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );    hPrimaryVsBBCE->Fill(BbcEastRate,nPrimary);
    hGlobalVsBBCE->Fill(BbcEastRate,nGlobal);    hPrimaryVsBBCsumE->Fill(BbcAdcSumEast,nPrimary);    hTowersVsBBCsumE->Fill(BbcAdcSumEast,nTowers);

    hChgVsNeuBG->Fill(neuRho,chgRho);            
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");


  hBG->Write();                             //  WRITE HISTOGRAMS
  hCHARGED->Write();
  hNEUTRAL->Write();
  hChgVsNeuBG->Write();
  hAllJetsPtRhoEta->Write();
  hPartPtEtaPhi->Write();
  hAllPtEtaPhi->Write();
  hPartPtDEtaDPhi->Write();
  hVertex->Write();
  hTowersPerEvent->Write();
  hTowersPerRun->Write();
  hPrimaryPerEvent->Write();
  hPrimaryPerRun->Write();
  hnPrimaryVSnTowers->Write();
  hPrimaryVsBBC->Write();
  hPrimaryVsGlobal->Write();
  hGlobalVsBBC->Write();
  hPrimaryVsBBCE->Write();
  hGlobalVsBBCE->Write();
  hLeadPtEtaPhi->Write();
  hLeadEtaPhi->Write();
  hAllJetsPtEtaPhi->Write();
  hSubEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hPt_UE_BBCsumE->Write();
  hTowersVsRho->Write();
  hLeadPtVsRho->Write();
  hPt_UE_RefMult->Write();
  hPrimaryVsBBCsumE->Write();
  hTowersVsBBCsumE->Write();
  hPrimaryVsRho->Write();
  hGlobalVsRho->Write();
  hCHARGED->Write();
  hNEUTRAL->Write();
  
  pAuFile->Close();
  
  return 0;
}
