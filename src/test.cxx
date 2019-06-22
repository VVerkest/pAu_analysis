#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;


int main(){
  
  int eID, rID, nEvents;                 string inFile, outFile;                 TString name, title;
  inFile = "production_pAu200_2015/HT/pAu_2015_200_HT*.root";    outFile = "out/pAuJets_test.root";    nEvents = 10000;

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TH3D *hBG = new TH3D("hBG","Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hCHARGED = new TH3D("hCHARGED","CHARGED: Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hNEUTRAL = new TH3D("hNEUTRAL","NEUTRAL: Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hPartPtEtaDPhi = new TH3D("hPartPtEtaDPhi","Background Particle p_{T} vs. #eta vs. #Delta#phi;Particle p_{T} (GeV);#eta;#Delta#phi", 200,0,50, 220,-1.1,1.1, 252,-pi,pi );
  
  double pmin1, pmax1, pmin2, pmax2, ptSum;
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons, towID, nHitsPoss, nHitsFit, Charge, dPhi, nCons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, BbcAdcSumEast, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, chgRho, neuRho, rho_avg;
  vector<double> partPt, partEta, partPhi, partEt, deltaPhi;         vector<int> partChg;

  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  TTree *sp[1];    // currently using all eta [-1,1]
  for ( int i=0; i<1; ++i ) {
    name = "sp" + etaBinName[i];          title = "Selected Particles: " + etaBinString[i];          sp[i] = new TTree( name, title );
    sp[i]->Branch("partPt",&partPt);        sp[i]->Branch("partEta",&partEta);        sp[i]->Branch("partPhi",&partPhi);
    sp[i]->Branch("partEt",&partEt);        sp[i]->Branch("chgRho",&chgRho);        sp[i]->Branch("neuRho",&neuRho);        sp[i]->Branch("rho",&rho);
    sp[i]->Branch("deltaPhi",&deltaPhi);    sp[i]->Branch("leadPt",&leadPt);        sp[i]->Branch("subPt",&subPt);
    sp[i]->Branch("nTowers",&nTowers);        sp[i]->Branch("partChg",&partChg);
  }
  
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

    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS  }
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }   //                                                              }   Move these to function
    if (!(header->HasTriggerId(500401) || header->HasTriggerId(500411))) {continue;}   //  ONLY SELECT JP2 TRIGGER EVENTS               }
    Vz = header->GetPrimaryVertexZ();           if ( abs(Vz) > vzCut ) { continue; }              //  VZ CUT                                                         }

    // GatherCharged( container, chgParticles);
    // GatherNeutral( container, neuParticles);
    for (int i=0; i<container->GetEntries(); ++i) {
      sv = container->Get(i);
      
      if ( abs(sv->Eta()) > etaCut  ||  sv->Pt() < partMinPt)      { continue; }  // removes particles with |eta|>1  OR   with pt<0.2GeV
      rawParticles.push_back( sv );
      if (sv->GetCharge()==0) {	neuParticles.push_back( sv ); }
      else if (sv->GetCharge()==1 || sv->GetCharge()==-1) { chgParticles.push_back( sv ); }
      else { cerr<<"charge error"<<endl; }

    }


    //   JET-FINDING
    //  GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT >= 0.2 GeV    and    |eta|<1.0
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS > 2.0 GEV
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    if ( rawJets.size()<2) { continue; }                                                       //  REQUIRE DIJET
    double phi1 = rawJets[0].phi();     double phi2 = rawJets[1].phi();    double dphi = fabs( fabs( phi1 - phi2 ) - pi );
    if ( dphi> R || rawJets[0].pt()<10.0 || rawJets[1].pt()<2.0 ) { continue; }
    else {                        //  CREATE DIJET PAIR  
      leadPt = rawJets[0].pt();      leadEta = rawJets[0].eta();      leadPhi = rawJets[0].phi();      leadEt = rawJets[0].Et();
      vector<PseudoJet> LeadCons= rawJets[0].constituents();      leadNcons = LeadCons.size();
      subPt = rawJets[1].pt();      subEta = rawJets[1].eta();      subPhi = rawJets[1].phi();      subEt = rawJets[1].Et();
      vector<PseudoJet> SubCons= rawJets[1].constituents();      subNcons = SubCons.size();
    }

    // JetDefinition bg_jet_def(kt_algorithm, 0.0001);     //  BACKGROUND ESTIMATION JET DEFINITION
    double phiRef1 = phi1;    double phiRef2 = phi1 - pi;
    Selector BGphiSelector = !( SelectorPhiRange( phi1-qpi, phi1+qpi ) || SelectorPhiRange( phi2-qpi, phi2+qpi ) );  //  EXCLUSIVE OR
    Selector bgSelector = BGphiSelector && etaSelector;

    chgBG = bgSelector( chgParticles );
    neuBG = bgSelector( neuParticles );

    double chgPtSum = 0;    double neuPtSum = 0;
    for (int i=0; i<chgBG.size(); ++i) {
      partPt.push_back( chgBG[i].pt() );
      partEta.push_back( chgBG[i].eta() );
      partPhi.push_back( chgBG[i].phi() );
      partEt.push_back( chgBG[i].Et() );
      deltaPhi.push_back( chgBG[i].delta_phi_to( rawJets[0] ) );
      partChg.push_back( chgBG[i].user_index() );
      dPhi = chgBG[i].delta_phi_to( rawJets[0] );
      hPartPtEtaDPhi->Fill( chgBG[i].pt(), chgBG[i].eta(), dPhi );
      Charge = chgBG[i].user_index();
      hCHARGED->Fill( leadPt, chgBG[i].pt(), chgBG[i].eta() );
      hBG->Fill( leadPt, chgBG[i].pt(), chgBG[i].eta() );
      chgPtSum+=chgBG[i].pt();
      ptSum+=chgBG[i].pt();
    }

    for (int i=0; i<neuBG.size(); ++i) {
      partPt.push_back( neuBG[i].pt() );
      partEta.push_back( neuBG[i].eta() );
      partPhi.push_back( neuBG[i].phi() );
      partEt.push_back( neuBG[i].Et() );
      deltaPhi.push_back( neuBG[i].delta_phi_to( rawJets[0] ) );
      partChg.push_back( neuBG[i].user_index() );
      dPhi = neuBG[i].delta_phi_to( rawJets[0] );
      hPartPtEtaDPhi->Fill( neuBG[i].pt(), neuBG[i].eta(), dPhi );
      Charge = neuBG[i].user_index();
      hNEUTRAL->Fill( leadPt, neuBG[i].pt(), neuBG[i].eta() );
      hBG->Fill( leadPt, neuBG[i].pt(), neuBG[i].eta() );
      neuPtSum+=neuBG[i].pt();
      ptSum+=neuBG[i].pt();
    }

    chgRho = chgPtSum / pi;
    neuRho = neuPtSum / pi;
    rho = (chgPtSum+neuPtSum) / pi;

    sp[0]->Fill();
    
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  for ( int e=0; e<1; ++e ) { sp[e]->Write(); }                                     //  WRITE PARTICLE TREE

  hBG->Write();
  hCHARGED->Write();
  hNEUTRAL->Write();
  hPartPtEtaDPhi->Write();
  
  pAuFile->Close();

  return 0;
}
