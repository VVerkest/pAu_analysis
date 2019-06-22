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
  
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Frequency;# of Towers", 700,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Frequency (per run);Run no.;# of Towers", 400,16124000,16164000, 140,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 4000,16124000,16164000, 40,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;# Primary Tracks", 140,0,700, 40,0,200);
  TH2D *hPrimaryVsGlobal = new TH2D("hPrimaryVsGlobal","# Primary Tracks vs. # Global Tracks;# Global Tracks;# Primary Tracks", 150,0,3000, 150,0,150 );
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 3500,0,3500000, 40,0,200 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 3500,0,3500000, 150,0,3000 );
  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 140,0,7000000, 150,0,3000 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 80,0,80000, 150,0,3000 );
  TH2D *hPrimaryVsBBCsumE = new TH2D("hPrimaryVsBBCsumE","# Primary Tracks vs. BBC ADC East Sum;BBC ADC East Sum;# Primary Tracks", 160,0,80000, 40,0,200 );
  TH2D *hTowersVsBBCsumE = new TH2D("hTowersVsBBCsumE","# Towers vs. BBC ADC East Sum;BBC ADC East Sum;# Towers", 160,0,80000, 140,0,700 );
  TH2D *hLeadEtaPhi = new TH2D("hLeadEtaPhi","Lead Jet #eta vs. #phi;#phi;#eta", 252,0,6.3, 160,-0.8,0.8);
  TH2D *hSubEtaPhi = new TH2D("hSubEtaPhi","Sub Jet #eta vs. #phi;#phi;#eta", 252,0,6.3, 160,-0.8,0.8);
  TH2D *hPrimaryVsRho = new TH2D("hPrimaryVsRho","# Primary Tracks vs. Underlying Event;#rho (GeV);# Primary Tracks", 100,0,25, 40,0,200);
  TH2D *hGlobalVsRho = new TH2D("hGlobalVsRho","# Global Tracks vs. Underlying Event;#rho (GeV);# Global Tracks", 100,0,25, 150,0,3000 );
  TH3D *hPt_UE_RefMult = new TH3D("hPt_UE_RefMult","UE vs. Ref. Mult;Lead Jet p_{T} (GeV);Underlying Event (GeV);Ref. Mult.", 500,0,125, 50,0,25, 300,0,150000 );  
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 500,0,125, 50,0,25, 140,0,7000000 );
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 160,0,80000 );
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 100,0,25, 140,0,700 );
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 140,0,35, 280,0,70);
  TH3D *hBG = new TH3D("hBG","Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hCHARGED = new TH3D("hCHARGED","CHARGED: Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hNEUTRAL = new TH3D("hNEUTRAL","NEUTRAL: Background Patricle #eta vs. p_{T};Lead Jet p_{T}(GeV);Particle p_{T} (GeV);Particle #eta", 280,0,70, 200,0,50, 220,-1.1,1.1 );
  TH3D *hPartPtEtaDPhi = new TH3D("hPartPtEtaDPhi","Background Particle p_{T} vs. #eta vs. #Delta#phi;Particle p_{T} (GeV);#eta;#Delta#phi", 200,0,50, 220,-1.1,1.1, 252,-pi,pi );
  
  double pmin1, pmax1, pmin2, pmax2, ptSum, dPhi;
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nJets, leadNcons, subNcons, towID, nHitsPoss, nHitsFit, Charge, nCons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, BbcAdcSumEast, vpdVz,  leadPt, leadEta, leadPhi, leadEt, subPt, subEta, subPhi, subEt, rho, chgRho, neuRho;
  vector<double> partPt, partEta, partPhi, partEt, deltaPhi;         vector<int> partChg;

  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  TTree *bgTree= new TTree( "bgTree", "Background Estimation Tree" );
  bgTree->Branch("partPt",&partPt);        bgTree->Branch("partEta",&partEta);        bgTree->Branch("partPhi",&partPhi);
  bgTree->Branch("partEt",&partEt);        bgTree->Branch("chgRho",&chgRho);        bgTree->Branch("neuRho",&neuRho);        bgTree->Branch("rho",&rho);
  bgTree->Branch("deltaPhi",&deltaPhi);    bgTree->Branch("leadPt",&leadPt);        bgTree->Branch("subPt",&subPt);
  bgTree->Branch("nTowers",&nTowers);        bgTree->Branch("partChg",&partChg);
  
  //  CREATE JET SELECTOR
  Selector etaSelector = SelectorAbsEtaMax( 1.0-R );    Selector ptMinSelector = SelectorPtMin(jetMinPt);  Selector etaPtSelector = etaSelector && ptMinSelector;
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  JetDefinition bg_jet_def(kt_algorithm, 0.0001);     //  BACKGROUND ESTIMATION JET DEFINITION
  
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
      
    //   JET-FINDING
    // GatherParticles( container, rawParticles);     //  GATHERS ALL PARTICLES WITH    pT >= 0.2 GeV    and    |eta|<1.0
    for (int i=0; i<container->GetEntries(); ++i) {
      sv = container->Get(i);
      
      if ( abs(sv->Eta()) > etaCut  ||  sv->Pt() < partMinPt)      { continue; }  // removes particles with |eta|>1  OR   with pt<0.2GeV
      PseudoJet current = PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );
      rawParticles.push_back( current );
      if (current.user_index()==0) {	neuParticles.push_back( current ); }
      else if (current.user_index()==1 || current.user_index()==-1) { chgParticles.push_back( current ); }
      else { cerr<<"charge error"<<endl; }
    }
    
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS > 2.0 GEV
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

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


    

double phiRef1 = phi1;    double phiRef2 = phi1 - pi;     //  BACKGROUND ESTIMATION 
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

    bgTree->Fill();

    
    eID = Reader.GetNOfCurrentEvent();          EventID = eID;
    rID = header->GetRunId();                        RunID = rID;
    nPrimary =  header->GetNOfPrimaryTracks();                           nTowers = header->GetNOfTowers();
    Vx = header->GetPrimaryVertexX();                                           Vy = header->GetPrimaryVertexY();
    nGlobal = header->GetNGlobalTracks();                                    nVertices = header->GetNumberOfVertices();
    refMult = header->GetReferenceMultiplicity();                           gRefMult = header->GetGReferenceMultiplicity();
    BbcCoincidenceRate = header->GetBbcCoincidenceRate();        vpdVz = header->GetVpdVz();
    BbcEastRate = header->GetBbcEastRate();                                  BbcWestRate = header->GetBbcWestRate();
    BbcAdcSumEast = header->GetBbcAdcSumEast();

    hLeadEtaPhi->Fill(leadPhi,leadEta);                                             hSubEtaPhi->Fill(subPhi,subEta);

    hPrimaryVsRho->Fill( rho, nPrimary );                                     hGlobalVsRho->Fill(rho, nGlobal );
    hPt_UE_BBCE->Fill(leadPt,rho,BbcEastRate);                            hPt_UE_BBCsumE->Fill(leadPt,rho,BbcAdcSumEast);
    hTowersVsRho->Fill(rho,nTowers);                                         hLeadPtVsRho->Fill(rho,leadPt);

    hVertex->Fill( Vx, Vy, Vz );                                        //  FILL HISTOGRAMS
    hTowersPerEvent->Fill( nTowers );    hTowersPerRun->Fill( RunID, nTowers );    hPrimaryPerEvent->Fill( nPrimary );    hPt_UE_RefMult->Fill( leadPt, rho, refMult);
    hPrimaryPerRun->Fill( RunID, nPrimary );    hnPrimaryVSnTowers->Fill( nTowers, nPrimary );    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );
    hPrimaryVsGlobal->Fill( nGlobal, nPrimary );    hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );    hPrimaryVsBBCE->Fill(BbcEastRate,nPrimary);
    hGlobalVsBBCE->Fill(BbcEastRate,nGlobal);    hPrimaryVsBBCsumE->Fill(BbcAdcSumEast,nPrimary);    hTowersVsBBCsumE->Fill(BbcAdcSumEast,nTowers);

  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  bgTree->Write();                                     //  WRITE PARTICLE TREE

  hBG->Write();                             //  WRITE HISTOGRAMS
  hCHARGED->Write();
  hNEUTRAL->Write();
  hPartPtEtaDPhi->Write();
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
  hLeadEtaPhi->Write();
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
  hPartPtEtaDPhi->Write();
  
  pAuFile->Close();
  
  return 0;
}
