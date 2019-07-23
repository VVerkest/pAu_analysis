//  pAu_analysis.cxx
//  Veronica Verkest		May 14, 2019

//for pAu, need to ask if the event has the trigger. The trigger IDs are:
//HT2*BBCMB : 500205, 500215		JP2 : 500401, 500411
//BBCMB : 500008, 500018			VPDMB :  500904

#include "pAuFunctions.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // tracks and towers have eta and DCA cuts specified in pAuFunctions.hh

  int nEvents;		string inFile, outFile;	TString name, title;
  TString BackgroundChargeBias = "allBG";     // options: chgBG, neuBG, allBG
  TString JetChargeBias = "allJets";     // options: chgJets, neuJets, allJets
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  6 ) {    inFile = arguments[0];    outFile = arguments[1];    nEvents = atoi(arguments[2].c_str());    BackgroundChargeBias = arguments[3];    JetChargeBias = arguments[4]; }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/pAuDijets_HT.root"; nEvents=100000; BackgroundChargeBias="allBG"; JetChargeBias="allJets"; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  TH3D *hVertex = new TH3D( "hVertex", "Event Vertex;v_{x};v_{y};v_{z}", 60,-0.3,0.3, 60,-0.3,0.3, 160,-40,40 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Frequency;# of Towers", 140,0,700 );
  TH2D *hTowersPerRun = new TH2D("hTowersPerRun","Tower Frequency (per run);Run no.;# of Towers", 400,16124000,16164000, 140,0,700 );
  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 700,0,700 );
  TH2D *hPrimaryPerRun = new TH2D("hPrimaryPerRun","Primary Track Multiplicity (per run);Run no.;# of Primary", 4000,16124000,16164000, 40,0,200 );
  TH2D *hnPrimaryVSnTowers = new TH2D("hnPrimaryVSnTowers","# of Primary Tracks vs. # of Towers;# Towers;# Primary Tracks", 140,0,700, 40,0,200);
  TH2D *hPrimaryVsGlobal = new TH2D("hPrimaryVsGlobal","# Primary Tracks vs. # Global Tracks;# Global Tracks;# Primary Tracks", 150,0,3000, 150,0,150 );
  TH2D *hPrimaryVsBBC = new TH2D("hPrimaryVsBBC","# Primary Tracks vs. BBC Coincidence Rate;BBC Rate;# Primary Tracks", 3500,0,3500000, 40,0,200 );
  TH2D *hGlobalVsBBC = new TH2D("hGlobalVsBBC","# Global Tracks vs. BBC Coincidence Rate;BBC Rate;# Global Tracks", 3500,0,3500000, 150,0,3000 );
  TH2D *hPrimaryVsBBCE = new TH2D("hPrimaryVsBBCE","# Primary Tracks vs. BBC East Rate;BBC East Rate;# Primary Tracks", 140,0,7000000, 40,0,200 );
  TH2D *hGlobalVsBBCE = new TH2D("hGlobalVsBBCE","# Global Tracks vs. BBC East Rate;BBC East Rate;# Global Tracks", 80,0,80000, 150,0,3000 );
  TH2D *hPrimaryVsBBCsumE = new TH2D("hPrimaryVsBBCsumE","# Primary Tracks vs. BBC ADC East Sum;BBC ADC East Sum;# Primary Tracks", 160,0,80000, 40,0,200 );
  TH2D *hTowersVsBBCsumE = new TH2D("hTowersVsBBCsumE","# Towers vs. BBC ADC East Sum;BBC ADC East Sum;# Towers", 160,0,80000, 140,0,700 );
  TH2D *hChgVsNeuBG = new TH2D("hChgVsNeuBG","Charged vs. Neutral Background;Neutral #rho;Charged #rho",100,0,25,100,0,25);
  TH3D *hAllJetsPtEtaPhi = new TH3D( "hAllJetsPtEtaPhi", "Inclusive Jets p_{T}, #eta, #phi;Jet p_{T} (GeV);Jet #eta;Jet #phi", 400,0.0,100.0, 40,-1.0,1.0, 120,0.0,2*pi  );
  TH3D *hAllJetsPtRhoEta = new TH3D( "hAllJetsPtRhoEta", "Inclusive Jets p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hSubPtEtaPhi = new TH3D("hSubPtEtaPhi","Sub Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH2D *hPrimaryVsRho = new TH2D("hPrimaryVsRho","# Primary Tracks vs. Underlying Event;#rho (GeV);# Primary Tracks", 100,0,25, 40,0,200);
  TH2D *hGlobalVsRho = new TH2D("hGlobalVsRho","# Global Tracks vs. Underlying Event;#rho (GeV);# Global Tracks", 100,0,25, 150,0,3000 );
  TH3D *hPt_UE_RefMult = new TH3D("hPt_UE_RefMult","UE vs. Ref. Mult;Lead Jet p_{T} (GeV);Underlying Event (GeV);Ref. Mult.", 500,0,125, 50,0,25, 100,0,1000 );  
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 500,0,125, 50,0,25, 140,0,7000000 );
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 160,0,80000 );
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 100,0,25, 140,0,700 );
  TH2D *hLeadPtVsRho = new TH2D("hLeadPtVsRho","Lead Jet p_{T} vs UE;#rho (GeV);p_{T}^{lead} (GeV)", 140,0,35, 280,0,70);
  TH3D *hBG = new TH3D("hBG","Background Particle #eta-#phi vs. p_{T};Background Particle p_{T}(GeV);Particle #eta;Particle #phi", 80,0,20, 40,-1,1, 120,0,2*pi );
  TH3D *hPartPtDEtaDPhi = new TH3D("hPartPtDEtaDPhi","Background Particle p_{T} vs. #Delta#eta vs. #Delta#phi;Particle p_{T} (GeV);#Delta#eta;#Delta#phi", 120,0,30, 80,-2.0,2.0, 120,-pi,pi );
  TH3D *hPartPtEtaPhi = new TH3D("hPartPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;Lead Jet p_{T} (GeV);Particle #eta;Particle #phi", 120,0,30, 40,-1.0,1.0, 120,0,2*pi );
  TH3D *hAllPtEtaPhi = new TH3D("hAllPtEtaPhi","All Particles p_{T} vs. #eta vs. #phi;Particle p_{T} (GeV);Particle #eta;Particle #phi", 120,0,30, 40,-1.0,1.0, 120,0,2*pi );
  
  int EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, towID, nHitsPoss, nHitsFit, nCons;
  double Vx, Vy, Vz, BbcCoincidenceRate, BbcEastRate, BbcWestRate, BbcAdcSumEast, vpdVz, rho, chgRho, neuRho;

  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  
  vector<PseudoJet> rawParticles, chgParticles, neuParticles, rawJets, chgBG, neuBG;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  
  
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    chgParticles.clear();    neuParticles.clear();    rawParticles.clear();    rawJets.clear();       //  CLEAR VECTORS
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept JP2 Trigger events
    
    //  GATHERS ALL PARTICLES WITH    pT >= 0.2 GeV    and    |eta|<1.0
    if ( JetChargeBias == "allJets" ) { GatherParticles( container, rawParticles ); }
    else if ( JetChargeBias == "chgJets" ) { GatherCharged( container, rawParticles ); }
    else if ( JetChargeBias == "neuJets" ) { GatherNeutral( container, rawParticles ); }
    else { cerr<<"INCORRECT ARGUMENT FOR JetChargeBias:  CHOOSE FROM { allJets, chgJets, neuJets }"<<endl; break; }
    
    //   JET-FINDING
    JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

    Selector etaSelector = SelectorAbsEtaMax( 1.0-R );        //  CREATE JET SELECTORS
    Selector ptMinSelector = SelectorPtMin(jetMinPt);
    Selector eastJetSelector = SelectorEtaMax( 0.0 );
    Selector etaPtSelector = etaSelector && ptMinSelector;
    // Selector etaPtSelector = etaSelector && ptMinSelector && eastJetSelector;		//  REQUIRE LEAD & SUB JETS BOTH TO BE IN THE GOLD-GOING ETA-RANGE
    
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS > 2.0 GEV
    vector<PseudoJet> rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS
    
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtEtaPhi->Fill( rawJets[i].pt(), rawJets[i].eta(), rawJets[i].phi() ); }
    
    if ( rawJets.size()<2) { continue; }                                                       //  REQUIRE DIJET
    if ( fabs( fabs( rawJets[0].phi() - rawJets[1].phi() ) - pi ) > R || rawJets[0].pt()<10.0 || rawJets[1].pt()<2.0 ) { continue; }    // Require Lead>=10; Sublead>=2;  pi-R < |dphi| < pi+R
    
    for ( int i=0; i<rawParticles.size(); ++i ) { hAllPtEtaPhi->Fill( rawParticles[i].pt(), rawParticles[i].eta(), rawParticles[i].phi() ); }
    
    if ( BackgroundChargeBias == "allBG" || BackgroundChargeBias == "chgBG" ) {  GatherChargedBG( rawJets[0], container, chgParticles);  }    //  Gather background particles based on
    if ( BackgroundChargeBias == "allBG" || BackgroundChargeBias == "neuBG" ) {  GatherNeutralBG( rawJets[0], container, neuParticles);  }     //   BackgroundChargeBias (arg)

    double chgPtSum = 0;    double neuPtSum = 0;					//  BACKGROUND ESTIMATION 
    BackGroundEstimationAndPlots( chgParticles, neuParticles, rawJets[0], hPartPtDEtaDPhi, hPartPtEtaPhi, hBG, chgPtSum, neuPtSum );
    chgRho = chgPtSum / AREA;		neuRho = neuPtSum / AREA;			rho = (chgPtSum+neuPtSum) / AREA;

    TList *SelectedTowers = Reader.GetListOfSelectedTowers();	nTowers = CountTowers( SelectedTowers );

    EventID = Reader.GetNOfCurrentEvent();
    GetHeaderInfo( header, nGlobal, nVertices, refMult, nPrimary, BbcCoincidenceRate, vpdVz, BbcEastRate, BbcWestRate, BbcAdcSumEast );
    
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtRhoEta->Fill( rawJets[i].pt(), rho,rawJets[i].eta() ); }                      //  FILL HISTOGRAMS
    
    hLeadPtEtaPhi->Fill(rawJets[0].pt(),rawJets[0].eta(),rawJets[0].phi());	//hSubPtEtaPhi->Fill(rawJets[1].pt(),rawJets[1].eta(),rawJets[1].phi());
    hPrimaryVsRho->Fill( rho, nPrimary );							hGlobalVsRho->Fill(rho, nGlobal );						hPt_UE_BBCE->Fill(rawJets[0].pt(),rho,BbcEastRate);
    hTowersVsRho->Fill(rho,nTowers);	       						hLeadPtVsRho->Fill(rho,rawJets[0].pt());					hPt_UE_BBCsumE->Fill(rawJets[0].pt(),rho,BbcAdcSumEast);
    hTowersPerEvent->Fill( nTowers );						       	hPrimaryPerEvent->Fill( nPrimary );						hPt_UE_RefMult->Fill( rawJets[0].pt(), rho, refMult);
    hTowersPerRun->Fill( header->GetRunId(), nTowers );			hPrimaryPerRun->Fill( header->GetRunId(), nPrimary );		hnPrimaryVSnTowers->Fill( nTowers, nPrimary );
    hPrimaryVsBBC->Fill( BbcCoincidenceRate, nPrimary );			hPrimaryVsGlobal->Fill( nGlobal, nPrimary );				hGlobalVsBBC->Fill( BbcCoincidenceRate, nGlobal );
    hPrimaryVsBBCE->Fill(BbcEastRate,nPrimary);					hGlobalVsBBCE->Fill(BbcEastRate,nGlobal);				hChgVsNeuBG->Fill(neuRho,chgRho);
    hPrimaryVsBBCsumE->Fill(BbcAdcSumEast,nPrimary);			hTowersVsBBCsumE->Fill(BbcAdcSumEast,nTowers);
    hVertex->Fill( header->GetPrimaryVertexX(), header->GetPrimaryVertexY(), Vz );									
  }
  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  //  ~ ~ ~ ~ ~ ~ ~ ~ WRITE HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~
  hBG->Write();			      	hChgVsNeuBG->Write();		hAllJetsPtRhoEta->Write();
  hPartPtEtaPhi->Write();			hAllPtEtaPhi->Write(); 		hPartPtDEtaDPhi->Write();      	hVertex->Write();
  hTowersPerEvent->Write();		hTowersPerRun->Write();	hPrimaryPerEvent->Write();	hPrimaryPerRun->Write();
  hnPrimaryVSnTowers->Write();	hPrimaryVsBBC->Write();		hPrimaryVsGlobal->Write();	hGlobalVsBBC->Write();
  hPrimaryVsBBCE->Write();		hGlobalVsBBCE->Write();		hLeadPtEtaPhi->Write();		
  hAllJetsPtEtaPhi->Write();		hSubPtEtaPhi->Write();	       	hPt_UE_BBCE->Write();		hPt_UE_BBCsumE->Write();
  hTowersVsRho->Write();		hLeadPtVsRho->Write();		hPt_UE_RefMult->Write();	hPrimaryVsBBCsumE->Write();
  hTowersVsBBCsumE->Write();	hPrimaryVsRho->Write();		hGlobalVsRho->Write();
  
  pAuFile->Close();
  
  return 0;
}
