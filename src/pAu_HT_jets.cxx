//  pAu_HT_jets.cxx
//  Veronica Verkest		September 27, 2019

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace jetreader;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/HT/pAuJets.root"; number_of_events=-1; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicty;# of Towers", 200,0,200 );
  TH1D *hTriggerTowerId = new TH1D("hTriggerTowerId","HT Trigger Tower ID;Tower ID", 4800,0,4800);

  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 100,0,25, 200,0,200 );
  
  TH3D *hAllJetsPtEtaPhi = new TH3D( "hAllJetsPtEtaPhi", "Inclusive Jets p_{T}, #eta, #phi;Jet p_{T} (GeV);Jet #eta;Jet #phi", 400,0.0,100.0, 40,-1.0,1.0, 120,0.0,2*pi  );
  TH3D *hLeadJetPtRhoEta = new TH3D( "hLeadJetPtRhoEta", "Lead Jet p_{T}, #rho, #eta;Jet p_{T} (GeV);#rho;Jet #eta", 400,0.0,100.0, 100,0,25, 40,-1.0,1.0 );  
  TH3D *hLeadPtEtaPhi = new TH3D("hLeadPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;p_{T} (GeV);#eta;#phi", 280,0,70, 40,-1.0,1.0, 120,0,6.3);
  TH3D *hPt_UE_BBCE = new TH3D("hPt_UE_BBCE","UE vs. BBC East Rate;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC East Rate", 500,0,125, 50,0,25, 140,0,7000000 );
  TH3D *hPt_UE_BBCsumE = new TH3D("hPt_UE_BBCsumE","UE vs. BBC ADC East Sum;Lead Jet p_{T} (GeV);Underlying Event (GeV);BBC ADC East Sum", 500,0,125, 50,0,25, 160,0,80000 );
  TH3D *hBG = new TH3D("hBG","Background Particle #eta-#phi vs. p_{T};Background Particle p_{T}(GeV);Particle #eta;Particle #phi", 80,0,20, 40,-1,1, 120,0,2*pi );
  TH3D *hPartPtDEtaDPhi = new TH3D("hPartPtDEtaDPhi","Background Particle p_{T} vs. #Delta#eta vs. #Delta#phi;Particle p_{T} (GeV);#Delta#eta;#Delta#phi", 120,0,30, 80,-2.0,2.0, 120,-pi,pi );
  TH3D *hPartPtEtaPhi = new TH3D("hPartPtEtaPhi","Lead Jet p_{T} vs. #eta vs. #phi;Lead Jet p_{T} (GeV);Particle #eta;Particle #phi", 120,0,30, 40,-1.0,1.0, 120,0,2*pi );

  TH3D *hBackground[nPtBins][nEtaBins][nChgBins];
  TH2D *hRhoByEta[nPtBins][nEtaBins][nChgBins];

  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {
	name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];	title = "";		hRhoByEta[p][e][c] = new TH2D( name, title, 3,-1.5,1.5, 100,0,25 );
	hRhoByEta[p][e][c]->SetLineColor( color[c] );	hRhoByEta[p][e][c]->SetMarkerColor( color[c] );	hRhoByEta[p][e][c]->SetMarkerStyle( marker[c] );

	name = "hBG" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];
	title = ptBinString[p] + " " + etaBinString[e] + " " + BackgroundChargeBias[c] +" Background Particles;Background Particle p_{T}(GeV);Particle #phi;Particle #eta";
	hBackground[p][e][c] = new TH3D(name, title, 80,0,20, 120,0,2*pi, 40,-1,1 );
      }
    }
  }

  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

  double Vz, chgRho, neuRho, rho, eastRho, midRho, westRho;

  int nTowers;
  int nEvents =0;
  
  PseudoJet leadJet;  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles, BGparticles;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    rawParticles.clear();            chgParticles.clear();            neuParticles.clear();            rawJets.clear();            //  CLEAR VECTORS
    nTowers = 0;
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, event, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept HT Trigger events


    
    GatherParticles( container, rawParticles );
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS

    Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
    Selector ptMinSelector = SelectorPtMin( jetMinPt );
    Selector allJetSelector = jetEtaSelector && ptMinSelector;

    Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
    Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >10GeV

    if ( rawJets.size()>0 ) { leadJet = rawJets[0]; }
    else { continue; }


    
    TList *SelectedTowers = Reader.GetListOfSelectedTowers();		nTowers = CountTowers( SelectedTowers );		hTowersPerEvent->Fill(nTowers);
    
    hPrimaryPerEvent->Fill( header->GetNOfPrimaryTracks() );

    GatherChargedBG( leadJet, container, chgParticles);
    GatherNeutralBG( leadJet, container, neuParticles);
    
    double chgPtSum = 0;    double neuPtSum = 0;					//  BACKGROUND ESTIMATION 
    BackGroundEstimationAndPlots( chgParticles, neuParticles, leadJet, hPartPtDEtaDPhi, hPartPtEtaPhi, hBG, chgPtSum, neuPtSum );
    chgRho = chgPtSum / AREA;		neuRho = neuPtSum / AREA;			rho = (chgPtSum+neuPtSum) / AREA;

    hTowersVsRho->Fill( rho, nTowers );
    hLeadJetPtRhoEta->Fill( leadJet.pt(), rho, leadJet.eta() );
    hLeadPtEtaPhi->Fill( leadJet.pt(), leadJet.eta(), leadJet.phi() );
    hPt_UE_BBCE->Fill( leadJet.pt(), rho, header->GetBbcEastRate() );
    hPt_UE_BBCsumE->Fill( leadJet.pt(), rho, header->GetBbcAdcSumEast() );

    rawJets = sorted_by_pt( allJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >2GeV
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtEtaPhi->Fill( rawJets[i].pt(), rawJets[i].eta(), rawJets[i].phi() ); }
    
    for ( int p=0; p<3; ++p ) {  // * * * * * * * * * * * * * * * * * * PT LOOP * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      for ( int e=0; e<3; ++e ) {  // * * * * * * * * * * * * * * * * ETA LOOP * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	BGparticles.clear();		rawJets.clear();
	
	Selector ptRangeSelector = SelectorPtRange( ptLo[p], ptHi[p] );          //  JET pT RANGE
	Selector etaRangeSelector = SelectorEtaRange( etaLo[e], etaHi[e] );          //  JET eta RANGE
	Selector etaPtSelector = etaRangeSelector && jetEtaSelector && ptRangeSelector;
	
	rawJets = sorted_by_pt( etaPtSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS
	
	if ( rawJets.size()>0 ) { leadJet = rawJets[0]; }
	else { continue; }
	
	for ( int c=0; c<3; ++c ) {	  //  BACKGROUND ESTIMATION 
	  BGparticles.clear();

	  if ( BackgroundChargeBias[c]=="_chgBG" ) { GatherChargedBG( leadJet, container, BGparticles); }      //  Gather background particles 
	  else if ( BackgroundChargeBias[c]=="_neuBG" ) { GatherNeutralBG( leadJet, container, BGparticles); }
	  else if ( BackgroundChargeBias[c]=="_allBG" ) { GatherBackground( leadJet, container, BGparticles); }
	  else { cerr<<"BOO"<<endl; }

	  double eastSum = 0;	  double midSum = 0;	  double westSum = 0;
	  for (int i=0; i<BGparticles.size(); ++i) {

	    hBackground[p][e][c]->Fill( BGparticles[i].pt(), BGparticles[i].phi(), BGparticles[i].eta() );

	    if ( BGparticles[i].eta() >=etaLo[0] && BGparticles[i].eta() <= etaHi[0]  ) { eastSum+=BGparticles[i].pt(); }
	    else if ( BGparticles[i].eta() >=etaLo[1] && BGparticles[i].eta() <= etaHi[1]  ) { midSum+=BGparticles[i].pt(); }
	    else if ( BGparticles[i].eta() >=etaLo[2] && BGparticles[i].eta() <= etaHi[2]  ) { westSum+=BGparticles[i].pt(); }
	    else { cout<<BGparticles[i].eta()<<endl;        continue; }
	  }

	  eastRho = eastSum/eastArea;			midRho = midSum/midArea;			westRho = westSum/westArea;
	  
	  hRhoByEta[p][e][c]->Fill(-1.0, eastRho );
	  hRhoByEta[p][e][c]->Fill( 0.0, midRho );
	  hRhoByEta[p][e][c]->Fill( 1.0, westRho );

	}
      }
    }

    nEvents += 1;
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  //  ~ ~ ~ ~ ~ ~ ~ ~ WRITE HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~
  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {
	// hRhoByEta[p][e][c]->Scale(1.0/hRhoByEta[p][e][c]->Integral());
	hRhoByEta[p][e][c]->Write();
	hBackground[p][e][c]->Write();
      }
    }
  }

  hPrimaryPerEvent->Write();
  hTowersPerEvent->Write();
  hTriggerTowerId->Write();
  hTowersVsRho->Write();
  hAllJetsPtEtaPhi->Write();
  hLeadJetPtRhoEta->Write();
  hLeadPtEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hPt_UE_BBCsumE->Write();
  hPartPtDEtaDPhi->Write();
  hPartPtEtaPhi->Write();
  hBG->Write();

  pAuFile->Close();

  return 0;
}
