//  pAu_HT_jets.cxx
//  Veronica Verkest		July 20, 2019

//HT2*BBCMB : 500205, 500215		JP2 : 500401, 500411
//BBCMB : 500008, 500018			VPDMB :  500904

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int nEvents;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    nEvents = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/pAuJets_HTJP2.root"; nEvents=100000; }
  else { cerr<< "incorrect number of command line arguments"; return -1; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TH3D *hBG = new TH3D("hBG","Background Particle p_{T} vs. #eta-#phi;Background Particle p_{T}(GeV);Particle #phi;Particle #eta", 80,0,20, 120,0,2*pi, 40,-1,1 );
  TH2D *hRhoByEta[nPtBins][nEtaBins][nChgBins];

  double eastArea = 1.4*(pi - 2);   // eta: [-1.0,-0.3]			(  etaMax - etaMin  ) X (  2*( pi-1 - 1 ) in phi  )
  double midArea = 1.2*(pi - 2);    //  eta: [-0.3,0.3]
  double westArea = 1.4*(pi - 2);    //  eta: [0.3,1.0]
  
  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {
	name = "hRho" + ptBinName[p] + etaBinName[e] + BackgroundChargeBias[c];	title = "";		hRhoByEta[p][e][c] = new TH2D( name, title, 7,0,7, 60,0,15 );
	hRhoByEta[p][e][c]->SetLineColor( color[c] );	hRhoByEta[p][e][c]->SetMarkerColor( color[c] );	hRhoByEta[p][e][c]->SetMarkerStyle( marker[c] );
      }
    }
  }

  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION

  vector<PseudoJet> rawParticles, chgParticles, BGparticles, neuParticles, rawJets;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = nEvents;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  double Vz, leadPt, leadEta, leadPhi, eastRho, midRho, westRho;		PseudoJet leadJet;


  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    rawParticles.clear();        //  CLEAR VECTORS
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept JP2 Trigger events

    GatherParticles( container, rawParticles );
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS

    for ( int p=0; p<3; ++p ) {  // * * * * * * * * * * * * * * * * * * PT LOOP * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      for ( int e=0; e<3; ++e ) {  // * * * * * * * * * * * * * * * * ETA LOOP * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	rawJets.clear();
	
	Selector ptRangeSelector = SelectorPtRange( ptLo[p], ptHi[p] );          //  JET pT RANGE
	Selector etaRangeSelector = SelectorEtaRange( etaLo[e], etaHi[e] );          //  JET eta RANGE
	Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
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

	    hBG->Fill( BGparticles[i].pt(), BGparticles[i].phi(), BGparticles[i].eta() );
	    
	    if ( etaLo[0] <= BGparticles[i].eta() <= etaHi[0]  ) { eastSum+=BGparticles[i].pt(); }
	    else if ( etaLo[1] < BGparticles[i].eta() < etaHi[1]  ) { midSum+=BGparticles[i].pt(); }
	    else if ( etaLo[2] <= BGparticles[i].eta() <= etaHi[2]  ) { westSum+=BGparticles[i].pt(); }
	    else { cerr<<BGparticles[i].eta()<<endl;        continue; }
	  }

	  eastRho = eastSum/eastArea;			midRho = midSum/midArea;			westRho = westSum/westArea;
	  
	  if ( eastRho != 0.0 ) { hRhoByEta[p][e][c]->Fill(1, eastRho ); }
	  if ( midRho != 0.0 ) { hRhoByEta[p][e][c]->Fill( 3, midRho ); }
	  if ( westRho != 0.0 ) { hRhoByEta[p][e][c]->Fill( 5, westRho ); }


	}
      }
    }

    



  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  for ( int p=0; p<3; ++p ) {	//  ~ ~ ~ ~ ~ ~ ~ ~ WRITE HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {
	// hRhoByEta[p][e][c]->Scale(1.0/hRhoByEta[p][e][c]->Integral());
	hRhoByEta[p][e][c]->Write();
      }
    }
  }
  hBG->Write();
  

  pAuFile->Close();

  return 0;
}
