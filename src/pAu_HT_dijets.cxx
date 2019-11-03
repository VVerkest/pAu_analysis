//  pAu_HT_dijets.cxx
//  Veronica Verkest		October 4, 2019

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

  bool monojet = false;
  // Sort by Event Activity per Dave's definitions (https://drupal.star.bnl.gov/STAR/system/files/QM2019_Stewart_10.pdf)
  double BBCEmin = 0;        double BBCEmax = 100000;
  bool HiEA = false;		if ( HiEA==true ){ BBCEmin = 38000;    BBCEmax = 100000; }
  bool MidEA = true;	if ( MidEA==true ){ BBCEmin = 16000;    BBCEmax = 38000; }
  bool LoEA = false;		if ( LoEA==true ){ BBCEmin = 8000;    BBCEmax = 16000; }
  
  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TH1D *hPrimaryPerEvent = new TH1D("hPrimaryPerEvent","Primary Track Multiplicity (per event);# of Primary", 200,0,200 );
  TH1D *hTowersPerEvent = new TH1D("hTowersPerEvent","Tower Multiplicty;# of Towers", 200,0,200 );
  TH1D *hRecoAbsDeltaPhi = new TH1D("hRecoAbsDeltaPhi","|#phi_{lead} - #phi_{reco}|;|#Delta#phi|", 120,0.0,2*pi );
  TH1D *hRho = new TH1D("hRho","Underlying Event;#rho (GeV)",120,0,30);
  TH1D *hLeadPhi = new TH1D("hLeadPhi","Lead Jet #phi;#phi_{lead}",12,0,2*pi);
  TH1D *hRecoPhi = new TH1D("hRecoPhi","Reco Jet #phi;#phi_{lead}",12,0,2*pi);
  
  TH2D *hTowersVsRho = new TH2D("hTowersVsRho","# of Towers vs. UE;#rho (GeV);# of Towers", 100,0,25, 200,0,200 );
  TH2D *hRecoVsLeadPt = new TH2D("hRecoVsLeadPt","Recoil Jet p_{T} vs. Lead Jet p_{T};Lead Jet p_{T} (GeV);Reco Jet p_{T} (GeV)", 400,0.0,100.0, 400,0.0,100.0 );
  TH2D *hRecoVsLeadEta = new TH2D("hRecoVsLeadEta","Recoil Jet #eta vs. Lead Jet #eta;Lead Jet #eta;Reco Jet #eta", 40,-1.0,1.0, 40,-1.0,1.0 );
  
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

  double Vz, chgRho, neuRho, rho, eastRho[nChgBins], midRho[nChgBins], westRho[nChgBins];

  int nTowers;
  int nEvents =0;
  
  PseudoJet leadJet, recoJet;  vector<PseudoJet> rawParticles, rawJets, chgParticles, neuParticles, BGparticles, recoCandidates;
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );          Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;                                int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );

  int pval, eval;

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {

    Reader.PrintStatus(10);

    rawParticles.clear();            chgParticles.clear();            neuParticles.clear();            rawJets.clear();            //  CLEAR VECTORS
    nTowers = 0;            pval = 99;            eval = 99;
    
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    
    Vz = header->GetPrimaryVertexZ();
    if ( UseEvent( header, event, vzCut, Vz ) == false ) { continue; }   //  Skip events based on: Run#, vz cut, BBCEastSum;    only accept HT Trigger events

    double EA = header->GetBbcAdcSumEast();         // make selections on event activity
    if ( EA<BBCEmin || EA>BBCEmax ) { continue; }
    
    GatherParticles( container, rawParticles );
    ClusterSequence jetCluster( rawParticles, jet_def );           //  CLUSTER ALL JETS

    Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
    Selector ptMinSelector = SelectorPtMin( jetMinPt );
    Selector allJetSelector = jetEtaSelector && ptMinSelector;

    Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
    Selector leadJetSelector = jetEtaSelector && leadPtMinSelector;

    Selector ptMaxSelector = SelectorPtMax( 30.0 );
    leadJetSelector = leadPtMinSelector && ptMaxSelector && jetEtaSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN 10-30 GeV RANGE
    if ( rawJets.size()>0 ) { leadJet = rawJets[0]; }
    else { continue; }
    
    for ( int p=0; p<3; ++p ) {
      if ( leadJet.pt() >= ptLo[p]  &&  leadJet.pt() <= ptHi[p] ) { pval = p; }
    }
    for ( int e=0; e<3; ++e ) {
      if ( leadJet.eta() >= etaLo[e]  &&  leadJet.eta() <= etaHi[e] ) { eval = e; }
    }
    if ( pval==99 || eval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl; }
    
    //  DIJET BACKGROUND INFORMATION ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
    BGparticles.clear();		recoCandidates.clear();

    //   REQUIRE RECOIL JET TO HAVE AT LEAST HALF OF LEAD PT AND BE IN THE SAME ETA RANGE!
    Selector recoPtRangeSelector = SelectorPtRange( leadJet.pt()/2, ptHi[pval] );          //  JET pT RANGE    { 10-15, 15-20, 20-30 }
    Selector etaRangeSelector = SelectorEtaRange( etaLo[eval], etaHi[eval] );          //  JET eta RANGE

    Selector recoJetSelector = recoPtRangeSelector && etaRangeSelector && jetEtaSelector;
    
    recoCandidates = sorted_by_pt( recoJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT SELECTED JETS

    if ( monojet == false ) {
    
      bool hasReco = false;
    
      if ( recoCandidates.size()>0 ) {

	for ( int i=0; i<recoCandidates.size(); ++i ) {
	  double deltaPhi_abs = fabs( recoCandidates[i].phi() - leadJet.phi() );
	  if ( fabs( deltaPhi_abs - pi ) <= R ) {
	    recoJet = recoCandidates[i];
	    hRecoVsLeadPt->Fill( leadJet.pt(), recoJet.pt() );
	    hRecoAbsDeltaPhi->Fill( deltaPhi_abs );
	    hRecoVsLeadEta->Fill( leadJet.eta(), recoJet.eta() );
	    hRecoPhi->Fill( recoJet.phi() );

	    hasReco = true;
	  }
	  if ( hasReco == true ) { continue; }  // exit loop with highest-pt jet in recoil range
	}

      }
      else { continue; }
      if ( hasReco == false ) { continue; }

    }

    double eastSum[nChgBins], midSum[nChgBins], westSum[nChgBins], nMidJets[nChgBins], nEastJets[nChgBins], nWestJets[nChgBins];
    
    for ( int c=0; c<3; ++c ) {	  //  BACKGROUND ESTIMATION (lead and reco jet must be in same eta-range)
      BGparticles.clear();
      eastSum[c] = 0;	   	   midSum[c] = 0;	   	   westSum[c] = 0;
      eastRho[c] = 0;	   	   midRho[c] = 0;	   	   westRho[c] = 0;
      nMidJets[c] = 0;	   	   nEastJets[c] = 0;	   	   nWestJets[c] = 0;
      
      if ( BackgroundChargeBias[c]=="_chgBG" ) {      //  Gather background particles 
	GatherChargedBG( leadJet, container, BGparticles);
	for (int i=0; i<BGparticles.size(); ++i) {
	  hBackground[pval][eval][c]->Fill( BGparticles[i].pt(), BGparticles[i].phi(), BGparticles[i].eta() );

	  if ( BGparticles[i].eta() >=etaLo[0] && BGparticles[i].eta() <= etaHi[0]  ) { eastSum[c]+=BGparticles[i].pt();          ++nEastJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[1] && BGparticles[i].eta() <= etaHi[1]  ) { midSum[c]+=BGparticles[i].pt();           ++nMidJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[2] && BGparticles[i].eta() <= etaHi[2]  ) { westSum[c]+=BGparticles[i].pt();           ++nWestJets[c]; }
	  else { cout<<BGparticles[i].eta()<<endl;        continue; }
	  
	}
	
      }
      else if ( BackgroundChargeBias[c]=="_neuBG" ) {
	GatherNeutralBG( leadJet, container, BGparticles);
	for (int i=0; i<BGparticles.size(); ++i) {
	  hBackground[pval][eval][c]->Fill( BGparticles[i].pt(), BGparticles[i].phi(), BGparticles[i].eta() );

	  if ( BGparticles[i].eta() >=etaLo[0] && BGparticles[i].eta() <= etaHi[0]  ) { eastSum[c]+=BGparticles[i].pt();          ++nEastJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[1] && BGparticles[i].eta() <= etaHi[1]  ) { midSum[c]+=BGparticles[i].pt();           ++nMidJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[2] && BGparticles[i].eta() <= etaHi[2]  ) { westSum[c]+=BGparticles[i].pt();           ++nWestJets[c]; }
	  else { cout<<BGparticles[i].eta()<<endl;        continue; }
	  
	}
	
      }
      else if ( BackgroundChargeBias[c]=="_allBG" ) {
	GatherBackground( leadJet, container, BGparticles);
	for (int i=0; i<BGparticles.size(); ++i) {
	  hBackground[pval][eval][c]->Fill( BGparticles[i].pt(), BGparticles[i].phi(), BGparticles[i].eta() );
	  
	  if ( BGparticles[i].eta() >=etaLo[0] && BGparticles[i].eta() <= etaHi[0]  ) { eastSum[c]+=BGparticles[i].pt();          ++nEastJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[1] && BGparticles[i].eta() <= etaHi[1]  ) { midSum[c]+=BGparticles[i].pt();           ++nMidJets[c]; }
	  else if ( BGparticles[i].eta() >=etaLo[2] && BGparticles[i].eta() <= etaHi[2]  ) { westSum[c]+=BGparticles[i].pt();           ++nWestJets[c]; }
	  else { cout<<BGparticles[i].eta()<<endl;        continue; }
	  
	}
	
      }
      else { cerr<<"Error in background estimation"<<endl; }

      rho = 0;
      for ( int c=0; c<nChgBins; ++c ) {
	eastRho[c] = eastSum[c]/eastArea;
	midRho[c] = midSum[c]/midArea;
	westRho[c] = westSum[c]/westArea;
	  
	hRhoByEta[pval][eval][c]->Fill(-1.0, eastRho[c] );
	hRhoByEta[pval][eval][c]->Fill( 0.0, midRho[c] );
	hRhoByEta[pval][eval][c]->Fill( 1.0, westRho[c] );
      }

      rho += (eastSum[2] + midSum[2] + westSum[2])/AREA;
      
    } // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 


    TList *SelectedTowers = Reader.GetListOfSelectedTowers();		nTowers = CountTowers( SelectedTowers );		hTowersPerEvent->Fill(nTowers);
    hRho->Fill( rho );
    hTowersVsRho->Fill( rho, nTowers );
    hLeadJetPtRhoEta->Fill( leadJet.pt(), rho, leadJet.eta() );
    hLeadPtEtaPhi->Fill( leadJet.pt(), leadJet.eta(), leadJet.phi() );
    hPt_UE_BBCE->Fill( leadJet.pt(), rho, header->GetBbcEastRate() );
    hPt_UE_BBCsumE->Fill( leadJet.pt(), rho, header->GetBbcAdcSumEast() );    
    hPrimaryPerEvent->Fill( header->GetNOfPrimaryTracks() );
    hLeadPhi->Fill( leadJet.phi() );

    rawJets = sorted_by_pt( allJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS >2GeV
    for ( int i=0; i<rawJets.size(); ++i ) { hAllJetsPtEtaPhi->Fill( rawJets[i].pt(), rawJets[i].eta(), rawJets[i].phi() ); }
    

  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  
  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  //  ~ ~ ~ ~ ~ ~ ~ ~ WRITE HISTOGRAMS ~ ~ ~ ~ ~ ~ ~ ~
  for ( int p=0; p<3; ++p ) {
    for ( int e=0; e<3; ++e ) {
      for ( int c=0; c<3; ++c ) {
	hRhoByEta[p][e][c]->Write();
	hBackground[p][e][c]->Write();
      }
    }
  }

  hRecoVsLeadPt->Write();
  hRecoVsLeadEta->Write();
  hRecoAbsDeltaPhi->Write();
  hPrimaryPerEvent->Write();
  hTowersPerEvent->Write();
  hTowersVsRho->Write();
  hAllJetsPtEtaPhi->Write();
  hLeadJetPtRhoEta->Write();
  hLeadPtEtaPhi->Write();
  hPt_UE_BBCE->Write();
  hPt_UE_BBCsumE->Write();
  hPartPtDEtaDPhi->Write();
  hPartPtEtaPhi->Write();
  hBG->Write();
  hRho->Write();
  hLeadPhi->Write();
  hRecoPhi->Write();

  pAuFile->Close();

  return 0;
}
