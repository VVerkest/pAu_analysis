//  pAuHTjetUE.cxx
//  Veronica Verkest		November 8, 2019

#include "pAuFunctions.hh"
#include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!!

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in pAuFunctions.hh

  int number_of_events;		string inFile, outFile;		TString name, title;
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  4 ) {    inFile = arguments[0];    outFile = arguments[1];    number_of_events = atoi(arguments[2].c_str()); }
  else if ( argc==1 ) { inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; outFile="out/UE/pAuHTjetUE.root"; number_of_events=-1; }
  else { cerr<< "incorrect number of command line arguments"; return 5000; }

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();

  TTree *HTjetTree = new TTree( "HTjetTree", "HT_JetTree" );

  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  int pval, jeval;
  
  //  Tree variables
  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nUEpart_chg, nUEpart_neu, nJetsAbove5, nHTtrig;
  double Vz, vpdVz, BbcAdcSumEast, leadPt, leadPtCorrected, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, chgEastRho_te,
    chgMidRho_te, chgWestRho_te, neuEastRho, neuMidRho, neuWestRho, leadArea, dPhiTrigLead, dRTrigLead, trigPhi;

  HTjetTree->Branch( "RunID", &RunID );
  HTjetTree->Branch( "EventID", &EventID );
  HTjetTree->Branch( "nTowers", &nTowers );
  HTjetTree->Branch( "nPrimary", &nPrimary );
  HTjetTree->Branch( "nGlobal", &nGlobal );
  HTjetTree->Branch( "nVertices", &nVertices );
  HTjetTree->Branch( "refMult", &refMult );
  HTjetTree->Branch( "gRefMult", &gRefMult );
  HTjetTree->Branch( "vpdVz", &vpdVz );
  HTjetTree->Branch( "Vz", &Vz );
  HTjetTree->Branch( "leadPt", &leadPt );
  HTjetTree->Branch( "leadPtCorrected", &leadPtCorrected );
  HTjetTree->Branch( "leadEta", &leadEta );
  HTjetTree->Branch( "BbcAdcSumEast", &BbcAdcSumEast );
  HTjetTree->Branch( "leadPhi", &leadPhi );
  HTjetTree->Branch( "chgEastRho", &chgEastRho );
  HTjetTree->Branch( "chgMidRho", &chgMidRho );
  HTjetTree->Branch( "chgWestRho_te", &chgWestRho_te );
  HTjetTree->Branch( "chgEastRho_te", &chgEastRho_te );
  HTjetTree->Branch( "chgMidRho_te", &chgMidRho_te );
  HTjetTree->Branch( "chgWestRho", &chgWestRho );
  HTjetTree->Branch( "neuEastRho", &neuEastRho );
  HTjetTree->Branch( "neuMidRho", &neuMidRho );
  HTjetTree->Branch( "neuWestRho", &neuWestRho );
  HTjetTree->Branch( "leadArea", &leadArea );
  HTjetTree->Branch( "nUEpart_chg", &nUEpart_chg );
  HTjetTree->Branch( "nUEpart_neu", &nUEpart_neu );
  HTjetTree->Branch( "nJetsAbove5", &nJetsAbove5 );
  HTjetTree->Branch( "nHTtrig", &nHTtrig );
  HTjetTree->Branch( "dPhiTrigLead", &dPhiTrigLead );
  HTjetTree->Branch( "dRTrigLead", &dRTrigLead );
  HTjetTree->Branch( "trigPhi", &trigPhi );

  // TH2D *hChgUePtEta[nPtBins];
  // TH2D *hNeuUePtEta[nPtBins];
  
  TH3D *hChgUE[nEtaBins];
  TH1D *hLeadPt[nEtaBins];
  const int xbins = 55;
  double xbinEdge[xbins+1];
  const int zbins = 20;
  double zbinEdge[zbins+1];  
  for (int i=0; i<xbins+1; ++i) {
    xbinEdge[i] = (double) i+4.0;
    if (i<zbins+1) { zbinEdge[i] = (double)(-10.0 + i)/10.0; }
  }
  // const int ybins = 15;
  // double ybinEdge[ybins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };

  const int ybins = 30;
  double ybinEdge[ybins+1];
  for (int i=0; i<=ybins; ++i) { ybinEdge[i] = 0.5*i; }

  TH2D *hLeadPtVsUEpT = new TH2D("hLeadPtVsUEpT",";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV)",xbins,xbinEdge,ybins,ybinEdge);
  TH1D *hLead = new TH1D("hLead", ";leading jet p_{T} (GeV)", 55,4.0,59.0);
  
  for (int e=0; e<nEtaBins; ++e) {
    
    name = "hChgUE" + etaBinName[e] + "Jet";
    title = etaBinString[e] + ";leading jet p_{T} (GeV);chg. UE part. p_{T} (GeV);chg. UE part. #eta";
    hChgUE[e] = new TH3D(name, title, xbins,xbinEdge,ybins,ybinEdge,zbins,zbinEdge);

    name = "hLeadPt" + etaBinName[e] + "Jet";
    title = etaBinString[e] + ";leading jet p_{T} (GeV)";
    hLeadPt[e] = new TH1D(name, title, 55,4.0,59.0);
    
  }

  for (int e=0; e<nEtaBins; ++e) {

    // name = "hChgUePtEta"; name += etaBinName[e];
    // title = "Charged Underlying Event #phi vs. #eta ("; title += etaBinString[e]; title += ") ;p_{T} (GeV);#eta";
    // hChgUePtEta[e] = new TH2D( name , title ,30,0,15,20,-1.0,1.0 );
    // name = "hNeuUePtEta"; name += etaBinName[e];
    // title = "Neutral Underlying Event #phi vs. #eta ("; title += etaBinString[e]; title += ") ;p_{T} (GeV);#eta";
    // hNeuUePtEta[e] = new TH2D( name, title, 30,0,15,20,-1.0,1.0 );
  }
  
  JetDefinition jet_def(antikt_algorithm, R);     //  JET DEFINITION
  Selector jetEtaSelector = SelectorAbsEtaMax( 1.0-R );
  Selector leadPtMinSelector = SelectorPtMin(leadJetMinPt);
  Selector leadJetSelector = jetEtaSelector;// && leadPtMinSelector;
  Selector ptMaxSelector = SelectorPtMax( 30.0 );

  PseudoJet leadJet, trigTowerPJ, towPJ;
  vector<PseudoJet> rawParticles, rawJets, allJets, chgParticles, neuParticles, UEparticles;
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TStarJetVectorContainer<TStarJetVector> * container;
  
  TChain* Chain = new TChain( "JetTree" );
  Chain->Add( inFile.c_str() );
  TStarJetPicoReader Reader;
  int numEvents = number_of_events;        // total events in HT: 152,007,032
  InitReader( Reader, Chain, numEvents );
  double deltaPhi, deltaR;       double trigTowEta, trigTowPhi;

  string efficFile = "src/trackeffic_loEA.root";
  // string efficFile = "src/trackeffic_hiEA.root";
  string UEcorrFile = "src/UEsubtractionPlots.root";

  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  while ( Reader.NextEvent() ) {
    Reader.PrintStatus(10);
    rawParticles.clear();    chgParticles.clear();    neuParticles.clear();    rawJets.clear();    allJets.clear();    //  CLEAR VECTORS
    event = Reader.GetEvent();
    header = event->GetHeader();
    container = Reader.GetOutputContainer();
    Vz = header->GetPrimaryVertexZ();
    if (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) { continue; }    //TEMPORARILY SKIPPING THESE RUNS
    if (header->GetRunId() == 16135031 || header->GetRunId() == 16135032) { continue; }
    if ( abs(Vz) > vzCut ) { continue; }
    if (!(header->HasTriggerId(500205) || header->HasTriggerId(500215))) {continue;}   //  ONLY SELECT HT TRIGGER EVENTS
    if ( header->GetBbcAdcSumEast() > 64000 ) { continue; }
    if ( header->GetBbcAdcSumEast() < 3559.12 ) { continue; }     //  neglect 90-100% event activity

    // //  HIGH EVENT ACTIVITY
    // if ( header->GetBbcAdcSumEast() < 26718.1 ) { continue; }  // LO: 3559.12-10126.1;  HI: 26718.1+// 

    //  LOW EVENT ACTIVITY
    if ( header->GetBbcAdcSumEast() > 10126.1 ) { continue; }  // LO: 3559.12-10126.1;  HI: 26718.1+
    
    TList *SelectedTowers = Reader.GetListOfSelectedTowers();
    nTowers = CountTowers( SelectedTowers );
		
    GatherParticles( container, rawParticles );
    GhostedAreaSpec gAreaSpec( 1.0, 1, 0.01 );
    AreaDefinition area_def(active_area_explicit_ghosts, gAreaSpec);
    ClusterSequenceArea jetCluster( rawParticles, jet_def, area_def); 
      
    leadJetSelector = jetEtaSelector; // && leadPtMinSelector;
    rawJets = sorted_by_pt( leadJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN >=10 GeV
    
    if ( rawJets.size()>0 ) {
      leadJet = rawJets[0]; 
      int trigTowId;
      TStarJetPicoTriggerInfo *trig;
      TStarJetPicoTower *tow, *triggerTower;
      double trigTowEt = 0.0;
      std::vector<int> trigTowers;
      for ( int i=0; i<event->GetTrigObjs()->GetEntries(); ++i ) {
	trig = (TStarJetPicoTriggerInfo *)event->GetTrigObj(i);
	if ( trig->isBHT2() && UseTriggerTower( trig->GetId()) ) { trigTowers.push_back( trig->GetId() ); }
      }
      sort(trigTowers.begin(), trigTowers.end());
      int nmatched = 0;
      for (int i=0; i<nTowers; ++i){				// loop throught selected towers in event
	tow = (TStarJetPicoTower*)SelectedTowers->At(i);
	if ( tow->GetEt()>=5.4 && count(trigTowers.begin(), trigTowers.end(), tow->GetId())) { // min 5.4 GeV tower and must be in list of HT towers
	  towPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
	  deltaPhi = fabs( leadJet.delta_phi_to( towPJ ) );
	  deltaR = leadJet.delta_R( towPJ );
	  if ( deltaR<=R || fabs(deltaPhi)>=(pi-R) ) {  // require trigger
	    if ( nmatched>0 && (tow->GetEt()<trigTowEt) ) { continue; }   // more than 1 trigger tower
	    else {							// first trigger tower
	      trigTowEt = tow->GetEt();
	      trigTowerPJ.reset_PtYPhiM( tow->GetEt(), tow->GetEta(), tow->GetPhi(), 0.0 ); //reset_PtYPhiM!!
	      nmatched += 1;
	    }
	  }
	  
	}
      }
      //if (nmatched==1) { // only accept events with 1 HT triggers
      if (nmatched>0) {
	nHTtrig = nmatched;
	trigPhi = trigTowerPJ.phi();
	dPhiTrigLead = fabs( leadJet.delta_phi_to( trigTowerPJ ) );
	dRTrigLead = leadJet.delta_R( trigTowerPJ );
	if ( dRTrigLead<=R || fabs(dPhiTrigLead)>=(pi-R) ) {
	  
	  Selector allJetSelector = SelectorPtMin(2.0) && ptMaxSelector && jetEtaSelector;
	  allJets = sorted_by_pt( allJetSelector( jetCluster.inclusive_jets() ) );     // EXTRACT ALL JETS IN >=5 GeV
	  
	  nJetsAbove5 = allJets.size();

	  vpdVz = header->GetvpdVz();
	  RunID = header->GetRunId();
	  EventID = Reader.GetNOfCurrentEvent();
	  nPrimary = header->GetNOfPrimaryTracks();
	  nGlobal = header->GetNGlobalTracks();
	  nVertices = header->GetNumberOfVertices();
	  refMult = header->GetReferenceMultiplicity();
	  gRefMult = header->GetGReferenceMultiplicity();
	  BbcAdcSumEast = header->GetBbcAdcSumEast();
	  leadPt = leadJet.pt();
	  leadPtCorrected = UEsubtraction( leadJet, UEcorrFile, BbcAdcSumEast );
	  leadEta = leadJet.eta();
	  leadPhi = leadJet.phi();
	  leadArea = leadJet.area();

	  //  UNDERLYING EVENT ESTIMATION
	  // GatherChargedUEwithEfficiency( leadJet, container, chgParticles, efficFile );
	  GatherChargedUE( leadJet, container, chgParticles );
	  GatherNeutralUE( leadJet, container, neuParticles );

	  nUEpart_chg = chgParticles.size();
	  nUEpart_neu = neuParticles.size();
    
	  chgEastSum = 0;  chgMidSum = 0;  chgWestSum = 0;  neuEastSum = 0;  neuMidSum = 0;  neuWestSum = 0;
	  CalculateRhoByChargeAndEta(chgParticles,neuParticles,chgEastSum,chgMidSum,chgWestSum,neuEastSum,neuMidSum,neuWestSum);
	  chgEastRho = chgEastSum/eastArea;
	  chgMidRho = chgMidSum/midArea;
	  chgWestRho = chgWestSum/westArea;
	  neuEastRho = neuEastSum/eastArea;
	  neuMidRho = neuMidSum/midArea;
	  neuWestRho = neuWestSum/westArea;
	  
	  pval = 99;    jeval = 99;

	  for ( int e=0; e<nEtaBins; ++e ) {
	    if ( leadEta >= etaLo[e]  &&  leadEta <= etaHi[e] ) { jeval = e; }
	  }
	  
	  if ( jeval==99 ) { cerr<<"UNABLE TO FIND PT OR ETA RANGE FOR LEAD JET"<<endl<<leadPt<<endl<<endl; }

	  
	  if (jeval==99) {continue;}
	  
	  for (int i=0; i<chgParticles.size(); ++i) {
	    hChgUE[jeval]->Fill( leadPt, chgParticles[i].pt(), chgParticles[i].eta() );
	    hLeadPtVsUEpT->Fill( leadPt, chgParticles[i].pt());
	    // hChgUE[jeval]->Fill( leadPtCorrected, chgParticles[i].pt(), chgParticles[i].eta() );
	    // hLeadPtVsUEpT->Fill( leadPtCorrected, chgParticles[i].pt());
	  }
	  // hLeadPt[jeval]->Fill( leadPtCorrected );
	  // hLead->Fill( leadPtCorrected );
	  hLeadPt[jeval]->Fill( leadPt );
	  hLead->Fill( leadPt );

	  chgParticles.clear(); // clear vector!

	  CalculateUEsubtractedChargedRho(chgParticles,chgEastSum,chgMidSum,chgWestSum);
	  chgEastRho_te = chgEastSum/eastArea;
	  chgMidRho_te = chgMidSum/midArea;
	  chgWestRho_te = chgWestSum/westArea;

	  
	  HTjetTree->Fill();
	}
      }
    }
  }  // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  END EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

  TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

  for ( int e=0; e<nEtaBins; ++e ) {  //  WRITE HISTOGRAMS & TREE
    hChgUE[e]->Write();
    hLeadPt[e]->Write();
  }
  HTjetTree->Write();  

  hLeadPtVsUEpT->Write();
  hLead->Write();
  
  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
