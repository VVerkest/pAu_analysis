//  pAuHTjetUE.cxx
//  Veronica Verkest		November 8, 2019 ~ Updated May 2022

#include "pAuFunctions.hh"
#include "pAuParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!!
#include "THnSparse.h"

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {         // funcions and cuts specified in src/pAuFunctions and src/pAuParameters

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  TString name, title;

  int number_of_events;
  string inFile, outFile; //, EA_selection
  bool UE_subtraction, require_dijet;
    
/*
 ARGUMENTS:
 [0]: inFile                -> "production_pAu200_2015/HT/pAu_2015_200_HT*.root"
 [1]: outFile               -> "out/UE/pAuHTjetUE.root"
 [2]: number_of_events      -> -1 for all events
 [3]: UE_subtraction        -> boolean option for subtracting <rho>A from the jet
 [4]: require_dijet         -> boolean option for requiring a dijet
 [5]: trigger_selection     -> TO DO!
 */
    
  
  vector<string> arguments( argv+1, argv+argc );
  if ( argc ==  6 ) {
      inFile = arguments[0];
      outFile = arguments[1];
      number_of_events = atoi(arguments[2].c_str());
      if (arguments[3]=="true"||arguments[3]=="1") { UE_subtraction = true; }
      else if (arguments[3]=="false"||arguments[3]=="0") { UE_subtraction = false; }
      else { cerr<<"4th argument must be 'true', 'false', 1, or 0"<<endl; return 9999; }
      if (arguments[4]=="true"||arguments[4]=="1") { require_dijet = true; }
      else if (arguments[4]=="false"||arguments[4]=="0") { require_dijet = false; }
      else { cerr<<"5th argument must be 'true', 'false', 1, or 0"<<endl; return 9999; }
  }
  else if ( argc==1 || (argc==2&&arguments[0]=="default") ) {
      inFile = "production_pAu200_2015/HT/pAu_2015_200_HT*.root";
      outFile = "out/UE/pAuHTjetUE.root";
      number_of_events = -1;
      UE_subtraction =  false;
      require_dijet = false;
  }
  else { cerr<< "incorrect number of command line arguments"; return 5000; }
  
    const int NspJetDims = 5; // leadPt;leadEta;leadPhi;leadNcons;iBBCEsum
    const int NspJetBins[NspJetDims] = { 50, 20, 30, 30, 10 };
//    spJetBins
    
    const double bin_leadPt[NspJetBins[0]+1] = {4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54.};
    const double bin_leadEta[NspJetBins[1]+1] = {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    const double bin_leadPhi[NspJetBins[2]+1] = {-3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.09439, -1.88495, -1.67551, -1.46608, -1.25664, -1.0472, -0.837757, -0.628318, -0.418879, -0.209439, 0., 0.209439, 0.418879, 0.628318, 0.837757, 1.0472, 1.25664, 1.46608, 1.67551, 1.88495, 2.09439, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159};
    const double bin_leadNcons[NspJetBins[3]+1] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30.};
    const double bin_iBBCEsum[NspJetBins[4]+1] = {0., 3559.12, 6735.12, 10126.1, 13752.1, 17669.1, 21948.1, 26718.1, 32283.1, 39473.1, 64000.};
        
    
    THnSparseD *spJet = new THnSparseD("spJet","Leading Jet;leadPt;leadEta;leadPhi;leadNcons;iBBCEsum",NspJetDims,NspJetBins,NULL,NULL);
    spJet->SetBinEdges( 0, bin_leadPt );
    spJet->SetBinEdges( 1, bin_leadEta );
    spJet->SetBinEdges( 2, bin_leadPhi );
    spJet->SetBinEdges( 3, bin_leadNcons );
    spJet->SetBinEdges( 4, bin_iBBCEsum );

    const int NspUEdims = 6;
    const int NspUEbins[NspUEdims] = { 50, 15, 20, 30, 20, 10 };
    
    const double bin_UEpt[NspUEbins[1]+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };
    const double bin_UEeta[NspUEbins[2]+1] = {-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
    const double bin_UEphi[NspUEbins[3]+1] = {-3.14159, -2.93215, -2.72271, -2.51327, -2.30383, -2.09439, -1.88495, -1.67551, -1.46608, -1.25664, -1.0472, -0.837757, -0.628318, -0.418879, -0.209439, 0., 0.209439, 0.418879, 0.628318, 0.837757, 1.0472, 1.25664, 1.46608, 1.67551, 1.88495, 2.09439, 2.30383, 2.51327, 2.72271, 2.93215, 3.14159};
    
    THnSparseD *spUE = new THnSparseD("spUE","Underlying Event;leadPt;UEpt;UEeta;UEphi;leadEta;iBBCEsum",NspUEdims,NspUEbins,NULL,NULL);
    spUE->SetBinEdges( 0, bin_leadPt );
    spUE->SetBinEdges( 1, bin_UEpt );
    spUE->SetBinEdges( 2, bin_UEeta );
    spUE->SetBinEdges( 4, bin_UEphi );
    spUE->SetBinEdges( 4, bin_leadEta );
    spUE->SetBinEdges( 5, bin_iBBCEsum );

    
  //  Tree variables
//  int RunID, EventID, nTowers, nPrimary, nGlobal, nVertices, refMult, gRefMult, nUEpart_chg, nUEpart_neu, nJetsAbove5, nHTtrig;
//  double Vz, vpdVz, BbcAdcSumEast, leadPt, leadPtCorrected, leadEta, leadPhi, chgEastRho, chgMidRho, chgWestRho, chgEastRho_te,
//    chgMidRho_te, chgWestRho_te, neuEastRho, neuMidRho, neuWestRho, leadArea, dPhiTrigLead, dRTrigLead, trigPhi;


  
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
//  TStarJetPicoReader Reader;
//  InitReader( Reader, Chain, number_of_events );
//  double deltaPhi, deltaR, trigTowEta, trigTowPhi;
  //  double chgEastSum, chgMidSum, chgWestSum, neuEastSum, neuMidSum, neuWestSum;
  //  int pval, jeval;

//  string efficFile = "src/trackeffic.root";
//  string UEcorrFile = "src/UEsubtractionPlots.root";

    
    
    /*
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
    
    */
    
    

    TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");
    spJet->Write();
    spUE->Write();
    
//  for ( int e=0; e<nEtaBins; ++e ) {  //  WRITE HISTOGRAMS & TREE
//    hChgUE[e]->Write();
//    hLeadPt[e]->Write();
//  }
//  HTjetTree->Write();
//
//  hLeadPtVsUEpT->Write();
//  hLead->Write();
  
  pAuFile->Write();
  pAuFile->Close();

  return 0;
}
