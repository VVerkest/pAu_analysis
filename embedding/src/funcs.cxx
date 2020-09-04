//!  functions.cxx
//!  Isaac Mooney, WSU - June 2019

#include "params.hh"
#include "funcs.hh"
//!#include "ktTrackEff.hh"

typedef fastjet::contrib::SoftDrop SD;

namespace Analysis {

  int CountTowers( TList *selectedtowers ) {

    TStarJetPicoTower *tow;
    int n_towers = 0;
    for (int i=0; i<selectedtowers->GetEntries(); ++i) {
      tow = (TStarJetPicoTower *) selectedtowers->At(i);
      if ( fabs(tow->GetEta())<=maxTrackEta ) { n_towers+=1; }
    }
    return n_towers;
  }

  
  //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from whence they came
  bool DiscardpAuEmbedEvent(const TString Filename, const std::vector<fastjet::PseudoJet> pJets, const std::vector<fastjet::PseudoJet> dJets) {
    bool bad_event = 0;

    std::string pt_hat[9] = { "pt-hat57", "pt-hat79", "pt-hat911", "pt-hat1115", "pt-hat1525", "pt-hat2535", "pt-hat3545", "pt-hat4555", "pt-hat5565" };
    double maxPtVal[9] = { 7.0, 9.0, 11.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0 };

    std::string name = (std::string) Filename;
    bool binFound = 0;

    // std::cout<<name<<" ==> "<<Filename<<std::endl;
    
    for ( int i=0; i<9; ++i ) {

      //std::cout<<pt_hat[i]<<"\t name.find(pt_hat[i]) = "<<name.find(pt_hat[i])<<std::endl;
      
      if ( name.find(pt_hat[i]) != std::string::npos ) {  // loop over and find pythia pt bin from file name

	binFound = 1;
	
	if (pJets.size() != 0) {
	  if (pJets[0].pt() > 2*maxPtVal[i]) { bad_event = 1; }  // if lead jet is over twice the upper range of pt bin, discard the event
	}

	if (dJets.size() != 0) {
	  if (dJets[0].pt() > 2*maxPtVal[i]) { bad_event = 1; }  // if lead jet is over twice the upper range of pt bin, discard the event
	}
	
      }
    }
    //std::cout<<"\n"<<std::endl;
    if ( binFound==0 ) { std::cerr<<"Pythia pT bin not found for this file: "<<Filename<<std::endl; }   
    return bad_event;
  }


  void DrawText(const char *text, float xp, float yp, int size){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    tex->SetTextSize(size);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    //tex->SetTextFont(42);
    tex->SetNDC();
    tex->Draw();
  }  


  std::vector<fastjet::PseudoJet> GatherParticles( TStarJetVectorContainer<TStarJetVector> *container, std::vector<fastjet::PseudoJet> &rawParticles ){
    for ( int i=0; i < container->GetEntries() ; ++i ) {
      TStarJetVector* sv = container->Get(i);
      fastjet::PseudoJet current = fastjet::PseudoJet( *sv );
      current.set_user_index( sv->GetCharge() );
      if ( std::abs(current.eta()) > maxTrackEta )      { continue; }  // removes particles with |eta|>1
      if ( current.pt() < partMinPt )      { continue; }  // removes particles with pt<0.2GeV

      rawParticles.push_back(current);
    }
    return rawParticles;
  }


  void GenerateWeightedPtResponse( TH1D *h_DetWt[3], TH1D *h_Det[21], TH1D *h_MissedJets, TString plot_dir ) {

    for (int p=0; p<nPtBins; ++p) {
      TString name = "hPtResponse"; name += ptBinName[p];
      h_DetWt[p] = new TH1D( name, ";det-level leading jet p_{T} (GeV)",55,4.5,59.5);
    }

    TCanvas * can5 = new TCanvas( "can5" , "" ,700 ,500 );              // CANVAS 5

    for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
      for (int i=0; i<21; ++i) {
	int ptVal = i + 5;
	int binno = i + 1;
	if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	  double wt = 1./h_MissedJets->GetBinContent(binno);  // weight according to misses at part-level
	  h_DetWt[p]->Add(h_Det[i],wt);
	}
      }    
      h_DetWt[p]->Scale(1./h_DetWt[p]->Integral());
      h_DetWt[p]->GetYaxis()->SetRangeUser(0.0,0.2);
      h_DetWt[p]->SetMarkerStyle(ptMarker[p]);
      h_DetWt[p]->SetStats(0);
      h_DetWt[p]->Draw("SAME PLC PMC");
    }
    can5->BuildLegend(0.4,0.68,0.9,0.9);
    can5->SaveAs("plots/test/WeightedPtResponse.pdf","PDF");

    for (int p=0; p<nPtBins; ++p) { h_DetWt[p]->SetStats(1); }
    
  }


  void GetEmbeddingHistograms( TFile *embedFile, TH2D *h_PtResponse, TH1D *h_FakeJets, TH1D *h_MissedJets, TString plot_dir ) {
    TH2D *hResponse[nEtaBins];
    TH1D *hFakes[nEtaBins];
    h_MissedJets = (TH1D*)embedFile->Get("hMisses");

    int nAccepted = 0; int nFakes = 0;
  
    for (int e=0; e<nEtaBins; ++e) {    
      TString name = "hPtResponse" + etaBinName[e] + "Jet";
      hResponse[e] = (TH2D*)embedFile->Get(name);
      name = "hFakes" + etaBinName[e] + "Jet";
      hFakes[e] = (TH1D*)embedFile->Get(name);

      nAccepted += hResponse[e]->GetEntries();
      nFakes += hFakes[e]->GetEntries();
    }

    for (int e=0; e<nEtaBins; ++e) {
      h_PtResponse->Add(hResponse[e]);
      h_FakeJets->Add(hFakes[e]);
    }

    TCanvas * can2 = new TCanvas( "can2" , "" ,700 ,500 );              // CANVAS 2
    can2->SetLogy();

    int nMissed = h_MissedJets->GetEntries();
    int nEvents = nAccepted + nMissed + nFakes;

    double scale = (double)nMissed/nEvents;
    h_MissedJets->Scale(scale/h_MissedJets->Integral());
    h_MissedJets->Draw();
    TString saveName = plot_dir + "MissedJets.pdf";
    can2->SaveAs(saveName,"PDF");
    
    h_FakeJets->Scale(scale/h_FakeJets->Integral());
    h_FakeJets->Draw();
    saveName = plot_dir + "fakeJets.pdf";
    can2->SaveAs(saveName,"PDF");
    
    TCanvas * can3 = new TCanvas( "can3" , "" ,700 ,500 );              // CANVAS 3
    can3->SetLogz();
    h_PtResponse->Draw("COLZ");
    saveName = plot_dir + "pTresponse.pdf";
    can3->SaveAs(saveName,"PDF");  // SAVE 2D PT RESPONSE

  }


  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  
  void InitReader( TStarJetPicoReader & reader, TChain* chain, int nEvents, double VertexZCut, double MaxEventPtCut, double MaxEventEtCut,
		   double MinEventEtCut, double VertexZDiffCut, double DCACut, int MinNFitPointsCut, double FitOverMaxPointsCut, double MaxPtCut,
		   double MaxEtCut, TString badTowerList ) {
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    // reader.SetApplyFractionHadronicCorrection( false );
    reader.SetFractionHadronicCorrection( 0.9999 );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetVertexZCut ( VertexZCut );
    evCuts->SetRefMultCut( RefMultCut );
    evCuts->SetMaxEventPtCut( MaxEventPtCut );
    evCuts->SetMaxEventEtCut( MaxEventEtCut );
    evCuts->SetMinEventEtCut( MinEventEtCut );
    evCuts->SetVertexZDiffCut( VertexZDiffCut );
    
    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( DCACut );
    trackCuts->SetMinNFitPointsCut( MinNFitPointsCut );
    trackCuts->SetFitOverMaxPointsCut( FitOverMaxPointsCut );
    trackCuts->SetMaxPtCut ( MaxPtCut );

    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( MaxEtCut );
    //towerCuts->AddBadTowers( "src/dummy_tower_list.txt" );
    //towerCuts->AddBadTowers( "/nfs/rhi/STAR/Data/P16id/resources/bad_towers_pAu2015.list" );
    towerCuts->AddBadTowers( badTowerList );  // "src/bad_towers_pAu2015_NEW.list"

    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1
  }
  

  //gets the cross section for each pT bin for run 6   
  double LookupRun6Xsec(TString currentfile ) {

    static const Double_t Xsec[12] = {
      1.0,      //!Placeholder for 2-3                                                                                                                      
      1.30E+09, //!3-4                                                                                                                                      
      3.15E+08, //!4-5                                                                                                                                      
      1.37E+08, //!5-7                                                                                                                                      
      2.30E+07, //!7-9                                                                                                                                      
      5.53E+06, //!9-11                                                                                                                                     
      2.22E+06, //!11-15                                                                                                                                    
      3.90E+05, //!15-25                                                                                                                                    
      1.02E+04, //!25-35                                                                                                                                    
      5.01E+02, //!35-45                                                                                                                                    
      2.86E+01, //!45-55                                                                                                                                    
      1.46E+00  //!55-65                                                                                                                                    
    };

    static const Double_t Nmc[12] = {
      1,      //!2-3                                                                                                                                        
      672518, //!3-4                                                                                                                                        
      672447, //!4-5                                                                                                                                        
      393498, //!5-7                                                                                                                                        
      417659, //!7-9                                                                                                                                        
      412652, //!9-11                                                                                                                                       
      419030, //!11-15                                                                                                                                      
      396744, //!15-25                                                                                                                                      
      399919, //!25-35                                                                                                                                      
      119995, //!35-45                                                                                                                                      
      117999, //!45-55                                                                                                                                      
      119999  //!55-65                                                                                                                                      
    };

    Double_t w[12];
    for ( int i=0; i<12 ; ++i ){
      w[i] = Xsec[i] / Nmc[i];
    }

    if ( currentfile.Contains("3_4") ) return w[1];
    if ( currentfile.Contains("4_5") ) return w[2];
    if ( currentfile.Contains("5_7") ) return w[3];
    if ( currentfile.Contains("7_9") ) return w[4];
    if ( currentfile.Contains("9_11") ) return w[5];
    if ( currentfile.Contains("11_15") ) return w[6];
    if ( currentfile.Contains("15_25") ) return w[7];
    if ( currentfile.Contains("25_35") ) return w[8];
    if ( currentfile.Contains("35_45") ) return w[9];
    if ( currentfile.Contains("45_55") ) return w[10];
    if ( currentfile.Contains("55_65") ) return w[11];
    return 1;
  }

  //gets the cross section for each pT bin for run 12
  double LookupRun12Xsec( TString filename ){
    const int NUMBEROFPT = 11;
    //! const char *PTBINS[NUMBEROFPT]={"2_3","3_4","4_5","5_7","7_9","9_11","11_15","15_20","20_25","25_35","35_-1"};                                      
    const static float XSEC[NUMBEROFPT] = {9.00581646, 1.461908221, 0.3544350863, 0.1513760388, 0.02488645725, 0.005845846143, 0.002304880181, 0.000342661835, 4.562988397e-05, 9.738041626e-06, 5.019978175e-07};
    const static float NUMBEROFEVENT[NUMBEROFPT] = {2100295, 600300, 600300, 300289, 300289, 300289, 160295, 100302, 80293, 76303, 23307};
    const static std::vector<std::string> vptbins={"pp12Pico_pt2_3","pp12Pico_pt3_4","pp12Pico_pt4_5","pp12Pico_pt5_7","pp12Pico_pt7_9","pp12Pico_pt9_11","pp12Pico_pt11_15","pp12Pico_pt15_20","pp12Pico_pt20_25","pp12Pico_pt25_35","35_-1"};
    
    for ( int i=0; i<vptbins.size(); ++i ){
      if ( filename.Contains(vptbins.at(i).data())) return XSEC[i] / NUMBEROFEVENT[i];
    }

    throw std::runtime_error("Not a valid filename");
    return -1;
  }


  double LookupRun15Xsec( TString filename ){
    const int NUMBEROFPT = 9;
    //! const char *PTBINS[NUMBEROFPT]={"5_7","7_9","9_11","11_15","15_25","25_35","35_45","45_55","55_65"};                                                             
    const static float XSEC[NUMBEROFPT] = {0.1604997783173907,
					   0.0279900193730690,
					   0.006924398431,
					   0.0028726057079642,
					   0.0005197051748372,
					   0.0000140447879818,
					   0.0000006505378525,
					   0.0000000345848665,
					   0.0000000016149182};
    const static float NUMBEROFEVENT[NUMBEROFPT] = {242090.0,
						    159181.0,
						    96283.0,
						    125463.0,
						    441145.0,
						    169818.0,
						    58406.0,
						    59431.0,
						    59973.0};
    const static std::vector<std::string> vptbins={"pt-hat57",
						   "pt-hat79",
						   "pt-hat911",
						   "pt-hat1115",
						   "pt-hat1525",
						   "pt-hat2535",
						   "pt-hat3545",
						   "pt-hat4555",
						   "pt-hat5565"};

    for ( int i=0; i<vptbins.size(); ++i ){
      if ( filename.Contains(vptbins.at(i).data())) return XSEC[i] / NUMBEROFEVENT[i];
    }

    throw std::runtime_error("Not a valid filename");
    return -1;
  }


  //  takes in part-level jets, det-level jets, and a container for matches; removes matched jets from original vectors and returns them
  void MissesFakesAndMatches(std::vector<fastjet::PseudoJet> &pJets, std::vector<fastjet::PseudoJet> &dJets,
			     std::vector<fastjet::PseudoJet> &pMatches, std::vector<fastjet::PseudoJet> &dMatches) {

    fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( R );

    std::vector<fastjet::PseudoJet> unmatchedPartJets = sorted_by_pt( pJets );
    std::vector<fastjet::PseudoJet> unmatchedDetJets = sorted_by_pt( dJets );
    pJets.clear();    dJets.clear();
    std::vector<fastjet::PseudoJet> matchedDetJets;
    std::vector<fastjet::PseudoJet> tempPartJets;

    while ( unmatchedDetJets.size()>0 ) {
      
      if ( unmatchedPartJets.size()==0 ) {  // when we run out of part-jets to match, put unmatched det jets in dJets and quit
	dJets = unmatchedDetJets;
	return;
      }
      
      tempPartJets.clear();      matchedDetJets.clear();
      selectMatchedJets.set_reference( unmatchedPartJets[0] );
      matchedDetJets = sorted_by_pt( selectMatchedJets( unmatchedDetJets ));
      
      if ( matchedDetJets.size()==0 ) {
	pJets.push_back(unmatchedPartJets[0]);
	unmatchedPartJets.erase(unmatchedPartJets.begin());
      }
      else {
	for ( int i=0; i<unmatchedDetJets.size(); ++i ) {
	  if ( matchedDetJets[0].delta_R( unmatchedDetJets[i] ) < 0.0001 ) {
	    pMatches.push_back( unmatchedPartJets[0] );
	    unmatchedPartJets.erase(unmatchedPartJets.begin());	
	    dMatches.push_back( unmatchedDetJets[i] );
	    // std::cout<<" \n";  std::cout<<" \n";
	    // std::cout<<unmatchedDetJets[i].pt()<<" \n"<<std::endl;
	    // for ( int j=0; j<unmatchedDetJets.size(); ++j ) { std::cout<<unmatchedDetJets[j].pt()<<", "; }
	    unmatchedDetJets.erase(unmatchedDetJets.begin()+i);
	    // std::cout<<" \n";
	    // for ( int j=0; j<unmatchedDetJets.size(); ++j ) { std::cout<<unmatchedDetJets[j].pt()<<", "; }
	  }
	}
      }
      
    }

    //  we get here if we run out of detector jets to match to the particle jets, so put unmatched part jets into pJets
    for ( int i=0; i<unmatchedPartJets.size(); ++i ) { pJets.push_back(unmatchedPartJets[i]); }
      
    return;
    
  }

  
  void ProjectAndSaveFinalUEPlots( TH1D *h_WtUEpt, TString suffix, TString plot_dir ){

    TString canvasName = "can8" + suffix;
    TCanvas * can8 = new TCanvas( canvasName , "" ,700 ,500 );              // CANVAS 8
    can8->SetLogy();
	
    double pt, effic, corrPt, corrErr;
    int ptBin;

    h_WtUEpt->SetStats(1);
    h_WtUEpt->SetAxisRange( 0.000001,10,"Y");
    h_WtUEpt->SetLineColor(kBlack);
    h_WtUEpt->SetMarkerColor(kBlack);
    h_WtUEpt->Draw();
	
    TString text = "#LT p_{T}^{ch}#GT = "; text += h_WtUEpt->GetMean(1);
    text = text(0,26);
    DrawText(text, 0.6, 0.7, 20);

    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += h_WtUEpt->Integral()/AREA;
    text = text(0,42);
    DrawText( text, 0.6, 0.55, 20 );
      
    TString name = plot_dir + "CorrectedWtUEpt"; name += suffix; name += ".pdf";
    can8->SaveAs(name,"PDF");

    can8->Destructor();
  }
  

  void ProjectAndScaleUEHistogramForAllPt( TH2D *h_ChgUE2D, TH1D *h_leadPt, TH1D *h_UEpt[55], TH1D *h_WtUEpt[nPtBins], TString dir_name ){
  
    TCanvas * can6 = new TCanvas( "can6" , "" ,700 ,500 );              // CANVAS 6

    can6->SetLogy();
    for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
      int ptVal = i + 5;
      int binno = i + 1;

      TString name = "hUEpt_"; name += ptVal; name += "GeV";
      h_UEpt[i] = (TH1D*)h_ChgUE2D->ProjectionY(name,binno,binno);

      if ( (int)h_leadPt->GetBinCenter(binno) != ptVal ) { std::cerr<<"failed pT matching"<<std::endl; }
      int nJets = h_leadPt->GetBinContent( binno );
      if (nJets>0){
	h_UEpt[i]->Scale(1./nJets);
	h_UEpt[i]->SetAxisRange( 0.00001,10, "Y");
	h_UEpt[i]->SetStats(0);
	h_UEpt[i]->Draw("SAME PLC PMC");
	h_UEpt[i]->SetStats(1);
	h_UEpt[i]->SetMarkerStyle(marker[i]);
	h_UEpt[i]->SetMarkerSize(1);
	h_UEpt[i]->SetStats(0);
	name = ""; name += ptVal; name += " GeV det. jet";
	h_UEpt[i]->SetNameTitle(name,";chg. UE particle p_{T} (GeV)");
      }
    }
    can6->BuildLegend(0.68,0.1,0.9,0.9);
    TString saveName = dir_name + "UEptByLeadPt.pdf";
    can6->SaveAs(saveName,"PDF");

  }
  

  void ProjectPartLevelJetPt( TH2D* h_PtResponse, TH1D *h_Det[21], TString plot_dir ) {

    
    for (int i=0; i<21; ++i) {  // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
      int ptVal = i + 5;
      int binno = i + 1;
    
      TString name = "hPtResponse_"; name += ptVal; name += "GeV";
      h_Det[i] = (TH1D*) h_PtResponse->ProjectionY(name,binno,binno);
      h_Det[i]->Scale(1./h_Det[i]->Integral());
      h_Det[i]->SetMarkerStyle(marker[i]);
      h_Det[i]->SetMarkerStyle(marker[i]);
      h_Det[i]->SetMarkerSize(1);
      h_Det[i]->SetStats(0);
      name = ""; name += ptVal; name += " GeV part. jet";
      h_Det[i]->SetNameTitle(name,";det-level leading jet p_{T} (GeV)");
    }


    TCanvas * can4 = new TCanvas( "can4" , "" ,700 ,500 );              // CANVAS 4
    for (int i=0; i<21; ++i) { h_Det[i]->Draw("SAME PLC PMC"); } // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
    can4->BuildLegend(0.68,0.1,0.9,0.9);
    TString saveName = plot_dir + "detPtResponses.pdf";
    can4->SaveAs(saveName,"PDF");
  
  }


  TH2D *ProjectUEHistograms( TH3D *h_ChgUE3D, TString plot_dir ) {
  
    TCanvas * can0 = new TCanvas( "can0" , "" ,700 ,500 );              // CANVAS 0
    can0->SetLogy(0);
    can0->SetLogz();
  
    // SAVE 2D HISTO OF UE PT vs. LEADING JET PT
    TH2D* h_ChgUE2D = (TH2D*)h_ChgUE3D->Project3D("YX");
    h_ChgUE2D->GetXaxis()->SetRangeUser(0.0,60.0);
    h_ChgUE2D->GetYaxis()->SetRangeUser(0.0,15.0);
    h_ChgUE2D->GetYaxis()->SetTitleOffset(1.25);
    h_ChgUE2D->Draw("COLZ");
    TString saveName = plot_dir + "ChgUE2D.pdf";
    can0->SaveAs(saveName,"PDF");
    h_ChgUE2D->GetXaxis()->SetRangeUser(4.5,59.5);
    h_ChgUE2D->GetYaxis()->SetRangeUser(0.0,30.0);

    return h_ChgUE2D;
  }

  
  void TrackingEfficiency2DCorrection( TH2D* h_ChgUE2D_te, TH2D* h_ChgUE2D_uncorr, TH1D *h_Effic ){
    
    double pt, effic, ptVal, corrVal, corrErr, jetPt; // eta, phi, e, px, py, pz, 
    int ptBin;
	
    for ( int iy=1; iy<=h_ChgUE2D_uncorr->GetNbinsY(); ++iy ) {  // loop over chg UE pT bins (y)
      pt = h_ChgUE2D_uncorr->GetYaxis()->GetBinCenter(iy);
      if ( pt > 3.0 ) { pt = 3.0; }
      ptBin = h_Effic->FindBin( pt );    // find histogram bin corresponding to track pt
      effic = h_Effic->GetBinContent( ptBin );
	
      for ( int ix=1; ix<=h_ChgUE2D_uncorr->GetNbinsX(); ++ix ) { // loop over all lead jet pt bins (x)

	ptVal = h_ChgUE2D_uncorr->GetYaxis()->GetBinCenter(iy);
	corrVal = h_ChgUE2D_uncorr->GetBinContent(ix,iy)/effic;      // calculate corrected bin content and error
	corrErr = h_ChgUE2D_uncorr->GetBinError(ix,iy)/effic;       // (divide bin content and error by efficiency)
	jetPt = h_ChgUE2D_uncorr->GetXaxis()->GetBinCenter(ix);
	ptBin = h_ChgUE2D_uncorr->GetYaxis()->FindBin( ptVal );
	  
	// h_ChgUE2D_te->Fill( jetPt, ptVal, corrVal );
	h_ChgUE2D_te->SetBinContent( ix, ptBin, corrVal );
	h_ChgUE2D_te->SetBinError( ix, ptBin, corrErr );

      }
    }
	
  }
  

  bool UseTriggerTower( int TriggerTowerId ) {
    
    int badTows[454] = { 31, 34, 59, 95, 96, 106, 113, 120, 134, 139, 157, 160, 175, 214, 224, 257, 266, 267, 275, 280, 282, 286, 287, 308, 360, 368, 385, 389, 395, 405, 410, 426, 433, 474, 479, 483, 484, 504, 509, 533, 541, 555, 561, 562, 563, 564, 585, 603, 615, 616, 627, 633, 638, 649, 650, 653, 657, 671, 673, 674, 680, 681, 693, 721, 722, 725, 740, 750, 753, 754, 755, 756, 758, 768, 773, 774, 775, 776, 779, 784, 790, 793, 794, 795, 796, 799, 801, 813, 814, 815, 816, 817, 822, 832, 835, 837, 840, 844, 846, 857, 860, 875, 880, 882, 893, 897, 899, 903, 916, 936, 939, 941, 946, 953, 954, 956, 986, 989, 993, 1005, 1012, 1020, 1023, 1026, 1027, 1028, 1045, 1046, 1048, 1057, 1063, 1080, 1081, 1085, 1100, 1104, 1120, 1125, 1128, 1130, 1132, 1142, 1154, 1158, 1159, 1160, 1171, 1180, 1183, 1184, 1187, 1189, 1190, 1197, 1200, 1202, 1204, 1207, 1214, 1217, 1219, 1220, 1221, 1224, 1232, 1233, 1237, 1238, 1244, 1256, 1257, 1263, 1280, 1283, 1294, 1301, 1306, 1312, 1313, 1318, 1329, 1341, 1348, 1353, 1354, 1369, 1375, 1378, 1388, 1400, 1401, 1405, 1407, 1409, 1434, 1439, 1440, 1448, 1475, 1486, 1537, 1563, 1564, 1567, 1574, 1575, 1588, 1592, 1597, 1599, 1602, 1612, 1654, 1668, 1679, 1705, 1709, 1720, 1728, 1753, 1762, 1765, 1766, 1773, 1776, 1789, 1807, 1823, 1840, 1856, 1866, 1877, 1878, 1879, 1880, 1921, 1945, 1952, 1976, 1983, 1984, 2027, 2032, 2043, 2051, 2066, 2073, 2077, 2092, 2093, 2097, 2104, 2107, 2111, 2141, 2160, 2162, 2168, 2175, 2176, 2177, 2190, 2193, 2194, 2195, 2196, 2197, 2213, 2214, 2215, 2216, 2217, 2223, 2233, 2234, 2235, 2236, 2253, 2254, 2255, 2256, 2299, 2300, 2303, 2305, 2313, 2340, 2390, 2391, 2392, 2414, 2415, 2417, 2439, 2459, 2476, 2493, 2520, 2569, 2580, 2589, 2590, 2633, 2697, 2737, 2749, 2822, 2834, 2863, 2865, 2874, 2929, 2953, 2954, 2955, 2961, 2969, 2973, 2974, 2975, 2976, 2981, 2993, 2994, 2995, 3005, 3020, 3063, 3070, 3071, 3146, 3160, 3167, 3186, 3233, 3234, 3235, 3236, 3253, 3254, 3255, 3256, 3263, 3273, 3274, 3275, 3276, 3293, 3294, 3295, 3296, 3299, 3300, 3337, 3354, 3355, 3356, 3360, 3362, 3385, 3407, 3436, 3451, 3473, 3481, 3492, 3493, 3494, 3495, 3498, 3504, 3513, 3515, 3544, 3584, 3588, 3611, 3666, 3668, 3670, 3678, 3679, 3682, 3690, 3692, 3702, 3718, 3720, 3725, 3737, 3738, 3739, 3741, 3769, 3777, 3822, 3831, 3834, 3838, 3840, 3859, 3861, 3984, 4006, 4013, 4017, 4018, 4019, 4047, 4053, 4057, 4079, 4099, 4104, 4130, 4169, 4177, 4217, 4223, 4302, 4331, 4350, 4355, 4357, 4405, 4440, 4458, 4460, 4469, 4480, 4495, 4496, 4500, 4518, 4534, 4563, 4595, 4659, 4677, 4678, 4712, 4742, 4743, 4744, 4762, 4763, 4764, 4766, 4768, 4778, 4781, 4782, 4783, 4784 };
    
    for( int i=0; i<454; ++i ) {
      if( badTows[i]==TriggerTowerId ) { return false; }
    }
    return true;
  }

  
  void WeightUEPtByLeadPtAndFakes( TH1D *h_WtUEpt[nPtBins],TH1D *h_UEpt[55],TH1D *h_DetWt[nPtBins],TH1D *h_leadPt,TH1D *h_FakeJets,TString plot_dir ){

    for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
      int ptVal = i + 5;
      int binno = i + 1;
      if (h_leadPt->GetBinContent( binno )>0){
        for (int p=0; p<nPtBins; ++p) {  // ADD UE DISTRIBUTIONS WEIGHTED BY DET-LEVEL FRACTIONAL CONTRIBUTION
	  double wt = h_DetWt[p]->GetBinContent(binno)*( 1.0 - h_FakeJets->GetBinContent(binno) );
	  h_WtUEpt[p]->Add( h_UEpt[i], wt );
	}
      }
    }
  
    TCanvas * can7 = new TCanvas( "can7" , "" ,700 ,500 );              // CANVAS 7
    for (int p=0; p<nPtBins; ++p) {
      h_WtUEpt[p]->SetMarkerStyle(ptMarker[p]);
      h_WtUEpt[p]->SetStats(0);
      h_WtUEpt[p]->Draw("PLC PMC SAME");    
    }
    can7->SetLogy();
    can7->BuildLegend(0.4,0.68,0.9,0.9);
    TString saveName = plot_dir + "weightedUEptByLeadPt.pdf";
    can7->SaveAs(saveName,"PDF");

  }
  

}
