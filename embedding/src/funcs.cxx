//!  functions.cxx
//!  Isaac Mooney, WSU - June 2019

#include "params.hh"
#include "funcs.hh"
//!#include "ktTrackEff.hh"

typedef fastjet::contrib::SoftDrop SD;

namespace Analysis {

  //! -------------------------                                                                                                                                                                            
  //! IO/OS Manip functionality                                                                                                                                                                            
  //! -------------------------                                                                                                                                                                             
  //! Used to understand which format of input file is being used                                                                                                                                          
  //! ( .root file, .txt, .list, etc )                                                                                                                                                                     
  //! ---------------------------------------------------------------------

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
  

    //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from which they came. Both the Py & Py+Ge event will be thrown out.                                                                                                                                                       
  //NOTE: I don't think the Picos save the pT-hat bin edges from the production, so I have to use the filenames instead                                         
  bool DiscardEmbedEvent(const TString Filename, const std::vector<fastjet::PseudoJet> p_Jets, const std::vector<fastjet::PseudoJet> g_Jets) {
    bool bad_event = 0;
    //have to do some manipulations to get the upper edge of the pT-hard bin from the file name, e.g. isolating "25" from "1525"                                
    //unfortunately this only works for the naming convention of the Pythia6 dijet embedding into Y15 pAu MB files on disk. Will need to be rewritten slightly if used on different files.                                                                                                                                     
    //upper-bin edges are either single or double digits, e.g. 7 or 15. So 10 positions from the end of the file name will either be the full 15 or the 5 that comes before the 7, in the previous two examples. There are only two sets of files with this single-digit problem: 57 and 79 (i.e. pT-hat in (7,9)). Therefore, we can check if the first two characters in this 10 character string are 57 or 79. Since these are not numbers that appear in the pT-hat ranges (5-7-9-11-15-25-35-45-55-65), we can separate into the case of single or double digit high-end by checking for these.                                                         
    std::string tail = ((std::string) Filename).substr(((std::string) Filename).size() - 10);
    std::string upstring = tail.substr(0,2);
    //if we found a "57" or "79" in the 2 digit string, we just take the last digit.                                                                            
    if (upstring.find("57") != std::string::npos || upstring.find("79") != std::string::npos) {
      upstring = upstring.substr(1,1);
    }
    int upbin = std::stoi(upstring);
    if (p_Jets.size() != 0) {
      if (p_Jets[0].pt() > 2*upbin) {
        std::cout << "DEBUG: 2upbin = " << 2*upbin << std::endl;
        std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << p_Jets[0].pt() << " " << p_Jets[0].eta() << " " << p_Jets[0].phi() << " " << p_Jets[0].m() << std::endl;
        bad_event = 1;
      }
      else {
        std::cout << "DEBUG: should be less than " << 2*upbin << ": " << p_Jets[0].pt() << std::endl;
      }
    }
    if (g_Jets.size() != 0) {
      //std::cout << "DEBUG: should never see this for toy_embedding since we use a dummy vector here which should have size == 0" << std::endl;                
      if (g_Jets[0].pt() > 2*upbin) {
        std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << g_Jets[0].pt() << " " << g_Jets[0].eta() << " " << g_Jets[0].phi() << " " << g_Jets[0].m() << std::endl;
        bad_event = 1;
      }
    }

    return bad_event;
  }

  bool HasEnding (std::string const &full_string, std::string const &ending) {
    if (full_string.length() >= ending.length()) {
      return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
    } else {
      return false;
    }
  }

  //! parse a CSV file to a set of unique entries.
  //! All comments must start on their own line, and be proceeded   
  //! by a pound sign (#)                                                                                                                                                              
  template <typename T>
  std::set<T> ParseCSV(std::string csv) {
    //! return set                                                                                                        
    std::set<T> ret;
    std::ifstream fs(csv);
    std::string line;
    //! first, split by line
    while (std::getline(fs, line)) {
      if (line.size() == 0) //!reject empty lines   
	continue;
      if (line[0] == '#') //!reject comments                                                                                                                                        
	continue;
      //! split the string by commas                                                                                                                                                   
      std::istringstream ss(line);
      while (ss) {
	std::string str_value;
	std::getline(ss, str_value, ',');
	if (CanCast<T>(str_value)) {
	  ret.insert(CastTo<T>(str_value));
	}
      }
    }
    return ret;
  }

  template<typename T>
  bool CanCast(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return iss && iss.eof();
  }

  template<typename T>
  T CastTo(std::string s) {
    std::istringstream iss(s);
    T dummy;
    iss >> std::skipws >> dummy;
    return dummy;
  }

  
  void drawText(const char *text, float xp, float yp, int size){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    tex->SetTextSize(size);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    //tex->SetTextFont(42);
    tex->SetNDC();
    tex->Draw();
  }
  
  
  //When we want e.g. R_jet = 0.4 for analysis, and take in "04" on the command line, this converts to (double) 0.4.
  double radius_str_to_double (std::string radius_str) {
    std::string radiusNum = radius_str.substr(0,1)+"."+radius_str.substr(1,1); //e.g. 0.4                                                               
    double radius_double = (double) std::stod(radiusNum); //converting the string 0.4 to a double
    std::cout << "DEBUG: jet radius = " << radius_double << std::endl;
    
    return radius_double;
  }

  /* //  INITIATE READER with some selections
  void InitReader( TStarJetPicoReader * reader, TChain* chain, int nEvents, const std::string trig, const double vZ, const double vZDiff, const double Pt, const double Et, const double Etmin, const double DCA, const double NFit, const double NFitRatio, const double maxEtTow, const double hc, const bool mip_correction, const std::string badTows, const std::string bad_run_list) {
    // set the chain
    if (chain != nullptr) {reader->SetInputChain( chain );}
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader->SetApplyFractionHadronicCorrection( true );
    reader->SetFractionHadronicCorrection( hc ); //0.9999
    reader->SetRejectTowerElectrons( kFALSE );
    reader->SetApplyMIPCorrection( mip_correction );
    // if bad run list is specified, add to reader
    if (!bad_run_list.empty()) {
      std::set<int> bad_runs = ParseCSV<int>(bad_run_list);
      for (auto run : bad_runs) {
	reader->AddMaskedRun(run);
      }
    }
    // Event and track selection - all explained in params.hh
    // -------------------------
    TStarJetPicoEventCuts* evCuts = reader->GetEventCuts();
    evCuts->SetTriggerSelection( trig.c_str() ); //All, MB, HT, pp, ppHT, ppJP2. For pAu I think trigger strings have now been hardcoded, but haven't tested
    evCuts->SetVertexZCut ( vZ );
    evCuts->SetVertexZDiffCut( vZDiff );
    evCuts->SetRefMultCut( refMultCut );
    evCuts->SetMaxEventPtCut( Pt );
    evCuts->SetMaxEventEtCut( Et );
    evCuts->SetMinEventEtCut( Etmin );
    
    std::cout << "Using these event cuts:" << std::endl;
    std::cout << " trigger: " << evCuts->GetTriggerSelection() << std::endl;
    std::cout << " vz: " << evCuts->GetVertexZCut() << std::endl;
    std::cout << " |vpdvz - tpcvz|: " << evCuts->GetVertexZDiffCut() << std::endl;
    std::cout << " event max pT: " << evCuts->GetMaxEventPtCut() << std::endl;
    std::cout << " event max Et: " << evCuts->GetMaxEventEtCut() << std::endl;
    std::cout << " event min Et: " << evCuts->GetMinEventEtCut() << std::endl;

    // Tracks cuts - all explained in params.hh
    TStarJetPicoTrackCuts* trackCuts = reader->GetTrackCuts();
    trackCuts->SetDCACut( DCA );
    trackCuts->SetMinNFitPointsCut( NFit );
    trackCuts->SetFitOverMaxPointsCut( NFitRatio );    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    
    // Towers - all explained in params.hh
    TStarJetPicoTowerCuts* towerCuts = reader->GetTowerCuts();
    towerCuts->SetMaxEtCut( maxEtTow );
    towerCuts->AddBadTowers( badTows );
    std::cout << badTows << std::endl;
    
    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader->SetProcessV0s(false);
    
    // Initialize the reader
    reader->Init( nEvents ); //runs through all events with -1
    }*/

  

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

  // //This function simply converts TStarJetVectors from the event into PseudoJets for later clustering into jets with FastJet.
  // //It also very importantly assigns a mass to each particle after the conversion.
  // //The assigned mass is dependent on whether the function call is for particle-level (e.g. Pythia) [where we know the rest masses] or detector-level (e.g. Geant, data) [where we don't know them].
  // void GatherParticles ( TStarJetVectorContainer<TStarJetVector> * container, TStarJetVector *sv, std::vector<fastjet::PseudoJet> & Particles, const bool full, const bool py){  //, TDatabasePDG *pdg){
  //   for ( int i = 0; i < container->GetEntries() ; ++i ) {
  //     sv = container->Get(i);
  //     //cout << "DEBUG: TStarJetVector info: " << 
  //     fastjet::PseudoJet current = fastjet::PseudoJet( *sv );


  //     if ((sv->GetCharge() == 0) && (full == 0)) { continue; } //!if we don't want full jets, skip neutrals                                                     
  //     current.set_user_index( sv->GetCharge() );
      
  //     //DEBUG:
  //     //std::cout << "ending with charge: " << current.user_index() << " for particle " << i << std::endl;
      
  //     Particles.push_back(current);
  //   }
  //   return;
  // }

  //this function takes in a jet sample, and keeps whichever ones pass a neutral energy fraction selection (i.e. have >= x% of their energy in tracks) 
  //!the name is a bit of a misnomer, we're actually using pT, not energy
  void ApplyNEFSelection(const std::vector<fastjet::PseudoJet> init, std::vector<fastjet::PseudoJet> & result) {
    //!Implementing a neutral energy fraction cut on inclusive jets 
    for (int i = 0; i < init.size(); ++ i) { //looping over the jet
      double towersum = 0; double ptsum = 0; //start with zero energy, neutral or otherwise
      for (int j = 0; j < init[i].constituents().size(); ++ j) { //loop over the jet constituents
        if (init[i].constituents()[j].user_index() == 0) { //if it's neutral, add its pT to the towers
          towersum += init[i].constituents()[j].pt();
        }
        ptsum += init[i].constituents()[j].pt(); //whether or not it's neutral, add its pT to the total
      }//constituent loop
      if (towersum / (double) ptsum < NEF_max) { //if the fraction in the towers is less than the threshold, accept it
        result.push_back(init[i]);
      }//if
    }//jet loop
    return;
  }
  
  //discards events on the grounds of them having jets of pT > double the high end of the pT-hat bin from which they came. Both the Py & Py+Ge event will be thrown out.
  //NOTE: I don't think the Picos save the pT-hat bin edges from the production, so I have to use the filenames instead
  bool DiscardEvent(const TString Filename, const std::vector<fastjet::PseudoJet> p_Jets, const std::vector<fastjet::PseudoJet> g_Jets) {
    bool bad_event = 0;
    //have to do some manipulations to get the upper edge of the pT-hard bin from the file name, e.g. isolating "25" from 20_25
    //unfortunately this only works for the naming convention of the Y12 Pythia6 and Pythia6+Geant files on disk. Will need to be rewritten slightly if
    //used on different files.
    std::string tail = ((std::string) Filename).substr(((std::string) Filename).size() - 10);
    std::string upstring = tail.substr(0,2);
    std::string upstring_copy = upstring;
    //since we search from the end of the string backward, there is a difference between files with two 2-digit pt-hat edges, with one, and with zero.
    if (upstring.find("_") != std::string::npos || upstring.find("-") != std::string::npos) {
      if (upstring.substr(1,1) != "_") {
	upstring = upstring.substr(1,1);
      }
      else {
	upstring = upstring.substr(0,1);
      }
    }
    int upbin = std::stoi(upstring);
    if (p_Jets.size() != 0) {
      if ((p_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << p_Jets[0].pt() << " " << p_Jets[0].eta() << " " << p_Jets[0].phi() << " " << p_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    if (g_Jets.size() != 0) {
      if ((g_Jets[0].pt() > 2*upbin) && upstring_copy != "-1") {
	std::cout << "Removing event from file " << Filename << /*<< " with weight " << mc_weight*/ " due to bad jet with pt, eta, phi, and m: " << g_Jets[0].pt() << " " << g_Jets[0].eta() << " " << g_Jets[0].phi() << " " << g_Jets[0].m() << std::endl;
	bad_event = 1;
      }
    }
    
    return bad_event;
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

  
  // //This function takes a vector of jets to be geometrically matched to another vector of "candidate" matches. Once matching occurs, a vector of indices is returned allowing one to index the original two vectors to fill responses etc. with the matched jet pairs. Redundantly (for debugging, etc.) the vectors of matches themselves are also updated for later use.
  // //Note: We need to be able to remove jets from the "candidates" vector after they've been matched, so we make a copy in the function. Also make a copy candidates vector for each iteration on toMatch since this vector has selections applied to it
  // //Note: In finding which jets were the matches, we know the toMatch jet match will be the 'i'th jet since we are iterating. The candidate_copy jet should be the highest pT match, so the first one in the candidate_copy list. Geometrically match the candidate_copy jet to the nearest candidate jet, since jets have been removed so they don't index to the same jet anymore
  // std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & c_matches, std::vector<fastjet::PseudoJet> & t_matches) {
  //   std::vector<int> match_indices;
  //   if (candidates_safe.size() == 0 || toMatch.size() == 0) {
  //     return match_indices; //later, match_indices being empty will tell us there were no matches
  //   }
  //   //define candidates outside the loop so list continually dwindles as we remove matched candidates
  //   std::vector<fastjet::PseudoJet> candidates = candidates_safe;
  //   for (int i = 0; i < toMatch.size(); ++ i) { //for each jet in toMatch, we try to find a match from candidates_copy
  //     //defined inside the loop so that for each toMatch jet there's a new set of candidates
  //     std::vector<fastjet::PseudoJet> candidates_copy = candidates;
  //     fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( R );
  //     selectMatchedJets.set_reference( toMatch[i] );
  //     //note: matchedToJet and candidates_copy are equivalent, assuming candidates_safe was already sorted by pT
  //     std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt( selectMatchedJets( candidates_copy ));
  //     if (matchedToJet.size() == 0) { continue; } //means no match to this jet. Remove none from candidates. Continuing on to the next one.
  //     else { //found at least one match. Need to remove the highest pT one from candidates and add the respective jets to the match vectors.
  // 	match_indices.push_back(i); //push back the toMatch match position
  // 	t_matches.push_back(toMatch[i]);
  // 	c_matches.push_back(matchedToJet[0]); //highest pT match
  // 	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
  // 	  if (matchedToJet[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet
  // 	    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
  // 	    match_indices.push_back(j); //push back the candidate match position
  // 	    break; //should exit only the c_matches loop.
  // 	  }
  // 	}
  //     }
  //   }
  //   return match_indices;
  // }
  
  //!this function is similar to the "MatchJets" function, but it instead takes a list of jets which have already been matched ("candidates_safe") and finds the jets to which they correspond in the list of unmatched jets ("toMatch") by (exact) geometrical matching. The remaining jets with no corresponding already-matched jets are either misses or fakes depending on the call.
  std::vector<int> FakesandMisses(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & unmatched) {
    std::vector<int> miss_fake_index;
    std::vector<fastjet::PseudoJet> candidates = candidates_safe;
    for (int i = 0; i < toMatch.size(); ++ i) {
      std::vector<fastjet::PseudoJet> candidates_copy = candidates;
      fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( 0.0001 ); //a "match" is now if we found the same exact jet.
      selectMatchedJets.set_reference( toMatch[i] );
      std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt( selectMatchedJets( candidates_copy ));
      if (matchedToJet.size() == 0) { //means no match to this jet. Remove none from candidates. Add it to unmatched & continue to next one.
	miss_fake_index.push_back(i); //the ith jet in toMatch is a miss or a fake
	unmatched.push_back(toMatch[i]);
	continue;
      }
      else { //found at least one match. Need to remove the highest pT one from candidates. [Should just be 1 match]
	for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
	  if (matchedToJet[0].delta_R(candidates[j]) < 0.0001) { //is probably the same jet
	    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
	    break; //should exit only the candidates loop.
	  }
	}
      }
    }
    return miss_fake_index;
  }
  
  //~~~~~~~FILL HISTS [not used currently, but might be useful in the future so saving it here]~~~~~~~//
  
  void FillHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "jet").c_str(), jet.m(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c2D.fill(("m_v_pt_rebin_" + flag + "_" + "jet").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "jet").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    for (int cons = 0; cons < jet.constituents().size(); ++ cons) {
      if (jet.constituents()[cons].pt() < partMinPt) {continue;} //ignores contributions from ghosts                       
      c3D.fill(("PtEtaPhi_" + flag + "_" + "cons").c_str(), jet.constituents()[cons].pt(), jet.constituents()[cons].eta(), jet.constituents()[cons].phi(), weight);
    }
    return;
  }
  
  void FillHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }
    /*
    //trigger & recoil
    std::vector<fastjet::PseudoJet> candidates;
    bool which_one = GetTriggerJet(candidates, jets);
    if (candidates.size() == 0 || jets.size() < 2) { // means there isn't a trigger or there isn't a recoil
      return;
    }
    if (candidates.size() == 1 && jets.size() > 1) { //potential trigger
      if (fabs(fabs(jets[which_one].delta_phi_to(jets[(which_one + 1) % 2])) - Pi) < R) { //found a recoil
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", jets[which_one], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", jets[(which_one + 1) % 2], weight); //filling hists for recoil
	return;
      }
    }
    if (candidates.size() == 2) {
      if (fabs(fabs(candidates[0].delta_phi_to(candidates[1])) - Pi) < R) { //trigger & recoil found!
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "trig", candidates[0], weight); //filling hists for trigger
	FillHistsHelper(c1D, c2D, c3D, flag1, flag2, "rec", candidates[1], weight); //filling hists for recoil
	return;
      }
    }
    */
    return;
  }

  void FillSDHistsHelper(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> &c2D, Collection<std::string, TH3D> & c3D, const std::string flag, const fastjet::PseudoJet jet, const double weight) {
    c1D.fill(("m_" + flag + "_" + "sd").c_str(), jet.m(), weight);
    c1D.fill(("zg_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().symmetry(), weight);
    c1D.fill(("thetag_" + flag + "_" + "sd").c_str(), jet.structure_of<fastjet::contrib::SoftDrop>().delta_R(), weight);
    c2D.fill(("m_v_pt_" + flag + "_" + "sd").c_str(), jet.m(), jet.pt(), weight);
    c3D.fill(("PtEtaPhi_" + flag + "_" + "sd").c_str(), jet.pt(), jet.eta(), jet.phi(), weight);
    return;
  }
   

  void FillSDHists(Collection<std::string, TH1D> & c1D, Collection<std::string, TH2D> & c2D, Collection<std::string, TH3D> & c3D, const std::vector<fastjet::PseudoJet> jets, const double weight) {
    //leading
    FillSDHistsHelper(c1D, c2D, c3D, "lead", jets[0], weight);
    //subleading
    if (jets.size() > 1) {
      FillSDHistsHelper(c1D, c2D, c3D, "sublead", jets[1], weight);
    }
    //inclusive
    for (int i = 0; i < jets.size(); ++ i) {
      FillSDHistsHelper(c1D, c2D, c3D, "incl", jets[i], weight);
    }

  }



  bool UseTriggerTower( int TriggerTowerId ) {
    
    int badTows[454] = { 31, 34, 59, 95, 96, 106, 113, 120, 134, 139, 157, 160, 175, 214, 224, 257, 266, 267, 275, 280, 282, 286, 287, 308, 360, 368, 385, 389, 395, 405, 410, 426, 433, 474, 479, 483, 484, 504, 509, 533, 541, 555, 561, 562, 563, 564, 585, 603, 615, 616, 627, 633, 638, 649, 650, 653, 657, 671, 673, 674, 680, 681, 693, 721, 722, 725, 740, 750, 753, 754, 755, 756, 758, 768, 773, 774, 775, 776, 779, 784, 790, 793, 794, 795, 796, 799, 801, 813, 814, 815, 816, 817, 822, 832, 835, 837, 840, 844, 846, 857, 860, 875, 880, 882, 893, 897, 899, 903, 916, 936, 939, 941, 946, 953, 954, 956, 986, 989, 993, 1005, 1012, 1020, 1023, 1026, 1027, 1028, 1045, 1046, 1048, 1057, 1063, 1080, 1081, 1085, 1100, 1104, 1120, 1125, 1128, 1130, 1132, 1142, 1154, 1158, 1159, 1160, 1171, 1180, 1183, 1184, 1187, 1189, 1190, 1197, 1200, 1202, 1204, 1207, 1214, 1217, 1219, 1220, 1221, 1224, 1232, 1233, 1237, 1238, 1244, 1256, 1257, 1263, 1280, 1283, 1294, 1301, 1306, 1312, 1313, 1318, 1329, 1341, 1348, 1353, 1354, 1369, 1375, 1378, 1388, 1400, 1401, 1405, 1407, 1409, 1434, 1439, 1440, 1448, 1475, 1486, 1537, 1563, 1564, 1567, 1574, 1575, 1588, 1592, 1597, 1599, 1602, 1612, 1654, 1668, 1679, 1705, 1709, 1720, 1728, 1753, 1762, 1765, 1766, 1773, 1776, 1789, 1807, 1823, 1840, 1856, 1866, 1877, 1878, 1879, 1880, 1921, 1945, 1952, 1976, 1983, 1984, 2027, 2032, 2043, 2051, 2066, 2073, 2077, 2092, 2093, 2097, 2104, 2107, 2111, 2141, 2160, 2162, 2168, 2175, 2176, 2177, 2190, 2193, 2194, 2195, 2196, 2197, 2213, 2214, 2215, 2216, 2217, 2223, 2233, 2234, 2235, 2236, 2253, 2254, 2255, 2256, 2299, 2300, 2303, 2305, 2313, 2340, 2390, 2391, 2392, 2414, 2415, 2417, 2439, 2459, 2476, 2493, 2520, 2569, 2580, 2589, 2590, 2633, 2697, 2737, 2749, 2822, 2834, 2863, 2865, 2874, 2929, 2953, 2954, 2955, 2961, 2969, 2973, 2974, 2975, 2976, 2981, 2993, 2994, 2995, 3005, 3020, 3063, 3070, 3071, 3146, 3160, 3167, 3186, 3233, 3234, 3235, 3236, 3253, 3254, 3255, 3256, 3263, 3273, 3274, 3275, 3276, 3293, 3294, 3295, 3296, 3299, 3300, 3337, 3354, 3355, 3356, 3360, 3362, 3385, 3407, 3436, 3451, 3473, 3481, 3492, 3493, 3494, 3495, 3498, 3504, 3513, 3515, 3544, 3584, 3588, 3611, 3666, 3668, 3670, 3678, 3679, 3682, 3690, 3692, 3702, 3718, 3720, 3725, 3737, 3738, 3739, 3741, 3769, 3777, 3822, 3831, 3834, 3838, 3840, 3859, 3861, 3984, 4006, 4013, 4017, 4018, 4019, 4047, 4053, 4057, 4079, 4099, 4104, 4130, 4169, 4177, 4217, 4223, 4302, 4331, 4350, 4355, 4357, 4405, 4440, 4458, 4460, 4469, 4480, 4495, 4496, 4500, 4518, 4534, 4563, 4595, 4659, 4677, 4678, 4712, 4742, 4743, 4744, 4762, 4763, 4764, 4766, 4768, 4778, 4781, 4782, 4783, 4784 };
    
    for( int i=0; i<454; ++i ) {
      if( badTows[i]==TriggerTowerId ) { return false; }
    }
    return true;
  }

}
