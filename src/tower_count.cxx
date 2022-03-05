// tower_count.cxx
// David Stewart  2022.03.04 -> modified from V. Verkest code

#include "pAuFunctions.hh"
/* #include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!! */
#include <iostream>
#include <fstream>

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

// funcions and cuts specified in pAuFunctions.hh
int main ( int argc, const char** argv ) {
    int number_of_events;		
    string inFile, outFile;		
    TString name, title;

    vector<string> arguments( argv+1, argv+argc );
    if ( argc ==  4 ) {    
        inFile = arguments[0];    
        outFile = arguments[1];    
        number_of_events = atoi(arguments[2].c_str()); 
    } else if ( argc==1 ) { 
        inFile="production_pAu200_2015/HT/pAu_2015_200_HT*.root"; 
        outFile="out/UE/pAuHTjetUE.root"; number_of_events=-1; 
    } else { cerr<< "incorrect number of command line arguments"; return 5000; }


    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;


    // read in a vector of good run list
    vector<int> good_runs;
    ifstream ifile;
    ifile.open("lists/good_run_id.list");
    string line;
    // skip the first line which is a comment line
    getline(ifile,line);
    while (getline(ifile,line)) {
        istringstream iline { line.c_str() };
        int id, runid, dur;
        iline >> id >> runid >> dur;
        good_runs.push_back(runid);
    } // note that good_runs are numerically sorted already


    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    TChain* Chain = new TChain( "JetTree" );
    Chain->Add( inFile.c_str() );
    TStarJetPicoReader Reader;
    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/dummy_tower.list");
    towerCuts->SetMaxEtCut( 1.0 );

    int numEvents = number_of_events;        // total events in HT: 152,007,032
    /* string bad_tower_option {"allTowers"}; */
    InitReader( Reader, Chain, numEvents );
    TFile *pAuFile = new TFile( outFile.c_str() ,"RECREATE");

    TProfile pr_Et_all { 
        "pr_Et_all","All runs; Tower ID;mean w/std-dev",4800, 0.5, 4800.5, "s" };
    TProfile pr_Et_bad { 
        "pr_Et_bad","Runs bad by Issac; Tower ID;mean w/std-dev",4800, 0.5, 4800.5, "s" };
    TProfile pr_Et_good{ 
        "pr_Et_good","Runs good by Issac; Tower ID;mean w/std-dev",4800, 0.5, 4800.5, "s" };

    TH1D hg_Et_all { "hg_Et_all","All runs; Tower ID;mean w/std-dev",4800, 0.5, 4800.5 };
    TH1D hg_Et_bad { "hg_Et_bad","Runs bad by Issac; Tower ID;mean w/std-dev",4800, 0.5, 4800.5 };
    TH1D hg_Et_good { "hg_Et_good","Runs good by Issac; Tower ID;mean w/std-dev",4800, 0.5, 4800.5 };

    TH1D nevents { "nevents","Number of events;0-all,1-bad,2-good,3-skipped (not on good file list)", 
        4, -0.5, 3.5 };

    while ( Reader.NextEvent() ) {
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();

        if (!(header->HasTriggerId(500001) || header->HasTriggerId(500004))) {continue;} 

        double Vz    { header->GetPrimaryVertexZ() };
        double vzVpd { header->GetvpdVz() };
        int    runid { header->GetRunId() };


        bool is_good_run { binary_search(good_runs.begin(), good_runs.end(), runid) };
        bool is_isaac_skip_run {
             (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001)
          || header->GetRunId() == 16135031 
          || header->GetRunId() == 16135032 };

        nevents.Fill(0);
        if (!is_good_run) nevents.Fill(3.);
        else nevents.Fill( is_isaac_skip_run ? 1. : 2. );

        for (int i{0}; i<header->GetNOfTowers(); ++i) {
            TStarJetPicoTower* tower { event->GetTower(i) };

            int    id { tower->GetId() };
            double Et { tower->GetEt() };

            if (Et<0.2) continue;
            pr_Et_all.Fill(id,Et);
            (is_isaac_skip_run ? pr_Et_bad : pr_Et_good).Fill(id,Et);

            hg_Et_all.Fill(id);
            (is_isaac_skip_run ? hg_Et_bad : hg_Et_good).Fill(id);
        }
    }

    hg_Et_all  .Write();
    hg_Et_bad  .Write();
    hg_Et_good .Write();

    pr_Et_all  .Write();
    pr_Et_bad  .Write();
    pr_Et_good .Write();
    nevents    .Write();

    pAuFile->Save();
    pAuFile->Close();

    return 0;
}
