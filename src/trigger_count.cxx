// tower_count.cxx
// David Stewart  2022.03.04 -> modified from V. Verkest code

#include "pAuFunctions.hh"
/* #include "pAu_HT_jetParameters.hh"      //  BbcAdcEastSum VS. BbcAdcSumEast!! */
#include <iostream>
#include <fstream>
#include <map>

using namespace std;
using namespace fastjet;
using namespace pAuAnalysis;

// functions and cuts specified in pAuFunctions.hh
// count the trigger's available on a per-run basis
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

    // write to an output text file:
    // runID  triggerID  n_events
    
    // there will be a variable number of line entries per run
    map<int,map<int,long int>> data;
    map<int,int> n_events;

    /* // read in a vector of good run list */
    /* vector<int> good_runs; */
    /* ifstream ifile; */
    /* ifile.open("lists/good_run_id.list"); */
    /* string line; */
    /* // skip the first line which is a comment line */
    /* getline(ifile,line); */
    /* while (getline(ifile,line)) { */
    /*     istringstream iline { line.c_str() }; */
    /*     int id, runid, dur; */
    /*     iline >> id >> runid >> dur; */
    /*     good_runs.push_back(runid); */
    /* } // note that good_runs are numerically sorted already */


    // ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~  BEGIN EVENT LOOP!  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    TChain* Chain = new TChain( "JetTree" );
    Chain->Add( inFile.c_str() );
    TStarJetPicoReader Reader;
    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/dummy_tower.list");
    towerCuts->SetMaxEtCut( 1.0 );
    /* cout << " a0 " << endl; */
    /* int numEvents = number_of_events;        // total events in HT: 152,007,032 */
    InitReader( Reader, Chain, number_of_events );

    while ( Reader.NextEvent() ) {
    /* cout << " a1 " << endl; */
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();
        int runid { header->GetRunId() };
        if (data.find(runid) == data.end()) {
            data[runid] = map<int, long int>{};
        }
        if (nevents.find(nevents.end())) nevents[runid] = 1;
        else                    kkkkkkkkkkk

        int n_triggers =  header->GetNOfTriggerIds();
        /* cout << "n_Triggers: " << n_triggers << endl; */

        const TArrayI& trig_array = *(header->GetTriggerIdArray());
        auto& l_dat = data[runid];

        for (int i{0}; i<n_triggers; ++i) {
            int trig = trig_array[i]; 
            if (l_dat.find(trig) == l_dat.end()) l_dat[trig] = 1;
            else                      ++l_dat[trig];
        }
            /* cout << "trigger " << i << " : " << (*trig_array)[i] << endl; */
            /* int k = (*trig_array)[i]; */
            /* cout << "trigger " << i << " : " << k << endl; */
        /* } */
        /* continue; */
    }

    
    for (auto& sub_map : data) {
        int runid = sub_map.first;
        cout << " run: " << runid << endl;
        for (auto& p : sub_map.second) cout << p.first << " " << p.second << endl;
    }

    ofstream f_out;
    f_out.open(outFile.c_str());
    f_out << Form("%10s %10s %10s","runId","triggerId","n_events") << endl;

    for (auto& sub_map : data) {
        int runid = sub_map.first;
        /* cout << " run: " << runid << endl; */
        for (auto& p : sub_map.second) {
            f_out << Form("%10i %10i %10i",runid, p.first, p.second) << endl;
        }
    }

        /* double Vz    { header->GetPrimaryVertexZ() }; */
        /* double vzVpd { header->GetvpdVz() }; */
        /* int    runid { header->GetRunId() }; */


        /* bool is_good_run { binary_search(good_runs.begin(), good_runs.end(), runid) }; */
        /* bool is_isaac_skip_run { */
        /*      (header->GetRunId() >= 16142059 && header->GetRunId() <= 16149001) */
        /*   || header->GetRunId() == 16135031 */ 
        /*   || header->GetRunId() == 16135032 }; */

        /* nevents.Fill(0); */
        /* if (!is_good_run) nevents.Fill(3.); */
        /* else nevents.Fill( is_isaac_skip_run ? 1. : 2. ); */

    /* hg_Et_all  .Write(); */
    /* hg_Et_bad  .Write(); */
    /* hg_Et_good .Write(); */

    /* pr_Et_all  .Write(); */
    /* pr_Et_bad  .Write(); */
    /* pr_Et_good .Write(); */
    /* nevents    .Write(); */

    /* pAuFile->Save(); */
    /* pAuFile->Close(); */

    return 0;
}
