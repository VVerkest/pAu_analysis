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

struct adc_entry {
    int   adc_0     {0}; 
    float energy_0  {0.};
    int   entries   {0};

    bool same_calc_ped   {true};
    bool same_calc_calib {true};
    bool same            {true}; // can be checked for adc constant values

    float calc_ped   {0.};
    float calc_calib {0.};

    adc_entry() {};
    adc_entry(int adc, float energy) 
    {
        if (adc != 0 || energy != 0) {
            adc_0 = adc;
            energy_0 = energy;
            entries = 1;
        }
    };

    void operator()(int adc, float energy) {
        if (adc == 0 || energy == 0) return;
        ++entries;
        if (adc_0 == 0 || energy_0 == 0) { 
            adc_0 = adc; 
            energy_0 = energy; 
            return;
        }
        if (adc_0 == adc) {
            if (energy != 0 && energy_0 != 0 && energy_0 != energy) {
            cout << " ADC same false ";
            cout <<  "acd: " << adc_0 << "<>" << adc 
                 << "  energy:"<<energy_0 <<"<>"<<energy << endl;
                same = false;
            } 
            return;
        } 
        float _calc_calib = (energy_0-energy) / (adc_0-adc);
        float _calc_ped   = adc-energy/_calc_calib;
        if (calc_calib == 0) {
            calc_calib = _calc_calib;
            calc_ped   = _calc_ped;
            return;
        }
        if (same_calc_calib && TMath::Abs(_calc_calib-calc_calib)>0.001) { 
            cout << " DELTA calc_calc: " << calc_calib << "<>" << _calc_calib  << endl;
            same_calc_calib = false; same = false; 
        }
        if (same_calc_ped   && TMath::Abs((_calc_ped-calc_ped)>0.001  )) { 
            cout << " DELTA calc_ped: " << calc_ped << "<>" << _calc_ped  << endl;
            same_calc_ped = false; 
            same = false; 
        }
        return;
    };
};

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
    } else { cerr<< "incorrect number of command line arguments"; return 5000; }


    TStarJetPicoEventHeader* header;
    TStarJetPicoEvent* event;

    // write to an output text file:
    // runID  triggerID  n_events
    
    // there will be a variable number of line entries per run
    array<map<int,adc_entry>,274> tower_data;
    vector<int> runids;
    ifstream ifile;
    ifile.open("lists/HT_runid.list");
    string line;
    // skip the first line which is a comment line
    while (getline(ifile,line)) {
        istringstream iline { line.c_str() };
        int runid;
        iline >> runid;
        runids.push_back(runid);
    } // note that good_runs are numerically sorted already
    cout << " size: " << runids.size() << endl;


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
    Reader.SetApplyFractionHadronicCorrection(false);
    TStarJetPicoTowerCuts* towerCuts = Reader.GetTowerCuts();
    towerCuts->AddBadTowers("lists/dummy_tower.list");
    towerCuts->SetMaxEtCut( 1.0 );
    towerCuts->SetMaxEtCut( 1.0 );
    Reader.SetApplyFractionHadronicCorrection(false);
    /* cout << " a0 " << endl; */
    /* int numEvents = number_of_events;        // total events in HT: 152,007,032 */
    InitReader( Reader, Chain, number_of_events );

    auto _data_begin = runids.begin();
    while ( Reader.NextEvent() ) {
    /* cout << " a1 " << endl; */
        Reader.PrintStatus(10);
        event = Reader.GetEvent();
        header = event->GetHeader();
        int _runid { header->GetRunId() };
        /* auto _loc = lower_bound(runids.begin(), runids.end(), _runid); */
        /* int pos = _data_begin */ 
        
        int pos  = (int)(lower_bound(runids.begin(),runids.end(),_runid)-_data_begin);
        /* cout << " pos: " << pos << endl; */
        auto& data = tower_data[pos];

        /* cout << " runid: " << _runid << endl; */
        // loop through the towers and wright them to data
        for (int i{0}; i<header->GetNOfTowers(); ++i) {
            TStarJetPicoTower* tower { event->GetTower(i) };

            int    id { tower->GetId() };
            float  energy { tower->GetEnergy() };
            int    adc{ tower->GetADC() };
            if (data.find(id)==data.end()) {
                data[id] = {adc, energy};
            } else {
                data[id](adc,energy);
            }
        }
    }

    // write out the results
    TFile fout { outFile.c_str(), "recreate" };
    TH2F ped  {"ped",  "ped;loc. run-id;Tower ID", 274, -0.5, 273.5, 4800, 0.5, 4800.5 };
    TH2F gain {"gain", "ped;loc. run-id;Tower ID", 274, -0.5, 273.5, 4800, 0.5, 4800.5 };

    for (int i{0}; i<274; ++i) {
        float runid = i;
        for (auto& tow : tower_data[i]) {
            float towerId = tow.first;
            ped.Fill(runid, towerId, tow.second.same_calc_ped   ? tow.second.calc_ped :  -100.);
            gain.Fill(runid,towerId, tow.second.same_calc_calib ? tow.second.calc_calib :-100.);
        };
    };
    ped.Write();
    gain.Write();
    fout.Save();

    return 0;
}
