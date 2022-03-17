// tower_count.cxx

#include "pAuFunctions.hh"
using namespace std;
using namespace pAuAnalysis;

int main ( int argc, const char** argv ) {

    if (true) {
        // test functionality of IntListSet
        IntListSet list ("lists/test.list", true);
        for (auto i : vector<int>{1,2,3,5,6,10}) {
            cout << " i: " << i << "  in list? " << list(i) << endl;
        }
    }
    return 1;
};

