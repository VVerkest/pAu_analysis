#ifndef myDiscreteLinFnc_h
#define myDiscreteLinFnc_h

#include <vector>
#include "TF1.h"
using std::vector;

struct myDiscreteLinFnc {
    vector<double> bins, b, m;
    int nbins;
    myDiscreteLinFnc (vector<double> _bins, vector<double> _b, vector<double> _m);
    TF1* makeTF1 (const char* name);
    double operator()(double x);
    double operator()(double* x, double *p);
};

#endif
