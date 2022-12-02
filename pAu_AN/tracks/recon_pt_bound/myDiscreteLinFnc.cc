#include "myDiscreteLinFnc.h"
#include <vector>
#include "TF1.h"
#include <iostream>

using std::cout;
using std::endl;

myDiscreteLinFnc::myDiscreteLinFnc (vector<double> _bins, vector<double> _b, vector<double> _m) :
    bins {_bins},  b{_b}, m{_m}, nbins {(int) bins.size()-1}
{};
TF1* myDiscreteLinFnc::makeTF1 (const char* name) {
    cout << "range: " << bins[0] << " -> " << bins[nbins] << endl;
    return new TF1(name, *this, bins[0], bins[nbins], nbins);
};
double myDiscreteLinFnc::operator()(double x) {
    if (x==0) x = 1.e-6;
    auto ptr = std::upper_bound(bins.begin(), bins.end(), x);
    if (ptr==bins.end()) throw std::runtime_error("FATAL: out of bounds in in fit_bounds.h::myDiscreteLinFnc");
    auto ibin = (int)(ptr-bins.begin())-1;
    if (ibin == -1) throw std::runtime_error("FATAL: out of bounds in in fit_bounds.h::myDiscreteLinFnc");
    return m[ibin]*x + b[ibin];
};
double myDiscreteLinFnc::operator()(double* x, double *p) {
    if (x[0]==0) x[0] = 1.e-6;
    auto ptr = std::upper_bound(bins.begin(), bins.end(), x[0]);
    if (ptr==bins.end()) throw std::runtime_error("FATAL: out of bounds in in fit_bounds.h::myDiscreteLinFnc");
    auto ibin = (int)(ptr-bins.begin())-1;
    if (ibin == -1) throw std::runtime_error("FATAL: out of bounds in in fit_bounds.h::myDiscreteLinFnc");
    return m[ibin]*x[0] + b[ibin];
};
