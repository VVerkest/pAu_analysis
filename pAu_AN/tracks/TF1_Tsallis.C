#include "TF1.h"
#include "../loc_lib/pAu_bins.h"
// requires .x tu_loadlibs.C to have been run first
TF1* TF1_Tsallis(const char* sys, const char* particle, const char* f_name="Tsallis_params.txt") {
    auto params = TsallisParams[Form("Tsallis_%s_%s",sys,particle)];
    double m0 = params[0];
    TF1* fn = new TF1( noiUniqueName(),
                [m0](double*x,double*p){
                return p[0] / TMath::Power( 
                        (1+(TMath::Sqrt(x[0]*x[0]+m0*m0)-m0)/(p[1]*p[2])), 
                        p[2]);},
                0.,15.,3);
    fn->SetParameters(params[1], params[2], params[3]);
    fn->SetParNames("A","T","n");
    return fn;
};
struct TsallisAll {
/* TF1* TF1_AllCh(const char* sys, const char* f_name="Tsallis_params.txt") { */
    array<TF1*,6> fn6;
    TsallisAll(const char* sys) {
        for (unsigned int i=0; i<6; ++i) {
            fn6[i] = TF1_Tsallis(sys, noi_geant05_ascii(i));
        }
    };
    double integral(int which_part, double x0, double x1) {
        return fn6[which_part]->Integral(x0,x1);
    };
    double integral_sum6(double x0, double x1) {
        double sum = 0;
        for (int i=0; i<6; ++i) {
            sum += fn6[i]->Integral(x0,x1);
        }
        return sum;
    };
    double int_ratio(int which_part, double x0, double x1) {
        return integral(which_part,x0,x1) / integral_sum6(x0, x1);
    };
};
