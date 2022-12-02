root -l <<EOF
#include "tu_fmt_loc.h"
#include "tuPads_loc.h"

    tuPadDim d {.1,.2};
    /* tuPads pads{3}; */

    /* TH1D hg {"hg","a;b;c",10,0.,10.}; */
    /* TRandom3 _rand; */
    /* for (int i=0;i<100;++i) { */
    /*     hg.Fill(_rand.Uniform(0.,10.)); */
    /* } */
    /* tu_fmt(&hg,{{"MarkerStyle",kFullSquare,"LineColorAlpha",kBlue,1.,"FillColorAlpha",kBlue,0.9,"FillColor",kGreen+2}}); */
    /* /1* hg.SetMarkerColorAlpha(kRed,0.3); *1/ */

    /* pads(0); */
    /* hg.Draw("PE3"); */
    /* pads(1); */
    /* hg.Draw("PE2"); */
    /* pads(2); */
    /* hg.Draw("PE1"); */


EOF
