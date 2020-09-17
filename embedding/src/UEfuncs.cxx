#include "params.hh"
#include "UEfuncs.hh"

namespace Analysis {

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

  TString GetEfficHistoName( TString eaString ) {
    if ( eaString=="loEA" ) { return "2_3"; }
    else if ( eaString=="hiEA" ) { return "8_10"; }
    else { std::cerr<<"EA STRING '"<<eaString<<"' IS UNACCEPTABLE!"<<std::endl; return eaString; }
  }
}
