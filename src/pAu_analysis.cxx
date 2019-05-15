//  pAu_analysis.cxx
//  Veronica Verkest       May 14, 2019


using namespace fastjet;
using namespace std;

void InitReader;


int main () {


  
  TChain* Chain = new TChain( "JetTree" );         Chain->Add( "pAu_2015_200_MB_156_160_2.root" );
  TStarJetPicoReader Reader;
  InitReader( Reader, Chain, -1 );
  TStarJetPicoEventHeader* header;    TStarJetPicoEvent* event;    TStarJetVector* sv;
  TStarJetVectorContainer<TStarJetVector> * container;



  

  return 0;
}
