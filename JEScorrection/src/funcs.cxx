// Veronica Verkest
// October 13, 2020

#include "params.hh"
#include "funcs.hh"

typedef fastjet::contrib::SoftDrop SD;

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


  double CalculateGeometricMean( TH1D* h ) {
    double geoMean = 0;
    double totalContent = 0;
    for (int i=1; i<=h->GetNbinsX(); ++i) {
      geoMean += h->GetBinCenter(i)*h->GetBinContent(i)*h->GetBinWidth(i);
      totalContent += h->GetBinContent(i)*h->GetBinWidth(i);
    }
    geoMean /= totalContent;
    return geoMean;
  }


  TH1D* GenerateFractionalContribution( TH2D* response, double pTlo, double pThi, TString dir, TString lohiEA ) {
    int lo = (int)pTlo;     int hi = (int)pThi;

    TCanvas *can = new TCanvas();
    can->SetLogy();

    response->GetXaxis()->SetRangeUser(pTlo,pThi);
    TString name = "hResponse_" + lohiEA + "EA_"; name+=lo; name+="_"; name+=hi;
    TH1D *response1D = (TH1D*)response->ProjectionY(name);
    response1D->Scale(1./response1D->Integral());
    TString title = "frac. contribution to a "; title+=lo; title+="-"; title+=hi; title+=" GeV part. jet";
    response1D->GetYaxis()->SetTitle(title);
    response1D->Draw();
    // response1D->Write();
    TString saveName = dir + "/response_" + lohiEA + "EA_"; saveName+=lo; saveName+="_"; saveName+=hi; saveName+=".pdf";
    can->SaveAs(saveName,"PDF");
    response->GetXaxis()->SetRange(1, -1);

    can->Destructor();
  
    return response1D;
  }


  void GenerateWeightedPtResponse( TH1D *h_DetWt[3], TH1D *h_Det[21], TH1D *h_MissedProb, TString plot_dir ) {

    for (int p=0; p<nPtBins; ++p) {
      TString name = "hPtResponse"; name += ptBinName[p];
      h_DetWt[p] = new TH1D( name, ";det-level leading jet p_{T} (GeV)",55,4.0,59.0);
    }

    TCanvas * can5 = new TCanvas( "can5" , "" ,700 ,500 );              // CANVAS 5

    for (int p=0; p<nPtBins; ++p) {  //  WEIGHT AND ADD PROJECTIONS TO GET WEIGHTED DETECTOR RESPONSE
      for (int i=0; i<21; ++i) {
	int ptVal = i + 5;
	int binno = i + 1;
	if ( ptVal >= ptLo[p] && ptVal <= ptHi[p] ) {
	  // std::cout<<h_MissedProb->GetBinContent(binno)<<std::endl;
	  double wt = 1+h_MissedProb->GetBinContent(binno);  // weight according to misses at part-level
	  h_DetWt[p]->Add(h_Det[i],wt);
	}
      }    
      h_DetWt[p]->Scale(1./h_DetWt[p]->Integral());
      h_DetWt[p]->GetYaxis()->SetRangeUser(0.0,0.2);
      h_DetWt[p]->SetMarkerStyle(ptMarker[p]);
      h_DetWt[p]->SetStats(0);
      h_DetWt[p]->Draw("SAME PLC PMC");
    }
    can5->BuildLegend(0.4,0.68,0.9,0.9);
    TString saveName = plot_dir + "WeightedPtResponse.pdf";
    can5->SaveAs(saveName,"PDF");

    for (int p=0; p<nPtBins; ++p) { h_DetWt[p]->SetStats(1); }
    
  }


  void GetEmbeddingHistograms( TFile *embedFile, TH2D *h_PtResponse, TH1D *h_FakeJets, TH1D *h_MissedJets, TString plot_dir ) {
    TH2D *hResponse[nEtaBins];
    TH1D *hFakes[nEtaBins];
    h_MissedJets = (TH1D*)embedFile->Get("hMisses");

    int nAccepted = 0; int nFakes = 0;
  
    for (int e=0; e<nEtaBins; ++e) {    
      TString name = "hPtResponse" + etaBinName[e] + "Jet";
      hResponse[e] = (TH2D*)embedFile->Get(name);
      name = "hFakes" + etaBinName[e] + "Jet";
      hFakes[e] = (TH1D*)embedFile->Get(name);

      nAccepted += hResponse[e]->GetEntries();
      nFakes += hFakes[e]->GetEntries();
    }

    for (int e=0; e<nEtaBins; ++e) {
      h_PtResponse->Add(hResponse[e]);
      h_FakeJets->Add(hFakes[e]);
    }

    TString title = "Fakes;fake det-level leading jet p_{T} (GeV);prob.";
    h_FakeJets->SetTitle(title);
    
    TCanvas * can2 = new TCanvas( "can2" , "" ,700 ,500 );              // CANVAS 2
    can2->SetLogy();

    int nMissed = h_MissedJets->GetEntries();
    int nEvents = nAccepted + nMissed + nFakes;

    double scale = (double)nMissed/nEvents;
    title = h_MissedJets->GetTitle(); title += ";prob.";
    h_MissedJets->SetTitle(title);
    h_MissedJets->Scale(scale/h_MissedJets->Integral());
    h_MissedJets->Draw();
    TString saveName = plot_dir + "MissedJets.pdf";
    can2->SaveAs(saveName,"PDF");
    
    h_FakeJets->Scale(scale/h_FakeJets->Integral());
    h_FakeJets->Draw();
    saveName = plot_dir + "fakeJets.pdf";
    can2->SaveAs(saveName,"PDF");
    
    TCanvas * can3 = new TCanvas( "can3" , "" ,700 ,500 );              // CANVAS 3
    can3->SetLogz();
    h_PtResponse->Draw("COLZ");
    saveName = plot_dir + "pTresponse.pdf";
    can3->SaveAs(saveName,"PDF");  // SAVE 2D PT RESPONSE

    can2->cd();
    TH1D *hDetPt = new TH1D("hDetPt",";det-level leading jet p_{T}",55,4.0,59.0);
    hDetPt = (TH1D*)h_PtResponse->ProjectionY();
    hDetPt->Scale(1./hDetPt->Integral());
    // h_FakeProb->Add(h_FakeJets);
    // h_FakeProb->Divide(hDetPt);
    // TString name = "hFakeProb";
    // title = ";missing part-level leading jet p_{T} (GeV);prob.";
    // h_FakeProb->SetNameTitle(name, title);

    // h_FakeProb->Draw();
    // saveName = plot_dir + "EmbeddingFakeJetProb.pdf";
    // can2->SaveAs(saveName,"PDF");

    
    TH1D *hPartPt = new TH1D("hPartPt",";part-level leading jet p_{T}",55,4.0,59.0);
    hPartPt = (TH1D*)h_PtResponse->ProjectionX();
    hPartPt->Scale(1./hPartPt->Integral());
    // h_MissedProb->Add(h_MissedJets);
    // h_MissedProb->Divide(hPartPt);
    // name = "hMissProb";
    // title = ";missing part-level leading jet p_{T} (GeV);prob.";
    // h_MissedProb->SetNameTitle( name, title );

    // h_MissedProb->Draw();
    // saveName = plot_dir + "EmbeddingMissedJetProb.pdf";
    // can2->SaveAs(saveName,"PDF");

    can2->SetLogy();
    hPartPt->Draw();
    saveName = plot_dir + "EmbeddingPartLevelPt.pdf";
    can2->SaveAs(saveName,"PDF");

    hDetPt->Draw();
    saveName = plot_dir + "EmbeddingDetLevelPt.pdf";
    can2->SaveAs(saveName,"PDF");


    can2->Destructor();
    can3->Destructor();
  }


  void ProjectAndPlotByEta( TH2D *h_AddedChgUE2D_corr, TH1D *h_PartJetUE, int eval, int pval, int aval, TString dir_name, TFile *outfile ){
    // h_AddedChgUE2D_corr->SetAxisRange(etaLo[eval],etaHi[eval],"Y");
    h_AddedChgUE2D_corr->GetYaxis()->SetRangeUser(etaLo[eval],etaHi[eval]);
    h_PartJetUE = (TH1D*)h_AddedChgUE2D_corr->ProjectionX();
    TString name = "hPartJetUE"; name += ptBinName[pval] + etaBinName[eval] + "_" + lohi[aval];
    h_PartJetUE->SetName(name);
    double mean = 0.;

    // std::cout<<h_PartJetUE->Integral("width")/area[eval]<<std::endl;
    
    for (int i=0; i<h_PartJetUE->GetNbinsX(); ++i) {
      int binno = i+1;

      if (h_PartJetUE->GetBinContent( binno )==0){ continue; }
      
      mean += h_PartJetUE->GetBinContent( binno )*h_PartJetUE->GetBinCenter( binno );///h_PartJetUE->Integral("width");

      double width = h_PartJetUE->GetXaxis()->GetBinWidth(binno);
      h_PartJetUE->SetBinContent( binno, h_PartJetUE->GetBinContent(binno)/width );
    }

    outfile->cd();
    h_PartJetUE->Write();

    std::cout<<mean<<std::endl;
    std::cout<<mean/h_PartJetUE->Integral("width")<<std::endl<<std::endl;
    mean /= h_PartJetUE->Integral("width");

    mean = h_PartJetUE->GetMean(1);
    
    TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
    can->SetLogy();
    gStyle->SetOptStat(0);
    h_PartJetUE->SetStats(kFALSE);
    h_PartJetUE->SetAxisRange(0.000001,4.,"Y");
    h_PartJetUE->Draw();

    TString text = "#LT p_{T}^{ch}#GT = "; text += mean;//h_PartJetUE->GetMean(1);
    text = text(0,26);
    DrawText(text, 0.6, 0.7, 20);
    
    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += h_PartJetUE->Integral("width")/area[eval];
    // text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += h_PartJetUE->Integral()/area[eval];
    text = text(0,42);
    DrawText( text, 0.6, 0.55, 20 );

    name = dir_name + "correctedUEpT" + ptBinName[pval] + etaBinName[eval] + ".pdf";
    can->SaveAs(name, "PDF");
  
    can->Destructor();
  }


  void ProjectAndSaveFinalUEPlots( TH1D *h_WtUEpt, TString suffix, double area, TString plot_dir ){

    TString canvasName = "can8" + suffix;
    TCanvas * can8 = new TCanvas( canvasName , "" ,700 ,500 );              // CANVAS 8
    can8->SetLogy();
	
    double pt, effic, corrPt, corrErr;
    int ptBin;

    h_WtUEpt->SetStats(1);
    h_WtUEpt->SetAxisRange( 0.000001,10,"Y");
    h_WtUEpt->SetLineColor(kBlack);
    h_WtUEpt->SetMarkerColor(kBlack);
    h_WtUEpt->Draw();
	
    TString text = "#LT p_{T}^{ch}#GT = "; text += h_WtUEpt->GetMean(1);
    text = text(0,26);
    DrawText(text, 0.6, 0.7, 20);

    // std::cout<<h_WtUEpt->Integral()<<std::endl;
    
    text = "#LT#frac{dN_{ch}}{d#eta d#phi}#GT = "; text += h_WtUEpt->Integral()/area;
    text = text(0,42);
    DrawText( text, 0.6, 0.55, 20 );
      
    TString name = plot_dir + "CorrectedWtUEpt"; name += suffix; name += ".pdf";
    can8->SaveAs(name,"PDF");

    can8->Destructor();
  }

  

  void ProjectPartLevelJetPt( TH2D* h_PtResponse, TH1D *h_Det[21], TString plot_dir ) {

    
    for (int i=0; i<21; ++i) {  // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
      int ptVal = i + 10;
      int binno = i + 6;

      // std::cout<<h_PtResponse->GetXaxis()->GetBinCenter(binno)<<std::endl;
    
      TString name = "hPtResponse_"; name += ptVal; name += "GeV";
      h_Det[i] = (TH1D*) h_PtResponse->ProjectionY(name,binno,binno);
      h_Det[i]->Scale(1./h_Det[i]->Integral());
      h_Det[i]->SetMarkerStyle(marker[i]);
      h_Det[i]->SetMarkerStyle(marker[i]);
      h_Det[i]->SetMarkerSize(1);
      h_Det[i]->SetStats(0);
      name = ""; name += ptVal; name += " GeV part. jet";
      h_Det[i]->SetNameTitle(name,";det-level leading jet p_{T} (GeV)");
    }


    TCanvas * can4 = new TCanvas( "can4" , "" ,700 ,500 );              // CANVAS 4
    for (int i=0; i<21; ++i) { h_Det[i]->Draw("SAME PLC PMC"); } // PROJECT FOR ALL PART-LEVEL BINS 10-30 GeV
    can4->BuildLegend(0.68,0.1,0.9,0.9);
    TString saveName = plot_dir + "detPtResponses.pdf";
    can4->SaveAs(saveName,"PDF");
  
  }


  void ProjectScaleAndSaveUE1D( TH2D* UE2D, TH1D &UE1D, TH1D* lead, TH1D &UE1D_dbw, double pTlo, double pThi, TString dir, TString lohiEA ){

    TString name, saveName;
  
    TCanvas *can = new TCanvas();
    can->SetLogy();
  
    int lo = (int)pTlo;     int hi = (int)pThi;
    
    //  UEpT_varBinning_pTlo_pThi
    //  UE pT distribution (scaled by Njets) for det-level jets pTlo-pThi GeV
    UE2D->GetXaxis()->SetRangeUser( pTlo, pThi );
    name = "hUE1D_" + lohiEA + "EA_"; name+=pTlo; name+="_"; name+=pThi; name +="GeV";
    TH1D* htemp = (TH1D*)UE2D->ProjectionY(name);
    UE1D = *htemp;
    int binlo = lo-3; int binhi = hi-4;
    UE1D.Scale(1./lead->Integral(lo-3,hi-4));
    UE1D.GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
    if (UE1D.GetEntries()>0) {
      UE1D.SetAxisRange(0.0000001,1.,"Y");
      UE1D.Draw();
      saveName = dir + "/UEpT_varBinning_" + lohiEA + "EA_"; saveName+=pTlo; saveName+="_"; saveName+=pThi; saveName +=".pdf";
      can->SaveAs(saveName,"PDF");
      std::cout<<"mean: "<<UE1D.GetMean()<<std::endl;
      std::cout<<"integral: "<<UE1D.Integral()<<std::endl<<std::endl;
    }

    //  UEpT_pTlo_pThi
    //  UE pT distribution (scaled by Njets, adjusted for bin width, and re-scaled to preserve integral) for det-level jets pTlo-pThi GeV
    name = "hUE1D_" + lohiEA + "EA_"; name+=pTlo; name+="_"; name+=pThi; name +="_dbw";

    UE1D_dbw = *htemp;
    binlo = lo-3; binhi = hi-4;
    UE1D_dbw.Scale(1./lead->Integral(lo-3,hi-4));
    UE1D_dbw.GetYaxis()->SetTitle("#frac{1}{N_{jets}} #frac{dN_{UE}}{dp_{T}^{UE}}");
  
    for (int i=1; i<=UE1D_dbw.GetNbinsX(); ++i) { UE1D_dbw.SetBinContent(i,UE1D_dbw.GetBinContent(i)/UE1D_dbw.GetBinWidth(i)); }
    UE1D_dbw.Scale(UE1D.Integral()/UE1D_dbw.Integral());
    if (UE1D_dbw.GetEntries()>0) {
      std::cout<<"geometric mean: "<< CalculateGeometricMean(&UE1D_dbw) <<std::endl;
      std::cout<<"integral: "<<UE1D_dbw.Integral()<<std::endl<<std::endl;
      UE1D_dbw.SetAxisRange(0.0000001,1.,"Y");
      UE1D_dbw.Draw();
      saveName = dir + "/UEpT_" + lohiEA + "EA_"; saveName+=pTlo; saveName+="_"; saveName+=pThi; saveName +="_dbw.pdf";
      can->SaveAs(saveName,"PDF");
    }

    // UE1D.Write();
    // UE1D_dbw.Write();
    can->Destructor();
  }

  
  TString RoundDecimal(double val, int nDigitsAfterDecimal ) {
    double n = val;
    n += 5/pow(10,nDigitsAfterDecimal+1);
    TString rounded = "";
    rounded += n;
    int nDigitsBeforeDecimal = 1 + rounded.First(".");  
    int input = nDigitsBeforeDecimal + nDigitsAfterDecimal;
    rounded = rounded(0,input);
    return rounded;
  }



  // void StackAndSaveNchPlots( TH1D *hPt[3][3], TH2D *hscale, TString save_name, TString dir_name ) {
  //   auto hs = new THStack("hs","Stacked 1D histograms colored using kOcean palette");

  //   const int etaColor[nEtaBins] = { 877, 596, 814 };
  //   int EAmarker[2] = { 24, 20 };

  //   const int n_bins = 3;
  //   double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
  //   double shiftedBins[n_bins+1];

  //   TH1D *hNch[nEtaBins]; 

  //   for (int e=0; e<nEtaBins; ++e) {
  //     hscale->GetXaxis()->SetBinLabel(e+1,ptBinString[e]);
  //     for (int i=0; i<=n_bins; ++i) { shiftedBins[i] = bin_edge[i] + 0.3*(e-1); }
  //     TString name = "hNch" + etaBinName[e];
  //     hNch[e] = new TH1D(name,";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,shiftedBins);
  //     hNch[e]->SetLineColor(etaColor[e]);
  //     hNch[e]->SetMarkerColor(etaColor[e]);
  //     hNch[e]->SetMarkerStyle(EAmarker[1]);
  //     hNch[e]->SetMarkerColor(etaColor[e]);
  //     for (int p=0; p<nPtBins; ++p) {
  // 	hPt[e][p]->SetMarkerColor(etaColor[e]);
  // 	hPt[e][p]->SetLineColor(etaColor[e]);
  // 	hPt[e][p]->SetMarkerStyle(EAmarker[1]);
  // 	hPt[e][p]->SetMarkerColor(etaColor[e]);
  // 	hNch[e]->SetBinContent(p+1,hPt[e][p]->Integral()/area[e]);
  // 	hNch[e]->SetBinError(p+1,hPt[e][p]->GetMeanError(1));
  //     }
  //     hNch[e]->GetYaxis()->SetRangeUser(0.5,1.8);
  //     hs->Add(hNch[e]);
  //   }
  //   TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );

  //   hscale->SetStats(0);
  //   hscale->Draw();
  //   hs->Draw("SAMEnostackEX0");   // draw the stack
  //   TString name = dir_name + save_name;
  //   can->SaveAs( name, "PDF");

  //   can->Destructor();
  //   // hNch->Destructor();
  // }



  void StackAndSavePtPlots( TH1D *hPt[3][3], TH2D *hscale, TString save_name, TString dir_name ) {
    auto hs = new THStack("hs","Stacked 1D histograms colored using kOcean palette");

    const int etaColor[nEtaBins] = { 877, 596, 814 };
    int EAmarker[2] = { 24, 20 };

    const int n_bins = 3;
    double bin_edge[n_bins+1] = { 10.0, 15.0, 20.0, 30.0 };
    double shiftedBins[n_bins+1];

    TH1D *hMeanPt[nEtaBins]; 

    for (int e=0; e<nEtaBins; ++e) {
      hscale->GetXaxis()->SetBinLabel(e+1,ptBinString[e]);
      for (int i=0; i<=n_bins; ++i) { shiftedBins[i] = bin_edge[i] + 0.3*(e-1); }
      TString name = "hMeanPt" + etaBinName[e];
      hMeanPt[e] = new TH1D(name,";leading jet p_{T} (GeV);#LT p_{T}^{ch}#GT (GeV)",n_bins,shiftedBins);
      hMeanPt[e]->SetAxisRange(0.55,0.85,"Y");
      hMeanPt[e]->SetLineColor(etaColor[e]);
      hMeanPt[e]->SetMarkerColor(etaColor[e]);
      hMeanPt[e]->SetMarkerStyle(EAmarker[1]);
      hMeanPt[e]->SetMarkerColor(etaColor[e]);
      for (int p=0; p<nPtBins; ++p) {
	hPt[e][p]->SetMarkerColor(etaColor[e]);
	hPt[e][p]->SetLineColor(etaColor[e]);
	hPt[e][p]->SetMarkerStyle(EAmarker[1]);
	hPt[e][p]->SetMarkerColor(etaColor[e]);
	hMeanPt[e]->SetBinContent(p+1,hPt[e][p]->GetMean(1));
	hMeanPt[e]->SetBinError(p+1,hPt[e][p]->GetMeanError(1));
      }
      hMeanPt[e]->SetAxisRange(0.55,0.85,"Y");
      hMeanPt[e]->GetYaxis()->SetRangeUser(0.55,0.85);
      hs->Add(hMeanPt[e]);
    }
    TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );

    hscale->SetStats(0);
    hscale->Draw();
    hs->Draw("SAMEnostackEX0");   // draw the stack
    TString name = dir_name + save_name;
    can->SaveAs( name, "PDF");

    can->Destructor();
    // hMeanPt->Destructor();
  }


  void TrackingEfficiencyByPtAndEta( TH2D* h_ChgUE2D[nPtBins], TH2D* h_ChgUE2D_corr[nPtBins], TFile *effic_File, TString ea_string, TString dir_name ){
    //  X: UE pT,   Y: UE eta
    TH1D *hEffic;
    TString name, saveName, bbcBins;
    double ptVal, oldVal, oldErr, effic, efficErr, relErr, newVal, newErr;

    TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
    can->SetLogz();
    gStyle->SetOptFit(1);

    TF1 *eff = new TF1("eff","((200.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
    // TF1 *eff = new TF1("eff","(([2]+log(x/100.0)))*exp([0]+[1]*x)",0.2,15.0);
    // TF1 *eff = new TF1("eff","(1+log(x/200.))*(exp([0]+[1]*x))+[2]",0.2,15.0);
    // TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  
    if ( ea_string=="lo" ) { bbcBins = "1_3"; }
    else if ( ea_string=="hi" ) { bbcBins = "8_10"; }
    else { std::cerr<< "invalid EA string provided!" <<std::endl; }
  
    for (int iy=1; iy<zbins+1; ++iy){ // loop over UE eta bins
    
      name = "eff_s_bin_" + bbcBins + "_bbc__"; name += iy; name += "_"; name += iy; name += "_eta";
      TH1D *hEffic = (TH1D*)effic_File->Get(name);
      
      hEffic->Fit( "eff", "EMR" );
      TF1* efficFit = (TF1*)hEffic->GetFunction("eff");

      hEffic->Draw();
    
      saveName = dir_name + name + ".pdf";
      can->SaveAs(saveName, "PDF");
      std::cout<<std::endl<<std::endl;
    
      for (int ix=1; ix<ybins+1; ++ix){ // loop over UE pt bins

	for (int i=0; i<nPtBins; ++i) {

	  ptVal = h_ChgUE2D[i]->GetXaxis()->GetBinCenter( ix );
	  effic = efficFit->Eval( ptVal );
	  // effic = hEffic->GetBinContent( ix );
	  // if (ptVal>=3.) { effic = hEffic->GetBinContent( hEffic->FindBin(3.) ); }

	  oldVal = h_ChgUE2D[i]->GetBinContent(ix,iy);
	  if ( oldVal==0.0 ) { continue; }
	  oldErr = h_ChgUE2D[i]->GetBinError(ix,iy);
	  
	  newVal = oldVal/effic;
	  efficErr = hEffic->GetBinError( hEffic->FindBin( ptVal) );
	  relErr = sqrt( oldErr*oldErr + efficErr*efficErr );
	  newErr = newVal*relErr;
      
	  h_ChgUE2D_corr[i]->SetBinContent(ix,iy,newVal);
	  h_ChgUE2D_corr[i]->SetBinError(ix,iy,newErr);
	}
      }
    }

    for (int i=0; i<nPtBins; ++i) {
      h_ChgUE2D_corr[i]->SetEntries(h_ChgUE2D[i]->GetEntries());
      if (h_ChgUE2D_corr[i]->GetEntries()==0){ continue; }
      h_ChgUE2D_corr[i]->Draw("COLZ");
      name = dir_name + h_ChgUE2D_corr[i]->GetName() + ".pdf";
      can->SaveAs(name, "PDF");

      // std::cout<<h_ChgUE2D_corr[i]->Integral()/AREA<<std::endl;

    }
    can->Destructor();
  }



  //    TrackingEfficiencyByPtAndEta55( hUE2D[a], hUE2D_detCorr[a], efficFile, lohi[a], name, efficShift );
  void TrackingEfficiencyByPtAndEta55( TH2D* h_ChgUE2D[55], TH2D* h_ChgUE2D_corr[55], TFile *effic_File, TString ea_string, TString dir_name, double EfficShift = 0. ){
    //  X: UE pT,   Y: UE eta
    TH1D *hEffic;
    TString name, saveName, bbcBins;
    double ptVal, oldVal, oldErr, effic, efficErr, relErr, newVal, newErr;

    TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
    can->SetLogz();
    gStyle->SetOptFit(1);

    // TF1 *eff = new TF1("eff","((200.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
    // TF1 *eff = new TF1("eff","(([2]+log(x/100.0)))*exp([0]+[1]*x)",0.2,15.0);
    // TF1 *eff = new TF1("eff","(1+log(x/200.))*(exp([0]+[1]*x))+[2]",0.2,15.0);
    TF1 *eff = new TF1("eff","((100.+log(x/200.0))/200.)*exp([0]+[1]*x)+[2]",0.2,15.0);
  
    if ( ea_string=="lo" ) { bbcBins = "1_3"; }
    else if ( ea_string=="hi" ) { bbcBins = "8_10"; }
    else { std::cerr<< "invalid EA string provided!" <<std::endl; }
    
    for (int iy=1; iy<zbins+1; ++iy){ // loop over UE eta bins

      name = "eff_s_bin_" + bbcBins + "_bbc__"; name += iy; name += "_"; name += iy; name += "_eta";
      TH1D *hEffic = (TH1D*)effic_File->Get(name);

      // hEffic->Fit( "eff", "EMR" );
      hEffic->Fit( "eff", "","",0.2,15.0 );

      TF1* efficFit = (TF1*)hEffic->GetFunction("eff");

      hEffic->Draw();

      saveName = dir_name + name + ".pdf";
      can->SaveAs(saveName, "PDF");
    
      for (int ix=1; ix<ybins+1; ++ix){ // loop over UE pt bins

      	for (int i=0; i<55; ++i) {
      	// // for (int i=0; i<3; ++i) {
      	  ptVal = h_ChgUE2D[i]->GetXaxis()->GetBinCenter( ix );
      	  effic = efficFit->Eval( ptVal );
      	  effic += EfficShift;
      	  // std::cout<< efficFit->Eval( ptVal ) <<" + "<< EfficShift <<" = "<< effic << std::endl;

      	  if ( effic > 1.0 ) { effic = 1.0; } 
      	  if ( effic < 0.0 ) { effic = 0.0; } 
	  
      	  // effic = hEffic->GetBinContent( ix );
      	  // if (ptVal>=3.) { effic = hEffic->GetBinContent( hEffic->FindBin(3.) ); }
      	  // if (ptVal>=3.) { effic = efficFit->Eval( 3. ); }

      	  oldVal = h_ChgUE2D[i]->GetBinContent(ix,iy);
      	  if ( oldVal==0.0 ) { continue; }
      	  oldErr = h_ChgUE2D[i]->GetBinError(ix,iy);
	  
      	  newVal = oldVal/effic;
      	  // efficErr = hEffic->GetBinError( hEffic->FindBin( ptVal) );
      	  // relErr = sqrt( oldErr*oldErr + efficErr*efficErr );
      	  // newErr = newVal*relErr;

      	  h_ChgUE2D_corr[i]->SetBinContent(ix,iy,newVal);
      	  h_ChgUE2D_corr[i]->SetBinError(ix,iy,oldErr);
      	  // h_ChgUE2D_corr[i]->SetBinContent(ix,iy,oldVal);  	  h_ChgUE2D_corr[i]->SetBinError(ix,iy,oldErr);  // NO TRACKING EFFICIENCY!
      	}
      }
    }

    for (int i=0; i<55; ++i) {
    // for (int i=0; i<3; ++i) {
      // if (h_ChgUE2D[i]->GetEntries()==0){ continue; }
      h_ChgUE2D_corr[i]->SetEntries(h_ChgUE2D[i]->GetEntries());
      if (h_ChgUE2D_corr[i]->GetEntries()==0){ continue; }
      h_ChgUE2D_corr[i]->Draw("COLZ");
      name = dir_name + h_ChgUE2D_corr[i]->GetName() + ".pdf";
      can->SaveAs(name, "PDF");

      // std::cout<<h_ChgUE2D_corr[i]->GetName()<<std::endl;
      // // std::cout<<h_ChgUE2D_corr[i]->Integral()/AREA<<std::endl;
      // std::cout<<h_ChgUE2D_corr[i]->GetMean(1)<<std::endl;
      // std::cout<<h_ChgUE2D[i]->GetName()<<std::endl;
      // // std::cout<<h_ChgUE2D[i]->Integral()/AREA<<std::endl<<std::endl;
      // std::cout<<h_ChgUE2D[i]->GetMean()<<std::endl<<std::endl;

    }
    can->Destructor();
  }


  void WeightAndAddCorrected2Ds( TH2D *h_AddedChgUE2D_corr, TH1D *h_DetWt, TH2D *h_ChgUE2D_corr[55], TString dir_name ){

    TCanvas * can = new TCanvas( "can" , "" ,700 ,500 );
    can->SetLogz();

    for (int i=0; i<55; ++i) {
      int ptLo = i + 4;
      int binno = i + 1;

      double weight = h_DetWt->GetBinContent(binno);
      h_AddedChgUE2D_corr->Add( h_ChgUE2D_corr[i], weight );
    }

    h_AddedChgUE2D_corr->Draw("COLZ");
    TString name = dir_name + "added" + h_AddedChgUE2D_corr->GetName() + ".pdf";
    // std::cout<<h_AddedChgUE2D_corr->Integral()/AREA<<std::endl;
    can->SaveAs(name, "PDF");

    can->Destructor();

  }


  TH1D* WeightAndSumByFC1D( TH1D* FC, TH1D* UE1D[55] ){

    const int bins = 15;
    double binEdge[bins+1] = { 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0 };
    TH1D* UE_part = new TH1D("UE_part","",bins,binEdge);

    for (int i=0; i<FC->GetNbinsX(); ++i) {
      int binno = i+4;
      double wt = FC->GetBinContent(binno);
      if (wt==0) { continue; }
      UE_part->Add(UE1D[i], wt);
    }

    UE_part->Write();
    return UE_part;
  }


  void WeightAndSumByFC2D( TH1D* FC, TH2D* UE2D[55], TH2D *UE2Dpart ){

    double sum = 0.;

    for (int i=0; i<FC->GetNbinsX(); ++i) {
    
      int binno = i+1;
    
      double wt = FC->GetBinContent(binno);
      // std::cout<<FC->GetBinCenter(binno)<<std::endl;
      if (wt==0) { continue; }
      UE2Dpart->Add(UE2D[i], wt);

      sum += wt;

    }

    // std::cout<<sum<<std::endl;
  
  }


  void WeightAndSumByFC2D_fakesCorrection( TH1D* FC, TH1D* FakeProb, TH2D* UE2D[55], TH2D *UE2Dpart ){

    double sum = 0.;
  
    for (int i=0; i<FC->GetNbinsX(); ++i) {
      int binno = i+1;
      FC->SetBinContent( binno, FC->GetBinContent(binno)*( 1. - FakeProb->GetBinContent(1+i) ) );
    }

    FC->Scale(1./FC->Integral());

    for (int i=0; i<FC->GetNbinsX(); ++i) {
      int binno = i+1;
      double wt = FC->GetBinContent(binno);
      if (wt==0) { continue; }
      UE2Dpart->Add(UE2D[i], wt);
      sum += wt;
    }
    // std::cout<<sum<<std::endl;
  }
  

  void WeightUEPtByLeadPtAndFakes( TH1D *h_WtUEpt[nPtBins],TH1D *h_UEpt[55],TH1D *h_DetWt[nPtBins],TH1D *h_leadPt,TH1D *h_FakeJets,TString plot_dir, TString suf ){

    for (int i=0; i<55; ++i) {  // det-level fractional pT contribution to part-level pT --> fc(pT_det)
      int ptVal = i + 5;
      int binno = i + 1;
      if (h_leadPt->GetBinContent( binno )>0){
        for (int p=0; p<nPtBins; ++p) {  // ADD UE DISTRIBUTIONS WEIGHTED BY DET-LEVEL FRACTIONAL CONTRIBUTION
	  double wt = h_DetWt[p]->GetBinContent(binno)*( 1.0 - h_FakeJets->GetBinContent(binno) );
	  h_WtUEpt[p]->Add( h_UEpt[i], wt );
	}
      }
    }
  
    TCanvas * can7 = new TCanvas( "can7" , "" ,700 ,500 );              // CANVAS 7
    for (int p=0; p<nPtBins; ++p) {
      h_WtUEpt[p]->SetMarkerStyle(ptMarker[p]);
      h_WtUEpt[p]->SetStats(0);
      h_WtUEpt[p]->Draw("PLC PMC SAME");    
    }
    can7->SetLogy();
    can7->BuildLegend(0.4,0.68,0.9,0.9);
    TString saveName = plot_dir + "weightedUEptByLeadPt_" + suf + "UE.pdf";
    can7->SaveAs(saveName,"PDF");
    can7->Destructor();
  }
  

}
