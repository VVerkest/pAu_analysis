// Veronica Verkest
// January 26, 2022

void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(63);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  //tex->SetTextFont(42);
  tex->SetNDC();
  tex->Draw();
}


TH1D* GenerateFractionalContribution( TH2D* response, double pTlo, double pThi, TString dir/*, TString lohiEA*/ ) {
  int lo = (int)pTlo;     int hi = (int)pThi;

  TCanvas *can = new TCanvas();
  can->SetLogy();

  response->GetXaxis()->SetRangeUser(pTlo,pThi);
  TString name = "hResponse_"; name+=lo; name+="_"; name+=hi;
  TH1D *response1D = (TH1D*)response->ProjectionY(name);
  response1D->Scale(1./response1D->Integral());
  TString title = "frac. contribution to a "; title+=lo; title+="-"; title+=hi; title+=" GeV part. jet";
  response1D->GetYaxis()->SetTitle(title);
  response1D->Draw();
  // response1D->Write();
  TString saveName = dir + "/response_"; saveName+=lo; saveName+="_"; saveName+=hi; saveName+=".pdf";
  can->SaveAs(saveName,"PDF");
  response->GetXaxis()->SetRange(1, -1);

  can->Destructor();
  
  return response1D;
}


TH1D* WeightAndSumByFC1D( TH1D* FC, TH1D* obs[55] ){

  TH1D* obs_part = new TH1D("obs_part","",12,3000.,70000.);

  for (int i=0; i<FC->GetNbinsX(); ++i) {
    int binno = i+4;
    double wt = FC->GetBinContent(binno);
    if (wt==0) { continue; }
    obs_part->Add(obs[i], wt);
  }

  // obs_part->Write();
  return obs_part;
}


// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
// ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


void LeadJetBBCdist(){

  TH1::SetDefaultSumw2();  TH2::SetDefaultSumw2();  TH3::SetDefaultSumw2();
  
  const double pi = 3.14159265;
  const double AREA = 4*(pi - 2);   // (  2 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  const int nPtBins = 3;
  const double ptLo[nPtBins] = { 10.0, 15.0, 20.0 };
  const double ptHi[nPtBins] = { 15.0, 20.0, 30.0 };
  const TString ptBinName[nPtBins] = { "_10_15GeV", "_15_20GeV", "_20_30GeV" };
  const TString ptBinString[nPtBins] = { "10 < p_{T,lead} < 15", "15 < p_{T,lead} < 20",  "20 < p_{T,lead} < 30" };
  const TString ptBinStringCorrected[nPtBins] = { "10 < p_{T,lead}^{corrected} < 15", "15 < p_{T,lead}^{corrected} < 20",  "20 < p_{T,lead}^{corrected} < 30" };
  const int ptColor[nPtBins] = { 797, 593, 892 };
  const int ptMarker[nPtBins] = { 20, 21, 29 };
  const int ptCorrectedMarker[nPtBins] = { 24, 25, 30 };

  const int nEAbins = 2;
  TString EAbinName[nEAbins] = { "Lo", "Hi" };
  TString EAbinString[nEAbins] = { "Low EA", "High EA" };
  TString BBCselection[nEAbins] = { "BbcAdcSumEast>3559.12 && BbcAdcSumEast<11503", "BbcAdcSumEast>26718.1" };

  const double eastArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double midArea = 2*(0.6)*(pi - 2);   // (  0.6 in eta  ) X (  2*( pi-1 - 1 ) in phi  )
  const double westArea = 2*(0.7)*(pi - 2);   // (  0.7 in eta  ) X (  2*( pi-1 - 1 ) in phi  )

  const int nEtaBins = 3;
  const double etaLo[nEtaBins] = { -1.0, -0.3, 0.3 };
  const double etaHi[nEtaBins] = { -0.3, 0.3, 1.0 };
  const TString etaBinName[nEtaBins] = { "_eastEta", "_midEta", "_westEta" };
  const TString emw[nEtaBins] = { "east", "mid", "west" };
  const TString etaBinString[nEtaBins] = { "-0.6<#eta_{jet}<-0.3", "-0.3<#eta_{jet}<0.3", "0.3<#eta_{jet}<0.6" };
  const int etaColor[nEtaBins] = { 877, 596, 814 };
  const double area[nEtaBins] = { eastArea, midArea, westArea };
  const TString UEetaBinString[nEtaBins] = { "-1.0 < UE #eta < -0.3", "-0.3 < UE #eta < 0.3", "0.3 < UE #eta < 1.0" };

  int jeval, ueeval, pval, eaval;
  TString name, saveName, title, avg, sigma, drawString;
  TString dirName = "plots/LeadJetBBCdist";

  double leadPt, leadPtCorrected, BbcAdcSumEast;
  int nUEpart_chg;
  
  TFile* inFile = new TFile("../out/UE/pAuHTjetUE_leadPtUncorrected.root", "READ");
  TH1D *hLead = (TH1D*)inFile->Get("hLead");
  
  TTree *jt = (TTree*)inFile->Get("HTjetTree");
  jt->SetBranchAddress("leadPt",&leadPt);
  jt->SetBranchAddress("leadPtCorrected",&leadPtCorrected);
  jt->SetBranchAddress("BbcAdcSumEast",&BbcAdcSumEast);
  jt->SetBranchAddress("nUEpart_chg",&nUEpart_chg);

  TH2D *hBBCEvsPt = new TH2D("hBBCEvsPt", ";Leading jet p_{T} (GeV);BBCE ADC sum", 55,4.,59., 12,3000.,70000. );
  jt->Draw("BbcAdcSumEast:leadPt>>hBBCEvsPt","","COLZ");
  TH2D *hBBCEvsEta = new TH2D("hBBCEvsEta", ";Leading jet #eta;BBCE ADC sum", 20,-1.,1., 12,3000.,70000. );
  jt->Draw("BbcAdcSumEast:leadEta>>hBBCEvsEta","","COLZ");
  TH2D *hBBCEvsPtCorrected = new TH2D("hBBCEvsPtCorrected", ";Corrected leading jet p_{T} (GeV);BBCE ADC sum", 55,4.,59., 12,3000.,70000. );
  jt->Draw("BbcAdcSumEast:leadPtCorrected>>hBBCEvsPtCorrected","","COLZ");
  TH2D *nUEpart_chgVsPt = new TH2D("nUEpart_chgVsPt", ";Leading jet p_{T} (GeV);N_{ch}", 55,4.,59., 20,0.,20. );
  jt->Draw("nUEpart_chg:leadPt>>nUEpart_chgVsPt","","COLZ");
  TH2D *nUEpart_chgVsEta = new TH2D("nUEpart_chgVsEta", ";Leading jet #eta;N_{ch}", 20,-1.,1., 20,0.,20. );
  jt->Draw("nUEpart_chg:leadEta>>nUEpart_chgVsEta","","COLZ");

  TH1D *hBBCE_byEta[nEtaBins];
  for (int e=0; e<nEtaBins; ++e) {
    name = "hBBCE" + etaBinName[e];
    hBBCEvsEta->GetXaxis()->SetRangeUser(etaLo[e],etaHi[e]);
    hBBCE_byEta[e] = (TH1D*)hBBCEvsEta->ProjectionY(name);  //,hBBCEvsPt->FindBin(ptLo[p]),hBBCEvsPt->FindBin(ptHi[p]));
    hBBCE_byEta[e]->Scale(1./hBBCE_byEta[e]->Integral());
    hBBCEvsEta->GetXaxis()->SetRange(1,-1);
  }

  
  TH1D *hBBCE_uncorrected[nPtBins];
  for (int p=0; p<nPtBins; ++p) {
    name = "hBBCE_uncorrected" + ptBinName[p];
    hBBCEvsPtCorrected->GetXaxis()->SetRangeUser(ptLo[p],ptHi[p]);
    hBBCE_uncorrected[p] = (TH1D*)hBBCEvsPtCorrected->ProjectionY(name);//,hBBCEvsPt->FindBin(ptLo[p]),hBBCEvsPt->FindBin(ptHi[p]));
    hBBCE_uncorrected[p]->Scale(1./hBBCE_uncorrected[p]->Integral());
    hBBCEvsPtCorrected->GetXaxis()->SetRange(1,-1);
  }
  TCanvas *can = new TCanvas();
  

  TH1D *hBBCE[55];
  for (int pp=0; pp<55; ++pp) {
    int binlo = pp + 1;  int binhi = pp + 2;  double p_lo = 4.0 + (1.0*pp);  double p_hi = 5.0 + (1.0*pp);  // project BBCE dist. for each 1GeV leading jet pT bin at det-level
    // cout<< p_lo << " \t" << hBBCEvsPt->GetXaxis()->GetBinLowEdge(binlo) <<endl;
    name = "hBBCE_"; name += (int)p_lo; name += "_"; name += (int)p_hi;
    hBBCE[pp] = (TH1D*)hBBCEvsPt->ProjectionY(name,binlo,binlo);
    // cout<< hBBCE[pp]->GetEntries() <<endl;
    hBBCE[pp]->SetName(name);
    hBBCE[pp]->Scale(1./hBBCE[pp]->Integral());  // normalize BBCE dist. in each pT bin to unity
    // hBBCE[pp]->Scale(1./hLead->GetBinContent(binlo));  // normalize BBCE dist. in each pT bin to nJets per bin
  }
  
  // TCanvas *can = new TCanvas();
  // can->SetLogz();
  // hBBCEvsPt->Draw("COLZ");
  // can->SaveAs("skrrrrt.pdf","PDF");

  TFile* embFile = new TFile("../embedding/out/sim/pAu2015embedding.root", "READ");
  TH2D *hResponse = (TH2D*)embFile->Get("hResponse");
  TH1D *hPart = (TH1D*)hResponse->ProjectionX("hPart");

  TH1D *FC_part[20], *hBBCE_part[20];
  for (int pp=0; pp<20; ++pp) {
    int plo = pp+7;  int phi = pp+8;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);
    
    FC_part[pp] = GenerateFractionalContribution( hResponse, p_lo, p_hi, dirName );   // project pT response to get det-level jet pT dist corresponding to a part-level pT bin
    name = "FC_part_"; name+=plo; name+="_"; name+=phi; name+="GeV";
    FC_part[pp]->SetName(name);

    hBBCE_part[pp] = WeightAndSumByFC1D( FC_part[pp], hBBCE );        // sum per FC to get BBCE dist per part-level leading jet pT
    name = "hBBCE_part_"; name+=plo; name+="_"; name+=phi; name+="GeV";
    hBBCE_part[pp]->SetName(name);
  }

  // sum BBCE dists at part-level per part-level pT spectrum from embedding
  TH1D *hBBCE_dist[nPtBins];
  for (int p=0; p<nPtBins; ++p) {
    name = "hBBCE_dist" + ptBinName[p];
    hBBCE_dist[p] = new TH1D(name,"BBCE ADC sum", 12,3000.,70000.);
  }
  
  for (int pp=0; pp<20; ++pp) {
    int plo = pp+7;  int phi = pp+8;  double p_lo = 10.0 + (1.0*pp);  double p_hi = 11.0 + (1.0*pp);

    int pval = 999;
    for (int p=0; p<nPtBins; ++p) {
      if ( p_lo>=ptLo[p] && p_hi<=ptHi[p] ) { pval = p; }
    }
    // cout<<pval<<endl;

    double weight = hPart->Integral(plo,plo) / hPart->Integral( hPart->FindBin(ptLo[pval]) , hPart->FindBin(ptHi[pval])-1  );
    hBBCE_dist[pval]->Add( hBBCE_part[pp], weight );
    
  }



  // set draw styles, plot, and make ratios here
  TCanvas *ratcan = new TCanvas( "ratcan" , "" ,700 ,700 );
  ratcan->cd();
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();

  // TH2D *hScale = new TH2D("hScale",";BBCE ADC sum;corrected / uncorrected", 12,3000.,70000., 5,0.,0.2 );
  // hScale->SetStats(0);
  TH2D *hRatioScale = new TH2D("hRatioScale",";BBCE ADC sum;corrected / uncorrected", 12,3000.,70000., 5,0.8,1.6 );
  hRatioScale->SetStats(0);
  TH1D *hBBCE_ratio[nPtBins];

  for (int p=0; p<nPtBins; ++p) {
    hBBCE_dist[p]->SetMarkerStyle(ptCorrectedMarker[p]);
    hBBCE_dist[p]->SetMarkerColor(ptColor[p]);
    hBBCE_dist[p]->SetLineColor(ptColor[p]);
    hBBCE_uncorrected[p]->SetMarkerStyle(ptMarker[p]);
    hBBCE_uncorrected[p]->SetMarkerColor(ptColor[p]);
    hBBCE_uncorrected[p]->SetLineColor(ptColor[p]);

    hBBCE_ratio[p] = (TH1D*)hBBCE_dist[p]->Clone();
    hBBCE_ratio[p]->Divide(hBBCE_uncorrected[p]);
    name = "hBBCE_ratio" + ptBinName[p];
    hBBCE_ratio[p]->SetName(name);
    // hBBCE_ratio[p]->SetMarkerStyle(ptMarker[p]);
  }
  // hBBCE_dist[2]->SetMarkerSize(  );  // the stars are sometimes smaller than the other markers
  // hBBCE_uncorrected[2]->SetMarkerSize(  );



    
  pad1->cd();
  // hScale->Draw();
  hBBCE_uncorrected[0]->SetStats(0);
  hBBCE_uncorrected[0]->GetYaxis()->SetRangeUser(0.,0.2);
  hBBCE_uncorrected[0]->Draw("lpf");
  title = ptBinString[0];
  hBBCE_uncorrected[0]->SetTitle(title);
  title = ptBinStringCorrected[0];
  hBBCE_dist[0]->SetTitle(title);
  hBBCE_dist[0]->Draw("SAMEnostacklpf");
  for (int p=1; p<nPtBins; ++p) {
    title = ptBinString[p];
    hBBCE_uncorrected[p]->GetYaxis()->SetRangeUser(0.,0.2);
    hBBCE_uncorrected[p]->SetTitle(title);
    hBBCE_uncorrected[p]->Draw("SAMEnostacklpf");
    title = ptBinStringCorrected[p];
    hBBCE_dist[p]->SetTitle(title);
    hBBCE_dist[p]->Draw("SAMEnostacklpf");
  }
  pad2->cd();
  hRatioScale->GetYaxis()->SetTitleOffset(0.5);
  hRatioScale->GetYaxis()->SetTitleSize(0.06);
  hRatioScale->Draw();
  for (int p=0; p<nPtBins; ++p) {
    hBBCE_ratio[p]->Draw("SAMEnostacklpf");
  }
  name = dirName + "/BBCEsumDistRatio.pdf";
  // ratcan->BuildLegend();
  // hScale->SetLineColorAlpha(1,0.);
  TLegend *leg = pad1->BuildLegend(0.45,0.65,0.89,0.89);
  leg->SetNColumns(2);
  leg->SetColumnSeparation(0.1);
  leg->SetLineColorAlpha(1,0.);
  hBBCE_uncorrected[0]->SetTitle("");
  ratcan->SaveAs( name, "PDF");




  for (int p=0; p<nPtBins; ++p) {
    hBBCE_dist[p]->SetMarkerStyle(ptMarker[p]);
  }
  
  TCanvas *can0 = new TCanvas( "can0" , "" ,700 ,500 );
  can0->SetLeftMargin(0.12);
  hBBCE_dist[0]->SetStats(0);
  hBBCE_dist[0]->GetYaxis()->SetRangeUser(0.,0.2);
  title = ptBinStringCorrected[0];
  hBBCE_dist[0]->SetTitle(title);
  hBBCE_dist[0]->Draw("lpf");
  for (int p=1; p<nPtBins; ++p) {
    title = ptBinString[p];
    title = ptBinStringCorrected[p];
    hBBCE_dist[p]->SetTitle(title);
    hBBCE_dist[p]->Draw("SAMEnostacklpf");
  }
  TLegend *leg0 = can0->BuildLegend(0.65,0.65,0.89,0.89);
  leg0->SetLineColorAlpha(1,0.);
  hBBCE_dist[0]->SetTitle(";EA_{BBC};#frac{1}{N_{events}} #frac{dN_{events}}{dEA_{BBC}}");

  TString plotString[5] = {"p+Au #sqrt{#it{s}_{NN}} = 200 GeV",
			   "E^{trig}_{T} > 5.4 GeV",
			   "anti-k_{T} R=0.4 jets",
			   "#left|#eta^{jets}#right| < 0.6",
			   "10 < p_{T,lead} < 30 [GeV/#it{c}]"};

  for (int q=0; q<5; ++q) {
    double xpos = 0.17;
    double ypos = 0.5 - 0.05*q;
    drawText( plotString[q], xpos, ypos, 16 );
  }
  
  name = dirName + "/EABBCEDist.pdf";
  can0->SaveAs( name, "PDF");
  // name = dirName + "/EABBCEDistRatio.png";
  // can0->SaveAs( name, "PNG");


  int etaMarker[nEtaBins] = { 29, 33, 34 };
  for (int e=0; e<nEtaBins; ++e) {
    hBBCE_byEta[e]->SetMarkerStyle(etaMarker[e]);
    // hBBCE_byEta[e]->SetMarkerSize(1.5);
    hBBCE_byEta[e]->SetMarkerColor(etaColor[e]);
    hBBCE_byEta[e]->SetLineColor(etaColor[e]);
  }
  
  // TCanvas *can0 = new TCanvas( "can0" , "" ,700 ,500 );
  can0->SetLeftMargin(0.12);
  hBBCE_byEta[0]->SetStats(0);
  hBBCE_byEta[0]->GetYaxis()->SetRangeUser(0.,0.2);
  title = etaBinString[0];
  hBBCE_byEta[0]->SetTitle(title);
  hBBCE_byEta[0]->Draw("lpf");
  for (int e=1; e<nEtaBins; ++e) {
    title = etaBinString[e];
    hBBCE_byEta[e]->SetTitle(title);
    hBBCE_byEta[e]->Draw("SAMEnostacklpf");
  }
  TLegend *leg1 = can0->BuildLegend(0.65,0.65,0.89,0.89);
  leg1->SetLineColorAlpha(1,0.);
  hBBCE_byEta[0]->SetTitle(";EA_{BBC};#frac{1}{N_{events}} #frac{dN_{events}}{dEA_{BBC}}");


  for (int q=0; q<5; ++q) {
    double xpos = 0.17;
    double ypos = 0.5 - 0.05*q;
    drawText( plotString[q], xpos, ypos, 16 );
  }
  
  name = dirName + "/EABBCEDist_byEta.pdf";
  can0->SaveAs( name, "PDF");
  // // name = dirName + "/EABBCEDist_byEta.png";
  // // can0->SaveAs( name, "PNG");


  


  TFile *outFile = new TFile("out/LeadJetBBCdist.root","RECREATE");
  // outFile->cd();
  hResponse->Write();
  hPart->Write();
  hBBCEvsPt->Write();
  nUEpart_chgVsPt->Write();
  hBBCEvsEta->Write();
  for (int pp=0; pp<55; ++pp) {
    hBBCE[pp]->Write();
    if (pp<20) { FC_part[pp]->Write(); }
    if (pp<nPtBins) {
      hBBCE_uncorrected[pp]->Write();
      hBBCE_dist[pp]->Write();
    }
  }

  for (int e=0; e<nEtaBins; ++e) { hBBCE_byEta[e]->Write(); }
  // outFile->Write();
  // // outFile->Close();
  
  // embFile->Close();
  // inFile->Close();
  
}
