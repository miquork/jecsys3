// Purpose: Derive MC truth JES + JER from dijet MC
//          1) Primary fit with low pT bias corrected JES+JER vs gen pT
//            1a) Arithmetic mean/RMS + low pT bias corrections
//            1b) Gaussian mean+sigma + total integral constraint (low pT bias)
//            1c) Median+CI68 + low pT bias corrections
//          2) Invert response vs reco pT for JEC
//            2a) Implement effective physical parameterization
//          3) Additional (small) JER+non-linearity correction to JEC 
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TH3D.h"

#include "tdrstyle_mod22.C"

bool debug = false;
bool correctLowPtBias = true;//false;

void evalJES(TH1D *h);
void evalJER(TH1D *h, TProfile *p);
void evalJET(TH1D *h, const TProfile *p);
void evalJEC(TH1D *h);

void fixJESandJER(TH1D *hjes, TH1D *hjec, TH1D *hjer, TH1D *heff,
		  const TF1 *f1jer);

void fixMedian(TH1D *hjes, TF1* jes_fit, TF1* jer_fit, double pT_threshold);

double evalXmin(double jes, double jer, double eff);
TH1D* cutHist(const TH1D *h, const double xmin);

//void MCTruth(string set="Winter25MG") {
//void MCTruth(string set="Summer24MG") {
void MCTruth(string set="Summer24MC_Flat") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/MCTruth");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/MCTruth");
  
  #include "Config.C"
  
  TFile *f(0);
  const char *cs = set.c_str();
  if (mfile.find(Form("JERC_%s_MC",cs))!=mfile.end()) {
    string file = mfile[Form("JERC_%s_MC",cs)];
    cout << "Reading JERC_" << set << "_MC from Config.C: " << file << endl;
    f = new TFile(file.c_str(),"READ");
  }
  else {
  if (set=="Summer24JME") f = new TFile("rootfiles/Prompt2024/v124_Jet/jmenano_mc_out_Summer24MC_FlatJMEN_v124.root","READ");
  if (set=="Summer24") f = new TFile("rootfiles/Prompt2024/v121_Jet/jmenano_mc_out_Summer24MG_v121_2024CDEFGHI.root","READ");
  // v4->v6: remove _Athens, count 3 leading gen (was 4), Jet_neMultiplicity for reco (was gen)
  // v6->v8: add p2eff_recEta hybrid map
  if (set=="Winter25") f = new TFile("rootfiles/Prompt2025/v123_v8_Jet/jmenano_mc_out_Winter25MC_Flat2022_v123_v8.root","READ");
  if (set=="Winter25MG") f = new TFile("rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root","READ");
  if (set=="Winter24") f = new TFile("rootfiles/Prompt2025/v123_v8_Jet/jmenano_mc_out_Winter24MCFlat_v123_v8.root","READ");
  }
  assert(f && !f->IsZombie());

  gDirectory->cd("HLT_MC");
  gDirectory->cd("MCtruth");
  //gDirectory->cd("LeadingJ");
  TDirectory *d = gDirectory;

  // <pTreco/pTgen> vs (eta_jet, pTgen); 41 bins in X=eta_jet
  TProfile2D *p2x = (TProfile2D*)d->Get("p2jes");   assert(p2x);
  //TProfile2D *p2e = (TProfile2D*)d->Get("p2eff");   assert(p2e);
  //TProfile2D *p2e = (TProfile2D*)d->Get("p2eff_recEta");   assert(p2e);
  TProfile2D *p2e = (TProfile2D*)d->Get("p2eff_recEta");
  if (!p2e) p2e = (TProfile2D*)d->Get("p2eff");   assert(p2e);
  TProfile2D *p2r = (TProfile2D*)d->Get("p2r_raw"); assert(p2r);
  TProfile2D *p2c = (TProfile2D*)d->Get("p2r");     assert(p2c);

  //if (set=="Winter25") {
  //p2r = (TProfile2D*)d->Get("LeadingJ/p2r_Athens"); assert(p2r);
  //}
    
  // 7*6=42>41 pads; last one for legend
  int size = 600;
  TCanvas *c1eff = new TCanvas("c1eff","c1eff",7*size,6*size);
  c1eff->Divide(7,6,0,0);
  TCanvas *c1jes = new TCanvas("c1jes","c1jes",7*size,6*size);
  c1jes->Divide(7,6,0,0);
  TCanvas *c1jesx = new TCanvas("c1jesx","c1jesx",7*size,6*size);
  c1jesx->Divide(7,6,0,0);
  TCanvas *c1jer = new TCanvas("c1jer","c1jer",7*size,6*size);
  c1jer->Divide(7,6,0,0);
  TCanvas *c1jerx = new TCanvas("c1jerx","c1jerx",7*size,6*size);
  c1jerx->Divide(7,6,0,0);
  TCanvas *c1jet = new TCanvas("c1jet","c1jet",7*size,6*size);
  c1jet->Divide(7,6,0,0);
  TCanvas *c1jec = new TCanvas("c1jec","c1jec",7*size,6*size);
  c1jec->Divide(7,6,0,0);
  TCanvas *c1jecx = new TCanvas("c1jecx","c1jecx",7*size,6*size);
  c1jecx->Divide(7,6,0,0);


  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);

  c1eff->cd(42);
  int nentry = (correctLowPtBias ? 6 : 3);
  TLegend *legeff = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legeff->SetTextSize(0.045*1.5);
  
  c1jes->cd(42);
  TLegend *legjes = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legjes->SetTextSize(0.045*1.5);

  c1jesx->cd(42);
  TLegend *legjesx = tdrLeg(0.05,0.90-0.05*1.5*9,0.50,0.90);
  legjesx->SetTextSize(0.045*1.5);
  
  c1jer->cd(42);
  TLegend *legjer = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legjer->SetTextSize(0.045*1.5);

  c1jerx->cd(42);
  TLegend *legjerx = tdrLeg(0.05,0.80-0.05*1.5*9,0.50,0.80);
  legjerx->SetTextSize(0.045*1.5);

  c1jec->cd(42);
  TLegend *legjec = tdrLeg(0.05,0.80-0.05*1.5*3,0.50,0.80);
  legjec->SetTextSize(0.045*1.5);

  c1jecx->cd(42);
  TLegend *legjecx = tdrLeg(0.05,0.90-0.05*1.5*9,0.50,0.90);
  legjecx->SetTextSize(0.045*1.5);


  double eps = 1e-4;

  vector<TF1*> vjer(p2r->GetNbinsX());
  vector<TF1*> vjes(p2r->GetNbinsX());
  vector<TH1D*> veff(p2r->GetNbinsX());
  for (int i = 1; i != p2r->GetNbinsX()+1; ++i) {
    //for (int j = 1; j != p2r->GetNbinsY()+1; ++j) {

    //////////////////////////////////////////////////////////////////////
    // Step 1. First round fitting low pT biased results at high pT end //
    //////////////////////////////////////////////////////////////////////

    // Prior MC truth JES as applied on MC
    TProfile *px = p2x->ProfileY(Form("px_%d",i),i,i);
    TH1D *hx = px->ProjectionX(Form("hx_%d",i));
    
    // Initial EFF, JES, JER and JET counts
    TProfile *pe = p2e->ProfileY(Form("pe_%d",i),i,i);
    TH1D *he = pe->ProjectionX(Form("he_%d",i));
    TProfile *pr = p2r->ProfileY(Form("pr_%d",i),i,i);
    TH1D *hr = pr->ProjectionX(Form("hr_%d",i));
    evalJES(hr);
    TProfile *pc = p2c->ProfileY(Form("pc_%d",i),i,i);
    TH1D *hc = pc->ProjectionX(Form("hc_%d",i));
    evalJEC(hc);
    pc->SetErrorOption("S");
    TH1D *hs = pc->ProjectionX(Form("hs_%d",i));
    evalJER(hs, pc);
    TH1D *hn = pc->ProjectionX(Form("hn_%d",i));
    evalJET(hn, pc);

    if (set=="Winter25") hn->Scale(1000.); // Flat2022QCD
    if (set=="Winter24") hn->Scale(1000.); // Flat2022QCD
    
    ////////////////////////////////////////////
    // Step 1.1 Jet reconstruction EFFiciency //
    ////////////////////////////////////////////
    
    c1eff->cd(i);
    gPad->SetLogx();

    double effmin(0.), effmax(1.2);
    TH1D *heff = tdrHist(Form("heff_%d",i),"Efficiency",effmin+eps,effmax-eps,
		      "p_{T,gen} (GeV)",5,5000);
    heff->GetXaxis()->SetMoreLogLabels(kFALSE);
    heff->Draw("AXIS");
    l->DrawLine(5,1,5000,1);

    double etamin = p2r->GetXaxis()->GetBinLowEdge(i);
    double etamax = p2r->GetXaxis()->GetBinLowEdge(i+1);
    double ptmin_eff = 15.;
    double ptmin = 30.;
    double ptmin_fwd = 30.;
    double ptmax = min(4000.,0.7*6800./cosh(etamax));
    TF1 *f1eff = new TF1(Form("f1eff_%d",i),
			 "0.5*(1+erf((x-[0])/([1]*[0])))",
			 5.,ptmax);
    f1eff->SetParameters(15.,0.30);
    he->Fit(f1eff,"QRN");
    
    tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(he,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1eff->Draw("SAME");

    if (i==1) {
      if (set=="Summer24JME") legeff->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legeff->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legeff->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legeff->SetHeader("Winter24, QCDFlat MC");
      else legeff->SetHeader(set.c_str());
      legeff->AddEntry(hr,"Data","PLE");
      legeff->AddEntry(f1eff,"Fit","L");
    }

    
    /////////////////////////////////////
    // Step 1.2 Jet Energy Scale (JES) //
    /////////////////////////////////////

    c1jes->cd(i);
    gPad->SetLogx();

    double jesmin(0.5), jesmax(1.5);
    // Range after low pT bias correction
    if (i<=7)       { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=14) { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=21) { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=28) { jesmin = 0.45; jesmax = 1.10; }
    else if (i<=35) { jesmin = 0.25; jesmax = 1.25; }
    else if (i<=42) { jesmin = 0.25; jesmax = 1.20; }
    TH1D *hjes = tdrHist(Form("hjes_%d",i),"JES",jesmin+eps,jesmax-eps,
		      "p_{T,gen} (GeV)",5,5000);
    hjes->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjes->Draw("AXIS");
    l->DrawLine(5,1,5000,1);
    
    TF1 *f1jes = new TF1(Form("f1jes_%d",i),
			 "[0]+[1]*pow(x,[2])+[3]/x+[4]*1e-5*x",
			 ptmin,ptmax);
    // Enforce physical ranges
    f1jes->SetParameters(1,-1,-0.3,0,1);
    f1jes->SetParLimits(0,0.8,1.2);
    f1jes->SetParLimits(1,-2,0.);
    f1jes->SetParLimits(2,-0.4,-0.2);
    f1jes->SetParLimits(3,-10.,+10.);
    f1jes->SetParLimits(4,-10,+10.);
    // Switch off trk parameters outside tracker
    if (etamin>2.65) {
      //f1jes->FixParameter(4, 0.);
      f1jes->SetParLimits(4, 0.,+10.);
      f1jes->SetRange(ptmin_fwd,ptmax);
    }
    hr->Fit(f1jes,"QRN");

    tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jes->Draw("SAME");

    if (i==1) {
      if (set=="Summer24JME") legjes->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjes->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjes->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjes->SetHeader("Winter24, QCDFlat MC");
      else legjes->SetHeader(set.c_str());
      legjes->AddEntry(hr,"Data","PLE");
      legjes->AddEntry(f1jes,"Fit","L");
    }

    
    c1jesx->cd(i);
    gPad->SetLogx();
    hjes->Draw("AXIS");
    l->DrawLine(5,1,5000,1);
    
    
    //////////////////////////////////////////
    // Step 1.3 Jet Energy Resolution (JER) //
    //////////////////////////////////////////
    
    c1jer->cd(i);
    gPad->SetLogx();

    double jermin(0.), jermax(0.5);
    // Range after low pT bias correction
    if (i<=7)       { jermin = 0.; jermax = 0.50; }
    else if (i<=14) { jermin = 0.; jermax = 0.50; }
    else if (i<=21) { jermin = 0.; jermax = 0.55; }
    else if (i<=28) { jermin = 0.; jermax = 1.00; }
    else if (i<=35) { jermin = 0.; jermax = 1.10; }
    else if (i<=42) { jermin = 0.; jermax = 0.75; }
    TH1D *hjer = tdrHist(Form("hjer_%d",i),"JER",jermin+eps,jermax-eps,
			 "p_{T,gen} (GeV)",5,5000);
    hjer->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjer->Draw("AXIS");

    TF1 *f1jer = new TF1(Form("f1jer_%d",i),
			 "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
			 ptmin,ptmax);
    f1jer->SetParameters(1,1,0.05);
    f1jer->SetParLimits(0,0,10.);
    f1jer->SetParLimits(1,0.5,1.5);
    f1jer->SetParLimits(2,0.03,0.12);
    if (etamin>3.139) {
      f1jer->FixParameter(2,0.08);
    }
    hs->Fit(f1jer,"QRN");

    tex->DrawLatex(i%7==1 ? 0.50 : 0.40, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jer->Draw("SAME");

    if (i==1) {
      if (set=="Summer24JME") legjer->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjer->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjer->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjer->SetHeader("Winter24, QCDFlat MC");
      else legjer->SetHeader(set.c_str());
      legjer->AddEntry(hr,"Data","PLE");
      legjer->AddEntry(f1jes,"Fit","L");
    }


    c1jerx->cd(i);
    gPad->SetLogx();
    hjer->Draw("AXIS");

    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jer->Draw("SAME");



    /////////////////////////
    // Step 1.4 JET counts //
    /////////////////////////
    
    c1jet->cd(i);
    gPad->SetLogx();
    gPad->SetLogy();
    
    TH1D *hjet = tdrHist(Form("hjet_%d",i),
			 "Weighted N(jet)/d#eta/dp_{T} (1/GeV)",
			 1e-10,1e10,"p_{T,gen} (GeV)",5,5000);
    hjet->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjet->GetYaxis()->SetMoreLogLabels(kFALSE);
    hjet->Draw("AXIS");

    TF1 *f1jet = new TF1("f1jet","[0]*pow(x,[1])*pow(1-2*x/[3],[2])",
			 ptmin,ptmax);
    f1jet->SetParameters(1e11,-5,5,13800./cosh(etamin));
    f1jet->FixParameter(3,13800./cosh(etamin));
    hn->Fit(f1jet,"QRN");

    tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hn,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jet->Draw("SAME");

    
    ///////////////////////////////////////////////////
    // Step 1.4 Jet Enercy Corrections (JEC) closure //
    ///////////////////////////////////////////////////
    
    c1jec->cd(i);
    gPad->SetLogx();

    double jecmin(-10), jecmax(10);
    /*
    // Range before low pT bias
    if (i<=7)       { jecmin = -5; jecmax = +10; }
    else if (i<=14) { jecmin = -5; jecmax = +10; }
    else if (i<=21) { jecmin = -5; jecmax = +10; }
    else if (i<=28) { jecmin = -8; jecmax = +16; }
    else if (i<=35) { jecmin = -10; jecmax = +16; }
    else if (i<=42) { jecmin = -12; jecmax = +16; }
    */
    // Range after low pT bias
    /*
    if (i<=7)       { jecmin = -5; jecmax = +10; }
    else if (i<=14) { jecmin = -5; jecmax = +10; }
    else if (i<=21) { jecmin = -5; jecmax = +10; }
    else if (i<=28) { jecmin = -8; jecmax = +16; }
    else if (i<=35) { jecmin = -10; jecmax = +16; }
    else if (i<=42) { jecmin = -12; jecmax = +16; }
    */
    TH1D *hjec = tdrHist(Form("hjec_%d",i),"JEC closure (%)", jecmin, jecmax,
			 "p_{T,gen} (GeV)",5,5000);
    hjec->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjec->Draw("AXIS");

    l->SetLineStyle(kDotted);
    l->DrawLine(5,-1,5000,-1);
    l->DrawLine(5,+1,5000,+1);
    l->SetLineStyle(kDashed);
    l->DrawLine(5,0,5000,0);

    tex->DrawLatex(i%7==1 ? 0.50 : 0.40, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hc,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

    
    if (i==1) {
      if (set=="Summer24JME") legjec->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjec->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjec->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjec->SetHeader("Winter24, QCDFlat MC");
      else legjec->SetHeader(set.c_str());
      //legjec->AddEntry(hc,"Data","PLE");
    }

    // Save fits and efficiency for Step 3.
    vjer[i-1] = f1jer;
    vjes[i-1] = f1jes;
    veff[i-1] = he;

    
    /////////////////////////////////////////////////////////////////////////
    // Step 2. Second round with low pT bias correction and extending down //
    /////////////////////////////////////////////////////////////////////////
    if (!correctLowPtBias) continue;

    TF1 *f1eff_raw = (TF1*)f1eff->Clone(Form("f1eff_raw_%d",i));
    TH1D *he_raw = (TH1D*)he->Clone(Form("he_raw_%d",i));
    TF1 *f1jes_raw = (TF1*)f1jes->Clone(Form("f1jes_raw_%d",i));
    TH1D *hr_raw = (TH1D*)hr->Clone(Form("hr_raw_%d",i));
    TF1 *f1jer_raw = (TF1*)f1jer->Clone(Form("f1jer_raw_%d",i));
    TH1D *hs_raw = (TH1D*)hs->Clone(Form("hs_raw_%d",i));

    TH1D *hc_raw = hc; hc_raw->SetName(Form("hc_raw_%d",i));
    pc->SetErrorOption("");
    hc = pc->ProjectionX(Form("hc_%d",i));

    fixJESandJER(hr, hc, hs, he, f1jer);
    
    evalJEC(hc);

    c1eff->cd(i);

    f1eff->SetRange(15.,ptmax);
    he->Fit(f1eff,"QRN");

    tdrDraw(he_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1eff_raw->Draw("SAME");
    tdrDraw(he,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1eff->SetLineColor(kGreen+2);
    f1eff->Draw("SAME");

    if (i==1) {
      legeff->AddEntry(he_raw,"Orig. data","PLE");
      legeff->AddEntry(f1eff_raw,"Orig. fit","L");
      legeff->AddEntry(f1eff_raw," "," ");
    }

    
    c1jer->cd(i);

    f1jer->SetRange(15.,ptmax);
    hs->Fit(f1jer,"QRN");

    tdrDraw(hs_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1jer_raw->Draw("SAME");
    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1jer->SetLineColor(kGreen+2);
    f1jer->Draw("SAME");

    if (i==1) {
      legjer->AddEntry(hs_raw,"Orig. data","PLE");
      legjer->AddEntry(f1jer_raw,"Orig. fit","L");
      legjer->AddEntry(f1jer_raw," ","");
    }


    c1jerx->cd(i);
	
    tdrDraw(hs_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1jer_raw->Draw("SAME");
    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1jer->Draw("SAME");

    if (i==1) {
      if (set=="Summer24JME") legjerx->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjerx->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjerx->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjerx->SetHeader("Winter24, QCDFlat MC");
      else legjerx->SetHeader(set.c_str());
      legjerx->AddEntry(hs_raw,"Raw RMS","PLE");
      legjerx->AddEntry(f1jer_raw,"Raw fit","L");
      legjerx->AddEntry(hs,"Corrected RMS","PLE");
      legjerx->AddEntry(f1jer_raw,"Corrected fit","L");
    }

        
    c1jes->cd(i);

    f1jes->SetRange(15.,ptmax);
    hr->Fit(f1jes,"QRN");

    //TH1D *hx_raw = (TH1D*)hx->Clone(Form("hx_raw_%d",i));
    //fixMedian(hx, f1jes, f1jer, 0.);
    
    //tdrDraw(hx_raw,"HIST][",kNone,kCyan+1,kSolid,-1,kNone,0);
    tdrDraw(hx,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
    
    tdrDraw(hr_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1jes_raw->Draw("SAME");
    tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1jes->SetLineColor(kGreen+2);
    f1jes->Draw("SAME");

    if (i==1) {
      legjes->AddEntry(hr_raw,"Orig. data","PLE");
      legjes->AddEntry(f1jes_raw,"Orig. fit","L");
      legjes->AddEntry(hx,"Previous JEC","F");
    }

    
    c1jesx->cd(i);
    gPad->SetLogx();

    tex->DrawLatex(i%7==1 ? 0.35 : 0.25, 0.90,
		   Form("%1.3f<|#eta_{rec}|<%1.3f",etamin,etamax));
    tdrDraw(hx,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
    f1jes_raw->Draw("SAME");
    f1jes->Draw("SAME");

    tdrDraw(hr_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);

    if (i==1) {
      if (set=="Summer24JME") legjesx->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjesx->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjesx->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjesx->SetHeader("Winter24, QCDFlat MC");
      else legjesx->SetHeader(set.c_str());
      legjesx->AddEntry(hx,"Official JEC","F");
      legjesx->AddEntry(hr_raw,"Raw Mean","PLE");
      legjesx->AddEntry(f1jes_raw,"Raw Mean fit (p_{T}>30 GeV)","L");
      legjesx->AddEntry(hr,"Corrected Mean","PLE");
      legjesx->AddEntry(f1jes,"Corrected Mean fit (p_{T}>15 GeV)","L");
    }


    c1jec->cd(i);

    TH1D *hr_ratio = (TH1D*)hr->Clone(Form("hr_ratio_%d",i));
    //hr_ratio->Divide(hx); evalJEC(hr_ratio);
    hr_ratio->Divide(f1jes); evalJEC(hr_ratio);
    
    tdrDraw(hc_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    //tdrDraw(hr_ratio,"Pz",kOpenCircle,kBlue,kSolid,-1,kNone,0,2.0,2);
    tdrDraw(hc,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    tdrDraw(hr_ratio,"Pz",kFullDotSmall,kGreen+2,kSolid,-1,kNone,0,2.0,1);

    if (i==1) {
      legjec->AddEntry(hc,"Data","PLE");
      legjec->AddEntry(hc_raw,"Orig. data","PLE");
    }


    c1jecx->cd(i);
    gPad->SetLogx();
    hjec->Draw("AXIS");
    l->DrawLine(5,0,5000,0);
    
    tex->DrawLatex(i%7==1 ? 0.35 : 0.25, 0.90,
		   Form("%1.3f<|#eta_{rec}|<%1.3f",etamin,etamax));


    tdrDraw(hc_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    tdrDraw(hc,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    //tdrDraw(hr_ratio,"Pz",kFullDiamond,kGreen+2,kSolid,-1,kNone,0,2.0,1);

    if (i==1) {
      if (set=="Summer24JME") legjecx->SetHeader("Summer24, JME QCDFlat MC");
      else if (set=="Summer24") legjecx->SetHeader("Summer24, MG QCD MC");
      else if (set=="Winter25") legjecx->SetHeader("Winter25, QCDFlat MC");
      else if (set=="Winter24") legjecx->SetHeader("Winter24, QCDFlat MC");
      else legjecx->SetHeader(set.c_str());
      legjecx->AddEntry(hc_raw,"Raw Mean","PLE");
      legjecx->AddEntry(hc,"Corrected Mean","PLE");
      //legjecx->AddEntry(hr_ratio,"Corr. Mean  / Fit","PLE");
    }
    
    //} // for j

  } // for i


  //const char *cs = set.c_str();
  c1eff->SaveAs(Form("pdf/MCTruth/MCTruth_EFF_%s.pdf",cs));
  c1jes->SaveAs(Form("pdf/MCTruth/MCTruth_JES_%s.pdf",cs));
  c1jer->SaveAs(Form("pdf/MCTruth/MCTruth_JER_%s.pdf",cs));
  c1jet->SaveAs(Form("pdf/MCTruth/MCTruth_JET_%s.pdf",cs));
  c1jec->SaveAs(Form("pdf/MCTruth/MCTruth_JEC_%s.pdf",cs));


  ///////////////////////////////////////////////////////////
  // Step 3. Add median/CI68 and Gaus mean/sigma from TH3D //
  ///////////////////////////////////////////////////////////

  //TF1 *fgaus = new TF1("fgaus","TMath::Gaus(x,[0],[1],0)",0,2);
  TF1 *fgaus = new TF1("fgaus","TMath::Gaus(x,[0],[1],1)",0,2);
  TH3D *h3r = (TH3D*)d->Get("Response3D_raw"); assert(h3r);
  TH3D *h3c = (TH3D*)d->Get("Response3D"); assert(h3c);

  // Make plots of Gaussian response in various pT intervals
  vector<pair<double,double> > ptbins;
  ptbins.push_back(make_pair<double,double>(15,20));
  ptbins.push_back(make_pair<double,double>(30,35));
  ptbins.push_back(make_pair<double,double>(86,110));

  
  TCanvas *c1gaus = new TCanvas("c1gaus","c1gaus",7*size,6*size);
  c1gaus->Divide(7,6,0,0);

  c1gaus->cd(42);
  TLegend *leggaus = tdrLeg(0.05,0.80-0.05*1.5*4,0.50,0.80);
  leggaus->SetTextSize(0.045*1.5);

  // x=eta_gen, y=pTgen, z=pTreco/pTgen
  for (int i = 1; i != h3r->GetNbinsX()+1; ++i) {

    double etamin = h3r->GetXaxis()->GetBinLowEdge(i);
    double etamax = h3r->GetXaxis()->GetBinLowEdge(i+1);
    
    TH1D *h1r = h3r->ProjectionY(Form("h1r_%d",i)); h1r->Reset();
    TH1D *h1c = h3c->ProjectionY(Form("h1c_%d",i)); h1c->Reset();
    TH1D *h1s = h3c->ProjectionY(Form("h1s_%d",i)); h1s->Reset();
    
    TH1D *h1r_median = h3r->ProjectionY(Form("h1r_median_%d",i));
    h1r_median->Reset();
    TH1D *h1c_median = h3c->ProjectionY(Form("h1c_median_%d",i));
    h1c_median->Reset();
    TH1D *h1s_CI68 = h3c->ProjectionY(Form("h1s_CI68_%d",i));
    h1s_CI68->Reset();

    TH1D *h1r_mu = h3r->ProjectionY(Form("h1r_mu_%d",i));
    h1r_mu->Reset();
    TH1D *h1c_mu = h3c->ProjectionY(Form("h1c_mu_%d",i));
    h1c_mu->Reset();
    TH1D *h1s_sigma = h3c->ProjectionY(Form("h1s_sigma_%d",i));
    h1s_sigma->Reset();
    
    for (int j = 1; j != h1r->GetNbinsX()+1; ++j) {
      
      TH1D *hr = h3r->ProjectionZ(Form("hzr_%d_%d",i,j),i,i,j,j);
      TH1D *hc = h3c->ProjectionZ(Form("hzc_%d_%d",i,j),i,i,j,j);

      // Arithmetic mean and RMS to compare to TProfile pr, pc
      double mean = hr->GetMean();
      double emean = hr->GetMeanError();
      double meanc = hc->GetMean();
      double emeanc = hc->GetMeanError();
      double rms = hc->GetRMS();
      double erms = hc->GetRMSError();
      h1r->SetBinContent(j, mean);
      h1r->SetBinError(j, emean);
      h1c->SetBinContent(j, meanc);
      h1c->SetBinError(j, emeanc);
      // This is weird, rms and mean from hr are closer to TProfile pc
      h1s->SetBinContent(j, meanc!=0 ? rms/meanc : rms);
      h1s->SetBinError(j, meanc!=0 ? erms/meanc : erms);

      // Median and CI68
      // probabilities where to evaluate the quantiles in [0,1]
      const int nq = 3;
      Double_t p[nq] = {0.16, 0.5, 0.84};
      Double_t xpr[nq], xpc[nq];
      if (hr->GetEntries()>5) {
	hr->GetQuantiles(nq,xpr,p);
	hc->GetQuantiles(nq,xpc,p);
	
	if (xpr[1]!=0) {
	  h1r_median->SetBinContent(j, xpr[1]);
	  h1r_median->SetBinError(j, xpr[1] * (mean!=0 ? emean / mean : 1));
	  h1c_median->SetBinContent(j, xpc[1]);
	  h1c_median->SetBinError(j, xpc[1] * (meanc!=0 ? emeanc / meanc : 1));
	  h1s_CI68->SetBinContent(j, 0.5*(xpc[2]-xpc[0]) / xpr[1]);
	  h1s_CI68->SetBinError(j, 0.5*(xpc[2]-xpc[0]) / xpr[1] *
				(meanc!=0 ? erms / meanc : 1));
	}
      }

      // Gaussian fits
      double pt = h3r->GetYaxis()->GetBinCenter(j);
      double ptmin = h3r->GetYaxis()->GetBinLowEdge(j);
      double ptmax = h3r->GetYaxis()->GetBinLowEdge(j+1);
      //double pt = h1r->GetBinCenter(j);
      //double ptmin = h1r->GetBinLowEdge(j);
      //double ptmax = h1r->GetBinLowEdge(j+1);
      if (hr->GetEntries()>5 && hr->Integral()>0 && hc->Integral()>0 && pt>20) {
	TF1 *f1jer = vjer[i-1];
	TF1 *f1jes = vjes[i-1];
	TH1D *he = veff[i-1];
	double eff = he->GetBinContent(j);
	hr->Scale(100./2.*eff/hr->Integral());
	hc->Scale(100./2.*eff/hc->Integral());

	double k = 2.5;//2.0;//1.5;
	double mu = f1jes->Eval(pt);
	double sigma = f1jer->Eval(pt) * mu;
	double xmin = evalXmin(mu, sigma/mu, eff);
	//double xmin = max(20.*mu/ptmin, mu-k*sigma);
	//fgaus->SetRange(max(0.1,mu-k*sigma), mu+k*sigma);
	fgaus->SetRange(max(0.1,xmin), mu+k*sigma);
	fgaus->SetParameters(max(0.1,min(1.9,mu)), max(0.01,min(1.0,sigma)));
	fgaus->SetParLimits(0, 0.1, 1.9);
	fgaus->SetParLimits(1, 0.01, 1.0);
	if (pt<100.)
	  fgaus->FixParameter(1, fgaus->GetParameter(1)); // fixed width
	else
	  fgaus->ReleaseParameter(1);
	hr->Fit(fgaus,"QRN"); // not fixed
	h1r_mu->SetBinContent(j, fgaus->GetParameter(0));
	h1r_mu->SetBinError(j, fgaus->GetParError(0));

	double ptmin_gaus(86), ptmax_gaus(110);
	if (pt>ptmin_gaus && pt<ptmax_gaus) {
	//if (pt>86 && pt<110) {
	  //if (pt>21 && pt<28) {
	  //if (pt>=10 && pt<15) {
	  c1gaus->cd(i);

	  if (set=="Summer24JME") { hr->Rebin(4); hr->Scale(0.25); }
	  hr->GetYaxis()->SetRangeUser(0,2.9);//3.);
	  if (ptmin_gaus>30) hr->GetYaxis()->SetRangeUser(0,4.9);
	  hr->UseCurrentStyle();
	  hr->SetYTitle("dN/dR");
	  
	  tdrDraw(hr,"HIST",kNone,kBlue,kSolid,-1,1001,kBlue-9);
	  hr->SetFillColorAlpha(kBlue-9, 0.3);
	  TH1D *hr_cut = cutHist(hr, xmin);
	  tdrDraw(hr_cut,"HIST",kNone,kBlue+1,kSolid,-1,1001,kBlue-9);
	  fgaus->SetLineColor(kCyan+2);
	  TF1 *f1r = (TF1*)fgaus->DrawClone("SAME");

	  tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
			 Form("%1.3f<|#eta_{rec}|<%1.3f",etamin,etamax));
	  if (i==1) {
	    leggaus->AddEntry(f1r,"Raw parameterization","L");
	    leggaus->AddEntry(hr_cut,"Raw response","F");

	    c1gaus->cd(42);
	    tex->DrawLatex(0.20, 0.90,
			   Form("%1.0f<p_{T,gen}<%1.0f GeV",
				ptmin_gaus,ptmax_gaus));
			   //Form("%1.0f<p_{T,gen}<%1.0f GeV",ptmin,ptmax));
	  }
	}

	mu = 1;
	sigma = f1jer->Eval(pt) * mu;
	double xminc = evalXmin(mu, sigma/mu, eff);
	//xmin = max(20.*mu/pt, mu-k*sigma);
	//fgaus->SetRange(max(0.1,mu-k*sigma), mu+k*sigma);
	fgaus->SetRange(max(0.1,xminc), mu+k*sigma);
	fgaus->SetParameters(1, max(0.01,min(1.0,sigma)));
	fgaus->SetParLimits(0, 0.1, 1.9);
	fgaus->SetParLimits(1, 0.01, 1.0);
	if (pt<100.)
	  fgaus->FixParameter(1, fgaus->GetParameter(1)); // fixed width
	else
	  fgaus->ReleaseParameter(1);
	hc->Fit(fgaus,"QRN"); // not fixed
	h1c_mu->SetBinContent(j, fgaus->GetParameter(0));
	h1c_mu->SetBinError(j, fgaus->GetParError(0));
	h1s_sigma->SetBinContent(j, fgaus->GetParameter(1) /
				 fgaus->GetParameter(0));
	h1s_sigma->SetBinError(j, fgaus->GetParError(1) /
			       fgaus->GetParameter(0));

	if (pt>ptmin_gaus && pt<ptmax_gaus) {
	  //if (pt>86 && pt<110) {
	//if (pt>21 && pt<28) {
	  //if (pt>15 && pt<20 && false) {
	  c1gaus->cd(i);
	  if (set=="Summer24JME") { hc->Rebin(4); hc->Scale(0.25); }
	  //if (set=="Winter25MG") hc->Scale(1.5);
	  tdrDraw(hc,"HIST",kNone,kMagenta+1,kSolid,-1,1001,kMagenta-9);
	  hc->SetFillColorAlpha(kMagenta-9,0.3);
	  TH1D *hc_cut = cutHist(hc, xminc);
	  tdrDraw(hc_cut,"HIST",kNone,kMagenta+2,kSolid,-1,1001,kMagenta-9);
	  hc_cut->SetFillColorAlpha(kMagenta-9,0.3);
	  fgaus->SetLineColor(kMagenta+2);
	  TF1 *f1c = (TF1*)fgaus->DrawClone("SAME");

	  if (i==1) {
	    leggaus->AddEntry(f1c,"Corrected parameterization","L");
	    leggaus->AddEntry(hc_cut,"Corrected response","F");
	  }
	  gPad->RedrawAxis();
	}

      }
      
    } // for j

    if (correctLowPtBias) {
      TF1 *f1jer = vjer[i-1];
      TH1D *he = veff[i-1];
      fixJESandJER(h1r, h1c, h1s, he, f1jer);
    }

    
    c1jesx->cd(i);

    tdrDraw(h1r_median,"Pz",kOpenDiamond,kMagenta+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1r,"Pz",kFullDiamond,kOrange+1,kSolid,-1,kNone,0,2.0);
    //tdrDraw(h1r_mu,"Pz",kOpenStar,kCyan+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1r_mu,"Pz",kFullDotSmall,kCyan+1,kSolid,-1,kNone,0,2.0);

    if (i==1) {
      legjesx->AddEntry(h1r_median,"Raw Median","PLE");
      legjesx->AddEntry(h1r,"Corrected Median","PLE");
      //legjesx->AddEntry(h1r_mu,"Gauss Mean (fixed)","PLE");
      legjesx->AddEntry(h1r_mu,"Gauss Mean","PLE");
    }

    
    c1jecx->cd(i);
	
    evalJEC(h1c);
    evalJEC(h1c_median);
    evalJEC(h1c_mu);
    tdrDraw(h1c,"Pz",kFullDiamond,kOrange+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1c_median,"Pz",kOpenDiamond,kMagenta+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1c_mu,"Pz",kFullDotSmall,kCyan+1,kSolid,-1,kNone,0,2.0);

    if (i==1) {
      legjecx->AddEntry(h1c_median,"Raw Median","PLE");
      legjecx->AddEntry(h1c,"Corrected Median","PLE");
      //legjecx->AddEntry(h1c_mu,"Gauss Mean (fixed)","PLE");
      legjecx->AddEntry(h1c_mu,"Gauss Mean","PLE");
    }

    
    c1jerx->cd(i);

    tdrDraw(h1s,"Pz",kFullDiamond,kOrange+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1s_CI68,"Pz",kOpenDiamond,kMagenta+1,kSolid,-1,kNone,0,2.0);
    tdrDraw(h1s_sigma,"Pz",kFullDotSmall,kCyan+1,kSolid,-1,kNone,0,2.0);

    if (i==1) {
      legjerx->AddEntry(h1s_CI68,"Raw CI68","PLE");
      legjerx->AddEntry(h1s,"Corrected CI68","PLE");
      //legjerx->AddEntry(h1s_sigma,"Gauss Sigma (fixed)","PLE");
      legjerx->AddEntry(h1s_sigma,"Gauss Sigma","PLE");
    }
  } // for i

  c1jesx->SaveAs(Form("pdf/MCTruth/MCTruth_JES_xtra_%s.pdf",cs));
  c1jecx->SaveAs(Form("pdf/MCTruth/MCTruth_JEC_xtra_%s.pdf",cs));
  c1jerx->SaveAs(Form("pdf/MCTruth/MCTruth_JER_xtra_%s.pdf",cs));
  c1gaus->SaveAs(Form("pdf/MCTruth/MCTruth_Gaus_%s.pdf",cs));
  
}


// Evaluate JES
// - add minimum uncertainty to avoid overfitting high statistics
void evalJES(TH1D *h) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double jes = h->GetBinContent(i);
    if (jes!=0) {
      double errmin = 0.0020;
      double ejes = sqrt(pow(h->GetBinError(i),2)+errmin*errmin);
      h->SetBinContent(i, jes);
      h->SetBinError(i, ejes);
    }
  }
  return;
}

// Evaluate JER
// - assume input histogram projected from JES profile after JEC applied
// - put bin content to RMS extracted from profile width divided by JES~1
// - update error to error on the mean over JER
void evalJER(TH1D *h, TProfile *p) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double jer = h->GetBinError(i); // RMS with option "S"
    p->SetErrorOption("");
    double ejer = (p ? p->GetBinError(i) : 0); // RMS/sqrt(N) with option ""
    double jes = (p ? p->GetBinContent(i) : 1);
    if (jer!=0 && jes!=0) {
      h->SetBinContent(i, jer/jes);
      h->SetBinError(i, ejer/(jer/jes));
    }
  }
  return;
}

// Evaluate JET counts
// - get jet counts from profile entries per bin
// - normalize by bin width in eta and pT
void evalJET(TH1D *h, const TProfile *p) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double n = p->GetBinEntries(i);
    if (n!=0) {
      double effn = p->GetBinEffectiveEntries(i);
      double en = (effn>0 ? n/sqrt(effn) : 0);
      double dpt = h->GetBinWidth(i);
      double deta = 2*0.087;
      h->SetBinContent(i, n / dpt / deta);
      h->SetBinError(i, en / dpt / deta);
    }
  }
  return;
}

// Evaluate JEC closure
// - turn ratio to percentage difference
void evalJEC(TH1D *h) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double jec = h->GetBinContent(i);
    if (jec!=0) {
      double ejec = h->GetBinError(i);
      h->SetBinContent(i, 100.*(jec-1));
      h->SetBinError(i, 100.*ejec);
    }
  }
  return;
}

// Low pT bias correction for JES
// - fill missing part with Gaussian extension
// - efficiency as profile to get accurate value per bin
// - JER as fit to extrapolate robustly to biased region
//   (also to have eventual JES consistent with JER)
void fixJESandJER(TH1D *hjes, TH1D *hjec, TH1D *hjer, TH1D *heff,
		  const TF1 *f1jer) {

  for (int i = 1; i != hjes->GetNbinsX()+1; ++i) {

    double pt = hjes->GetBinCenter(i);
    double ptmin = hjes->GetBinLowEdge(i);
    double ptmax = hjes->GetBinLowEdge(i+1);
    double jes_raw = hjes->GetBinContent(i);
    double ejes = hjes->GetBinError(i);
    double jer_raw = hjer->GetBinContent(i);
    double ejer = hjer->GetBinError(i);
    double eff = heff->GetBinContent(i);
    double eeff = heff->GetBinError(i);
    // Cap jer at 100%
    double jer_fit = min(1.0, f1jer->Eval(pt));

    double jec_raw = hjec->GetBinContent(i);
    double ejec = hjec->GetBinError(i);
    
    // Only run the extra correction if efficiency is less than 99%
    double jes(jes_raw), jec(jec_raw), jer(jer_raw), eff_cut(eff);
    //if (eff<0.16) { // 1-sided 1-sigma cut
    //if (eff<0.5) { // cut at peak
    if (eff<0.25) { // 
      jes = ejes = jec = ejec = jer = ejer = eff_cut = eeff = 0;
    }
    else if (eff<0.999) {
      
      // Starting points:
      // 1) we know the relative width sigma/mu of the Gaussian from JER (f1jer)
      // 2) we know what fraction of the Gaussian is lost from EFF (heff)
      // 3) we don't know unbiased mean from JES (hjes), but can start with mu=1
      // => estimate zcut that integrates to eff, 
      //    then solve for the mean and sigma of the truncated Gaussian
      // Math v1: https://chatgpt.com/share/67efebc4-0d30-800f-a459-71796a13bdb5
      double mu(1), sigma(1);
      double zcut = TMath::ErfInverse(1. - 2.*eff) * sqrt(2.);
      double phi = TMath::Gaus(zcut, 0, 1, true); // normalized
      double mu_bias = 1 + jer_fit * phi / eff;
      double sigma_bias = sqrt(1 + zcut * phi / eff - pow(phi / eff, 2))
	/ mu_bias;

      // Correct jes_raw based on mu/mu_raw
      jes = (mu_bias>0 ? jes_raw * mu / mu_bias : jes_raw);
      jec = (mu_bias>0 ? jec_raw * mu / mu_bias : jec_raw);
      jer = (sigma_bias>0  && mu_bias>0 ?
	     (jer_raw * sigma / sigma_bias) : jer_raw);

      // Rough place-holder errors for now
      ejes = sqrt(pow(ejes,2)+pow(jes_raw*(mu/mu_bias-1)*0.5,2));
      ejec = sqrt(pow(ejec,2)+pow(jec_raw*(mu/mu_bias-1)*0.5,2));
      ejer = sqrt(pow(ejer,2)+pow(jer_raw*(sigma/sigma_bias-1)*0.5,2));

      if (debug)
	cout << Form("pt=%1.0f [%1.0f,%1.0f] GeV, eff=%1.3f, jer=%1.3f, "
		     "zcut=%1.3f, phi=%1.3f, mu_bias=%1.3f\n"
		     "* jes_raw=%1.3f, jes=%1.3f, "
		     "* jer_raw=%1.3f, jer=%1.3f\n",
		     pt, ptmin, ptmax, eff, jer, zcut, phi, mu_bias,
		     jes_raw, jes, jer_raw, jer);
    }
    heff->SetBinContent(i, eff_cut);
    heff->SetBinError(i, eeff);
    hjes->SetBinContent(i, jes);
    hjes->SetBinError(i, ejes);
    hjec->SetBinContent(i, jec);
    hjec->SetBinError(i, ejec);
    hjer->SetBinContent(i, jer);
    hjer->SetBinError(i, ejer);
  } // for i
  return;
}


void fixMedian(TH1D *hjes, TF1* jes_fit, TF1* jer_fit, double pT_threshold) {

  //double GetMedianBias(double pT, TF1* jes_fit, TF1* jer_fit, double pT_threshold) {
  for (int i = 1; i != hjes->GetNbinsX()+1; ++i) {

    double pT = hjes->GetBinCenter(i);
    
    double jes_mean = jes_fit->Eval(pT);
    double jer = jer_fit->Eval(pT);
    double jes_sigma = jes_mean * jer;

    double jes_min = jes_mean * (pT_threshold / pT);
    double x = (jes_min - jes_mean) / jes_sigma;

    double Fmin = 0.5 * (1.0 + TMath::Erf(x / sqrt(2)));
    double F_median = 0.5 * (1.0 + Fmin);
    double quantile = jes_mean + jes_sigma * sqrt(2) * TMath::ErfInverse(2 * F_median - 1);

    double bias = (quantile - jes_mean) / jes_mean;
    //return bias;
    
    double jes_median_raw = hjes->GetBinContent(i);
    double jes_median_corr = jes_median_raw - bias;
    if (jes_median_raw>0) {
      hjes->SetBinContent(i, jes_median_corr);
    }
  }
  return;
} // fixMedian

double evalXmin(double jes, double jer, double eff) {
  double zcut = TMath::ErfInverse(1. - 2.*eff) * sqrt(2.);
  return (jes*(1 + zcut*jer));
}

TH1D* cutHist(const TH1D *h, const double xmin_out) {
  TH1D *ho = (TH1D*)h->Clone(Form("%s_cut",h->GetName()));
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double xmin = h->GetBinLowEdge(i);
    double xmax = h->GetBinLowEdge(i+1);
    if (xmax<=xmin_out) {
      ho->SetBinContent(i, 0.);
      ho->SetBinError(i, 0.);
    }
  }
  return ho;
}
