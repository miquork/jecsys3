// Purpose: Draw pairs (or multiplets) of time scans made by drawTimeStability.C
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include <set>

//#include "drawTimeStability.C" // BuildRanges, addLines
#include "../tdrstyle_mod22.C"

// pr50pr110 for RCF application
bool useRCFstyle = true;

// Clean out bins from hd that are empty for ha or hb
void cleanEmpty(TH1D *hd, TH1D *ha, TH1D *hb);
// Clean out bins outside given range
void cleanRange(TH1D *h, double xmin, double xmax);
// Find scale from hscales
double getScale(TH1D *hscales, string name);
// Draw Gamma scale vs Zmm scale
void drawGamVsZmm(string mode);
// Draw Central vs Outer Barrel photon scale
void drawGamVsGam();
// Draw TTBar vs Gam
void drawGamVsTTBar(string mode);
// Draw PFcomposition + JES (+ reference run if not 0)
void drawPFcomp(string ref, double refrun=0);

// To-do:
// drawDBvsMPF to check stability of HDM inputs wrt each other

// Main call for drawing various pairs
void drawTimeStabilityPairs() {

  drawGamVsZmm("MPF");
  drawGamVsZmm("DB");

  drawPFcomp("pr50pr110");
  /*
  drawPFcomp("pr50");
  drawPFcomp("pr110");
  drawPFcomp("pr230");
  */
  /*
  drawPFcomp("zpt30");
  drawPFcomp("zpt50");
  drawPFcomp("zpt110");
  */
  //drawGamVsGam();

  //drawGamVsTTBar("prof_top_inWindow");
  //drawGamVsTTBar("prof_W_inWindow");
}

void drawGamVsZmm(string mode) {

  const char *cm = mode.c_str();
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  //double xmin(0), xmax(175), ymin(-2.4), ymax(+3.4), ymind(-1.5), ymaxd(+1.9);
  double xmin(0), xmax(280), ymin(-2.8), ymax(+5.8), ymind(-2.5), ymaxd(+2.5);
  TH1D *h = tdrHist("h_u","JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  TH1D *h_d = tdrHist("h_d","#gamma+j - Z+j (%)",ymind,ymaxd,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  //lumi_136TeV = "Run 3, 2022-24";
  lumi_136TeV = "Run 3, 2022-25";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h,h_d,8,11);

  TCanvas *c1b(0);
  TLegend *leg1b(0);
  if (useRCFstyle) {
    //TH1D *h1b = tdrHist("h1b","(#gamma+j - Z+j) vs late 2024 [V2.0] (%)",ymind+1e-3,ymaxd+0.5-1e-3,
    TH1D *h1b = tdrHist("h1b","(#gamma+j - Z+j) vs 2024 (%)",ymind+1e-3,ymaxd+0.5-1e-3,
			"Cumulative luminosity (fb^{-1})",xmin,xmax);
    c1b = tdrCanvas("c1b",h1b,8,11,kSquare);

    l->SetLineStyle(kDotted);
    l->DrawLine(xmin,-0.5,xmax,-0.5);
    l->DrawLine(xmin,+0.5,xmax,+0.5);
    l->SetLineStyle(kDashed);
    l->DrawLine(xmin,0,xmax,0);

    leg1b = tdrLeg(0.58,0.90-3*0.05,0.83,0.90);
  }

  
  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);
  //map<string,<pair,pair> ml = BuildRanges("rootfiles/brilcalc/fibs.txt");
  //addLines(ml, hcumlum2, 1.0, xmin, xmax, ymin, ymax, 0);

  TH1D *ha(0), *hb(0);
  if (mode=="MPF") {
    ha = (TH1D*)f->Get("pr50m"); assert(ha);
    hb = (TH1D*)f->Get("mpf_run_zpt50"); assert(hb);
  }
  if (mode=="DB") {
    ha = (TH1D*)f->Get("pr50b"); assert(ha);
    hb = (TH1D*)f->Get("db_run_zpt50"); assert(hb);
  }
  TH1D *hd = (TH1D*)ha->Clone("hd");
  hd->Add(hb,-1);
  cleanEmpty(hd,ha,hb);

  //cleanRange(ha,63.2,200); // noisy
  // //cleanRange(hb,63.2,200); // ok
  //cleanRange(hd,63.2,200); // noisy
  
  // Extra pairs
  TH1D *ha2(0), *hb2(0);
  if (mode=="MPF") {
    ha2 = (TH1D*)f->Get("pr110m"); assert(ha2);
    hb2 = (TH1D*)f->Get("mpf_run_zpt110"); assert(hb2);
  }
  if (mode=="DB") {
  ha2 = (TH1D*)f->Get("pr110b"); assert(ha2);
  hb2 = (TH1D*)f->Get("db_run_zpt110"); assert(hb2);
  }
  TH1D *hd2 = (TH1D*)ha2->Clone("hd2");
  hd2->Add(hb2,-1);
  cleanEmpty(hd2,ha2,hb2);

  //cleanRange(ha2,0,63.2);
  //cleanRange(hb2,0,63.2);
  //cleanRange(hd2,0,63.2);
  
  c1->cd(1);

  TLine *ll = new TLine();
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();
  tex->DrawLatex(0.40,0.80,cm);
  tex->SetNDC(kFALSE);
  tex->SetTextSize(0.03);
  set<double> breakset;
  breakset.insert(xmin);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));
    /*
    ll->SetLineStyle(kSolid);
    ll->SetLineColor(kGray);
    if (ts.Contains("V")) { ll->SetLineColor(kGreen+1); } //kRed-9); //continue; }
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    if (ts.Contains("TR")) { ll->SetLineColor(kOrange+1); continue; }
    if (ll->GetLineColor()==kBlue-9)
      ll->DrawLine(cumlum,ymin,cumlum,ymax);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kGreen+1);//kRed-9);
    if (ts.Contains("IC")) tex->SetTextColor(kBlue-9);
    if (ts.Contains("TR")) tex->SetTextColor(kOrange+1);
    if (tex->GetTextColor()==kBlue-9)
      tex->DrawLatex(cumlum+1,ymin+0.4-((i-1)%3)*0.15,ts.Data());

    if (useRCFstyle) {
      c1b->cd();
      if (ts.Contains("2") && ll->GetLineColor()==kGray || ts=="V2.0") {
	ll->DrawLine(cumlum,ymind,cumlum,ymaxd+0.5);
	tex->DrawLatex(cumlum+1,ymind+0.2,ts.Data());
      }
      if (ll->GetLineColor()==kGray) {
	ll->SetLineStyle(kDotted);
	ll->DrawLine(cumlum,ymind,cumlum,ymaxd+0.5);
      }
      c1->cd(1);
    }
    */
    
    //if (ts.Contains("IC"))
      breakset.insert(cumlum);
  }
  breakset.insert(xmax);
  
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+1,xmax,+1);
  l->DrawLine(xmin,-1,xmax,-1);
  
  tdrDraw(hb,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.6);
  tdrDraw(hb2,"Pz",kOpenSquare,kMagenta-9,kSolid,-1,kNone,0,0.6);
  tdrDraw(ha2,"Pz",kOpenCircle,kCyan+1-9,kSolid,-1,kNone,0,0.6);
  tdrDraw(ha,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  TLegend *leg = tdrLeg(0.65,0.89-0.05*4,0.90,0.89);
  leg->SetFillStyle(1001);
  leg->AddEntry(hb,"Z(#mu#mu)+jet 50","PLE");
  leg->AddEntry(hb2,"Z(#mu#mu)+jet 110","PLE");
  leg->AddEntry(ha,"#gamma+jet 50","PLE");
  leg->AddEntry(ha2,"#gamma+jet 110","PLE");

  gPad->RedrawAxis();
  gPad->Update();
  
  c1->cd(2);

  tex->SetTextSize(0.03*1.5);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    ll->SetLineColor(kGray);
    if (ts.Contains("V")) { ll->SetLineColor(kRed-9); continue; }
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    if (ts.Contains("TR")) { ll->SetLineColor(kOrange+1); continue; }
    ll->DrawLine(cumlum,ymind,cumlum,ymaxd);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kRed-9);
    if (ts.Contains("IC")) {
      tex->SetTextColor(kBlue-9);
      //tex->DrawLatex(cumlum+1,ymind+0.5,ts.Data());
      tex->DrawLatex(cumlum+1,ymind+0.5-((i-1)%3)*0.20,ts.Data());
    }
    if (ts.Contains("TR")) {
      tex->SetTextColor(kOrange+1);
      tex->DrawLatex(cumlum+1,ymind+0.5,ts.Data());
    }
    //tex->DrawLatex(cumlum+1,ymind+0.4-((i-1)%3)*0.15,ts.Data());
  }

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+0.5,xmax,+0.5);
  l->DrawLine(xmin,-0.5,xmax,-0.5);

  tdrDraw(hd2,"Pz",kOpenSquare,kMagenta-9,kSolid,-1,kNone,0,0.6);
  tdrDraw(hd,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  if (useRCFstyle) {
    c1b->cd();
    TLatex *tex = new TLatex(); tex->SetNDC(); tex->SetTextSize(0.045);
    tex->DrawLatex(0.40,0.80,cm);
    tdrDraw(hd2,"Pz",kOpenSquare,kMagenta-9,kSolid,-1,kNone,0,0.6);
    tdrDraw(hd,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);
    leg1b->AddEntry(hd2,"p_{T} > 110 GeV","PLE");
    leg1b->AddEntry(hd,"p_{T} > 50 GeV","PLE");
    gPad->RedrawAxis();
    c1->cd(2);
  }
  
  // Fit nibs
  TF1 *f1 = new TF1("f1","[0]",xmin,xmax);
  f1->SetLineColor(kBlue);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(new TGraphErrors(hd));
  mg->Add(new TGraphErrors(hd2));
  //mg->Draw("SAMEHIST");
  typedef set<double>::const_iterator IT;
  double xfit1(0),xfit2(0);
  double vx[breakset.size()+2];
  double vy[breakset.size()+1];
  double vye[breakset.size()+1];
  int nib(0);
  vx[0] = xmin;
  for (IT it = breakset.begin(); it != breakset.end(); ++it) {
    xfit1 = xfit2;
    xfit2 = (*it);
    f1->SetRange(xfit1,xfit2);
    mg->Fit(f1,"QRN");
    f1->DrawClone("SAME");
    if (xfit2!=xfit1) {
      vx[nib+1] = xfit2;
      vy[nib] = f1->GetParameter(0);
      vye[nib] = f1->GetParError(0);
      ++nib;
    }
    if (useRCFstyle) {
      c1b->cd();
      f1->DrawClone("SAME");
      c1->cd(2);
    }
  }
  TH1D *hnib = new TH1D(Form("GamJetMinusZmmJet_%s",cm),"",nib,&vx[0]);
  for (int i = 0; i != nib; ++i) {
    hnib->SetBinContent(i+1, vy[i]);
    hnib->SetBinError(i+1, vye[i]);
  }
  tdrDraw(hnib,"E2",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
  hnib->SetFillColorAlpha(kCyan+1,0.5);
  
  gPad->RedrawAxis();

  if (useRCFstyle) {
      c1b->cd();
      tdrDraw(hnib,"E2",kNone,kBlue,kSolid,-1,1001,kCyan+1);
      hnib->SetFillColorAlpha(kCyan+1,0.5);
      leg1b->AddEntry(hnib,"Fit","FL");
      c1->cd(2);
      gPad->RedrawAxis();
  }

  c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStabilityPairs_GamVsZmm_%s.pdf",cm));

  // Save photon correction
  cout << "Save "<<Form("GamJetVsZmmJet_%s",cm)<<" to rootfiles/drawTimeStabilityPairs.root" << endl;
  TFile *fout = new TFile("rootfiles/drawTimeStabilityPairs.root","UPDATE");
  hnib->Write(Form("GamJetMinusZmmJet_%s",cm),TObject::kOverwrite);
  fout->Close();

  // RCF2025 plots
  if (useRCFstyle) {
    c1b->SaveAs(Form("pdf/drawTimeStability/drawTimeStabilityPairs_RCF2025_GamVsZmm_%s.pdf",cm));
  }
} // drawGamVsZmm

void drawGamVsGam() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  double xmin(0), xmax(175), ymin(-2.4), ymax(+3.4), ymind(-1.5), ymaxd(+1.9);
  TH1D *h = tdrHist("h_u","JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  TH1D *h_d = tdrHist("h_d","#gamma+j - Z+j (%)",ymind,ymaxd,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  lumi_136TeV = "Run3, 2022-24";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h,h_d,8,11);

  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);
  TH1D *hscales = (TH1D*)f->Get("hscales"); assert(hscales);

  double jeshi(1), jeslo(1);
  for (int i = 1; i != hscales->GetNbinsX()+1; ++i) {
    string s = hscales->GetXaxis()->GetBinLabel(i);
    if (s=="pr50m_eta08hi") jeshi = hscales->GetBinContent(i);
    if (s=="pr50m_eta08lo") jeslo = hscales->GetBinContent(i);
  }
  
  TH1D *ha = (TH1D*)f->Get("pr50m_eta08lo"); assert(ha);
  ha->Scale(jeslo);
  TH1D *hb = (TH1D*)f->Get("pr50m_eta08hi"); assert(hb);
  hb->Scale(jeshi);
  TH1D *hd = (TH1D*)ha->Clone("hd");
  hd->Add(hb,-1);
  cleanEmpty(hd,ha,hb);

  c1->cd(1);

  TLine *ll = new TLine();
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.03);
  set<double> breakset;
  breakset.insert(xmin);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    ll->SetLineColor(kGray);
    if (ts.Contains("V")) ll->SetLineColor(kRed-9);
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    ll->DrawLine(cumlum,ymin,cumlum,ymax);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kRed-9);
    if (ts.Contains("IC")) tex->SetTextColor(kBlue-9);
    tex->DrawLatex(cumlum+1,ymin+0.4-((i-1)%3)*0.15,ts.Data());

    breakset.insert(cumlum);
  }
  breakset.insert(xmax);
  
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+1,xmax,+1);
  l->DrawLine(xmin,-1,xmax,-1);
  
  tdrDraw(hb,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.6);
  tdrDraw(ha,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  TLegend *leg = tdrLeg(0.60,0.89-0.05*2,0.85,0.89);
  leg->SetFillStyle(1001);
  leg->AddEntry(hb,"#gamma+jet 50 |#eta_{#gamma}|>0.8","PLE");
  leg->AddEntry(ha,"#gamma+jet 50 |#eta_{#gamma}|<0.8","PLE");

  gPad->RedrawAxis();
  gPad->Update();
  
  c1->cd(2);

  tex->SetTextSize(0.03*1.5);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    ll->SetLineColor(kGray);
    if (ts.Contains("V")) ll->SetLineColor(kRed-9);
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    ll->DrawLine(cumlum,ymind,cumlum,ymaxd);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kRed-9);
    if (ts.Contains("IC")) {
      tex->SetTextColor(kBlue-9);
      tex->DrawLatex(cumlum+1,ymind+0.5,ts.Data());
    }
    //tex->DrawLatex(cumlum+1,ymind+0.4-((i-1)%3)*0.15,ts.Data());
  }

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+0.5,xmax,+0.5);
  l->DrawLine(xmin,-0.5,xmax,-0.5);

  tdrDraw(hd,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  // Fit nibs
  TF1 *f1 = new TF1("f1","[0]",xmin,xmax);
  f1->SetLineColor(kBlue);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(new TGraphErrors(hd));
  //mg->Draw("SAMEHIST");
  typedef set<double>::const_iterator IT;
  double xfit1(0),xfit2(0);
  double vx[breakset.size()+2];
  double vy[breakset.size()+1];
  double vye[breakset.size()+1];
  int nib(0);
  vx[0] = xmin;
  for (IT it = breakset.begin(); it != breakset.end(); ++it) {
    xfit1 = xfit2;
    xfit2 = (*it);
    f1->SetRange(xfit1,xfit2);
    mg->Fit(f1,"QRN");
    f1->DrawClone("SAME");
    if (xfit2!=xfit1) {
      vx[nib+1] = xfit2;
      vy[nib] = f1->GetParameter(0);
      vye[nib] = f1->GetParError(0);
      ++nib;
    }
  }
  TH1D *hnib = new TH1D("hnib","",nib,&vx[0]);
  for (int i = 0; i != nib; ++i) {
    hnib->SetBinContent(i+1, vy[i]);
    hnib->SetBinError(i+1, vye[i]);
  }
  tdrDraw(hnib,"E2",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
  hnib->SetFillColorAlpha(kCyan+1,0.5);
  
  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/drawTimeStability/drawTimeStabilityPairs_GamVsGam.pdf");
} // drawTimeStabilityPairs()

void drawGamVsTTBar(string mode) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  //double xmin(0), xmax(175), ymin(-2.4), ymax(+3.4), ymind(-1.5), ymaxd(+1.9);
  double xmin(0), xmax(175), ymin(-2.4), ymax(+3.4), ymind(-1.9), ymaxd(+2.4);
  TH1D *h = tdrHist("h_u","JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  TH1D *h_d = tdrHist("h_d","#gamma+j - t#bar{t} (%)",ymind,ymaxd,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  lumi_136TeV = "Run 3, 2022-24";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1",h,h_d,8,11);

  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);

  TH1D *ha(0), *hb(0);
  ha = (TH1D*)f->Get("pr50m"); assert(ha);
  //hb = (TH1D*)f->Get("prof_top_inWindow"); assert(hb);
  //hb = (TH1D*)f->Get("prof_W_inWindow"); assert(hb);
  const char *cm = mode.c_str();
  TString tm(cm);
  hb = (TH1D*)f->Get(cm); assert(hb);
  TH1D *hd = (TH1D*)ha->Clone("hd");
  hd->Add(hb,-1);
  cleanEmpty(hd,ha,hb);

  //cleanRange(ha,63.2,200); // noisy
  //cleanRange(hb,63.2,200); // ok
  //cleanRange(hd,63.2,200); // noisy
  
  c1->cd(1);

  TLine *ll = new TLine();
  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC(kFALSE);
  tex->SetTextSize(0.03);
  set<double> breakset;
  breakset.insert(xmin);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    ll->SetLineStyle(kSolid);
    ll->SetLineColor(kGray);
    if (ts.Contains("V")) { ll->SetLineColor(kGreen+1); } //kRed-9); //continue; }
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    if (ts.Contains("TR")) { ll->SetLineColor(kOrange+1); continue; }
    if (ll->GetLineColor()==kBlue-9)
      ll->DrawLine(cumlum,ymin,cumlum,ymax);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kGreen+1);//kRed-9);
    if (ts.Contains("IC")) tex->SetTextColor(kBlue-9);
    if (ts.Contains("TR")) tex->SetTextColor(kOrange+1);
    if (tex->GetTextColor()==kBlue-9)
      tex->DrawLatex(cumlum+1,ymin+0.4-((i-1)%3)*0.15,ts.Data());

    if (ts.Contains("IC"))
      breakset.insert(cumlum);
  }
  breakset.insert(xmax);
  
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+1,xmax,+1);
  l->DrawLine(xmin,-1,xmax,-1);
  
  tdrDraw(hb,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.6);
  tdrDraw(ha,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  TLegend *leg = tdrLeg(0.65,0.89-0.05*4,0.90,0.89);
  leg->SetFillStyle(1001);

  if (tm.Contains("top"))
    leg->AddEntry(hb,"t#bar{t} m_{t}","PLE");
  else
    leg->AddEntry(hb,"t#bar{t} m_{W}","PLE");
  leg->AddEntry(ha,"#gamma+jet 50","PLE");

  gPad->RedrawAxis();
  gPad->Update();
  
  c1->cd(2);

  tex->SetTextSize(0.03*1.5);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    ll->SetLineColor(kGray);
    if (ts.Contains("V")) { ll->SetLineColor(kRed-9); continue; }
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    if (ts.Contains("TR")) { ll->SetLineColor(kOrange+1); continue; }
    ll->DrawLine(cumlum,ymind,cumlum,ymaxd);

    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kRed-9);
    if (ts.Contains("IC")) {
      tex->SetTextColor(kBlue-9);
      //tex->DrawLatex(cumlum+1,ymind+0.5,ts.Data());
      tex->DrawLatex(cumlum+1,ymind+0.5-((i-1)%3)*0.20,ts.Data());
    }
    if (ts.Contains("TR")) {
      tex->SetTextColor(kOrange+1);
      tex->DrawLatex(cumlum+1,ymind+0.5,ts.Data());
    }
    //tex->DrawLatex(cumlum+1,ymind+0.4-((i-1)%3)*0.15,ts.Data());
  }

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+0.5,xmax,+0.5);
  l->DrawLine(xmin,-0.5,xmax,-0.5);

  tdrDraw(hd,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

  // Fit nibs
  TF1 *f1 = new TF1("f1","[0]",xmin,xmax);
  f1->SetLineColor(kBlue);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(new TGraphErrors(hd));
  //mg->Draw("SAMEHIST");
  typedef set<double>::const_iterator IT;
  double xfit1(0),xfit2(0);
  double vx[breakset.size()+2];
  double vy[breakset.size()+1];
  double vye[breakset.size()+1];
  int nib(0);
  vx[0] = xmin;
  for (IT it = breakset.begin(); it != breakset.end(); ++it) {
    xfit1 = xfit2;
    xfit2 = (*it);
    f1->SetRange(xfit1,xfit2);
    mg->Fit(f1,"QRN");
    f1->DrawClone("SAME");
    if (xfit2!=xfit1) {
      vx[nib+1] = xfit2;
      vy[nib] = f1->GetParameter(0);
      vye[nib] = f1->GetParError(0);
      ++nib;
    }
  }
  TH1D *hnib = new TH1D("GamJetMinusTTBar","",nib,&vx[0]);
  for (int i = 0; i != nib; ++i) {
    hnib->SetBinContent(i+1, vy[i]);
    hnib->SetBinError(i+1, vye[i]);
  }
  tdrDraw(hnib,"E2",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
  hnib->SetFillColorAlpha(kCyan+1,0.5);
  
  gPad->RedrawAxis();

  if (tm.Contains("top"))
    c1->SaveAs("pdf/drawTimeStability/drawTimeStabilityPairs_GamVsTTBar_mt.pdf");
  else
    c1->SaveAs("pdf/drawTimeStability/drawTimeStabilityPairs_GamVsTTBar_mW.pdf");

} // drawGamVsTTBar


void drawPFcomp(string ref, double refrun) {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  bool isRCFstyle = (useRCFstyle && ref=="pr50pr110");
  const char *cr = (ref=="pr50pr110" ? "pr50" : ref.c_str());
  TString tr(cr);  

  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  // Photon scale from drawGamVsZmm()
  TFile *fs = new TFile("rootfiles/drawTimeStabilityPairs.root","READ");
  assert(fs && !fs->IsZombie());
  
  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  //double xmin(0), xmax(175), ymin(-3.7), ymax(+6.3);
  double xmin(0), xmax(280), ymin(-9.0), ymax(+9.5);
  TH1D *h = tdrHist(Form("hpf_%s",cr),"PFJES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  //lumi_136TeV = "Run3, 2022-24";
  lumi_136TeV = "Run3, 2022-25";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1pf_%s",cr),h,8,11,kRectangular);

  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);
  TH1D *hscales = (TH1D*)f->Get("hscales"); assert(hscales);

  // DB has bigger effect than MPF in 2022-23. Good?
  //TH1D *hgamvsz = (TH1D*)fs->Get("GamJetMinusZmmJet_DB"); assert(hgamvsz);
  // MPF seems more continous for pr50m. Is any correction needed?
  TH1D *hgamvsz = (TH1D*)fs->Get("GamJetMinusZmmJet_MPF"); assert(hgamvsz);
  // NONE with no correction at all. MPF looks better.
  //hgamvsz->Reset();

  // For MPF use MPF correction
  TH1D *hgamvszm = (TH1D*)fs->Get("GamJetMinusZmmJet_MPF"); assert(hgamvszm);

  // Format string for Z+jet
  TH1D *hr(0), *hc(0), *hn(0), *he(0), *hm(0);
  double kr(1), fc(0.65), fe(0.25), fn(0.10), km(1);
  if (tr.Contains("pr")) {
    //hr = (TH1D*)f->Get(Form("%sm_jes",cr)); assert(hr); // MPF
    hr = (TH1D*)f->Get(Form("%sb_jes",cr)); assert(hr); // DB
    hc = (TH1D*)f->Get(Form("%schf",cr));  assert(hc);
    hn = (TH1D*)f->Get(Form("%snhf",cr));  assert(hn);
    he = (TH1D*)f->Get(Form("%snef",cr));  assert(he);
    hm = (TH1D*)f->Get(Form("%sm_jes",cr)); assert(hr); // MPF

    //kr = getScale(hscales,Form("%sm_jes",cr)); // MPF
    kr = getScale(hscales,Form("%sb_jes",cr)); // DB
    fc = getScale(hscales,Form("%schf",cr));
    fe = getScale(hscales,Form("%snef",cr));
    fn = getScale(hscales,Form("%snhf",cr));
    km = getScale(hscales,Form("%sm_jes",cr)); // MPF
  }
  if (tr.Contains("zpt")) {
    hr = (TH1D*)f->Get(Form("mpf_run_%s_jes",cr)); assert(hr);
    hc = (TH1D*)f->Get(Form("chf_run_%s",cr));  assert(hc);
    hn = (TH1D*)f->Get(Form("nhf_run_%s",cr));  assert(hn);
    he = (TH1D*)f->Get(Form("nef_run_%s",cr));  assert(he);

    kr = getScale(hscales,Form("mpf_run_%s",cr));
    fc = getScale(hscales,Form("chf_run_%s",cr));
    fe = getScale(hscales,Form("nef_run_%s",cr));
    fn = getScale(hscales,Form("nhf_run_%s",cr));
  }
  assert(hr);

  // Get second photon trigger to patch full range
  TH1D *hr2(0), *hc2(0), *hn2(0), *he2(0), *hm2(0);
  double kr2(1), fc2(0.65), fe2(0.25), fn2(0.10), km2(0);
  if (ref=="pr50pr110") {
    const char *cr2 = "pr110";
    //hr2 = (TH1D*)f->Get(Form("%sm_jes",cr2)); assert(hr2); // MPF
    hr2 = (TH1D*)f->Get(Form("%sb_jes",cr2)); assert(hr2); // DB
    hc2 = (TH1D*)f->Get(Form("%schf",cr2));  assert(hc2);
    hn2 = (TH1D*)f->Get(Form("%snhf",cr2));  assert(hn2);
    he2 = (TH1D*)f->Get(Form("%snef",cr2));  assert(he2);
    hm2 = (TH1D*)f->Get(Form("%sm_jes",cr2)); assert(hm2); // MPF

    //kr2 = getScale(hscales,Form("%sm_jes",cr2)); // MPF
    kr2 = getScale(hscales,Form("%sb_jes",cr2)); // DB
    fc2 = getScale(hscales,Form("%schf",cr2));
    fe2 = getScale(hscales,Form("%snef",cr2));
    fn2 = getScale(hscales,Form("%snhf",cr2));
    km2 = getScale(hscales,Form("%sm_jes",cr2)); // MPF
  }

  // Calculate JES and scale with per-fib JES
  TH1D *hjesr = (TH1D*)hr->Clone(Form("hjesr_%s",cr));
  TH1D *hjesc = (TH1D*)hc->Clone(Form("hjesc_%s",cr));
  TH1D *hjesn = (TH1D*)hn->Clone(Form("hjesn_%s",cr));
  TH1D *hjese = (TH1D*)he->Clone(Form("hjese_%s",cr));
  TH1D *hmpfr = (TH1D*)hm->Clone(Form("hmpfr_%s",cr));
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {

    TH1D *hrs(hr), *hcs(hc), *hns(hn), *hes(he), *hms(hm);
    double krs(kr), fcs(fc), fns(fn), fes(fe), kms(km);
    if (ref=="pr50pr110" && hr->GetBinLowEdge(i)<63.2) {
      assert(hr2); assert(hc2); assert(hn2); assert(he2); assert(hm2);
      hrs = hr2; hcs = hc2; hns = hn2; hes = he2; hms = hm2;
      krs = kr2; fcs = fc2; fns = fn2; fes = fe2; kms = km2;
    }
    
    double s = hrs->GetBinContent(i); // s = (jes-1)*100
    double es = hrs->GetBinError(i);
    double jes = s*kr*0.01+1;

    double m = hms->GetBinContent(i); // m = (mpf-1)*100
    double em = hms->GetBinError(i);
    double mpf = m*km*0.01+1;

    // Patch photon scale before 2024F from GamVsZmm_DB
    double cumlum = hrs->GetBinCenter(i);
    double j = hgamvsz->GetXaxis()->FindBin(cumlum);
    double gamvsz = (1+0.01*hgamvsz->GetBinContent(j));
    jes *= 1./gamvsz;

    double jm = hgamvszm->GetXaxis()->FindBin(cumlum);
    double gamvszm = (1+0.01*hgamvszm->GetBinContent(jm));
    mpf *= 1./gamvszm;
    
    double jesr = (jes-1)*100.;
    double ejesr = es*krs;
    hjesr->SetBinContent(i, jesr);
    hjesr->SetBinError(i, ejesr);
    double c = hcs->GetBinContent(i);
    double ec = hcs->GetBinError(i);
    double jesc = ((1+c*0.01)*fc*jes-fc)*100;
    double ejesc = ec*fcs*jes;
    hjesc->SetBinContent(i, jesc);
    hjesc->SetBinError(i, ejesc);
    double n = hns->GetBinContent(i);
    double en = hns->GetBinError(i);
    double jesn = ((1+n*0.01)*fn*jes-fn)*100;
    double ejesn = en*fns*jes;
    hjesn->SetBinContent(i, jesn);
    hjesn->SetBinError(i, ejesn);
    double e = hes->GetBinContent(i);
    double ee = hes->GetBinError(i);
    double jese = ((1+e*0.01)*fe*jes-fe)*100;
    double ejese = ee*fes*jes;
    hjese->SetBinContent(i, jese);
    hjese->SetBinError(i, ejese);
    double mpfr = (mpf-1)*100.;
    double empfr = em*kms;
    hmpfr->SetBinContent(i, mpfr);
    hmpfr->SetBinError(i, empfr);
  }

  // Add known breaks
  TLine *ll = new TLine();
  TLatex *tex = new TLatex();
  //tex->SetTextSize(0.03);
  tex->SetTextSize(0.015);
  tex->SetTextAngle(270);
  tex->SetTextAlign(31);
  
  set<double> breakset;
  breakset.insert(xmin);
  double prevlum(0);
  for (int i = 1; i != hbreaks->GetNbinsX()+1; ++i) {
    double cumlum = hbreaks->GetBinContent(i);
    TString ts(hbreaks->GetXaxis()->GetBinLabel(i));

    /*
    ll->SetLineColor(kGray);
    ll->SetLineStyle(kSolid);
    if (ts.Contains("V")) ll->SetLineColor(kGreen+1);
    if (ts.Contains("IC")) ll->SetLineColor(kBlue-9);
    if (ts.Contains("TR")) ll->SetLineColor(kRed-9);
    if (ll->GetLineColor()==kGray || refrun>=0) {
      //if (!isRCFstyle || (ll->GetLineColor()==kGreen+1 || ll->GetLineColor()==kGray))
      //if (!isRCFstyle || (ll->GetLineColor()==kGray && ts.Contains("2") || ts=="V2.0"))
      if (isRCFstyle && !ts.Contains("2")) ll->SetLineStyle(kDotted);
      if (!isRCFstyle || (ll->GetLineColor()==kGray || ts=="V2.0"))
	ll->DrawLine(cumlum,ymin,cumlum,ymax);
    }
	  
    tex->SetTextColor(kGray);
    if (ts.Contains("V")) tex->SetTextColor(kGreen+2);
    if (ts.Contains("IC")) tex->SetTextColor(kBlue);
    if (ts.Contains("TR")) tex->SetTextColor(kRed);
    if (ll->GetLineColor()==kGray || refrun>=0) {
      //if (!isRCFstyle) // || (ll->GetLineColor()==kGreen+1 || ll->GetLineColor()==kGray))
      if (!isRCFstyle || (ll->GetLineColor()==kGray && ts.Contains("2") || ts=="V2.0"))
	tex->DrawLatex(cumlum+1,ymin+0.4-((i-1)%3)*0.15,ts.Data());
    }
    */

    // Avoid label overlaps
    bool isYear = (ts.Length()==4);
    if (cumlum-prevlum>3.0 || isYear) {
      //if (hbreaks) hbreaks->GetXaxis()->SetBinLabel(i,tag.c_str());
      tex->SetTextColorAlpha(kGray,isYear ? 1.0 : 0.3);
      ll->SetLineColorAlpha(kGray,isYear ? 1.0 : 0.1);
      ll->DrawLine(cumlum,ymin,cumlum,ymax);
      tex->DrawLatex(cumlum+0.01,ymin,ts.Data());
      prevlum = cumlum;
    }
    else {
      ll->SetLineColorAlpha(kGray,0.075);
      ll->DrawLine(cumlum,ymin+3.0,cumlum,ymax);
    }
			  
    breakset.insert(cumlum);
  }
  breakset.insert(xmax);

  // Add new break
  if (refrun!=0 && fabs(refrun)<10000) {
    ll->SetLineColor(kBlack);
    ll->DrawLine(refrun,ymin,refrun,ymax);
  }
  
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,+1,xmax,+1);
  l->DrawLine(xmin,-1,xmax,-1);
  
  //tdrDraw(hjesc,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.6);
  //tdrDraw(hjese,"Pz",kOpenSquare,kBlue,kSolid,-1,kNone,0,0.6);
  //tdrDraw(hjesn,"Pz",kOpenCircle,kGreen+2,kSolid,-1,kNone,0,0.6);
  //tdrDraw(hjesr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);

  tdrDraw(hmpfr,"Pz",kFullDiamond,kOrange+2,kSolid,-1,kNone,0,0.6);
  tdrDraw(hjesc,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.5);
  tdrDraw(hjese,"Pz",kOpenSquare,kBlue,kSolid,-1,kNone,0,0.5);
  tdrDraw(hjesn,"Pz",kOpenCircle,kGreen+2,kSolid,-1,kNone,0,0.5);
  tdrDraw(hjesr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.5);

  //TLegend *leg = tdrLeg(0.65,0.89-0.05*4,0.90,0.89);
  TLegend *leg = tdrLeg(0.58,0.89-0.05*5,0.83,0.89);
  leg->SetFillStyle(1001);
  if (isRCFstyle) {
    //h->GetYaxis()->SetTitle("Difference to late 2024 [V2.0] (%)");
    h->GetYaxis()->SetTitle("Difference to 2024 (%)");
    h->GetYaxis()->SetTitleOffset(0.8);
    gPad->SetLeftMargin(0.10);
    leg->SetX1(0.43); leg->SetX2(0.68);
    //leg->AddEntry(hjesr,"Jet Energy Scale, split by:","PLE");
    //leg->AddEntry(hjesc,"  Tracking (charged hadrons)","PLE");
    //leg->AddEntry(hjesn,"  HCAL (neutral hadrons)","PLE");
    //leg->AddEntry(hjese,"  ECAL (photons)","PLE");
    leg->AddEntry(hjesr,"Jet Energy Scale, split by PF:","PLE");
    leg->AddEntry(hjesc,"  Charged hadrons (tracks)","PLE");
    leg->AddEntry(hjesn,"  Neutral hadrons (HCAL)","PLE");
    leg->AddEntry(hjese,"  Photons (ECAL)","PLE");
    leg->AddEntry(hmpfr,"Missing E_{T} (MPF)","PLE");
  }
  else if (ref=="pr50pr110") {
    //leg->AddEntry(hjesr,"#gamma+jet MPF 50/110","PLE");
    leg->AddEntry(hjesr,"#gamma+jet DB 50/110","PLE");
    leg->AddEntry(hjesc,"#gamma+jet CHF 50/110","PLE");
    leg->AddEntry(hjesn,"#gamma+jet NHF 50/110","PLE");
    leg->AddEntry(hjese,"#gamma+jet NEF 50/110","PLE");
    leg->AddEntry(hmpfr,"#gamma+jet MPF 50/110","PLE");
  }
  if (ref=="pr50" || ref=="pr110" || ref=="pr230") {
    int ipt; sscanf(cr,"pr%d",&ipt);
    //leg->AddEntry(hjesr,Form("#gamma+jet MPF %d",ipt),"PLE");
    leg->AddEntry(hjesr,Form("#gamma+jet DB %d",ipt),"PLE");
    leg->AddEntry(hjesc,Form("#gamma+jet CHF %d",ipt),"PLE");
    leg->AddEntry(hjesn,Form("#gamma+jet NHF %d",ipt),"PLE");
    leg->AddEntry(hjese,Form("#gamma+jet NEF %d",ipt),"PLE");
    leg->AddEntry(hmpfr,Form("#gamma+jet MPF %d",ipt),"PLE");
  }
  if (ref=="zpt30" || ref=="zpt50" || ref=="zpt110") {
    int ipt; sscanf(cr,"zpt%d",&ipt);
    //leg->AddEntry(hjesr,Form("Z_{#mu#mu}+jet MPF %d",ipt),"PLE");
    leg->AddEntry(hjesr,Form("Z_{#mu#mu}+jet DB %d",ipt),"PLE");
    leg->AddEntry(hjesc,Form("Z_{#mu#mu}+jet CHF %d",ipt),"PLE");
    leg->AddEntry(hjesn,Form("Z_{#mu#mu}+jet NHF %d",ipt),"PLE");
    leg->AddEntry(hjese,Form("Z_{#mu#mu}+jet NEF %d",ipt),"PLE");
    leg->AddEntry(hmpfr,Form("Z_{#mu#mu}+jet MPF %d",ipt),"PLE");
  }
    
  gPad->RedrawAxis();
  gPad->Update();

  if (isRCFstyle) {
    c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStabilityPairs_RCF2025_%s.pdf",ref.c_str()));
  }
  else if (ref=="pr50pr110")
    c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStabilityPairs_PFcomp_%s.pdf",ref.c_str()));
  else
    c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStabilityPairs_PFcomp_%s.pdf",cr));
} // drawTimeStabilityPairs()



void cleanEmpty(TH1D *hd, TH1D *ha, TH1D *hb) {

  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    if ((ha->GetBinContent(i)==0 && ha->GetBinError(i)==0) ||
	(hb->GetBinContent(i)==0 && hb->GetBinError(i)==0)) {
      hd->SetBinContent(i,0);
      hd->SetBinError(i,0);
    }
  }
} // cleanEmpty

void cleanRange(TH1D *h, double xmin, double xmax) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinLowEdge(i)<xmin || h->GetBinLowEdge(i+1)>=xmax) {
      h->SetBinContent(i,0);
      h->SetBinError(i,0);
    }
  }
} // cleanRange


double getScale(TH1D *hscales, string name) {

  for (int i = 1; i != hscales->GetNbinsX()+1; ++i) {
    if (string(hscales->GetXaxis()->GetBinLabel(i))==name) return hscales->GetBinContent(i);
  } // for i
  cout << "getScale not found for " << name << endl << flush;
  return 1;
} // getScale
