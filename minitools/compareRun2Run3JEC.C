// Purpose: compare the JES and the residual actually applied to data and to MC
//          between a Run2 and a Run3 epoch, channel by channel.
//          Run2 data carries Summer20 L2L3 plus Summer19 CHS residuals while
//          the MC carries Summer24, so data/MC is not the same correction on
//          both sides and the difference shows up as a JES offset.
//
//          Reads the pjes*_data, pjes*_mc, pjes*_ratio (and pres*) objects
//          written by reprocess.C into the data/, mc/ and ratio/ directories.
//
// Usage:   root -l -b -q 'minitools/compareRun2Run3JEC.C("2018ABCD","2024_nib")'
#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <iostream>
#include <string>
#include <vector>

#include "../tdrstyle_mod22.C"

bool debug = false;
using namespace std;

// Channels, in the order they are drawn
struct chan_t { const char *id; const char *label; int color;
                int mfull; int mopen; };
const chan_t chans[] = {
  {"",  "#gamma+jet #oplus Z+jet", kGray+2,  kFullSquare,      kOpenSquare},
  {"p", "#gamma+jet",              kBlue+1,   kFullSquare,      kOpenSquare},
  {"z", "Z+jet",                   kRed+1,    kFullCircle,      kOpenCircle},
  {"m", "Multijet",                kBlack, kFullDiamond,    kOpenDiamond},
  {"w", "W#rightarrowqq",          kGreen+2, kFullTriangleDown, kOpenTriangleDown}
};
const int nchan = sizeof(chans)/sizeof(chans[0]);

TH1D *getHist(TFile *f, const char *dir, const char *name, const char *tag);
TH1D *ratioOfHists(TH1D *h1, TH1D *h2, const char *name);
void drawSet(TFile *f1, TFile *f2, string era1, string era2,
	     const char *dir, const char *base, const char *sfx,
	     const char *ytitle, double ymin, double ymax,
	     double rmin, double rmax);

void compareRun2Run3JEC(string era1 = "2018ABCD", string era2 = "2024_nib",
			string etabin = "eta00-13") {

  setTDRStyle();
  gSystem->mkdir("pdf/compareRun2Run3JEC", kTRUE);

  TFile *f1 = new TFile(Form("rootfiles/jecdata%s.root",era1.c_str()),"READ");
  TFile *f2 = new TFile(Form("rootfiles/jecdata%s.root",era2.c_str()),"READ");
  if (!f1 || f1->IsZombie()) {
    cout << "Could not open jecdata" << era1 << ".root" << endl << flush;
    return;
  }
  if (!f2 || f2->IsZombie()) {
    cout << "Could not open jecdata" << era2 << ".root" << endl << flush;
    return;
  }
  cout << "Comparing " << era1 << " (solid) to " << era2 << " (open)"
       << ", " << etabin << endl << flush;

  extraText = "Private work";

  // NB: keep these as strings, Form() reuses a rotating buffer
  string sd = Form("data/%s",etabin.c_str());  const char *cd = sd.c_str();
  string sm = Form("mc/%s",etabin.c_str());    const char *cm = sm.c_str();
  string sr = Form("ratio/%s",etabin.c_str()); const char *cr = sr.c_str();

  // JES applied to data, and to MC
  drawSet(f1, f2, era1, era2, cd, "pjes", "_data",
	  "JES applied to data", 0.83, 1.07, 0.95, 1.05);
  drawSet(f1, f2, era1, era2, cm, "pjes", "_mc",
	  "JES applied to MC",   0.83, 1.07, 0.95, 1.05);
  // Residual applied to data (the MC one is unity by construction)
  drawSet(f1, f2, era1, era2, cd, "pres", "_data",
	  "Residual applied to data", 0.95, 1.10, 0.95, 1.05);
  // Data over MC, the quantity the global fit sees
  drawSet(f1, f2, era1, era2, cr, "pjes", "_ratio",
	  "JES data/MC", 0.95, 1.10, 0.95, 1.05);
  drawSet(f1, f2, era1, era2, cr, "pres", "_ratio",
	  "Residual data/MC", 0.95, 1.10, 0.95, 1.05);

  cout << "Wrote pdf/compareRun2Run3JEC/*.pdf" << endl << flush;
} // compareRun2Run3JEC


// One canvas: all channels for both epochs on top, epoch ratio at the bottom
void drawSet(TFile *f1, TFile *f2, string era1, string era2,
	     const char *dir, const char *base, const char *sfx,
	     const char *ytitle, double ymin, double ymax,
	     double rmin, double rmax) {

  string scn = Form("%s%s",base,sfx); const char *cn = scn.c_str();

  TH1D *hup = new TH1D(Form("hup_%s",cn),
		       Form(";p_{T} (GeV);%s",ytitle), 100, 15, 3500);
  hup->SetMinimum(ymin);
  hup->SetMaximum(ymax);
  hup->GetXaxis()->SetMoreLogLabels();
  hup->GetXaxis()->SetNoExponent();

  string sera1(era1), sera2(era2);
  if (era1=="2018ABCD") sera1 = "2018";
  if (era2=="2024_nib") sera2 = "2024";
  if (era1=="2026C_w87") sera1 = "w87";
  if (era2=="2026C_w90") sera2 = "w90";
  
  TH1D *hdw = new TH1D(Form("hdw_%s",cn),
		       Form(";p_{T} (GeV);%s / %s",
			    sera1.c_str(),sera2.c_str()), 100, 15, 3500);
  hdw->SetMinimum(rmin);
  hdw->SetMaximum(rmax);
  hdw->GetXaxis()->SetMoreLogLabels();
  hdw->GetXaxis()->SetNoExponent();

  TCanvas *c = tdrDiCanvas(Form("c_%s",cn), hup, hdw, 0, 11);

  c->cd(1);
  gPad->SetLogx();
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->DrawLine(15, 1, 3500, 1);

  TLegend *leg = tdrLeg(0.55, 0.90-0.045*nchan, 0.90, 0.90);//, 0.040);
  leg->SetHeader(Form("#bf{%s} solid, #bf{%s} open",
		      sera1.c_str(), sera2.c_str()));

  int ndrawn(0);
  for (int i = 0; i != nchan; ++i) {

    const chan_t &ch = chans[i];
    string scname = Form("%s%s%s",base,ch.id,sfx);
    const char *cname = scname.c_str();

    TH1D *h1 = getHist(f1, dir, cname, era1.c_str());
    TH1D *h2 = getHist(f2, dir, cname, era2.c_str());
    if (!h1 && !h2) continue;

    c->cd(1);
    if (h1) tdrDraw(h1, "Pz", ch.mfull, ch.color, kSolid, -1, kNone);
    if (h2) tdrDraw(h2, "Pz", ch.mopen, ch.color, kSolid, -1, kNone);
    if (h1 && h2) leg->AddEntry(h1, ch.label, "PL");
    else if (h1)  leg->AddEntry(h1, Form("%s (%s only)",
					 ch.label, era1.c_str()), "PL");
    else          leg->AddEntry(h2, Form("%s (%s only)",
					 ch.label, era2.c_str()), "PL");
    ++ndrawn;

    // Bottom panel: epoch over epoch, point by point
    if (h1 && h2) {
      TH1D *hr = ratioOfHists(h1, h2, Form("hr_%s",cname));
      if (hr) {
	c->cd(2);
	gPad->SetLogx();
	tdrDraw(hr, "Pz", ch.mfull, ch.color, kSolid, -1, kNone);
      }
    }
  } // for i

  if (ndrawn==0) {
    cout << "  " << dir << "/" << cn << ": nothing found, skipping canvas"
	 << endl << flush;
    return;
  }

  c->cd(2);
  l->DrawLine(15, 1, 3500, 1);

  c->cd(1);
  gPad->RedrawAxis();
  c->cd(2);
  gPad->RedrawAxis();

  c->SaveAs(Form("pdf/compareRun2Run3JEC/compareRun2Run3JEC_%s_%s_vs_%s.pdf",
		 cn, era1.c_str(), era2.c_str()));
} // drawSet


// Fetch an object as a histogram, whether it was stored as a TProfile or not
TH1D *getHist(TFile *f, const char *dir, const char *name, const char *tag) {

  if (!f) return 0;
  string s = Form("%s/%s",dir,name);
  if (debug) cout << "Fetching " << f->GetName() << ":" << s << endl;
  TObject *o = f->Get(s.c_str());
  if (!o) {
    cout << "  missing " << f->GetName() << ":" << dir << "/" << name
	 << endl << flush;
    return 0;
  }

  if (o->InheritsFrom("TProfile"))
    return ((TProfile*)o)->ProjectionX(Form("%s_%s_px", name, tag));
  if (o->InheritsFrom("TH1D"))
    return (TH1D*)((TH1D*)o)->Clone(Form("%s_%s_cl", name, tag));

  cout << "  " << dir << "/" << name << " is a " << o->ClassName()
       << ", not drawn" << endl << flush;

  return 0;
} // getHist


// Point-by-point ratio, only if the binning matches
TH1D *ratioOfHists(TH1D *h1, TH1D *h2, const char *name) {

  if (!h1 || !h2) return 0;
  if (h1->GetNbinsX()!=h2->GetNbinsX() ||
      h1->GetXaxis()->GetXmin()!=h2->GetXaxis()->GetXmin() ||
      h1->GetXaxis()->GetXmax()!=h2->GetXaxis()->GetXmax()) {
    cout << "  binning mismatch for " << name << " ("
	 << h1->GetNbinsX() << " vs " << h2->GetNbinsX()
	 << " bins), no epoch ratio drawn" << endl << flush;
    return 0;
  }

  TH1D *h = (TH1D*)h1->Clone(name);
  h->Divide(h1, h2);

  return h;
} // ratioOfHists
