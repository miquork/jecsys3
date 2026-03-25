// L3Res time dependence with inclusive jets
// Use Incjet/h2pteta_lumi
// Projection sanity check with Incjet/p2pteta vs hpt13
#include "TFile.h"
#include "TH2D.h"

#include <fstream>
#include <cstdio>

#include "../tdrstyle_mod22.C"

string triggerString(double pt) {
  // Inclusive jet thresholds, mi[]
  //if (pt>686.) return "HLT_PFJet550";
  if (pt>686.) return "HLT_PFJet500";
  if (pt>548.) return "HLT_PFJet450";
  if (pt>468.) return "HLT_PFJet400";
  if (pt>395.) return "HLT_PFJet320";
  if (pt>330.) return "HLT_PFJet260";
  if (pt>272.) return "HLT_PFJet200";
  if (pt>196.) return "HLT_PFJet140";
  //if (pt>196.) return "HLT_PFJet110";
  if (pt>114.) return "HLT_PFJet80";
  if (pt>84.)  return "HLT_PFJet60";
  if (pt>49.)  return "HLT_PFJet40";
  return "HLT_PFJet40"; // "HLT_ZeroBias"
}

map<string, map<string, double> > _trgLumCache;
double triggerLumi(string run, string trg) {
  if (_trgLumCache[run][trg]!=0) return _trgLumCache[run][trg];

  const char *cr = run.c_str();
  string sf = Form("rootfiles/Prompt2024/lumiText/lumi_HLT_CombinedJSONAug01_%s_PFJet.csv",cr);
  cout << "Opening " << sf << endl << flush;
  ifstream f(sf.c_str());
  assert(f.is_open());

  double lumi(0);
  float rec(0), del(0); // in /fb
  int ver,nfill,nrun,ncms,nruntot(0);
  string s;
  string sRegExp = Form("#%s_v%%d,%%d,%%d,%%d,%%f,%%f",trg.c_str());
  const char *csr = sRegExp.c_str();
  cout << "Reading regExp " << csr << endl << flush;
  while (f >> s) {
    const char *cs = s.c_str();
    if (sscanf(cs,csr,&ver,&nfill,&nrun,&ncms,&del,&rec)==6) {
      cout << "Parsing string " << cs << endl << flush;
      cout << trg << "_v" << ver<<": " << rec << endl << flush;
      lumi += rec;
      nruntot += nrun;
    }
  } // while f

  cout << trg << ": " << lumi << " ("<<nruntot<<" runs)" << endl << flush;
  if (lumi==0) lumi = -1;
  _trgLumCache[run][trg] = lumi;

  return lumi;
}

void normalizeHistByLumi(TH1D *h, string run) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double pt = h->GetBinCenter(i);
    string trg = triggerString(pt);
    double lumi = triggerLumi(run, trg);
    
    h->SetBinContent(i, h->GetBinContent(i)/lumi);
    h->SetBinError(i, h->GetBinError(i)/lumi);
  }

} // normalizeHistByLumi

TH1D* timeDeps(string run="2024BCD") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TFile *f = new TFile(Form("rootfiles/Prompt2024/v86_2024/jmenano_data_cmb_%s_JME_v86_2024.root",cr),"READ");
  if (!f || f->IsZombie()) 
    f = new TFile(Form("rootfiles/Prompt2024/v83_2024/jmenano_data_cmb_%s_JME_v83_2024.root",cr),"READ");
  assert(f && !f->IsZombie());

  TH1D *hpt13 = (TH1D*)f->Get("Incjet/hpt13"); assert(hpt13);
  TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta"); assert(h2pteta);
  TH2D *h2pteta_lum = (TH2D*)f->Get("Incjet/h2pteta_lumi"); assert(h2pteta_lum);

  int i1 = h2pteta->GetXaxis()->FindBin(-1.305+1e-2);
  int i2 = h2pteta->GetXaxis()->FindBin(+1.305-1e-2);
  TH1D *hpt13p = h2pteta->ProjectionY(Form("hpt13p_%s",cr),i1,i2);
  hpt13->Scale(1000./14.674,"width");
  hpt13p->Scale(1000.,"width");
  TH1D *hpt13l = h2pteta_lum->ProjectionY(Form("hpt13l_%s",cr),i1,i2);
  hpt13l->Scale(1000./91,"width");
  
  //TH1D *h = tdrHist("h","N_{jet} / TeV, Xsec (fb / TeV)",0.5,2e17,
  TH1D *h = tdrHist(Form("h_%s",cr),"Cross section (fb / TeV)",0.2,2e17,
		    "p_{T,jet} (GeV)",15,4500.);
  TH1D *hd = tdrHist(Form("hd_%s",cr),"Ratio",0.95,1.05,
		     "p_{T,jet} (GeV)",15,4500.);
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cr),h,hd,8,11);

  normalizeHistByLumi(hpt13p,cr);
  
  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(hpt13,"Pz",kOpenCircle);
  tdrDraw(hpt13p,"Pz",kFullDiamond);
  tdrDraw(hpt13l,"Pz",kOpenDiamond);

  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();
  TH1D *hpt13r = (TH1D*)hpt13p->Clone("hpt13r");
  hpt13r->Divide(hpt13);
  tdrDraw(hpt13r,"Pz",kFullDiamond);

  TH1D *hpt13lr = (TH1D*)hpt13l->Clone("hpt13lr");
  hpt13lr->Divide(hpt13p);
  hpt13lr->Scale(1.04);//0.22);
  tdrDraw(hpt13lr,"Pz",kOpenDiamond);
  
  gPad->RedrawAxis();

  c1->Close();
  return hpt13p;
} // TH1D* timeDeps

void timeDep() {

  TH1D *h24bcd = timeDeps("2024BCD");
  TH1D *h24e   = timeDeps("2024E");
  TH1D *h24f   = timeDeps("2024F");

  TH1D *h24 = (TH1D*)h24bcd->Clone("h24"); h24->Reset();
  double lum24bcd = _trgLumCache["2024BCD"]["HLT_PFJet500"];
  double lum24e = _trgLumCache["2024E"]["HLT_PFJet500"];
  double lum24f = _trgLumCache["2024F"]["HLT_PFJet500"];
  double lumsum = lum24bcd + lum24e + lum24f;
  h24->Add(h24bcd, lum24bcd);
  h24->Add(h24e,   lum24e);
  h24->Add(h24f,   lum24f);
  h24->Scale(1./lumsum);

  TH1D *h = tdrHist("h","Cross section (fb / TeV)",0.2,2e17,
		    "p_{T,jet} (GeV)",15,4500.);
  TH1D *hd = tdrHist("hd","Ratio to 2024",0.75,1.25,
		     "p_{T,jet} (GeV)",15,4500.);
  lumi_136TeV = Form("2024, %1.1f fb^{-1}",lumsum);
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(h24,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(h24bcd,"Pz",kFullDiamond,kRed);
  tdrDraw(h24e,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(h24f,"Pz",kFullDiamond,kGreen+2);

  TLegend *leg = tdrLeg(0.60,0.85-4*0.05,0.85,0.85);
  leg->AddEntry(h24,"2024","L");
  leg->AddEntry(h24bcd,Form("BCD, %1.1f fb^{-1}",lum24bcd),"PLE");
  leg->AddEntry(h24e,Form("E, %1.1f fb^{-1}",lum24e),"PLE");
  leg->AddEntry(h24f,Form("F, %1.1f fb^{-1}",lum24f),"PLE");

  TLatex *tex = new TLatex(0.35,0.80);
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex("|#eta| < 1.305");
  
  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();

  TH1D *hr24bcd = (TH1D*)h24bcd->Clone("hr24bcd");
  TH1D *hr24e = (TH1D*)h24e->Clone("hr24e");
  TH1D *hr24f = (TH1D*)h24f->Clone("hr24f");
  TH1D *hr24 = (TH1D*)h24->Clone("hr24");
  hr24bcd->Divide(h24);
  hr24e->Divide(h24);
  hr24f->Divide(h24);
  hr24->Divide(h24);

  tdrDraw(hr24,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hr24bcd,"Pz",kFullDiamond,kRed);
  tdrDraw(hr24e,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(hr24f,"Pz",kFullDiamond,kGreen+2);
    
  gPad->RedrawAxis();

  c1->SaveAs("pdf/timeDep/timeDep.pdf");
}
