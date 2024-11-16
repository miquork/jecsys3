// Purpose: Draw pairs (or multiplets) of time scans made by drawTimeStability.C
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include <set>

#include "../tdrstyle_mod22.C"

// Clean out bins from hd that are empty for ha or hb
void cleanEmpty(TH1D *hd, TH1D *ha, TH1D *hb);
// Draw Gamma scale vs Zmm scale
void drawGamVsZmm();
void drawPFcomp();

// Main call for drawing various pairs
void drawTimeStabilityPairs() {

  //drawGamVsZmm();
  drawPFcomp();
}

void drawGamVsZmm() {
  
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
  
  TH1D *ha = (TH1D*)f->Get("pr50m"); assert(ha);
  TH1D *hb = (TH1D*)f->Get("mpf_run_zpt50"); assert(hb);
  TH1D *hd = (TH1D*)ha->Clone("hd");
  hd->Add(hb,-1);
  cleanEmpty(hd,ha,hb);

  // Extra pairs
  TH1D *ha2 = (TH1D*)f->Get("pr110m"); assert(ha2);
  TH1D *hb2 = (TH1D*)f->Get("mpf_run_zpt110"); assert(hb2);
  TH1D *hd2 = (TH1D*)ha2->Clone("hd2");
  hd2->Add(hb2,-1);
  cleanEmpty(hd2,ha2,hb2);
  
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

  tdrDraw(hd2,"Pz",kOpenSquare,kMagenta-9,kSolid,-1,kNone,0,0.6);
  tdrDraw(hd,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.6);

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
  }
  TH1D *hnib = new TH1D("hnib","",nib,&vx[0]);
  for (int i = 0; i != nib; ++i) {
    hnib->SetBinContent(i+1, vy[i]);
    hnib->SetBinError(i+1, vye[i]);
  }
  tdrDraw(hnib,"E2",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
  hnib->SetFillColorAlpha(kCyan+1,0.5);
  
  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/drawTimeStability/drawTimeStabilityPairs_GamVsZmm.pdf");
} // drawTimeStabilityPairs()



void drawPFcomp() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  double xmin(0), xmax(175), ymin(-3.6), ymax(+5.4);//, ymind(-1.5), ymaxd(+1.9);
  TH1D *h = tdrHist("hpf_u","PFJES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  //TH1D *h_d = tdrHist("hpf_d","#gamma+j - Z+j (%)",ymind,ymaxd,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  lumi_136TeV = "Run3, 2022-24";
  extraText = "Private";
  //TCanvas *c1 = tdrDiCanvas("c1",h,h_d,8,11);
  TCanvas *c1 = tdrCanvas("c1pf",h,8,11,kSquare);

  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);

  /*
  TH1D *hr = (TH1D*)f->Get("pr50m_jes"); assert(hr);
  TH1D *hc = (TH1D*)f->Get("pr50chf");   assert(hc);
  TH1D *hn = (TH1D*)f->Get("pr50nhf");   assert(hn);
  TH1D *he = (TH1D*)f->Get("pr50nef");   assert(he);
  */
  TH1D *hr = (TH1D*)f->Get("pr110m_jes"); assert(hr);
  TH1D *hc = (TH1D*)f->Get("pr110chf");   assert(hc);
  TH1D *hn = (TH1D*)f->Get("pr110nhf");   assert(hn);
  TH1D *he = (TH1D*)f->Get("pr110nef");   assert(he);
  
  // Approximate PF fractions
  //hc->Scale(0.65);
  //he->Scale(0.25);
  //hn->Scale(0.10);

  // Calculate JES and scale with per-fib JES
  TH1D *hjesr = (TH1D*)hr->Clone("hjesr");
  TH1D *hjesc = (TH1D*)hc->Clone("hjesc");
  TH1D *hjesn = (TH1D*)hn->Clone("hjesn");
  TH1D *hjese = (TH1D*)he->Clone("hjese");
  double kr(1.02), fc(0.65), fe(0.25), fn(0.10);
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    double s = hr->GetBinContent(i); // s = (jes-1)*100
    double jes = s*kr*0.01+1;
    double jesr = (jes-1)*100.;
    hjesr->SetBinContent(i, jesr);
    double c = hc->GetBinContent(i);
    double jesc = ((1+c*0.01)*fc*jes-fc)*100;
    hjesc->SetBinContent(i, jesc);
    double n = hn->GetBinContent(i);
    double jesn = ((1+n*0.01)*fn*jes-fn)*100;
    hjesn->SetBinContent(i, jesn);
    double e = he->GetBinContent(i);
    double jese = ((1+e*0.01)*fe*jes-fe)*100;
    hjese->SetBinContent(i, jese);
  }
  
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
  
  tdrDraw(hjesc,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.6);
  tdrDraw(hjese,"Pz",kOpenSquare,kBlue,kSolid,-1,kNone,0,0.6);
  tdrDraw(hjesn,"Pz",kOpenCircle,kGreen+2,kSolid,-1,kNone,0,0.6);
  tdrDraw(hjesr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);

  //TLegend *leg = tdrLeg(0.65,0.89-0.05*4,0.90,0.89);
  TLegend *leg = tdrLeg(0.60,0.89-0.05*4,0.85,0.89);
  leg->SetFillStyle(1001);
  /*
  leg->AddEntry(hjesr,"#gamma+jet MPF 50","PLE");
  leg->AddEntry(hjesc,"#gamma+jet CHF 50","PLE");
  leg->AddEntry(hjesn,"#gamma+jet NHF 50","PLE");
  leg->AddEntry(hjese,"#gamma+jet NEF 50","PLE");
  */
  leg->AddEntry(hjesr,"#gamma+jet MPF 110","PLE");
  leg->AddEntry(hjesc,"#gamma+jet CHF 110","PLE");
  leg->AddEntry(hjesn,"#gamma+jet NHF 110","PLE");
  leg->AddEntry(hjese,"#gamma+jet NEF 110","PLE");
  
  gPad->RedrawAxis();
  gPad->Update();
  
  c1->SaveAs("pdf/drawTimeStability/drawTimeStabilityPairs_PFcomp.pdf");
} // drawTimeStabilityPairs()



void cleanEmpty(TH1D *hd, TH1D *ha, TH1D *hb) {

  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    if ((ha->GetBinContent(i)==0 && ha->GetBinError(i)==0) ||
	(hb->GetBinContent(i)==0 && hb->GetBinError(i)==0)) {
      hd->SetBinContent(i,0);
      hd->SetBinError(i,0);
    }
  }
}
