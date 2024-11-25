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
// Find scale from hscales
double getScale(TH1D *hscales, string name);
// Draw Gamma scale vs Zmm scale
void drawGamVsZmm();
void drawGamVsGam();
void drawPFcomp(string ref);

// Main call for drawing various pairs
void drawTimeStabilityPairs() {

  //drawGamVsZmm();

  drawPFcomp("pr50");
  drawPFcomp("pr110");
  drawPFcomp("pr230");
  /*
  drawPFcomp("zpt30");
  drawPFcomp("zpt50");
  drawPFcomp("zpt110");
  */
  //drawGamVsGam();
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


void drawPFcomp(string ref) {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = ref.c_str();
  TString tr(cr);  

  TFile *f = new TFile("rootfiles/drawTimeStability.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  //double xmin(0), xmax(175), ymin(-3.6), ymax(+5.4);//, ymind(-1.5), ymaxd(+1.9);
  double xmin(0), xmax(175), ymin(-3.7), ymax(+6.3);//, ymind(-1.5), ymaxd(+1.9);
  TH1D *h = tdrHist(Form("hpf_%s",cr),"PFJES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
  lumi_136TeV = "Run3, 2022-24";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1pf_%s",cr),h,8,11,kSquare);

  TH1D *hbreaks = (TH1D*)f->Get("hbreaks"); assert(hbreaks);
  TH1D *hscales = (TH1D*)f->Get("hscales"); assert(hscales);

  /*
  TH1D *hr = (TH1D*)f->Get("pr50m_jes"); assert(hr);
  TH1D *hc = (TH1D*)f->Get("pr50chf");   assert(hc);
  TH1D *hn = (TH1D*)f->Get("pr50nhf");   assert(hn);
  TH1D *he = (TH1D*)f->Get("pr50nef");   assert(he);
  */
  /*
  TH1D *hr = (TH1D*)f->Get("pr110m_jes"); assert(hr);
  TH1D *hc = (TH1D*)f->Get("pr110chf");   assert(hc);
  TH1D *hn = (TH1D*)f->Get("pr110nhf");   assert(hn);
  TH1D *he = (TH1D*)f->Get("pr110nef");   assert(he);
  */
  /*
  TH1D *hr = (TH1D*)f->Get(Form("%sm_jes",cr)); assert(hr);
  TH1D *hc = (TH1D*)f->Get(Form("%schf",cr));  assert(hc);
  TH1D *hn = (TH1D*)f->Get(Form("%snhf",cr));  assert(hn);
  TH1D *he = (TH1D*)f->Get(Form("%snef",cr));  assert(he);
  */

  // Format string for Z+jet
  TH1D *hr(0), *hc(0), *hn(0), *he(0);
  double kr(1), fc(0.65), fe(0.25), fn(0.10);
  if (tr.Contains("pr")) {
    hr = (TH1D*)f->Get(Form("%sm_jes",cr)); assert(hr);
    hc = (TH1D*)f->Get(Form("%schf",cr));  assert(hc);
    hn = (TH1D*)f->Get(Form("%snhf",cr));  assert(hn);
    he = (TH1D*)f->Get(Form("%snef",cr));  assert(he);

    kr = getScale(hscales,Form("%sm_jes",cr));
    fc = getScale(hscales,Form("%schf",cr));
    fe = getScale(hscales,Form("%snef",cr));
    fn = getScale(hscales,Form("%snhf",cr));
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

  /*
  double kr = getScale(hscales,Form("%sm_jes",cr));//hr->GetName());
  double fc = getScale(hscales,Form("%schf",cr));//hc->GetName());
  double fe = getScale(hscales,Form("%snef",cr));//he->GetName());
  double fn = getScale(hscales,Form("%snhf",cr));//hn->GetName());
  */
  // Calculate JES and scale with per-fib JES
  TH1D *hjesr = (TH1D*)hr->Clone(Form("hjesr_%s",cr));
  TH1D *hjesc = (TH1D*)hc->Clone(Form("hjesc_%s",cr));
  TH1D *hjesn = (TH1D*)hn->Clone(Form("hjesn_%s",cr));
  TH1D *hjese = (TH1D*)he->Clone(Form("hjese_%s",cr));
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    double s = hr->GetBinContent(i); // s = (jes-1)*100
    double es = hr->GetBinError(i);
    double jes = s*kr*0.01+1;
    double jesr = (jes-1)*100.;
    double ejesr = es*kr;
    hjesr->SetBinContent(i, jesr);
    hjesr->SetBinError(i, ejesr);
    double c = hc->GetBinContent(i);
    double ec = hc->GetBinError(i);
    double jesc = ((1+c*0.01)*fc*jes-fc)*100;
    double ejesc = ec*fc*jes;
    hjesc->SetBinContent(i, jesc);
    hjesc->SetBinError(i, ejesc);
    double n = hn->GetBinContent(i);
    double en = hn->GetBinError(i);
    double jesn = ((1+n*0.01)*fn*jes-fn)*100;
    double ejesn = en*fn*jes;
    hjesn->SetBinContent(i, jesn);
    hjesn->SetBinError(i, ejesn);
    double e = he->GetBinContent(i);
    double ee = he->GetBinError(i);
    double jese = ((1+e*0.01)*fe*jes-fe)*100;
    double ejese = ee*fe*jes;
    hjese->SetBinContent(i, jese);
    hjese->SetBinError(i, ejese);
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
  TLegend *leg = tdrLeg(0.58,0.89-0.05*4,0.83,0.89);
  leg->SetFillStyle(1001);
  if (ref=="pr50" || ref=="pr110" || ref=="pr230") {
    int ipt; sscanf(cr,"pr%d",&ipt);
    leg->AddEntry(hjesr,Form("#gamma+jet MPF %d",ipt),"PLE");
    leg->AddEntry(hjesc,Form("#gamma+jet CHF %d",ipt),"PLE");
    leg->AddEntry(hjesn,Form("#gamma+jet NHF %d",ipt),"PLE");
    leg->AddEntry(hjese,Form("#gamma+jet NEF %d",ipt),"PLE");
  }
  if (ref=="zpt30" || ref=="zpt50" || ref=="zpt110") {
    int ipt; sscanf(cr,"zpt%d",&ipt);
    leg->AddEntry(hjesr,Form("Z_{#mu#mu}+jet MPF %d",ipt),"PLE");
    leg->AddEntry(hjesc,Form("Z_{#mu#mu}+jet CHF %d",ipt),"PLE");
    leg->AddEntry(hjesn,Form("Z_{#mu#mu}+jet NHF %d",ipt),"PLE");
    leg->AddEntry(hjese,Form("Z_{#mu#mu}+jet NEF %d",ipt),"PLE");
  }
    
  gPad->RedrawAxis();
  gPad->Update();
  
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
}

double getScale(TH1D *hscales, string name) {

  for (int i = 1; i != hscales->GetNbinsX()+1; ++i) {
    if (string(hscales->GetXaxis()->GetBinLabel(i))==name) return hscales->GetBinContent(i);
  } // for i
  cout << "getScale not found for " << name << endl << flush;
  return 1;
} // getScale
