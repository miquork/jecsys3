// Purpose: test JER SF text file produced by JERSF.C
#include "../CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "../CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "TMath.h"

#include "../tdrstyle_mod22.C"

#include <iostream>
#include <vector>

using namespace std;

const bool debug = true;

void testMCJER(string filename, string mcname="") {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/testMCJER");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/testMCJER");

  setTDRStyle();

  cout << "Open " << filename << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(filename.c_str());
  vector<JetCorrectorParameters> v;
  v.push_back(*p);
  FactorizedJetCorrector *jer = new FactorizedJetCorrector(v);

  if (debug) {
    jer->setJetEta(0.);
    jer->setJetPt(100.);
    jer->setRho(20.85); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=100,rho=20.85)="<<jer->getCorrection()<<endl;
    jer->setJetEta(0.);
    jer->setJetPt(100.);
    jer->setRho(5.); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=100,rho=5)="<<jer->getCorrection()<<endl;
  }

  const int ndiv = 1;//5;
  const int nbins = 60*ndiv;
  const double deta = TMath::TwoPi()/72./ndiv;
  const double maxeta = nbins*deta; // 5.236
  const double rho = 20.85;
  TH1D *h10 = new TH1D("h10",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h15 = new TH1D("h15",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h30 = new TH1D("h30",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h100 = new TH1D("h100",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h1000 = new TH1D("h1000",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *hxmax = new TH1D("hxmax",";#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h10m = new TH1D("h10m",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h15m = new TH1D("h15m",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h30m = new TH1D("h30m",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h100m = new TH1D("h100m",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *h1000m = new TH1D("h1000m",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  TH1D *hxmaxm = new TH1D("hxmaxm",";-#eta_{jet};MC JER",nbins,0.,maxeta);
  for (int i = 0; i != h100->GetNbinsX()+1; ++i) {

    double eta = h100->GetBinCenter(i);
    double etamin = h100->GetBinLowEdge(i);

    jer->setJetEta(eta);
    jer->setJetPt(10.);
    jer->setRho(rho);
    h10->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(-eta);
    jer->setJetPt(10.);
    jer->setRho(rho);
    h10m->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(eta);
    jer->setJetPt(15.);
    jer->setRho(rho);
    h15->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(-eta);
    jer->setJetPt(15.);
    jer->setRho(rho);
    h15m->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(eta);
    jer->setJetPt(30.);
    jer->setRho(rho);
    h30->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(-eta);
    jer->setJetPt(30.);
    jer->setRho(rho);
    h30m->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(eta);
    jer->setJetPt(100.);
    jer->setRho(rho);
    h100->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(-eta);
    jer->setJetPt(100.);
    jer->setRho(rho);
    h100m->SetBinContent(i, jer->getCorrection());

    if (1000.*cosh(etamin)<6800.) {
      jer->setJetEta(eta);
      jer->setJetPt(1000.);
      jer->setRho(rho);
      h1000->SetBinContent(i, jer->getCorrection());

      jer->setJetEta(-eta);
      jer->setJetPt(1000.);
      jer->setRho(rho);
      h1000m->SetBinContent(i, jer->getCorrection());
    }
      
    double ptmax = min(4000.,6800./cosh(etamin));
    jer->setJetEta(eta);
    jer->setJetPt(ptmax);
    jer->setRho(rho);
    hxmax->SetBinContent(i, jer->getCorrection());

    jer->setJetEta(-eta);
    jer->setJetPt(ptmax);
    jer->setRho(rho);
    hxmaxm->SetBinContent(i, jer->getCorrection());
  }

  TString name = filename.c_str();
  Ssiz_t lastSlashPos = name.Last('/'); // Remove path before file name
  if (lastSlashPos != kNPOS) name = name(lastSlashPos+1, name.Length());
  name.ReplaceAll(".txt",""); // Remove .txt ending
  TString name_short = name;
  name_short.ReplaceAll("_MC_SF_AK4PFPuppi","");
  
  double miny = 0;
  double maxy = 1.20-1e-4;
  TH1D *h = tdrHist("h","MC JER",miny,maxy,"|#eta_{jet}|",0,5.2);
  lumi_136TeV = mcname;//name_short;
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0,1,5.2,1);
  l->SetLineStyle(kDotted);
  l->SetLineColor(kGray+1);
  l->DrawLine(0,0.5,5.2,0.5);
  l->SetLineColor(kGray);
  l->DrawLine(0,0.1,5.2,0.1);
  l->DrawLine(1.305,0.,1.305,0.65);
  l->DrawLine(2.5,0.,2.5,1.0);
  l->DrawLine(2.964,0.,2.964,1.0);
  
  //tdrDraw(h10,"HPz",kOpenDiamond,kMagenta+2,kSolid,-1,kNone);
  tdrDraw(h15,"HPz",kOpenDiamond,kMagenta+1,kSolid,-1,kNone);
  tdrDraw(h30,"HPz",kOpenDiamond,kMagenta+3,kSolid,-1,kNone);
  tdrDraw(h100,"HPz",kFullCircle,kGreen+2,kSolid,-1,kNone);
  tdrDraw(h1000,"H][Pz",kFullDiamond,kRed,kSolid,-1,kNone);
  //tdrDraw(hxmax,"HPz",kOpenDiamond,kBlack,kSolid,-1,kNone);

  //tdrDraw(h10m,"HPz",kOpenDiamond,kMagenta+2-9,kDotted,-1,kNone);
  tdrDraw(h15m,"HPz",kOpenDiamond,kMagenta+1-9,kDotted,-1,kNone);
  tdrDraw(h30m,"HPz",kOpenDiamond,kMagenta+3-9,kDotted,-1,kNone);
  tdrDraw(h100m,"HPz",kOpenCircle,kGreen+1,kDotted,-1,kNone);
  tdrDraw(h1000m,"H][Pz",kFullDiamond,kRed-9,kDotted,-1,kNone);
  //tdrDraw(hxmaxm,"HPz",kOpenDiamond,kGray+2,kDotted,-1,kNone);
  //h10m->SetMarkerSize(0.6);
  h15m->SetMarkerSize(0.6);
  h30m->SetMarkerSize(0.6);
  h100m->SetMarkerSize(0.6);
  h1000m->SetMarkerSize(0.6);
  hxmaxm->SetMarkerSize(0.6);

  TLegend *leg = tdrLeg(0.24,0.78-0.04*5,0.49,0.78);
  leg->SetHeader("+|#eta|");
  leg->AddEntry(h15,"15 GeV","PLE");
  leg->AddEntry(h30,"30 GeV","PLE");
  leg->AddEntry(h100,"100 GeV","PLE");
  leg->AddEntry(h1000,"1000 GeV","PLE");

  double dx = -0.06;
  TLegend *legm = tdrLeg(0.24+dx,0.78-0.04*5,0.49+dx,0.78);
  legm->SetHeader("-|#eta|");
  legm->AddEntry(h15m," ","PLE");
  legm->AddEntry(h30m," ","PLE");
  legm->AddEntry(h100m," ","PLE");
  legm->AddEntry(h1000m," ","PLE");
  
  c1->SaveAs(Form("pdf/testMCJER/testMCJER_%s.pdf",name.Data()));
} // testMCJER
