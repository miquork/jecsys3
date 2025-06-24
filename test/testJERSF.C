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

void testJERSF(string filename) {

  setTDRStyle();

  cout << "Open " << filename << endl << flush;
  JetCorrectorParameters *p = new JetCorrectorParameters(filename.c_str());
  vector<JetCorrectorParameters> v;
  v.push_back(*p);
  FactorizedJetCorrector *jer = new FactorizedJetCorrector(v);

  if (debug) {
    jer->setJetEta(0.);
    jer->setJetPt(100.);
    //jer->setRho(20.85); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=100,rho=20.85)="<<jer->getCorrection()<<endl;
    jer->setJetEta(0.);
    jer->setJetPt(100.);
    //jer->setRho(5.); // UL18 Z+jet vs pT at [60,130] GeV, |eta|<1.3
    cout << "JER(eta=0,pt=100,rho=5)="<<jer->getCorrection()<<endl;
  }

  const int ndiv = 1;//5;
  const int nbins = 60*ndiv;
  const double deta = TMath::TwoPi()/72./ndiv;
  const double maxeta = nbins*deta; // 5.236
  const double rho = 20.85;
  TH1D *h10 = new TH1D("h10",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h100 = new TH1D("h100",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h1000 = new TH1D("h1000",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *hxmax = new TH1D("hxmax",";#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h10m = new TH1D("h10m",";-#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h100m = new TH1D("h100m",";-#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *h1000m = new TH1D("h1000m",";-#eta_{jet};JER SF",nbins,0.,maxeta);
  TH1D *hxmaxm = new TH1D("hxmaxm",";-#eta_{jet};JER SF",nbins,0.,maxeta);
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
  
  double miny = 0.9; // 0.85
  double maxy = 2.8; // 1.75
  TH1D *h = tdrHist("h","JER SF",miny,maxy,"|#eta_{jet}|",0,5.2);
  //lumi_13TeV = "UL2018, 59.9 fb^{-1}";
  //TCanvas *c1 = tdrCanvas("c1",h,4,11,kSquare);
  lumi_136TeV = name_short; //"Summer23, Prompt23";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(0,1,5.2,1);
  
  tdrDraw(h10,"HPz",kOpenDiamond,kMagenta+2,kSolid,-1,kNone);
  tdrDraw(h100,"HPz",kFullCircle,kGreen+2,kSolid,-1,kNone);
  tdrDraw(h1000,"H][Pz",kFullDiamond,kRed,kSolid,-1,kNone);
  tdrDraw(hxmax,"HPz",kOpenDiamond,kBlack,kSolid,-1,kNone);

  tdrDraw(h10m,"HPz",kOpenDiamond,kMagenta+2-9,kDotted,-1,kNone);
  tdrDraw(h100m,"HPz",kOpenCircle,kGreen+1,kDotted,-1,kNone);
  tdrDraw(h1000m,"H][Pz",kFullDiamond,kRed-9,kDotted,-1,kNone);
  tdrDraw(hxmaxm,"HPz",kOpenDiamond,kGray+2,kDotted,-1,kNone);
  h10m->SetMarkerSize(0.6);
  h100m->SetMarkerSize(0.6);
  h1000m->SetMarkerSize(0.6);
  hxmaxm->SetMarkerSize(0.6);

  TLegend *leg = tdrLeg(0.62,0.90-0.05*7,0.87,0.90);
  leg->AddEntry(hxmax,"p_{T,max}","PLE");
  leg->AddEntry(h1000,"1000 GeV","PLE");
  leg->AddEntry(h100,"100 GeV","PLE");
  leg->AddEntry(h10,"10 GeV","PLE");
  leg->AddEntry(hxmaxm,"p_{T,max}, -|#eta|","PLE");
  leg->AddEntry(h1000m,"1000 GeV, -|#eta|","PLE");
  leg->AddEntry(h100m,"100 GeV, -|#eta|","PLE");
  leg->AddEntry(h10m,"10 GeV, -|#eta|","PLE");
  
  c1->SaveAs(Form("pdf/JERSF/testJERSF_%s.pdf",name.Data()));
} // testJERSF
