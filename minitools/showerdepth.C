// Purpose: estimate shower depth in ECAL+HCAL, to turn HCAL scale per layer
//          into HCAL scale vs pT
// Assumptions: 1) ECAL 1 lambda, each HCAL layer (out of 13) 0.5 Lambda
//              2) Each nuclear interaction produces 5 secondaries (5-6)
//              3) Each secondary deposits MIP energy of 1 GeV
//              4) Each layer causes probability branching, so 2^14=16384 paths
// Run3 : Turn this into fraction of energy not in depth #4.
//        Still need to map from particle to jet
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"

#include "../tdrstyle_mod22.C"

#include <iostream>
using namespace std;

const int Nlayer = 15; // two f
double Elayer[Nlayer];
const double khcal = 1;//1./0.9;
//const double w[Nlayer+2] = 
//  {0.92,0.92, 0.79,0.87,0.90,0.89, 0.93,0.96,0.98,0.97, 0.99,0.99,1.00,1.01, 1,1,1}; // model of front layer radiation damage
//{1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1};
  //  {1,1, 0.76,0.81,0.82,0.79, 0.92,0.94,0.97,...} // slide13
// For calculating energy fraction not in depth4
const double w[Nlayer+2] = 
  {1,1, 1, 1,1,1, 1,1,1, 0,0,0,0, 1, 1,1,1}; // model of front layer radiation damage
void propagate(double pn, double E, int layer) {

  // stop once going past last layer
  if (layer>=Nlayer || E<=0) return;

  const int NsecMax = 5; // Typical number of secondaries in nuclear interaction
  const int NsecMin = 2; // Minimum number of secondaries
  const double Emin = 1; // Minimum energy per secondary
  const double EMIP = 1; // MIP energy loss/deposition per layer
  const double P = 1-exp(-0.5); // Probability of interaction per layer, ~40%
  const double fE = 0.27; // Fraction of energy going to EM branch through pi0s
  const double eoh = 1.8; // Efficiency of hadronic energy deposition over EM

  double Nsec(0), Esec(0), Emip(0), Edep(0);

  // Case 1: number of secondaries not limited, split energy evenly
  if (E>NsecMax*Emin) {
    Nsec = NsecMax;
    Esec = E/Nsec - EMIP; // nuclear interaction
    Emip = E - EMIP;      // MIP through
    Edep = (1-P)*EMIP/eoh + P*Nsec*((1-fE)*EMIP/eoh+fE*(Esec+EMIP));
  }
  // Case 2: number of secondaries limited, produce maximum at Emin
  else if (E>NsecMin*Emin) {
    Nsec = int(E/Emin);
    Esec = E/Nsec - EMIP; // nuclear interaction, Esec may be negative
    Emip = E - EMIP; // MIP through, Emip may be negative
    Edep = (1-P)*(EMIP+min(0.,Emip))/eoh
      + P*Nsec*((1-fE)*(EMIP+min(0.,Esec))/eoh+fE*(Esec+EMIP));
  }
  // Case 3: no more energy for secondary production, just deposit MIP
  else if (E>0) {
    Nsec = 0;
    Esec = 0;
    Emip = E - EMIP; // Emip may be negative
    Edep = (EMIP+min(0.,Emip))/eoh;
  }

  // Add deposited energy to this layer
  Elayer[layer] += pn*Edep;

  // Continue to next layer;
  propagate(pn*(1-P), Emip, layer+1);
  propagate(pn*Nsec*(1-fE)*P, Esec, layer+1);

  return;
}

TGraph* showerdepths(double E, int marker, int color,
		     double &fw, double &fhw, double &fe) {

  //double E = 100;
  for (int i = 0; i != Nlayer; ++i) Elayer[i]=0;
  propagate(1,E,0);

  double Esum(0), Esumw(0), Hsum(0), Hsumw(0);
  TGraph *g = (marker!=-1 ? new TGraph(Nlayer) : 0);
  for (int i = 0; i != Nlayer; ++i) {
    Esum += Elayer[i]/E;
    Hsum += (i<2 ? 0 : Elayer[i]/E);
    //Esumw += Elayer[i]/E * w[i];
    Esumw += Elayer[i]/E * (i<2 ? w[i] : w[i]*khcal);
    Hsumw += Elayer[i]/E * (i<2 ? 0 : w[i]*khcal);
    // Two first layers are effectively ECAL
    if (marker!=-1) g->SetPoint(i, i-1, Elayer[i]/E);
  }
  fw = Esumw / Esum;
  fhw = Hsumw / Hsum;
  fe = (Esum-Hsum) / Esum;
  cout << "E = " << E << ", Esum/E = " << Esum << ", fw = " << fw << endl;

  if (marker!=-1) {
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->Draw("SAMEPL");
  }
   
 return g;
}

void showerdepth() {

  setTDRStyle();

  TH1D *h1 = new TH1D("h",";Layer (#);Fraction of energy",15,-1.5,13.5);
  h1->SetMaximum(0.35);

  extraText = "Private";
  lumi_13TeV = "Run3 toy model";
  TCanvas *c1 = tdrCanvas("c1",h1,4,0,kSquare);

  double f10,f50,f100,f500, fwh, fe;
  TGraph *g10 = showerdepths(10, kFullSquare, kBlue, f10, fwh,fe);
  TGraph *g50 = showerdepths(50, kOpenSquare, kGreen+2, f50, fwh,fe);
  TGraph *g100 = showerdepths(100, kFullCircle, kOrange+2, f100, fwh,fe);
  TGraph *g500 = showerdepths(500, kOpenCircle, kRed, f500, fwh,fe);

  TLegend *leg = tdrLeg(0.40,0.70,0.60,0.90);
  leg->AddEntry(g10,Form("E=10 GeV (%1.2f, %1.2f)",f10,fwh),"P");
  leg->AddEntry(g50,Form("E=50 GeV (%1.2f, %1.2f)",f50,fwh),"P");
  leg->AddEntry(g100,Form("E=100 GeV (%1.2f, %1.2f)",f100,fwh),"P");
  leg->AddEntry(g500,Form("E=500 GeV (%1.2f, %1.2f)",f500,fwh),"P");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(0.5,0,0.5,0.35);
  l->SetLineStyle(kDotted);
  l->DrawLine(1.5,0,1.5,0.20); // depth #1 end (0.5 lambda, 1 layer)
  l->DrawLine(4.5,0,4.5,0.22); // depth #2 end (1.5 lambda, 3 layers)
  l->SetLineStyle(kSolid);
  l->DrawLine(7.5,0,7.5,0.22); // depth #3 end (1.5 lambda, 3 layers)
  l->DrawLine(11.5,0,11.5,0.22); // depth #4 end (2.0 lambda, 4 layers)

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.035);
  tex->SetNDC();
  tex->DrawLatex(0.165,0.85,"ECAL");
  tex->DrawLatex(0.265,0.85,"HCAL");
  
  tex->SetNDC(kFALSE);
  tex->DrawLatex(0.75,0.20,"depth #1");
  tex->DrawLatex(1.75,0.18,"depth #2");
  tex->DrawLatex(4.75,0.20,"depth #3");
  tex->DrawLatex(8.25,0.20,"depth #4");
  tex->SetNDC();


  double ve[] = {3, 3.4, 4, 5, 6, 7, 8, 10, 12, 15, 18, 22, 26, 32, 38, 46, 55, 66, 79, 95, 115, 138, 165, 198, 238, 286, 343, 412, 494, 593, 712, 854, 1025};
  const int ne = sizeof(ve)/sizeof(ve[0]);
  TH1D *h = new TH1D("h",";p_{T} (GeV);Energy fraction",ne-1,ve);
  TGraph *g = new TGraph(ne);
  TGraph *gh = new TGraph(ne);
  TGraph *ge = new TGraph(ne);
  for (int i = 0; i != ne; ++i) {
    double fw, fwh, fe;
    showerdepths(ve[i], -1,-1, fw, fwh, fe);
    h->SetBinContent(i+1, ve[i], fw);
    g->SetPoint(i, ve[i], fw);
    gh->SetPoint(i, ve[i], fwh);
    ge->SetPoint(i, ve[i], fe);
  }

  TH1D *h3 = new TH1D("h3",";p_{T} (GeV);ECAL fraction",997,3,1000);
  h3->SetMaximum(1.);
  h3->SetMinimum(0.);
  h3->GetXaxis()->SetMoreLogLabels();
  h3->GetXaxis()->SetNoExponent();

  TCanvas *c3 = tdrCanvas("c3",h3,4,0,kSquare);
  gPad->SetLogx();

  ge->SetMarkerStyle(kFullCircle);
  ge->Draw("SAMEP");


  //TH1D *h2 = new TH1D("h2",";p_{T} (GeV);Scaled hadron response",997,3,1000);
  TH1D *h2 = new TH1D("h2",";p_{T} (GeV);Energy fraction",997,3,3000);
  //h2->SetMaximum(1.15);
  //h2->SetMinimum(0.75);
  h2->SetMaximum(1.05);
  h2->SetMinimum(0.85);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();

  TCanvas *c2 = tdrCanvas("c2",h2,4,0,kSquare);
  gPad->SetLogx();

  h->Draw("SAMEHIST");

  //g->SetMarkerStyle(kFullCircle);
  g->SetMarkerStyle(kOpenCircle);
  g->Draw("SAMEP");
  //gh->SetMarkerStyle(kFullCircle);
  //gh->Draw("SAMEP");

  //TF1 *fIsoTrack = new TF1("fIsoTrack","[0]",40,60);
  //gh->Fit(fIsoTrack,"QRN");
  //fIsoTrack->SetLineWidth(3);
  //fIsoTrack->Draw("SAME");

  //TGraph *ghc = new TGraph(ne);
  //for (int i = 0; i != ne; ++i) {
  //double E = gh->GetX()[i];
  //ghc->SetPoint(i, E, gh->GetY()[i]/fIsoTrack->Eval(E));
  //}

  //ghc->SetMarkerStyle(kFullCircle);
  //ghc->SetMarkerColor(kRed);
  //ghc->Draw("SAMEP");

  //TF1 *fspr = new TF1("fspr","[2]*(1+[0]*pow(x,[1]))",3,1000);
  //fspr->SetParameters(1.1,-0.5,-0.3);
  //ghc->Fit(fspr,"QNR");
  //fspr->SetLineWidth(2);
  //fspr->Draw("SAME");

  //tex->SetTextColor(kRed);
  //tex->DrawLatex(0.35,0.20,Form("f_{SPR}(p_{T}) = %1.3f #times "
  //				"(1%+1.3f(p_{T}/GeV)^{%+1.3f})",
  //				fspr->GetParameter(2), fspr->GetParameter(0),
  //				fspr->GetParameter(1)));

  //TLegend *leg2 = tdrLeg(0.40,0.75,0.60,0.90);
  //leg2->AddEntry(g,"ECAL#times0.92+HCAL","P");
  //leg2->AddEntry(gh,"HCAL only","P");
  //leg2->AddEntry(ghc,"HCAL+IsoTrack","P");

  l->DrawLine(3,1,1000,1);
  //l->SetLineStyle(kDotted);
  //l->SetLineWidth(2);
  //l->SetLineColor(kGreen+2);
  //l->DrawLine(3,1+0.03,1000,1+0.03);
  //l->DrawLine(3,1-0.03,1000,1-0.03);
  //l->SetLineColor(kRed);
  //l->DrawLine(3,1+0.03*sqrt(2),1000,1+0.03*sqrt(2));
  //l->DrawLine(3,1-0.03*sqrt(2),1000,1-0.03*sqrt(2));


  // Try mapping from single particle to jet under the assumption that
  // leading particle carries 33% of jet energy, and depth #4 dominated by that
  // So taking jet pt=particle pt / 0.33, fraction x0.33,
  // then adding points from 11% of jet pT for another 33%
  // and finally from 3.7% for final 33%
  TGraphErrors *gjet = new TGraphErrors(g->GetN());
  for (int i = 0; i != g->GetN(); ++i) {
    double ptjet = g->GetX()[i]/0.33;
    double pt1 = ptjet/3.;  int i1 = max(h->GetXaxis()->FindBin(pt1),1);
    double pt3 = ptjet/9.;  int i3 = max(h->GetXaxis()->FindBin(pt3),1);
    double pt9 = ptjet/27.; int i9 = max(h->GetXaxis()->FindBin(pt9),1);
    double f4_1 = (1-h->GetBinContent(i1)) / 3.;
    double f4_3 = (1-h->GetBinContent(i3)) / 3.;
    double f4_9 = (1-h->GetBinContent(i9)) / 3.;
    gjet->SetPoint(i, ptjet, 1 - f4_1 - f4_3 - f4_9);
  }
  gjet->SetMarkerStyle(kFullCircle);
  gjet->Draw("SAMEP");

  TF1 *f4 = new TF1("f4","[0]+log(x)*([1]+log(x)*([2]+log(x)*[3]))",9,3000);
  f4->SetParameters(1,0,0,0);
  gjet->Fit(f4,"QRNW");
  f4->Draw("SAME");
  f4->Print();
  cout << "f4->SetParameters(";
  for (int i = 0; i != f4->GetNpar(); ++i) {
    cout << Form("%s%1.4g",(i==0 ? "" : " ,"),f4->GetParameter(i));
  }
  cout << ");" << endl;

  TLegend *leg2 = tdrLeg(0.20,0.80,0.40,0.90);
  leg2->AddEntry(g,"Particle not in depth #4","PL");
  leg2->AddEntry(gjet,"Jet not in depth #4","P");
  
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.20,0.29,"Jet model:");
  tex->DrawLatex(0.20,0.25,"1/3 of energy in 1 particle of p_{T}/3");
  tex->DrawLatex(0.20,0.21,"1/3 of energy in 3 particles of p_{T}/9");
  tex->DrawLatex(0.20,0.17,"1/3 of energy in 9 particles of p_{T}/27");

  c1->SaveAs("pdf/showerdepth_depositions.pdf");
  c2->SaveAs("pdf/showerdepth_hbtime.pdf");
}
