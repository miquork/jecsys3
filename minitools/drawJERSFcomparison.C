// Purpose: draw comparison of JER SF for |eta|<0.7xx in 23Cv123,23Cv4,23D
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"
#include "../tools.C"

TH1D *profileY(const char *name, TH2D *h2, int ix1, int ix2) {
  TH1D *h1 = h2->ProjectionY(name,ix1,ix2); h1->Scale(1./3.);
  // Do better merging
  for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
    double sumw(0), sum(0), esum2(0);
    for (int i = ix1; i != ix2; ++i) {
      double y = h2->GetBinContent(i,j);
      double ey = h2->GetBinError(i,j);
      double w = (ey!=0 ? 1./pow(ey,2) : 0);
      sumw += w;
      sum += w*y;
      esum2 += pow(w*ey,2);
    } // for i
    double ynew = (sumw!=0 ? sum/sumw : 0);
    double eynew = (sumw!=0 ? sqrt(esum2/sumw)/sqrt(sumw) : 0);
    h1->SetBinContent(j, ynew);
    h1->SetBinError(j, eynew);
  } //for j

  return h1;
} // profileY

void drawJERSFcomparison() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/JERSF_Summer23_L2ResOnly.root","READ");
  assert(f && !f->IsZombie());

  TH2D *h2d(0), *h2dr(0);
  h2d = (TH2D*)f->Get(Form("h2jersf_%s_%s","2023D","Summer23MGBPix"));
  h2dr = (TH2D*)f->Get(Form("h2jersfRaw_%s_%s","2023D","Summer23MGBPix"));
  assert(h2d);
  assert(h2dr);
  
  //TH1D *h1d = h2d->ProjectionY("h1d",1,3); h1d->Scale(1./3.);
  //TH1D *h1dr = h2dr->ProjectionY("h1dr",1,3); h1dr->Scale(1./3.);
  TH1D *h1d = profileY("h1d",h2d,1,3);
  TH1D *h1dr = profileY("h1dr",h2dr,1,3);
  
  TH2D *h2c(0), *h2cr(0);
  h2c = (TH2D*)f->Get(Form("h2jersf_%s_%s","2023Cv123","Summer23MG"));
  h2cr = (TH2D*)f->Get(Form("h2jersfRaw_%s_%s","2023Cv123","Summer23MG"));
  assert(h2c);
  assert(h2cr);

  //TH1D *h1c = h2c->ProjectionY("h1c",1,3); h1c->Scale(1./3.);
  //TH1D *h1cr = h2cr->ProjectionY("h1cr",1,3); h1cr->Scale(1./3.);
  TH1D *h1c = profileY("h1c",h2c,1,3);
  TH1D *h1cr = profileY("h1cr",h2cr,1,3);
  
  TH2D *h2c4(0), *h2c4r(0);
  h2c4 = (TH2D*)f->Get(Form("h2jersf_%s_%s","2023Cv4","Summer23MG"));
  h2c4r = (TH2D*)f->Get(Form("h2jersfRaw_%s_%s","2023Cv4","Summer23MG"));
  assert(h2c4);
  assert(h2c4r);

  //TH1D *h1c4 = h2c4->ProjectionY("h1c4",1,3); h1c4->Scale(1./3.);
  //TH1D *h1c4r = h2c4r->ProjectionY("h1c4r",1,3); h1c4r->Scale(1./3.);
  TH1D *h1c4 = profileY("h1c4",h2c4,1,3);
  TH1D *h1c4r = profileY("h1c4r",h2c4r,1,3);

  
  TFile *f24 = new TFile("rootfiles/JERSF.root","READ");
  assert(f24 && !f24->IsZombie());

  TH2D *h224bc(0), *h224bcr(0);
  h224bc = (TH2D*)f24->Get(Form("Fits/h2jersf_%s_%s","2024BC","Summer23MGBPix"));
  h224bcr = (TH2D*)f24->Get(Form("Dijet/h2jersfRaw_%s_%s","2024BC","Summer23MGBPix"));
  assert(h224bc);
  assert(h224bcr);
  /*
  TH2D *h224c(0), *h224cr(0);
  h224c = (TH2D*)f24->Get(Form("Fits/h2jersf_%s_%s","2024C","Summer23MGBPix"));
  h224cr = (TH2D*)f24->Get(Form("Dijet/h2jersfRaw_%s_%s","2024C","Summer23MGBPix"));
  assert(h224c);
  assert(h224cr);
  */
  int ix2 = h224bc->GetXaxis()->FindBin(h2d->GetXaxis()->GetBinLowEdge(3+1))-1;
  TH1D *h124bc = profileY("h124bc",h224bc,1,ix2);
  TH1D *h124bcr = profileY("h124bcr",h224bcr,1,ix2);
  /*
  int ix2 = h224c->GetXaxis()->FindBin(h2d->GetXaxis()->GetBinLowEdge(3+1))-1;
  TH1D *h124c = profileY("h124c",h224c,1,ix2);
  TH1D *h124cr = profileY("h124cr",h224cr,1,ix2);
  */

  // Produce estimate for 2023D based on 2023Cv123 (+) 5% constant factor
  TH1D *h1d_est = (TH1D*)h1c->Clone("h1d_est");
  TH2D *h2_mc = (TH2D*)f24->Get("Dijet/h2jermcRaw_2024BC_Summer23MGBPix"); assert(h2_mc);
  int i2 = h2_mc->GetXaxis()->FindBin(0.783-0.05);
  TH1D *h1_mc = profileY("h1c_mc",h2_mc,1,i2);
  for (int i = 1; i != h1c->GetNbinsX()+1; ++i) {
    double pt = h1c->GetBinCenter(i);
    double sf = h1c->GetBinContent(i);
    double esf = h1c->GetBinError(i);
    double j = h1_mc->GetXaxis()->FindBin(pt);
    double jer = h1_mc->GetBinContent(j);
    //double c = 0.05*min(max(log(pt)/log(1800.)-log(555)/log(1800.),0.),0.);
    double pt2 = 1300.;
    double pt1 = 50.;
    double c = 0.028*min(max((log(pt)/log(pt2)-log(pt1)/log(pt2))/(1-log(pt1)/log(pt2)),-1.),1.);
    h1d_est->SetBinContent(i, jer>0 ? sqrt(pow(sf,2)+c*fabs(c)/(jer*jer)) : 0.);
    h1d_est->SetBinError(i, jer>0 ? esf : 0);
  } // for i
  
  //TH1D *h = tdrHist("h","JER SF",0.9,1.4);
  TH1D *h = tdrHist("h","JER SF",0.5,3.0);
  //lumi_136TeV = "2023Cv123+2023D";
  lumi_136TeV = "2023-24";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,1,3500,1);
  
  tdrDraw(h1c,"E3",kNone,kBlue,kSolid,-1,1001,kBlue-9);
  tdrDraw(new TGraph(h1c),"L",kNone,kBlue,kSolid,-1);
  h1c->SetFillColorAlpha(kBlue-9,0.70);

  //tdrDraw(h1c4,"E3",kNone,kMagenta,kSolid,-1,1001,kMagenta-9);
  //tdrDraw(new TGraph(h1c4),"L",kNone,kMagenta+2,kSolid,-1);
  //h1c4->SetFillColorAlpha(kMagenta-9,0.70);

  tdrDraw(h1d,"E3",kNone,kBlue,kSolid,-1,1001,kRed-9);
  tdrDraw(new TGraph(h1d),"L",kNone,kRed,kSolid,-1);
  h1d->SetFillColorAlpha(kRed-9,0.70);

  tdrDraw(h1d_est,"E3",kNone,kMagenta+2,kSolid,-1,1001,kMagenta-9);
  tdrDraw(new TGraph(h1d_est),"L",kNone,kMagenta+2,kSolid,-1);
  h1d_est->SetFillColorAlpha(kMagenta-9,0.70);

  tdrDraw(h124bc,"E3",kNone,kGray+2,kSolid,-1,1001,kGray);
  tdrDraw(new TGraph(h124bc),"L",kNone,kBlack,kSolid,-1);
  h124bc->SetFillColorAlpha(kGray,0.70);
  /*
  tdrDraw(h124c,"E3",kNone,kGray+2,kSolid,-1,1001,kGray);
  tdrDraw(new TGraph(h124c),"L",kNone,kBlack,kSolid,-1);
  h124c->SetFillColorAlpha(kGray,0.70);
  */
  tdrDraw(h1cr,"Pz",kFullSquare,kBlue,kSolid,-1,1001,kBlue-9);
  //tdrDraw(h1c4r,"Pz",kOpenSquare,kMagenta+2,kSolid,-1,1001,kMagenta-9);
  tdrDraw(h1dr,"Pz",kFullCircle,kRed,kSolid,-1,1001,kRed-9);
  tdrDraw(h124bcr,"Pz",kFullDiamond,kBlack,kSolid,-1,1001,kGray);
  //tdrDraw(h124cr,"Pz",kFullDiamond,kBlack,kSolid,-1,1001,kGray);

  TLegend *leg = tdrLeg(0.35,0.90-5*0.05,0.60,0.90);
  leg->SetHeader(Form("|#eta| < %1.3f",h2d->GetXaxis()->GetBinLowEdge(3+1)));
  leg->AddEntry(h124bcr,"2024BC (per-depth)","PLEF");
  //leg->AddEntry(h124cr,"2024C (per-depth)","PLEF");
  leg->AddEntry(h1dr,"2023D (per-depth)","PLEF");
  //leg->AddEntry(h1c4r,"2023Cv4 (per-depth)","PLEF");
  leg->AddEntry(h1cr,"2023Cv123","PLEF");
  leg->AddEntry(h1d_est,"2023Cv123 #oplus c_{depth}","LF");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparison24BC.pdf");
  c1->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparison24BC.png");

  /*
  TFile *fa = new TFile("rootfiles/dijet_balance_jer_Summer23BPixPrompt23_AK4Puppi.root","READ");
  assert(fa && !fa->IsZombie());

  TGraphErrors *gd1(0), *gd2(0), *gd3(0), *gm1(0), *gm2(0), *gm3(0);
  TGraphErrors *gd(0), *gm(0), *gr(0), *gr1(0);
  // SM or FE?
  gd1 = (TGraphErrors*)fa->Get("dijet_balance_jer_Data_0p0_0p261_FE_nominal");
  gd2 = (TGraphErrors*)fa->Get("dijet_balance_jer_Data_0p261_0p522_FE_nominal");
  gd3 = (TGraphErrors*)fa->Get("dijet_balance_jer_Data_0p522_0p783_FE_nominal");
  gm1 = (TGraphErrors*)fa->Get("dijet_balance_jer_MC_0p0_0p261_FE_nominal");
  gm2 = (TGraphErrors*)fa->Get("dijet_balance_jer_MC_0p261_0p522_FE_nominal");
  gm3 = (TGraphErrors*)fa->Get("dijet_balance_jer_MC_0p522_0p783_FE_nominal");
  assert(gd1);
  assert(gd2);
  assert(gd3);
  assert(gm1);
  assert(gm2);
  assert(gm3);

  gd = tools::diffGraphs(gd1,gd2,1./3.,-1./3.);
  gd = tools::diffGraphs(gd,gd3,2./3.,-1./3.);
  gm = tools::diffGraphs(gm1,gm2,1./3.,-1./3.);
  gm = tools::diffGraphs(gm,gm3,2./3.,-1./3.);
  
  gr = ratioGraphs(gd,gm);
  gr1 = ratioGraphs(gd1,gm1);

  tdrDraw(gr,"Pz",kFullDiamond,kBlack);
  //tdrDraw(gr1,"Pz",kOpenDiamond,kBlack);
  gr->SetLineWidth(2);

  TF1 *f1 = new TF1("f1","[0]",15,3500);
  f1->SetLineColor(kBlack);
  f1->SetLineStyle(kDotted);
  gr->Fit(f1,"QRN");
  f1->Draw("SAME");

  TF1 *f2 = new TF1("f2","sqrt([0]*[0]/x+[1]*[1]*0.05*0.05)/"
		    "sqrt(1./x+0.05*0.05)",15,3500);
  f2->SetLineColor(kBlack);
  f2->SetLineStyle(kDashed);
  gr->Fit(f2,"QRN");
  f2->Draw("SAME");

  TF1 *f3 = new TF1("f3","sqrt([0]*[0]/x+[1]*[1]*0.05*0.05)/"
		    "sqrt(1./x+0.05*0.05)",300,3500);
  f3->SetParLimits(0,1,3);
  f3->SetParameters(1,1);
  f3->SetLineColor(kBlack);
  f3->SetLineStyle(kSolid);
  gr->Fit(f3,"QRN");
  f3->SetRange(15.,3500.);
  f3->Draw("SAME");
  
  c1->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparison24_withDB.pdf");
  */
} // void drawJERSFcomparison
