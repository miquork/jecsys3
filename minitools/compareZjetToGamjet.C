// Purpose: compare Z+jet to Gam+jet at different stages of bias corrections:
//          1) QCD background for gamma+jet
//          2) +Gluon JES from data for gamma+jet
//          [3) Gluon JES from  data for Z+jet]
#include "TFile.h"

#include "../tdrstyle_mod22.C"

void compareZjetToGamjet() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/compareZjetToGamjet");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/compareZjetToGamjet");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f0 = new TFile("rootfiles/jecdata20205CDEFG_scaleEMtoOne.root","READ");
  assert(f0 && !f0->IsZombie());

  TFile *f1 = new TFile("rootfiles/jecdata20205CDEFG_scaleEMforQCDBkgInHDM.root","READ");
  assert(f1 && !f1->IsZombie());

  TFile *f2 = new TFile("rootfiles/jecdata20205CDEFG_scaleEMforDataGluonJES.root","READ");
  assert(f2 && !f2->IsZombie());

  curdir->cd();

  TH1D *hz = (TH1D*)f0->Get("ratio/eta00-13/hdm_mpfchs1_zjet"); assert(hz);
  TH1D *hg = (TH1D*)f0->Get("ratio/eta00-13/hdm_mpfchs1_gamjet"); assert(hg);
  TH1D *hg1 = (TH1D*)f1->Get("ratio/eta00-13/hdm_mpfchs1_gamjet"); assert(hg1);
  TH1D *hg2 = (TH1D*)f2->Get("ratio/eta00-13/hdm_mpfchs1_gamjet"); assert(hg2);

  TH1D *hzr = (TH1D*)hz->Clone("hzr"); hzr->Reset();//Divide(hz);
  for (int i = 1; i != hzr->GetNbinsX()+1; ++i) {
    double pt = hzr->GetBinCenter(i);
    double rz = hz->Interpolate(pt);
    if (rz!=0) {
      hzr->SetBinContent(i, (hz->GetBinContent(i)/rz-1)*100);
      hzr->SetBinError(i, hz->GetBinError(i)/rz*100);
    }
  } // for i
  
  TH1D *hgr = (TH1D*)hg->Clone("hgr"); hgr->Reset();
  TH1D *hg1r = (TH1D*)hg1->Clone("hg1r"); hg1r->Reset();
  TH1D *hg2r = (TH1D*)hg2->Clone("hg2r"); hg2r->Reset();
  for (int i = 1; i != hgr->GetNbinsX()+1; ++i) {
    double pt = hgr->GetBinCenter(i);
    double rz = hz->Interpolate(pt);
    if (rz!=0) {
      hgr->SetBinContent(i, (hg->GetBinContent(i)/rz-1)*100);
      hgr->SetBinError(i, hg->GetBinError(i)/rz*100);
      hg1r->SetBinContent(i, (hg1->GetBinContent(i)/rz-1)*100);
      hg1r->SetBinError(i, hg1->GetBinError(i)/rz*100);
      hg2r->SetBinContent(i, (hg2->GetBinContent(i)/rz-1)*100);
      hg2r->SetBinError(i, hg2->GetBinError(i)/rz*100);
    }
  } // for i

  
  double xmin(30), xmax(1500), eps(1e-4), ydmin(-1.0+eps), ydmax(1.0+eps);
  TH1D *h = tdrHist("h","HDM response",0.97+eps,1.03-eps,
		    "p_{T,ref} (GeV)",xmin,xmax);
  TH1D *hd = tdrHist("hd","#DeltaHDM (%)",ydmin,ydmax,"p_{T,ref} (GeV)",xmin,xmax);
  lumi_136TeV = "2025CDEFG, 109 fb^{-1}";
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  l->DrawLine(xmin,1,xmax,1);
  
  hz->GetXaxis()->SetRangeUser(xmin,xmax);
  
  tdrDraw(hz,"HISTE][",kNone,kRed,kSolid,-1,kNone,0);
  tdrDraw(hg,"HISTE][",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hg1,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hg2,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

  TLegend *leg = tdrLeg(0.45,0.85-0.05*4,0.70,0.85);
  leg->AddEntry(hz,"Z+jet","LE");
  leg->AddEntry(hg,"#gamma+jet (uncorrected)","LE");
  leg->AddEntry(hg1,"#gamma+jet + QCD bkg corr.","PLE");
  leg->AddEntry(hg2,"         + gluonJES corr.","PLE");

  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(xmin,0,xmax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,-0.1,xmax,-0.1);
  l->DrawLine(xmin,+0.1,xmax,+0.1);
  l->DrawLine(90,ydmin,90,ydmax);
  l->DrawLine(400,ydmin,400,ydmax);
  
  hzr->GetXaxis()->SetRangeUser(xmin,xmax);

  tdrDraw(hzr,"HISTE][",kNone,kRed,kSolid,-1,kNone,0);
  tdrDraw(hgr,"HISTE][",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hg1r,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0,0.7);
  tdrDraw(hg2r,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.2);
  tex->SetTextColor(kGray+1);
  tex->DrawLatex(0.19,0.85,"Z gluon grows");
  tex->DrawLatex(0.70,0.85,"#gamma gain1 grows");

  tex->SetTextColor(kGreen+2);
  tex->DrawLatex(0.405,0.85,"Golden region 90-400 GeV");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.40,0.78,"!!!Agreement within 0.1%!!!");

  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/compareZjetToGamjet/compareZjetToGamjet.pdf");
} // compareZjetToGamjet
