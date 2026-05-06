// Purpose: Use Ravindra's GJet+QCD vs GJet to parameterize background
//          correction for gamma+jet results from Bettina.
//          Stage 2: use GJet+QCD vs GJet flavor fractions to propagete
//          gluonJES from data to background correction
#include "TFile.h"
#include "TProfile.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

void photonBkgCorrections(string method = "MPF");
TH1D *getHDM(TFile *f, string sample);
void photonBkgCorrectionForGluonJES();
  
void photonBkgCorrection() {
  /*
  photonBkgCorrections("MPF");
  photonBkgCorrections("DB");
  photonBkgCorrections("HDM");
  */
  photonBkgCorrectionForGluonJES();
}

void photonBkgCorrections(string method) {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/photonBkgCorrection");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/photonBkgCorrection");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/Ravindra/GamJet_vs_GamJetFake_2024.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();
  
  if (method=="MPF" || method=="HDM") {
    f->cd("plot_006_passL3Residual_MpfResp");
  }
  else if (method=="DB") {
    f->cd("plot_005_passL3Residual_DbResp");
  }
  gDirectory->cd("top_pad");
  TDirectory *d = gDirectory;

  TProfile *pmix = (TProfile*)d->Get("hMCGJetsPlusMCQCD"); assert(pmix);
  TProfile *pref = (TProfile*)d->Get("hMCGJets"); assert(pref);

  TH1D *hmix = pmix->ProjectionX("hmix");
  TH1D *href = pref->ProjectionX("href");

  if (method=="HDM") {
    hmix = getHDM(f,"hMCGJetsPlusMCQCD");
    href = getHDM(f,"hMCGJets");
  }
  
  TH1D *hr = (TH1D*)hmix->Clone("hr");
  hr->Divide(href);

  // Binomial statistics should give more correct uncertainty
  //TH1D *hrinv = (TH1D*)href->Clone("hrinv");
  //hrinv->Divide(href,hmix,1,1,"B");
  //TH1D *hrb = (TH1D*)hr->Clone("hrb");
  //for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
  //hrb->SetBinError(i, hrinv->GetBinError(i));
  //}

  const char *cm = method.c_str();
  double xmin(30), xmax(2000);
  TH1D *hu = tdrHist("hu",cm,0.985,1.085,"p_{T,#gamma} (GeV)",xmin,xmax);
  TH1D *h = tdrHist("hd","Ratio",0.99,1.005,"p_{T,#gamma} (GeV)",xmin,xmax);
  lumi_136TeV = "Summer24 MC";
  TCanvas *c1 = tdrDiCanvas("c1",hu,h,8,11);

  if (method=="DB") {
    hu->GetYaxis()->SetRangeUser(0.86,1.26);
    h->GetYaxis()->SetRangeUser(0.98,1.03);
  }
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  
  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(xmin,1,xmax,1);

  tdrDraw(href,"HIST",kNone,kRed,kSolid,-1,kNone,0);
  tdrDraw(hmix,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
  hmix->GetXaxis()->SetRangeUser(xmin,xmax);

  TLegend *leg = tdrLeg(0.45,0.80-0.05*2,0.70,0.80);
  leg->AddEntry(href,"Pure #gamma+jet","FL");
  leg->AddEntry(hmix,"Mixed #gamma+jet and QCD","PLE");

  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(xmin,1,xmax,1);
    
  //tdrDraw(hr,"Pz",kNone,kRed,kSolid,-1,kNone,0);
  //tdrDraw(hrb,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

  TF1 *f1 = new TF1("f1","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  if (method=="DB") f1->SetRange(40.,xmax);
  f1->SetParameters(1,-0.01,-0.3,0);
  hr->Fit(f1,"RN");
  //f1->FixParameter(3, f1->GetParameter(3));
  //hr->Fit(f1,"RN");
  //f1->FixParameter(2, f1->GetParameter(2));
  //hr->Fit(f1,"RN");
  if (method=="DB") f1->SetLineColor(kRed);
  if (method=="MPF") f1->SetLineColor(kBlue);
  if (method=="HDM") f1->SetLineColor(kGreen+2);
  f1->Draw("SAME");

  // DB
  TF1 *f1d = new TF1("f1d","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1d->SetParameters(0.9481, 0.097, -0.067, -1.082);
  f1d->SetLineStyle(kDotted);
  f1d->SetLineColor(kRed);
  f1d->Draw("SAME");

  // MPF
  TF1 *f1m = new TF1("f1m","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1m->SetParameters(1.0034, -0.128, -0.482, 0.683);
  f1m->SetLineStyle(kDotted);
  f1m->SetLineColor(kBlue);
  f1m->Draw("SAME");

  // HDM
  TF1 *f1h = new TF1("f1h","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1h->SetParameters(1.0030, -0.366, -0.624, 1.326);
  f1h->SetLineStyle(kDotted);
  f1h->SetLineColor(kGreen+2);
  f1h->Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);
  tex->DrawLatex(0.40,0.80,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f)",
				f1->GetChisquare(), f1->GetNDF(),
				f1->GetChisquare() / f1->GetNDF()));
  tex->DrawLatex(0.40,0.42,Form("%s",f1->GetExpFormula().Data()));
  tex->DrawLatex(0.40,0.35,Form("p0: %1.4f, p1: %1.3f, p2: %1.3f, p4: %1.3f",
				f1->GetParameter(0), f1->GetParameter(1),
				f1->GetParameter(2), f1->GetParameter(3)));

  gPad->RedrawAxis();


  c1->SaveAs(Form("pdf/photonBkgCorrection/photonBkgCorrection_%s.pdf",cm));
} // photonBkgCorrection


TH1D *getHDM(TFile *f, string sample) {

  const char *c = sample.c_str();
  TProfile *p0, *p1, *pn, *pu;
  p0 = (TProfile*)f->Get(Form("plot_006_passL3Residual_MpfResp/top_pad/%s",c));
  p1 = (TProfile*)f->Get(Form("plot_005_passL3Residual_DbResp/top_pad/%s",c));
  pn = (TProfile*)f->Get(Form("plot_007_passL3Residual_MpfnResp/top_pad/%s",c));
  pu = (TProfile*)f->Get(Form("plot_008_passL3Residual_MpfuResp/top_pad/%s",c));
  assert(p0);
  assert(p1);
  assert(pn);
  assert(pu);

  TH1D *h = p0->ProjectionX(Form("hdm_%s",c)); h->Reset();

  // Following HDM code is mostly copied from softrad3.C
  
  // Helper function for solving jet response x from MPF decomposition
  TF1 *fm = new TF1("fm","x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])",0,13000);
  TF1 *fp = new TF1("fp","x - ([0] + (x/[3]  )*[1] + (x/[4]  )*[2])",0,13000);

  double xmin(0.1), xmax(1.9);
  double Rn_m(1.000), Ru_m(0.92);
  for (int i = 1; i != p0->GetNbinsX()+1; ++i) {

    double r0 = p0->GetBinContent(i);
    double er0 = p0->GetBinError(i);
    double ru = pu->GetBinContent(i);
    double rn = pn->GetBinContent(i);
    if (r0!=0 && er0>0) {
      
      // Solve master equation numerically
      // R1 = r0 + [(R1-Rn)/Rn] * rn + [(R1-Ru)/Rn] * ru
      // => fm = "x - ([0] + (x/[3]-1)*[1] + (x/[4]-1)*[2])"
      fm->SetParameters(r0, rn, ru, Rn_m, Ru_m);
      double R1 = fm->GetX(0,xmin,xmax);

      h->SetBinContent(i, R1);
      h->SetBinError(i, er0);
    }
  } // for i
  
  // Delete to avoid errors from same name later
  delete fm;
  delete fp;

  return h;
} // getHDM


void photonBkgCorrectionForGluonJES() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *fg = new TFile("rootfiles/Ravindra/flavour_fraction_EMEnrichedQCD.root","READ");
  assert(fg && !fg->IsZombie());

  TFile *fq = new TFile("rootfiles/Ravindra/GamJet_vs_GamJetFake_2024.root","READ");
  assert(fq && !fq->IsZombie());

  curdir->cd();

  TH1D *hgq = (TH1D*)fg->Get("passDPhiTagProbe_WhenPho/HistFlavor_ProbeInTagPt/hFrac_Flavor_g"); assert(hgq);
  TH1D *hgp = (TH1D*)hgq->Clone("hgp"); hgp->Clear();
  for (int i = 1; i != hgp->GetNbinsX()+1; ++i) {
    hgp->SetBinContent(i, 0.20);
    hgp->SetBinError(i, 0.2*hgq->GetBinError(i));
  }

  // 004 seems too high QCD fraction (flat 40%)? And GJet drops very steep
  // 003 is higher on GJet, 002 low again, 001 very high
  // Picking 003 that has better behaved QCD fraction
  // 001 has almost 100% QCD at low pT, only 10% at high pT
  // 002 fraction similar to 003, bit lower at high pT
  TH1D *hqq = (TH1D*)fq->Get("plot_003_passL3Residual_Events/top_pad/hMCQCD");
  TH1D *hqb = (TH1D*)fq->Get("plot_003_passL3Residual_Events/top_pad/hMCGJetsPlusMCQCD");
  assert(hqq);
  assert(hqb);
  TH1D *hq = (TH1D*)hqq->Clone("hq");
  hq->Divide(hqb);

  // Flavor fit parameters from minitools/Zflavor.C (Eta13)
  TF1 *f1q = new TF1("f1q","[0]+[1]*(pow(0.01*x,[2])-1)",45,300);
  TF1 *f1g = new TF1("f1g","[0]+[1]*(pow(0.01*x,[2])-1)",45,300);
  TF1 *f1c = new TF1("f1c","[0]",45,300);
  TF1 *f1b = new TF1("f1b","[0]",45,300);
  TF1 *f1z = new TF1("f1z","[0]+[1]*(pow(0.01*x,[2])-1)",45,300);
  f1q->SetParameters(0.1933, -0.4463, -1); // chi2/NDF=3.5/7
  f1g->SetParameters(-0.8837, -0.6869, -1); // chi2/NDF=27.3/7
  f1c->SetParameter(0,0.09416); // chi2/NDF=13.4/8
  f1b->SetParameter(0,1.75); // chi2/NDF=170.2/8
  f1z->SetParameters(0, 0.5, -1); // chi2/NDF=0.0/0
  // pdf/Zflavor/Zflavor_datamc_Eta13_2025CDEFG.pdf has been created

  TH1D *hd = (TH1D*)hq->Clone("hd"); hd->Reset();
  for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
    double pt = hd->GetBinCenter(i);
    int j = hgq->GetXaxis()->FindBin(pt);
    double f_gluon_qcd = hgq->GetBinContent(i);
    double ef_gluon = hgq->GetBinError(i);
    double f_gluon_ref = 0.20;
    double f_qcd = hq->GetBinContent(i);
    double ef_qcd = hq->GetBinError(i);
    double r_gluon = 1+0.01*f1g->Eval(pt);
    double r_quark = 1+0.01*f1q->Eval(pt);
    double er_gluon = 0.01*0.1;

    double df_gluon = f_gluon_qcd-f_gluon_ref;
    double dr_gluon = r_gluon-r_quark;

    double d_jes = df_gluon * f_qcd * dr_gluon;
    double ed_jes = sqrt(pow(ef_gluon*f_qcd*dr_gluon,2) +
			 pow(df_gluon*ef_qcd*dr_gluon,2) +
			 pow(df_gluon*f_qcd*er_gluon,2));
    hd->SetBinContent(i, 100.*d_jes);
    hd->SetBinError(i, 100.*ed_jes);
  }

  double xmin(15), xmax(3500), eps(1e-4);
  TH1D *hu = tdrHist("hu2","Fraction",0+eps,1-eps,"p_{T} (GeV)",xmin,xmax);
  TH1D *h = tdrHist("hd2","#DeltaJES (%)",-0.5,0.3,"p_{T} (GeV)",xmin,xmax);
  lumi_136TeV = "Summer24 MC";
  TCanvas *c1 = tdrDiCanvas("c2",hu,h,8,11);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  c1->cd(1);
  gPad->SetLogx();

  tdrDraw(hq,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0);
  tdrDraw(hgq,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hgp,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0);

  TLegend *leg = tdrLeg(0.37,0.90-0.05*3,0.62,0.90);
  leg->AddEntry(hgq,"Gluon fraction in QCD MC","PLE");
  leg->AddEntry(hgp,"Gluon fraction in #gamma+jet MC","PLE");
  leg->AddEntry(hq,"Fraction of QCD in mixed MC","PLE");

  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(xmin,0,xmax,0);

  tdrDraw(hd,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);

  TF1 *f1 = new TF1("f2","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1->SetParameters(-0.277, 0.0760, 0.177, -4.57);
  hd->Fit(f1,"RN");
  f1->SetLineColor(kMagenta+1);
  f1->Draw("SAME");
  
  // HDM
  TF1 *f1h = new TF1("f2h","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1h->SetParameters((1.0030-1)*100., -0.366*100., -0.624, 1.326*100);
  f1h->SetLineStyle(kDotted);
  f1h->SetLineColor(kGreen+2);
  f1h->Draw("SAME");

  // Data gluon JES
  TF1 *f1d = new TF1("f1d","[0]+[1]*pow(x,[2])+[3]/x",xmin,xmax);
  f1d->SetParameters(-0.277, 0.0760, 0.177, -4.57);
  f1d->SetLineStyle(kDotted);
  f1d->SetLineColor(kMagenta+2);
  f1h->SetLineWidth(2);
  f1d->Draw("SAME");

  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);
  tex->DrawLatex(0.25,0.80,"Data gluon JES induced bias on #gamma+jet HDM");

  tex->DrawLatex(0.37,0.73,Form("#chi^{2} / NDF = %1.1f / %d (%1.1f)",
				f1->GetChisquare(), f1->GetNDF(),
				f1->GetChisquare() / f1->GetNDF()));
  tex->DrawLatex(0.37,0.42,Form("%s",f1->GetExpFormula().Data()));
  tex->DrawLatex(0.37,0.35,Form("p0: %1.3f, p1: %1.4f, p2: %1.3f, p4: %1.2f",
				f1->GetParameter(0), f1->GetParameter(1),
				f1->GetParameter(2), f1->GetParameter(3)));
  
  c1->SaveAs("pdf/photonBkgCorrection/photonBkgCorrectionForGluonJES.pdf");
} // photonBkgCorrectionForGluonJES
