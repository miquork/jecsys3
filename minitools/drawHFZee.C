// Purpose: Draw HF Zee data to demonstrate pT threshold bias at high |eta|
#include "TFile.h"

#include "../tdrstyle_mod22.C"

void drawHFZee() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Read in histograms from James Vincent Natoli
  TFile *fm = new TFile("rootfiles/ZeeFromJamesNatoli/output_MC_DYto2L-4Jets_MLL-50_Winter23_PLOTS_mc_noPU.root","READ");
  assert(fm && !fm->IsZombie());

  TFile *fd = new TFile(" rootfiles/ZeeFromJamesNatoli/output_data_EGamma_Run2023D_BPix_PLOTS_data_noPU.root","READ");
  assert(fd && !fd->IsZombie());

  curdir->cd();

  // Merge plus and minus sides for better statistics
  TH1D *hm1p = (TH1D*)fm->Get("etaPlus31"); assert(hm1p);
  TH1D *hm1m = (TH1D*)fm->Get("etaMinus31"); assert(hm1m);
  TH1D *hm1 = (TH1D*)hm1p->Clone("hm1");
  hm1->Add(hm1m);

  TH1D *hm2p = (TH1D*)fm->Get("etaPlus37"); assert(hm2p);
  TH1D *hm2m = (TH1D*)fm->Get("etaMinus37"); assert(hm2m);
  TH1D *hm2 = (TH1D*)hm2p->Clone("hm2");
  hm2->Add(hm2m);

  TH1D *hd1p = (TH1D*)fd->Get("etaPlus31"); assert(hd1p);
  TH1D *hd1m = (TH1D*)fd->Get("etaMinus31"); assert(hd1m);
  TH1D *hd1 = (TH1D*)hd1p->Clone("hd1");
  hd1->Add(hd1m);

  TH1D *hd2p = (TH1D*)fd->Get("etaPlus37"); assert(hd2p);
  TH1D *hd2m = (TH1D*)fd->Get("etaMinus37"); assert(hd2m);
  TH1D *hd2 = (TH1D*)hd2p->Clone("hd2");
  hd2->Add(hd2m);

  // Rebin for better statistics
  hm1->Rebin(2);
  hm2->Rebin(2);
  hd1->Rebin(2);
  hd2->Rebin(2);

  // Normalize by upper half times HF scale
  // (to account for scale times resolution compressing peak higher)
  // Use posterior HF scales from minitools/drawHF.C estimates
  double mz = 91.2;
  double mzmax = mz+30;//120.;
  double km1 = 82./91.2;
  int im1_1 = hm1->FindBin(km1*mz);
  int im1_2 = hm1->FindBin(km1*mzmax);
  hm1->Scale(0.5/hm1->Integral(im1_1,im1_2) * km1);

  double km2 = 82./91.2;
  int im2_1 = hm2->FindBin(km2*mz);
  int im2_2 = hm2->FindBin(km2*mzmax);
  hm2->Scale(0.5/hm2->Integral(im2_1,im2_2) * km2);

  double kd1 = km1*sqrt(0.93);
  int id1_1 = hd1->FindBin(kd1*mz);
  int id1_2 = hd1->FindBin(kd1*mzmax);
  hd1->Scale(0.5/hd1->Integral(id1_1,id1_2) * kd1);

  double kd2 = km2*sqrt(0.60);
  int id2_1 = hd2->FindBin(kd2*mz);
  int id2_2 = hd2->FindBin(kd2*mzmax);
  hd2->Scale(0.5/hd2->Integral(id2_1,id2_2) * kd2);
  
  
  TH1D *h = tdrHist("h","Zee counts (a.u.)",0,0.12,"M_{ee} (GeV)",40,140);
  lumi_136TeV = "2023D";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  tdrDraw(hm1,"HISTE",kNone,kBlue+1,kSolid,-1,1001,kBlue-9);
  tdrDraw(hm2,"HISTE",kNone,kRed+1,kSolid,-1,1001,kRed-9);
  hm2->SetFillColorAlpha(kRed-9,0.7);
  
  tdrDraw(hd1,"Pz",kOpenCircle,kBlue+1,kSolid,-1,1001,kBlue-9);
  tdrDraw(hd2,"Pz",kFullCircle,kRed+1,kSolid,-1,1001,kRed-9);

  // Initial MC fit
  TF1 *fm2 = new TF1("fm2","gaus",km2*(mz-15.),km2*(mz+15.));
  hm2->Fit(fm2,"QRN");
  fm2->SetLineWidth(2);
  fm2->SetLineColor(kRed);
  fm2->SetLineStyle(kDashed);
  fm2->DrawClone("SAME");

  // Reference MC fit using reliable range at Mee>72 GeV
  fm2->SetRange(72.,km2*(mz+15.));
  hm2->Fit(fm2,"QRN");
  fm2->SetLineStyle(kSolid);
  fm2->SetRange(km2*(mz-15.),km2*(mz+15.));
  fm2->DrawClone("SAME");
  fm2->SetLineStyle(kDashed);

  // Biased data fit in full range mZ+/-15 GeV
  TF1 *fd2 = new TF1("fd2","gaus",kd2*(mz-15.),kd2*(mz+15.));
  hd2->Fit(fd2,"QRN");
  fd2->SetLineWidth(2);
  fd2->SetLineColor(kRed);
  fd2->SetLineStyle(kDashed);
  fd2->DrawClone("SAME");

  // Unbiased data shape estimated from reference MC times data HF scale
  fd2->SetParameters(fm2->GetParameter(0),
		     fm2->GetParameter(1)*kd2/km2,
		     fm2->GetParameter(2)*kd2/km2);
  fd2->SetLineStyle(kSolid);
  fd2->SetRange(kd2*(mz-15.),kd2*(mz+15.));
  fd2->DrawClone("SAME");

  TLegend *leg = tdrLeg(0.55,0.90-6*0.045,0.80,0.90);
  leg->AddEntry(hd2,"Data |ieta|=37","PLE");
  leg->AddEntry(hd1,"Data |ieta|=31","PLE");
  leg->AddEntry(hm2,"MC |ieta|=37","F");
  leg->AddEntry(hm1,"MC |ieta|=31","F");
  leg->AddEntry(fm2,"Fit [m_{Z}-15,m_{Z}+15]","L");
  leg->AddEntry(fd2,"Fit [72,m_{Z}+15]","L");

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(72,0,72,0.12);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.60,0.55,"Normalized by");
  tex->DrawLatex(0.60,0.50,"[m_{Z},m_{Z}+30 GeV]");
  tex->DrawLatex(0.60,0.45,"#times HF scale #times 0.5");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawHFZee.pdf");
} // draw HFZee
