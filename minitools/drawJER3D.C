// Purpose: draw JER distribution from TH3D *h3m0 and *h3m2
#include "TFile.h"
#include "TH3D.h"
#include "TF1.h"
#include "TMath.h"

#include "../tdrstyle_mod22.C"

void drawJER3D() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("~/Downloads/jmenano_data_out_2024C_v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF.root","READ");
  //TFile *f = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF/jmenano_data_out_2024C_v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF.root");
  //TFile *f = new TFile("rootfiles/Prompt2024/jmenano_data_out_2024C_JME_v39_2024_Prompt_Golden_29April.root","READ");
  //TFile *f = new TFile("rootfiles/Prompt2024/jmenano_data_out_2024BC_JME_v39_2024_Prompt_Golden_29April.root","READ"); // golden 0.74/fb
  //TFile *f = new TFile("rootfiles/Prompt2024/v41_2024_Golden/jmenano_data_out_2024BC_JME_v41_2024_Golden.root","READ"); // golden 3/fb
  TFile *f = new TFile("rootfiles/Prompt2024/v50_2024/jmenano_data_out_2024BCD_JME_v50_2024.root","READ"); // May 16 golden, 12.3/fb
  assert(f && !f->IsZombie());

  TFile *fm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF/jmenano_mc_cmb_Summer23MG_Cv4_v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF.root","READ");
  assert(fm && !fm->IsZombie());

  //TFile *f3 = new TFile("rootfiles/jmenano_data_cmb_2023D_JME_v46_2023_prompt.root","READ"); // Wrong errors?
  TFile *f3 = new TFile("rootfiles/jmenano_data_out_2023D_JME_v46_2023_prompt.root","READ");
  assert(f3 && !f3->IsZombie());


  curdir->cd();
  
  TH3D *h3m0(0), *h3m2(0), *h3m0b(0);
  h3m0 = (TH3D*)f->Get("HLT_PFJet500/Dijet2/h3m0"); assert(h3m0);
  h3m0b = (TH3D*)f->Get("HLT_PFJet200/Dijet2/h3m0"); assert(h3m0b);
  h3m2 = (TH3D*)f->Get("HLT_PFJet500/Dijet2/h3m2"); assert(h3m2);

  
  h3m0 = (TH3D*)h3m0->Clone("h3m0");
  h3m0b = (TH3D*)h3m0b->Clone("h3m0b");
  h3m2 = (TH3D*)h3m2->Clone("h3m2");

  TH3D *h3m0m(0), *h3m2m(0);
  h3m0m = (TH3D*)fm->Get("Dijet2/h3m0"); assert(h3m0m);
  h3m2m = (TH3D*)fm->Get("Dijet2/h3m2"); assert(h3m2m);

  
  h3m0m = (TH3D*)h3m0m->Clone("h3m0m");
  h3m2m = (TH3D*)h3m2m->Clone("h3m2m");

  TH3D *h3m03(0), *h3m23(0), *h3m0b3(0);
  h3m03 = (TH3D*)f3->Get("HLT_PFJet500/Dijet2/h3m0"); assert(h3m03);
  h3m0b3 = (TH3D*)f3->Get("HLT_PFJet200/Dijet2/h3m0"); assert(h3m0b3);
  h3m23 = (TH3D*)f3->Get("HLT_PFJet500/Dijet2/h3m2"); assert(h3m23);

  h3m03 = (TH3D*)h3m03->Clone("h3m03");
  h3m0b3 = (TH3D*)h3m0b3->Clone("h3m0b3");
  h3m23 = (TH3D*)h3m23->Clone("h3m23");
  
  // x=eta, y=pT, z=MPF/DB
  /*
  // 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783,
  // 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
  int ix1 = h3m0->GetXaxis()->FindBin(0.);
  int ix2 = h3m0->GetXaxis()->FindBin(0.783-0.05);
  //int ix2 = h3m0->GetXaxis()->FindBin(1.305-0.05);
  //97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
  //507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
  //1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
  //2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
  int iy1 = h3m0->GetYaxis()->FindBin(1032.);
  int iy2 = h3m0->GetYaxis()->FindBin(1248.);
  //h3m0->GetXaxis()->SetRange(ix1,ix2);
  //h3m0->GetYaxis()->SetRange(iy1,iy2);
  */

  h3m0->GetXaxis()->SetRangeUser(0,0.783);
  h3m0->GetYaxis()->SetRangeUser(548,592);
  TH1D *h1m0_500 = (TH1D*)h3m0->Project3D("z");
  h1m0_500->SetName("h1m0_500");
  h1m0_500->Scale(1./h1m0_500->Integral());
  
  h3m0->GetXaxis()->SetRangeUser(0,0.783);
  h3m0->GetYaxis()->SetRangeUser(1032,1248);
  TH1D *h1m0_1000 = (TH1D*)h3m0->Project3D("z");
  h1m0_1000->SetName("h1m0_1000");
  h1m0_1000->Scale(1./h1m0_1000->Integral());

  h3m0b->GetXaxis()->SetRangeUser(0,0.783);
  h3m0b->GetYaxis()->SetRangeUser(245,300);
  TH1D *h1m0_250 = (TH1D*)h3m0b->Project3D("z");
  h1m0_250->SetName("h1m0_250");
  h1m0_250->Scale(1./h1m0_250->Integral());


  h3m0m->GetXaxis()->SetRangeUser(0,0.783);
  h3m0m->GetYaxis()->SetRangeUser(548,592);
  TH1D *h1m0_500m = (TH1D*)h3m0m->Project3D("z");
  h1m0_500m->SetName("h1m0_500m");
  h1m0_500m->Scale(1./h1m0_500m->Integral());
  
  h3m0m->GetXaxis()->SetRangeUser(0,0.783);
  h3m0m->GetYaxis()->SetRangeUser(1032,1248);
  TH1D *h1m0_1000m = (TH1D*)h3m0m->Project3D("z");
  h1m0_1000m->SetName("h1m0_1000m");
  h1m0_1000m->Scale(1./h1m0_1000m->Integral());

  h3m0m->GetXaxis()->SetRangeUser(0,0.783);
  h3m0m->GetYaxis()->SetRangeUser(245,300);
  TH1D *h1m0_250m = (TH1D*)h3m0m->Project3D("z");
  h1m0_250m->SetName("h1m0_250m");
  h1m0_250m->Scale(1./h1m0_250m->Integral());
  
  h3m03->GetXaxis()->SetRangeUser(0,0.783);
  h3m03->GetYaxis()->SetRangeUser(548,592);
  TH1D *h1m03_500 = (TH1D*)h3m03->Project3D("z");
  h1m03_500->SetName("h1m03_500");
  h1m03_500->Scale(1./h1m03_500->Integral());
  
  h3m03->GetXaxis()->SetRangeUser(0,0.783);
  h3m03->GetYaxis()->SetRangeUser(1032,1248);
  TH1D *h1m03_1000 = (TH1D*)h3m03->Project3D("z");
  h1m03_1000->SetName("h1m03_1000");
  h1m03_1000->Scale(1./h1m03_1000->Integral());

  h3m0b3->GetXaxis()->SetRangeUser(0,0.783);
  h3m0b3->GetYaxis()->SetRangeUser(245,300);
  TH1D *h1m03_250 = (TH1D*)h3m0b3->Project3D("z");
  h1m03_250->SetName("h1m03_250");
  h1m03_250->Scale(1./h1m03_250->Integral());
  
  
  //h3m2->GetXaxis()->SetRangeUser(0,0.783);
  //h3m2->GetYaxis()->SetRangeUser(1032,1248);
  //TH1D *h1m2 = (TH1D*)h3m2->Project3D("z");

  //h1m0_1000->Draw();
  //h1m0_500->Draw("SAME");
  //h1m2->Draw("SAME");

  TH1D *h = tdrHist("h","Fraction of N_{events}",1.5e-4,7e-2,
		    "MPF",0.2,2.3);
  lumi_136TeV = "2024BCD, 12.3 fb^{-1}";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  tdrDraw(h1m0_500m,"HIST",kOpenDiamond,kBlue,kSolid,-1,1001,kBlue-9);
  tdrDraw(h1m0_250m,"HIST",kOpenDiamond,kGreen+2,kSolid,-1,1001,kGreen+2-9);
  h1m0_250m->SetFillColorAlpha(kGreen+2-9,0.5);
  //tdrDraw(h1m0_1000m,"Pz",kOpenCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(h1m0_1000m,"HIST",kOpenCircle,kRed,kSolid,-1,1001,kRed-9);
  h1m0_1000m->SetFillColorAlpha(kRed-9,0.5);
  h1m0_1000m->SetMarkerSize(0.5);

  tdrDraw(h1m0_250,"Pz",kFullDiamond,kGreen+2,kSolid,-1,kNone);
  tdrDraw(h1m0_500,"Pz",kFullDiamond,kBlue,kSolid,-1,kNone);
  tdrDraw(h1m0_1000,"Pz",kFullCircle,kRed,kSolid,-1,kNone);
  h1m0_250->SetMarkerSize(0.5);
  h1m0_500->SetMarkerSize(0.5);
  h1m0_1000->SetMarkerSize(0.5);

  tdrDraw(h1m03_250,"Pz",kOpenDiamond,kGreen+2,kSolid,-1,kNone);
  tdrDraw(h1m03_500,"Pz",kOpenDiamond,kBlue,kSolid,-1,kNone);
  tdrDraw(h1m03_1000,"Pz",kOpenCircle,kRed,kSolid,-1,kNone);
  h1m03_250->SetMarkerSize(0.5);
  h1m03_500->SetMarkerSize(0.5);
  h1m03_1000->SetMarkerSize(0.5);
  
  /*
  tdrDraw(h1m0_500,"HIST",kNone,kBlue,kSolid,-1,1001,kBlue-9);
  tdrDraw(h1m0_250,"HIST",kNone,kGreen+2,kSolid,-1,1001,kGreen+2-9);
  h1m0_250->SetFillColorAlpha(kGreen+2-9,0.5);
  tdrDraw(h1m0_1000,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  h1m0_1000->SetMarkerSize(0.5);
  */
  
  double rms250 = h1m0_250->GetRMS();
  TF1 *f250 = new TF1("f250","gaus",1-1.5*rms250,1+1.5*rms250);
  h1m0_250->Fit(f250,"QRN");
  f250->SetLineColor(kGreen+2);
  double sigma250 = f250->GetParameter(2);
  f250->SetRange(1-3.2*sigma250,1+3.2*sigma250);
  f250->SetNpx(400);
  f250->Draw("SAME");
  
  double rms500 = h1m0_500->GetRMS();
  TF1 *f500 = new TF1("f500","gaus",1-1.5*rms500,1+1.5*rms500);
  h1m0_500->Fit(f500,"QRN");
  f500->SetLineColor(kBlue);
  double sigma500 = f500->GetParameter(2);
  f500->SetRange(1-3.2*sigma500,1+3.2*sigma500);
  f500->SetNpx(400);
  f500->Draw("SAME");
  
  double rms1000 = h1m0_1000->GetRMS();
  TF1 *f1000 = new TF1("f1000","gaus",1-1.5*rms1000,1+1.5*rms1000);
  h1m0_1000->Fit(f1000,"QRN");
  f1000->SetLineColor(kRed);
  double sigma1000 = f1000->GetParameter(2);
  f1000->SetRange(1-3.2*sigma1000,1+3.2*sigma1000);
  f1000->SetNpx(400);
  f1000->Draw("SAME");

  TF1 *f1000x2 = new TF1("f1000x2","[0]*TMath::Gaus(x,[1],[2],{kTRUE})+[3]*TMath::Gaus(x,[4],[5],{kTRUE})",1-6.4*sigma1000,1+6.4*sigma1000);
  f1000x2->SetParameters(f1000->GetParameter(0),f1000->GetParameter(1),
			 f1000->GetParameter(2),
			 0.1*f1000->GetParameter(0),f1000->GetParameter(0),
			 2.*f1000->GetParameter(0));
  h1m0_1000->Fit(f1000x2,"QRN");
  f1000x2->SetLineColor(kRed+1);
  f1000x2->SetLineWidth(2);
  f1000x2->SetNpx(400);
  f1000x2->Draw("SAME");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->SetTextColor(kRed+1);
  tex->DrawLatex(0.75,0.30,Form("N_{2}=%1.2f",f1000x2->GetParameter(3)/(f1000x2->GetParameter(0)+f1000x2->GetParameter(3))));
  tex->DrawLatex(0.75,0.25,Form("#sigma_{2}=%1.2f",f1000x2->GetParameter(5)/f1000x2->GetParameter(2)));
  
  TLegend *leg = tdrLeg(0.55,0.90-0.05*6,0.80,0.90);
  leg->SetHeader("|#eta| < 0.783");
  leg->AddEntry(h1m0_1000,"2024BCD Data","PLE");
  leg->AddEntry(h1m03_1000,"2023D Data","PLE");
  //leg->AddEntry(h1m0_1000m,"[1032,1248] GeV","PLEF");
  leg->AddEntry(h1m0_1000m,"[1032,1248] GeV","F");
  leg->AddEntry(h1m0_500m,"[548,592] GeV","F");
  leg->AddEntry(h1m0_250m,"[245,300] GeV","F");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawJER3D/drawJER3D_lin_2024BCD.pdf");
  gPad->SetLogy();
  c1->SaveAs("pdf/drawJER3D/drawJER3D_log_2024BCD.pdf");
} // void drawJER3D
