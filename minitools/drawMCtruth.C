// Purpose: draw MC truth from HLT_MC/MCtruth folder
//          compare Winter24 to Summer23 for response
#include "TFile.h"
#include "TProfile2D.h"

#include "../tdrstyle_mod22.C"

void drawMCtruth() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fref = new TFile("rootfiles/Summer23_L2L3Res/jmenano_mc_out_Summer23MG_v39_noRwPU_SmearJets_L2Res_v1_SF_Cv4.root","READ");
  assert(fref && !fref->IsZombie());
  TDirectory *dref = fref->GetDirectory("HLT_MC/MCtruth");
  
  TFile *fw24 = new TFile("rootfiles/Winter24Flat/v39_Winter24_noJetveto_SmearJets/jmenano_mc_out_Winter24MCFlat_v39_Winter24_noJetveto_SmearJets.root","READ");
  assert(fw24 && !fw24->IsZombie());
  TDirectory *dw24 = fw24->GetDirectory("HLT_MC/MCtruth");

  TFile *fs23 = new TFile("rootfiles/Winter24Flat/v39_Summer23Flat_noJetveto_SmearJets/jmenano_mc_out_Summer23MCFlat_v39_Summer23Flat_noJetveto_SmearJets.root","READ");
  assert(fs23 && !fs23->IsZombie());
  TDirectory *ds23 = fs23->GetDirectory("HLT_MC/MCtruth");

  curdir->cd();
  
  TProfile2D *p2ref(0), *p2w24(0), *p2s23(0), *p2w24c(0), *p2s23c(0);
  p2ref = (TProfile2D*)dref->Get("p2r"); assert(p2ref);
  p2w24 = (TProfile2D*)dw24->Get("p2r_raw"); assert(p2w24);
  p2s23 = (TProfile2D*)ds23->Get("p2r_raw"); assert(p2s23);
  p2w24c = (TProfile2D*)dw24->Get("p2r"); assert(p2w24c);
  p2s23c = (TProfile2D*)ds23->Get("p2r"); assert(p2s23c);

  int i1ref = p2ref->GetXaxis()->FindBin(0.);
  int i2ref = p2ref->GetXaxis()->FindBin(1.305-0.05);
  TProfile *p1ref = p2ref->ProfileY("p1ref",i1ref,i2ref);
  
  int i1 = p2w24->GetXaxis()->FindBin(0.);
  int i2 = p2w24->GetXaxis()->FindBin(1.305-0.05);
  TProfile *p1w24 = p2w24->ProfileY("p1w24",i1,i2);
  TProfile *p1s23 = p2s23->ProfileY("p1s23",i1,i2);
  TProfile *p1w24c = p2w24c->ProfileY("p1w24c",i1,i2);
  TProfile *p1s23c = p2s23c->ProfileY("p1s23c",i1,i2);

  TH1D *h = tdrHist("h","Response",0.85+1e-4,1.1-1e-4);
  lumi_136TeV = "Winter24 vs Summer23";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,1,3500,1);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.35,0.86,"|#eta|<1.3");

  tdrDraw(p1ref,"HIST",kNone,kGreen+2,kSolid,-1,kNone);
  tdrDraw(p1s23c,"HIST",kNone,kBlack,kDashed,-1,kNone);
  tdrDraw(p1w24c,"HIST",kNone,kRed,kDashed,-1,kNone);
  tdrDraw(p1s23,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(p1w24,"HIST",kNone,kRed,kSolid,-1,kNone);

  TLegend *leg = tdrLeg(0.50,0.90-0.05*5,0.75,0.90);
  leg->AddEntry(p1ref,"QCD-MG-S23 (ref.)","L");
  leg->AddEntry(p1w24c,"Winter24 (nom.)","L");
  leg->AddEntry(p1s23c,"Summer23 (nom.)","L");
  leg->AddEntry(p1w24,"Winter24 (raw)","L");
  leg->AddEntry(p1s23,"Summer23 (raw)","L");
    
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawMCtruth/drawMCtruth_W24vsS23.pdf");
} // drawMCtruth
