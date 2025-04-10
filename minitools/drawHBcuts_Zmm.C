// Purpose: Draw HE scale variations to data from Fikri and Sami
#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

// Minus 1, times 100
void shiftAndScale(TH1D *h, double dy = 1, double k = 100);

void drawHBcuts_Zmm_sets(string set);
void drawHBcuts_Zmm() {
  drawHBcuts_Zmm_sets("MPF");
  drawHBcuts_Zmm_sets("DB");
  drawHBcuts_Zmm_sets("NHF");
} // drawHBcuts_Zmm

void drawHBcuts_Zmm_sets(string set="MPF") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cs = set.c_str();
  
  string vf[] = {"Base","MC","Prompt25","1x","1p5x","2x","4x"};
  const int nf = sizeof(vf)/sizeof(vf[0]);

  map<string,const char*> mleg;
  mleg["Base"] = "";
  mleg["Prompt25"] = "Prompt25 (DD)";
  mleg["MC"] = "Summer24 MC";
  mleg["1x"] = "Prompt24 (DI)";
  mleg["1p5x"] = "#times 1.5 (DI)";
  mleg["2x"] = "#times 2    (DI)";
  mleg["4x"] = "#times 4    (DI)";

  map<string, int> mcol;
  mcol["Base"] = kNone;
  mcol["Prompt25"] = kGray+1;
  mcol["MC"] = kBlack;
  mcol["1x"] = kBlack;
  mcol["1p5x"] = kGreen+2;
  mcol["2x"] = kBlue;
  mcol["4x"] = kRed;

  map<string, int> mmar;
  mmar["Base"] = 0;
  mmar["Prompt25"] = kNone;
  mmar["MC"] = kNone;
  mmar["1x"] = kFullDiamond;
  mmar["1p5x"] = kFullDiamond;
  mmar["2x"] = kFullDiamond;
  mmar["4x"] = kFullDiamond;
  
  TH1D *h = tdrHist(Form("h_%s",cs),Form("%s / L2L3Res",cs),0.83,1.57,
		    "|#eta_{jet}|",0,5.2);
  extraText = "Private";
  lumi_136TeV = "HB cut variations, 2024I Z(#mu#mu)+jet";
  TCanvas *c1 = tdrCanvas(Form("c1_%s",cs),h,8,11,kSquare);

  TH1D *h_2 = tdrHist(Form("h_2_%s",cs),Form("%s / L2L3Res",cs),0.87,1.12,
		      "p_{T,Z} (GeV)",15,500);
  extraText = "Private";
  lumi_136TeV = "HB cut variations, 2024I Z(#mu#mu)+jet";
  TCanvas *c2 = tdrCanvas(Form("c2_%s",cs),h_2,8,11,kSquare);

  TH1D *h_3 = tdrHist(Form("h_3_%s",cs),Form("%s (100/Base)",cs),0.001,100.,
		      "NHF",0,1);
  extraText = "Private";
  lumi_136TeV = "HB cut variations, 2024I Z(#mu#mu)+jet";
  TCanvas *c3 = tdrCanvas(Form("c3_%s",cs),h_3,8,11,kSquare);

  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  c1->cd();
  l->DrawLine(0,1,5.2,1);
  l->DrawLine(0.522,0.83,0.522,1.35);

  if (set=="MPF") tex->DrawLatex(0.18,0.75,"15<p_{T,Z}<28 GeV");
  if (set=="NHF") tex->DrawLatex(0.18,0.75,"15<p_{T,Z}<28 GeV");
  if (set=="DB")  tex->DrawLatex(0.18,0.75,"28<p_{T,Z}<49 GeV");
  //if (set=="DB")  tex->DrawLatex(0.18,0.75,"28<p_{T,Z}<37 GeV");
  
  TLegend *leg = tdrLeg(0.60,0.90-0.045*6,0.85,0.90);

  c2->cd();
  gPad->SetLogx();
  l->DrawLine(15,1,500,1);  
  if (set=="MPF") {
    l->DrawLine(15,0.87,15,1.05);
    l->DrawLine(28,0.87,28,1.05);
  }
  if (set=="NHF") {
    l->DrawLine(15,0.87,15,1.05);
    l->DrawLine(28,0.87,28,1.05);
  }
  if (set=="DB") {
    l->DrawLine(28,0.87,28,1.05);
    l->DrawLine(37,0.87,37,1.05);
  }

  //tex->DrawLatex(0.18,0.75,"|#eta_{jet}|<1.305");
  tex->DrawLatex(0.18,0.75,"|#eta_{jet}|<0.522");

  TLegend *leg2 = tdrLeg(0.55,0.90-0.045*6,0.80,0.90);

  c3->cd();
  gPad->SetLogy();
  tex->DrawLatex(0.25,0.34,"|#eta_{jet}|<1.305");
  tex->DrawLatex(0.25,0.27,"28<p_{T,Z}<49 GeV");
  //tex->DrawLatex(0.25,0.20,"NHF = p_{T,nh} / p_{T,Z}"); 
  TLegend *leg3 = tdrLeg(0.55,0.90-0.045*6,0.80,0.90);
  
  
  TH1D *hbase(0), *hbase_pt(0), *hbase_f(0);
  for (int i = 0; i != nf; ++i) {

    string sf = vf[i];
    const char *cf = sf.c_str();
    TFile *f(0);
    if (sf=="Base")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_nib1_Zmm_NoPU_V8M_v96.root","READ");
    else if (sf=="1x" || sf=="MC")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_nib1_Zmm_NoPU_V8M_v96.root","READ");
    else if (sf=="Prompt25")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RAWToPFNANO_V8M_v96.root","READ");
    else 
      f = new TFile(Form("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_HcalPFCutsHB%s_RAWToPFNANO_V8M_v96.root",cf),"READ");
    assert(f && !f->IsZombie());

    if (sf=="MC") f->cd("mc");
    else f->cd("data");
    TDirectory *d = gDirectory;

    TProfile2D *p2(0);
    if (set=="MPF") { p2 = (TProfile2D*)d->Get("l2res/p2m0tc"); assert(p2); }
    if (set=="NHF") { p2 = (TProfile2D*)d->Get("l2res/p2m0tc"); assert(p2); }
    if (set=="DB")  { p2 = (TProfile2D*)d->Get("l2res/p2m2tc"); assert(p2); }
    assert(p2);
    TProfile2D *p2res = (TProfile2D*)d->Get("l2res/p2res"); assert(p2res);

    TH2D *h2f(0);
    if (set=="NHF") { h2f = (TH2D*)d->Get("eta_00_13/h_Zpt_neHEF_alpha100"); assert(h2f); }
    
    c1->cd();

    int j1(-1), j2(-1);
    if (set=="MPF" || set=="NHF") {
      j1 = p2->GetYaxis()->FindBin(15.);
      j2 = p2->GetYaxis()->FindBin(28.)-1;
    }
    if (set=="DB") {
      j1 = p2->GetYaxis()->FindBin(28.);
      j2 = p2->GetYaxis()->FindBin(49.)-1;
      //j2 = p2->GetYaxis()->FindBin(37.)-1;
    }
    assert(j1!=-1 && j2!=-1);
    
    TProfile *p = p2->ProfileX(Form("p_%s_%s",cs,cf),j1,j2);
    TProfile *pres = p2res->ProfileX(Form("pres_%s_%s",cs,cf),j1,j2);
    TH1D *h = p->ProjectionX(Form("h_%s_%s",cs,cf));
    TH1D *hres = pres->ProjectionX(Form("hres_%s_%s",cs,cf));

    // Undo JEC residuals (only leaves MC JEC)
    h->Multiply(hres);

    if (i==0) cout << Form("pt=[%1.0f,%1.0f] GeV\n",
			   p2->GetYaxis()->GetBinLowEdge(j1),
			   p2->GetYaxis()->GetBinLowEdge(j2+1));
    if (i==0) hbase = (TH1D*)h->Clone(Form("hbase_%s",cs)); 
    
    // Divide by baseline scale
    assert(hbase);
    h->Divide(h,hbase,1,1);//,"B");
    
    if (i!=0) {
      tdrDraw(h,mmar[sf]==kNone ? "HIST" : "Pz",mmar[sf],mcol[sf],kSolid,-1,kNone,0);
      leg->AddEntry(h,mleg[sf],mmar[sf]==kNone ? "F" : "PLE");
    }

    
    c2->cd();

    int i1 = p2->GetXaxis()->FindBin(0.);
    //int i2 = p2->GetXaxis()->FindBin(1.305)-1;
    int i2 = p2->GetXaxis()->FindBin(0.522)-1;
    TProfile *p_pt = p2->ProfileY(Form("p_pt_%s_%s",cs,cf),i1,i2);
    TProfile *pres_pt = p2res->ProfileY(Form("pres_pt_%s_%s",cs,cf),i1,i2);
    TH1D *h_pt = p_pt->ProjectionX(Form("h_pt_%s_%s",cs,cf));
    TH1D *hres_pt = pres_pt->ProjectionX(Form("hres_pt_%s_%s",cs,cf));

    // Undo JEC residuals (only leaves MC JEC)
    h_pt->Multiply(hres_pt);

    if (i==0) cout << Form("|eta|=[%1.3f,%1.3f]\n",
			   p2->GetXaxis()->GetBinLowEdge(i1),
			   p2->GetXaxis()->GetBinLowEdge(i2+1));
    if (i==0) hbase_pt = (TH1D*)h_pt->Clone(Form("hbase_pt_%s",cs));

    // Divide by baseline scale
    assert(hbase_pt);
    h_pt->Divide(h_pt,hbase_pt,1,1);//,"B");

    if (i!=0) {
      tdrDraw(h_pt,mmar[sf]==kNone ? "HIST" : "Pz",mmar[sf],mcol[sf],kSolid,-1,kNone,0);
      leg2->AddEntry(h_pt,mleg[cf],mmar[sf]==kNone ? "F" : "PLE");
    }


    if (set!="NHF") continue;
    
    c3->cd();
    int k1 = h2f->GetXaxis()->FindBin(28.);
    int k2 = h2f->GetXaxis()->FindBin(49.)-1;
    TH1D *h_f = h2f->ProjectionY(Form("hf_%s_%s",cs,cf),k1,k2);

    if (i==0) cout << Form("NHF p_{T,Z}=[%1.3f,%1.3f] (GeV)\n",
			   h2f->GetXaxis()->GetBinLowEdge(k1),
			   h2f->GetXaxis()->GetBinLowEdge(k2+1));
    if (i==0) hbase_f = (TH1D*)h_f->Clone(Form("hbase_f_%s",cs));

    // Divide by baseline scale
    assert(hbase_f);
    //h_f->Divide(h_f,hbase_f,1,1);//,"B");
    h_f->Rebin(2); h_f->Scale(0.5);
    //h_f->Rebin(4); h_f->Scale(0.25);
    //if (sf!="MC") h_f->Scale(100./hbase_f->Integral());
    //if (sf=="MC") h_f->Scale(100./h_f->Integral());
    h_f->Scale(100./hbase_f->Integral());
    if (sf=="MC") h_f->Scale(hbase_f->GetEntries()/h_f->GetEntries());
    
    if (i!=0) {
      tdrDraw(h_f,mmar[sf]==kNone ? "HIST" : "Pz",mmar[sf],mcol[sf],kSolid,-1,kNone,0);
      leg3->AddEntry(h_f,mleg[cf],mmar[sf]==kNone ? "F" : "PLE");
    }
  } // for i


  if (set=="MPF" || set=="DB") {
    c1->cd();
    gPad->RedrawAxis();
    c1->SaveAs(Form("pdf/drawHBcuts/drawHBcuts_Zmm_%s_vsEta.pdf",cs));
    
    c2->cd();
    gPad->RedrawAxis();
    c2->SaveAs(Form("pdf/drawHBcuts/drawHBcuts_Zmm_%s_vsPt.pdf",cs));
  }
  if (set=="NHF") {
    c3->cd();
    gPad->RedrawAxis();
    c3->SaveAs(Form("pdf/drawHBcuts/drawHBcuts_Zmm_%s_vsNHF.pdf",cs));
  }
    
} // drawHEscale_Zmm


void shiftAndScale(TH1D *h, double dy, double k) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinContent(i)!=0) {
      h->SetBinContent(i, k*(h->GetBinContent(i)-dy));
      h->SetBinError(i, k*h->GetBinError(i));
    }
  } // for i
} // shiftAndScale
