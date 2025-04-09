// Purpose: Draw HE scale variations to data from Fikri and Sami
#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

// Minus 1, times 100
void shiftAndScale(TH1D *h, double dy = 1, double k = 100);

void drawHEscale_Zmm() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  string vf[] = {"x1","Prompt","MC","x1p5a","x2p0a","x1p5b","x2p0b"};
  const int nf = sizeof(vf)/sizeof(vf[0]);
  int color[] = {kBlack,kGray+1,kBlack,kGreen+2,kOrange+1,kBlue,kRed};
  int marker[] = {kFullDiamond,kOpenDiamond,kNone,kOpenDiamond,kOpenDiamond,kFullDiamond,kFullDiamond};
  map<string,string> mleg;
  mleg["Prompt"] = "Prompt24";
  mleg["x1"] = "Prompt25";
  mleg["MC"] = "Summer24 MC";
  mleg["x1p5a"] = "#times 1.5 (26-29)";
  mleg["x2p0a"] = "#times 2.0 (26-29)";
  mleg["x1p5b"] = "#times 1.5 (28-29)";
  mleg["x2p0b"] = "#times 2.0 (28-29)";
  
  //TH1D *h = tdrHist("h","MPF-1 (%)",-50,20,
  //TH1D *h = tdrHist("h","MPF",0.9,1.5, // Corrected MPF
  TH1D *h = tdrHist("h","MPF",0.65,1.7, // Raw-level MPF
  //TH1D *h = tdrHist("h","MPF",0.95,1.15, // Raw-level MPF over base
  //TH1D *h = tdrHist("h","MPF",0.90,1.46, // Raw-level MPF over base, with MC
		    "|#eta_{jet}|",0,5.2);
  extraText = "Private";
  lumi_136TeV = "HE scale variations, 2024I Z(#mu#mu)+jet";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0,1,5.2,1);
  
  //tex->DrawLatex(0.18,0.75,"|#eta| < 1.3");
  tex->DrawLatex(0.18,0.75,"28<p_{T,Z}<49 GeV");
  
  TLegend *leg = tdrLeg(0.60,0.90-0.045*6,0.85,0.90);
  //leg->SetNColumns(2);

  TH1D *hbase(0);
  for (int i = 0; i != nf; ++i) {

    string sf = vf[i];
    const char *cf = sf.c_str();
    TFile *f(0);
    if (sf=="Prompt" || sf=="MC")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_nib1_Zmm_NoPU_V8M_v96.root","READ");
    else if (sf=="x1")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RAWToPFNANO_V8M_v96.root","READ");
    else if (sf=="x1p5a")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RHAbsIEta26to29D1x1p5_RAWToPFNANO_V8M_v96.root","READ");
    else if (sf=="x2p0a")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RHAbsIEta26to29D1x2p0_RAWToPFNANO_V8M_v96.root","READ");
    else if (sf=="x1p5b")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RHAbsIEta28to29D1x1p5_RAWToPFNANO_V8M_v96.root","READ");
    else if (sf=="x2p0b")
      f = new TFile("rootfiles/Prompt2024/v96_Zmm/jme_Zj_2024I_UpdatedGainsRespCorrs_RHAbsIEta28to29D1x2p0_RAWToPFNANO_V8M_v96.root","READ");
    assert(f && !f->IsZombie());

    if (sf=="MC") f->cd("mc");
    else f->cd("data");
    TDirectory *d = gDirectory;
    
    TProfile2D *p2 = (TProfile2D*)d->Get("l2res/p2m0tc"); assert(p2);
    TProfile2D *p2res = (TProfile2D*)d->Get("l2res/p2res"); assert(p2res);
    int j1 = p2->GetYaxis()->FindBin(28.);
    int j2 = p2->GetYaxis()->FindBin(49.)-1;
    TProfile *p = p2->ProfileX(Form("p_%s",cf),j1,j2);
    TProfile *pres = p2res->ProfileX(Form("pres_%s",cf),j1,j2);
    TH1D *h = p->ProjectionX(Form("h_%s",cf));
    TH1D *hres = pres->ProjectionX(Form("hres_%s",cf));
    h->Multiply(hres);

    if (i==0) hbase = (TH1D*)h->Clone("hbase");
    //h = (TH1D*)h->Clone(Form("h_%s",cf));
    assert(hbase);
    //h->Divide(h,hbase,1,1,"B");
    //shiftAndScale(h);

    tdrDraw(h,sf=="MC" ? "HIST" : "Pz",marker[i],color[i],kSolid,-1,kNone,0);
    leg->AddEntry(h,mleg[cf].c_str(),sf=="MC" ? "F" : "PLE");

    //tdrDraw(p,"Pz",marker[i],color[i],kSolid,-1,kNone,0);
    //leg->AddEntry(p,mleg[cf].c_str(),"PLE");

  } // for i


  c1->cd();
  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawHEscale/drawHEscale_Zmm_MPF.pdf");
  
} // drawHEscale_Zmm


void shiftAndScale(TH1D *h, double dy, double k) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinContent(i)!=0) {
      h->SetBinContent(i, k*(h->GetBinContent(i)-dy));
      h->SetBinError(i, k*h->GetBinError(i));
    }
  } // for i
} // shiftAndScale
