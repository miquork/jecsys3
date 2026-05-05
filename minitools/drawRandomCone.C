// Purpose: draw comparison of Hirak's RandomCone results from various eras
#include "TFile.h"

#include <map>

#include "../tdrstyle_mod22.C"

void drawRandomConeStack();

void drawRandomCone() {

  drawRandomConeStack(); exit(0);

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/drawRandomCone");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/drawRandomCone");
  
  TFile *f24 = new TFile("rootfiles/Hirak_compare_sf_Run2024_Summer24MC-merge3.root","READ");
  assert(f24 && !f24->IsZombie());
  
  TFile *f25 = new TFile("rootfiles/Hirak_2025CDEFG/compare_sf_Run2025-Summer24MC.root","READ");
  assert(f25 && !f25->IsZombie());

  curdir->cd();

  double miny(0.5+1e-4), maxy(2.0-1e-4);
  TH1D *h = tdrHist("h","Offset Scale Factor",miny,maxy,"#eta",-5.2,5.2);
  //lumi_136TeV = "2024-2025";
  lumi_136TeV = "2025, 109 fb^{-1}";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);//kRectangular);

  // 2025 vs 2024
  /*
  string vera[] = {
    "2024C","2024D","2024E-v1","2024C-reReco",
    "2024F-nib2","2024H","2024I",
    "2025Cv1","2025D","2025G"};
  */
  // 2025 detailed breakdown
  string vera[] = {
    "25Cv1","25Cv2","25D","25E","25Fv1","25Fv2","25G"};
  const int nera = sizeof(vera)/sizeof(vera[0]);

  map<string, int> color;
  color["2024C"] = kBlue;
  color["2024D"] = kBlue-9;
  color["2024E-v1"] = kGreen+2;//kCyan+1;
  color["2024E-reReco"] = kGreen+2;
  color["2024C-reReco"] = kCyan+1;
  color["2024F-nib2"] = kOrange+1;
  color["2024H"] = kOrange+2;
  color["2024I"] = kBlack;
  color["2025Cv1"] = kRed-9;
  color["2025D"] = kRed;
  color["2025G"] = kRed+1;
  //
  color["25Cv1"] = kBlue;
  color["25Cv2"] = kBlue-9;
  color["25D"] = kCyan+1;
  color["25E"] = kGreen+2;
  color["25Fv1"] = kOrange+1;
  color["25Fv2"] = kOrange+2;
  color["25G"] = kRed;

  map<string, const char*> label;
  label["2024C"] = "24C (start'24)";
  label["2024D"] = "24D";
  label["2024E-v1"] = "24E (1st HF)";
  label["2024E-reReco"] = "24E (CDE re-reco)";
  label["2024C-reReco"] = "24C (CDE re-reco)";
  label["2024F-nib2"] = "24F (HCAL DI)";
  label["2024H"] = "24H";
  label["2024I"] = "24I (2nd HF)";
  label["2025Cv1"] = "25C (start'25)";
  label["2025D"] = "25D (HCAL timing)";
  label["2025G"] = "25G (end'25)";
  //
  label["25Cv1"] = "25Cv1 (start'25)";
  label["25Cv2"] = "25Cv2";
  label["25D"] = "25D (HCAL timing)";
  label["25Fv1"] = "25Fv1";
  label["25Fv2"] = "25Fv2";
  label["25G"] = "25G (end'25)";
  
  TLegend *leg = tdrLeg(0.30,0.90-0.035*7,0.50,0.90);
  leg->SetTextSize(0.035);
  TLegend *leg2 = tdrLeg(0.60,0.90-0.035*3,0.80,0.90);
  leg2->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.65,0.70,"#sigma_{MB} = 75.3 mb");
  tex->DrawLatex(0.65,0.65,"Summer24 MC");

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);
  l->DrawLine(-5.2,1,5.2,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(-5.2,1.2,5.2,1.2);
  l->DrawLine(-5.2,0.8,5.2,0.8);
  
  for (int i = 0; i != nera; ++i) {
    string se = vera[i];
    const char *ce = se.c_str();
    TH1D *hrc(0);
    hrc = (TH1D*)f24->Get(Form("Run%s",ce));
    if (!hrc) hrc = (TH1D*)f25->Get(Form("Run%s",ce));
    if (!hrc) hrc = (TH1D*)f25->Get(Form("Run20%s",ce));
    if (!hrc) { cout << "Missing Run" << se << endl << flush; exit(0); }

    if (hrc->GetBinContent(42)>2) hrc->Scale(0.5); // merge3 duplicates
    hrc->Scale(69.3/75.3);
    
    tdrDraw(hrc,"HIST",kNone,color[se],kSolid,-1,kNone,0);
    if (i<7) leg->AddEntry(hrc,label[se],"L");
    if (i>=7) leg2->AddEntry(hrc,label[se],"L");
  }

  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawRandomCone/drawRandomCone.pdf");
} // void drawRandomCone


void drawRandomConeStack() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f = new TFile("rootfiles/Hirak_2025CDEFG_stack/stack_mikko_Run3-Run2025CDEFG-Summer24_DataMC_R4_all_nPU60.root","READ");
  assert(f && !f->IsZombie());

  curdir->cd();

  double miny(0.+1e-4), maxy(0.65-1e-4), minyd(-0.55+1e-4), maxyd(0.55-1e-4);
  TH1D *h = tdrHist("hs","p_{T} offset / #mu",miny,maxy,"#eta",-5.2,5.2);
  //TH1D *hd = tdrHist("hds","Data - MC",minyd,maxyd,"#eta",-5.2,5.2);
  TH1D *hd = tdrHist("hds","(D-M)/totM",minyd,maxyd,"#eta",-5.2,5.2);
  lumi_136TeV = "Summer24 MC + 2025, 109 fb^{-1}";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas("c1s",h,hd,8,11);

  string vs[] = {"Assoc_CH","Unassoc_CH","Photons","Neutral_Had","HF_EM","HF_Had"};
  const int ns = sizeof(vs)/sizeof(vs[0]);

  map<string, int> mcolor;
  mcolor["Assoc_CH"] = kRed;
  mcolor["Unassoc_CH"] = kRed-9;
  mcolor["Photons"] = kBlue;
  mcolor["Neutral_Had"] = kGreen+2;
  mcolor["HF_EM"] = kCyan+1;
  mcolor["HF_Had"] = kMagenta+1;
  mcolor["S"] = kBlack;
  
  map<string, int> mmarker;
  mmarker["Assoc_CH"] = kFullSquare;
  mmarker["Unassoc_CH"] = kOpenSquare;
  mmarker["Photons"] = kFullCircle;
  mmarker["Neutral_Had"] = kOpenDiamond;
  mmarker["HF_EM"] = kFullStar;
  mmarker["HF_Had"] = kOpenStar;
  mmarker["S"] = kFullDiamond;

  //map<string, const char*> label;
  //label["25G"] = "25G (end'25)";

  c1->cd(1);
  
  TLegend *leg = tdrLeg(0.40,0.90-0.045*ns,0.65,0.90);
  //leg->SetTextSize(0.035);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.40,0.54,"#sigma_{MB} = 75.3 mb");
  //tex->DrawLatex(0.40,0.55,"Summer24 MC");

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);
  //l->DrawLine(-5.2,1,5.2,1);
  //l->SetLineStyle(kDotted);
  //l->DrawLine(-5.2,1.2,5.2,1.2);
  //l->DrawLine(-5.2,0.8,5.2,0.8);

  TH1D *hsumd(0), *hsumm(0);
  vector<TH1D*> vh(ns);
  for (int i = 0; i != ns; ++i) {
    string ss = vs[i];
    const char *cs = ss.c_str();
    TH1D *h(0), *hm(0), *hr;
    cout << "Reading " << cs << endl << flush;
    hd = (TH1D*)f->Get(Form("%s_Data",cs)); assert(hd);
    hm = (TH1D*)f->Get(Form("%s_MC",cs));   assert(hm);
    hr = (TH1D*)hd->Clone(Form("%s_diff",cs)); assert(hr);
    hr->Add(hm,-1);

    if (hsumd) hsumd->Add(hd);
    else       hsumd = (TH1D*)hd->Clone("hsumd");
    if (hsumm) hsumm->Add(hm);
    else       hsumm = (TH1D*)hm->Clone("hsumm");
    
    c1->cd(1);
    tdrDraw(hm,"HIST",kNone,mcolor[ss],kSolid,-1,1001,mcolor[ss],0.6,1,0.2);
    //hm->SetFillColorAlpha(mcolor[ss]-9,0.3);
    tdrDraw(hd,"Pz",mmarker[ss],mcolor[ss],kSolid,-1,1001,mcolor[ss],0.6,1,0.2);
	
    leg->AddEntry(hd,cs,"PFLE");

    //c1->cd(2);
    //tdrDraw(hr,"Pz",mmarker[ss],mcolor[ss],kSolid,-1,kNone,0,0.8);
    
    vh[i] = hr;
  }

  c1->cd(1);
  tdrDraw(hsumm,"HIST",kNone,mcolor["S"],kSolid,-1,kNone,mcolor["S"],0.6,1,0.2);
  tdrDraw(hsumd,"Pz",mmarker["S"],mcolor["S"],kSolid,-1,1001,kWhite,0.6,1,0.2);
  
  leg->SetY1(leg->GetY1()-0.05);
  leg->AddEntry(hsumd,"Total","PFLE");
  
  c1->cd(2);
  TH1D *hsumr = (TH1D*)hsumd->Clone("hsumr");
  hsumr->Add(hsumm,-1);
  //tdrDraw(hsumr,"Pz",mmarker["S"],mcolor["S"],kSolid,-1,kNone,0,0.8);

  for (int i = 0; i != ns; ++i) {
    string ss = vs[i];
    TH1D *hr = vh[i];
    hr->Divide(hsumm);
    tdrDraw(hr,"Pz",mmarker[ss],mcolor[ss],kSolid,-1,kNone,0,0.8);
  }
  hsumr->Divide(hsumm);
  tdrDraw(hsumr,"Pz",mmarker["S"],mcolor["S"],kSolid,-1,kNone,0,0.8);
  
  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawRandomCone/drawRandomConeStack.pdf");
} // void drawRandomCone
