// Purpose: Draw comparison of L3Res correction vs time for Run3
//          Use global fit results from jecdata*.root
//          (bonus: consider also using raw merged data from X+jet)
#include "TFile.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"
bool doRawData = false; // not yet implemented/tested

void drawL3ResVsTimes(string set);

void drawL3ResVsTime() {

  drawL3ResVsTimes("Prompt24");
  drawL3ResVsTimes("Prompt25");
  drawL3ResVsTimes("Prompt26");
  drawL3ResVsTimes("Prompt24to26");
}
  
void drawL3ResVsTimes(string set) {

  const char *cs = set.c_str();
  
  gROOT->ProcessLine(".! mkdir pdf/drawL3ResVsTime");
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! touch pdf/drawL3ResVsTime");
  gROOT->ProcessLine(".! touch pdf");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  double eps(1e-4);
  double ymin(0.93+eps), ymax(1.15-eps), ymind(0.95-eps), ymaxd(1.05+eps);
  if (set=="Prompt24") {
    ymin = 0.93+eps; ymax = 1.15-eps; ymind = 0.95-eps; ymaxd = 1.05+eps;
  }
  else if (set=="Prompt25") {
    ymin = 0.93+eps; ymax = 1.15-eps; ymind = 0.97-eps; ymaxd = 1.03+eps;
  }
  else if (set=="Prompt26") {
    ymin = 0.93+eps; ymax = 1.15-eps; ymind = 0.985-eps; ymaxd = 1.015+eps;
  }
  else if (set=="Prompt24to26") {
    ymin = 0.93+eps; ymax = 1.15-eps; ymind = 0.92-eps; ymaxd = 1.08+eps;
  }

  double xmin(15), xmax(4500);
  TH1D *hu = tdrHist(Form("hu_%s",cs),
		     "p_{T}-dependent JES",ymin,ymax,
		     "p_{T,ref} (GeV)",xmin,xmax);
  TH1D *hd = tdrHist(Form("hd_%s",cs),
		     "Ratio to 26C",ymind,ymaxd,
		     "p_{T,ref} (GeV)",xmin,xmax);
  
#include "../Config.C"
  if (set=="Prompt24")
    lumi_136TeV = Form("%s, %s","Prompt2024",mlum["2024"].c_str());
  else if (set=="Prompt25")
    lumi_136TeV = Form("%s, %s","Prompt2025",mlum["2025"].c_str());
  else if (set=="Prompt26")
    lumi_136TeV = Form("%s, >%s","Prompt2026",mlum["2026"].c_str());
  else if (set=="Prompt24to26")
    lumi_136TeV = Form("%s, >%s","2024 to 2026",mlum["24to26C"].c_str());
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cs),hu,hd,8,11);

  vector<string> vrun;
  if (set=="Prompt24to26") {
    //string arun[] = {"2026C","2026B","2025CDEFG",/*"2025C",*/
    //		     "2024FGHI_nib","2024F_nib1","2024CDE_nib"};
    string arun[] = {"2024CDE_nib","2024F_nib1","2024FGHI_nib","2025CDEFG",
		     "2026B","2026C","2026D"};

    const int nrun = sizeof(arun)/sizeof(arun[0]);
    vrun.resize(nrun);
    for (int i = 0; i != nrun; ++i) vrun[i] = arun[i];
  }
  else if (set=="Prompt26") {
    string arun[] = {"2026B","2026C","2026D"};
    const int nrun = sizeof(arun)/sizeof(arun[0]);
    vrun.resize(nrun);
    for (int i = 0; i != nrun; ++i) vrun[i] = arun[i];
  }
  else if (set=="Prompt25") {
    string arun[] = {"2025CDEFG","2025C","2025D","2025E","2025F","2025G"};
    const int nrun = sizeof(arun)/sizeof(arun[0]);
    vrun.resize(nrun);
    for (int i = 0; i != nrun; ++i) vrun[i] = arun[i];
  }
  else if (set=="Prompt24") {
    string arun[] = {"2024_nib","2024C_nib1","2024D_nib1","2024E_nib1",
		     "2024F_nib1","2024F_nib2","2024F_nib3",
		     "2024G_nib1","2024G_nib2","2024H_nib1","2024I_nib1"};
    const int nrun = sizeof(arun)/sizeof(arun[0]);
    vrun.resize(nrun);
    for (int i = 0; i != nrun; ++i) vrun[i] = arun[i];
  }
  else
    assert(false);
  const int nrun = int(vrun.size());
  
  map<string, int> mcolor;
  map<string, int> fcolor;
  if (set=="Prompt24to26") {
    mcolor["2024CDE_nib"] = kBlue;
    mcolor["2024F_nib1"] = kCyan+1;
    mcolor["2024FGHI_nib"] = kCyan+3;//kMagenta+2;
    mcolor["2025C"] = kSpring+2;
    mcolor["2025CDEFG"] = kGreen+2;
    mcolor["2026B"] = kOrange+2;
    fcolor["2026B"] = kOrange;
    mcolor["2026C"] = kRed;
    mcolor["2026D"] = kMagenta+2;
  }
  else if (set=="Prompt26") {
    mcolor["2026B"] = kBlack;
    mcolor["2026C"] = kRed;
    mcolor["2026D"] = kBlue;

  }
  else if (set=="Prompt25") {
    mcolor["2025CDEFG"] = kBlack;
    mcolor["2025C"] = kBlue;
    mcolor["2025D"] = kCyan+2;
    mcolor["2025E"] = kGreen+2;
    mcolor["2025F"] = kOrange+2;
    mcolor["2025G"] = kRed;
  }
  else if (set=="Prompt24") {
    mcolor["2024CDE_nib"] = kBlue;
    mcolor["2024C_nib1"] = kBlue;
    mcolor["2024D_nib1"] = kMagenta+1;
    mcolor["2024E_nib1"] = kMagenta+2;
    mcolor["2024F_nib1"] = kCyan+1;
    mcolor["2024F_nib2"] = kGreen+1;
    mcolor["2024F_nib3"] = kGreen+2;
    mcolor["2024G_nib1"] = kSpring+2;
    mcolor["2024G_nib2"] = kSpring+3;
    mcolor["2024H_nib1"] = kSpring+4;//kRed;
    mcolor["2024I_nib1"] = kBlack;
    mcolor["2024FGHI_nib"] = kMagenta+2;
    mcolor["2024_nib"] = kGray+2;

  }
  else
    assert(false);

  string ref("");
  map<string, const char*> mlabel;
  if (set=="Prompt24to26") {
    ref = "24CDE";//"26C";
    mlabel["2024_nib"] = "2024CDEFGHI";
    mlabel["2024CDE_nib"] = "24CDE (re-reco DI)";
    mlabel["2024F_nib1"] = "24F_nib1 (prompt DD)";
    mlabel["2024FGHI_nib"] = "24FGHI (prompt DI)";
    mlabel["2025C"] = "25C (prompt pre-timing)";
    mlabel["2025CDEFG"] = "25CDEFG (prompt DD)";
    mlabel["2026B"] = "26B (high PU)";
    mlabel["2026C"] = "26C (low PU)";
    mlabel["2026D"] = "26D (high PU)";
  }
  else if (set=="Prompt26") {
    ref = "26B";
    mlabel["2026B"] = "2026B (high PU)";
    mlabel["2026C"] = "2026C (low PU)";
    mlabel["2026D"] = "2026D (high PU)";
  }
  else if (set=="Prompt25") {
    ref = "2025";
    mlabel["2025C"] = "2025C (before HCAL time alignment)";
  }
  else if (set=="Prompt24") {
    ref = "2024";
    mlabel["2024_nib"] = "2024CDEFGHI";
    mlabel["2024C_nib1"] = "2024C (re-reco DI, HFv3)";
    mlabel["2024D_nib1"] = "2024D (re-reco)";
    mlabel["2024E_nib1"] = "2024E (re-reco)";
    mlabel["2024F_nib1"] = "2024F_nib1 (prompt DD)";
    mlabel["2024F_nib2"] = "2024F_nib2 (prompt DI)";
    mlabel["2024F_nib3"] = "2024F_nib3 (DI)";
    mlabel["2024G_nib1"] = "2024G_nib1 (DI)";
    mlabel["2024G_nib2"] = "2024G_nib2 (DI)";
    mlabel["2024H_nib1"] = "2024H (DI)";
    mlabel["2024I_nib1"] = "2024I (DI, HFv3)";
  }

  hd->SetYTitle(Form("Ratio to %s",ref!="" ? ref.c_str() : mlabel[vrun[0]] ? mlabel[vrun[0]] : vrun[0].c_str()));
  
  
  TLine *l = new TLine();
  l->SetLineColor(kGray+1);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  c1->cd(1);
  gPad->SetLogx();
  drawCustomLogXLabels(hu);
  
  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.03,xmax,1.03);
  l->DrawLine(xmin,0.95,xmax,0.95);

  TLegend *leg = tdrLeg(0.37,0.76-0.03*(nrun-0),0.62,0.76); // 2024
  leg->SetTextSize(0.035); // 2024
  
  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();
  drawCustomLogXLabels(hd);

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.01,xmax,1.01);
  l->DrawLine(xmin,0.99,xmax,0.99);
  l->SetLineStyle(kDashDotted);
  if (ymind<0.96) {
    l->DrawLine(xmin,1.05,xmax,1.05);
    l->DrawLine(xmin,0.95,xmax,0.95);
  }
  
  TH1D *h(0), *href(0);
  double pt1(0), pt2(0);
  for (int irun = 0; irun != nrun; ++irun) {

    string run = vrun[irun];
    const char *cr = run.c_str();

    string s = Form("rootfiles/jecdata%s.root",cr);
    TFile *f = new TFile(s.c_str(),"READ");
    if (!f || f->IsZombie()) cout << "Missing file "<<s<<endl<<flush;
    assert(f && !f->IsZombie());
    curdir->cd();
    
    TH1D *hl3 = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet"); assert(hl3);
    TH1D *hdt = (TH1D*)f->Get("ratio/eta00-13/hdm_cmb_mj_res"); assert(hdt);
    
    if (doRawData) h = (TH1D*)hdt->Clone(Form("h_%s_%s",cs,cr));
    else           h = (TH1D*)hl3->Clone(Form("h_%s_%s",cs,cr));

    if (!href) href = (TH1D*)h->Clone(Form("href_%s",cs));
    TH1D *hr = (TH1D*)h->Clone(Form("hr_%s_%s",cr,cs));
    hr->Divide(href);
    
    c1->cd(1);
    //tdrDraw(h,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);
    int fillcolor(mcolor[run]-9);
    if (mcolor[run]==kBlack) fillcolor = kGray+1;
    if (fcolor[run]!=0) fillcolor = fcolor[run];
    tdrDraw(h,"E3",kNone,mcolor[run],kSolid,-1,1001,fillcolor);
    h->SetFillColorAlpha(fillcolor,0.7);
    tdrDraw(new TGraph(h),"L",kNone,mcolor[run],kSolid,-1,kNone,0);

    if (mlabel[cr]) leg->AddEntry(h,mlabel[cr],"FLE");
    else            leg->AddEntry(h,cr,"FLE");
    
    c1->cd(2);
    //tdrDraw(hr,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);
    tdrDraw(hr,"E3",kNone,mcolor[run],kSolid,-1,1001,fillcolor);
    hr->SetFillColorAlpha(fillcolor,0.7);
    tdrDraw(new TGraph(hr),"L",kNone,mcolor[run],kSolid,-1,kNone,0);

  }

  c1->cd(1);
  gPad->RedrawAxis();

  if (doRawData)
    tex->DrawLatex(0.36,0.85,"Merged Z+jet, #gamma+jet and dijet data");
  else
    tex->DrawLatex(0.36,0.85,"Global fit of Z/#gamma+jet and multijet");
  tex->DrawLatex(0.36,0.80,"|#eta_{jet}| < 1.3");

  c1->cd(2);
  gPad->RedrawAxis();
  
  c1->SaveAs(Form("pdf/drawL3ResVsTime/drawL3ResVsTime_%s.pdf",cs));
} // drawL2ResVsTime
