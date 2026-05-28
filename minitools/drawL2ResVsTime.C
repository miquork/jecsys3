// Purpose: Draw comparison of L2Res correction vs time for Run3
//          Use raw merged data from X+jet
#include "TFile.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"
bool doRawData = false;

void drawL2ResVsTimes(string set);

void drawL2ResVsTime() {

  drawL2ResVsTimes("Prompt24");
  drawL2ResVsTimes("Prompt25");
  drawL2ResVsTimes("Prompt26");
  drawL2ResVsTimes("Prompt24to26");
}
  
void drawL2ResVsTimes(string set) {

  const char *cs = set.c_str();
  
  gROOT->ProcessLine(".! mkdir pdf/drawL2ResVsTime");
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! touch pdf/drawL2ResVsTime");
  gROOT->ProcessLine(".! touch pdf");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f(0);
  if (set=="Prompt26" || set=="Prompt25" || set=="Prompt24" ||
      set=="Prompt24to26") {
    //f = new TFile("rootfiles/L2Res_20250615_v2.root");
    f = new TFile("rootfiles/L2Res_Prompt24to26_V11MV5MV2M_draft.root");
  }
  else
    f = new TFile("rootfiles/L2Res.root");
  assert(f && !f->IsZombie());

  double eps(1e-4);
  //double ymin(0.8+1e-4), ymax(1.2-1e-4), ymind(0.945+1e-4), ymaxd(1.055-1e-4);
  //double ymin(0.70+1e-4), ymax(1.35-1e-4), ymind(0.73+1e-4), ymaxd(1.22-1e-4);
  double ymin(0.70+1e-4), ymax(1.45-1e-4), ymind(0.85+1e-4), ymaxd(1.15-1e-4);
  if (set=="Prompt24") {
    ymin = 0.70+eps; ymax = 1.45-eps; ymind = 0.85+eps; ymaxd = 1.15-eps;
  }
  else if (set=="Prompt25") {
    ymin = 0.74+eps; ymax = 1.20-eps; ymind = 0.95-eps; ymaxd = 1.05+eps;
  }
  else if (set=="Prompt26") {
    ymin = 0.64+eps; ymax = 1.20-eps; ymind = 0.85-eps; ymaxd = 1.10+eps;
  }
  else if (set=="Prompt24to26") {
    ymin = 0.70+eps; ymax = 1.35-eps; ymind = 0.73+eps; ymaxd = 1.22-eps;
  }
  TH1D *h = tdrHist(Form("h_%s",cs),
		    "#eta-dependent JES",ymin,ymax,"#eta",-5.2,5.2);
  TH1D *hd = tdrHist(Form("hd_%s",cs),
		     "Ratio to 26C",ymind,ymaxd,"#eta",-5.2,5.2);
  
#include "../Config.C"
  //lumi_136TeV = Form("%s, %s","2025",mlum["2025"].c_str());
  if (set=="Prompt24")
    lumi_136TeV = Form("%s, %s","Prompt2024",mlum["2024"].c_str());
  else if (set=="Prompt25")
    lumi_136TeV = Form("%s, %s","Prompt2025",mlum["2025"].c_str());
  else if (set=="Prompt26")
    lumi_136TeV = Form("%s, >%s","Prompt2026",mlum["2026"].c_str());
  else if (set=="Prompt24to26")
    lumi_136TeV = Form("%s, >%s","2024 to 2026",mlum["24to26C"].c_str());
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cs),h,hd,8,11);

  vector<string> vrun;
  if (set=="Prompt24to26") {
    string arun[] = {"2026C","2026B","2025CDEFG",/*"2025C",*/
		     "2024FGHI_nib","2024F_nib1","2024CDE_nib"};
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
  //string vrun[] = {"2025CDEFG","2025C","2025D","2025E","2025F","2025G"};
  //const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  const int nrun = int(vrun.size());
  
  map<string, int> mcolor;
  if (set=="Prompt24to26") {
    mcolor["2024CDE_nib"] = kBlue;
    mcolor["2024F_nib1"] = kCyan+1;
    mcolor["2024FGHI_nib"] = kMagenta+2;
    mcolor["2025C"] = kSpring+2;
    mcolor["2025CDEFG"] = kGreen+2;
    mcolor["2026B"] = kRed;
    mcolor["2026C"] = kBlack;
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
  /*
  mcolor["2025CDEFG"] = kBlack;
  mcolor["2025DEFG"] = kGray+2;
  mcolor["2025C"] = kBlue;
  mcolor["2025D"] = kCyan+2;
  mcolor["2025E"] = kGreen+2;
  mcolor["2025F"] = kOrange+1;
  mcolor["2025G"] = kRed;
  */
  string ref("");
  map<string, const char*> mlabel;
  if (set=="Prompt24to26") {
    ref = "26C";
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

  l->SetLineStyle(kDashed);
  l->DrawLine(-5.2,1,+5.2,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(-5.2,1.05,+5.2,1.05);
  //l->DrawLine(-5.2,0.95,+5.2,0.95);
  l->DrawLine(-5.2,0.90,+5.2,0.90);
  l->DrawLine(-1.305,ymin,-1.305,ymax);
  l->DrawLine(+1.305,ymin,+1.305,ymax);
  l->DrawLine(-2.500,ymin,-2.500,ymax);
  l->DrawLine(+2.500,ymin,+2.500,ymax);
  l->DrawLine(-2.964,ymin,-2.964,ymax);
  l->DrawLine(+2.964,ymin,+2.964,ymax);

  //TLegend *leg = tdrLeg(0.37,0.75-0.04*(nrun-0),0.62,0.75);
  TLegend *leg = tdrLeg(0.37,0.76-0.03*(nrun-0),0.62,0.76); // 2024
  leg->SetTextSize(0.035); // 2024
  
  gPad->RedrawAxis();
  
  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->DrawLine(-5.2,1,+5.2,1);
  l->SetLineStyle(kDotted);
  //l->DrawLine(-5.2,1.01,+5.2,1.01);
  //l->DrawLine(-5.2,0.99,+5.2,0.99);
  l->DrawLine(-5.2,1.10,+5.2,1.10);
  l->DrawLine(-5.2,0.90,+5.2,0.90);
  l->DrawLine(-1.305,ymind,-1.305,ymaxd);
  l->DrawLine(+1.305,ymind,+1.305,ymaxd);
  l->DrawLine(-2.500,ymind,-2.500,ymaxd);
  l->DrawLine(+2.500,ymind,+2.500,ymaxd);
  l->DrawLine(-2.964,ymind,-2.964,ymaxd);
  l->DrawLine(+2.964,ymind,+2.964,ymaxd);
  
  TH1D *href(0);
  double pt1(0), pt2(0);
  for (int irun = 0; irun != nrun; ++irun) {

    string run = vrun[irun];
    const char *cr = run.c_str();

    //TH2D *h2m = (TH2D*)f->Get(Form("hmerged1_%s",cr)); assert(h2m);
    //TH2D *h2m = (TH2D*)f->Get(Form("h2merged1_%s",cr)); assert(h2m);
    TH2D *h2m = (TH2D*)f->Get(Form("h2jes1_%s",cr)); assert(h2m);
    double pt = 100;
    int i = h2m->GetYaxis()->FindBin(pt);
    pt1 = h2m->GetYaxis()->GetBinLowEdge(i);
    pt2 = h2m->GetYaxis()->GetBinLowEdge(i+1);
    TH1D *hm = h2m->ProjectionX(Form("hm_%s_%s",cr,cs),i,i);

    if (!href) href = (TH1D*)hm->Clone("href");
    TH1D *hr = (TH1D*)hm->Clone(Form("hr_%s_%s",cr,cs));
    hr->Divide(href);
    
    c1->cd(1);
    tdrDraw(hm,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);

    //if (run!="2025CDEFG") leg->AddEntry(hm,cr,"LE");
    //if (run!="2026C") leg->AddEntry(hm,cr,"LE");
    //if (run=="2026C") leg->AddEntry(hm,"2026C (low PU)","LE");
    if (mlabel[cr]) leg->AddEntry(hm,mlabel[cr],"LE");
    else            leg->AddEntry(hm,cr,"LE");
    
    c1->cd(2);
    tdrDraw(hr,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);
  }

  c1->cd(1);
  gPad->RedrawAxis();

  if (doRawData)
    tex->DrawLatex(0.36,0.85,"Merged Z+jet, #gamma+jet and dijet data");
  else
    tex->DrawLatex(0.36,0.85,"Fit Z/#gamma+jet, dijet and Random Cone");
  tex->DrawLatex(0.36,0.80,Form("%1.0f < p_{T,ref} < %1.0f GeV",pt1,pt2));

  c1->cd(2);
  gPad->RedrawAxis();
  
  c1->SaveAs(Form("pdf/drawL2ResVsTime/drawL2ResVsTime_%s.pdf",cs));
} // drawL2ResVsTime
