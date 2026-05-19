// Purpose: Draw comparison of JER SF vs time for Run3
//          Use global fit results from JERSF.root
//          (bonus: consider also using raw merged data from X+jet)
#include "TFile.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"
bool doRawData = false; // not yet implemented/tested

void drawJERSFvsTimes(string set);

void drawJERSFvsTime() {

  drawJERSFvsTimes("Prompt24");
  drawJERSFvsTimes("Prompt25");
  drawJERSFvsTimes("Prompt26");
  drawJERSFvsTimes("Prompt24to26");
}
  
void drawJERSFvsTimes(string set) {

  const char *cs = set.c_str();
  
  gROOT->ProcessLine(".! mkdir pdf/drawJERSFvsTime");
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! touch pdf/drawJERSFvsTime");
  gROOT->ProcessLine(".! touch pdf");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  string s = "rootfiles/JERSF.root";
  TFile *f = new TFile(s.c_str(),"READ");
  if (!f || f->IsZombie()) cout << "Missing file "<<s<<endl<<flush;
  assert(f && !f->IsZombie());
  curdir->cd();

  double eps(1e-4);
  double ymin(0.95+eps), ymax(1.9-eps), ymind(0.75-eps), ymaxd(1.35+eps);
  if (set=="Prompt24") {
    ymin = 0.95+eps; ymax = 1.9-eps; ymind = 0.75-eps; ymaxd = 1.25+eps;
  }
  else if (set=="Prompt25") {
    ymin = 0.95+eps; ymax = 1.9-eps; ymind = 0.80-eps; ymaxd = 1.20+eps;
  }
  else if (set=="Prompt26") {
    ymin = 0.95+eps; ymax = 1.9-eps; ymind = 0.80-eps; ymaxd = 1.20+eps;
  }
  else if (set=="Prompt24to26") {
    ymin = 0.95+eps; ymax = 1.9-eps; ymind = 0.75-eps; ymaxd = 1.35+eps;
  }
  
  double xmin(0), xmax(5.2);
  TH1D *hu = tdrHist(Form("hu_%s",cs),
		     "JER scale factor",ymin,ymax,"|#eta|",xmin,xmax);
  TH1D *hd = tdrHist(Form("hd_%s",cs),
		     "Ratio to 26C",ymind,ymaxd,"|#eta|",xmin,xmax);
  
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
  map<string, int> lwid;
  if (set=="Prompt24to26") {
    ref = "24CDE";//"26C";
    mlabel["2024_nib"] = "2024CDEFGHI";
    mlabel["2024CDE_nib"] = "24CDE (re-reco DI)";
    mlabel["2024F_nib1"] = "24F_nib1 (prompt DD)";
    mlabel["2024FGHI_nib"] = "24FGHI (prompt DI)";
    mlabel["2025C"] = "25C (prompt pre-timing)";
    mlabel["2025CDEFG"] = "25CDEFG (prompt DD)";
    lwid["2025CDEFG"] = 3;
    mlabel["2026B"] = "26B (high PU)";
    mlabel["2026C"] = "26C (low PU)";
    lwid["2026C"] = 2;
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
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.15,xmax,1.15);
  l->SetLineStyle(kDashDotted);
  l->DrawLine(xmin,1.50,xmax,1.50);

  TLegend *leg = tdrLeg(0.17,0.76-0.03*(nrun-0),0.42,0.76); // 2024
  leg->SetTextSize(0.035); // 2024
  
  gPad->RedrawAxis();
  
  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->DrawLine(xmin,1,xmax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(xmin,1.10,xmax,1.10);
  l->DrawLine(xmin,0.90,xmax,0.90);

  
  TH1D *h(0), *href(0);
  double pt1(0), pt2(0);
  for (int irun = 0; irun != nrun; ++irun) {

    string run = vrun[irun];
    const char *cr = run.c_str();

    curdir->cd();

    const char *cm = "Summer24MG_NOJERSF";
    string s = Form("Fits/h2jersf_%s_%s",cr,cm);
    TH2D *h2 = (TH2D*)f->Get(s.c_str());
    if (!h2) cout << "Missing histogram " << s << endl << flush;
    assert(h2);
    double pt = 100;
    int i = h2->GetYaxis()->FindBin(pt);
    pt1 = h2->GetYaxis()->GetBinLowEdge(i);
    pt2 = h2->GetYaxis()->GetBinLowEdge(i+1);
    h = h2->ProjectionX(Form("h_%s_%s",cr,cs),i,i);

    if (!href) href = (TH1D*)h->Clone(Form("href_%s",cs));
    TH1D *hr = (TH1D*)h->Clone(Form("hr_%s_%s",cr,cs));
    hr->Divide(href);
    
    c1->cd(1);
    tdrDraw(h,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);
    if (lwid[run]!=0) h->SetLineWidth(lwid[run]);
    /*
    int fillcolor(mcolor[run]-9);
    if (mcolor[run]==kBlack) fillcolor = kGray+1;
    if (fcolor[run]!=0) fillcolor = fcolor[run];
    tdrDraw(h,"E3",kNone,mcolor[run],kSolid,-1,1001,fillcolor);
    h->SetFillColorAlpha(fillcolor,0.7);
    tdrDraw(new TGraph(h),"L",kNone,mcolor[run],kSolid,-1,kNone,0);
    */
    
    if (mlabel[cr]) leg->AddEntry(h,mlabel[cr],"LE");
    else            leg->AddEntry(h,cr,"LE");
    
    c1->cd(2);
    tdrDraw(hr,"HISTE",kNone,mcolor[run],kSolid,-1,kNone,0);
    if (lwid[run]!=0) hr->SetLineWidth(lwid[run]);
    //tdrDraw(hr,"E3",kNone,mcolor[run],kSolid,-1,1001,fillcolor);
    //hr->SetFillColorAlpha(fillcolor,0.7);
    //tdrDraw(new TGraph(hr),"L",kNone,mcolor[run],kSolid,-1,kNone,0);

  }

  c1->cd(1);
  gPad->RedrawAxis();

  if (doRawData)
    tex->DrawLatex(0.36,0.85,"Merged Z+jet, #gamma+jet and dijet data");
  else
    //tex->DrawLatex(0.36,0.85,"Global fit of Z/#gamma+jet and multijet");
    tex->DrawLatex(0.36,0.85,"Fit of dijet data vs Summer24");
  tex->DrawLatex(0.36,0.80,"p_{T,jet} = 100 GeV, no PU reweighing");

  c1->cd(2);
  gPad->RedrawAxis();
  
  c1->SaveAs(Form("pdf/drawJERSFvsTime/drawJERSFvsTime_%s.pdf",cs));
} // drawL2ResVsTime
