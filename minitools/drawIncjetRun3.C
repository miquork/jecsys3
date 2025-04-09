// Purpose: Draw inclusive jet spectrum in Run3 to estimate time stability
#include "TFile.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

TGraphErrors* remapHist(TH1D *h) {
  TGraphErrors *g = new TGraphErrors(0);
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    double x = h->GetBinCenter(i);
    double ex = 0.5*h->GetBinWidth(i);
    if (y!=0 || ey!=0) {
      int n = g->GetN();
      g->SetPoint(n, 0.5*x, y);
      g->SetPointError(n, 0.5*ex, ey);
    }
  }
  return g;
} // remapHist

void drawIncjetRun3() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  string files[] =
    {"jmenano_data_cmb_2223_JME_v1.root",
     "jmenano_data_cmb_2022CD_JME_v1.root",
     //jmenano_data_cmb_2022C_JME_v1.root
     //jmenano_data_cmb_2022D_JME_v1.root
     //"jmenano_data_cmb_2022EFG_JME_v1.root",
     "jmenano_data_cmb_2022E_JME_v1.root",
     "jmenano_data_cmb_2022FG_JME_v1.root",
     //jmenano_data_cmb_2022F_JME_v1.root
     //jmenano_data_cmb_2022G_JME_v1.root
     //jmenano_data_cmb_2223_JME_v1.root,
     "nano_data_cmb_2023BCv123_JME_v1.root",
     //nano_data_cmb_2023B_JME_v1.root
     //nano_data_cmb_2023Cv123_JME_v1.root
     "nano_data_cmb_2023Cv4D_JME_v1.root",
     //nano_data_cmb_2023Cv4_JME_v1.root
     //nano_data_cmb_2023D_JME_v1.root
     //"jmenano_mc_cmb_Summer22Both_v1.root",
     "jmenano_mc_out_Summer22Both_v1.root",
     //"jmenano_mc_cmb_Summer22_v1.root",
     //"jmenano_mc_cmb_Summer22EE_v1.root"
    };
  int nfiles = sizeof(files)/sizeof(files[0]);
  //double k1(1.65), k2(0.47), k3(1.4), k4(1.48);
  //double k1(1), k1b(1), k2(1), k2b(1), k3(1), k4(1);
  //double k1(1.70), k1b(1.00), k2b(0.33), k3(1.42), k4(1.5), kmc(1);
  double k1(1.70*0.95), k1b(1.00*1.05*0.85), k2b(0.33*1.05), k3(1.42*0.8), k4(1.5*1.05), kmc(1);
  double lumi[] =
    {//(5071+3006)*k1+(5878+18007+3122)*k2+(589+8728)*k3+(4271+9525)*k4,
      (5071+3006)*k1+(5878)*k1b+(18007+3122)*k2b+(589+8728)*k3+(4271+9525)*k4,
     (5071+3006)*k1,
     (5878)*k1b,
     (18007+3122)*k2b,
     (589+8728)*k3,
     (4271+9525)*k4,
     1e-3/1e5*7.4e4,
    };
  int nlumi = sizeof(lumi)/sizeof(lumi[0]);
  assert(nfiles==nlumi);
  int color[] =
    /*(
    {kBlack,
     kBlue,
     kGreen+2,
     kOrange+2,
     kRed};
    */
    //{kGray, kBlack, kYellow+1, kRed, kGreen+2, kBlue, kOrange+2};
    {kGray, kBlue, kYellow+1, kBlack, kRed, kGreen+2, kOrange+2};
  string legs[] =
    {"Run3", "22CD", "22E", "22FG", "23BC123", "23C4D", "22(EE) MC"};

  double ptmin = 700; // 650
  double ptmax = 3800;//5000;
  TH1D *hu = tdrHist("hu","Xsec (pb/GeV)",//1e-7,1e3,
		     1e-8,1e2,"p_{T} (GeV)",ptmin,ptmax);
  TH1D *hd = tdrHist("hu","Ratio",0.5,1.5,//0,3,
		     "p_{T} (GeV)",ptmin,ptmax);
  TH1D *h2 = tdrHist("h2","Era / Run3",0+1e-4,2.0-1e-4,
		     "p_{T} or M_{jj}/2 (GeV)",ptmin,ptmax);
  lumi_136TeV = "Run3";
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,8,11);
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);
  
  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  
  TLegend *leg = tdrLeg(0.60,0.85,0.90,0.85);
  
  c1->cd(2);
  gPad->SetLogx();
  
  c2->cd();
  gPad->SetLogx();
  TLegend *leg2 = tdrLeg(0.17,0.15,0.47,0.15);
  
  TH1D *href(0);
  const char *cd = "rootfiles/Iita_20230824_jetveto";
  for (int i = 0; i != nfiles; ++i) {
    
    TFile *f = new TFile(Form("%s/%s",cd,files[i].c_str()),"READ");
    assert(f && !f->IsZombie());
    const char *run;
    
    TH1D *h = (TH1D*)f->Get("Incjet/hpt13");
    if (!h) h = (TH1D*)f->Get("HLT_MC/Incjet/hpt13");
    assert(h);
    h = (TH1D*)h->Clone(Form("h%d",i));
    h->Scale(1./lumi[i],"width");
    h->GetXaxis()->SetRangeUser(ptmin,ptmax);
    
    c1->cd(1);
    tdrDraw(h,"Pz",kFullCircle,color[i]);
    
    leg->AddEntry(h,legs[i].c_str(),"PLE");
    leg->SetY1(leg->GetY1()-0.05);

    c1->cd(2);
    TH1D *hr = (TH1D*)h->Clone(Form("hr%d",i));
    if (!href) href = (TH1D*)h->Clone("href");
    hr->Divide(href);
    if (TString(files[i].c_str()).Contains("mc"))
      tdrDraw(hr,"HIST",kNone,color[i],kSolid,-1,kNone);
    else
      tdrDraw(hr,"Pz",kFullCircle,color[i]);
    
    c2->cd();
    
    if (legs[i]=="Run3") {}
    else if (TString(files[i].c_str()).Contains("mc"))
      tdrDraw(hr,"HIST",kNone,color[i],kSolid,-1,kNone);
    else
      tdrDraw(hr,"Pz",kFullCircle,color[i]);
    
    if (legs[i]!="Run3") {
      leg2->AddEntry(h,legs[i].c_str(),"PLE");
      leg2->SetY2(leg2->GetY2()+0.04);
    }
    
  }
  
  c1->cd(1);
  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/drawIncJetRun3.pdf");
  
  //TFile *f2 = new TFile("rootfiles/DijetMass_Ratios_inRun3.root");
  TFile *f2 = new TFile("rootfiles/DijetMass_Ratios_inRun3-1.root");
  assert(f2 && !f2->IsZombie());
  
  string hists[] =
    //{"2022CDE_vs_fullRun3", "", "2022FG_vs_fullRun3",
    {"2022CD_vs_fullRun3", "2022E_vs_fullRun3", "2022FG_vs_fullRun3",
     "2023BC123_vs_fullRun3", "2023C4D_vs_fullRun3"};
  const int nh = sizeof(hists)/sizeof(hists[0]);
  
  for (int i = 0; i != nh; ++i) {
   
    if (hists[i]=="") continue;
    TH1D *h = (TH1D*)f2->Get(hists[i].c_str()); assert(h);
    
    TGraphErrors *g = remapHist(h);

    c2->cd();
    tdrDraw(g,"Pz",kOpenSquare,color[i+1]);
  }

  c2->SaveAs("pdf/drawIncJetRun3_withMjj.pdf");
}
