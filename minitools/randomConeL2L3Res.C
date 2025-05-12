// Purpose: compare Random Cone data/MC scale factor to L2L3Res data,
//          specifically gamma+jet, Z+jet, dijet(ZeroBias) and their fit
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

void randomConeL2L3Res_era(string era = "2024I");

void randomConeL2L3Res() {
  randomConeL2L3Res_era("2024B");
  randomConeL2L3Res_era("2024C");
  randomConeL2L3Res_era("2024D");
  randomConeL2L3Res_era("2024E-v1");
  randomConeL2L3Res_era("2024E-v2");
  randomConeL2L3Res_era("2024F");
  randomConeL2L3Res_era("2024G");
  randomConeL2L3Res_era("2024H");
  randomConeL2L3Res_era("2024I");

  exit(0);
}

void randomConeL2L3Res_era(string era) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // Get RC SF
  TFile *frc = new TFile("rootfiles/Hirak_RC_compare_sf_2024_Winter24MC.root","READ");
  assert(frc && !frc->IsZombie());

  curdir->cd();

  const char *cr = era.c_str();
  //TH1D *hrc = (TH1D*)frc->Get("Run2024I"); assert(hrc);
  TH1D *hrc = (TH1D*)frc->Get(Form("Run%s",cr)); assert(hrc);

  // Correct for sigmaMB=69.2mb->75.3mb
  TH1D *hrc753 = (TH1D*)hrc->Clone(Form("hrc753_%s",cr));
  hrc753->Scale(69.2/75.3);

  // Turn to L2Res by scaling out |eta|<1.3
  TH1D *hrc13 = (TH1D*)hrc753->Clone(Form("hrc13_%s",cr));
  TF1 *f13 = new TF1(Form("f13_%s",cr),"[0]",-1.305,1.305);
  hrc753->Fit(f13,"QRNW");
  hrc13->Scale(1./f13->GetParameter(0));
  // Remove RC error for later ratio
  for (int i = 1; i != hrc13->GetNbinsX()+1; ++i) {
    hrc13->SetBinError(i, 0.);
  } // for i
  
  // Get L2Res
  //TFile *fl2 = new TFile("rootfiles/L2Res.root","READ");
  TFile *fl2 = new TFile("rootfiles/L2Res_Prompt2024_V8M.root","READ");
  assert(fl2 && !fl2->IsZombie());

  map<string, const char*> m;
  m["2024B"] = "2024B_nib1";
  m["2024C"] = "2024C_nib1";
  m["2024D"] = "2024D_nib1";
  m["2024E-v1"] = "2024Ev1_nib1";
  m["2024E-v2"] = "2024Ev2_nib1";
  m["2024F"] = "2024F_nib2";
  m["2024G"] = "2024G_nib1";
  m["2024H"] = "2024H_nib1";
  m["2024I"] = "2024I_nib1";
  TFile *fl3 = new TFile(Form("rootfiles/jecdata_V8M/jecdata%s.root",m[era]),
			 "READ");
  assert(fl3 && !fl3->IsZombie());

  curdir->cd();
  
  string era2 = TString(cr).ReplaceAll("-","").Data();
  const char *cr2 = era2.c_str();
  
  //TH2D *h2l2 = (TH2D*)fl2->Get("h2jes_2024I_nib1"); assert(h2l2);
  //TH2D *h2l2 = (TH2D*)fl2->Get(Form("h2jes_%s_nib1",cr2)); assert(h2l2);
  TH2D *h2l2 = (TH2D*)fl2->Get(Form("h2jes_%s",m[era])); assert(h2l2);
  double pt1 = 10;
  int ipt1 = h2l2->GetYaxis()->FindBin(pt1);
  TH1D *hl2_1 = h2l2->ProjectionX(Form("hl2_%s_1",cr),ipt1,ipt1);
  double pt2 = 20;//15;
  int ipt2 = h2l2->GetYaxis()->FindBin(pt2);
  TH1D *hl2_2 = h2l2->ProjectionX(Form("hl2_%s_2",cr),ipt2,ipt2);

  TH1D *hl3 = (TH1D*)fl3->Get("ratio/eta00-13/run3/hFit_Rjet");
  assert(hl3);
  
  // Copy L2Res*L3Res to RC bins, especially to +eta and -eta
  TH1D *hl2_1_rc = (TH1D*)hrc13->Clone(Form("hl2_%s_1_rc",cr));
  hl2_1_rc->Clear();
  TH1D *hl2_2_rc = (TH1D*)hrc13->Clone(Form("hl2_%s_2_rc",cr));
  hl2_2_rc->Clear();
  for (int i = 1; i != hl2_1_rc->GetNbinsX()+1; ++i) {

    int k10 = hl3->FindBin(10.);
    double r10 = 1;//;hl3->GetBinContent(k10);
    double er10 = 0;//hl3->GetBinError(k10);

    int k15 = hl3->FindBin(15.);
    double r15 = 1;//hl3->GetBinContent(k15);
    double er15 = 0;//hl3->GetBinError(k15);

    int k20 = hl3->FindBin(20.);
    double r20 = 1;//hl3->GetBinContent(k20);
    double er20 = 0;//hl3->GetBinError(k20);
    
    double eta = hl2_1_rc->GetBinCenter(i);
    int j = hl2_1->FindBin(fabs(eta));
    hl2_1_rc->SetBinContent(i, r10*hl2_1->GetBinContent(j));
    hl2_1_rc->SetBinError(i, r10*hl2_1->GetBinError(j));
    hl2_2_rc->SetBinContent(i, r20*hl2_2->GetBinContent(j));
    hl2_2_rc->SetBinError(i, r20*hl2_2->GetBinError(j));
  } // for i
  
  curdir->cd();

  // Get L3Res
  //TFile *fl3 = new TFile(Form("rootfiles/jecdata%s_nib1.root",cr),"READ");
  //assert(fl3 && !fl3->IsZombie());
  
  // Draw comparison
  #include "../Config.C"
  double eps = 1e-4;
  TH1D *h = tdrHist(Form("h_%s",cr),"JES or RC SF ",0.4+eps,1.8-eps,
		    "#eta",-5.2,5.2);
  TH1D *hd = tdrHist(Form("hd_%s",cr),"JES / RC SF ",0.8+eps,1.8-eps,
		     "#eta",-5.2,5.2);
  //lumi_136TeV = Form("%s, %s",cr,mlum[Form("%s_nib1",cr2)].c_str());
  lumi_136TeV = Form("%s, %s",cr,mlum[m[era]].c_str());
  extraText = "Private";
  //TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cr),h,hd,8,11);

  TLine *l = new TLine();

  c1->cd(1);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(-5.2,1,5.2,1);
  
  //tdrDraw(hrc,"HIST",kNone,kGray,kSolid,-1,kNone);
  //tdrDraw(hrc753,"HIST",kNone,kGray+1,kSolid,-1,kNone);
  tdrDraw(hrc13,"HISTE",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hl2_1_rc,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_1_rc,"Pz",kFullCircle,kGray+1,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_2_rc,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_1,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_2,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0,0.7);

  //TLegend *leg = tdrLeg(0.50,0.90-0.05*5,0.75,0.90);
  TLegend *leg = tdrLeg(0.50,0.90-0.05*3,0.75,0.90);
  //leg->AddEntry(hrc,"RC@69.2 mb","L");
  //leg->AddEntry(hrc753,"RC@75.3 mb","L");
  leg->AddEntry(hrc13,"RC/(|#eta|<1.3)","L");
  leg->AddEntry(hl2_1,"JES(L2Res)@10 GeV","L");
  //leg->AddEntry(hl2_2,"JES@15 GeV","L");
  //leg->AddEntry(hl2_2,"JES@20 GeV","L");

  gPad->RedrawAxis();
  
  c1->cd(2);

  l->SetLineStyle(kSolid);
  l->SetLineColor(kBlack);
  l->DrawLine(-5.2,1,5.2,1);

  TH1D *hr_rc = (TH1D*)hrc->Clone(Form("hr_rc_%s",cr));
  hr_rc->Divide(hrc13);
  TH1D *hr_rc753 = (TH1D*)hrc753->Clone(Form("hr753_rc753_%s",cr));
  hr_rc753->Divide(hrc13);
  TH1D *hr_1 = (TH1D*)hl2_1_rc->Clone(Form("hr_1_%s",cr));
  hr_1->Divide(hrc13);
  TH1D *hr_2 = (TH1D*)hl2_2_rc->Clone(Form("hr_2_%s",cr));
  hr_2->Divide(hrc13);

  //tdrDraw(hr_rc,"HIST",kNone,kGray,kSolid,-1,kNone);
  //tdrDraw(hr_rc753,"HIST",kNone,kGray+1,kSolid,-1,kNone);
  //tdrDraw(hr_1,"Pz",kFullCircle,kGray+1,kSolid,-1,kNone,0,0.5);
  tdrDraw(hr_1,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.5);
  //tdrDraw(hr_2,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone,0,0.5);

  gPad->RedrawAxis();
  
  //c1->SaveAs("pdf/randomConeL2L3Res/randomConeL2L3Res.pdf");
  c1->SaveAs(Form("pdf/randomConeL2L3Res/randomConeL2L3Res_%s.pdf",cr));

  frc->Close();
  fl2->Close();
  fl3->Close();
  
} // randomConeL2L3Res
