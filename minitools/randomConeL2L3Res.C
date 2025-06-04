// Purpose: compare Random Cone data/MC scale factor to L2L3Res data,
//          specifically gamma+jet, Z+jet, dijet(ZeroBias) and their fit
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TPRegexp.h"

#include "../tdrstyle_mod22.C"

void randomConeL2L3Res_era(string era = "2024I");

void randomConeL2L3Res() {
  /*
  randomConeL2L3Res_era("2024B");
  randomConeL2L3Res_era("2024C");
  randomConeL2L3Res_era("2024D");
  randomConeL2L3Res_era("2024E-v1");
  randomConeL2L3Res_era("2024E-v2");
  randomConeL2L3Res_era("2024F");
  randomConeL2L3Res_era("2024G");
  randomConeL2L3Res_era("2024H");
  randomConeL2L3Res_era("2024I");
  */

  randomConeL2L3Res_era("2024C_nib1");
  randomConeL2L3Res_era("2024D_nib1");
  randomConeL2L3Res_era("2024Ev1_nib1");
  randomConeL2L3Res_era("2024Ev2_nib1");
  randomConeL2L3Res_era("2024F_nib1");
  randomConeL2L3Res_era("2024F_nib2");
  randomConeL2L3Res_era("2024F_nib3");
  randomConeL2L3Res_era("2024G_nib1");
  randomConeL2L3Res_era("2024G_nib2");
  randomConeL2L3Res_era("2024H_nib1");
  randomConeL2L3Res_era("2024I_nib1");

  randomConeL2L3Res_era("2025B");
  randomConeL2L3Res_era("2025C");
  
  //exit(0);
}

void randomConeL2L3Res_era(string era) {

  cout << "Processing era " << era << endl << flush;
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  // Get RC SF
  //TFile *frc = new TFile("rootfiles/Hirak_RC_compare_sf_2024_Winter24MC.root","READ");
  //TFile *frc = new TFile("rootfiles/Hirak_RC_compare_sf_2024_Summer24MC-v2.root","READ");
  // NB: this file has doubled 24C,D,E-v1,E-v2,F,G; others ok
  TFile *frc = new TFile("rootfiles/Hirak_compare_sf_Run2024_Summer24MC-merge3.root","READ");
  assert(frc && !frc->IsZombie());

  TFile *frc25 = new TFile("rootfiles/Hirak_2025B/compare_sf_Run2024GHI-Summer24MC_Run2025B-Winter25MC.root","READ");
  assert(frc25 && !frc25->IsZombie());

  
  curdir->cd();

  const char *cr = era.c_str();
  TString s = cr; TPRegexp r("_nib[0-9]+");
  if (s.Contains("2024H") || s.Contains("2024I"))
      r.Substitute(s, "", "g");
  s.ReplaceAll("Ev","E-v");
  s.ReplaceAll("_nib","-nib");
  // Patch for missing re-reco RC
  if (s.Contains("2024C") || s.Contains("2024D") || s.Contains("2024E"))
    //s = "2024F";
    s.ReplaceAll("nib1","reReco");
  if (s.Contains("2024E")) {
    s.ReplaceAll("E-v1","E"); s.ReplaceAll("E-v2","E");
  }
    
  const char *cr2 = s.Data();
  cout << "Reading RC for era " << cr2 << endl << flush;
  //TH1D *hrc = (TH1D*)frc->Get("Run2024I"); assert(hrc);
  //TH1D *hrc = (TH1D*)frc->Get(Form("Run%s",cr)); assert(hrc);
  //TH1D *hrc = (TH1D*)frc->Get(Form("Run%s",cr2)); assert(hrc);
  TH1D *hrc(0);
  if (s.Contains("2024")) hrc = (TH1D*)frc->Get(Form("Run%s",cr2));
  if (s.Contains("2025")) hrc = (TH1D*)frc25->Get(Form("Run%s",cr2));
  if (s.Contains("2025C") && !hrc) hrc = (TH1D*)frc25->Get("Run2025B");
  assert(hrc);
  
  // Correct for sigmaMB=69.2mb->75.3mb
  TH1D *hrc753 = (TH1D*)hrc->Clone(Form("hrc753_%s",cr));
  hrc753->Scale(69.2/75.3);

  // For 2024CDE re-reco, apply extra scaling to match HF
  if (s.Contains("reReco") &&
      (s.Contains("2024C") || s.Contains("2024D") || s.Contains("2024E"))) {
    cout << "Apply extra scaling to reReco" << endl;
    hrc753->Scale(pow(69.2/75.3,2));
  }

  // Scale tracker-covered region by 0.4 to emulate difference from
  // primarily having differences from calorimeter side only
  TH1D *hrc753jes = (TH1D*)hrc753->Clone(Form("hrc753jes_%s",cr));
  for (int i = 1; i != hrc753->GetNbinsX()+1; ++i) {
    double eta = hrc753->GetBinCenter(i);
    double r = hrc753->GetBinContent(i);
    // For 2024CDE re-reco, apply scaling to EB only, not EE
    if (s.Contains("reReco") &&
	(s.Contains("2024C") || s.Contains("2024D") || s.Contains("2024E"))) {
      if (fabs(eta)<1.218) // HE-1
	hrc753jes->SetBinContent(i, (r-1)*0.4+1);
      else if (fabs(eta)<1.305) // HB-HE
	hrc753jes->SetBinContent(i, (r-1)*0.5+1);
      else if (fabs(eta)<1.395) // HB+1,EE-1
	hrc753jes->SetBinContent(i, (r-1)*0.6+1);
      else if (fabs(eta)<1.479) // EB-EE
	hrc753jes->SetBinContent(i, (r-1)*0.7+1);
    }
    else if (s.Contains("2024F-nib1")) {
      // do nothing
    }
    else if (s.Contains("2025")) {
      if (fabs(eta)<1.044)
	hrc753jes->SetBinContent(i, (r-1)*0.8+1);
      else if (fabs(eta)<1.305)
	hrc753jes->SetBinContent(i, (r-1)*0.7+1);
      else if (fabs(eta)<2.65)
	hrc753jes->SetBinContent(i, (r-1)*0.6+1);
      // do nothing
    }
    else if (fabs(eta)<2.65)
      hrc753jes->SetBinContent(i, (r-1)*0.4+1);
  } // for i

  // Turn to L2Res by scaling out |eta|<1.3
  TH1D *hrc13 = (TH1D*)hrc753->Clone(Form("hrc13_%s",cr));
  TF1 *f13 = new TF1(Form("f13_%s",cr),"[0]",-1.305,1.305);
  hrc753->Fit(f13,"QRNW");
  hrc13->Scale(1./f13->GetParameter(0));
  TH1D *hrc13jes = (TH1D*)hrc753jes->Clone(Form("hrc13jes_%s",cr));
  TF1 *f13jes = new TF1(Form("f13jes_%s",cr),"[0]",-1.305,1.305);
  hrc753jes->Fit(f13jes,"QRNW");
  hrc13jes->Scale(1./f13jes->GetParameter(0));
  // Remove RC error for later ratio, and symmetrize vs eta
  for (int i = 1; i != hrc13->GetNbinsX()+1; ++i) {
    hrc13->SetBinError(i, 0.);
    hrc13jes->SetBinError(i, 0.);
    double eta = hrc13jes->GetBinCenter(i);
    if (eta<0) {
      int j = hrc13jes->GetXaxis()->FindBin(-eta);
      hrc13jes->SetBinContent(i, 0.5*(hrc13jes->GetBinContent(i) +
				      hrc13jes->GetBinContent(j)));
    }
    else {
      int j = hrc13jes->GetXaxis()->FindBin(-eta);
      hrc13jes->SetBinContent(i, hrc13jes->GetBinContent(j));
    }
  } // for i

  // Add some reasonable uncertainties for L2Res fits later
  for (int i = 1; i != hrc13jes->GetNbinsX()+1; ++i) {
    double eta = fabs(hrc13jes->GetBinCenter(i));
    double dr = fabs(hrc13jes->GetBinContent(i)-1);
    if      (eta>3.139) hrc13jes->SetBinError(i, max(0.25*dr,0.02));
    else if (eta>2.65)  hrc13jes->SetBinError(i, 0.100);
    else if (eta>2.5)   hrc13jes->SetBinError(i, 0.070);
    else if (eta>2.322) hrc13jes->SetBinError(i, 0.050);
    else if (eta>1.305) hrc13jes->SetBinError(i, 0.030);
    else if (eta>=0.00) hrc13jes->SetBinError(i, 0.010);
  }
  

  
  
  // Get L2Res
  //TFile *fl2 = new TFile("rootfiles/L2Res.root","READ");
  //TFile *fl2 = new TFile("rootfiles/L2Res_Prompt2024_V8M.root","READ");
  TFile *fl2(0);
  if (s.Contains("24")) fl2=new TFile("rootfiles/L2Res_V9M_draft1.root","READ");
  if (s.Contains("25")) fl2=new TFile("rootfiles/L2Res.root","READ");
  assert(fl2 && !fl2->IsZombie());

  //map<string, const char*> m;
  //m["2024B"] = "2024B_nib1";
  //m["2024C"] = "2024C_nib1";
  //m["2024D"] = "2024D_nib1";
  //m["2024E-v1"] = "2024Ev1_nib1";
  //m["2024E-v2"] = "2024Ev2_nib1";
  //m["2024F"] = "2024F_nib2";
  //m["2024G"] = "2024G_nib1";
  //m["2024H"] = "2024H_nib1";
  //m["2024I"] = "2024I_nib1";
  //TFile *fl3 = new TFile(Form("rootfiles/jecdata_V8M/jecdata%s.root",m[era]),
  //TFile *fl3 = new TFile(Form("rootfiles/jecdata%s.root",m[era]),
  TFile *fl3 = new TFile(Form("rootfiles/jecdata%s.root",cr),
			 "READ");
  assert(fl3 && !fl3->IsZombie());

  curdir->cd();
  
  //string era2 = TString(cr).ReplaceAll("-","").Data();
  //const char *cr2 = era2.c_str();
  
  //TH2D *h2l2 = (TH2D*)fl2->Get("h2jes_2024I_nib1"); assert(h2l2);
  //TH2D *h2l2 = (TH2D*)fl2->Get(Form("h2jes_%s_nib1",cr2)); assert(h2l2);
  //TH2D *h2l2 = (TH2D*)fl2->Get(Form("h2jes_%s",m[era])); assert(h2l2);
  TH2D *h2l2 = (TH2D*)fl2->Get(Form("h2jes_%s",cr)); assert(h2l2);
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
  //TH1D *h = tdrHist(Form("h_%s",cr),"JES or RC SF ",0.4+eps,1.8-eps,
  TH1D *h = tdrHist(Form("h_%s",cr),"JES or RC",0.25+eps,1.7-eps,
		    "#eta",-5.2,5.2);
  //TH1D *hd = tdrHist(Form("hd_%s",cr),"JES / RC SF ",0.8+eps,1.8-eps,
  TH1D *hd = tdrHist(Form("hd_%s",cr),"JES / RCF",0.4+eps,1.6-eps,
		     "#eta",-5.2,5.2);
  //lumi_136TeV = Form("%s, %s",cr,mlum[Form("%s_nib1",cr2)].c_str());
  //lumi_136TeV = Form("%s, %s",cr,mlum[m[era]].c_str());
  lumi_136TeV = Form("%s, %s",cr,mlum[cr].c_str());
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
  tdrDraw(hrc13jes,"E2",kNone,kBlack,kSolid,-1,1001,kYellow+1);
  tdrDraw(hrc13,"HISTE",kNone,kGray+1,kSolid,-1,kNone);
  TH1D *hrc13jes_clone = (TH1D*)hrc13jes->Clone(Form("hrc13jes_clone_%s",cr));
  tdrDraw(hrc13jes_clone,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hl2_1_rc,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_1_rc,"Pz",kFullCircle,kGray+1,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_2_rc,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_1,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hl2_2,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0,0.7);

  //TLegend *leg = tdrLeg(0.50,0.90-0.05*5,0.75,0.90);
  TLegend *leg = tdrLeg(0.45,0.90-0.05*3,0.70,0.90);
  //leg->AddEntry(hrc,"RC@69.2 mb","L");
  //leg->AddEntry(hrc753,"RC@75.3 mb","L");
  leg->AddEntry(hrc13,"RC SF / (|#eta|<1.3)","L");
  leg->AddEntry(hrc13jes,"RC#otimesf(neutral) / (|#eta|<1.3)","FL");
  leg->AddEntry(hl2_1,"JES(L2Res)@10 GeV","PLE");
  //leg->AddEntry(hl2_2,"JES@15 GeV","L");
  //leg->AddEntry(hl2_2,"JES@20 GeV","L");

  gPad->RedrawAxis();
  
  c1->cd(2);

  TH1D *hr_rc = (TH1D*)hrc->Clone(Form("hr_rc_%s",cr));
  hr_rc->Divide(hrc13jes);
  TH1D *hr_rc753 = (TH1D*)hrc753->Clone(Form("hr753_rc753_%s",cr));
  hr_rc753->Divide(hrc13jes);
  TH1D *hr_1 = (TH1D*)hl2_1_rc->Clone(Form("hr_1_%s",cr));
  hr_1->Divide(hrc13jes);
  TH1D *hr_2 = (TH1D*)hl2_2_rc->Clone(Form("hr_2_%s",cr));
  hr_2->Divide(hrc13jes);

  TH1D *hr_1_band = (TH1D*)hr_1->Clone(Form("hr_1_band_%s",cr));
  for (int i = 1; i != hr_1_band->GetNbinsX()+1; ++i) {
    hr_1_band->SetBinContent(i, 1);
    hr_1_band->SetBinError(i, hrc13jes->GetBinError(i));
  }
  tdrDraw(hr_1_band,"E2",kFullCircle,kBlack,kSolid,-1,1001,kYellow+1,0.5);

  l->SetLineStyle(kSolid);
  l->SetLineColor(kBlack);
  l->DrawLine(-5.2,1,5.2,1);
  
  //tdrDraw(hr_rc,"HIST",kNone,kGray,kSolid,-1,kNone);
  //tdrDraw(hr_rc753,"HIST",kNone,kGray+1,kSolid,-1,kNone);
  //tdrDraw(hr_1,"Pz",kFullCircle,kGray+1,kSolid,-1,kNone,0,0.5);
  tdrDraw(hr_1,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.5);
  //tdrDraw(hr_2,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone,0,0.5);
  
  gPad->RedrawAxis();

  gROOT->ProcessLine(".! touch pdf/randomConeL2L3Res");
  //c1->SaveAs("pdf/randomConeL2L3Res/randomConeL2L3Res.pdf");
  c1->SaveAs(Form("pdf/randomConeL2L3Res/randomConeL2L3Res_%s.pdf",cr));

  TFile *fout = new TFile("rootfiles/randomConeL2L3Res.root","UPDATE");
  hrc13jes->Write(Form("hrc13jes_%s",cr), TObject::kOverwrite);
  fout->Close();
  
  frc->Close();
  frc25->Close();
  fl2->Close();
  fl3->Close();
  
} // randomConeL2L3Res
