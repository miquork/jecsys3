// Purpose: Compare two sets of MC truths stored in rootfiles/MCTruth.root
//          Currently setup for EEZS vs nominal Winter25
#include "TFile.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"

void compareMCTruth() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/compareMCTruth");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/compareMCTruth");
  
  TFile *f = new TFile("rootfiles/MCTruth.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  ////////////////////////
  // JES comparisons    //
  ////////////////////////
  
  TH2D *h2jes_var = (TH2D*)f->Get("h2jes_Winter25MC_Flat_EEZS9p5");
  assert(h2jes_var);
  //TH2D *h2jes_ref = (TH2D*)f->Get("h2jes_Winter25MG");
  TH2D *h2jes_ref = (TH2D*)f->Get("h2jes_Winter25MC_Flat22");
  assert(h2jes_ref);

  TH2D *h2jesFit_var = (TH2D*)f->Get("h2jesFit_Winter25MC_Flat_EEZS9p5");
  assert(h2jesFit_var);
  //TH2D *h2jesFit_ref = (TH2D*)f->Get("h2jesFit_Winter25MG");
  TH2D *h2jesFit_ref = (TH2D*)f->Get("h2jesFit_Winter25MC_Flat22");
  assert(h2jesFit_ref);

  double eps = 1e-4;
  double yup1(0.50+eps), yup2(1.25-eps);
  TH1D *hup = tdrHist("hup","JES",yup1,yup2,"|#eta|",0,5.2);
  double ydw1(0.95+eps), ydw2(1.07-eps);
  TH1D *hdw = tdrHist("hdw","Ratio",ydw1,ydw2,"|#eta|",0,5.2);
  
  extraText = "Private";
  //lumi_136TeV = "EEZS 9.5 vs Winter25, Simulation";
  lumi_136TeV = "Simulation";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,8,11);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  c1->cd(1);
  TLegend *leg = tdrLeg(0.32,0.89-5*0.05,0.57,0.89);
  leg->SetHeader("EEZS 9.5 vs Winter25");

  //double vpt[] = {30,100,300};
  //int color[] = {kBlue, kGreen+2, kRed};
  double vpt[] = {30,300,100};
  int color[] = {kBlue, kRed, kGreen+2};
  const int npt = sizeof(vpt)/sizeof(vpt[0]);

  for (int jpt = 0; jpt != npt; ++jpt) {
    
    double pt = vpt[jpt];//100.;
    int ipt = h2jes_ref->GetYaxis()->FindBin(pt);
    TH1D *h1jes_var = h2jes_var->ProjectionX(Form("h1jes_var_%d",jpt),ipt,ipt);
    TH1D *h1jes_ref = h2jes_ref->ProjectionX(Form("h1jes_ref_%d",jpt),ipt,ipt);
    TH1D *hr = (TH1D*)h1jes_var->Clone(Form("hr_%d",jpt));
    hr->Divide(h1jes_var,h1jes_ref,1,1,"");//"B");

    TH1D *h1jesFit_var = h2jesFit_var->ProjectionX(Form("h1jesFit_var_%d",jpt),ipt,ipt);
    TH1D *h1jesFit_ref = h2jesFit_ref->ProjectionX(Form("h1jesFit_ref_%d",jpt),ipt,ipt);
    TH1D *hrFit = (TH1D*)h1jesFit_var->Clone(Form("hrFit_%d",jpt));
    hrFit->Divide(h1jesFit_var,h1jesFit_ref,1,1,"");//"B");
    
    c1->cd(1);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,yup1,1.305,1.00);
    l->DrawLine(2.964,yup1,2.964,1.00);
    
    tdrDraw(h1jesFit_ref,"HIST][",kNone,color[jpt]-9,kSolid,-1,kNone,0);
    tdrDraw(h1jesFit_var,"HIST][",kNone,color[jpt],kSolid,-1,kNone,0);
    tdrDraw(h1jes_ref,"Pz",kOpenCircle,color[jpt]-9,kSolid,-1,kNone,0,0.7);
    tdrDraw(h1jes_var,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);

    leg->AddEntry(h1jes_var,Form("p_{T} = %1.0f GeV",pt),"PLE");
    if (jpt == npt-1) {
      leg->AddEntry(h1jes_ref,"Winter25");
      leg->AddEntry(h1jesFit_var,"JES fits");
    }
    
    c1->cd(2);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,ydw1,1.305,ydw2);
    l->DrawLine(2.964,ydw1,2.964,ydw2);
    tdrDraw(hrFit,"HIST][",kNone,color[jpt],kSolid,-1,kNone,0);
    //if (pt==100) {
    tdrDraw(hr,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);
    //}
    //else {
    //TGraph *gr = new TGraph(hr);
    //tdrDraw(gr,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);
    //}
      
  }
  
  c1->cd(1);
  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->RedrawAxis();

  c1->SaveAs("pdf/compareMCTruth/compareMCTruth_JES_EEZS9p6.pdf");



  ////////////////////////
  // JER comparisons    //
  ////////////////////////
  
  TH2D *h2jer_var = (TH2D*)f->Get("h2jer_Winter25MC_Flat_EEZS9p5");
  assert(h2jer_var);
  //TH2D *h2jer_ref = (TH2D*)f->Get("h2jer_Winter25MG");
  TH2D *h2jer_ref = (TH2D*)f->Get("h2jer_Winter25MC_Flat22");
  assert(h2jer_ref);

  TH2D *h2jerFit_var = (TH2D*)f->Get("h2jerFit_Winter25MC_Flat_EEZS9p5");
  assert(h2jerFit_var);
  //TH2D *h2jerFit_ref = (TH2D*)f->Get("h2jerFit_Winter25MG");
  TH2D *h2jerFit_ref = (TH2D*)f->Get("h2jerFit_Winter25MC_Flat22");
  assert(h2jerFit_ref);

  //double eps = 1e-4;
  double y2up1(0.0+eps), y2up2(1.25-eps);
  TH1D *hup2 = tdrHist("hup2","JER",y2up1,y2up2,"|#eta|",0,5.2);
  double y2dw1(0.80+eps), y2dw2(1.20-eps);
  TH1D *hdw2 = tdrHist("hdw2","Ratio",y2dw1,y2dw2,"|#eta|",0,5.2);
  
  TCanvas *c2 = tdrDiCanvas("c2",hup2,hdw2,8,11);

  c2->cd(1);
  TLegend *leg2 = tdrLeg(0.32,0.89-5*0.05,0.57,0.89);
  leg2->SetHeader("EEZS 9.5 vs Winter25");

  //double vpt[] = {30,300,100};
  //int color[] = {kBlue, kRed, kGreen+2};
  //const int npt = sizeof(vpt)/sizeof(vpt[0]);

  for (int jpt = 0; jpt != npt; ++jpt) {
    
    double pt = vpt[jpt];//100.;
    int ipt = h2jer_ref->GetYaxis()->FindBin(pt);
    TH1D *h1jer_var = h2jer_var->ProjectionX(Form("h1jer_var_%d",jpt),ipt,ipt);
    TH1D *h1jer_ref = h2jer_ref->ProjectionX(Form("h1jer_ref_%d",jpt),ipt,ipt);
    TH1D *hr = (TH1D*)h1jer_var->Clone(Form("hr2_%d",jpt));
    hr->Divide(h1jer_var,h1jer_ref,1,1,"");//"B");

    TH1D *h1jerFit_var = h2jerFit_var->ProjectionX(Form("h1jerFit_var_%d",jpt),ipt,ipt);
    TH1D *h1jerFit_ref = h2jerFit_ref->ProjectionX(Form("h1jerFit_ref_%d",jpt),ipt,ipt);
    TH1D *hrFit = (TH1D*)h1jerFit_var->Clone(Form("hr2Fit_%d",jpt));
    hrFit->Divide(h1jerFit_var,h1jerFit_ref,1,1,"");//"B");
    
    c2->cd(1);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,y2up1,1.305,1.00);
    l->DrawLine(2.964,y2up1,2.964,1.00);
    
    tdrDraw(h1jerFit_ref,"HIST][",kNone,color[jpt]-9,kSolid,-1,kNone,0);
    tdrDraw(h1jerFit_var,"HIST][",kNone,color[jpt],kSolid,-1,kNone,0);
    tdrDraw(h1jer_ref,"Pz",kOpenCircle,color[jpt]-9,kSolid,-1,kNone,0,0.7);
    tdrDraw(h1jer_var,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);

    leg2->AddEntry(h1jer_var,Form("p_{T} = %1.0f GeV",pt),"PLE");
    if (jpt == npt-1) {
      leg2->AddEntry(h1jer_ref,"Winter25");
      leg2->AddEntry(h1jerFit_var,"JER fits");
    }
    
    c2->cd(2);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,y2dw1,1.305,y2dw2);
    l->DrawLine(2.964,y2dw1,2.964,y2dw2);
    tdrDraw(hrFit,"HIST][",kNone,color[jpt],kSolid,-1,kNone,0);
    //if (pt==100) {
    //tdrDraw(hr,"HPz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);
    //}
    //else {
    TGraph *gr = new TGraph(hr);
    tdrDraw(gr,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);
    //}
  } // jpt
  
  c2->cd(1);
  gPad->RedrawAxis();
  
  c2->cd(2);
  gPad->RedrawAxis();

  c2->SaveAs("pdf/compareMCTruth/compareMCTruth_JER_EEZS9p6.pdf");



  ////////////////////////
  // EFF comparisons    //
  ////////////////////////
  
  TH2D *h2eff_var = (TH2D*)f->Get("h2effRaw_Winter25MC_Flat_EEZS9p5");
  assert(h2eff_var);
  TH2D *h2eff_ref = (TH2D*)f->Get("h2effRaw_Winter25MC_Flat22");
  assert(h2eff_ref);

  TH2D *h2effFit_var = (TH2D*)f->Get("h2effFit_Winter25MC_Flat_EEZS9p5");
  assert(h2effFit_var);
  TH2D *h2effFit_ref = (TH2D*)f->Get("h2effFit_Winter25MC_Flat22");
  assert(h2effFit_ref);

  double y3up1(0.0+eps), y3up2(1.50-eps);
  TH1D *hup3 = tdrHist("hup3","Efficiency",y3up1,y3up2,"|#eta|",0,5.2);
  //double ydw1(0.95+eps), ydw2(1.07-eps); // JES
  double y3dw1(0.80+eps), y3dw2(1.20-eps); // EFF
  TH1D *hdw3 = tdrHist("hdw3","Ratio",y3dw1,y3dw2,"|#eta|",0,5.2);
  
  TCanvas *c3 = tdrDiCanvas("c3",hup3,hdw3,8,11);

  c3->cd(1);
  TLegend *leg3 = tdrLeg(0.32,0.89-5*0.05,0.57,0.89);
  leg3->SetHeader("EEZS 9.5 vs Winter25");

  double vpt3[] = {15,30,100};
  int color3[] = {kBlue, kRed, kGreen+2};
  const int npt3 = sizeof(vpt3)/sizeof(vpt3[0]);

  for (int jpt = 0; jpt != npt3; ++jpt) {
    
    double pt = vpt3[jpt];
    int ipt = h2eff_ref->GetYaxis()->FindBin(pt);
    TH1D *h1eff_var = h2eff_var->ProjectionX(Form("h1eff_var_%d",jpt),ipt,ipt);
    TH1D *h1eff_ref = h2eff_ref->ProjectionX(Form("h1eff_ref_%d",jpt),ipt,ipt);
    TH1D *hr = (TH1D*)h1eff_var->Clone(Form("hr3_%d",jpt));
    hr->Divide(h1eff_var,h1eff_ref,1,1,"");//"B");

    TH1D *h1effFit_var = h2effFit_var->ProjectionX(Form("h1effFit_var_%d",jpt),ipt,ipt);
    TH1D *h1effFit_ref = h2effFit_ref->ProjectionX(Form("h1effFit_ref_%d",jpt),ipt,ipt);
    TH1D *hrFit = (TH1D*)h1effFit_var->Clone(Form("hr3Fit_%d",jpt));
    hrFit->Divide(h1effFit_var,h1effFit_ref,1,1,"");//"B");
    
    c3->cd(1);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,y2up1,1.305,1.00);
    l->DrawLine(2.964,y2up1,2.964,1.00);
    
    tdrDraw(h1effFit_ref,"HIST][",kNone,color[jpt]-9,kSolid,-1,kNone,0);
    tdrDraw(h1effFit_var,"HIST][",kNone,color[jpt],kSolid,-1,kNone,0);
    tdrDraw(h1eff_ref,"Pz",kOpenCircle,color[jpt]-9,kSolid,-1,kNone,0,0.7);
    tdrDraw(h1eff_var,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);

    leg3->AddEntry(h1eff_var,Form("p_{T} = %1.0f GeV",pt),"PLE");
    if (jpt == npt-1) {
      leg3->AddEntry(h1eff_ref,"Winter25");
      leg3->AddEntry(h1effFit_var,"EFF fits");
    }
    
    c3->cd(2);
    
    l->DrawLine(0,1,5.2,1);
    l->DrawLine(1.305,y2dw1,1.305,y2dw2);
    l->DrawLine(2.964,y2dw1,2.964,y2dw2);
    tdrDraw(hrFit,"HIST][",kNone,color3[jpt],kSolid,-1,kNone,0);
    //TGraph *gr = new TGraph(hr);
    //tdrDraw(gr,"Pz",kFullCircle,color[jpt],kSolid,-1,kNone,0,0.7);
    tdrDraw(hr,"Pz",kFullCircle,color3[jpt],kSolid,-1,kNone,0,0.7);
  } // jpt
  
  c3->cd(1);
  gPad->RedrawAxis();
  
  c3->cd(2);
  gPad->RedrawAxis();

  c3->SaveAs("pdf/compareMCTruth/compareMCTruth_EFF_EEZS9p6.pdf");
} // void compareMCTruth
