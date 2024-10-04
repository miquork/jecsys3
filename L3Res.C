// Purpose: Derive L3Res from Z+jet, gamma+jet, multijet based on new TProfile2D
//          This is to migrate away from GlobalFit.C to a more unified input
// To-do:   Add PF composition

#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#include <iostream>
#include <fstream>

#include "tdrstyle_mod22.C"

TLegend *_leg(0);
bool fitZ = true; // Z+jet
bool fitG = true; // gamma+jet
bool fitD = true; // Multijet (pT,ave)
bool fitP = true; // Multijet (pT,probe)
bool fitJ = true; // Multijet (pT,tag)

bool doClosure = false; // Do not undo L2L3Res for closure test

// Global variable for renaming histograms
string _run;


// Step 1. Slice 1D profile out of 2D in given range and draw it
TProfile* drawPt(TProfile2D *p2, double etamin, double etamax,
		 string draw, int marker, int color, string label="",
		 TProfile2D *p2x = 0) {
  
  assert(p2);
  int ix1 = p2->GetXaxis()->FindBin(etamin);
  int ix2 = p2->GetXaxis()->FindBin(etamax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(),
		   1000*etamin, 1000*etamax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileY(Form("p%s",id.c_str()),ix1,ix2);

  if (p2x) {
    TProfile *py = p2x->ProfileY(Form("px%s",id.c_str()),ix1,ix2);
    for (int i = 0; i != p->GetNbinsX()+1; ++i) {
      int j = py->GetXaxis()->FindBin(p->GetBinCenter(i));
      double mean_value = p->GetBinContent(i) * py->GetBinContent(j);
      double entries = p->GetBinEntries(i);
      p->SetBinEntries(i, entries);
      p->SetBinContent(i, entries * mean_value);
    }
  }
  
  tdrDraw(p,draw,marker,color,kSolid,-1,kNone);

  if (_leg && label!="") {
    double eta1 = p2->GetXaxis()->GetBinLowEdge(ix1);
    double eta2 = p2->GetXaxis()->GetBinLowEdge(ix2+1);
    _leg->AddEntry(p, Form("%s+jet",label.c_str()), "PLE");
    _leg->SetY1(_leg->GetY1()-0.05);
  }

  return p;
} // drawPt

// Step 3,4. Take a ratio and draw it
TH1D *drawRatio(TH1D *h, TH1D *hm, string draw, int marker, int color) {
  assert(h);
  assert(hm);
  string id = Form("%s_%s_%s_%d_%d", h->GetName(), hm->GetName(),
		   draw.c_str(), marker, color);
  TH1D *hr = (TH1D*)h->Clone(Form("hr_%s_%s",id.c_str(),_run.c_str()));
  hr->Divide(hm);

  tdrDraw(hr,draw,marker,color,kSolid,-1,kNone);
  
  return hr;
}

// Step 5. Clean up data (for fitting) and draw it
TH1D *drawCleaned(TH1D *h, string data, string draw,
		  int marker, int color) {
  assert(h);
  TH1D *hc = (TH1D*)h->Clone(Form("hc_%s_%s",h->GetName(),_run.c_str()));
  for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
    double pt = hc->GetBinCenter(i);
    double ptmin = hc->GetBinLowEdge(i);
    double ptmax = hc->GetBinLowEdge(i+1);
    double emax = ptmax * cosh(1.3);
    double sqrts = 13600.;
    double keep(false);

    // L2Res pT binning (central+forward hybrid)
    //{15, 21, 28, 37, 49, 59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575, 638, 737, 846, 967, 1101, 1248, 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
    // L2Res eta binning
    //{0., 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.839, 4.013, 4.583, 5.191};
    // New L2Res eta binning (same as L2Relative)
    //{0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

    // Prompt24
    if (data=="Z") {
      if      (ptmin>=15. && ptmax<700) keep = true;
    }
    if (data=="G") {
      if      (ptmin>=30. && ptmax<1300.) keep = true;
    }
    if (data=="J") {
      if       (ptmin>=15. && ptmax<2500.) keep = true;
    }
    if (data=="P") {
      if       (ptmin>=15. && ptmax<2500.) keep = true;
    }
    if (data=="D") {
      if       (ptmin>=15. && ptmax<2500.) keep = true;

    }
    
    // Check that no points out of bound
    if (h->GetBinContent(i)>1.3 || h->GetBinContent(i)<0.3) {
      keep = false;
    }

    // Remove points we don't want to keep
    if (!keep) {
      hc->SetBinContent(i, 0.);
      hc->SetBinError(i, 0.);
    }
    // Set minimum uncertainty for others
    else {
      double errmin = 0.002;
      hc->SetBinError(i, sqrt(pow(h->GetBinError(i),2)+pow(errmin,2)));
    }
  } // for i

  tdrDraw(hc,draw,marker,color,kSolid,-1,kNone);
  
  return hc;
}
// Helper for TMultiGraphErrors
TGraphErrors *cleanGraph(TGraphErrors* g, bool clone = false) {
  assert(g);
  assert(!clone);
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0 && g->GetEY()[i]==0)
      g->RemovePoint(i);
  }
  return g;
}
// Step 7. Draw h2jes
TH1D *drawH2JES(TH2D *h2, double pt, string draw, int marker, int color) {
  int ipt = h2->GetYaxis()->FindBin(pt);
  string id = Form("h%s_%s_%d_%d",h2->GetName(),draw.c_str(),marker,color);
  TH1D *h = h2->ProjectionX(Form("h%s",id.c_str()),ipt,ipt);

  tdrDraw(h,draw.c_str(),marker,color,kSolid,-1,kNone);

  return h;
} // drawH2JES


void L3Res() {

  // Set graphical styles
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Set output directory;
  TFile *fout = new TFile("rootfiles/L3Res.root","RECREATE");

  // Make sure graphics output directories exists
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/L3Res");

  map<string, string> mlum;
  mlum["2024BCD"] = "12.3 fb^{-1}";
  mlum["2024E"] = "X fb^{-1}";
  mlum["2024F"] = "X fb^{-1}";
  
  string vrun[] = {"2024BCD","2024E","2024F"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  //string vmc[] = {"Summer23BPix","Summer23BPix"};
  string vmc[] = {"Winter24","Winter24"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);
  assert(nmc==nrun);

  string sc = (doClosure ? "_closure" : "");
  const char *cc = sc.c_str();
  
  for (int irun = 0; irun != nrun; ++irun) {
    string run = vrun[irun];
    const char *cr = run.c_str();
    _run = run;
    string mc = vmc[irun];
    const char *cm = mc.c_str();
    string lum = (mlum[run]!="" ? Form(", %s",mlum[run].c_str()) : "");
    const char *cl = lum.c_str();
    
    //gROOT->ProcessLine(Form(".! mkdir pdf/L3Res/%s",cr));

  // (No indent here for the resf of the loop, maybe function call later)


  // Mikko's temporary code to create L2Res.root file
  // Load Z+jet
  TFile *fz(0);
  if (run=="2024CR") {
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ");
  }
  else if (TString(cr).Contains("2024")) {
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root",cr),"READ");
  }
  assert(fz && !fz->IsZombie());

  TDirectory *dz = fz->GetDirectory("data/l2res"); assert(dz);
  TDirectory *dzm = fz->GetDirectory("mc/l2res"); assert(dzm);
  
  // Load G+jet
  TFile *fg(0), *fgm(0);
  if (run=="2024CR") {
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w29.root","2024C"),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
  }
  else if (TString(cr).Contains("2024")) {
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w29.root",cr),"READ"); // June 6 hybrid
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
  }
  assert(fg && !fg->IsZombie());
  assert(fgm && !fgm->IsZombie());
  
  TDirectory *dg = fg->GetDirectory("Gamjet2"); assert(dg);
  TDirectory *dgm = fgm->GetDirectory("Gamjet2"); assert(dgm);

  // Load dijet
  TFile *fd(0), *fdm(0);
  if (run=="2024CR") {
    fd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_ECAL_JME_v76_2024.root","2024C"),"READ");
    fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (TString(cr).Contains("2024")) {
    fd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root",cr),"READ");
    fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  assert(fd && !fd->IsZombie());
  assert(fdm && !fdm->IsZombie());

  TDirectory *dd = fd->GetDirectory("Dijet2"); assert(dd);
  TDirectory *ddm = fdm->GetDirectory("HLT_MC/Dijet2"); assert(ddm);
  
  // Get the TProfile2D inputs. Z+jet, gam+jet, 3x dijet
  
  // Z+jet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2z = (TProfile2D*)dz->Get("p2m0tc"); assert(p2z);
  TProfile2D *p2zm = (TProfile2D*)dzm->Get("p2m0tc"); assert(p2zm);
  TProfile2D *p2zx = (TProfile2D*)dz->Get("p2restc"); assert(p2zx);
  
  // G+jet: x:eta, y:pT,avp (no other variants yet)
  TProfile2D *p2g = (TProfile2D*)dg->Get("p2m0"); assert(p2g);
  TProfile2D *p2gm = (TProfile2D*)dgm->Get("p2m0"); assert(p2gm);
  TProfile2D *p2gx = (TProfile2D*)dg->Get("p2res"); assert(p2gx);

  // Dijet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2j = (TProfile2D*)dd->Get("p2m0tc"); assert(p2j);
  TProfile2D *p2jm = (TProfile2D*)ddm->Get("p2m0tc"); assert(p2jm);
  TProfile2D *p2jx = (TProfile2D*)dd->Get("p2restc"); assert(p2jx);
  TProfile2D *p2p = (TProfile2D*)dd->Get("p2m0pf"); assert(p2p);
  TProfile2D *p2pm = (TProfile2D*)ddm->Get("p2m0pf"); assert(p2pm);
  TProfile2D *p2px = (TProfile2D*)dd->Get("p2respf"); assert(p2px);
  TProfile2D *p2d = (TProfile2D*)dd->Get("p2m0"); assert(p2d);
  TProfile2D *p2dm = (TProfile2D*)ddm->Get("p2m0"); assert(p2dm);
  TProfile2D *p2dx = (TProfile2D*)dd->Get("p2res"); assert(p2dx);

  // Reset L2L3Res to zero for closure tests so it's not undone later
  if (doClosure) {
    p2zx = p2gx = p2jx = p2px = p2dx = 0;
  }
  
  /*
  // Use common file instead
  TFile *f = new TFile("rootfiles/L2Res_2023_Summer22.root","READ");
  //TFile *f = new TFile("rootfiles/L2Res_2023_Winter23.root","READ");
  assert(f && !f->IsZombie());

  TProfile2D *p2z, *p2zm, *p2g, *p2gm;
  p2z = (TProfile2D*)f->Get(Form("p2m0tc_zjet_data_%s",cr)); assert(p2z);
  p2zm = (TProfile2D*)f->Get(Form("p2m0tc_zjet_mc_%s",cm)); assert(p2zm);
  p2g = (TProfile2D*)f->Get(Form("p2m0tc_gjet_data_%s",cr)); assert(p2g);
  p2gm = (TProfile2D*)f->Get(Form("p2m0tc_gjet_mc_%s",cm)); assert(p2gm);

  TProfile2D *p2j, *p2jm, *p2p, *p2pm, *p2d, *p2dm;
  p2j = (TProfile2D*)f->Get(Form("p2m0tc_dijet_data_%s",cr)); assert(p2j);
  p2jm = (TProfile2D*)f->Get(Form("p2m0tc_dijet_mc_%s",cm)); assert(p2jm);
  p2p = (TProfile2D*)f->Get(Form("p2m0pf_dijet_data_%s",cr)); assert(p2p);
  p2pm = (TProfile2D*)f->Get(Form("p2m0pf_dijet_mc_%s",cm)); assert(p2pm);
  p2d = (TProfile2D*)f->Get(Form("p2m0ab_dijet_data_%s",cr)); assert(p2d);
  p2dm = (TProfile2D*)f->Get(Form("p2m0ab_dijet_mc_%s",cm)); assert(p2dm);  
  */
  // Step 1. Extract barrel reference
  TProfile *p13zm(0), *p13gm(0), *p13dm(0), *p13jm(0), *p13pm(0);
  TProfile *p13z(0), *p13g(0), *p13d(0), *p13j(0), *p13p(0);

  TH1D *h13 = tdrHist("h13","JES",0.70,1.75);
  lumi_136TeV = Form("%s - %s%s",cr,cm,cl);
  extraText = "Private";
  TCanvas *c13 = tdrCanvas("c13",h13,8,33,kSquare);
  c13->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,1,3500,1);
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  TLegend *leg13 = tdrLeg(0.20,0.90,0.45,0.90);
  _leg = leg13;

  double etamin(0.), etamax(1.30);
  int ieta1 = p2d->GetXaxis()->FindBin(etamin);
  int ieta2 = p2d->GetXaxis()->FindBin(etamax);
  double eta1 = p2d->GetXaxis()->GetBinLowEdge(ieta1);
  double eta2 = p2d->GetXaxis()->GetBinLowEdge(ieta2+1);
  p13zm = drawPt(p2zm,etamin,etamax,"HISTE",kNone,kRed);
  p13gm = drawPt(p2gm,etamin,etamax,"HISTE",kNone,kBlue);
  p13jm = drawPt(p2jm,etamin,etamax,"HISTE",kNone,kGreen+2);
  p13pm = drawPt(p2pm,etamin,etamax,"HISTE",kNone,kOrange+2);
  p13dm = drawPt(p2dm,etamin,etamax,"HISTE",kNone,kBlack);
  
  p13z = drawPt(p2z,etamin,etamax,"Pz",kFullSquare,kRed,"Z",p2zx);
  p13g = drawPt(p2g,etamin,etamax,"Pz",kFullCircle,kBlue,"#gamma",p2gx);
  p13j = drawPt(p2j,etamin,etamax,"Pz",kFullDiamond,kGreen+2,"Tag",p2jx);
  p13p = drawPt(p2p,etamin,etamax,"Pz",kFullDiamond,kOrange+2,"Probe",p2px);
  p13d = drawPt(p2d,etamin,etamax,"Pz",kOpenDiamond,kBlack,"Dijet",p2dx);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c13->SaveAs(Form("pdf/L3Res/L3Res_%s%s_datamc.pdf",cr,cc));
  
  // Step 2. Draw data/MC ratio
  TH1D *h3 = tdrHist(Form("h3_%s",cr),"JES Data/MC",0.90,1.25);
  TCanvas *c3 = tdrCanvas(Form("c3_%s",cr),h3,8,33,kSquare);
  c3->SetLogx();

  l->DrawLine(15,1,3500,1);

  leg13->Draw("SAME");
  
  TH1D *hzr, *hgr, *hdr, *hjr, *hpr;
  hzr = drawRatio(p13z->ProjectionX(),p13zm,"Pz",kFullSquare,kRed);
  hgr = drawRatio(p13g->ProjectionX(),p13gm,"Pz",kFullCircle,kBlue);
  hjr = drawRatio(p13j->ProjectionX(),p13jm,"Pz",kFullDiamond,kGreen+2);
  hpr = drawRatio(p13p->ProjectionX(),p13pm,"Pz",kFullDiamond,kOrange+2);
  hdr = drawRatio(p13d->ProjectionX(),p13dm,"Pz",kOpenDiamond,kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c3->SaveAs(Form("pdf/L3Res/L3Res_%s%s_ratio.pdf",cr,cc));

  // Step 5. Curate data and do final fit of response
  // Options: flat, log-linear, quadratic, +1/x
  TH1D *h5 = tdrHist(Form("h5_%s",cr),"Data/MC",0.50,1.35);
  TCanvas *c5 = tdrCanvas(Form("c5_%s",cr),h5,8,33,kSquare);
  c5->SetLogx();

  double eta = 0.5*(eta1+eta2);
  double ptmaxe = 0.5*13600./cosh(eta1);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmaxe,0.50,ptmaxe,1.35);
  l->SetLineStyle(kDashed);
  l->DrawLine(15.,1,3500.,1);

  leg13->Draw("SAME");

  // Draw full range at the back
  tdrDraw(hzr,"Pz",kOpenSquare,kRed-9);
  tdrDraw(hgr,"Pz",kOpenCircle,kBlue-9);
  tdrDraw(hjr,"Pz",kOpenDiamond,kGreen+2-9);
  tdrDraw(hpr,"Pz",kOpenDiamond,kOrange+2-9);
  tdrDraw(hdr,"Pz",kOpenDiamond,kGray);

  TH1D *hzrf, *hgrf, *hdrf, *hjrf, *hprf;
  hzrf = drawCleaned(hzr,"Z","Pz",kFullSquare,kRed);
  hgrf = drawCleaned(hgr,"G","Pz",kFullCircle,kBlue);
  hjrf = drawCleaned(hjr,"J","Pz",kFullDiamond,kGreen+3);
  hprf = drawCleaned(hpr,"P","Pz",kFullDiamond,kOrange+2);
  hdrf = drawCleaned(hdr,"D","Pz",kFullDiamond,kBlack);

  TMultiGraph *mg = new TMultiGraph(Form("mg_%s",cr),"mg");
  if (fitZ) mg->Add(cleanGraph(new TGraphErrors(hzrf)),"SAMEP");
  if (fitG) mg->Add(cleanGraph(new TGraphErrors(hgrf)),"SAMEP");
  if (fitD) mg->Add(cleanGraph(new TGraphErrors(hdrf)),"SAMEP");
  if (fitP) mg->Add(cleanGraph(new TGraphErrors(hprf)),"SAMEP");
  if (fitJ) mg->Add(cleanGraph(new TGraphErrors(hjrf)),"SAMEP");
  mg->Draw("Pz");//"SAME Pz");

  TF1 *f0 = new TF1(Form("f0_%s",cr),"[0]",15.,3500.);
  TF1 *f1 = new TF1(Form("f1_%s",cr),
		    "[0]+[1]*log10(0.01*x)",15.,3500.);
  TF1 *f2 = new TF1(Form("f2_%s",cr),
		    "[0]+[1]*log10(0.01*x)+[2]/(x/10.)",15.,3500.);
  TF1 *f3 = new TF1(Form("f3_%s",cr),
		    "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)",15,3500);
  TF1 *f4 = new TF1(Form("f4_%s",cr),
		    "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)"
		    "+[3]/(x/10.)",15.,3500.);
  // Reference JES
  TF1 *fref = new TF1(Form("fref_%s",cr),
		      "[0]+[1]*log10(0.01*x)+[2]/(x/10.)",15.,3500.);
  
  // Sequential fitting with more and more parameters
  // Extra parameters are of size of expected prior uncertainty
  f0->SetParameter(0,1.000);
  mg->Fit(f0,"QRN");
  f1->SetParameters(f0->GetParameter(0),-0.01);
  mg->Fit(f1,"QRN");
  f2->SetParameters(f1->GetParameter(0),f1->GetParameter(1),+0.02);
  f2->SetParLimits(2,-0.5,0.5); // offset no more than 50% at 10 GeV
  mg->Fit(f2,"QRN");
  f3->SetParameters(f1->GetParameter(0),f1->GetParameter(1),+0.005);
  mg->Fit(f3,"QRN");
  f4->SetParameters(f3->GetParameter(0),f3->GetParameter(1),
		    f3->GetParameter(2),f2->GetParameter(2));
  f4->SetParLimits(2,-0.5,0.5); // offset no more than 50% at 10 GeV
  mg->Fit(f4,"QRN");

  fref->SetParameters(f1->GetParameter(0),f1->GetParameter(1),+0.02);
  fref->SetParLimits(2,-0.5,0.5); // offset no more than 50% at 10 GeV
  mg->Fit(fref,"QRN");

  f0->Draw("SAME"); f0->SetLineColor(kMagenta+2);
  f1->Draw("SAME"); f1->SetLineColor(kBlue);
  f2->Draw("SAME"); f2->SetLineColor(kGreen+1);
  f3->Draw("SAME"); f3->SetLineColor(kOrange+2);
  f4->Draw("SAME"); f4->SetLineColor(kRed);

  fref->Draw("SAME"); fref->SetLineColor(kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c5->SaveAs(Form("pdf/L3Res/L3Res_%s%s_ratioc.pdf",cr,cc));
  

  // Step 9. Print out text files
  //ofstream ftxt(Form("textfiles/Prompt24/Prompt24_Run%s_V3M_DATA_L2Residual_AK4PFPuppi.txt",cr));
  //ftxt << Form("{ 1 JetEta 1 JetPt 1./(%s) Correction L2Relative}",
  //	       vf1[0]->GetExpFormula().Data()) << endl;
  if (fout) {
    fout->cd();
    //h2jes->Write(Form("h2jes_%s",cr),TObject::kOverwrite);
    p2z->Write(Form("p2m0tc_zjet_data_%s",cr),TObject::kOverwrite);
    p2zm->Write(Form("p2m0tc_zjet_mc_%s",cm),TObject::kOverwrite);
    p2g->Write(Form("p2m0tc_gjet_data_%s",cr),TObject::kOverwrite);
    p2gm->Write(Form("p2m0tc_gjet_mc_%s",cm),TObject::kOverwrite);
    p2j->Write(Form("p2m0tc_dijet_data_%s",cr),TObject::kOverwrite);
    p2jm->Write(Form("p2m0tc_dijet_mc_%s",cm),TObject::kOverwrite);
    p2p->Write(Form("p2m0pf_dijet_data_%s",cr),TObject::kOverwrite);
    p2pm->Write(Form("p2m0pf_dijet_mc_%s",cm),TObject::kOverwrite);
    p2d->Write(Form("p2m0ab_dijet_data_%s",cr),TObject::kOverwrite);
    p2dm->Write(Form("p2m0ab_dijet_mc_%s",cm),TObject::kOverwrite);
    curdir->cd();
  }

  } // for irun

  if (fout) {
    fout->Write();
    fout->Close();
  }
} // L3Res
