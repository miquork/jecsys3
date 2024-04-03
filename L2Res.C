// Purpose: derive L2Res from Z+jet, gamma+jet, dijet based on new TProfile2D
// To-do: update to finer L2Relative bins to fully capture detector structure,
//        especially |eta|=1.5 with EB-EE transition, but also HF slope

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
bool fitD = true; // Dijet (pT,ave)
bool fitP = true; // Dijet (pT,probe)
bool fitJ = true; // Dijet (pT,tag)

// Step 1. Slice 1D profile out of 2D in given range and draw it
TProfile* drawEta(TProfile2D *p2, double ptmin, double ptmax,
		  string draw, int marker, int color, string label="") {
  
  assert(p2);
  int iy1 = p2->GetYaxis()->FindBin(ptmin);
  int iy2 = p2->GetYaxis()->FindBin(ptmax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(), ptmin, ptmax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileX(Form("p%s",id.c_str()),iy1,iy2);
  tdrDraw(p,draw,marker,color,kSolid,-1,kNone);

  if (_leg && label!="") {
    int pt1 = int(p2->GetYaxis()->GetBinLowEdge(iy1));
    int pt2 = int(p2->GetYaxis()->GetBinLowEdge(iy2+1));
    if (label == "Z")
      _leg->AddEntry(p, Form("%s+jet [%d,%d]",label.c_str(),pt1,pt2),"PLE");
    else
      _leg->AddEntry(p, Form("%s+jet",label.c_str()), "PLE");
    _leg->SetY1(_leg->GetY1()-0.05);
  }

  return p;
} // drawEta
TProfile* drawPt(TProfile2D *p2, double etamin, double etamax,
		 string draw, int marker, int color, string label="") {
  
  assert(p2);
  int ix1 = p2->GetXaxis()->FindBin(etamin);
  int ix2 = p2->GetXaxis()->FindBin(etamax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(),
		   1000*etamin, 1000*etamax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileY(Form("p%s",id.c_str()),ix1,ix2);
  tdrDraw(p,draw,marker,color,kSolid,-1,kNone);

  if (_leg && label!="") {
    double eta1 = p2->GetXaxis()->GetBinLowEdge(ix1);
    double eta2 = p2->GetXaxis()->GetBinLowEdge(ix2+1);
    //if (label == "Z")
    //_leg->AddEntry(p, Form("%s+jet [%1.3f,%1.3f]",label.c_str(),eta1,eta2),"PLE");
    //else
    _leg->AddEntry(p, Form("%s+jet",label.c_str()), "PLE");
    _leg->SetY1(_leg->GetY1()-0.05);
  }

  return p;
} // drawPt

// Step 2. Project profile to histogram, normalize by |eta|<1.3 and draw it
TH1D *drawNormEta(TProfile *p, string draw, int marker, int color) {

  assert(p);
  string id = Form("%s_%s_%d_%d",p->GetName(),draw.c_str(),marker,color);
  TH1D *h = p->ProjectionX(Form("h%s",id.c_str()));
  TF1 *f1 = new TF1(Form("f1%s_",id.c_str()),"[0]",0,1.305);
  h->Fit(f1,"QRN");
  f1->SetRange(0,5.2);
  h->Divide(f1);
  tdrDraw(h,draw,marker,color,kSolid,-1,kNone);
  
  return h;
}
TH1D *drawNormPt(TProfile *p, TProfile *p13,
		 string draw, int marker, int color) {

  assert(p);
  assert(p13);
  string id = Form("%s_%s_%d_%d",p->GetName(),draw.c_str(),marker,color);
  TH1D *h = p->ProjectionX(Form("h%s",id.c_str()));
  h->Divide(p13);
  tdrDraw(h,draw,marker,color,kSolid,-1,kNone);
  
  return h;
}

// Step 3,4. Take a ratio and draw it
TH1D *drawRatio(TH1D *h, TH1D *hm, string draw, int marker, int color) {
  assert(h);
  assert(hm);
  string id = Form("%s_%s_%s_%d_%d", h->GetName(), hm->GetName(),
		   draw.c_str(), marker, color);
  TH1D *hr = (TH1D*)h->Clone(Form("hr_%s",id.c_str()));
  hr->Divide(hm);

  tdrDraw(hr,draw,marker,color,kSolid,-1,kNone);
  
  return hr;
}

// Step 5. Clean up data (for fitting) and draw it
TH1D *drawCleaned(TH1D *h, double eta, string data, string draw,
		  int marker, int color) {
  assert(h);
  TH1D *hc = (TH1D*)h->Clone(Form("hc_%s",h->GetName()));
  for (int i = 1; i != hc->GetNbinsX()+1; ++i) {
    double pt = hc->GetBinCenter(i);
    double ptmin = hc->GetBinLowEdge(i);
    double ptmax = hc->GetBinLowEdge(i+1);
    double emax = ptmax * cosh(eta);
    double sqrts = 13600.;
    double keep(false);

    // L2Res pT binning (central+forward hybrid)
    //{15, 21, 28, 37, 49, 59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575, 638, 737, 846, 967, 1101, 1248, 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
    // L2Res eta binning
    //{0., 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.839, 4.013, 4.583, 5.191};

    if (data=="Z") {
      // Summer23
      if      (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.500) keep = true;
      else if (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.964) keep = true;
      else if (ptmin>=59. && ptmax<200 && eta>2.964  && eta<3.839) keep = true;
    }
    if (data=="G") {
      // Summer23
      if      (ptmin>=110. && ptmax<1300. && emax<2500.) keep = true;
      else if (ptmin>=110. && ptmax<200.  && emax<0.5*sqrts) keep = true;
    }
    if (data=="J") {
      // Summer23_v2
      if       (ptmin>=15. && ptmax<2500. && eta<1.044) keep = true;
      else if  (ptmin>=15. && ptmax<2000. && eta<2.500) keep = true;
      else if  (ptmin>=15. && ptmax<967.  && eta<2.964) keep = true;
      else if  (ptmin>=15. && ptmax<846.  && eta<3.839) keep = true;
      else if  (ptmin>=15. && ptmax<460.  && eta<5.191) keep = true;
    }
    if (data=="P") {
      // Summer23
      if (ptmin>=15. && ptmax<=110. && eta<2.322) keep = true;
    }
    if (data=="D") {
      // Summer23
      if (ptmin>=59. && ptmax<=110. && emax<3100. && eta<2.853) keep = true;
    }
    
    // Check that no points out of bound
    if (h->GetBinContent(i)>1.3 || h->GetBinContent(i)<0.3) {
      keep = false;
    }
    // Remove BPIX bad point, HF bad point
    if (data=="J" || data=="P" || data=="D") {
      if (eta>1.653 && eta < 1.930 && pt>86 && pt<110) keep = false;
      if (eta>4.018 && eta < 4.583 && pt>279 && pt<302) keep = false;
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


void L2Res() {

  // Set graphical styles
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Set output directory;
  TFile *fout = new TFile("rootfiles/L2Res.root","RECREATE");

  // Make sure graphics output directories exists
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/L2Res");
  gROOT->ProcessLine(".! mkdir pdf/L2Res/vsEta");
  gROOT->ProcessLine(".! mkdir pdf/L2Res/vsPt");

  string vrun[] = {"2023Cv123","2023Cv4","2023D"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  string vmc[] = {"Summer23","Summer23","Summer23BPIX"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);
  assert(nmc==nrun);
  for (int irun = 0; irun != nrun; ++irun) {
    string run = vrun[irun];
    const char *cr = run.c_str();
    string mc = vmc[irun];
    const char *cm = mc.c_str();
  // (No indent here for the resf of the loop, maybe function call later)


  // Mikko's temporary code to create L2Res.root file
  // Load Z+jet
  TFile *fz(0);
  fz = new TFile(Form("rootfiles/Summer23_L2ResOnly/jme_bplusZ_%s_Zmm_sync_v70.root",cr),"READ"); // Summer23 L2Res_V1
  assert(fz && !fz->IsZombie());

  TDirectory *dz = fz->GetDirectory("data/l2res");
  TDirectory *dzm = fz->GetDirectory("mc/l2res");
  
  // Load G+jet
  TFile *fg(0), *fgm(0);
  fg = new TFile(Form("rootfiles/Summer23_L2ResOnly/gamjet_all/GamHistosFill_data_%s_w4.root",cr),"READ"); // Summer23 with L2Res
  assert(fg && !fg->IsZombie());
  //
  fgm = new TFile(run=="2023D" ? "rootfiles/Summer23_L2ResOnly/gamjet_all/GamHistosFill_mc_2023P8-BPix_w4.root" : "rootfiles/Summer23_L2ResOnly/gamjet_all/GamHistosFill_mc_2023P8_w4.root","READ"); // Summer23 L2Res_V1
  assert(fgm && !fgm->IsZombie());
  
  TDirectory *dg = fg->GetDirectory("Gamjet2");
  TDirectory *dgm = fgm->GetDirectory("Gamjet2");

  // Load dijet
  TFile *fd(0), *fdm(0);
  fd = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_data_cmb_%s_JME_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cr),"READ"); // Summer23 L2Res_V1
  assert(fd && !fd->IsZombie());
  //
  fdm = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MG%s_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",run=="2023D" ? "BPix" : ""),"READ");
  assert(fdm && !fdm->IsZombie());

  TDirectory *dd = fd->GetDirectory("Dijet2");
  TDirectory *ddm = fdm->GetDirectory("Dijet2");
  
  // Get the TProfile2D inputs. Z+jet, gam+jet, 3x dijet
  
  // Z+jet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2z = (TProfile2D*)dz->Get("p2m0tc"); assert(p2z);
  TProfile2D *p2zm = (TProfile2D*)dzm->Get("p2m0tc"); assert(p2zm);

  // G+jet: x:eta, y:pT,avp (no other variants yet)
  TProfile2D *p2g = (TProfile2D*)dg->Get("p2m0"); assert(p2g);
  TProfile2D *p2gm = (TProfile2D*)dgm->Get("p2m0"); assert(p2gm);

  // Dijet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2j = (TProfile2D*)dd->Get("p2m0tc"); assert(p2j);
  TProfile2D *p2jm = (TProfile2D*)ddm->Get("p2m0tc"); assert(p2jm);
  TProfile2D *p2p = (TProfile2D*)dd->Get("p2m0pf"); assert(p2p);
  TProfile2D *p2pm = (TProfile2D*)ddm->Get("p2m0pf"); assert(p2pm);
  TProfile2D *p2d = (TProfile2D*)dd->Get("p2m0"); assert(p2d);
  TProfile2D *p2dm = (TProfile2D*)ddm->Get("p2m0"); assert(p2dm);

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
  // Step 0. Extract barrel reference to normalize it out later
  TProfile *p13zm(0), *p13gm(0), *p13dm(0), *p13jm(0), *p13pm(0);
  TProfile *p13z(0), *p13g(0), *p13d(0), *p13j(0), *p13p(0);

  TH1D *h13 = tdrHist("h13","JES",0.65,1.85);//0.65,1.45);
  lumi_136TeV = Form("%s - %s",cr,cm);
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
  
  p13z = drawPt(p2z,etamin,etamax,"Pz",kFullSquare,kRed,"Z");
  p13g = drawPt(p2g,etamin,etamax,"Pz",kFullCircle,kBlue,"#gamma");
  p13j = drawPt(p2j,etamin,etamax,"Pz",kFullDiamond,kGreen+2,"Tag");
  p13p = drawPt(p2p,etamin,etamax,"Pz",kFullDiamond,kOrange+2,"Probe");
  p13d = drawPt(p2d,etamin,etamax,"Pz",kOpenDiamond,kBlack,"Dijet");

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c13->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000.*eta2),cr,"c13"));
  
  // Create giant canvas for all eta bins (7*3=21; more in the future)
  int nxy = p2d->GetNbinsX();
  TCanvas *cx = new TCanvas("cx","cx",7*250,3*250);
  cx->Divide(7,3,0,0);
  TH2D *h2jes = p2d->ProjectionXY("h2jes"); h2jes->Reset();
  
  // Loop over the ieta bins
  vector<TF1*> vf1(p2d->GetNbinsX()+1);
  TH1D *hmin = p2d->ProjectionX(Form("hmin_%s",cr)); hmin->Reset();
  TH1D *hmax = p2d->ProjectionX(Form("hmax_%s",cr)); hmax->Reset();
  for (int ieta = 1; ieta != p2d->GetNbinsX()+1; ++ieta) {
    //int ieta = p2d->GetXaxis()->FindBin(2.7);
    double etamin = p2d->GetXaxis()->GetBinCenter(ieta);
    double etamax = etamin;
    double eta1 = p2d->GetXaxis()->GetBinLowEdge(ieta);
    double eta2 = p2d->GetXaxis()->GetBinLowEdge(ieta+1);
  // (No indent here for the resf of the loop, maybe function call later)
    
  TH1D *h1 = tdrHist("h1","JES",0.65,1.85);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,8,33,kSquare);
  c1->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,1,3500,1);

  TLegend *leg1 = tdrLeg(0.20,0.90,0.45,0.90);
  _leg = leg1;

  TProfile *pzm, *pgm, *pdm, *pjm, *ppm;
  pzm = drawPt(p2zm,etamin,etamax,"HISTE",kNone,kRed);
  pgm = drawPt(p2gm,etamin,etamax,"HISTE",kNone,kBlue);
  pjm = drawPt(p2jm,etamin,etamax,"HISTE",kNone,kGreen+2);
  ppm = drawPt(p2pm,etamin,etamax,"HISTE",kNone,kOrange+2);
  pdm = drawPt(p2dm,etamin,etamax,"HISTE",kNone,kBlack);

  TProfile *pz, *pg, *pd, *pj, *pp;
  pz = drawPt(p2z,etamin,etamax,"Pz",kFullSquare,kRed,"Z");
  pg = drawPt(p2g,etamin,etamax,"Pz",kFullCircle,kBlue,"#gamma");
  pj = drawPt(p2j,etamin,etamax,"Pz",kFullDiamond,kGreen+2,"Tag");
  pp = drawPt(p2p,etamin,etamax,"Pz",kFullDiamond,kOrange+2,"Probe");
  pd = drawPt(p2d,etamin,etamax,"Pz",kOpenDiamond,kBlack,"Dijet");

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c1->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000.*eta2),cr,"c1"));

  
  // Step 2. Project profile to histogram, normalize by |eta|<1.3
  TH1D *h2 = tdrHist("h2","Rel. JES",0.65,1.85);
  TCanvas *c2 = tdrCanvas("c2",h2,8,33,kSquare);
  c2->SetLogx();

  l->DrawLine(15,1,3500,1);

  leg1->Draw("SAME");
  
  TH1D *hzm, *hgm, *hdm, *hjm, *hpm;
  hzm = drawNormPt(pzm,p13zm,"HISTE",kNone,kRed);
  hgm = drawNormPt(pgm,p13gm,"HISTE",kNone,kBlue);
  hjm = drawNormPt(pjm,p13jm,"HISTE",kNone,kGreen+2);
  hpm = drawNormPt(ppm,p13pm,"HISTE",kNone,kOrange+2);
  hdm = drawNormPt(pdm,p13dm,"HISTE",kNone,kBlack);

  TH1D *hz, *hg, *hd, *hj, *hp;
  hz = drawNormPt(pz,p13z,"Pz",kFullSquare,kRed);
  hg = drawNormPt(pg,p13g,"Pz",kFullCircle,kBlue);
  hj = drawNormPt(pj,p13j,"Pz",kFullDiamond,kGreen+2);
  hp = drawNormPt(pp,p13p,"Pz",kFullDiamond,kOrange+2);
  hd = drawNormPt(pd,p13d,"Pz",kOpenDiamond,kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));

  c2->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000*eta2),cr,"c2"));

  // Step 3. Draw data/MC ratio before normalization
  TH1D *h3 = tdrHist("h3","JES Data/MC",0.5,1.3);//0.80,1.15);
  TCanvas *c3 = tdrCanvas("c3",h3,8,33,kSquare);
  c3->SetLogx();

  l->DrawLine(15,1,3500,1);

  leg1->Draw("SAME");
  
  TH1D *hzr, *hgr, *hdr, *hjr, *hpr;
  hzr = drawRatio(pz->ProjectionX(),pzm,"Pz",kFullSquare,kRed);
  hgr = drawRatio(pg->ProjectionX(),pgm,"Pz",kFullCircle,kBlue);
  hjr = drawRatio(pj->ProjectionX(),pjm,"Pz",kFullDiamond,kGreen+2);
  hpr = drawRatio(pp->ProjectionX(),ppm,"Pz",kFullDiamond,kOrange+2);
  hdr = drawRatio(pd->ProjectionX(),pdm,"Pz",kOpenDiamond,kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c3->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000.*eta2),cr,"c3"));

  // Step 4. Draw data/MC ratio of normalized JES
  TH1D *h4 = tdrHist("h4","Rel. JES Data/MC",0.50,1.35);
  TCanvas *c4 = tdrCanvas("c4",h4,8,33,kSquare);
  c4->SetLogx();

  l->DrawLine(15.,1,3500.,1);

  leg1->Draw("SAME");
	     
  TH1D *hzrn, *hgrn, *hdrn, *hjrn, *hprn;
  hzrn = drawRatio(hz,hzm,"Pz",kFullSquare,kRed);
  hgrn = drawRatio(hg,hgm,"Pz",kFullCircle,kBlue);
  hjrn = drawRatio(hj,hjm,"Pz",kFullDiamond,kGreen+2);
  hprn = drawRatio(hp,hpm,"Pz",kFullDiamond,kOrange+2);
  hdrn = drawRatio(hd,hdm,"Pz",kOpenDiamond,kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));

  c4->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000.*eta2),cr,"c4"));


  // Step 5. Curate data and do final fit of response
  // Options: flat, log-linear, quadratic, +1/x
  TH1D *h5 = tdrHist("h5","Rel. JES Data/MC",0.50,1.35);
  TCanvas *c5 = tdrCanvas("c5",h5,8,33,kSquare);
  c5->SetLogx();

  double eta = 0.5*(eta1+eta2);
  double ptmaxe = 0.5*13600./cosh(eta1);
  l->SetLineStyle(kDotted);
  l->DrawLine(ptmaxe,0.50,ptmaxe,1.35);
  l->SetLineStyle(kDashed);
  l->DrawLine(15.,1,3500.,1);

  leg1->Draw("SAME");

  // Draw full range at the back
  tdrDraw(hzrn,"Pz",kOpenSquare,kRed-9);
  tdrDraw(hgrn,"Pz",kOpenCircle,kBlue-9);
  tdrDraw(hjrn,"Pz",kOpenDiamond,kGreen+2-9);
  tdrDraw(hprn,"Pz",kOpenDiamond,kOrange+2-9);
  tdrDraw(hdrn,"Pz",kOpenDiamond,kGray);

  TH1D *hzrf, *hgrf, *hdrf, *hjrf, *hprf;
  hzrf = drawCleaned(hzrn,eta,"Z","Pz",kFullSquare,kRed);
  hgrf = drawCleaned(hgrn,eta,"G","Pz",kFullCircle,kBlue);
  hjrf = drawCleaned(hjrn,eta,"J","Pz",kFullDiamond,kGreen+3);
  hprf = drawCleaned(hprn,eta,"P","Pz",kFullDiamond,kOrange+2);
  hdrf = drawCleaned(hdrn,eta,"D","Pz",kFullDiamond,kBlack);

  TMultiGraph *mg = new TMultiGraph();
  if (fitZ) mg->Add(cleanGraph(new TGraphErrors(hzrf)),"SAMEP");
  if (fitG) mg->Add(cleanGraph(new TGraphErrors(hgrf)),"SAMEP");
  if (fitD) mg->Add(cleanGraph(new TGraphErrors(hdrf)),"SAMEP");
  if (fitP) mg->Add(cleanGraph(new TGraphErrors(hprf)),"SAMEP");
  if (fitJ) mg->Add(cleanGraph(new TGraphErrors(hjrf)),"SAMEP");
  mg->Draw();//"SAME Pz");

  TF1 *f0 = new TF1(Form("f0_%d_%s",ieta,cr),"[0]",15.,3500.);
  TF1 *f1 = new TF1(Form("f1_%d_%s",ieta,cr),
		    "[0]+[1]*log10(0.01*x)",15.,3500.);
  TF1 *f2 = new TF1(Form("f2_%d_%s",ieta,cr),
		    "[0]+[1]*log10(0.01*x)+[2]/(x/10.)",15.,3500.);
  TF1 *f3 = new TF1(Form("f3_%d_%s",ieta,cr),
		    "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)",15,3500);
  TF1 *f4 = new TF1(Form("f4_%d_%s",ieta,cr),
		    "[0]+[1]*log10(0.01*x)+[2]*pow(log10(0.01*x),2)"
		    "+[3]/(x/10.)",15.,3500.);

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

  // Bonus: constrain barrel f2 fit to flat if not significant
  if (eta<1.3 && false) {
    if (fabs(f2->GetParameter(1))<2.*f2->GetParError(1)) {
      f2->FixParameter(1,0.);
      mg->Fit(f2,"QRN");

      if (fabs(f2->GetParameter(2))<2.*f2->GetParError(2)) {
	f2->FixParameter(2,0.);
	mg->Fit(f2,"QRN");
      }
    }
  }
  
  // Keep track of log-lin+1/x for text files
  vf1[ieta-1] = f2;
  
  f0->Draw("SAME"); f0->SetLineColor(kMagenta+2); f0->SetLineStyle(kDashed);
  f1->Draw("SAME"); f1->SetLineColor(kBlue);
  f2->Draw("SAME"); f2->SetLineColor(kGreen+1);
  f3->Draw("SAME"); f3->SetLineColor(kOrange+2);
  f4->Draw("SAME"); f4->SetLineColor(kRed);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c5->SaveAs(Form("pdf/L2Res/vsPt/L2Res_vsPt_%04d_%04d_%s_%s.pdf",
		  int(1000.*eta1),int(1000.*eta2),cr,"c5"));
  

  // Step 6. Also draw final results into a giant canvas
  cx->cd(ieta);
  TH1D *h6 = tdrHist("h6","Rel. JES Data/MC",0.65,1.35);
  if      (eta<1.653) h6->GetYaxis()->SetRangeUser(0.8,1.2);
  else if (eta<2.964) h6->GetYaxis()->SetRangeUser(0.7,1.3);
  else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.3,1.3);
  h6->Draw();
  gPad->SetLogx();

  l->SetLineStyle(kDotted);
  l->DrawLine(ptmaxe,0.30,ptmaxe,1.35);
  l->SetLineStyle(kDashed);
  l->DrawLine(15.,1,3500.,1);

  // Draw full range at the back
  tdrDraw(hzrn,"Pz",kOpenSquare,kRed-9);
  tdrDraw(hgrn,"Pz",kOpenCircle,kBlue-9);
  tdrDraw(hdrn,"Pz",kOpenDiamond,kGray);
  tdrDraw(hprn,"Pz",kOpenDiamond,kOrange-9);
  tdrDraw(hjrn,"Pz",kOpenDiamond,kGreen-9);

  //f0->Draw("SAME");
  f1->Draw("SAME");
  //f2->Draw("SAME");
  f3->Draw("SAME");
  f4->Draw("SAME");
  f2->Draw("SAME");
  
  mg->Draw();

  if (ieta==1 || ieta==nxy) leg1->Draw("SAME");
  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));


  // Store results with uncertainty to output histogram
  for (int ipt = 1; ipt != h2jes->GetNbinsY()+1; ++ipt) {
    double pt = h2jes->GetYaxis()->GetBinCenter(ipt);
    double pt1 = h2jes->GetYaxis()->GetBinLowEdge(ipt);
    double emax1 = pt1*cosh(eta1);
    double jes4 = f4->Eval(pt); //quadlog+1/x
    double jes3 = f3->Eval(pt); //quadlog
    double jes2 = f2->Eval(pt); //loglin+1/x
    double jes1 = f2->Eval(pt); //loglin
    //double jes = jes4; // quadlog+1/x
    double jes = jes2; // loglin+1/x
    double ejes = sqrt(pow(jes1-jes,2) + pow(jes2-jes,2) + 
		       pow(jes3-jes,2) + pow(jes4-jes,2));
    if (emax1 < 13600.*0.5) {
      h2jes->SetBinContent(ieta, ipt, jes);
      h2jes->SetBinError(ieta, ipt, ejes); 
    }
  }
  hmin->SetBinContent(ieta, f2->Eval(10.));
  hmax->SetBinContent(ieta, f2->Eval(6800./cosh(eta1)));

  
  // Rename to avoid loop leakage and errors
  h1->SetName(Form("h1pt_%s_%d",cr,ieta));
  h2->SetName(Form("h2pt_%s_%d",cr,ieta));
  h3->SetName(Form("h3pt_%s_%d",cr,ieta));
  h4->SetName(Form("h4pt_%s_%d",cr,ieta));
  h5->SetName(Form("h5t_%s_%d",cr,ieta));
  h6->SetName(Form("h6t_%s_%d",cr,ieta));
  
  f0->SetName(Form("f0pt_%s_%d",cr,ieta));
  f1->SetName(Form("f1pt_%s_%d",cr,ieta));
  f2->SetName(Form("f2pt_%s_%d",cr,ieta));
  f3->SetName(Form("f3pt_%s_%d",cr,ieta));
  f4->SetName(Form("f4pt_%s_%d",cr,ieta));
  
  c1->SetName(Form("c1pt_%s_%d",cr,ieta));
  c2->SetName(Form("c2pt_%s_%d",cr,ieta));
  c3->SetName(Form("c3pt_%s_%d",cr,ieta));
  c4->SetName(Form("c4pt_%s_%d",cr,ieta));
  c5->SetName(Form("c5pt_%s_%d",cr,ieta));

  } // for ieta

  cx->SaveAs(Form("pdf/L2res/L2Res_AllEta_%s.pdf",cr));
  cx->SetName(Form("cx_%s",cr));

  // Step 7. Draw summary of final results in a single plot
  TH1D *hy = tdrHist("hy","JES",0.6,1.35,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *cy = tdrCanvas("cy",hy,8,33,kSquare);

  l->SetLineStyle(kDotted);
  l->DrawLine(1.3,0.6,1.3,1.1);
  l->DrawLine(2.5,0.6,2.5,1.1);
  l->DrawLine(2.964,0.6,2.964,1.35);
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);

  tdrDraw(hmin,"HIST",kNone,kMagenta+2,kSolid,-1,kNone,kMagenta+2);
  tdrDraw(hmax,"HIST",kNone,kBlack,kSolid,-1,kNone,kBlack);

  TLegend *legy0 = tdrLeg(0.20,0.20,0.45,0.20+2*0.045);
  legy0->AddEntry(hmax,"E < #sqrt{s}/2","FL");
  legy0->AddEntry(hmin,"p_{T} > 10 GeV","FL");
  
  TH1D *hy15, *hy30, *hy100, *hy300, *hy1000, *hy3000;
  hy100  = drawH2JES(h2jes,100., "HISTE][",kNone,kGreen+2);
  hy100->SetLineWidth(3);

  hy15   = drawH2JES(h2jes,15.,  "HISTE][",kNone,kMagenta+2);
  hy30   = drawH2JES(h2jes,30.,  "HISTE][",kNone,kBlue);
  hy300  = drawH2JES(h2jes,300., "HISTE][",kNone,kOrange+2);
  hy1000 = drawH2JES(h2jes,1000.,"HISTE][",kNone,kRed);
  hy3000 = drawH2JES(h2jes,3000.,"HISTE][",kNone,kBlack);

  TLegend *legy = tdrLeg(0.20,0.90-6*0.045,0.45,0.90);
  legy->AddEntry(hy15,  "p_{T} = 15 GeV",  "PLE");
  legy->AddEntry(hy30,  "p_{T} = 30 GeV",  "PLE");
  legy->AddEntry(hy100, "p_{T} = 100 GeV", "PLE");
  legy->AddEntry(hy300, "p_{T} = 300 GeV", "PLE");
  legy->AddEntry(hy1000,"p_{T} = 1000 GeV","PLE");
  legy->AddEntry(hy3000,"p_{T} = 3000 GeV","PLE");
  
  cy->SaveAs(Form("pdf/L2res/L2Res_Summary_%s.pdf",cr));
  cy->SetName(Form("cy_%s",cr));

  // Step 8. Print out text files
  ofstream ftxt(Form("textfiles/Summer23_L2ResClosure/Summer23Prompt23_Run%s_V1_DATA_L2Residual_AK4PFPuppi.txt",cr));
  ftxt << Form("{ 1 JetEta 1 JetPt 1./(%s) Correction L2Relative}",
	       vf1[0]->GetExpFormula().Data()) << endl;
  for (int ieta = p2d->GetNbinsX(); ieta != 0; --ieta) {
    double eta1 = -p2d->GetXaxis()->GetBinLowEdge(ieta+1);
    double eta2 = -p2d->GetXaxis()->GetBinLowEdge(ieta);
    TF1 *f1 = vf1[ieta-1];
    ftxt << Form("  %+1.3f %+1.3f  %d  %d %4d  ", eta1, eta2,
		 2 + f1->GetNpar(),  10, int(6800. / cosh(eta2)));
    for (int i = 0; i != f1->GetNpar(); ++i) {
      ftxt << Form(" %7.4f",f1->GetParameter(i));
    }
    ftxt << endl;
  }
  for (int ieta = 1; ieta != p2d->GetNbinsX()+1; ++ieta) {
    double eta1 = p2d->GetXaxis()->GetBinLowEdge(ieta);
    double eta2 = p2d->GetXaxis()->GetBinLowEdge(ieta+1);
    TF1 *f1 = vf1[ieta-1];
    ftxt << Form("  %+1.3f %+1.3f  %d  %d %4d  ", eta1, eta2,
		 2 + f1->GetNpar(),  10, int(6800. / cosh(eta1)));
    for (int i = 0; i != f1->GetNpar(); ++i) {
      ftxt << Form(" %7.4f",f1->GetParameter(i));
    }
    ftxt << endl;
  }

  
  // Loop over gamjet pT bins for plotting
  for (int ipt = 1; ipt != p2g->GetNbinsY()+1; ++ipt) {
    double ptmin = p2g->GetYaxis()->GetBinCenter(ipt);
    double ptmax = ptmin;
    double pt1 = p2g->GetYaxis()->GetBinLowEdge(ipt);
    double pt2 = p2g->GetYaxis()->GetBinLowEdge(ipt+1);
  // (No indent here for the resf of the loop, maybe function call later)

    
  // Step 1. Slice pT, draw response vs eta. No other manipulation yet
  TH1D *h1 = tdrHist("h1","JES",0.3,1.5,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,8,33,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0,1,5.2,1);

  TLegend *leg1 = tdrLeg(0.20,0.90,0.45,0.90);
  _leg = leg1;
  
  TProfile *pzm, *pgm, *pdm, *pjm, *ppm;
  pzm = drawEta(p2zm,ptmin,ptmax,"HISTE",kNone,kRed);
  pgm = drawEta(p2gm,ptmin,ptmax,"HISTE",kNone,kBlue);
  pjm = drawEta(p2jm,ptmin,ptmax,"HISTE",kNone,kGreen+2);
  ppm = drawEta(p2pm,ptmin,ptmax,"HISTE",kNone,kOrange+2);
  pdm = drawEta(p2dm,ptmin,ptmax,"HISTE",kNone,kBlack);

  TProfile *pz, *pg, *pd, *pj, *pp;
  pz = drawEta(p2z,ptmin,ptmax,"Pz",kFullSquare,kRed,"Z");
  pg = drawEta(p2g,ptmin,ptmax,"Pz",kFullCircle,kBlue,"#gamma");
  pj = drawEta(p2j,ptmin,ptmax,"Pz",kFullDiamond,kGreen+2,"Tag");
  pp = drawEta(p2p,ptmin,ptmax,"Pz",kFullDiamond,kOrange+2,"Probe");
  pd = drawEta(p2d,ptmin,ptmax,"Pz",kOpenDiamond,kBlack,"Dijet");

  c1->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c1"));
    

  // Step 2. Project profile to histogram, normalize by |eta|<1.3
  TH1D *h2 = tdrHist("h2","Rel. JES",0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c2 = tdrCanvas("c2",h2,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
  
  TH1D *hzm, *hgm, *hdm, *hjm, *hpm;
  hzm = drawNormEta(pzm,"HISTE",kNone,kRed);
  hgm = drawNormEta(pgm,"HISTE",kNone,kBlue);
  hjm = drawNormEta(pjm,"HISTE",kNone,kGreen+2);
  hpm = drawNormEta(ppm,"HISTE",kNone,kOrange+2);
  hdm = drawNormEta(pdm,"HISTE",kNone,kBlack);

  TH1D *hz, *hg, *hd, *hj, *hp;
  hz = drawNormEta(pz,"Pz",kFullSquare,kRed);
  hg = drawNormEta(pg,"Pz",kFullCircle,kBlue);
  hj = drawNormEta(pj,"Pz",kFullDiamond,kGreen+2);
  hp = drawNormEta(pp,"Pz",kFullDiamond,kOrange+2);
  hd = drawNormEta(pd,"Pz",kOpenDiamond,kBlack);

  c2->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c2"));
    
  
  // Step 3. Draw data/MC ratio before normalization
  TH1D *h3 = tdrHist("h3","JES Data/MC",0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c3 = tdrCanvas("c3",h3,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
  
  TH1D *hzr, *hgr, *hdr, *hjr, *hpr;
  hzr = drawRatio(pz->ProjectionX(),pzm,"Pz",kFullSquare,kRed);
  hgr = drawRatio(pg->ProjectionX(),pgm,"Pz",kFullCircle,kBlue);
  hjr = drawRatio(pj->ProjectionX(),pjm,"Pz",kFullDiamond,kGreen+2);
  hpr = drawRatio(pp->ProjectionX(),ppm,"Pz",kFullDiamond,kOrange+2);
  hdr = drawRatio(pd->ProjectionX(),pdm,"Pz",kOpenDiamond,kBlack);

  c3->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c3"));
  
  
  // Step 4. Draw data/MC ratio of normalized JES
  TH1D *h4 = tdrHist("h4","Rel. JES Data/MC",0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c4 = tdrCanvas("c4",h4,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
	     
  TH1D *hzrn, *hgrn, *hdrn, *hjrn, *hprn;
  hzrn = drawRatio(hz,hzm,"Pz",kFullSquare,kRed);
  hgrn = drawRatio(hg,hgm,"Pz",kFullCircle,kBlue);
  hjrn = drawRatio(hj,hjm,"Pz",kFullDiamond,kGreen+2);
  hprn = drawRatio(hp,hpm,"Pz",kFullDiamond,kOrange+2);
  hdrn = drawRatio(hd,hdm,"Pz",kOpenDiamond,kBlack);

  c4->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c4"));

  // Rename to avoid loop leakage and errors
  h1->SetName(Form("h1_%s_%d",cr,ipt));
  h2->SetName(Form("h2_%s_%d",cr,ipt));
  h3->SetName(Form("h3_%s_%d",cr,ipt));
  h4->SetName(Form("h4_%s_%d",cr,ipt));
  
  c1->SetName(Form("c1_%s_%d",cr,ipt));
  c2->SetName(Form("c2_%s_%d",cr,ipt));
  c3->SetName(Form("c3_%s_%d",cr,ipt));
  c4->SetName(Form("c4_%s_%d",cr,ipt));

  } // for ipt


  // Save results for easier book-keeping
  // Notes: HF needs correction to map from pT,tag bin to <pT,probe>
  if (fout) {
    fout->cd();
    h2jes->Write(Form("h2jes_%s",cr),TObject::kOverwrite);
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
} // L2Res
