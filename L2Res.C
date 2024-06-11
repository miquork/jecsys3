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

bool doClosure = false; // Do not undo L2L3Res for closure test

bool flattenHF = false;//true; // Const vs pT for HF
double flatHFptmin = 80;//60;//50;//80; // If flat, use only high pT
double flatHFetamin = 3.139;//2.964;
bool posOffHF = false;//true; // Positive offset for HF (2024BC,3/fb patch)
bool posOffEC2 = false;//true; // Positive offset for EC2 (2024BCD, 12.3/fb patch)
bool negLogHF = false;//true; // Positive offset for HF (2024BC,3/fb patch)

// Global variable for renaming histograms
string _run;


// Step 1. Slice 1D profile out of 2D in given range and draw it
TProfile* drawEta(TProfile2D *p2, double ptmin, double ptmax,
		  string draw, int marker, int color, string label="",
		  TProfile2D *p2x = 0, TProfile2D *p2xx = 0) {
  
  assert(p2);
  int iy1 = p2->GetYaxis()->FindBin(ptmin);
  int iy2 = p2->GetYaxis()->FindBin(ptmax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(), ptmin, ptmax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileX(Form("p%s_%s",id.c_str(),_run.c_str()),iy1,iy2);

  if (p2x) {
    TProfile *px = p2x->ProfileX(Form("px%s",id.c_str()),iy1,iy2);
    for (int i = 1; i != p->GetNbinsX()+1; ++i) {
      int j = px->GetXaxis()->FindBin(p->GetBinCenter(i));
      double mean_value = p->GetBinContent(i) * px->GetBinContent(j);
      //if (label=="Z") // inconsistent definition of p2res
      //mean_value = p->GetBinContent(i) / px->GetBinContent(j);
      //double calculated_error = p->GetBinError(i) * px->GetBinContent(j);
      double entries = p->GetBinEntries(i);
      p->SetBinEntries(i, entries);
      p->SetBinContent(i, entries * mean_value);
      //p->SetBinError(i, calculated_error);
    }
  }

  // Patch for re-reco. NB: not yet averaging over pT properly
  if (p2xx) {
    for (int i = 1; i != p->GetNbinsX()+1; ++i) {
      double eta = p->GetXaxis()->GetBinCenter(i);
      double pt = max(50., p2->GetYaxis()->GetBinCenter(iy1));
      int ieta = p2xx->GetYaxis()->FindBin(eta);
      int ipt = p2xx->GetXaxis()->FindBin(pt);
      double d = p2xx->GetBinContent(ipt,ieta); // 0.5*(B-A)/tag
      double k = 1 + 2*d;
      double mean_value = p->GetBinContent(i) * k;
      double entries = p->GetBinEntries(i);
      p->SetBinEntries(i, entries);
      p->SetBinContent(i, entries * mean_value);
    }
  }

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
		 string draw, int marker, int color, string label="",
		 TProfile2D *p2x = 0, TProfile2D *p2xx = 0) {
  
  assert(p2);
  int ix1 = p2->GetXaxis()->FindBin(etamin);
  int ix2 = p2->GetXaxis()->FindBin(etamax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(),
		   1000*etamin, 1000*etamax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileY(Form("p%s",id.c_str()),ix1,ix2);

  if (p2x) {
    TProfile *py = p2x->ProfileY(Form("px%s",id.c_str()),ix1,ix2);
    for (int i = 1; i != p->GetNbinsX()+1; ++i) {
      int j = py->GetXaxis()->FindBin(p->GetBinCenter(i));
      double mean_value = p->GetBinContent(i) * py->GetBinContent(j);
      //if (label=="Z") // inconsistent definition of p2res
      //mean_value = p->GetBinContent(i) / py->GetBinContent(j);
      //double calculated_error = p->GetBinError(i) * py->GetBinContent(j);
      double entries = p->GetBinEntries(i);
      p->SetBinEntries(i, entries);
      p->SetBinContent(i, entries * mean_value);
      //p->SetBinError(i, calculated_error);
    }
  }

  // Patch for re-reco. NB: not yet averaging over |eta| properly
  if (p2xx) {
    for (int i = 1; i != p->GetNbinsX()+1; ++i) {
      double pt = max(50., p->GetXaxis()->GetBinCenter(i));
      double eta = p2->GetXaxis()->GetBinCenter(ix1);
      int ieta = p2xx->GetYaxis()->FindBin(eta);
      int ipt = p2xx->GetXaxis()->FindBin(pt);
      double d = p2xx->GetBinContent(ipt,ieta); // 0.5*(B-A)/tag
      double k = 1 + 2*d;
      double mean_value = p->GetBinContent(i) * k;
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

// Step 2. Project profile to histogram, normalize by |eta|<1.3 and draw it
TH1D *drawNormEta(TProfile *p, string draw, int marker, int color) {

  assert(p);
  string id = Form("%s_%s_%d_%d",p->GetName(),draw.c_str(),marker,color);
  TH1D *h = p->ProjectionX(Form("h%s",id.c_str()));
  TF1 *f1 = new TF1(Form("f1%s_%s",id.c_str(),_run.c_str()),"[0]",0,1.305);
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
  TH1D *h = p->ProjectionX(Form("h%s_%s",id.c_str(),_run.c_str()));
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
  TH1D *hr = (TH1D*)h->Clone(Form("hr_%s_%s",id.c_str(),_run.c_str()));
  hr->Divide(hm);

  tdrDraw(hr,draw,marker,color,kSolid,-1,kNone);
  
  return hr;
}

// Step 5. Clean up data (for fitting) and draw it
TH1D *drawCleaned(TH1D *h, double eta, string data, string draw,
		  int marker, int color) {
  assert(h);
  TH1D *hc = (TH1D*)h->Clone(Form("hc_%s_%s",h->GetName(),_run.c_str()));
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
    // New L2Res eta binning (same as L2Relative)
    //{0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};

    if (TString(_run).Contains("2024")) {
	
    // Prompt24
    if (data=="Z") {
      if      (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.500) keep = true;
      else if (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.964) keep = true;
      else if (ptmin>=30. && ptmax<300 && eta>2.964  && eta<3.489) keep = true;
      else if (ptmin>=15. && ptmax<200 && eta>3.489  && eta<3.839) keep = true;
      else if (ptmin>=15. && ptmax<100 && eta>3.839  && eta<5.191) keep = true;
      // Additional veto
      if (eta>2.322 && eta<3.139) keep = false;
      // Additional veto for 2024B, 3/fb
      //if (eta>4.538 && eta<4.716 && ptmax<30) keep = false;
      // Additional veto for 2024B,C,D, 12.3/fb
      if (eta>3.139 && ptmax<70) keep = false;
    }
    if (data=="G") { // HLT_Photon50EB_TightId_TightIso, HLT_Photon30EB...
      if      (ptmin>=30. && ptmax<1300. && emax<2500.) keep = true;
      else if (ptmin>=30. && ptmax<200.  && emax<0.5*sqrts) keep = true;
      // Additional veto for 2024C 0.7/fb golden special
      if (fabs(eta)>2.322 && eta<2.964 && ptmin<60) keep = false;
      // Additional veto for 2024BC 3.3/fb golden closure (errors bad)
      if (doClosure && ptmin>500) keep = false;
    }
    if (data=="J") {
      if       (ptmin>=15. && ptmax<2500. && eta<1.044) keep = true;
      else if  (ptmin>=15. && ptmax<2000. && eta<2.500) keep = true;
      else if  (ptmin>=15. && ptmax<967.  && eta<2.964) keep = true;
      else if  (ptmin>=15. && ptmax<846.  && eta<3.839) keep = true;
      else if  (ptmin>=15. && ptmax<460.  && eta<5.191) keep = true;
      // Additional veto for 2024B, 3/fb
      //if (eta>4.538 && eta<4.716 && ptmax<30) keep = false;
      // Additional veto for 2024BC 3.3/fb golden closure (bad errors?)
      if (doClosure && emax>0.5*sqrts) keep = false;
    }
    if (data=="P") {
      if (ptmin>=15. && ptmax<=110. && eta<2.322) keep = true;
    }
    if (data=="D") {
      if (ptmin>=59. && ptmax<=110. && emax<3100. && eta<2.853) keep = true;
      // Additional veto
      if (eta>2.322) keep = false;
    }
    } // Prompt24
    else {
    // Summer23
    if (data=="Z") {
      if      (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.500) keep = true;
      else if (ptmin>=15. && ptmax<550 && emax<1550. && eta<2.964) keep = true;
      else if (ptmin>=59. && ptmax<200 && eta>2.964  && eta<3.839) keep = true;
    }
    if (data=="G") {
      if      (ptmin>=110. && ptmax<1300. && emax<2500.) keep = true;
      else if (ptmin>=110. && ptmax<200.  && emax<0.5*sqrts) keep = true;
    }
    if (data=="J") {
      if       (ptmin>=15. && ptmax<2500. && eta<1.044) keep = true;
      else if  (ptmin>=15. && ptmax<2000. && eta<2.500) keep = true;
      else if  (ptmin>=15. && ptmax<967.  && eta<2.964) keep = true;
      else if  (ptmin>=15. && ptmax<846.  && eta<3.839) keep = true;
      else if  (ptmin>=15. && ptmax<460.  && eta<5.191) keep = true;
    }
    if (data=="P") {
      if (ptmin>=15. && ptmax<=110. && eta<2.322) keep = true;
    }
    if (data=="D") {
      if (ptmin>=59. && ptmax<=110. && emax<3100. && eta<2.853) keep = true;
    }
    } // !Prompt24
    
    // Check that no points out of bound
    if (h->GetBinContent(i)>1.3 || h->GetBinContent(i)<0.3) {
      keep = false;
    }
    // Remove BPIX bad point, HF bad point
    if (data=="J" || data=="P" || data=="D") {
      if (eta>1.653 && eta < 1.930 && pt>86 && pt<110) keep = false;
      if (eta>4.018 && eta < 4.583 && pt>279 && pt<302) keep = false;
    }
    // Remove low pT for flat HF fits
    if (flattenHF && eta>flatHFetamin && pt<flatHFptmin) keep = false;

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
  //gROOT->ProcessLine(".! mkdir pdf/L2Res/vsEta");
  //gROOT->ProcessLine(".! mkdir pdf/L2Res/vsPt");

  map<string, string> mlum;
  mlum["2024B"] = "0.13 fb^{-1}";
  mlum["2024C"] = "7.5 fb^{-1}";
  mlum["2024D"] = "8.3 fb^{-1}"; //partial, 4.6->7.9->8.3
  mlum["2024BC"] = "7.7 fb^{-1}";
  mlum["2024BCD"] = "16.0 fb^{-1}"; // partial, 12.3->15.6->16.0
  mlum["2024CP"] = "7.5 fb^{-1}";
  mlum["2024CR"] = "7.5 fb^{-1}";
  mlum["2024CS"] = "7.5 fb^{-1}";
  mlum["2024E"] = "11.0 fb^{-1}"; // partial, 9.3->11.0
  
  //string vrun[] = {"2023Cv123","2023Cv4","2023D"};
  //string vrun[] = {"2022CD","2022E","2022F","2022G"};
  //string vrun[] = {"2024B","2024C","2024D"};
  //string vrun[] = {"2024BCD","2024CR"};
  //string vrun[] = {"2024E","2024C","2024D","2024BCD"};
  //string vrun[] = {"2024E"};
  string vrun[] = {"2024CS"};
  //string vrun[] = {"2024BCD","2024E","2024C","2024CR","2024CS"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  //string vmc[] = {"Summer23","Summer23","Summer23BPIX"};
  //string vmc[] = {"Summer22","Summer22EE","Summer22EE","Summer22EE"};
  //string vmc[] = {"Summer22","Summer22","Summer22","Summer22"}; // 22EE buggy?
  //string vmc[] = {"Summer23BPix","Summer23BPix","Summer23BPix"};
  //string vmc[] = {"Summer23BPix","Summer23BPix"};
  //string vmc[] = {"Summer23BPix","Summer23BPix","Summer23BPix","Summer23BPix"};
  string vmc[] = {"Summer23BPix"};
  //string vmc[] = {"Summer23BPix","Summer23BPix","Summer23BPix",
  //		  "Summer23BPix","Summer23BPix"};
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
    
    gROOT->ProcessLine(Form(".! mkdir pdf/L2Res/%s",cr));
    gROOT->ProcessLine(Form(".! mkdir pdf/L2Res/%s/vsEta",cr));
    gROOT->ProcessLine(Form(".! mkdir pdf/L2Res/%s/vsPt",cr));

  // (No indent here for the resf of the loop, maybe function call later)


  // Mikko's temporary code to create L2Res.root file
  // Load Z+jet
  TFile *fz(0), *fsz(0);
  TProfile2D *p2zxx(0);
  if (run=="2024CS") {
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ"); // June 5 hybrid, 15.6/fb
    // HCALDI/ECALRATIO ~ HCALDI/Prompt for Z+jet
    fsz = new TFile("rootfiles/compareLite_2024CR.root","READ");
    assert(fsz && !fsz->IsZombie());
    p2zxx = (TProfile2D*)fsz->Get("2D/p2d_tp"); assert(p2zxx);
  }
  else if (run=="2024CR") {
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ"); // June 5 hybrid, 15.6/fb
  }
  else if (TString(cr).Contains("2024")) {
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v78.root",cr),"READ"); // v77: Golden 2024B, DCSOnly 2024C
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v78golden.root",cr),"READ"); // 0.74/fb
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v79golden.root",cr),"READ"); // 3/fb
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v80golden.root",cr),"READ"); // 3/fb closure
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v81.root",cr),"READ"); // May 16 golden, 12.3/fb
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v82.root",cr),"READ"); // May 27 golden + DCSOnly, X/fb
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root",cr),"READ"); // June 5 golden + DCSOnly, 15.6/fb
  }
  else if (TString(cr).Contains("2023")) {
    //fz = new TFile(Form("rootfiles/Summer23_L2ResOnly/jme_bplusZ_%s_Zmm_sync_v70.root",cr),"READ"); // Summer23 L2Res_V1
    //fz = new TFile(Form("rootfiles/Summer23_L2L3Res/jme_bplusZ_%s_Zmm_sync_v72noPUrw.root",cr),"READ"); // MC botched?
    //fz = new TFile(Form("rootfiles/Summer23_L2L3Res/jme_bplusZ_%s_Zmm_sync_v73smearoff.root",cr),"READ");
    fz = new TFile(Form("rootfiles/Summer23_L2L3Res/jme_bplusZ_%s_Zmm_sync_v77.root",cr),"READ");
  }
  else if (TString(cr).Contains("2022")) {
    //fz = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v66.root",cr),"READ");
    if (run=="2022F" || run=="2022G" || run=="2022FG")
      fz = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v78.root",cr),"READ"); // 22Sep2023 instead of 19Dec2023
    else
      fz = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v76.root",cr),"READ");
  }
  assert(fz && !fz->IsZombie());

  TDirectory *dz = fz->GetDirectory("data/l2res"); assert(dz);
  TDirectory *dzm = fz->GetDirectory("mc/l2res"); assert(dzm);
  
  // Load G+jet
  TFile *fg(0), *fgm(0);
  if (run=="2024CS") {
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALR-HCALDI_w29.root","2024C"),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
  }
  else if (run=="2024CR") {
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w29.root","2024C"),"READ"); // w26->w29
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
  }
  else if (TString(cr).Contains("2024")) {
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w12.root",cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w13.root",cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w14.root",cr),"READ"); // Golden JSON 0.74/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w16.root",cr),"READ"); // Golden JSON 3/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w18.root",cr),"READ"); // Golden JSON 3/fb closure
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w22.root",cr),"READ"); // May 16 golden, 12.3/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w27.root",cr),"READ"); // DCSOnly
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w29.root",cr),"READ"); // June 5 hybrid, 15.6/fb
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w12.root","READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
  }
  else if (TString(cr).Contains("2023")) {
    fg = new TFile(Form("rootfiles/Summer23_L2L3Res/GamHistosFill_data_%s_w8.root",cr),"READ");
    fgm = new TFile(run=="2023D" ? "rootfiles/Summer23_L2L3Res/GamHistosFill_mc_2023P8-BPix_w8.root" : "rootfiles/Summer23_L2L3Res/GamHistosFill_mc_2023P8_w8.root","READ");
  }
  else if (TString(cr).Contains("2022")) {
    fg = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s_v32.root",cr),"READ");
    const char *cmg = TString(cm).ReplaceAll("Summer","20").Data();
    fgm = new TFile(Form("../gamjet/rootfiles/GamHistosFill_mc_%sP8_v32.root",cmg),"READ");
  }
  assert(fg && !fg->IsZombie());
  assert(fgm && !fgm->IsZombie());
  
  TDirectory *dg = fg->GetDirectory("Gamjet2"); assert(dg);
  TDirectory *dgm = fgm->GetDirectory("Gamjet2"); assert(dgm);

  // Load dijet
  TFile *fd(0), *fdm(0);
  //if (run=="2024E") {
  //fd = new TFile(Form("rootfiles/Prompt2024/v68_2024/jmenano_data_cmb_%s_JME_v68_2024.root",cr),"READ");
  //fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  //}
  if (run=="2024CS") {
    fd = new TFile(Form("rootfiles/Prompt2024/v77_2024/jmenano_data_cmb_%sS_JME_v77_2024.root","2024C"),"READ"); // was v76/*2024Crs*
    fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (run=="2024CR") {
    fd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_ECAL_JME_v76_2024.root","2024C"),"READ"); // One half v66->v76 full
    fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (TString(cr).Contains("2024")) {
    //fd = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_data_cmb_%s_v39_2024_Prompt_etabin_DCSOnly.root",cr),"READ");
    //fd = new TFile(Form("rootfiles/Prompt2024/jmenano_data_cmb_%s_JME_v39_2024_Prompt_Golden_29April.root",cr),"READ"); // golden 0.74/fb
    //fd = new TFile(Form("rootfiles/Prompt2024/v41_2024_Golden/jmenano_data_cmb_%s_JME_v41_2024_Golden.root",cr),"READ"); // golden 3/fb
    //fd = new TFile(Form("rootfiles/Prompt2024/v43_2024_Golden/jmenano_data_cmb_%s_JME_v43_2024_Golden.root",cr),"READ"); // golden 3/fb closure
    //fd = new TFile(Form("rootfiles/Prompt2024/v50_2024/jmenano_data_cmb_%s_JME_v50_2024.root",cr),"READ"); // May 16 golden, 12.3/fb
    //fd = new TFile(Form("rootfiles/Prompt2024/v67_2024/jmenano_data_cmb_%s_JME_v67_2024.root",cr),"READ"); // May 27 golden, X/fb
    fd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root",cr),"READ"); // June 5 hybrid, 15.6/fb
    //if (run=="2024B")
    //fd = new TFile(Form(" rootfiles/Prompt2024/v39_2024_Prompt_etabin_SFD_Golden/jmenano_data_cmb_%s_v39_2024_Prompt_etabin_SFD_Golden.root",cr),"READ");
    fdm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_out_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (TString(cr).Contains("2023")) {
    fd = new TFile(Form("rootfiles/Summer23_L2L3ResJERSF/v39_2023_etabin_SFv2/jmenano_data_cmb_%s_JME_v39_2023_etabin_SFv2.root",cr),"READ");
    fdm = new TFile(Form("rootfiles/Summer23_L2L3ResJERSF/v39_2023_etabin_SFv2/jmenano_mc_out_Summer23MG%s_v39_2023_etabin_SFv2.root",run=="2023D" ? "BPix" : "_Cv123"),"READ");
  }
  else if (TString(cr).Contains("2022")) {
    //fd = new TFile(Form("rootfiles/Summer22_L2L3ResOnly/v39_2022_JV_noSJets_noRwPU/jmenano_data_cmb_%s_JME_v39_2022_JV_noSJets_noRwPU.root",cr),"READ");
    //fdm = new TFile(Form("rootfiles/Summer22_L2L3ResOnly/v39_2022_JV_noSJets_noRwPU/jmenano_mc_out_%sMG_full_v39_2022_JV_noSJets_noRwPU.root",cm),"READ");
    fd = new TFile(Form("rootfiles/Summer22_L2L3ResJERSF/v39_2022_etabin_SFv/jmenano_data_cmb_%s_JME_v39_2022_etabin_SFv.root",cr),"READ");
    fdm = new TFile(Form("rootfiles/Summer22_L2L3ResJERSF/v39_2022_etabin_SFv/jmenano_mc_out_%sMG_full_v39_2022_etabin_SFv.root",cm),"READ");
  }
  assert(fd && !fd->IsZombie());
  assert(fdm && !fdm->IsZombie());

  TDirectory *dd = fd->GetDirectory("Dijet2"); assert(dd);
  //TDirectory *ddm = fdm->GetDirectory("Dijet2");
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
  // Step 0. Extract barrel reference to normalize it out later
  TProfile *p13zm(0), *p13gm(0), *p13dm(0), *p13jm(0), *p13pm(0);
  TProfile *p13z(0), *p13g(0), *p13d(0), *p13j(0), *p13p(0);

  TH1D *h13 = tdrHist("h13","JES",0.65,1.85);//0.65,1.45);
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
  
  p13z = drawPt(p2z,etamin,etamax,"Pz",kFullSquare,kRed,"Z",p2zx,p2zxx);
  p13g = drawPt(p2g,etamin,etamax,"Pz",kFullCircle,kBlue,"#gamma",p2gx);
  p13j = drawPt(p2j,etamin,etamax,"Pz",kFullDiamond,kGreen+2,"Tag",p2jx);
  p13p = drawPt(p2p,etamin,etamax,"Pz",kFullDiamond,kOrange+2,"Probe",p2px);
  p13d = drawPt(p2d,etamin,etamax,"Pz",kOpenDiamond,kBlack,"Dijet",p2dx);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c13->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		   cr,int(1000.*eta1),int(1000.*eta2),cr,cc,"c13"));
  
  // Create giant canvas for all eta bins (7*3=21; more in the future)
  int nxy = p2d->GetNbinsX();
  //TCanvas *cx = new TCanvas("cx","cx",9*250,5*250);
  //TCanvas *cx = new TCanvas("cx","cx",9*400,5*400);
  TCanvas *cx = new TCanvas(Form("cx_%s",cr),"cx",9*300,5*300);
  cx->Divide(9,5,0,0);
  TH2D *h2jes = p2d->ProjectionXY(Form("h2jes_%s",cr)); h2jes->Reset();
  
  // Loop over the ieta bins
  vector<TF1*> vf1(p2d->GetNbinsX()+1);
  TH1D *hmin = p2d->ProjectionX(Form("hmin_%s",cr)); hmin->Reset();
  TH1D *hmax = p2d->ProjectionX(Form("hmax_%s",cr)); hmax->Reset();
  TH1D *hchi2 = p2d->ProjectionX(Form("hchi2_%s",cr)); hchi2->Reset();
  for (int ieta = 1; ieta != p2d->GetNbinsX()+1; ++ieta) {
    //int ieta = p2d->GetXaxis()->FindBin(2.7);
    double etamin = p2d->GetXaxis()->GetBinCenter(ieta);
    double etamax = etamin;
    double eta1 = p2d->GetXaxis()->GetBinLowEdge(ieta);
    double eta2 = p2d->GetXaxis()->GetBinLowEdge(ieta+1);
  // (No indent here for the resf of the loop, maybe function call later)
    
  TH1D *h1 = tdrHist(Form("h1_%s",cr),"JES",0.65,1.85);
  lumi_136TeV = Form("%s - %s%s",cr,cm,cl);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1_%s",cr),h1,8,33,kSquare);
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
  pz = drawPt(p2z,etamin,etamax,"Pz",kFullSquare,kRed,"Z",p2zx,p2zxx);
  pg = drawPt(p2g,etamin,etamax,"Pz",kFullCircle,kBlue,"#gamma",p2gx);
  pj = drawPt(p2j,etamin,etamax,"Pz",kFullDiamond,kGreen+2,"Tag",p2jx);
  pp = drawPt(p2p,etamin,etamax,"Pz",kFullDiamond,kOrange+2,"Probe",p2px);
  pd = drawPt(p2d,etamin,etamax,"Pz",kOpenDiamond,kBlack,"Dijet",p2dx);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c1->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		  cr,int(1000.*eta1),int(1000.*eta2),cr,cc,"c1"));

  
  // Step 2. Project profile to histogram, normalize by |eta|<1.3
  TH1D *h2 = tdrHist(Form("h2_%s",cr),"Rel. JES",0.65,1.85);
  TCanvas *c2 = tdrCanvas(Form("c2_%s",cr),h2,8,33,kSquare);
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

  c2->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		  cr,int(1000.*eta1),int(1000*eta2),cr,cc,"c2"));

  // Step 3. Draw data/MC ratio before normalization
  TH1D *h3 = tdrHist(Form("h3_%s",cr),"JES Data/MC",0.5,1.3);//0.80,1.15);
  TCanvas *c3 = tdrCanvas(Form("c3_%s",cr),h3,8,33,kSquare);
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
  
  c3->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		  cr,int(1000.*eta1),int(1000.*eta2),cr,cc,"c3"));

  // Step 4. Draw data/MC ratio of normalized JES
  TH1D *h4 = tdrHist(Form("h4_%s",cr),"Rel. JES Data/MC",0.50,1.35);
  TCanvas *c4 = tdrCanvas(Form("c4_%s",cr),h4,8,33,kSquare);
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

  c4->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		  cr,int(1000.*eta1),int(1000.*eta2),cr,cc,"c4"));


  // Step 5. Curate data and do final fit of response
  // Options: flat, log-linear, quadratic, +1/x
  TH1D *h5 = tdrHist(Form("h5_%s",cr),"Rel. JES Data/MC",0.50,1.35);
  TCanvas *c5 = tdrCanvas(Form("c5_%s",cr),h5,8,33,kSquare);
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

  TMultiGraph *mg = new TMultiGraph(Form("mg_%s",cr),"mg");
  if (fitZ) mg->Add(cleanGraph(new TGraphErrors(hzrf)),"SAMEP");
  if (fitG) mg->Add(cleanGraph(new TGraphErrors(hgrf)),"SAMEP");
  if (fitD) mg->Add(cleanGraph(new TGraphErrors(hdrf)),"SAMEP");
  if (fitP) mg->Add(cleanGraph(new TGraphErrors(hprf)),"SAMEP");
  if (fitJ) mg->Add(cleanGraph(new TGraphErrors(hjrf)),"SAMEP");
  mg->Draw("Pz");//"SAME Pz");

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
  // Reference JES
  TF1 *fref = new TF1(Form("fref_%d_%s",ieta,cr),
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
  if (eta>2.964 && posOffHF && !doClosure) {
    fref->SetParLimits(2,0.3,0.5); // 2024BC,3/fb special
  }
  if (eta>2.5 && eta<2.650 && posOffEC2 && !doClosure) {
    fref->SetParLimits(2,0.3,0.5); // 2024BCD, 12.3/fb special
  }
  //if (eta>3.139 && negLogHF) {
  if (eta>4.538 && eta<4.716 && negLogHF && !doClosure) {
    fref->SetParLimits(1,-1,-0.1); // 2024BC,3/fb special
  }
  if (eta>flatHFetamin && flattenHF) {
    fref->FixParameter(1,0.);
    fref->FixParameter(2,0.);
  }
  mg->Fit(fref,"QRN");

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
  // Bonus2: constrain HF |eta|>4.013 to constant
  if (eta>4 && false) {
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
  
  f0->Draw("SAME"); f0->SetLineColor(kMagenta+2); //f0->SetLineStyle(kDashed);
  f1->Draw("SAME"); f1->SetLineColor(kBlue);
  f2->Draw("SAME"); f2->SetLineColor(kGreen+1);
  f3->Draw("SAME"); f3->SetLineColor(kOrange+2);
  f4->Draw("SAME"); f4->SetLineColor(kRed);

  fref->Draw("SAME"); fref->SetLineColor(kBlack);

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  c5->SaveAs(Form("pdf/L2Res/%s/vsPt/L2Res_vsPt_%04d_%04d_%s%s_%s.pdf",
		  cr,int(1000.*eta1),int(1000.*eta2),cr,cc,"c5"));
  

  // Step 6. Also draw final results into a giant canvas
  cx->cd(ieta);
  double eps = 1e-4;
  TH1D *h6 = tdrHist(Form("h6_%s",cr),"Rel. JES Data/MC",0.65+eps,1.35-eps);
  //if      (eta<1.653) h6->GetYaxis()->SetRangeUser(0.8+eps,1.2-eps);
  //else if (eta<2.964) h6->GetYaxis()->SetRangeUser(0.7+eps,1.3-eps);
  //else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.3+eps,1.3-eps);
  if      (eta<1.466) h6->GetYaxis()->SetRangeUser(0.8+eps,1.2-eps);
  else if (eta<2.650) h6->GetYaxis()->SetRangeUser(0.65+eps,1.3-eps);
  else if (eta<4.013) h6->GetYaxis()->SetRangeUser(0.55+eps,1.35-eps);
  else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.30+eps,1.35-eps);
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

  if (eta>flatHFetamin && flattenHF) {
    f0->Draw("SAME");
  }
  f1->Draw("SAME");
  f2->Draw("SAME");
  f3->Draw("SAME");
  f4->Draw("SAME");
  //f2->Draw("SAME");

  fref->Draw("SAME");
  
  mg->Draw("Pz");//"SAME Pz");

  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  //if (ieta==1) leg1->Draw("SAME");
  //if (ieta==nxy) { cx->cd(ieta+1); leg1->Draw("SAME"); }
  if (ieta==nxy) {
    cx->cd(ieta+1);
    leg1->SetTextSize(2.5*0.045);
    leg1->SetY1NDC(leg1->GetY2NDC()-2.5*fabs(leg1->GetY2NDC()-leg1->GetY1NDC()));
    leg1->Draw("SAME");
    //leg1->SetTextSize(0.045);
    //leg1->SetY2(leg1->GetY1()-0.5*fabs(leg1->GetY1()-leg1->GetY2()));
    cx->cd(ieta+2);
    double siz = tex->GetTextSize();
    tex->SetTextSize(siz*2.5);
    //tex->DrawLatex(0.20,0.50,Form("%s - %s%s",cr,cm,cl));
    tex->DrawLatex(0.05,0.80,Form("%s vs",cr));
    tex->DrawLatex(0.05,0.65,cm);
    tex->DrawLatex(0.05,0.50,mlum[run].c_str());
    tex->SetTextSize(siz);
  }


  // Store results with uncertainty to output histogram
  for (int ipt = 1; ipt != h2jes->GetNbinsY()+1; ++ipt) {
    double pt = h2jes->GetYaxis()->GetBinCenter(ipt);
    double pt1 = h2jes->GetYaxis()->GetBinLowEdge(ipt);
    double emax1 = pt1*cosh(eta1);
    double jes4 = f4->Eval(pt); //quadlog+1/x
    double jes3 = f3->Eval(pt); //quadlog
    double jes2 = f2->Eval(pt); //loglin+1/x
    double jes1 = f1->Eval(pt); //loglin
    double jes0 = f0->Eval(pt); //const
    double jesref = fref->Eval(pt); //loglin+1/x, loglin or const (vs eta)
    //double jes = jes4; // quadlog+1/x
    //double jes = jes2; // loglin+1/x
    double jes = jesref; // loglin+1/x
    double ejes = sqrt(pow(jes1-jes,2) + pow(jes2-jes,2) + 
		       pow(jes3-jes,2) + pow(jes4-jes,2));
    if (eta>flatHFetamin && flattenHF) {
      ejes = sqrt(pow(jes0-jes,2) + pow(jes1-jes,2));
    }
    if (emax1 < 13600.*0.5) {
      h2jes->SetBinContent(ieta, ipt, jes);
      h2jes->SetBinError(ieta, ipt, ejes); 
    }
  }
  //hmin->SetBinContent(ieta, f2->Eval(10.));
  //hmax->SetBinContent(ieta, f2->Eval(6800./cosh(eta1)));
  hmin->SetBinContent(ieta, fref->Eval(10.));
  hmax->SetBinContent(ieta, fref->Eval(6800./cosh(eta1)));
  hchi2->SetBinContent(ieta, fref->GetChisquare() / max(fref->GetNDF(),1));
  hchi2->SetBinError(ieta, 1. / sqrt(max(fref->GetNDF(),1)));
  
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

  cx->SaveAs(Form("pdf/L2res/L2Res_AllEta_%s%s.pdf",cr,cc));
  cx->SetName(Form("cx_%s",cr));

  // Step 7. Draw summary of final results in a single plot
  TH1D *hy = tdrHist(Form("hy_%s",cr),"JES",0.5,1.45,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s%s",cr,cm,cl);
  extraText = "Private";
  TCanvas *cy = tdrCanvas(Form("cy_%s",cr),hy,8,33,kSquare);

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
  
  cy->SaveAs(Form("pdf/L2res/L2Res_Summary_%s%s.pdf",cr,cc));
  cy->SetName(Form("cy_%s",cr));

  // Step 8. Print out chi2/NDF vs eta
  TH1D *hc = tdrHist(Form("hc_%s",cr),"#chi^{2} / NDF",0.,10.,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s%s",cr,cm,cl);
  extraText = "Private";
  TCanvas *cchi2 = tdrCanvas(Form("cchi2_%s",cr),hc,8,33,kSquare);

  l->SetLineColor(kGray+2);
  l->SetLineStyle(kDotted);
  l->DrawLine(1.3,0.,1.3,10);
  l->DrawLine(2.5,0.,2.5,10);
  l->DrawLine(2.964,0.,2.964,10);
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);

  tdrDraw(hchi2,"HE",kNone,kBlack,kSolid,-1,kNone,0);

  cchi2->SaveAs(Form("pdf/L2res/L2Res_Chi2_%s%s.pdf",cr,cc));
  cchi2->SetName(Form("cc_%s",cr));
  
  
  // Step 9. Print out text files
  //ofstream ftxt(Form("textfiles/Summer23_L2ResClosure/Summer23Prompt23_Run%s_V1_DATA_L2Residual_AK4PFPuppi.txt",cr));
  //ofstream ftxt(Form("textfiles/Prompt24/Prompt24_Run%s_V1M_DATA_L2Residual_AK4PFPuppi.txt",cr));
  //ofstream ftxt(Form("textfiles/Prompt24/Prompt24_Run%s_V2M_DATA_L2Residual_AK4PFPuppi.txt",cr));
  //ofstream ftxt(Form("textfiles/Prompt24/Prompt24_Run%s_V3M_DATA_L2Residual_AK4PFPuppi.txt",cr));
  ofstream ftxt(Form("textfiles/Prompt24/Prompt24_Run%s_V4M_DATA_L2Residual_AK4PFPuppi.txt",cr));
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
  TH1D *h1 = tdrHist(Form("h1p_%s",cr),"JES",0.3,1.5,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s%s",cr,cm,cl);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas(Form("c1p_%s",cr),h1,8,33,kSquare);

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
  pz = drawEta(p2z,ptmin,ptmax,"Pz",kFullSquare,kRed,"Z",p2zx,p2zxx);
  pg = drawEta(p2g,ptmin,ptmax,"Pz",kFullCircle,kBlue,"#gamma",p2gx);
  pj = drawEta(p2j,ptmin,ptmax,"Pz",kFullDiamond,kGreen+2,"Tag",p2jx);
  pp = drawEta(p2p,ptmin,ptmax,"Pz",kFullDiamond,kOrange+2,"Probe",p2px);
  pd = drawEta(p2d,ptmin,ptmax,"Pz",kOpenDiamond,kBlack,"Dijet",p2dx);

  c1->SaveAs(Form("pdf/L2Res/%s/vsEta/L2Res_vsEta_%04d_%04d_%s%s_%s.pdf",
		  cr,int(pt1),int(pt2),cr,cc,"c1"));
    

  // Step 2. Project profile to histogram, normalize by |eta|<1.3
  TH1D *h2 = tdrHist(Form("h2p_%s",cr),"Rel. JES",0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c2 = tdrCanvas(Form("c2p_%s",cr),h2,8,33,kSquare);

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

  c2->SaveAs(Form("pdf/L2Res/%s/vsEta/L2Res_vsEta_%04d_%04d_%s%s_%s.pdf",
		  cr,int(pt1),int(pt2),cr,cc,"c2"));
    
  
  // Step 3. Draw data/MC ratio before normalization
  TH1D *h3 = tdrHist(Form("h3p_%s",cr),"JES Data/MC",0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c3 = tdrCanvas(Form("c3p_%s",cr),h3,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
  
  TH1D *hzr, *hgr, *hdr, *hjr, *hpr;
  hzr = drawRatio(pz->ProjectionX(),pzm,"Pz",kFullSquare,kRed);
  hgr = drawRatio(pg->ProjectionX(),pgm,"Pz",kFullCircle,kBlue);
  hjr = drawRatio(pj->ProjectionX(),pjm,"Pz",kFullDiamond,kGreen+2);
  hpr = drawRatio(pp->ProjectionX(),ppm,"Pz",kFullDiamond,kOrange+2);
  hdr = drawRatio(pd->ProjectionX(),pdm,"Pz",kOpenDiamond,kBlack);

  c3->SaveAs(Form("pdf/L2Res/%s/vsEta/L2Res_vsEta_%04d_%04d_%s%s_%s.pdf",
		  cr,int(pt1),int(pt2),cr,cc,"c3"));
  
  
  // Step 4. Draw data/MC ratio of normalized JES
  TH1D *h4 = tdrHist(Form("h4p_%s",cr),"Rel. JES Data/MC",
		     0.3,1.5,"|#eta|",0,5.2);
  TCanvas *c4 = tdrCanvas(Form("c4p_%s",cr),h4,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
	     
  TH1D *hzrn, *hgrn, *hdrn, *hjrn, *hprn;
  hzrn = drawRatio(hz,hzm,"Pz",kFullSquare,kRed);
  hgrn = drawRatio(hg,hgm,"Pz",kFullCircle,kBlue);
  hjrn = drawRatio(hj,hjm,"Pz",kFullDiamond,kGreen+2);
  hprn = drawRatio(hp,hpm,"Pz",kFullDiamond,kOrange+2);
  hdrn = drawRatio(hd,hdm,"Pz",kOpenDiamond,kBlack);

  c4->SaveAs(Form("pdf/L2Res/%s/vsEta/L2Res_vsEta_%04d_%04d_%s%s_%s.pdf",
		  cr,int(pt1),int(pt2),cr,cc,"c4"));

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
