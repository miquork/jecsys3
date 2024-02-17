// Purpose: derive JER SF with MPFX method
//          combine multiple channels, start from TProfile2D *p2m0
#include "TFile.h"
#include "TProfile2D.h"

#include "tdrstyle_mod22.C"


// Calculate JER with MPFX. If JER13 not given, calculate it
TH1D *getJER(TProfile2D *p2s, TProfile2D *p2x, TH1D *hjer13,
	     double etamin, double etamax, const char *name = 0) {
  assert(p2s);
  assert(p2x);

  // Calculate reference JER if not given
  bool newJER13(false);
  if (!hjer13) {
    newJER13 = true;
    int ieta00 = p2s->GetXaxis()->FindBin(0.);
    int ieta13 = p2s->GetXaxis()->FindBin(1.3); // bin edge 1.305
    TProfile *p13s, *p13x;
    p13s = p2s->ProfileY("p13s",ieta00,ieta13,"S");
    p13x = p2x->ProfileY("p13x",ieta00,ieta13,"S");
    hjer13 = p13s->ProjectionX("hjer13_");
    for (int i = 1; i != hjer13->GetNbinsX()+1; ++i) {
      double s = p13s->GetBinError(i);
      double x = p13x->GetBinError(i);
      double jer = sqrt(max(s*s - x*x,0.))/sqrt(2.);
      hjer13->SetBinContent(i, jer);
      hjer13->SetBinError(i, 0.); // do properly later
    }
    delete p13s;
    delete p13x;
  }
  assert(hjer13);

  // Calculate JER2 with MPFX, then subtract reference JER
  int ieta1 = p2s->GetXaxis()->FindBin(etamin);
  int ieta2 = p2s->GetXaxis()->FindBin(etamax);
  TProfile *ps, *px;
  ps = p2s->ProfileY("ps",ieta1,ieta2,"S");
  px = p2x->ProfileY("px",ieta2,ieta2,"S");
  TH1D *hjer = ps->ProjectionX(name!=0 ? name : "hjer_");
  for (int i = 1; i != hjer->GetNbinsX()+1; ++i) {
    double s = ps->GetBinError(i);
    double x = px->GetBinError(i);
    double jer2 = sqrt(max(s*s - x*x,0.));
    double jer13 = hjer13->GetBinContent(i);
    double jer = sqrt(max(jer2*jer2 - jer13*jer13, 0.));
    hjer->SetBinContent(i, jer);
    //hjer->SetBinError(i, 0.); // do properly later

    // Error propagation
    // y = sqrt(a^2 + b^2)
    // => dy = 2*a*0.5/y (oplus) 2*b*0.5/y
    double ns = ps->GetBinEffectiveEntries(i);
    double nx = px->GetBinEffectiveEntries(i);
    double es = (ns>0 ? max(s,0.05)/sqrt(ns) : 0);
    double ex = (nx>0 ? max(x,0.05)/sqrt(nx) : 0);
    double ejer2 = (jer2>0 ? sqrt(pow(s*es,2) + pow(x*ex,2)) / jer2 : 0);
    double ejer13 = hjer13->GetBinError(i);
    double ejer = (jer>0 ? sqrt(pow(jer2*ejer2,2)+pow(jer13*ejer13,2))/jer : 0);
    //ejer = max(ejer,0.005); // Minimum error
    hjer->SetBinError(i, ejer);

    // Cleanup for bad data or MC
    if (jer<0.025) {
      hjer->SetBinContent(i, 0.);
      hjer->SetBinError(i, 0.);
    }
  }
  delete ps;
  delete px;
  if (newJER13) delete hjer13;

  return hjer;
} // getJER

TF1 *fitJER(TH1D *hjer, double ptmin, double ptmax, double eta,
	    const char *name, int color = kBlack, TF1 *f1ref = 0) {

  TF1 *f1 = new TF1(name,"sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])",
		    ptmin,ptmax);
  f1->SetParameters(0,1,0.05);
  // Fix sensible range for each parameter
  if (!f1ref) {
    f1->SetParLimits(0,  -1,  10); // Noise
    f1->SetParLimits(1, 0.5, 1.5); // Stochastic
    f1->SetParLimits(2,0.03,0.25); // Constant
    if (eta>3.139) f1->FixParameter(2,0.065);
  }
  if (f1ref) {
    f1->SetParameters(f1ref->GetParameter(0), f1ref->GetParameter(1),
		      f1ref->GetParameter(2));
    double N = f1ref->GetParameter(0);
    double S = f1ref->GetParameter(1);
    double C = f1ref->GetParameter(2);
    //f1->SetParLimits(0, N-0.2*fabs(N), N+0.2*fabs(N));
    //f1->SetParLimits(0, N-0.05*max(fabs(N),1.), N+0.05*max(fabs(N),1.));
    f1->SetParLimits(0, N, N+0.05*max(fabs(N),1.));
    //f1->SetParLimits(1, S, 1.5*S);
    //f1->SetParLimits(2, C, 5.0*C);
    //f1->SetParLimits(0, N,  10); // Noise
    f1->SetParLimits(1, S, 1.5); // Stochastic
    f1->SetParLimits(2, C,0.25); // Constant

  }
  hjer->Fit(f1,"QRN");
  //f1->SetRange(15,3500);
  f1->SetLineColor(color);
  
  return f1;
} // fitJER
TF1 *ratioJER(TF1 *f1, TF1 *f1m, const char *name, int color) {
  TF1 *f1r = new TF1(name,"sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])/"
		     "sqrt([3]*fabs([3])/(x*x)+[4]*[4]/x+[5]*[5])",
		     15,3500);
  f1r->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		     f1->GetParameter(2),f1m->GetParameter(0),
		     f1m->GetParameter(1),f1m->GetParameter(2));
  f1r->SetLineColor(color);
  
  return f1r;
} // ratioJER

// Draw h2jersf
TH1D *drawH2JERSF(TH2D *h2, double pt, string draw, int marker, int color) {
  int ipt = h2->GetYaxis()->FindBin(pt);
  string id = Form("h%s_%s_%d_%d",h2->GetName(),draw.c_str(),marker,color);
  TH1D *h = h2->ProjectionX(Form("h%s",id.c_str()),ipt,ipt);

  tdrDraw(h,draw.c_str(),marker,color,kSolid,-1,kNone);

  return h;
} // drawH2JERSF

void JERSF() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  string run = "2023D";
  const char *cr = run.c_str();
  string mc = "Summer23MGBPix";
  const char *cm = mc.c_str();
    
  TFile *f = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_data_cmb_%s_JME_v36_Summer23.root",cr),"READ");
  assert(f && !f->IsZombie());

  TFile *fm = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_mc_cmb_%s_v36_Summer23.root",cm),"READ");
  assert(fm && !fm->IsZombie());
  

  curdir->cd();

  TProfile2D *p2s, *p2x, *p2sm, *p2xm;
  p2s = (TProfile2D*)f->Get("Dijet2/p2m0");  assert(p2s);
  p2x = (TProfile2D*)f->Get("Dijet2/p2m0x"); assert(p2x);
  p2sm = (TProfile2D*)fm->Get("Dijet2/p2m0");  assert(p2sm);
  p2xm = (TProfile2D*)fm->Get("Dijet2/p2m0x"); assert(p2xm);

  TH1D *hjer13  = getJER(p2s, p2x, 0, 0,1.3,Form("hjer13_%s",cr));
  TH1D *hjer13m = getJER(p2sm,p2xm,0, 0,1.3,Form("hjer13m_%s",cr));
  TH1D *hsf13 = (TH1D*)hjer13->Clone(Form("hsf13_%s",cr));
  hsf13->Divide(hjer13m);

  TF1 *f13m = fitJER(hjer13m,50,1500,0,Form("f13m_%s",cr),kBlue-9);
  TF1 *f13  = fitJER(hjer13, 50,1500,0,Form("f13_%s",cr),kBlue,f13m);
  TF1 *f13r = ratioJER(f13,f13m,Form("f13r_%s",cr),kBlue);

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);
  
  TCanvas *cx = new TCanvas("cx","cx",7*250,3*250);
  cx->Divide(7,3,0,0);
  TH2D *h2jersf = p2s->ProjectionXY(Form("h2jersf_%s",cr)); h2jersf->Reset();
  
  // Loop over the ieta bins
  vector<TF1*> vf1(p2s->GetNbinsX()+1);
  TH1D *hmin = p2s->ProjectionX(Form("hmin_%s",cr)); hmin->Reset();
  TH1D *hmax = p2s->ProjectionX(Form("hmax_%s",cr)); hmax->Reset();
  for (int ieta = 1; ieta != p2s->GetNbinsX()+1; ++ieta) {
    double etamin = p2s->GetXaxis()->GetBinCenter(ieta);
    double etamax = etamin;
    double eta = p2s->GetXaxis()->GetBinCenter(ieta);
    double eta1 = p2s->GetXaxis()->GetBinLowEdge(ieta);
    double eta2 = p2s->GetXaxis()->GetBinLowEdge(ieta+1);
  // (No indent here for the resf of the loop, maybe function call later)
    
  //double etamin(2.6), etamax(2.6);
  //double etamin(0.0), etamax(1.3);
  TH1D *hjer  = getJER(p2s,p2x,hjer13,etamin,etamax,Form("hjer_%d_%s",ieta,cr));
  TH1D *hjerm = getJER(p2sm,p2xm,hjer13m,etamin,etamax,Form("hjerm_%d_%s",ieta,cr));
  TH1D *hsf = (TH1D*)hjer->Clone(Form("hsf_%d_%s",ieta,cr));
  hsf->Divide(hjerm);

  double ptmin = 50;
  //double ptmax = 1500./cosh(eta1);
  double ptmax = min(1500.,4890./cosh(eta1));
  TF1 *f1m = fitJER(hjerm,ptmin,ptmax,eta,Form("f1m_%d_%s",ieta,cr),kGray+1);
  TF1 *f1  = fitJER(hjer, ptmin,ptmax,eta,Form("f1_%d_%s",ieta,cr),kBlack,f1m);
  TF1 *f1r = ratioJER(f1,f1m,Form("f1r_%d_%s",ieta,cr),kBlack);
  TF1 *f1rb = ratioJER(f1,f1m,Form("f1rb_%d_%s",ieta,cr),kBlack);
  f1rb->SetRange(f1->GetXmin(),f1->GetXmax());
  f1rb->SetLineWidth(2);
  
  TH1D *h = tdrHist(Form("h_%d_%s",ieta,cr),"JER",0,0.45);
  TH1D *hd = tdrHist(Form("hd_%d_%s",ieta,cr),"Data/MC",0.5,2.5);
  lumi_136TeV = cr;//"2023D";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c11_%d_%s",ieta,cr),h,hd,8,33);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  c1->cd(1);
  gPad->SetLogx();

  tex->DrawLatex(0.35,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));

  //tdrDraw(hjer13m,"HISTE",kNone,kBlue-9,kSolid,-1,kNone);
  //tdrDraw(hjer13,"Pz",kOpenCircle,kBlue);
  tdrDraw(hjerm,"HISTE",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hjer,"Pz",kFullCircle,kBlack);

  //f13->Draw("SAME");
  //f13m->Draw("SAME");
  f1->Draw("SAME");
  f1m->Draw("SAME");

  // Print out parameters
  if (true) {
    tex->DrawLatex(0.35,0.75,Form("N=%3.1f#pm%3.1f, S=%4.2f#pm%4.2f",f1->GetParameter(0),f1->GetParError(0),f1->GetParameter(1),f1->GetParError(1)));
    tex->DrawLatex(0.35,0.70,Form("100C=%4.2f#pm%3.2f",100.*f1->GetParameter(2),100.*f1->GetParError(2)));
    tex->DrawLatex(0.35,0.65,Form("n=%3.1f#pm%3.1f, s=%4.2f#pm%4.2f",f1m->GetParameter(0),f1m->GetParError(0),f1m->GetParameter(1),f1m->GetParError(1)));
    tex->DrawLatex(0.35,0.60,Form("100c=%4.2f#pm%3.2f",100.*f1m->GetParameter(2),100.*f1m->GetParError(2)));
    tex->DrawLatex(0.35,0.55,Form("#chi^{2}/NDF=%1.1f/%d, #chi^{2}/ndf=%1.1f/%d",f1->GetChisquare(),f1->GetNDF(),f1m->GetChisquare(),f1m->GetNDF()));
  }
  
  c1->cd(2);
  gPad->SetLogx();

  l->DrawLine(15,1,3500,1);
  
  //tdrDraw(hsf13,"Pz",kOpenCircle,kBlue);
  tdrDraw(hsf,"Pz",kFullCircle,kBlack);

  //f13r->Draw("SAME");
  f1r->Draw("SAME");
  f1rb->Draw("SAME");
  
  
  c1->SaveAs(Form("pdf/JERSF/vsEta/JERSF_eta_%04d_%04d.pdf",
		  int(1000.*eta1),int(1000.*eta2)));

  // Also draw final results into a giant canvas
  cx->cd(ieta);
  TH1D *h6 = tdrHist(Form("h6_%d_%s",ieta,cr),"JER SF",0.5,2.5);
  double ptmaxe = 13600.*0.5/cosh(eta1);
  if      (eta<1.653) h6->GetYaxis()->SetRangeUser(0.8,1.6);
  else if (eta<2.964) h6->GetYaxis()->SetRangeUser(0.5,3.0);
  else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.5,2.5);
  h6->Draw();
  gPad->SetLogx();

  l->SetLineStyle(kDotted);
  l->DrawLine(ptmaxe,0.50,ptmaxe,2.5);
  l->SetLineStyle(kDashed);
  l->DrawLine(15.,1,3500.,1);

  // Draw full range at the back
  tdrDraw(hsf,"Pz",kFullCircle,kBlack);

  f1r->Draw("SAME");
  f1rb->Draw("SAME");
  
  //if (ieta==1 || ieta==nxy) leg1->Draw("SAME");
  tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  
  // Store results with uncertainty to output histogram
  for (int ipt = 1; ipt != h2jersf->GetNbinsY()+1; ++ipt) {
    double pt = h2jersf->GetYaxis()->GetBinCenter(ipt);
    double pt1 = h2jersf->GetYaxis()->GetBinLowEdge(ipt);
    double emax1 = pt1*cosh(eta1);
    double jersf = f1r->Eval(pt);
    double ejersf = 0.; // do properly later

    if (emax1 < 13600.*0.5) {
      h2jersf->SetBinContent(ieta, ipt, jersf);
      h2jersf->SetBinError(ieta, ipt, ejersf); 
    }
  }
  hmin->SetBinContent(ieta, f1r->Eval(10.));
  hmax->SetBinContent(ieta, f1r->Eval(6800./cosh(eta1)));
  
  } // for ieta

  cx->SaveAs(Form("pdf/JERSF/JERSF_AllEta_%s.pdf",cr));
  //cx->SetName(Form("cx_%s",cr));

  // Draw summary of final results in a single plot
  double ymin(0.9), ymax(2.8);
  TH1D *hy = tdrHist("hy","JER SF",ymin,ymax,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *cy = tdrCanvas("cy",hy,8,33,kSquare);

  l->SetLineStyle(kDotted);
  l->DrawLine(1.3,ymin,1.3,2.0);
  l->DrawLine(2.5,ymin,2.5,ymax);
  l->DrawLine(2.964,ymin,2.964,ymax);
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);

  tdrDraw(hmin,"HIST",kNone,kMagenta-9,kSolid,-1,kNone,kMagenta-9);
  tdrDraw(hmax,"HIST",kNone,kGray+1,kSolid,-1,kNone,kGray+1);

  //TLegend *legy0 = tdrLeg(0.18,0.20,0.43,0.20+2*0.045);
  TLegend *legy0 = tdrLeg(0.65,0.63,0.90,0.63+2*0.045);
  legy0->AddEntry(hmax,"E < #sqrt{s}/2","FL");
  legy0->AddEntry(hmin,"p_{T} > 10 GeV","FL");
  
  TH1D *hy15, *hy30, *hy100, *hy300, *hy1000, *hy3000;
  hy100  = drawH2JERSF(h2jersf,100., "HISTE][",kNone,kGreen+2);
  hy100->SetLineWidth(3);

  hy15   = drawH2JERSF(h2jersf,15.,  "HISTE][",kNone,kMagenta+2);
  hy30   = drawH2JERSF(h2jersf,30.,  "HISTE][",kNone,kBlue);
  hy300  = drawH2JERSF(h2jersf,300., "HISTE][",kNone,kOrange+2);
  hy1000 = drawH2JERSF(h2jersf,1000.,"HISTE][",kNone,kRed);
  hy3000 = drawH2JERSF(h2jersf,3000.,"HISTE][",kNone,kBlack);

  TLegend *legy = tdrLeg(0.18,0.90-6*0.045,0.43,0.90);
  legy->AddEntry(hy15,  "p_{T} = 15 GeV",  "PLE");
  legy->AddEntry(hy30,  "p_{T} = 30 GeV",  "PLE");
  legy->AddEntry(hy100, "p_{T} = 100 GeV", "PLE");
  legy->AddEntry(hy300, "p_{T} = 300 GeV", "PLE");
  legy->AddEntry(hy1000,"p_{T} = 1000 GeV","PLE");
  legy->AddEntry(hy3000,"p_{T} = 3000 GeV","PLE");
  
  cy->SaveAs(Form("pdf/JERSF/JERSF_Summary_%s.pdf",cr));
  cy->SetName(Form("cy_%s",cr));
  
} // JERSF
