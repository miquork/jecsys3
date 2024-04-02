// Purpose: derive JER SF with MPFX method
//          combine multiple channels, start from TProfile2D *p2m0
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFitResultPtr.h"

#include "tdrstyle_mod22.C"
#include "tools.C"


bool plotEtaBins = true; // lots of plots


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

// Calculate JER with MPFX for Z+jet. JER13 not needed
TH1D *getJERZ(TProfile2D *p2s, TProfile2D *p2x, 
	      double etamin, double etamax, const char *name = 0) {
  assert(p2s);
  assert(p2x);

  // Calculate JER with MPFX
  int ieta1 = p2s->GetXaxis()->FindBin(etamin);
  int ieta2 = p2s->GetXaxis()->FindBin(etamax);
  TProfile *ps, *px;
  ps = p2s->ProfileY("ps",ieta1,ieta2,"S");
  px = p2x->ProfileY("px",ieta2,ieta2,"S");
  TH1D *hjer = ps->ProjectionX(name!=0 ? name : "hjer_");
  for (int i = 1; i != hjer->GetNbinsX()+1; ++i) {
    double s = ps->GetBinError(i);
    double x = px->GetBinError(i);
    double jer = sqrt(max(s*s - x*x,0.));
    hjer->SetBinContent(i, jer);
    //hjer->SetBinError(i, 0.); // do properly later

    // Error propagation
    // y = sqrt(a^2 + b^2)
    // => dy = 2*a*0.5/y (oplus) 2*b*0.5/y
    double ns = ps->GetBinEffectiveEntries(i);
    double nx = px->GetBinEffectiveEntries(i);
    double es = (ns>0 ? max(s,0.05)/sqrt(ns) : 0);
    double ex = (nx>0 ? max(x,0.05)/sqrt(nx) : 0);
    double ejer = (jer>0 ? sqrt(pow(s*es,2) + pow(x*ex,2)) / jer : 0);
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

  return hjer;
} // getJERZ

TF1 *fitJER(TH1D *hjer, double ptmin, double ptmax, double eta,
	    const char *name, int color = kBlack, TF1 *f1ref = 0,
	    TFitResultPtr *f1ptr = 0) {

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
  if (f1ptr) (*f1ptr) = hjer->Fit(f1,"QRNS");
  else                  hjer->Fit(f1,"QRN");
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

  // Set graphical styles
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf/JERSF/vsEta");
  //gROOT->ProcessLine(".! mkdir pdf/JERSF/vsPt");
  
  // Set output directory;
  TFile *fout = new TFile("rootfiles/JERSF.root","RECREATE");
  fout->mkdir("Fits");
  fout->mkdir("Dijet");
  fout->mkdir("Extras");
  fout->mkdir("Raw");

  string vrun[] = {"2023Cv123","2023Cv4","2023D"};//,"2023Cv4D"};
  //string vrun[] = {"2023D"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  //string vmc[] = {"Summer23","Summer23","Summer23BPIX"};//,"Summer23BPIX"};
  string vmc[] = {"Summer23MG","Summer23MG","Summer23MGBPix"};
  //string vmc[] = {"Summer23MG_Cv123","Summer23MG_Cv4","Summer23MGBPix_D"};
  //string vmc[] = {"Summer23BPIX"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);
  assert(nmc==nrun);
  for (int irun = 0; irun != nrun; ++irun) {
    string run = vrun[irun];
    const char *cr = run.c_str();
    string mc = vmc[irun];
    const char *cm = mc.c_str();
  // (No indent here for the resf of the loop, maybe function call later)
    
  //string run = "2023D";
  //const char *cr = run.c_str();
  //string mc = "Summer23MGBPix";
  //const char *cm = mc.c_str();
    
  //TFile *f = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_data_cmb_%s_JME_v36_Summer23.root",cr),"READ");
    TFile *f = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_data_cmb_%s_JME_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cr),"READ");
    //TFile *f = new TFile(Form("rootfiles/Summer23_L2L3Res/jmenano_data_cmb_%s_JME_v39_L2Rel_L2L3Res_v2_SF.root",cr),"READ");
  assert(f && !f->IsZombie());

  //TFile *fm = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_mc_cmb_%s_v36_Summer23.root",cm),"READ");
  //TFile *fm = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_mc_cmb_%s_v36_Summer23.root","Summer23MGBPix"),"READ"); // Summer23 patch
  TFile *fm = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_%s_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cm),"READ");
  //TFile *fm = new TFile(Form("rootfiles/Summer23_L2L3Res/jmenano_mc_cmb_%s_RwPU_v39_SmearJets_L2Res_v1_SF.root",cm),"READ");
  assert(fm && !fm->IsZombie());

  //TFile *fz = new TFile(Form("rootfiles/Summer23_noL2L3Res/jme_bplusZ_%s_Zmm_sync_v69.root",cr),"READ");
  TFile *fz = new TFile(Form("rootfiles/Summer23_L2ResOnly/jme_bplusZ_%s_Zmm_sync_v70.root",cr),"READ");
  assert(fz && !fz->IsZombie());

  //TFile *fg = new TFile(Form("rootfiles/Summer23_noL2L3Res/GamHistosFill_data_%s_w2.root",cr),"READ"); // Summer23 (w3->w2)
  TFile *fg = new TFile(Form("rootfiles/Summer23_L2ResOnly/GamHistosFill_data_%s_w6.root",cr),"READ");
  //TFile *fg = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s_w4.root",cr),"READ"); // Summer23 with L2Res
  assert(fg && !fg->IsZombie());
  //
  //TFile *fgm = new TFile(run=="2023D" ? "rootfiles/Summer23_noL2L3Res/GamHistosFill_mc_2023P8-BPix_w2.root" : "rootfiles/Summer23_noL2L3Res/GamHistosFill_mc_2023P8_w2.root","READ"); // Summer23 (w3->w2)
  TFile *fgm = new TFile(run=="2023D" ? "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8-BPix_w6.root" : "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8_w6.root","READ"); // Summer23 (w3->w2)
  assert(fgm && !fgm->IsZombie());

  curdir->cd();

  TProfile2D *p2s, *p2x, *p2sm, *p2xm;
  p2s = (TProfile2D*)f->Get("Dijet2/p2m0");  assert(p2s);
  p2x = (TProfile2D*)f->Get("Dijet2/p2m0x"); assert(p2x);
  p2sm = (TProfile2D*)fm->Get("Dijet2/p2m0");  assert(p2sm);
  p2xm = (TProfile2D*)fm->Get("Dijet2/p2m0x"); assert(p2xm);

  TProfile2D *p2zs, *p2zx, *p2zsm, *p2zxm;
  p2zs = (TProfile2D*)fz->Get("data/l2res/p2m0");  assert(p2zs);
  p2zx = (TProfile2D*)fz->Get("data/l2res/p2m0x"); assert(p2zx);
  p2zsm = (TProfile2D*)fz->Get("mc/l2res/p2m0");  assert(p2zsm);
  p2zxm = (TProfile2D*)fz->Get("mc/l2res/p2m0x"); assert(p2zxm);

  TProfile2D *p2gs, *p2gx, *p2gsm, *p2gxm;
  p2gs = (TProfile2D*)fg->Get("Gamjet2/p2m0");  assert(p2gs);
  p2gx = (TProfile2D*)fg->Get("Gamjet2/p2m0x"); assert(p2gx);
  p2gsm = (TProfile2D*)fgm->Get("Gamjet2/p2m0");  assert(p2gsm);
  p2gxm = (TProfile2D*)fgm->Get("Gamjet2/p2m0x"); assert(p2gxm);
  
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
  TH2D *h2jersf0 = p2s->ProjectionXY(Form("h2jersf0_%s",cr)); h2jersf0->Reset();
  //
  TH2D *h2jerdt = p2s->ProjectionXY(Form("h2jerdt_%s",cr)); h2jerdt->Reset();
  TH2D *h2jermc = p2s->ProjectionXY(Form("h2jermc_%s",cr)); h2jermc->Reset();
  //
  TH2D *h2zjer = p2zs->ProjectionXY(Form("h2zjer_%s",cr)); h2zjer->Reset();
  TH2D *h2zjerdt =p2zs->ProjectionXY(Form("h2zjerdt_%s",cr)); h2zjerdt->Reset();
  TH2D *h2zjermc =p2zs->ProjectionXY(Form("h2zjermc_%s",cr)); h2zjermc->Reset();
  //
  TH2D *h2gjer = p2gs->ProjectionXY(Form("h2gjer_%s",cr)); h2gjer->Reset();
  TH2D *h2gjerdt =p2gs->ProjectionXY(Form("h2gjerdt_%s",cr)); h2gjerdt->Reset();
  TH2D *h2gjermc =p2gs->ProjectionXY(Form("h2gjermc_%s",cr)); h2gjermc->Reset();
  
  // Loop over the ieta bins
  vector<TF1*> vf1(p2s->GetNbinsX()+1);
  TH1D *hchi2 = p2s->ProjectionX(Form("hchi2_%s",cr)); hchi2->Reset();
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

  TH1D *hjerz  = getJERZ(p2zs,p2zx,etamin,etamax,Form("hjerz_%d_%s",ieta,cr));
  TH1D *hjerzm = getJERZ(p2zsm,p2zxm,etamin,etamax,Form("hjerzm_%d_%s",ieta,cr));
  TH1D *hsfz = (TH1D*)hjerz->Clone(Form("hsfz_%d_%s",ieta,cr));
  hsfz->Divide(hjerzm);

  TH1D *hjerg  = getJERZ(p2gs,p2gx,etamin,etamax,Form("hjerg_%d_%s",ieta,cr));
  TH1D *hjergm = getJERZ(p2gsm,p2gxm,etamin,etamax,Form("hjergm_%d_%s",ieta,cr));
  TH1D *hsfg = (TH1D*)hjerg->Clone(Form("hsfg_%d_%s",ieta,cr));
  hsfg->Divide(hjergm);
  
  double ptmin = 50;
  //double ptmax = 1500./cosh(eta1);
  double ptmax = min(1500.,4890./cosh(eta1));
  // Clean up some high chi2 bins
  if (eta>2.964 && eta<3.139) ptmax = 302;
  if (eta>3.139 && eta<3.314) ptmax = 236;

  TFitResultPtr f1mptr, f1ptr;
  TF1 *f1m = fitJER(hjerm,ptmin,ptmax,eta,Form("f1m_%d_%s",ieta,cr),kGray+1,
		    0, &f1mptr);
  TF1 *f1  = fitJER(hjer, ptmin,ptmax,eta,Form("f1_%d_%s",ieta,cr),kBlack,
		    f1m, &f1ptr);
  TF1 *f1r = ratioJER(f1,f1m,Form("f1r_%d_%s",ieta,cr),kBlack);
  TF1 *f1rb = ratioJER(f1,f1m,Form("f1rb_%d_%s",ieta,cr),kBlack);
  f1rb->SetRange(f1->GetXmin(),f1->GetXmax());
  f1rb->SetLineWidth(3);//2);
  vf1[ieta-1] = f1r; // save for printout
  
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

  tdrDraw(hjerzm,"HISTE",kNone,kRed,kSolid,-1,kNone);
  tdrDraw(hjerz,"Pz",kOpenCircle,kRed);

  tdrDraw(hjergm,"HISTE",kNone,kBlue,kSolid,-1,kNone);
  tdrDraw(hjerg,"Pz",kOpenDiamond,kBlue);

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

  tdrDraw(hsfz,"Pz",kOpenCircle,kRed);
  tdrDraw(hsfg,"Pz",kOpenDiamond,kBlue);
  
  //tdrDraw(hsf13,"Pz",kOpenCircle,kBlue);
  tdrDraw(hsf,"Pz",kFullCircle,kBlack);

  //f13r->Draw("SAME");
  f1r->Draw("SAME");
  f1rb->Draw("SAME");
  
  if (plotEtaBins)
    c1->SaveAs(Form("pdf/JERSF/vsEta/JERSF_eta_%04d_%04d_%s.pdf",
		    int(1000.*eta1),int(1000.*eta2),cr));

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
    //double ejersf = 0.; // do properly later
    double d = f1->Eval(pt);
    double ed = tools::getFitErr(f1, f1ptr, pt);
    double m = f1m->Eval(pt);
    double em = tools::getFitErr(f1m, f1mptr, pt);
    // r = d/m => dr = dd/m (+) d*dm/m^2 => dr = d/m*(dd/d (+) dm/m)
    double ejersf = jersf * sqrt(pow(ed/d,2) + pow(em/m,2));
    
    if (emax1 < 13600.*0.5) {
      h2jersf->SetBinContent(ieta, ipt, jersf);
      h2jersf->SetBinError(ieta, ipt, ejersf);
    }
    h2jersf0->SetBinContent(ieta, ipt, hsf->GetBinContent(ipt));
    h2jersf0->SetBinError(ieta, ipt, hsf->GetBinError(ipt));
    //
    h2jerdt->SetBinContent(ieta, ipt, hjer->GetBinContent(ipt));
    h2jerdt->SetBinError(ieta, ipt, hjer->GetBinError(ipt));
    h2jermc->SetBinContent(ieta, ipt, hjerm->GetBinContent(ipt));
    h2jermc->SetBinError(ieta, ipt, hjerm->GetBinError(ipt)); 

  }
  // Store results with uncertainty to output histogram for Z/gamma+jet
  for (int ipt = 1; ipt != h2zjer->GetNbinsY()+1; ++ipt) {
    h2zjer->SetBinContent(ieta, ipt, hsfz->GetBinContent(ipt));
    h2zjer->SetBinError(ieta, ipt, hsfz->GetBinError(ipt));
    h2zjerdt->SetBinContent(ieta, ipt, hjerz->GetBinContent(ipt));
    h2zjerdt->SetBinError(ieta, ipt, hjerz->GetBinError(ipt));
    h2zjermc->SetBinContent(ieta, ipt, hjerzm->GetBinContent(ipt));
    h2zjermc->SetBinError(ieta, ipt, hjerzm->GetBinError(ipt)); 
  }
  for (int ipt = 1; ipt != h2gjer->GetNbinsY()+1; ++ipt) {
    h2gjer->SetBinContent(ieta, ipt, hsfg->GetBinContent(ipt));
    h2gjer->SetBinError(ieta, ipt, hsfg->GetBinError(ipt));
    h2gjerdt->SetBinContent(ieta, ipt, hjerg->GetBinContent(ipt));
    h2gjerdt->SetBinError(ieta, ipt, hjerg->GetBinError(ipt));
    h2gjermc->SetBinContent(ieta, ipt, hjergm->GetBinContent(ipt));
    h2gjermc->SetBinError(ieta, ipt, hjergm->GetBinError(ipt)); 
  }

  
  hmin->SetBinContent(ieta, f1r->Eval(10.));
  hmax->SetBinContent(ieta, f1r->Eval(6800./cosh(eta1)));

  hchi2->SetBinContent(ieta, f1->GetChisquare()/max(1,f1->GetNDF()));
  hchi2->SetBinError(ieta, 1./sqrt(max(1,f1->GetNDF())));

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

  
  // Draw chi2 to have quick control of fit quality
  TH1D *hz = tdrHist("hz","Fit #chi^{2} / NDF",0,10,"|#eta|",0,5.2);
  TCanvas *cz = tdrCanvas("cz",hz,8,33,kSquare);

  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);
  
  tdrDraw(hchi2,"HPz",kFullCircle,kBlack,kSolid,1,kNone);
  
  cz->SaveAs(Form("pdf/JERSF/JERSF_Chi2_%s.pdf",cr));


  // Produce output text file (FactorizedJetCorrecter Style
  ofstream txt(Form("pdf/JERSF/Summer23_%s_JRV1_MC_SF_AK4PFPuppi.txt",cr));
  txt << "{1 JetEta 1 JetPt "
      << "sqrt([0]*fabs([0])/(x*x)+[1]*[1]/x+[2]*[2])/"
      << "sqrt([3]*fabs([3])/(x*x)+[4]*[4]/x+[5]*[5])"
      << " Correction L2Relative}" << endl; // awkward name, but runs
  // Run2 reference header:
  //txt << "{1 JetEta 2 JetPt Rho  "
  //  << "sqrt(([0]*[0]*[3]*[3]*y/[6])/(x*x)+[1]*[1]*[4]*[4]/x+[2]*[2]*[5]*[5])"
  //  << "/sqrt(([0]*[0]*y/[6])/(x*x)+[1]*[1]/x+[2]*[2])"
  //  << " Correction L2Relative}" << endl; // test: runs

  int neta = p2s->GetNbinsX();
  for (int i = neta-1; i != -1; --i) {
    TF1 *f1 = vf1[i];
    double abseta1 = p2s->GetXaxis()->GetBinLowEdge((i+1));
    double abseta2 = p2s->GetXaxis()->GetBinLowEdge((i+1)+1);
    txt << Form("%5.3f %5.3f %3d %4d %4d %8.3f %6.4f %7.5f %8.3f %6.4f %7.5f\n",
		-abseta2, -abseta1,
		2+6, 15, int(6800. / cosh(abseta1)),
		f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),
		f1->GetParameter(3),f1->GetParameter(4),f1->GetParameter(5));
    // Run2 reference header: 
    //JER &jer = vjer[i];
    //txt << Form("%5.3f %6.3f %2d %2d %4d %2d %2d %5.2f %5.3f %6.4f "
    //		" %5.3f %5.3f %5.3f %5.2f\n",
    //		-jer.eta-jer.deta, -jer.eta+jer.deta, 11, 8, 4000, 0, 70,
    //		  jer.n, jer.s, jer.c, jer.kn, jer.ks, jer.kc, rho);
  } // for i in -neta
  for (int i = 0; i != neta; ++i) {
    TF1 *f1 = vf1[i];
    double abseta1 = p2s->GetXaxis()->GetBinLowEdge((i+1));
    double abseta2 = p2s->GetXaxis()->GetBinLowEdge((i+1)+1);
    txt << Form("%5.3f %5.3f %3d %4d %4d %8.3f %6.4f %7.5f %8.3f %6.4f %7.5f\n",
		abseta1, abseta2,
		2+6, 15, int(6800. / cosh(abseta1)),
		f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),
		f1->GetParameter(3),f1->GetParameter(4),f1->GetParameter(5));
  } // for i in +neta

  // Store results to output file
  fout->cd();

  //if (!fout->FindObject("Fits")) fout->mkdir("Fits");
  fout->cd("Fits");
  h2jersf->Write(Form("h2jersf_%s_%s",cr,cm),TObject::kOverwrite);

  //if (!fout->FindObject("Dijet")) fout->mkdir("Dijet");
  fout->cd("Dijet");
  h2jersf0->Write(Form("h2jersfRaw_%s_%s",cr,cm),TObject::kOverwrite);
  h2jerdt->Write(Form("h2jerdtRaw_%s_%s",cr,cm),TObject::kOverwrite);
  h2jermc->Write(Form("h2jermcRaw_%s_%s",cr,cm),TObject::kOverwrite);

  //if (!fout->FindObject("Extras")) fout->mkdir("Extras");
  fout->cd("Extras");
  h2zjer->Write(Form("h2jersfZjet_%s_%s",cr,cm),TObject::kOverwrite);
  h2zjerdt->Write(Form("h2jerdtZjet_%s_%s",cr,cm),TObject::kOverwrite);
  h2zjermc->Write(Form("h2jermcZjet_%s_%s",cr,cm),TObject::kOverwrite);
  //
  h2gjer->Write(Form("h2jersfGamjet_%s_%s",cr,cm),TObject::kOverwrite);
  h2gjerdt->Write(Form("h2jerdtGamjet_%s_%s",cr,cm),TObject::kOverwrite);
  h2gjermc->Write(Form("h2jermcGamjet_%s_%s",cr,cm),TObject::kOverwrite);

  //if (!fout->FindObject("Raw")) fout->mkdir("Raw");
  fout->cd("Raw");
  p2s->Write(Form("p2m0dt_%s_%s",cr,cm),TObject::kOverwrite);
  p2sm->Write(Form("p2m0mc_%s_%s",cr,cm),TObject::kOverwrite);
  p2x->Write(Form("p2m0xdt_%s_%s",cr,cm),TObject::kOverwrite);
  p2xm->Write(Form("p2m0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  //
  p2zs->Write(Form("p2zm0dt_%s_%s",cr,cm),TObject::kOverwrite);
  p2zsm->Write(Form("pzgm0mc_%s_%s",cr,cm),TObject::kOverwrite);
  p2zx->Write(Form("p2zm0xdt_%s_%s",cr,cm),TObject::kOverwrite);
  p2zxm->Write(Form("p2zm0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  //
  p2gs->Write(Form("p2gm0dt_%s_%s",cr,cm),TObject::kOverwrite);
  p2gsm->Write(Form("p2gm0mc_%s_%s",cr,cm),TObject::kOverwrite);
  p2gx->Write(Form("p2gm0xdt_%s_%s",cr,cm),TObject::kOverwrite);
  p2gxm->Write(Form("p2gm0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  curdir->cd();

  } // for irunx

  fout->Close();
} // JERSF
