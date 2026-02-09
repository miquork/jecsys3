// Purpose: Draw comparison of pT,ave to pT,Z and pT,jet1 from jecdata.root
//          Compare MPF, MPF1, MPFn, MPFu, and also HDM
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include "../tdrstyle_mod22.C"
#include "../tools.C"

#include <map>
#include <vector>
#include <string>
using namespace std;

const bool debug = true;//false;

const bool addZ = true;//false; // add Z+jet
const bool addGam = true; // add photon+jet


void drawZJetPtAveV2_jer();
void drawZJetPtAveV2_fsr();
void drawZJetPtAveV2_bin();
void drawZJetPtAveV2_bins(bool _addZ = addZ, bool _addGam = addGam);

void drawZJetPtAveV2() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/drawZJetPtAveV2");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/drawZJetPtAveV2");
  
  //drawZJetPtAveV2_jer();
  //drawZJetPtAveV2_fsr();
  //drawZJetPtAveV2_bin();
  drawZJetPtAveV2_bins(true,false);
  drawZJetPtAveV2_bins(false,true);
  drawZJetPtAveV2_bins();
}

void cleanGraph(TGraph *g, double xmin, double xmax) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetX()[i]<xmin || g->GetX()[i]>xmax) {
      g->RemovePoint(i);
    }
  }
} // cleanGraph

TGraphErrors* addGraphs(TGraphErrors *g1, TGraphErrors *g2,
			double c1=1, double c2=1) {
  TGraphErrors *g = new TGraphErrors(g1->GetN());
  for (int i = 0; i != min(g->GetN(),g2->GetN()); ++i) {
    g->SetPoint(i, 0.5*(g1->GetX()[i]+g2->GetX()[i]),
		(c1*g1->GetY()[i]+c2*g2->GetY()[i]));
    g->SetPointError(i, sqrt(pow(0.5*fabs(g1->GetX()[i]-g2->GetX()[i]),2) +
			     pow(g1->GetEX()[i],2)+pow(g2->GetEX()[i],2)),
		     sqrt(pow(c1*g1->GetEY()[i],2)+pow(c2*g2->GetEY()[i],2)));
  }
  return g;
} // addGraphs

// SF vs noSF
void drawZJetPtAveV2_jer() {

  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  string tag = "";
  
  //TGraphErrors *gm(0), *gd(0), *gr(0);
  TH1D *hu(0), *hd(0);
  lumi_13TeV = tag;
  TCanvas *c1(0);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);

  hu = tdrHist("hu","HDM response",0.985,1.045);
  hd = tdrHist("hd","JERSF-NoSF (%)",-1.8,+0.8);
  c1 = tdrDiCanvas("c1",hu,hd,4,11);
  double ptmin = hu->GetXaxis()->GetXmin();
  double ptmax = hu->GetXaxis()->GetXmax();
    
  TLegend *leg = tdrLeg(0.45,0.92-2*0.04,0.75,0.92);
  
  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(ptmin,1,ptmax,1);
  l->DrawLine(ptmin,0,ptmax,0);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(ptmin,0,ptmax,0);

  //TFile *fs = new TFile("rootfiles/ZJetPtAveV2_zpt_20230824_v57/jecdataRun2Test.root","READ");
  TFile *fs = new TFile("rootfiles/ZJetPtAveV2_ptave_20230824_v57/jecdataRun2Test.root","READ");
  assert(fs && !fs->IsZombie());
  //TFile *fn = new TFile("rootfiles/ZJetPtAveV2_zpt_20230825_v53/jecdataRun2Test.root","READ");
  TFile *fn = new TFile("rootfiles/ZJetPtAveV2_ptave_20230825_v53/jecdataRun2Test.root","READ");
  assert(fn && !fn->IsZombie());
  curdir->cd();

  string s = "ratio/eta00-13/hdm_mpfchs1_zjet";
  TH1D *hs = (TH1D*)fs->Get(s.c_str()); assert(hs);
  TH1D *hn = (TH1D*)fn->Get(s.c_str()); assert(hn);

  c1->cd(1);
  tdrDraw(hs,"Pz",kFullCircle,kRed);
  tdrDraw(hn,"Pz",kOpenCircle,kBlue);
  leg->AddEntry(hs,"JER SF","PE");
  leg->AddEntry(hn,"No SF","PE");

  c1->cd(2);

  //TH1D *hr = (TH1D*)hs->Clone("hr");
  //hr->Divide(hn);
  TGraphErrors *gr = tools::diffGraphs(new TGraphErrors(hs),
				       new TGraphErrors(hn),
				       100,100);

  //tdrDraw(hr,"Pz",kFullCircle,kRed);
  tdrDraw(gr,"Pz",kFullCircle,kRed);

  c1->SaveAs("pdf/drawZJetPtAveV2/drawZJetPtAveV2_jer.pdf");

} // void drawZJetPtAve_jer

// FSR = 1 - (DB / HDM)
void drawZJetPtAveV2_fsr() {

  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  string tag = "";
  
  TH1D *hu(0), *hd(0);
  lumi_13TeV = tag;
  TCanvas *c1(0);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);

  hu = tdrHist("hu","FSR",-0.15,0.30);
  hd = tdrHist("hd","Data/MC-1 (%)",-3,+3);
  c1 = tdrDiCanvas("c1",hu,hd,4,11);
  double ptmin = hu->GetXaxis()->GetXmin();
  double ptmax = hu->GetXaxis()->GetXmax();
    
  TLegend *leg = tdrLeg(0.45,0.92-3*0.04,0.75,0.92);
  
  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(ptmin,1,ptmax,1);
  l->DrawLine(ptmin,0,ptmax,0);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(ptmin,0,ptmax,0);

  TFile *fs = new TFile("rootfiles/ZJetPtAveV2_zpt_20230824_v57/jecdataRun2Test.root","READ");
  assert(fs && !fs->IsZombie());
  TFile *fa = new TFile("rootfiles/ZJetPtAveV2_ptave_20230824_v57/jecdataRun2Test.root","READ");
  assert(fa && !fa->IsZombie());
  TFile *fj = new TFile("rootfiles/ZJetPtAveV2_jetpt_20230824_v57/jecdataRun2Test.root","READ");
  assert(fj && !fj->IsZombie());
  curdir->cd();

  string shd = "data/eta00-13/hdm_mpfchs1_zjet";
  TH1D *hhd = (TH1D*)fs->Get(shd.c_str()); assert(hhd);
  string shm = "mc/eta00-13/hdm_mpfchs1_zjet";
  TH1D *hhm = (TH1D*)fs->Get(shm.c_str()); assert(hhm);

  TH1D *hhad = (TH1D*)fa->Get(shd.c_str()); assert(hhad);
  TH1D *hham = (TH1D*)fa->Get(shm.c_str()); assert(hham);
  TH1D *hhjd = (TH1D*)fj->Get(shd.c_str()); assert(hhjd);
  TH1D *hhjm = (TH1D*)fj->Get(shm.c_str()); assert(hhjm);

  TGraphErrors *gdd(0), *gdm(0), *gdad(0), *gdam(0), *gdjd(0), *gdjm(0);
  string sdd = "data/eta00-13/mpf1_zjet_a100";
  gdd = (TGraphErrors*)fs->Get(sdd.c_str()); assert(gdd);
  string sdm = "mc/eta00-13/mpf1_zjet_a100";
  gdm = (TGraphErrors*)fs->Get(sdm.c_str()); assert(gdm);

  gdad = (TGraphErrors*)fa->Get(sdd.c_str()); assert(gdad);
  gdam = (TGraphErrors*)fa->Get(sdm.c_str()); assert(gdam);
  gdjd = (TGraphErrors*)fj->Get(sdd.c_str()); assert(gdjd);
  gdjm = (TGraphErrors*)fj->Get(sdm.c_str()); assert(gdjm);

  TGraphErrors *gfd = (TGraphErrors*)gdd->Clone("gfd");
  for (int i = 0; i != gfd->GetN(); ++i) {
    double pt = gdd->GetX()[i];
    double j = hhd->GetXaxis()->FindBin(pt);
    double r = hhd->GetBinContent(j);
    gfd->SetPoint(i, pt, r!=0 ? 1. - (gdd->GetY()[i] / r) : 0);
  }
  //
  TGraphErrors *gfm = (TGraphErrors*)gdm->Clone("gfm");
  for (int i = 0; i != gfm->GetN(); ++i) {
    double pt = gdm->GetX()[i];
    double j = hhm->GetXaxis()->FindBin(pt);
    double r = hhm->GetBinContent(j);
    gfm->SetPoint(i, pt, r!=0 ? 1. - (gdm->GetY()[i] / r) : 0);
  }
  //
  TGraphErrors *gfad = (TGraphErrors*)gdad->Clone("gfad");
  for (int i = 0; i != gfad->GetN(); ++i) {
    double pt = gdad->GetX()[i];
    double j = hhad->GetXaxis()->FindBin(pt);
    double r = hhad->GetBinContent(j);
    gfad->SetPoint(i, pt, r!=0 ? 1. - (gdad->GetY()[i] / r) : 0);
  }
  //
  TGraphErrors *gfam = (TGraphErrors*)gdam->Clone("gfam");
  for (int i = 0; i != gfam->GetN(); ++i) {
    double pt = gdam->GetX()[i];
    double j = hham->GetXaxis()->FindBin(pt);
    double r = hham->GetBinContent(j);
    gfam->SetPoint(i, pt, r!=0 ? 1. - (gdam->GetY()[i] / r) : 0);
  }
  //
  TGraphErrors *gfjd = (TGraphErrors*)gdjd->Clone("gfjd");
  for (int i = 0; i != gfjd->GetN(); ++i) {
    double pt = gdjd->GetX()[i];
    double j = hhjd->GetXaxis()->FindBin(pt);
    double r = hhjd->GetBinContent(j);
    gfjd->SetPoint(i, pt, r!=0 ? 1. - (gdjd->GetY()[i] / r) : 0);
  }
  //
  TGraphErrors *gfjm = (TGraphErrors*)gdjm->Clone("gfjm");
  for (int i = 0; i != gfjm->GetN(); ++i) {
    double pt = gdjm->GetX()[i];
    double j = hhjm->GetXaxis()->FindBin(pt);
    double r = hhjm->GetBinContent(j);
    gfjm->SetPoint(i, pt, r!=0 ? 1. - (gdjm->GetY()[i] / r) : 0);
  }


  c1->cd(1);
  tdrDraw(gfad,"Pz",kFullSquare,kGreen+2);
  tdrDraw(gfam,"Pz",kOpenSquare,kGreen+2);
  tdrDraw(gfd,"Pz",kFullCircle,kRed);
  tdrDraw(gfm,"Pz",kOpenCircle,kRed);
  tdrDraw(gfjd,"Pz",kFullDiamond,kBlue);
  tdrDraw(gfjm,"Pz",kOpenDiamond,kBlue);
  leg->AddEntry(gfd,"Data ZPt","PE");
  //leg->AddEntry(gfm,"MC","PE");
  leg->AddEntry(gfad,"Data PtAve","PE");
  //leg->AddEntry(gfam,"MC PtAve","PE");
  leg->AddEntry(gfjd,"Data JetPt","PE");
  //leg->AddEntry(gfjm,"MC JetPt","PE");

  c1->cd(2);

  TGraphErrors *gfr = tools::diffGraphs(gfd,gfm,100,100);
  TGraphErrors *gfar = tools::diffGraphs(gfad,gfam,100,100);
  TGraphErrors *gfjr = tools::diffGraphs(gfjd,gfjm,100,100);
  tdrDraw(gfar,"Pz",kFullSquare,kGreen+2);
  tdrDraw(gfr,"Pz",kFullCircle,kRed);
  tdrDraw(gfjr,"Pz",kFullDiamond,kBlue);

  c1->SaveAs("pdf/drawZJetPtAveV2/drawZJetPtAveV2_fsr.pdf");

} // void drawZJetPtAve_fsr

// pTave vs zpt
void drawZJetPtAveV2_bin() {

  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  string tag = "";
  
  TH1D *hu(0), *hd(0);
  lumi_13TeV = tag;
  TCanvas *c1(0);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);

  hu = tdrHist("hu","HDM response",0.985,1.045);
  hd = tdrHist("hd","pTave-Zpt (%)",-2.2,2.2);
  c1 = tdrDiCanvas("c1",hu,hd,4,11);
  double ptmin = hu->GetXaxis()->GetXmin();
  double ptmax = hu->GetXaxis()->GetXmax();
    
  TLegend *leg = tdrLeg(0.45,0.92-2*0.04,0.75,0.92);
  
  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(ptmin,1,ptmax,1);
  l->DrawLine(ptmin,0,ptmax,0);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(ptmin,0,ptmax,0);

  TFile *fs = new TFile("rootfiles/ZJetPtAveV2_ptave_20230824_v57/jecdataRun2Test.root","READ");
  //TFile *fs = new TFile("rootfiles/ZJetPtAveV2_ptave_20230825_v53/jecdataRun2Test.root","READ");
  assert(fs && !fs->IsZombie());
  TFile *fn = new TFile("rootfiles/ZJetPtAveV2_zpt_20230824_v57/jecdataRun2Test.root","READ");
  //TFile *fn = new TFile("rootfiles/ZJetPtAveV2_zpt_20230825_v53/jecdataRun2Test.root","READ");
  assert(fn && !fn->IsZombie());
  curdir->cd();

  string s = "ratio/eta00-13/hdm_mpfchs1_zjet";
  TH1D *hs = (TH1D*)fs->Get(s.c_str()); assert(hs);
  TH1D *hn = (TH1D*)fn->Get(s.c_str()); assert(hn);

  c1->cd(1);
  tdrDraw(hs,"Pz",kFullCircle,kRed);
  tdrDraw(hn,"Pz",kOpenCircle,kBlue);
  leg->AddEntry(hs,"PtAve","PE");
  leg->AddEntry(hn,"Zpt","PE");

  c1->cd(2);

  TGraphErrors *gr = tools::diffGraphs(new TGraphErrors(hs),
				       new TGraphErrors(hn),
				       100,100);

  tdrDraw(gr,"Pz",kFullCircle,kRed);

  c1->SaveAs("pdf/drawZJetPtAveV2/drawZJetPtAveV2_bin.pdf");

} // void drawZJetPtAve_bin

void drawZJetPtAveV2_bins(bool addGam, bool addZ) {

  double Ru_d(0.92*0.80), Ru_m(0.92);
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  //string tag = "2025CDEFG, 98.1 fb^{-1}";
  #include "../Config.C"
  string tag = Form("2025CDEFG, %s",mlum["2025"].c_str());
  
  vector<string> vobs;
  vobs.push_back("hdm_mpfchs1");
  vobs.push_back("mpfchs1");
  vobs.push_back("mpf1");
  vobs.push_back("mpfn");
  vobs.push_back("mpfu");
  vobs.push_back("mpfnu");
  vobs.push_back("counts");

  vector<string> vbin;
  // Run 2 convention
  //vbin.push_back("zpt");
  //vbin.push_back("ptave");
  //vbin.push_back("jetpt");
  // Run 3 convention
  if (addZ) {
    vbin.push_back("zjet");
    vbin.push_back("zjav");
    vbin.push_back("jetz");
  }
  // Run 3 extra
  if (addGam) {
    vbin.push_back("pjet");
    vbin.push_back("pjav");
    vbin.push_back("jetp");
  }

  map<string,const char*> label;
  label["hdm_mpfchs1"] = "HDM response";
  label["mpfchs1"] = "MPF response";
  label["mpf1"] = "MPF jet1";
  label["mpfn"] = "MPF jetN";
  label["mpfu"] = "MPF unclustered";
  label["mpfnu"] = "MPF jetN+k#timesU";
  label["zpt"] = "Z p_{T}";
  label["jetpt"] = "Jet p_{T}";
  label["ptave"] = "Avg. p_{T}";
  label["zjet"] = "Z p_{T}";
  label["jetz"] = "Jet p_{T}";
  label["zjav"] = "Avg. p_{T}";
  label["pjet"] = "#gamma p_{T}";
  label["jetp"] = "Jet p_{T}";
  label["pjav"] = "Avg. p_{T}";
  //label["counts"] = "Counts (%) per GeV";
  label["counts"] = "Counts (a.u.) per GeV";

  map<string,double> ymin;
  ymin["hdm_mpfchs1"] = 0.75;//0.95;
  ymin["mpfchs1"] = 0.75;//0.95;
  ymin["mpf1"] = 0.75;
  ymin["mpfn"] = -0.16;//-0.15;
  ymin["mpfu"] = -0.04;//-0.15;//-0.05;
  ymin["mpfnu"] = -0.22;//-0.16;
  ymin["counts"] = 2e-5;

  map<string,double> ymax;
  ymax["hdm_mpfchs1"] = 1.40;//1.20;
  ymax["mpfchs1"] = 1.40;//1.20;
  ymax["mpf1"] = 1.40;
  ymax["mpfn"] = +0.25;
  ymax["mpfu"] = +0.25;//+0.20;
  ymax["mpfnu"] = +0.32;//+0.25;
  ymax["counts"] = 3;//2e4;

  map<string,int> color;
  color["zpt"] = kRed;
  color["jetpt"] = kBlue;
  color["ptave"] = kGreen+2;
  color["zjet"] = kRed;
  color["jetz"] = kBlue;
  color["zjav"] = kGreen+2;
  color["pjet"] = kOrange+2;
  color["jetp"] = kMagenta+1;
  color["pjav"] = kCyan+1;//2;

  map<string,int> dmarker;
  dmarker["zpt"] = kFullTriangleUp;
  dmarker["jetpt"] = kFullTriangleDown;
  dmarker["ptave"] = kFullCircle;
  dmarker["zjet"] = kFullTriangleUp;
  dmarker["jetz"] = kFullTriangleDown;
  dmarker["zjav"] = kFullCircle;
  dmarker["pjet"] = kFullTriangleUp;
  dmarker["jetp"] = kFullTriangleDown;
  dmarker["pjav"] = kFullDiamond;//kFullCircle;

  map<string,int> mmarker;
  mmarker["zpt"] = kOpenTriangleUp;
  mmarker["jetpt"] = kOpenTriangleDown;
  mmarker["ptave"] = kOpenSquare;//Circle;
  mmarker["zjet"] = kOpenTriangleUp;
  mmarker["jetz"] = kOpenTriangleDown;
  mmarker["zjav"] = kOpenSquare;
  mmarker["pjet"] = kOpenTriangleUp;
  mmarker["jetp"] = kOpenTriangleDown;
  mmarker["pjav"] = kOpenSquare;

  map<string, map<string, map<string, TGraphErrors*>>> mg;
  TGraphErrors *gm(0), *gd(0), *gr(0);
  TH1D *hm(0), *hd(0), *hr(0);
  TH1D *hu_(0), *hd_(0);
  lumi_136TeV = tag;
  TCanvas *c1(0);
  TLine *l = new TLine(); l->SetLineStyle(kDashed);
  for (unsigned int i = 0; i != vobs.size(); ++i) {

    string so = vobs[i];
    const char *co = so.c_str();
    hu_ = tdrHist(Form("hu_%s",co),label[co],ymin[co],ymax[co]);
    hd_ =(so=="counts" ? tdrHist(Form("hd_%s",co),"Data/MC",0.3,1.3) :
	  so=="mpfu" ? tdrHist(Form("hd_%s",co),"Data/MC",0.,1.5) :
	   tdrHist(Form("hd_%s",co),"Data-MC (%)",-6,+5));
    c1 = tdrDiCanvas(Form("c1_%s",co),hu_,hd_,8,11);
    double ptmin = hu_->GetXaxis()->GetXmin();
    double ptmax = hu_->GetXaxis()->GetXmax();

    double lsize = (addGam && addZ ? 0.030 : 0.045);
    TLegend *legd = tdrLeg(0.45,0.92-(vbin.size()+1)*(lsize-0.005),0.75,0.92);
    legd->SetTextSize(lsize);
    legd->SetHeader("  Data");
    TLegend *legm = tdrLeg(0.40,0.92-(vbin.size()+1)*(lsize-0.005),0.70,0.92);
    legm->SetTextSize(lsize);
    legm->SetHeader("MC");

    c1->cd(1);
    gPad->SetLogx();
    l->DrawLine(ptmin,1,ptmax,1);
    l->DrawLine(ptmin,0,ptmax,0);

    c1->cd(2);
    gPad->SetLogx();
    if (so=="counts")    l->DrawLine(ptmin,1,ptmax,1);
    else if (so=="mpfu") l->DrawLine(ptmin,1,ptmax,1);
    else                 l->DrawLine(ptmin,0,ptmax,0);

    for (unsigned int j = 0; j != vbin.size(); ++j) {

      string sb = vbin[j];
      const char *cb = sb.c_str();

      // Separate file for each binning
      //TFile *f = new TFile(Form("rootfiles/ZJetPtAveV2_%s_20230824_v57/jecdataRun2Test.root",cb),"READ"); tag="Run2Test_v57";
      //TFile *f = new TFile(Form("rootfiles/ZJetPtAveV2_%s_20230825_v53/jecdataRun2Test.root",cb),"READ"); tag="Run2Test_v53";
      // Same file for each binning
      TFile *f = new TFile("rootfiles/jecdata2025CDEFG.root","READ"); tag="2025";
      assert(f && !f->IsZombie());
      curdir->cd();

      string sd = Form("data/eta00-13/%s_%s",co,cb);
      //string sd = Form("data/eta_00_13/%s_%s_a100",co,cb);
      //string sd = Form("data/eta00-13/%s_zjet_a100",co);
      //if (so=="hdm_mpfchs1") sd = Form("data/eta00-13/%s_zjet",co);
      if (debug) cout << sd << endl << flush;
      //gd = (TGraphErrors*)f->Get(sd.c_str()); assert(gd);
      TObject *od = f->Get(sd.c_str());
      if (!od) {
	sd = Form("data/eta00-13/%s_%s_a100",co,cb);
	if (debug) cout << " => " << sd << endl << flush;
	od = f->Get(sd.c_str());
      }
      if (!od && so=="mpfnu") {
	if (debug) cout << " => => merge mpfn+mpfu" << endl << flush;
	TGraphErrors *gn = mg["data"][sb]["mpfn"]; assert(gn);
	TGraphErrors *gu = mg["data"][sb]["mpfu"]; assert(gu);
	TGraphErrors *gnu = addGraphs(gn,gu,1.,1./Ru_d);
	od = (TObject*)gnu;
      }
      assert(od);
      if (od->InheritsFrom("TGraphErrors")) {
	gd = (TGraphErrors*)od;
	hd = 0;
      }
      else if (od->InheritsFrom("TH1D")) {
	hd = (TH1D*)od;
	gd = new TGraphErrors(hd);//(TH1D*)od);
      }
      else
	assert(false);

      string sm = Form("mc/eta00-13/%s_%s",co,cb);
      //string sm = Form("mc/eta_00_13/%s_%s_a100",co,cb);
      //string sm = Form("mc/eta00-13/%s_zjet_a100",co);
      //if (so=="hdm_mpfchs1") sm = Form("mc/eta00-13/%s_zjet",co);
      if (debug) cout << sm << endl << flush;
      //gm = (TGraphErrors*)f->Get(sm.c_str()); assert(gm);
      TObject *om = f->Get(sm.c_str()); //assert(om);
      if (!om) {
	sm = Form("mc/eta00-13/%s_%s_a100",co,cb);
	if (debug) cout << " => " << sm << endl << flush;
	om = f->Get(sm.c_str());
      }
      if (!om && so=="mpfnu") {
	if (debug) cout << " => => merge mpfn+mpfu" << endl << flush;
	TGraphErrors *gn = mg["mc"][sb]["mpfn"]; assert(gn);
	TGraphErrors *gu = mg["mc"][sb]["mpfu"]; assert(gu);
	TGraphErrors *gnu = addGraphs(gn,gu,1.,1./Ru_m);
	om = (TObject*)gnu;
      }
      assert(om);
      if (om->InheritsFrom("TGraphErrors")) {
	gm = (TGraphErrors*)om;
	hm = 0;
      }
      else if (om->InheritsFrom("TH1D")) {
	hm = (TH1D*)om;
	gm = new TGraphErrors(hm);//(TH1D*)om);
      }
      else
	assert(false);

      if (so=="counts") {
	//TH1D *hm = (TH1D*)gm;
	//TH1D *hd = (TH1D*)gd;
	assert(hd);
	assert(hm);
	hm = (TH1D*)hm->Clone(Form("mc_%s_%s",co,cb));
	hd = (TH1D*)hd->Clone(Form("data_%s_%s",co,cb));
	int i1 = hm->GetXaxis()->FindBin(180.);//50.);
	int i2 = hm->GetXaxis()->FindBin(1000.)-1;
	double k = 0.01;
	hd->Scale(k*100./hd->Integral(i1,i2),"width");
	hm->Scale(k*100./hm->Integral(i1,i2),"width");
	gd = new TGraphErrors(hd);
	gm = new TGraphErrors(hm);	
      }
      if (sb=="pjet" || sb=="jetp" || sb=="pjav") {
	cleanGraph(gd,40,3600);
	cleanGraph(gm,40,3600);
      }

      mg["data"][sb][so] = gd;
      mg["mc"][sb][so] = gm;
      
      /*
      if (so=="hdm_mpfchs1") {
	TH1D *hm = (TH1D*)gm;
	TH1D *hd = (TH1D*)gd;
	hm = (TH1D*)hm->Clone(Form("mc_%s_%s",co,cb));
	hd = (TH1D*)hd->Clone(Form("data_%s_%s",co,cb));
	gd = new TGraphErrors(hd);
	gm = new TGraphErrors(hm);	
      }
      */
      
      c1->cd(1);
      tdrDraw(gm,"Pz",mmarker[cb],color[cb]);
      tdrDraw(gd,"Pz",dmarker[cb],color[cb]);
      legd->AddEntry(gd,label[cb],"PE");
      legm->AddEntry(gm," ","PE");

      // Add function to estimate FSR and UE contributions
      // Normally rhoUE=2.0 GeV at reco, so AK4 with area 0.5 gets 1 GeV
      if (so=="mpfu") {
	TF1 *f1 = new TF1(Form("f1_%s_%s",co,cb),"[0]*(1+[1]*log(x))/x-1./x",60,1000);
	f1->SetParameters(15.*0.130,(0.130-0.118)/log(91.2));
	gd->Fit(f1,"QRN");
	f1->SetLineColor(color[cb]);
	f1->SetRange(ptmin,ptmax);
	f1->Draw("SAME");

	// Reference UE-only line
	TF1 *f2 = new TF1(Form("f2_%s_%s",co,cb),"-1./x",ptmin,ptmax);
	f2->SetLineColor(kBlack);
	f2->Draw("SAME");
      }
      
      c1->cd(2);
      if (so=="counts" || so=="mpfu")
	gr = tools::ratioGraphs(gd,gm);
      else
	gr = tools::diffGraphs(gd,gm,100,100);
      
      if (so=="statistics_rmpf") {
	gr = tools::ratioGraphs(gd,gm);
	hd_->SetYTitle("Data / MC");
	hd_->GetYaxis()->SetRangeUser(0.5,1.5);
	l->DrawLine(ptmin,1,ptmax,1);
      }

      tdrDraw(gr,"Pz",dmarker[cb],color[cb]);
    } // for j in vbin

    string stag = (tag!="" ? ("_"+tag) : "");
    const char *ct = stag.c_str();

    string s = "";
    if (addGam  && !addZ) s = "_Gam";
    if (!addGam &&  addZ) s = "_Zmm";
    const char *cs = s.c_str();
    
    c1->SaveAs(Form("pdf/drawZJetPtAveV2/drawZJetPtAveV2_%s%s%s.pdf",co,ct,cs));

    if (so=="statistics_rmpf" || so=="counts") {
      c1->cd(1);
      gPad->SetLogy();
      hu_->GetYaxis()->SetRangeUser(2e-8,9e3);
      c1->SaveAs(Form("pdf/drawZJetPtAveV2/drawZJetPtAveV2_%s%s%s_log.pdf",co,ct,cs));
      gPad->SetLogy(kFALSE);
    }
  } // for i in vobs


} // void drawZJetPtAve_bins
