// Purpose: Perform Run 3 or Run 2 Legacy global fit
//
// Pre-requisites:
// - reprocess.C : produce rootfiles/jecdata[X].root input file
// - softrad3.C : produce FSR+ISR corrections for HDM (MPF+DB)
// - globalFitSyst.C : produce uncertainty shapes
// - globalFitSettings.h : input definitions and configurable settings
// Post-processing:
// - globalFitPulls.C : plot pull distributions
// - createL2L3ResTextFile.C : produce simplified (L2)L3Res text file
// [- minitools/mergerL2L3ResTextFiles.C : combine L3Res with L2Res]
//
// Author: Mikko Voutilainen
//
// Notes: enable systematic source offsetting?
#include "TFile.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TLine.h"
#include "TProfile.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <set>

#include "tools.C"
#include "tdrstyle_mod22.C"
#include "globalFitSettings.h"

//using namespace globalFit;
using namespace std;
const bool debug = true;
const bool saveROOT = false;

// Helper functions to draw fit uncertainty band for arbitrary TF1 
Double_t fitError(Double_t *x, Double_t *p);
Double_t jesFit(Double_t *x, Double_t *p);
void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag);
void cleanGraph(TGraphErrors *g);
void globalFitEtaBin(double etamin, double etamax, string run, string version);
void globalFitDraw(string run, string version);

// Define global variables used in fitError
TF1 *_fitError_func(0);                 // fitError uncertainty function and
TMatrixD *_fitError_emat(0);            // error matrix

// Define global variables used in jesFit and jesFitter
// The structs are defined in globalFitSettings.h
vector<fitData> _vdt;                   // jesFitter input data,
set<string> sources;
map<string, vector<fitSyst> > _msrc;    // data->sources,
map<string, vector<fitShape> > _mshape; // obs->shapes,
string _obs;                            // data type switch,
TH1D *_hjesref(0);                      // and reference JES
int cnt(0), Nk(0);
TF1 *_jesFit(0);                        // JES fit used in jesFitter

// More small helper functions for plotting
void scaleGraph(TGraph *g, double k) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], k * g->GetY()[i]);
    if (g->InheritsFrom("TGraphErrors"))
      ((TGraphErrors*)g)->SetPointError(i, g->GetEX()[i], k * g->GetEY()[i]);
  }
} // scaleGraph
void scaleGraph(TGraph *g, TH1D *h) {
  for (int i = 0; i != g->GetN(); ++i) {
    double k = h->Interpolate(g->GetX()[i]);
    g->SetPoint(i, g->GetX()[i], k * g->GetY()[i]);
    if (g->InheritsFrom("TGraphErrors"))
      ((TGraphErrors*)g)->SetPointError(i, g->GetEX()[i], k * g->GetEY()[i]);
  }
} // scaleGraph
void shiftGraph(TGraph* g, double dy) {
  for (int i = 0; i != g->GetN(); ++i) {
    g->SetPoint(i, g->GetX()[i], g->GetY()[i]+dy);
  }
} // shiftGraph
void scaleHist(TH1D *h, TH1D *hs) {
  assert(h); assert(hs);
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    int j = hs->FindBin(h->GetBinCenter(i));
    h->SetBinContent(i, h->GetBinContent(i) * hs->GetBinContent(j));
  }
} // scaleHist
void shiftHist(TH1D* h, double dy) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    h->SetBinContent(i, h->GetBinContent(i) + dy);
  }
} // shiftHist

// Call global fit for each eta bin separately
void globalFit(string run = "All", string version = "vX") {

  globalFitEtaBin(0.0, 1.3, run, version);
} // globalFit

void globalFitEtaBin(double etamin, double etamax, string run, string version) {

  // Set fancy plotting style (CMS TDR style)
  setTDRStyle();

  // Keep track of current working directory so plots don't disappear
  TDirectory *curdir = gDirectory;
  

  // 1. Load input data sets
  ///////////////////////////
  if (debug) cout << "Reading in active datasets" << endl << flush;

  // Open file created by minitools/runAllIOVs.py and recombine.C
  const char *crun = run.c_str();
  const char *cv = version.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",crun),"UPDATE");
  assert(f && !f->IsZombie());

  // Prepare links to all relevant subdirectories
  f->cd("ratio");
  TDirectory *dratio = gDirectory;
  dratio->cd(Form("eta%02.0f-%02.0f", 10.*etamin, 10*etamax));
  TDirectory *deta = gDirectory;
  deta->cd(Form("sys"));
  TDirectory *dsys = gDirectory;
  deta->cd(Form("fsr"));
  TDirectory *dfsr = gDirectory;
  curdir->cd();

  // Load reference JES and uncertainty
  _hjesref = (TH1D*)deta->Get("herr_l2l3res"); assert(_hjesref);
  _hjesref = (TH1D*)_hjesref->Clone("hjesref");
  TH1D *herr = (TH1D*)deta->Get("herr"); assert(herr);

  // Set whitelists for quickly selecting only subset of datasets or shapes
  set<string> whitelist;
  for (unsigned int i = 0; i != _gf_datasets_whitelist.size(); ++i) {
    if (_gf_datasets_whitelist[i]!="")
      whitelist.insert(_gf_datasets_whitelist[i]);
  }
  set<string> whitelistshape;
  for (unsigned int i = 0; i != _gf_shapes_whitelist.size(); ++i) {
    if (_gf_shapes_whitelist[i]!="")
      whitelistshape.insert(_gf_shapes_whitelist[i]);
  }

  // Create listing of all active datasets
  set<string> datasets;
  for (unsigned int i = 0; i != _gf_datasets.size(); ++i) {

    // Read in information from globalFitSettings.h
    const char *name  = _gf_datasets[i][0].c_str();
    const char *type  = _gf_datasets[i][1].c_str();
    const char *name2 = _gf_datasets[i][2].c_str();

    // Check if input dataset is whitelisted
    if (!whitelist.empty() && whitelist.find(name)==whitelist.end()) continue;
    if (string(name)=="") continue; // missing elements in dataset array
    
    // Retrieve graph for input data
    TGraphErrors *g = (TGraphErrors*)deta->Get(name);
    if (!g) cout << "Input " << name << " not found." << endl << flush;
    assert(g);

    if (debug) cout << "..." << name << ": " << type << endl << flush;

    // Multijet special for crecoil
    TGraphErrors *g2(0);
    if (string(name2)!="") {
      g2 = (TGraphErrors*)deta->Get(name2);
      assert(g2);
    }
    
    // Patch HDM input from TH1D to graph [temporary]
    // => fix in softrad3.C and globalFitSyst.C (for hadw)
    if (TString(name).Contains("hdm_")) {
      TH1D *h = (TH1D*)deta->Get(name);
      assert(h);
      g = new TGraphErrors(h);
      cleanGraph(g);
    }

    // Undo previous L2L3Res, if so requested
    if (string(type)=="Rjet" && _gf_undoJESref) {
      if (debug) cout << "...undoing JES for " << type << endl << flush;

      // Special treatment for multijet
      if (TString(name).Contains("multijet")) {

	for (int i = 0; i != g->GetN(); ++i) {
	  double pt = g->GetX()[i];
	  assert(g2);
	  double ptref = g2->GetY()[i] * pt;
	  double ijes = _hjesref->FindBin(pt);
	  double jes = _hjesref->GetBinContent(ijes);
	  double ijesref = _hjesref->FindBin(ptref);
	  double jesref = _hjesref->GetBinContent(ijesref);
	  double k = jes / jesref;
	  g->SetPoint(i, g->GetX()[i], k * g->GetY()[i]);
	  g->SetPointError(i, g->GetEX()[i], k * g->GetEY()[i]);
	}
      }
      else  {
	scaleGraph(g, _hjesref);
      }
    }

    // Create fitData object
    fitData data;
    data.name = name;
    data.type = type;
    data.input = (TGraphErrors*)g->Clone(Form("%s_in",name));
    data.output = (TGraphErrors*)g->Clone(Form("%s_out",name));

    // Multijet special
    if (g2) { //string(name2)!="") {
      data.input2 = (TGraphErrors*)g2->Clone(Form("%s_in2",name2));
      data.output2 = (TGraphErrors*)g2->Clone(Form("%s_out2",name2));
    }

    // Save fitData for list and vector
    datasets.insert(name);
    _vdt.push_back(data);
  } // for i in datasets


  // 2. Load input uncertainty sources
  ////////////////////////////////////

  if (debug) cout << "Loading systematics" << endl << flush; 

  // Create listing and mapping of active uncertainty sources
  //set<string> sources;
  map<string, int> msrcidx;
  for (unsigned int i = 0; i != _gf_sources.size(); ++i) {

    // Read in information from globalFitSettings.h
    const char *name      = _gf_sources[i][0].c_str();
    const char *appliesTo = _gf_sources[i][1].c_str();
    const char *histname  = _gf_sources[i][2].c_str();

    // Activate only sources that apply to input data sets
    if (datasets.find(appliesTo)==datasets.end()) continue;

    // Use running indexing for active sources
    if (sources.find(name)==sources.end()) {
      msrcidx[name] = sources.size();
    }

    if (debug) cout << "..." << name << "->" << appliesTo
		    << "(" << histname << ")" << endl << flush;

    // Retrieve histogram for shape
    TH1D *h = (TH1D*)dsys->Get(histname);
    if (!h) h = (TH1D*)dfsr->Get(histname);
    assert(h);
  
    // Create fitSyst object
    fitSyst syst;
    syst.idx = msrcidx[name];
    syst.name = name;
    syst.appliesTo = appliesTo;
    syst.hist = (TH1D*)h->Clone(Form("%s_sys",histname));

    // Save fitSyst to list and map
    sources.insert(name);
    _msrc[syst.appliesTo].push_back(syst);
  } // for i in sources  


  // 3. Load input fit shapes
  ///////////////////////////
  if (debug) cout << "Loading input fit shapes" << endl << flush;

  // Set list of positive-definite sources
  set<string> posdeflist;
  for (unsigned int i = 0; i != _gf_posdef.size(); ++i) {
      posdeflist.insert(_gf_posdef[i]);
  }

  // Create listing and mapping of active response/composition shapes
  set<string> shapes;
  map<string, int> mshapeidx;
  for (unsigned int i = 0; i != _gf_shapes.size(); ++i) {
    
    string name       = _gf_shapes[i][0];
    string appliesTo  = _gf_shapes[i][1];
    string funcstring = _gf_shapes[i][2];

    if (name=="") continue; // ignore extra empty (commented out) elements
    // Check if input dataset is whitelisted
    if (!whitelistshape.empty() &&
	whitelistshape.find(name)==whitelistshape.end()) continue;

    // Use running indexing for active shapes
    if (shapes.find(name)==shapes.end())
      mshapeidx[name] = shapes.size();

    if (debug) cout << "..." << name << "->" << appliesTo << endl << flush;
    //<< "(" << funcstring << ")" << endl << flush;

    // Create fitShape object
    fitShape shape;
    shape.idx = mshapeidx[name];
    shape.name = name;
    shape.appliesTo = appliesTo;
    shape.ispos = (posdeflist.find(name)!=posdeflist.end());
    shape.func = new TF1(Form("f1_%s_%s",name.c_str(),appliesTo.c_str()),
			 funcstring.c_str(),15.,6500.);

    // Save fitShapes to list and map
    shapes.insert(name);
    _mshape[appliesTo].push_back(shape);
  } // for i in shapes

  // Create function to plot for JES (or composition)
  double minpt = 15.;
  double maxpt = 2000;//1500.;
  int njesFit = shapes.size();
  _jesFit = new TF1("jesFit",jesFit,minpt,maxpt,njesFit);


  // 4. Perform chi2 fit
  ///////////////////////

  //const int npar = _jesFit->GetNpar();
  const int npar = shapes.size();
  const int nsrc = sources.size();//_msrc->size();
  Int_t ntot = npar+nsrc;

  cout << endl;
  cout << "Global fit has " << npar << " fit parameters and "
       << nsrc << " nuisance parameters." << endl;
  if (_gf_penalizeFitPars)
    cout << "Fit parameters have Gaussian prior" << endl;
  cout << endl;

  // Setup global chi2 fit (jesFitter is our function)
  TFitter *fitter = new TFitter(ntot);
  fitter->SetFCN(jesFitter);

  // Set parameters
  vector<double> a(ntot, 0);
  vector<string> parnames(ntot); // empty for now, to fill with shapes/sources
  for (int i = 0; i != ntot; ++i)
    fitter->SetParameter(i, parnames[i].c_str(), a[i], (i<npar ? 0.01 : 1),
			 -100, 100);

  // Run fitter (multiple times if needed)
  const int nfit = 1;
  cnt = 0;
  for (int i = 0; i != nfit; ++i)
    fitter->ExecuteCommand("MINI", 0, 0);
  TMatrixD emat(ntot, ntot);
  gMinuit->mnemat(emat.GetMatrixArray(), ntot);

  // Retrieve the chi2 the hard way
  Double_t tmp_par[ntot], tmp_err[ntot];
  TVectorD vpar(ntot);
  TVectorD verr(ntot);
  Double_t chi2_gbl(0), chi2_src(0), chi2_par(0), chi2_data(0);
  int npar_true(0), nsrc_true(0), ndt(0);
  Double_t grad[ntot];
  Int_t flag = 1;

  for (int i = 0; i != ntot; ++i) {
    tmp_par[i] = fitter->GetParameter(i);
    tmp_err[i] = fitter->GetParError(i);
    vpar[i] = fitter->GetParameter(i);
    verr[i] = fitter->GetParError(i);
  }
  jesFitter(ntot, grad, chi2_gbl, tmp_par, flag);

  for (int i = 0; i != ntot; ++i) {
    if (fabs(tmp_par[i])!=0 || fabs(tmp_err[i]-1)>1e-2) {
      if (i < npar) ++npar_true;
      else          ++nsrc_true;
    }
    if (i < npar) chi2_par += pow(tmp_par[i],2);
    else          chi2_src += pow(tmp_par[i],2);
  }
  for (unsigned int i = 0; i != _vdt.size(); ++i) {
    TGraphErrors *gout = _vdt[i].output;
    _obs = _vdt[i].type; // for _jesFit
    for (int j = 0; j != gout->GetN(); ++j) {
      double x = gout->GetX()[j];
      double y = gout->GetY()[j];
      double ey = gout->GetEY()[j];
      chi2_data += pow((y - _jesFit->Eval(x)) / ey, 2);
      ++ndt;
    }
  }
  
  cout << endl;
  cout << "Listing global fit results for " << run << endl;
  cout << Form("Used %d data points, %d fit parameters and %d nuisances.\n"
	       "Data chi2/NDF = %1.1f / %d [%1.0f,%1.0f]\n"
	       "Nuisance chi2/Nsrc = %1.1f / %d\n"
	       "Parameter chi2/Npar = %1.1f / %d\n"
	       "Total chi2/NDF = %1.1f / %d\n",
	       ndt, npar, nsrc_true, chi2_data, ndt - npar,
	       _jesFit->GetXmin(), _jesFit->GetXmax(),
	       chi2_src, nsrc_true,
	       chi2_par, npar_true,
	       chi2_gbl, ndt - npar);

  cout << "Listing shapes (for Rjet):" << endl;
  vector<fitShape> &v = _mshape["Rjet"];
  assert(int(v.size())==njesFit);
  for (unsigned int i = 0; i != v.size(); ++i) {
    cout << Form("  %5s : %+5.2f +/- %5.2f", v[i].name.c_str(),
		 vpar[v[i].idx], verr[v[i].idx]) << endl;
  } // for i in _mshape

  cout << "Listing sources" << endl;
  map<string, vector<fitSyst> >::const_iterator it;
  for (it = _msrc.begin(); it != _msrc.end(); ++it) {
    string name = it->first;
    cout << " for " << name << ": " << endl;
    vector<fitSyst> &vs = _msrc[name];
    //assert(int(vs.size())==sources.size());
    for (unsigned int i = 0; i != vs.size(); ++i) {
      cout << Form("  %5s : %+5.2f +/- %5.2f", vs[i].name.c_str(),
		   vpar[vs[i].idx], verr[vs[i].idx]) << endl;
    }
  } // for it in _msrc
  

  // 5. Save results
  //////////////////
  deta->mkdir("run3");
  deta->cd("run3");
  
  // Fit function(s) and error matrix
  _jesFit->SetRange(10.,6500.); // nice range
  _jesFit->SetNpx(6490.); // dense binning for log scale
  const int nobs = 4;
  const char *obs[nobs] = {"Rjet", "chf",    "nhf", "nef"};
  const int color[nobs] = {kBlack,  kRed, kGreen+2, kBlue};
  for (int i = 0; i != nobs; ++i) {
    
    _obs = obs[i]; _jesFit->SetLineColor(color[i]);
    _jesFit->Write(Form("jesFit_%s",_obs.c_str()),TObject::kOverwrite);
  }
  emat.Write("emat",TObject::kOverwrite);

  // Fit functions as histograms with error band
  double vx[] = { /*1, 5, 6, 8,*/ 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const double nx = sizeof(vx)/sizeof(vx[0])-1;

  for (int i = 0; i != nobs; ++i) {

    _obs = obs[i];
    TH1D *h = new TH1D(Form("hFit_%s",_obs.c_str()),
		       Form(";p_{T} (GeV);%s",_obs.c_str()),nx,vx);
    _fitError_func = _jesFit;
    _fitError_emat = &emat;
    double k = 1;

    for (int j = 1; j != h->GetNbinsX()+1; ++j) {
      double pt = h->GetBinCenter(j);
      h->SetBinContent(j, pt, _jesFit->Eval(pt));
      h->SetBinError(j, 0., fitError(&pt, &k) - _jesFit->Eval(pt));
    } // for j in h

    h->SetLineColor(color[i]);
    h->SetMarkerStyle(kFullCircle);
    h->Write(Form("hFit_%s",_obs.c_str()),TObject::kOverwrite);
  } // for i in nobs
  
  // Fit functions as graphs with error band
  for (int i = 0; i != nobs; ++i) {
    
    _obs = obs[i];
    TGraph *gr = new TGraph(_jesFit); // use graph to keep '_obs' setting
    TGraphErrors *gre = new TGraphErrors(gr->GetN());
    _fitError_func = _jesFit;
    _fitError_emat = &emat;
    double k = 1;
    for (int i = 0; i != gr->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gre->SetPoint(i, pt, gr->GetY()[i]);
      gre->SetPointError(i, 0., fitError(&pt, &k) - gr->GetY()[i]);
    }
    delete gr;
    gre->SetLineColor(color[i]);
    gre->SetMarkerStyle(kNone);
    gre->Write(Form("gFit_%s",_obs.c_str()),TObject::kOverwrite);
  }

  f->Close();
  curdir->cd();

  globalFitDraw(run, version);
  if (debug) cout << "Finishing code" << endl << flush;
  exit(0); // avoid pageful of THashList::Delete errors
} // globalFit

// Separate drawing to factorize code and enable direct redrawing from file
void globalFitDraw(string run, string version) {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Load globalFit fit results
  const char *crun = run.c_str();
  const char *cv = version.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",run.c_str()),"READ");
  assert(f && !f->IsZombie());

  curdir->cd();
  
  const char *p = "ratio/eta00-13/run3";
  TH1D *herr = (TH1D*)f->Get("ratio/eta00-13/herr"); assert(herr);
  TMatrixD *pemat = (TMatrixD*)f->Get(Form("%s/emat",p)); assert(pemat);


  // 5. Draw results
  //////////////////
  bool drawResults = true;
  bool usingMu(false);
  if (drawResults) {

    // Graphical settings in globalFitStyles.h
    #include "globalFitStyles.h"
    curdir->cd();

    // Create canvas
    lumi_136TeV = Form("globalFit.C(\"%s\")",run.c_str());
    if (run=="Run3") lumi_136TeV = "Run3, 64 fb^{-1}";
    //TH1D *h = tdrHist("h","JES",0.982+1e-5,1.025-1e-5); // ratio (hdm)
    TH1D *h = tdrHist("h","JES",
		      0.75+1e-5,1.20-1e-5, // Summer23
		      //0.75+1e-5,1.20-1e-5, // Summer22
		      //0.83+1e-5,1.15-1e-5,
		      //0.88+1e-5,1.15-1e-5,
		      //0.88+1e-5,1.05-1e-5,
		      "p_{T} (GeV)",15,4500); // ratio (hdm)
    string epoch = string(crun);
    //if (epoch=="Run22C" || epoch=="Run22D" ||
    //	epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G") {
    //h->SetMaximum(1.12-1e-5);
    //h->SetMinimum(0.88+1e-5);
    //}
    //if (epoch=="Run22F" || epoch=="Run22G") {
    //h->SetMaximum(1.15-1e-5);
    //h->SetMinimum(0.91+1e-5);
    //}
    TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
    gPad->SetLogx();

    TLine *l = new TLine();
    //TLegend *leg = tdrLeg(0.60,0.90,0.80,0.90);
    TLegend *leg = tdrLeg(0.50,0.90,0.70,0.90);
    //TLegend *leg2 = tdrLeg(0.45,0.15,0.65,0.30);
    TLegend *leg2 = tdrLeg(0.45,0.15,0.65,0.25);

    // Draw fit on the back
    assert(_jesFit);
    _jesFit->SetRange(15.,4500.); // nice range
    _jesFit->SetNpx(4485.); // dense binning for log scale
    _obs = "Rjet";
    TGraph *gr = new TGraph(_jesFit); // use graph to keep '_obs' setting
    TGraphErrors *gre = new TGraphErrors(gr->GetN());
    _fitError_func = _jesFit;
    _fitError_emat = pemat;
    double k = 1;
    for (int i = 0; i != gr->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gre->SetPoint(i, pt, gr->GetY()[i]);
      gre->SetPointError(i, 0., fitError(&pt, &k) - gr->GetY()[i]);
    }
    if (debug) cout << "Draw JES" << endl << flush;
    tdrDraw(herr,"E3",kFullCircle,kCyan+2,kSolid,-1,1001,kCyan+1);

    // Add reference JEC for Run 2 UL reruns
    if (TString(epoch).Contains("UL")) {
      assert(_hjesref);
      tdrDraw(_hjesref,"E3",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
      _hjesref->SetFillColorAlpha(kOrange+1,0.7);
      TH1D *hjesrefl = (TH1D*)_hjesref->Clone("hjesrefl");
      tdrDraw(hjesrefl,"HIST",kNone,kOrange+2,kSolid,-1,kNone);
      //leg2->SetY2(leg2->GetY2()+0.05);
      leg2->AddEntry(_hjesref,"Ref. JES.","LF");

      TFile *f2(0);
      if (epoch=="2016ULAPV")
	f2 = new TFile("../jecsys2020/rootfiles/jecdata2016BCDEF.root","READ");
      if (epoch=="2016UL")
	f2 = new TFile("../jecsys2020/rootfiles/jecdata2016GH.root","READ");
      if (epoch=="2017UL")
	f2 = new TFile("../jecsys2020/rootfiles/jecdata2017BCDEF.root","READ");
      if (epoch=="2018UL")
	f2 = new TFile("../jecsys2020/rootfiles/jecdata2018ABCD.root","READ");
      assert(f2);
      TH1D *hrefit = (TH1D*)f2->Get("ratio/eta00-13/sys/hjesfit2"); // UL refit
      assert(hrefit);
      TH1D *hrefitl = (TH1D*)hrefit->Clone("hrefitl");
      tdrDraw(hrefit,"E3",kNone,kRed+1,kSolid,-1,1001,kRed);
      tdrDraw(hrefitl,"HIST",kNone,kRed+1,kSolid,-1,kNone);
      hrefit->SetFillColorAlpha(kRed,0.7);
      leg2->SetY2(leg2->GetY2()+0.05);
      leg2->AddEntry(hrefit,"UL refit","LF");
    }

    //if (epoch!="Run22F") {
    tdrDraw(gre,"E3",kNone,kBlack,kSolid,-1,1001,kYellow+1);
    tdrDraw(gr,"Lz",kNone,kBlack);
    //}
    l->SetLineStyle(kDashed);
    l->DrawLine(15,1,4500,1);
    
    //if (!TString(epoch).Contains("UL"))
    //leg2->AddEntry(l,"Run 3 avg.","L");
    //leg2->AddEntry(herr,"Total unc.","F");
    leg2->AddEntry(herr,"Run2 total unc.","F");
    leg2->AddEntry(gre,"Fit unc.","FL");

    // Separate canvas for CHF, NHF, NEF
    //TH1D *hc = tdrHist("hc","PF composition (0.01)",-2+1e-5,2-1e-5);
    //TH1D *hc = tdrHist("hc","PF composition (0.01)",-10+1e-4,10-1e-4,
    TH1D *hc = tdrHist("hc","PF composition (0.01)",-15+1e-4,15-1e-4,
		       "p_{T} (GeV)",15,4500);
    //hc->GetYaxis()->SetRangeUser(-4.5+1e-5,4.5-1e-5);
    TCanvas *c1c = tdrCanvas("c1c",hc,8,11,kSquare);
    gPad->SetLogx();
    l->DrawLine(15,0,4500,0);
    TLegend *legc = tdrLeg(0.45,0.90,0.65,0.90);

    // Separate canvas for CEF, MUF
    TH1D *hl = tdrHist("hl","PF composition (0.01)",-2e-3+1e-7,2.5e-3-1e-7);
    TCanvas *c1l = tdrCanvas("c1l",hl,8,11,kSquare);
    gPad->SetLogx();
    l->DrawLine(15,0,3500,0);

    // Sanity check PF composition sums
    TGraphErrors *gpfjet(0);
    TGraphErrors *gzjet(0);
    TGraphErrors *gzljet(0);
    TGraphErrors *gzmjet(0);
    TGraphErrors *ggjet(0);
  
    if (debug) cout << "Draw data" << endl << flush;
    for (unsigned int i = 0; i != _vdt.size(); ++i) {
      
      TGraphErrors *gi = (TGraphErrors*)_vdt[i].input->Clone(Form("gi%d",i));
      TGraphErrors *go = (TGraphErrors*)_vdt[i].output->Clone(Form("go%d",i));
      string name = _vdt[i].name;
      string type = _vdt[i].type;
      if (debug) cout << "..." << name << ": " << type << endl << flush;

      if (type=="Rjet") {
        c1->cd();
	leg->AddEntry(go,_gf_label[name],"PLE");
	leg->SetY1NDC(leg->GetY1NDC()-0.05);
      }
      else if (type=="cef" || type=="muf") {
	usingMu = true;
	c1l->cd();
      }
      else {
	c1c->cd();
	if (_gf_label[name] && string(_gf_label[name])!="X") {
	  legc->AddEntry(go,_gf_label[name],"PLE");
	  legc->SetY1NDC(legc->GetY1NDC()-0.05);
	}
	// Turn into percentages
	scaleGraph(gi,100);
	scaleGraph(go,100);
      }
      
      // Default settings
      if (_gf_color[name]==0)  _gf_color[name] = kBlack;
      if (_gf_marker[name]==0) _gf_marker[name] = kFullCircle;
      if (_gf_label[name]==0)  _gf_label[name] = name.c_str();
      if (_gf_size[name]==0)   _gf_size[name] = 1.0;
      
      tdrDraw(go,"Pz",_gf_marker[name],_gf_color[name]);
      if (name=="hdm_mpfchs1_multijet" || name=="mpfchs1_multijet_a100" ||
	  name=="ptchs_multijet_a100") {
	//if (TString(name.c_str()).Contains("multijet")) {
	tdrDraw(gi,"Pz",kOpenTriangleUp,_gf_color[name]);
      }
      //if (name=="hdm_cmb_mj" || (run=="2017H" && name=="hdm_cmb")) {
      if (name=="mpfchs1_zjet_a100" || name=="hdm_mpfchs1_zjet" ||
	  name=="ptchs_zjet_a100" || //) {
	  name=="mpfchs1_gamjet_a100" || name=="hdm_mpfchs1_gamjet" ||
	  name=="ptchs_gamjet_a100") {
	c1c->cd();
	TGraphErrors *gr = (TGraphErrors*)go->Clone("gr");
	assert(_hjesref);
	assert(gr);
	if (_gf_useJESref) scaleGraph(gr,_hjesref);
	shiftGraph(gr,-1);
	scaleGraph(gr,100);
	//tdrDraw(gr,"Pz",kFullCircle,kBlack);
	tdrDraw(gr,"Pz",
		(_gf_marker[name+"2"] ? _gf_marker[name+"2"] : kFullCircle),
		(_gf_color[name+"2"] ? _gf_color[name+"2"] : kBlack));
	if (string(_gf_label[name+"2"])!="X") {
	  legc->AddEntry(gr,_gf_label[name+"2"],"PLE");
	  legc->SetY1NDC(legc->GetY1NDC()-0.05);
	}
      }

      go->SetMarkerSize(_gf_size[name]);
      gi->SetMarkerSize(_gf_size[name]);      
    } // for i in _vdt   

    if (debug) cout << "Draw EFL, MUF compositions" << endl << flush;

    c1l->cd();
    _obs = "cef";
    TGraph *gcef = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    tdrDraw(gcef,"Lz",kNone,kCyan+2,kSolid);

    _obs = "muf";
    TGraph *gmuf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    tdrDraw(gmuf,"Lz",kNone,kMagenta+2,kSolid);

    if (debug) cout << "Draw CHF, NHF, NEF compositions" << endl << flush;

    c1c->cd();

    _obs = "chf";
    TGraph *gchf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gchfe = new TGraphErrors(gchf->GetN());
    for (int i = 0; i != gchf->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gchfe->SetPoint(i, pt, gchf->GetY()[i]);
      gchfe->SetPointError(i, 0., fitError(&pt, &k) - gchf->GetY()[i]);
    }

    scaleGraph(gchf,100);
    tdrDraw(gchf,"Lz",kNone,kRed,_gf_fitPFcomp ? kSolid : kDotted);

    _obs = "nhf";
    TGraph *gnhf = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gnhfe = new TGraphErrors(gnhf->GetN());
    for (int i = 0; i != gnhf->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gnhfe->SetPoint(i, pt, gnhf->GetY()[i]);
      gnhfe->SetPointError(i, 0., fitError(&pt, &k) - gnhf->GetY()[i]);
    }

    scaleGraph(gnhf,100);
    tdrDraw(gnhf,"Lz",kNone,kGreen+2,_gf_fitPFcomp ? kSolid : kDotted);

    _obs = "nef";
    TGraph *gnef = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    TGraphErrors *gnefe = new TGraphErrors(gnef->GetN());
    for (int i = 0; i != gnef->GetN(); ++i) {
      double pt = gr->GetX()[i];
      gnefe->SetPoint(i, pt, gnef->GetY()[i]);
      gnefe->SetPointError(i, 0., fitError(&pt, &k) - gnef->GetY()[i]);
    }
    //tdrDraw(gnefe,"E3",kNone,kBlue+1,kSolid,-1,1001,kBlue-9);
    scaleGraph(gnef,100);
    tdrDraw(gnef,"Lz",kNone,kBlue,_gf_fitPFcomp ? kSolid : kDotted);

    _obs = "Rjet";
    TGraph *grjt = new TGraph(_jesFit);  // use graph to keep '_obs' setting
    if (_gf_useJESref) scaleGraph(grjt,_hjesref);
    shiftGraph(grjt,-1);
    scaleGraph(grjt,100);
    tdrDraw(grjt,"Lz",kNone,kBlack,_gf_fitPFcomp ? kSolid : kDotted);


    // test case: gamma+jet true response in MC
    if (false) {
      c1->cd();
      TFile *fg = new TFile("../gamjet/files/GamHistosFill_mc_2018P8.root",
			    "READ");
      assert(fg && !fg->IsZombie());
      TProfile *p = (TProfile*)fg->Get("control/prgen"); assert(p);
      TProfile *pr = (TProfile*)fg->Get("control/prjet"); assert(pr);
      TProfile *pg = (TProfile*)fg->Get("control/pgjet"); assert(pg);
      TH1D *hp = pr->ProjectionX("hp");
      hp->Divide(pg);
      tdrDraw(hp,"HIST",kNone,kCyan+2,kSolid,-1,kNone);
    }
    // test case: Z+jet true response in MC
    if (false) {
      TFile *fz = new TFile("rootfiles/jme_bplusZ_merged_v37_2016FH_mu.root","READ");
      assert(fz && !fz->IsZombie());
      // rz_zmmjet_a100 : (Zpt,genZpt/Zpt) => Z resp. log-lin slope +1% to -2%
      // rpt_zmmjet_a100 : (Zpt,RpT) => pTreco/pTZ vs pTZ
      // rbal_zmmjet_a100 : (Zpt,ljet_pt/Zpt) => pTrgen/pTZ vs pTZ
      // *rgenjet1_zmmjet_a100 : (Zpt,...) => pTgen*cos(dphi)/pTZ vs pTZ
      // *rmpfjet1_zmmjet_a100 : (Zpt,RMPFjet1) => pTreco*cos(dphi)/pTZ vs pTZ
      TGraphErrors *gr2 = (TGraphErrors*)fz->Get("mc/eta_00_13/rmpfjet1_zmmjet_a100"); assert(gr2);
      TGraphErrors *gg2 = (TGraphErrors*)fz->Get("mc/eta_00_13/rgenjet1_zmmjet_a100"); assert(gg2);
      TGraphErrors *g2 = tools::ratioGraphs(gr2,gg2);
      TGraphErrors *gr = (TGraphErrors*)fz->Get("mc/eta_00_13/rpt_zmmjet_a100"); assert(gr);
      TGraphErrors *gg = (TGraphErrors*)fz->Get("mc/eta_00_13/rbal_zmmjet_a100"); assert(gg);
      TGraphErrors *g = tools::ratioGraphs(gr,gg);
      tdrDraw(g,"Lz",kNone,kRed+2,kSolid,-1,kNone);
      tdrDraw(g2,"Lz",kNone,kRed+2,kDashed,-1,kNone);
      // rz may be relevant, HDM is above MC truth at high pT
      TGraphErrors *gz = (TGraphErrors*)fz->Get("mc/eta_00_13/rz_zmmjet_a100");
      assert(gz);
      TGraphErrors *gz1 = tools::ratioGraphs(g,gz);
      tdrDraw(gz1,"Lz",kNone,kRed+2,kDotted,-1,kNone);
      // Much too large fix at high pT. Different in 2016FH? Something else?
    }

    c1->cd(); gPad->RedrawAxis();
    c1c->cd(); gPad->RedrawAxis();
    c1l->cd(); gPad->RedrawAxis();

    if (debug) cout << "Draw plots" << endl << flush;
   
    c1->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_rjet.pdf",crun,cv));
    c1c->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_pf.pdf",crun,cv));
    if (usingMu) c1l->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_mu.pdf",crun,cv));
    //
    if (saveROOT) c1->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_rjet.root",crun,cv));

    // Test
    c1->cd();
    //TGraphErrors *gnhf0 = (TGraphErrors*)f->Get("ratio/eta00-13/nhf_multijet_a100");
    TGraphErrors *gnhf0 = (TGraphErrors*)f->Get("ratio/eta00-13/nhf_incjet_a100");
    assert(gnhf0);
    gnhf0 = (TGraphErrors*)gnhf0->Clone("gnhf0");
    TGraphErrors *gnhf0z = (TGraphErrors*)f->Get("ratio/eta00-13/nhf_zjet_a100");
    assert(gnhf0z);
    gnhf0z = (TGraphErrors*)gnhf0z->Clone("gnhf0z");
    double nhf_off(1);
    /*
    // 22Sep2023
    if (run=="Run22CD") nhf_off = 1.0;
    if (run=="Run22E")  nhf_off = 2.0;
    if (run=="Run22FG")  nhf_off = 0.0;
    if (run=="Run23C123")  nhf_off = 1.0;
    if (run=="Run23C4")  nhf_off = 4.5;
    if (run=="Run23D")  nhf_off = 4.5;
    if (run=="Run23C4D")  nhf_off = 4.5;
    if (run=="Run3")  {
      nhf_off = ((5.1+3.0)*1.0 + 5.9*2.0 + (18.0+3.1)*0.0 +
		 8.7*1.0 + (9.8+9.5)*4.5) /
	((5.1+3.0) + 5.9 + (18.0+3.1) + 8.7 + (9.8+9.5));
    }
    */
    // 22Sep2023
    if (run=="Run22CD") nhf_off = 2.0;
    if (run=="Run22E")  nhf_off = 2.5;
    if (run=="Run22FG")  nhf_off = 3.0;
    if (run=="Run23C123")  nhf_off = -1.0;//1.0;
    if (run=="Run23C4")  nhf_off = 7.0;//4.0;
    if (run=="Run23D")  nhf_off = 9.0;//4.5;
    if (run=="Run3")  {
      nhf_off = ((5.1+3.0)*2.0 + 5.9*3.0 + (18.0+3.1)*3.0 +
		 //8.7*1.0 + 9.8*4.0 + 9.5*4.5) / //Summer22
		 8.7*-1.0 + 9.8*7.0 + 9.5*9.0) /  //Summer23
	((5.1+3.0) + 5.9 + (18.0+3.1) + 8.7 + (9.8+9.5));
    }
    
    for (int i = 0; i != gnhf0->GetN(); ++i) {
      double x = gnhf0->GetX()[i];
      //double ex = gnhf0->GetEX()[i];
      double y = gnhf0->GetY()[i];
      //double ey = gnhf0->GetEY()[i];
      //gnhf0->SetPoint(i, x, 1+y);
      gnhf0->SetPoint(i, x, (1+y)-nhf_off*0.01); // small offset
      //gnhf0->SetPointError(i, ex, ey);
    }
    for (int i = gnhf0z->GetN()-1; i != -1; --i) {
      double x = gnhf0z->GetX()[i];
      double y = gnhf0z->GetY()[i];
      gnhf0z->SetPoint(i, x, (1+y)-nhf_off*0.01); // small offset
      if (x>50.) gnhf0z->RemovePoint(i); // only range with no incjet
    }
    
    //leg->AddEntry(gnhf0,"Not fit: Dijet NHF","PLE");
    //leg->AddEntry(gnhf0,"Not fit: Incjet NHF-1%","PLE");
    //leg->AddEntry(gnhf0,Form("Not fit: jet NHF-%1.1f%%",nhf_off),"PLE");
      leg->AddEntry(gnhf0,Form("Not fit: NHF%+1.1f%%",nhf_off),"PLE");
    leg->SetY1NDC(leg->GetY1NDC()-0.05);
    //tdrDraw(gnhf0,"Pz",kFullSquare,kGreen+2);
    tdrDraw(gnhf0,"Pz",kOpenSquare,kGreen+2);
    tdrDraw(gnhf0z,"Pz",kOpenSquare,kGreen+3);

    // Add average JEC also on the plot
    if (run=="Run3") {
      TFile *f2 = new TFile("rootfiles/jecdataRun3Data.root","READ");
      assert(f2 && !f2->IsZombie());
      TH1D *hrjet = (TH1D*)f2->Get("ratio/eta00-13/run3/hFit_Rjet");
      assert(hrjet);
      tdrDraw(hrjet,"E3",kNone,kOrange+2,kSolid,-1,1001,kOrange+1);
      hrjet->SetFillColorAlpha(kOrange+1,0.7);
      leg->AddEntry(hrjet,"Avg. JES","F");
      leg->SetY1NDC(leg->GetY1NDC()-0.05);
    }

    if (string(crun)=="Run3")
      c1->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_rjet_wNHF.pdf",crun,cv));
    else
      c1->SaveAs(Form("pdf/globalFit/globalFit_%s_%s_rjet.pdf",crun,cv));
  } // drawResults
} // globalFitEtaBin


// Generic fit function for JES and composition
Double_t jesFit(Double_t *x, Double_t *p) {

  // Choose JES (var~1) or PF composition (var~0) as baseline
  double var = (_obs=="Rjet" ? 1. : 0.);
  double pt = x[0];

  // Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
  double jesref = 1;
  if (_gf_useJESref && _obs=="Rjet") {
    assert(_hjesref);
    jesref = _hjesref->Interpolate(pt);
  }

  // Load available shapes for this observable
  vector<fitShape> &v = _mshape[_obs];
  for (unsigned int i = 0; i != v.size(); ++i) {

    // Calculate variation and add it to total
    int idx = v[i].idx;   assert(idx>=0);
    TF1 *f1 = v[i].func;  assert(f1);
    double par = p[idx];
    if (v[i].ispos) par = max(par,0.);
    //if (v[i].ispos) par = min(1.,max(par,0.));

    var += par * f1->Eval(pt) * 0.01; // fullSimShapes in %'s
  } // for i in mshape

  return (var / jesref);
} // jesFit


// General chi2 fit function
void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par,
	       Int_t flag) {

  // Basic checks
  assert(_jesFit);

  // Parametes for nuisances (sources)
  Double_t *ps = &par[_jesFit->GetNpar()];
  int ns = npar - _jesFit->GetNpar();

  if (flag) {
    
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    Nk = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (unsigned int ig = 0; ig != _vdt.size(); ++ig) {

      string name         = _vdt[ig].name;
      string type         = _vdt[ig].type;
      TGraphErrors *gin   = _vdt[ig].input;  assert(gin);
      TGraphErrors *gin2  = _vdt[ig].input2;
      TGraphErrors *gout  = _vdt[ig].output; assert(gout);
      TGraphErrors *gout2 = _vdt[ig].output2;
      
      for (int i = 0; i != gin->GetN(); ++i) {

	// Retrieve central value and uncertainty for this point
	double pt = gin->GetX()[i];
	double data = gin->GetY()[i];
	double sigma = gin->GetEY()[i];

	// Calculate fit value at this point
	for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	  _jesFit->SetParameter(ipar, par[ipar]);
	}
	_obs = type;
	double fit = _jesFit->EvalPar(&pt,par);

	// For multijet balancing, multiply data by reference JES
	//if (TString(name.c_str()).Contains("multijet")) {
	if (name=="hdm_mpfchs1_multijet" || name=="mpfchs1_multijet_a100" ||
	    name=="ptchs_multijet_a100") {
	  assert(gin2);
	  double ptref = gin2->GetY()[i] * pt;
	  double fitRef = _jesFit->EvalPar(&ptref,par);
	  data *= fitRef;
	  sigma *= fitRef;
	} // multijet

	// For PF composition, check if we want to add this to chi2
	bool addChi2 = true;
	bool isPF = (TString(name.c_str()).Contains("chf") ||
		     TString(name.c_str()).Contains("nhf") ||
		     TString(name.c_str()).Contains("nef"));
	if (!_gf_fitPFcomp && isPF) addChi2 = false;
	    
	// For Z+jet HDM, check if we want to add this to chi2
	if (!_gf_fitZjetHDM && name=="hdm_mpfchs1_zjet") addChi2 = false;
	
	// Calculate total shift caused by all nuisance parameters
	double shifts = 0;
	vector<fitSyst> &_vsrc = _msrc[name];
	for (unsigned int is = 0; is != _vsrc.size(); ++is) {

	  TH1D *hsrc = _vsrc[is].hist; assert(hsrc);
	  int idx = _vsrc[is].idx;

	  int ipt = hsrc->FindBin(pt);
	  shifts += ps[idx] * hsrc->GetBinContent(ipt);
	}

	
	// Add chi2 from residual
	double chi = (data + shifts - fit) / oplus(sigma,globalErrMin);
	if (isPF) { // quick hack to de-emphasize PF composition a bit
	  chi = max(fabs(data+shifts-fit)-_gf_fitPFcomp_minErr,0.)
	    / oplus(sigma,oplus(globalErrMin,_gf_fitPFcomp_minErr));
	}
	if (addChi2) {
	  chi2 += chi * chi;
	  ++Nk;
	}

	// Store shifted data
	assert(gout->GetN()==gin->GetN() && gout->GetX()[i]==pt);
	gout->SetPoint(i, pt, data + shifts);

	// For multijets, store also downward extrapolation
	if (TString(name.c_str()).Contains("multijet")) {
	  double jes = _jesFit->EvalPar(&pt, par);
	  double ptref = pt * gin2->GetY()[i];
	  double jesref = _jesFit->EvalPar(&ptref, par);
	  // MJB = jes / jesref
	  // data = MJB*jesref => "jesref" = jes / MJB = jes * jesref/data
	  // NB: double check the logic here
	  gout2->SetPoint(i, ptref, jes * jesref / (data + shifts));
	  gout2->SetPointError(i, gout->GetEX()[i], gout2->GetEY()[i]);
	} // multijet
      } // for ipt
    } // for ig
  
    // Add chi2 from nuisance parameters
    for (int is = 0; is != ns; ++is) {
      
      chi2 += ps[is]*ps[is];
    } // for ipar

    // Add penalty for fit parameters (Bayesian prior, essentially)
    if (_gf_penalizeFitPars) {
      for (int ipar = 0; ipar != _jesFit->GetNpar(); ++ipar) {
	chi2 += par[ipar] * par[ipar];
	++Nk;
      }
    } // penalizeFitPars
    
    // Give some feedback on progress in case loop gets stuck
    if ((++cnt)%100==0) cout << "." << flush;
  } // if flag
  else {
    if (grad) {}; // suppress warning;
    return;
  }

} // jesFitter


// Calculate fit uncertainty for an arbitrary function                          
// using standard error propagation with differentials                          
Double_t fitError(Double_t *xx, Double_t *pp) {

  // Sanity checks on inputs                                                    
  assert(_fitError_func);
  assert(_fitError_emat);
  double x = *xx;
  double k = pp[0];
  TF1 *f = _fitError_func;
  int n = f->GetNpar();
  int nsrc = sources.size();
  TMatrixD &emat = (*_fitError_emat);
  assert(emat.GetNrows()==n+nsrc);
  assert(emat.GetNcols()==n+nsrc);

  // Partial derivatives as differentials with 10% step size                    
  vector<double> df(n);
  for (int i = 0; i != n; ++i) {

    double p = f->GetParameter(i);
    double dp = 0.1*sqrt(emat[i][i]); // step size 10% of uncertainty           
    f->SetParameter(i, p+dp);
    double fup = f->Eval(x);
    f->SetParameter(i, p-dp);
    double fdw = f->Eval(x);
    f->SetParameter(i, p);

    // Calculate partial derivative as a simple differential                    
    df[i] = (dp ? (fup - fdw) / (2.*dp) : 0);
  }

  // Perform standard error propagation                                         
  double sumerr2(0);
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
      sumerr2 += emat[i][j]*df[i]*df[j];
    }
  }

  double err = sqrt(sumerr2);

  return (f->Eval(x) + k*err);
}

void cleanGraph(TGraphErrors *g) {
  for (int i = g->GetN()-1; i != -1; --i) {
    if (g->GetY()[i]==0 && g->GetEY()[i]==0)
      g->RemovePoint(i);
  }
} // cleanGraph
