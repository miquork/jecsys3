// File: reprocess.C
// Created by Mikko Voutilainen, on Sep 6th, 2012
// Originally taken from jecsys package, heavily stripped down
// Purpose: Combine graphs from difference channels for simpler JEC analysis
//           Macro examplePlot() shows how to create plots with "CMS JEC style"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TLine.h"
#include "TProfile.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "tools.h"
#include "tdrstyle_mod22.C"

#include <string>
#include <iostream>
#include <vector>
#include <map>

using namespace std;

const bool rp_debug = true; // verbose messages
bool correctZMass = false; // pT with Run 2 Z+jet mass (def:true)

// PS weight variations, 0,1,2,3 and "" for reference
string sPSWgtZ = ""; // use 'PSWeight[X]/' variant of Z+jet MC (def:"") [PSWeight0 has very poor effective statistics]

//"L1L2": "MCTruth corrections" applied, in reality L1Data/MC,L2
//"L1L2Res": MCTruth + L2Res applied
//"L1L2L3Res": MCTruth + L2L3Res applied; reference JES central values set to 1.0 (affects plotting as well)
//string CorLevel = "L1L2";    // Rarely used
//string CorLevel = "L1L2Res"; // Input to (L2)L3Res
string CorLevel = "L1L2L3Res"; // Closure test for L2L3Res


// Settings for cleaned up global fit
/////////////////////////////////////

// Z+jet
double fzptmin(40);//15.);//40.);   // Z+jet pTmin
double fzptmax(700.);//1300.);  // Z+jet pTmax
double fzmpfptmin(15); // Z+jet MPF pTmin
double fzmpfptmax(700.);//1300);// Z+jet MPF pTmax
double fzbalptmin(15); // Z+jet DB pTmin
double fzbalptmax(700.);//1300);// Z+jet DB pTmax
double fzbptmax(300.); // Z+b pTmax

// Helper functions to handle JEC
FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
FactorizedJetCorrector *getFJC(string l1, string l2="", string res="",
			       string path="");
void setEtaPtRho(FactorizedJetCorrector *jec, double eta,double pt,double rho);
Double_t funcCorrPt(Double_t *x, Double_t *p);
double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho = 0);


// put all the different methods in a single file for easy access for everybody
void reprocess(string epoch="") {

  // Set TDR style to have correct graphical settings when storing graphs
  setTDRStyle();
 
  TDirectory *curdir = gDirectory;

  const char *cep = epoch.c_str();
  TFile *fout = new TFile(Form("rootfiles/jecdata%s.root", cep), "RECREATE");

  ////////////////////////////
  // Z+jet                  //
  ////////////////////////////

  TFile *fz(0), *fzjes(0);
  TH1D *hzjes(0);
  const char *cdz = "../JERCProtoLab/Winter22Run3/L3Residual_Z";
  if (epoch=="RunF") {
    //fz = new TFile(Form("%s/jme_ZplusJet_muon_EOS_v2p1_Run2022F_v1.root",cdz),"READ"); fz->cd("Data_2022F"); gDirectory->cd("NoCorrection"); fz = (TFile*)gDirectory; // RunF_noL2L3Res
    //fz = new TFile(Form("%s/jme_ZplusJet_muon_EOS_v2p1_Run2022F_v2.root",cdz),"READ");  // RunF_defaultMET
    fz = new TFile(Form("%s/jme_ZplusJet_Muon_Run2022F_NoCorrection_v3.root",cdz),"READ");
    //fz->cd("Data_2022F"); 
    fz->cd("Run2022F"); 
    gDirectory->cd("NoCorrection"); fz = (TFile*)gDirectory;
  }
  if (epoch=="RunE") {
    // what's 2022E_v10?
    //fz = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9.root",cdz),"READ"); // RunE_V2_L2L3Res input
    fz = new TFile(Form("%s/jme_ZplusJet_Muon_Run2022E_NoCorrection_v3.root",cdz),"READ");
    //fz->cd("Data_2022E"); 
    fz->cd("Run2022E");
    gDirectory->cd("NoCorrection"); fz = (TFile*)gDirectory;
  }
  if (epoch=="RunCD") {
    fz = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p0_Run2022CD_v6.root",cdz),"READ"); // RunC_V2_L2L3Res input
    //fz = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p0_Run2022CD_v7.root",cdz),"READ"); // RunC_V2_L2L3Res closure test
  }
  const char *cdz22 = "../JERCProtoLab/Summer22Run3/L3Residual_Z";
  if (epoch=="Run22C") {
    fz = new TFile(Form("%s/jme_bplusZ_2022C_Zmm_sync_v53.root",cdz22),"READ");
  }
  if (epoch=="Run22D") {
    fz = new TFile(Form("%s/jme_bplusZ_2022D_Zmm_sync_v53.root",cdz22),"READ");
  }
  const char *cdz22e = "../JERCProtoLab/Summer22EERun3/L3Residual_Z";
  if (epoch=="Run22E") {
    fz = new TFile(Form("%s/jme_bplusZ_2022E_Zmm_sync_v53.root",cdz22e),"READ");
  }
  if (epoch=="Run22F") {
    fz = new TFile(Form("%s/jme_bplusZ_2022F_Zmm_sync_v53.root",cdz22e),"READ");
  }
  if (epoch=="Run22G") {
    fz = new TFile(Form("%s/jme_bplusZ_2022G_Zmm_sync_v53.root",cdz22e),"READ");
  }
  const char *cdz23 = "../JERCProtoLab/Summer23Run3/L3Residual_Z";
  const char *cdzx = "rootfiles";
  if (epoch=="Run23B") {
    fz = new TFile(Form("%s/jme_bplusZ_2023B_Zmm_sync_v53.root",cdz23),"READ");
    //fz = new TFile(Form("%s/jme_bplusZ_2023B_Zmm_sync_v53-2.root",cdzx),"READ");
  }
  if (epoch=="Run23BC123") {
    fz = new TFile(Form("%s/jme_bplusZ_2023BC123_Zmm_sync_v53.root",cdzx),"READ");
  }
  if (epoch=="Run23C") {
    fz = new TFile(Form("%s/jme_bplusZ_2023C_Zmm_sync_v53.root",cdz23),"READ");
  }
  if (epoch=="Run23C1") {
    fz = new TFile(Form("%s/jme_bplusZ_2023C1_Zmm_sync_v53.root",cdzx),"READ");
  }
  if (epoch=="Run23C2") {
    fz = new TFile(Form("%s/jme_bplusZ_2023C2_Zmm_sync_v53.root",cdzx),"READ");
  }
  if (epoch=="Run23C3") {
    fz = new TFile(Form("%s/jme_bplusZ_2023C3_Zmm_sync_v53.root",cdzx),"READ");
  }
  if (epoch=="Run23C4") {
    fz = new TFile(Form("%s/jme_bplusZ_2023C4_Zmm_sync_v53.root",cdzx),"READ");
  }
  if (epoch=="Run3") {
    fz = new TFile("rootfiles/jecdataRun3Data.root","READ");
  }

  assert(fz && !fz->IsZombie());
  
  TFile *fmz = fz; 
  assert(fmz && !fmz->IsZombie());
  
  string sr = "eta_00_13";
  const char *cr = sr.c_str();
  const char *cl = CorLevel.c_str();

  TH1D *hcounts(0);
  TH1D *hmz_dt(0), *hmz_mc(0);
  if (epoch=="RunE" || epoch=="RunF") {
    //TH2D *hmz_dt2 = (TH2D*)fmz->Get("2022F/NoCorrection/DATA_ZpT_ZMass_a10_eta_00_03");
    TH2D *hmz_dt2 = (TH2D*)fmz->Get("data/eta_00_13/zpt_mass_zmmjet_a100");
    assert(hmz_dt2);
    hmz_dt = (hmz_dt2 ? hmz_dt2->ProfileX()->ProjectionX("hmz_dt") : 0);
    hcounts = hmz_dt2->ProjectionX("hcounts");
    assert(hmz_dt);

    //TH2D *hmz_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100",cr));
    TH2D *hmz_mc2 = (TH2D*)fmz->Get("mc/eta_00_13/zpt_mass_zmmjet_a100");
    assert(hmz_mc2);
    hmz_mc = (hmz_mc2 ? hmz_mc2->ProfileX()->ProjectionX("hmz_mc") : 0);
    //hmz_mc = (TH1D*)hmz_dt->Clone("hmz_mc"); hmz_mc->Divide(hmz_dt);
    //hmz_mc->Scale(91.2);
    assert(hmz_mc);
  }
  //if (epoch=="RunE" || epoch=="RunCD") {
  if (epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
      epoch=="Run23C1" ||epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4"
      ) {
    TH2D *hmz_dt2 = (TH2D*)fmz->Get(Form("data/%s/h_Zpt_mZ_alpha100",cr));
    assert(hmz_dt2);
    hmz_dt = (hmz_dt2 ? hmz_dt2->ProfileX()->ProjectionX("hmz_dt") : 0);
    assert(hmz_dt);

    TH2D *hmz_mc2 = (TH2D*)fmz->Get(Form("mc/%s/h_Zpt_mZ_alpha100",cr));
    assert(hmz_mc2);
    hmz_mc = (hmz_mc2 ? hmz_mc2->ProfileX()->ProjectionX("hmz_mc") : 0);
    assert(hmz_mc);
  }
  if (epoch=="Run3") {
    hmz_dt = (TH1D*)fz->Get("data/eta00-13/mass_zjet_a100"); assert(hmz_dt);
    hmz_mc = (TH1D*)fz->Get("mc/eta00-13/mass_zjet_a100"); assert(hmz_mc);
  }
  assert(hmz_dt);
  assert(hmz_mc);

  TH1D *hmz = (TH1D*)hmz_dt->Clone("hmz");
  hmz->Divide(hmz_mc);

  // \BEGIN copy-paste from minitools/drawZmass.C
  TF1 *f1mz = new TF1("f1mz","[0]+[1]*log(0.01*x)+[2]*pow(log(0.01*x),2)",
                      fzptmin, fzptmax);
  //TF1 *f1ez = new TF1("f1ez","sqrt([0]+pow(log(0.01*x),2)*[1]"
  //                  "+pow(log(0.01*x),4)*[2]"
  //                  "+2*log(0.01*x)*[3]+2*pow(log(0.01*x),2)*[4]"
  //                  "+2*pow(log(0.01*x),3)*[5])",
  //                  fzptmin, fzptmax);
  if (correctZMass) {
      // Run2Test fit with minitools/drawZmass.C
    f1mz->SetParameters(0.99875, 0.00118, 0.00059);
    //f1ez->SetParameters(+2.51e-09, +8.73e-09, +3.99e-09,
    //                  +2.15e-09, +7.07e-11, +4.95e-09);
  }

  // \END copy-paste from minitools/drawZmass.C

  // Run 2 2D pT-eta distribution down to 5 GeV based on P8 dijet sample
  TFile *feta = new TFile("rootfiles/P8_dijet_45M_TH2D_correct_weighting.root",
			  "READ");
  assert(feta && !feta->IsZombie());
  TH2D *h2pteta = (TH2D*)feta->Get("pTjet_etajet"); assert(h2pteta);

  // Store pointers to all files in a map for easy handling later
  map<string, TFile*> files;
  files["zjet"] = fz;

  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Results from Sami's Z+b analysis
  //if (epoch=="RunE" || epoch=="RunCD") {
  if (epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
      epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4"
      ) {
    rename["zjet"]["ratio"] = "data"; // missing => PATCH
    rename["zjet"]["data"] = "data";
    rename["zjet"]["mc"] = "mc";
    rename["zjet"]["mpfchs"] = "rmpf";
    rename["zjet"]["mpfchs1"] = "rmpf";
    rename["zjet"]["ptchs"] = "rmpfjet1";
    rename["zjet"]["counts"] = "statistics_rmpf";
    rename["zjet"]["chf"] = "chHEF";
    rename["zjet"]["nef"] = "neEmEF";
    rename["zjet"]["nhf"] = "neHEF";
    rename["zjet"]["cef"] = "chEmEF";
    rename["zjet"]["muf"] = "muEF";
    //
    rename["zjet"]["mpf1"] = "rmpfjet1";
    rename["zjet"]["mpfn"] = "rmpfjetn";
    rename["zjet"]["mpfu"] = "rmpfuncl";
    rename["zjet"]["rho"] = "rho";
    
    rename["zjet"]["rjet"] = "rbal";
    rename["zjet"]["gjet"] = "rgenjet1";
  }
  /*
  if (epoch=="RunF") { // v2
    rename["zjet"]["ratio"] = "2022F/NoCorrection/";
    rename["zjet"]["data"] = "2022F/NoCorrection/";
    rename["zjet"]["mc"] = "MC/NoCorrection/";
    rename["zjet"]["mpfchs"] = "RMPF";
    rename["zjet"]["mpfchs1"] = "RMPF";
    rename["zjet"]["ptchs"] = "RMPFjet1";
    rename["zjet"]["counts"] = "statistics_rmpf";
    rename["zjet"]["chf"] = "chHEF";
    rename["zjet"]["nef"] = "neEmEF";
    rename["zjet"]["nhf"] = "neHEF";
    rename["zjet"]["cef"] = "chEmEF";
    rename["zjet"]["muf"] = "muEF";
    //
    rename["zjet"]["mpf1"] = "RMPFjet1";
    rename["zjet"]["mpfn"] = "RMPFjetn";
    rename["zjet"]["mpfu"] = "RMPFuncl";
  }
  */
  if (epoch=="RunE" || epoch=="RunF") { // v3
    rename["zjet"]["ratio"] = "data"; // patch
    rename["zjet"]["data"] = "data";
    rename["zjet"]["mc"] = "mc";
    rename["zjet"]["mpfchs"] = "mpf";
    rename["zjet"]["mpfchs1"] = "mpf";
    rename["zjet"]["ptchs"] = "db";
    rename["zjet"]["counts"] = "RawNEvents_zpt_a100";
    rename["zjet"]["chf"] = "chf";
    rename["zjet"]["nef"] = "nef";
    rename["zjet"]["nhf"] = "nhf";
    rename["zjet"]["cef"] = "cef";
    rename["zjet"]["muf"] = "muf";
    //
    rename["zjet"]["mpf1"] = "mpf1";
    rename["zjet"]["mpfn"] = "mpfn";
    rename["zjet"]["mpfu"] = "mpfu";
    //
    rename["zjet"]["rho"] = "rho";
    rename["zjet"]["npv"] = "npv";
  }

  // color and style codes
  map<string, map<string, int> > style;

  style["zjet"]["mpfchs1"] = kFullDiamond;
  style["zjet"]["ptchs"] = kOpenDiamond;
  style["zjet"]["chf"] = kFullCircle;
  style["zjet"]["nhf"] = kFullDiamond;
  style["zjet"]["nef"] = kFullSquare;
  style["zjet"]["cef"] = kFullDiamond;
  style["zjet"]["muf"] = kFullDiamond;
  style["zjet_mc"]["chf"] = kOpenCircle;
  style["zjet_mc"]["nhf"] = kOpenDiamond;
  style["zjet_mc"]["nef"] = kOpenSquare;
  style["zjet_mc"]["cef"] = kOpenDiamond;
  style["zjet_mc"]["muf"] = kOpenDiamond;
  style["zjet"]["mpf1"] = kFullTriangleUp;
  style["zjet"]["mpfn"] = kFullTriangleUp;
  style["zjet"]["mpfu"] = kFullTriangleDown;
  style["zjet"]["rho"] = kFullTriangleDown;
  style["zjet"]["npv"] = kFullTriangleDown;
  style["zjet_mc"]["mpf1"] = kOpenTriangleUp;
  style["zjet_mc"]["mpfn"] = kOpenTriangleUp;
  style["zjet_mc"]["mpfu"] = kOpenTriangleDown;
  style["zjet_mc"]["rho"] = kOpenTriangleDown;
  style["zjet_mc"]["npv"] = kOpenTriangleDown;

  map<string, int> color;
  color["zjet"] = kRed+1;
  color["zjet_muf"] = kMagenta+1;
  color["zjet_mpf1"] = kRed;
  color["zjet_mpfn"] = kGreen+2;
  color["zjet_mpfu"] = kBlue;
  color["zjet_rho"] = kBlack;
  color["zjet_chf"] = kRed;
  color["zjet_nhf"] = kGreen+2;
  color["zjet_nef"] = kBlue;
  color["zjet_cef"] = kCyan+1;
  color["zjet_muf"] = kMagenta+1;

  // Select which subsets of data to process:
  // datamc x method x sample x etabin (x alphacut
  //    3   x    3   x   4    x   6-8  (x    4)      = 864-1152 graphs
  vector<string> dirs;
  dirs.push_back("mc");
  dirs.push_back("data");
  dirs.push_back("ratio");

  vector<string> types;
  types.push_back("counts");
  types.push_back("crecoil");
  types.push_back("mpfchs1"); // Type-I MET
  types.push_back("ptchs");
  // for pfjet only (activate puf, cef, muf later?)
  if (epoch=="RunCD" || epoch=="RunE" || epoch=="RunF" ||
      epoch=="Run22C" || epoch=="Run22D" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
      epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" ||
      epoch=="Run3") {
    types.push_back("chf");
    types.push_back("nef");
    types.push_back("nhf");
    types.push_back("cef");
    types.push_back("muf");
  }
  //types.push_back("puf");
  types.push_back("mpf1");
  types.push_back("mpfn");
  types.push_back("mpfu");
  if (epoch=="RunCD" || epoch=="RunE" || epoch=="RunF" ||
      epoch=="Run22C" || epoch=="Run22D" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
      epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" ||
      epoch=="Run3") {
    types.push_back("rho");
  }
  if (epoch=="RunE" || epoch=="RunF") {
    types.push_back("npv");
  }
  // <pT,reco> and <pT,gen> vs ref pT (MC only)
  //if (epoch=="RunE" || epoch=="RunCD") {
  if (epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
      epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" ||
      epoch=="Run3") {
    types.push_back("rjet");
    types.push_back("gjet");
  }

  vector<string> sets;
  sets.push_back("zjet");

  vector<pair<double,double> > etas;
  etas.push_back(make_pair<double,double>(0,1.305));
  //etas.push_back(make_pair<double,double>(0,2.500));
  // Narrow eta bins for L2Res
  
  /*
  etas.push_back(make_pair<double,double>(0.000,0.261)); 
  etas.push_back(make_pair<double,double>(0.261,0.522)); 
  // PATCH!! V15 gamjet combines 0-0.5 bins
  etas.push_back(make_pair<double,double>(0.522,0.783)); 
  etas.push_back(make_pair<double,double>(0.783,1.044)); 
  etas.push_back(make_pair<double,double>(1.044,1.305)); 
  // PATCH!! V15 gamjet has 0.8-1.1? 
  etas.push_back(make_pair<double,double>(1.305,1.479)); 
  etas.push_back(make_pair<double,double>(1.479,1.653)); 
  // PATCH!! V15 gamjet has 1.3-1.7? 
  etas.push_back(make_pair<double,double>(1.653,1.930)); 
  etas.push_back(make_pair<double,double>(1.930,2.172)); 
  etas.push_back(make_pair<double,double>(2.172,2.322)); 
  etas.push_back(make_pair<double,double>(2.322,2.500)); 
  etas.push_back(make_pair<double,double>(2.500,2.650)); 
  etas.push_back(make_pair<double,double>(2.650,2.853)); 
  etas.push_back(make_pair<double,double>(2.853,2.964)); 
  etas.push_back(make_pair<double,double>(2.964,3.139)); 
  etas.push_back(make_pair<double,double>(3.139,3.489)); 
  etas.push_back(make_pair<double,double>(3.489,3.839)); 
  etas.push_back(make_pair<double,double>(3.839,5.191));

  // Wide eta bins for L2L3Res closure
  etas.push_back(make_pair<double,double>(1.305,1.93));
  etas.push_back(make_pair<double,double>(1.93,2.5));
  etas.push_back(make_pair<double,double>(2.5,2.964));
  etas.push_back(make_pair<double,double>(2.964,3.2));
  etas.push_back(make_pair<double,double>(3.2,5.191));
  */


  ///////////////////////////////////////////
  // Rename selected graphs and store them //
  ///////////////////////////////////////////

  map<string, map<string, map<string, map<int, TGraphErrors*> > > > grs;
  map<string, map<string, map<int, TH1D*> > > counts;

  // Loop over data, MC, ratio
  for (unsigned int idir = 0; idir != dirs.size(); ++idir) {

    string d = dirs[idir];
    const char *dd = d.c_str();

    fout->mkdir(dd);
    assert(fout->cd(dd));
    TDirectory *dout0 = fout->GetDirectory(dd); assert(dout0);

    // Loop over eta bins
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
      
      double eta1 = etas[ieta].first;
      double eta2 = etas[ieta].second;
      const char *dd0 = Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.);
      cout << dd0 << endl << flush;
      
      dout0->mkdir(dd0);
      assert(dout0->cd(dd0));
      TDirectory *dout = dout0->GetDirectory(dd0); assert(dout);
      dout->mkdir("orig");
      dout->cd();

      // Save background histogram for easy drawing from root file
      TH1D *h = new TH1D("h",Form(";p_{T} (GeV);Response (%s)",dd),
			 int(4000-10),10,4000);
      h->GetXaxis()->SetMoreLogLabels();
      h->GetXaxis()->SetNoExponent();
      h->SetMinimum(0.50);
      h->SetMaximum(1.50);
      h->Write();
      
      // Loop over methods (MPF, pT balance)
      for (unsigned int itype = 0; itype != types.size(); ++itype) {
	
	string t = types[itype];
	const char* tt = t.c_str();

	// Loop over samples (Z+jet, gamma+jet, multijet, W>qq')
	for (unsigned int iset = 0; iset != sets.size(); ++iset) {
	
	  string s = sets[iset];
	  const char* ss = s.c_str();
	  string sp = (rename[s]["parent"]!=0 ? rename[s]["parent"] : "");

	  TFile *f = files[s];

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

	  bool ismpfc = (t=="mpf1" || t=="mpfn" || t=="mpfu");// || t=="rho");
	  bool isfrac = (t=="chf"||t=="nef"||t=="nhf"||
			 t=="cef"||t=="muf"||t=="puf");
	  bool isflavor = 
	    (s=="zi"  || s=="zb"  || s=="zc"  || s=="zq"  || s=="zg"||s=="zn")||
	    (s=="gi"  || s=="gb"  || s=="gc"  || s=="gq"  || s=="gg"||s=="gn");
	  bool isflavormc = 
	    (s=="zii" || s=="zbi" || s=="zci" || s=="zqi"||s=="zgi"||s=="zni"||
	     s=="zib" || s=="zbb" || s=="zcb" || s=="zqb"||s=="zgb"||s=="znb"||
	     s=="zic" || s=="zbc" || s=="zcc" || s=="zqc"||s=="zgc"||s=="znc"||
	     s=="ziq" || s=="zbq" || s=="zcq" || s=="zqq"||s=="zgq"||s=="znq"||
	     s=="zig" || s=="zbg" || s=="zcg" || s=="zqg"||s=="zgg"||s=="zng"||
	     s=="zin" || s=="zbn" || s=="zcn" || s=="zqn"||s=="zgn"||s=="znn"||
	     s=="ziu" || s=="zis")||
	    (s=="gii" || s=="gbi" || s=="gci" || s=="gqi"||s=="ggi"||s=="gni"||
	     s=="gib" || s=="gbb" || s=="gcb" || s=="gqb"||s=="ggb"||s=="gnb"||
	     s=="gic" || s=="gbc" || s=="gcc" || s=="gqc"||s=="ggc"||s=="gnc"||
	     s=="giq" || s=="gbq" || s=="gcq" || s=="gqq"||s=="ggq"||s=="gnq"||
	     s=="gig" || s=="gbg" || s=="gcg" || s=="gqg"||s=="ggg"||s=="gng"||
	     s=="gin" || s=="gbn" || s=="gcn" || s=="gqn"||s=="ggn"||s=="gnn");
	  if (isflavormc && d!="mc") continue;

	  // true responses only for gamjet for now => add to Z+jet
	  //if ((t=="rjet" || t=="gjet") && sp!="gamjet") continue;
	  if ((t=="rjet" || t=="gjet") && d!="mc") continue;

	  double alpha = 1.00;
	  eta1 = etas[ieta].first; // reset to avoid trouble with below
	  eta2 = etas[ieta].second; // reset to avoid trouble with below
	  
	  //int ieta1 = int(10.*eta1+0.5);
	  //int ieta2 = int(10.*eta2+0.5);
	  
	  // Reconstruct naming scheme used in each of the files
	  // If non-conventional naming schemes, patch here
	  const char *c(0);
	  if (s=="zjet") {
	    //if (epoch=="RunE" || epoch=="RunCD") {
	    if (epoch=="RunCD" ||
		epoch=="Run22C" || epoch=="Run22D" ||
		epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
		epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
		epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
		epoch=="Run23C4"
		) {
	      c = Form("%s/eta_%02.0f_%02.0f/%s%s_zmmjet_a%1.0f",
		       rename[s][d],10*eta1,10*eta2,
		       ((d=="data"||d=="ratio"||t=="counts") ? "" : 
			sPSWgtZ.c_str()),
		       rename[s][t],100.*alpha);
	      if (t=="rho")
		c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f",
			 rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
	      if (isfrac)
		c = Form("%s/eta_%02.0f_%02.0f/h_Zpt_%s_alpha%1.0f",
			 rename[s][d],10*eta1,10*eta2,rename[s][t],100.*alpha);
	    }
	    //if (epoch=="RunF") {
	    //c = Form("%sh_ZpT_%s_alpha100_eta13",
	    //	       rename[s][d],rename[s][t]);
	    //}
	    if (epoch=="RunE" || epoch=="RunF") {
	      c = Form("%s/eta_00_13/zpt_%s_zmmjet_a100",
		       rename[s][d],rename[s][t]);
	    }
	    if (epoch=="Run3") {
	      c = Form("%s/eta00-13/%s_%s_a100",d.c_str(),t.c_str(),s.c_str());
	    }
	  } // "zjet"
    
	  
	  TObject *obj = f->Get(c);
	  if (!obj) {
	    cout << "Graph " << c << " not found for "
		 << s << " " << t << " " << d << "!" <<endl << flush;
	    cout << "File: " << f->GetName() << endl << flush;
	    cout << "Eta " << eta1 << " - " << eta2 <<  endl << flush;
	  }
	  if (t=="counts" && !obj) obj = hcounts;
	  
	  // Calculate data/MC ratio
	  if ((s=="zjet" || sp=="zjet") &&
	      d=="ratio" && (t=="mpfchs1"||t=="ptchs")) {
	    TGraphErrors *gd = grs["data"][t][s][ieta];
	    TGraphErrors *gm = grs["mc"][t][s][ieta];
	    assert(gd);
	    assert(gm);
	    
	    TGraphErrors *g = tools::ratioGraphs(gd, gm);
	    obj = (TObject*)g;
	  }
	  if (s=="zjet" && d=="ratio" &&
	      (isfrac || t=="rmpf1" || t=="rmpfn" || t=="rpmfn")) {
	    TGraphErrors *gd = grs["data"][t][s][ieta];
	    TGraphErrors *gm = grs["mc"][t][s][ieta];
	    assert(gd);
	    assert(gm);
	    
	    TGraphErrors *g = tools::diffGraphs(gd, gm);
	    obj = (TObject*)g;
	  }

	  if (!obj) {
	    cout << "Missing " << c << endl << flush;
	    cout << "s="<<s<<", d="<<d<<", t="<<t<<endl<<flush;
	  }
	  assert(obj);
	  
	  // write out counts to jecdata.root (as TH1F)
	  if (t=="counts") {
	      
	    assert(obj->InheritsFrom("TH1D") ||obj->InheritsFrom("TH1F"));
	    dout->cd();
	    TH1D *h = (TH1D*)obj;
	    h->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alpha));

	    h->Write();
	    counts[d][s][ieta] = h;
	    continue;
	  }

	  // If data stored in TH2D instead of TGrapherrors, patch here
	  // Patch for zjet that has PF fractions in TH2D (v24)
	  if (obj->InheritsFrom("TH2")) {
	    obj = new TGraphErrors(((TH2D*)obj)->ProfileX()->ProjectionX());
	  }

	  // If data stored in TH1D instead of TGraphErrors, patch here
	  // Patch for dijet file that has TH1D's instead of graphs
	  if (obj->InheritsFrom("TH1")) {
	    obj = new TGraphErrors((TH1D*)obj);
	  }
	  
	  // If data stored in TProfile instead of TGraphErrors, patch here
	  // Patch for pfjet file that has TProfiles's instead of graphs
	  // Patch for zjet file that has TProfiles instead of graphs
	  if (obj->InheritsFrom("TProfile")) {
	    obj = new TGraphErrors(((TProfile*)obj)->ProjectionX());
	  }

	  assert(obj->InheritsFrom("TGraphErrors"));
	  TGraphErrors *g = (TGraphErrors*)obj;
	    
	  // Clean out empty points from TH1D->TGraphErrors conversion
	  for (int i = g->GetN()-1; i != -1; --i) {
	    assert(i<=g->GetN()-1);
	    // Clean out spurious empty pooints
	    if (g->GetY()[i]==0 && g->GetEY()[i]==0) g->RemovePoint(i);
	  } // for i

	  // Make a copy of raw data before cleaning cuts and corrections
	  // for documenting it
	  TGraphErrors *g_orig = (TGraphErrors*)g->Clone();
	  g_orig->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alpha));
	    
	  // Select stable range and good statistics points for global fit
	  for (int i = g->GetN()-1; i != -1; --i) {
	    assert(i<=g->GetN()-1);
	    // remove points with large error (more than 0.2 right now)
	    if (g->GetEY()[i]>0.2 && !isflavormc)  g->RemovePoint(i);
	    // Clean out point outside good ranges
	    else if (s=="zjet" &&
		     (g->GetX()[i]<fzptmin || g->GetX()[i]>fzptmax))
	      g->RemovePoint(i);
	  } // for i
	  
	  // Mass corrections for Z+jet
	  if (correctZMass && s=="zjet" && (d=="data" || d=="ratio") &&
	      //(t=="mpfchs1" || t=="mpf1" || t=="ptchs")) {
	      (ismpfc || t=="ptchs"))  {
	    for (int i = 0; i != g->GetN(); ++i) {
	      double pt = g->GetX()[i];
	      //double ek = f1ez->Eval(pt);
	      double k = f1mz->Eval(pt);
	      g->SetPoint(i, g->GetX()[i], g->GetY()[i]*k);
	      //if (correctUncert)
	      //g->SetPointError(i, g->GetEX()[i], 
	      //		 sqrt(pow(g->GetEY()[i]*k,2) + ek*ek));
	    }
	  }

	  // Patch mpf1, mpfn and mpfu
	  if (epoch=="RunE" || epoch=="RunF") {

	    if (s=="zjet" && (d=="data" || d=="mc") && t=="mpf1") {
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]+0.10);
	      }
	    } // patch mpf1
	    if (s=="zjet" && (d=="data" || d=="mc") && t=="mpfn") {
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]+1-0.20);
	      }
	    } // patch mpfu
	    if (s=="zjet" && (d=="data" || d=="mc") && t=="mpfu") {
	      for (int i = 0; i != g->GetN(); ++i) {
		g->SetPoint(i, g->GetX()[i], g->GetY()[i]-1+0.10);
	      }
	    } // patch mpfn
	  } // RunE, RunF patch

	  dout->cd("orig");
	  g_orig->Write();

	  dout->cd();
	  
	  // Set uniform naming scheme and graphical style
	  g->SetName(Form("%s_%s_a%1.0f",tt,ss,100.*alpha));
	  g->UseCurrentStyle(); // Basic TDR style
	  g->SetMarkerStyle(style[s][t]);
	  g->SetMarkerColor(color[s]);
	  g->SetLineColor(color[s]);
	  if (isfrac || ismpfc) {
	    g->SetMarkerStyle((d=="mc" && style[s+"_"+d][t]!=0) ?
			      style[s+"_"+d][t] : style[s][t]);
	    g->SetMarkerColor(color[s+"_"+t]!=0 ? color[s+"_"+t] : color[s]);
	    g->SetLineColor(color[s+"_"+t]!=0 ? color[s+"_"+t] : color[s]);
	  }
	  g->SetDrawOption("SAMEP");

	  g->Write();

	  grs[d][t][s][ieta] = g;
	} // for iset
      } // for itype
    } // for ieta
  } // for itier

  // Save mass histos as well
  cout << "Save mass histos" << endl << flush;
  fout->cd("ratio/eta00-13");
  hmz->Write("mass_zjet_a100");

  fout->cd("data/eta00-13");
  hmz_dt->Write("mass_zjet_a100");

  fout->cd("mc/eta00-13");
  hmz_mc->Write("mass_zjet_a100");

  if (fz) fz->Close();
  curdir->cd();


  /////////////////////////////////////////////////////
  // Calculate JEC central values and uncertainties, //
  // and store them for easy visualialization later  //
  /////////////////////////////////////////////////////

  if (rp_debug) cout << "Instantiating JECs..." << endl << flush;

  map<string,const char*> mera;
  mera["2018ABCD"] = "ABCD";

  // Calculate L2L3Res with JEC uncertainty
  {

    const char *s, *s2;
    // Usual directory for text files
    const char *cd = "CondFormats/JetMETObjects/data";
    
    // New JEC for plotting on the back and mcjec for softrad3.C
    // ** Also used as reference for CorLevel=="L1L2L3Res" **
    // So be careful when running minitools/createL2L3Res.C
    FactorizedJetCorrector *mcjec(0), *jec(0);
    //jec = getFJC("","",Form("Winter22Run3_Run%s_V1_DATA_L2L3Residual","A"));
    //mcjec = getFJC("",Form("Winter22Run3_Run%s_V1_DATA_L2Relative","A"));
    if (epoch=="Run22C" || epoch=="Run22D") {
      jec = getFJC("","","Winter22Run3_RunC_V2_DATA_L2L3Residual");
      mcjec = getFJC("","Winter22Run3_RunC_V2_DATA_L2Relative");
    }
    if (epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G") {
      jec = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
      mcjec = getFJC("","Winter23Prompt23_RunA_V1_DATA_L2Relative");
    }
    if (epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
	epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
	epoch=="Run23C4" ||
	epoch=="Run3") {
      jec = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
      mcjec = getFJC("","Winter23Prompt23_RunA_V1_DATA_L2Relative");
    }
    assert(jec);
    assert(mcjec);

    //#define PAIR(a,b) (make_pair<double,FactorizedJetCorrector*>((a),getFJC("","",(b))))
    //add vjec here later to combine IOVs

    if (rp_debug) cout << "Loading reference JECs..." << endl << flush;
    
    FactorizedJetCorrector *jecrun1, *jecrun2, *jecold;

    // Reference Run I and Run II JEC for plotting on the back
    jecrun1 = getFJC("","","Winter14_V8_DATA_L2L3Residual_AK5PFchs"); 
    jecrun2 = getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs"); 
    
    // Store old JEC for undoing it in global fit (JEC from closure files)
    //jecold = getFJC("","",Form("Winter22Run3_Run%s_V1_DATA_L2L3Residual","A"));
    if (epoch=="Run22C" || epoch=="Run22D") {
      jecold = getFJC("","","Winter22Run3_RunC_V2_DATA_L2L3Residual");
    }
    if (epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G") {
      jecold = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
    }
    if (epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C" ||
	epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
	epoch=="Run23C4" ||
	epoch=="Run3") {
      jecold = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
    }
    assert(jecold);

    if (rp_debug) cout << "Loading uncertainty sources..." << endl << flush;

    // Run I uncertainty
    s = Form("%s/Winter14_V8_DATA_UncertaintySources_AK5PFchs.txt",cd);
    s2 = "SubTotalAbsolute";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref1 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref1 = new JetCorrectionUncertainty(*p_ref1);

    // Run II uncertainty
    s = Form("%s/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.txt",cd);
    s2 = "SubTotalAbsolute";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref2 = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref2 = new JetCorrectionUncertainty(*p_ref1);

    // Total uncertainty, excluding Flavor and Time
    s = Form("%s/Winter22Run3_RunA_V1_DATA_UncertaintySources_AK4PFPuppi.txt",cd);
    s2 = "TotalNoFlavorNoTime";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_unc = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p_unc);
    //
    s2 = "SubTotalAbsolute";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ref = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ref = new JetCorrectionUncertainty(*p_ref);
    //
    s2 = "SubTotalPt";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_pt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_pt = new JetCorrectionUncertainty(*p_pt);
    //
    s2 = "SinglePionHCAL";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_hcal = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_hcal = new JetCorrectionUncertainty(*p_hcal);
    //
    s2 = "SinglePionECAL";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_ecal = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_ecal = new JetCorrectionUncertainty(*p_ecal);
    //
    s2 = "FlavorPureGluon";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_glu = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_glu = new JetCorrectionUncertainty(*p_glu);
    //
    s2 = "FlavorPureQuark";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_uds = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_uds = new JetCorrectionUncertainty(*p_uds);
    //
    s2 = "FlavorPureCharm";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_cha = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_cha = new JetCorrectionUncertainty(*p_cha);
    //
    s2 = "FlavorPureBottom";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_bot = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_bot = new JetCorrectionUncertainty(*p_bot);
    //
    s2 = "FlavorZJet";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_zjt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_zjt = new JetCorrectionUncertainty(*p_zjt);
    //
    s2 = "FlavorPhotonJet";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_gjt = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_gjt = new JetCorrectionUncertainty(*p_gjt);
    //
    s2 = "FlavorQCD";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_qcd = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_qcd = new JetCorrectionUncertainty(*p_qcd);
    //
    s2 = "SubTotalPileUp";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_pu = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_pu = new JetCorrectionUncertainty(*p_pu);
    //
    s2 = "TotalNoFlavor";
    cout << s << ":" << s2 << endl << flush;
    JetCorrectorParameters *p_noflv = new JetCorrectorParameters(s,s2);
    JetCorrectionUncertainty *unc_noflv = new JetCorrectionUncertainty(*p_noflv);

    // Loop over eta bins, but do JEC for data/MC ratio only
    for (unsigned int ieta = 0; ieta != etas.size(); ++ieta) {
      
      assert(fout->cd("ratio"));
      TDirectory *dout1 = fout->GetDirectory("ratio"); assert(dout1);
      double eta1 = etas[ieta].first; double eta2 = etas[ieta].second;
      const char *dd1 = Form("eta%02.0f-%02.0f",eta1*10.,eta2*10.);
      assert(dout1->cd(dd1));
      TDirectory *dout2 = dout1->GetDirectory(dd1); assert(dout2);
      dout2->cd();
      cout << "Eta bin:" << eta1 <<"-"<< eta2 << endl;

      const double ptbins[] = {15, 16, 18, 20, 22, 25,
			       30, 35, 40, 50, 60, 70, 85, 100, 125, 155, 180,
			       210, 250, 300, 350, 400, 500, 600, 800, 1000,
			       1200, 1500,
			       1800, 2100, 2400, 2700, 3000, 3300, 3600,
			       3900, 4200, 4500, 4501};
      const int npt = sizeof(ptbins)/sizeof(ptbins[0])-1;

      // Uncertainty bands
      TH1D *herr = new TH1D("herr",";p_{T} (GeV);JEC uncertainty;",
			    npt, &ptbins[0]);
      TH1D *herr_l2l3res = new TH1D("herr_l2l3res",";p_{T} (GeV);L2L3Res;",
				npt, &ptbins[0]);
      TH1D *herr_ref = new TH1D("herr_ref",";p_{T} (GeV);TotalNoFlavorNoTime;",
				npt, &ptbins[0]);
      TH1D *herr_spr = new TH1D("herr_spr",";p_{T} (GeV);SPR uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_glu = new TH1D("herr_glu",";p_{T} (GeV);Gluon uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_uds = new TH1D("herr_uds",";p_{T} (GeV);Quark uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_cha = new TH1D("herr_cha",";p_{T} (GeV);Charm uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_bot = new TH1D("herr_bot",";p_{T} (GeV);Bottom uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_zjt = new TH1D("herr_zjt",";p_{T} (GeV);Z+jet uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_gjt = new TH1D("herr_gjt",";p_{T} (GeV);#gamma+jet uncert.;",
				npt, &ptbins[0]);
      TH1D *herr_qcd = new TH1D("herr_qcd",";p_{T} (GeV);QCD uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_pu = new TH1D("herr_pu",";p_{T} (GeV);PU uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_noflv = new TH1D("herr_noflv",";p_{T} (GeV);TotalNoFlavor;",
				npt, &ptbins[0]);
      TH1D *herr_mpf = new TH1D("herr_mpf",";p_{T} (GeV);MPF uncertainty;",
				npt, &ptbins[0]);
      TH1D *herr_pt = new TH1D("herr_pt",";p_{T} (GeV);p_{T} uncertainty;",
			       npt, &ptbins[0]);

      // Run I JEC central value
      TH1D *hrun1 = new TH1D("hrun1",";p_{T} (GeV);Run I JES;",
			     npt, &ptbins[0]);

      // Run II JEC central value
      TH1D *hrun2 = new TH1D("hrun2",";p_{T} (GeV);Run II JES;",
			     npt, &ptbins[0]);

      // JEC central value
      TH1D *hjes = new TH1D("hjes",";p_{T} (GeV);Old JES;",npt, &ptbins[0]);


      if(rp_debug) cout << "calculate reference R_uncl (hruncl)" << endl<<flush;

      TH1D *hruncl = new TH1D("hruncl",";mode;Runcl",4,0.5,4.5);
      if (CorLevel=="L1L2L3Res") { // R_uncl
	double ptref = 15.;
	int jref = h2pteta->GetXaxis()->FindBin(ptref);
	jref = max(1,min(jref, h2pteta->GetXaxis()->GetNbins()));
	TH1D *hpteta = h2pteta->ProjectionY("hpteta",jref,jref);

	double suml2l3(0), suml2l3res(0), sumw5(0), sumw8(0);
	for (int ieta = 1; ieta != hpteta->GetNbinsX()+1; ++ieta) {

	  double abseta = hpteta->GetBinCenter(ieta);
	  double w = hpteta->GetBinContent(ieta);
	  sumw8 += w;
	  if (abseta<5.191) {
	    sumw5 += w;
	    double l2l3p = 1./getJEC(mcjec,+abseta,ptref);
	    double l2l3m = 1./getJEC(mcjec,-abseta,ptref);
	    double l2l3 = 0.5*(l2l3p+l2l3m);
	    suml2l3 += w*l2l3;
	    double l2l3resp = 1./getJEC(jec,+abseta,ptref);
	    double l2l3resm = 1./getJEC(jec,-abseta,ptref);
	    double l2l3res = 0.5*(l2l3resp+l2l3resm);
	    suml2l3res += w*l2l3res;
	  }
	} // for ieta
	hruncl->SetBinContent(1, sumw5/sumw8);
	hruncl->GetXaxis()->SetBinLabel(1,"|#eta|<5.2");
	hruncl->SetBinContent(2, suml2l3/sumw5);
	hruncl->GetXaxis()->SetBinLabel(2,"L2L3");
	hruncl->SetBinContent(3, suml2l3res/sumw5);
	hruncl->GetXaxis()->SetBinLabel(3,"L2L3Res");
	hruncl->SetBinContent(4, (suml2l3/sumw5)*(suml2l3res/sumw5));
	hruncl->GetXaxis()->SetBinLabel(4,"DATA");
      } // R_uncl

      if(rp_debug) cout << "create reference JES bands" << endl << flush;

      for (int ipt = 1; ipt != herr->GetNbinsX()+1; ++ipt) {

	double pt = herr->GetBinCenter(ipt);

	int jpt = h2pteta->GetXaxis()->FindBin(pt);
	jpt = max(1,min(jpt, h2pteta->GetXaxis()->GetNbins()));

	// Take Z+jet eta,pT distribution for correctly averaging JEC
	TH1D *heta = h2pteta->ProjectionY(Form("heta_%d",ipt), jpt, jpt);
	const int ieta2 = heta->FindBin(eta2);
	const int ieta1 = heta->FindBin(eta1);
	const int intw = heta->Integral(ieta1,ieta2-1);
	const int neta = 2*(ieta2 - ieta1);

	// Loop over fine eta bins to average JEC and uncertainties
	double sumval(0), sumvall2l3res(0), sumerr2(0), sumw(0);
	double sumrun1(0), sumrun2(0), sumjes(0);
	//
	double sumerr2_pt(0), sumerr2_hcal(0), sumerr2_ecal(0), sumerr2_pu(0);
	double sumerr2_glu(0), sumerr2_uds(0), sumerr2_cha(0), sumerr2_bot(0);
	double sumerr2_zjt(0), sumerr2_gjt(0), sumerr2_qcd(0);
	double sumerr2_noflv(0), sumerr2_ref(0);
	double sumerr2_ref1(0), sumerr2_ref2(0);

	for (int jeta = 0; jeta != neta; ++jeta) {

	  assert(eta2 > eta1);
	  assert(eta1 >= 0);

	  // average over both plus and minus sides
	  int keta = ieta1 + (jeta % (neta/2));
	  double eta = (jeta<neta/2 ? -1 : +1) * heta->GetBinCenter(keta);
	  double w = (intw ? heta->GetBinContent(keta) / (2*intw) : 1./neta);
	  if (!(w<=1)) {
	    cout << "pt =" << pt << " w="<<w << endl << flush;
	    assert(w<=1);
	  }

	  // JEC central values first
	  double val = 1./getJEC(jec,eta,pt);

	  // Run I and Run II reference JECs
	  double jesrun1 = 1./getJEC(jecrun1,eta,pt);
	  double jesrun2 = 1./getJEC(jecrun2,eta,pt);

	  // old JEC
	  double jes = 1./getJEC(jecold,eta,pt);

	  double vall2l3res = val; // For closure test 
	  if(CorLevel=="L1L2L3Res") { //to get proper bands during closure test
	    val /= jes; 
	    jesrun1 /= jes;
	    jesrun2 /= jes;
	    jes /= jes;
	  }

	  sumvall2l3res += w*vall2l3res;
	  sumrun1 += w*jesrun1;
	  sumrun2 += w*jesrun2;
	  sumval += w*val;
	  sumjes += w*jes;
	  sumw += w; // sum weights only once
	  
	  // JEC uncertainties
	  unc->setJetEta(eta);
	  unc->setJetPt(pt);
	  double err = unc->getUncertainty(true);
	  sumerr2 += w*err*err; // sum squared weights only once

	  unc_ref->setJetEta(eta);
	  unc_ref->setJetPt(pt);
	  double err_ref = unc_ref->getUncertainty(true);
	  sumerr2_ref += w*err_ref*err_ref;

	  unc_ref1->setJetEta(eta);
	  unc_ref1->setJetPt(pt);
	  double err_ref1 = unc_ref1->getUncertainty(true);
	  sumerr2_ref1 += w*err_ref1*err_ref1;

	  unc_ref2->setJetEta(eta);
	  unc_ref2->setJetPt(pt);
	  double err_ref2 = unc_ref2->getUncertainty(true);
	  sumerr2_ref2 += w*err_ref2*err_ref2;

	  unc_pt->setJetEta(eta);
	  unc_pt->setJetPt(pt);
	  double err_pt = unc_pt->getUncertainty(true);
	  sumerr2_pt += w*err_pt*err_pt;

	  unc_hcal->setJetEta(eta);
	  unc_hcal->setJetPt(pt);
	  double err_hcal = unc_hcal->getUncertainty(true);
	  sumerr2_hcal += w*err_hcal*err_hcal;

	  unc_ecal->setJetEta(eta);
	  unc_ecal->setJetPt(pt);
	  double err_ecal = unc_ecal->getUncertainty(true);
	  sumerr2_ecal += w*err_ecal*err_ecal;

	  unc_glu->setJetEta(eta);
	  unc_glu->setJetPt(pt);
	  double err_glu = unc_glu->getUncertainty(true);
	  sumerr2_glu += w*err_glu*err_glu;

	  unc_uds->setJetEta(eta);
	  unc_uds->setJetPt(pt);
	  double err_uds = unc_uds->getUncertainty(true);
	  sumerr2_uds += w*err_uds*err_uds;

	  unc_cha->setJetEta(eta);
	  unc_cha->setJetPt(pt);
	  double err_cha = unc_cha->getUncertainty(true);
	  sumerr2_cha += w*err_cha*err_cha;

	  unc_bot->setJetEta(eta);
	  unc_bot->setJetPt(pt);
	  double err_bot = unc_bot->getUncertainty(true);
	  sumerr2_bot += w*err_bot*err_bot;

	  unc_zjt->setJetEta(eta);
	  unc_zjt->setJetPt(pt);
	  double err_zjt = unc_zjt->getUncertainty(true);
	  sumerr2_zjt += w*err_zjt*err_zjt;

	  unc_gjt->setJetEta(eta);
	  unc_gjt->setJetPt(pt);
	  double err_gjt = unc_gjt->getUncertainty(true);
	  sumerr2_gjt += w*err_gjt*err_gjt;

	  unc_qcd->setJetEta(eta);
	  unc_qcd->setJetPt(pt);
	  double err_qcd = unc_qcd->getUncertainty(true);
	  sumerr2_qcd += w*err_qcd*err_qcd;

	  unc_pu->setJetEta(eta);
	  unc_pu->setJetPt(pt);
	  double err_pu = unc_pu->getUncertainty(true);
	  sumerr2_pu += w*err_pu*err_pu;

	  unc_noflv->setJetEta(eta);
	  unc_noflv->setJetPt(pt);
	  double err_noflv = unc_noflv->getUncertainty(true);
	  sumerr2_noflv += w*err_noflv*err_noflv;
	} // for jeta

	// normalize by total weight for correct average
	double vall2l3res = sumvall2l3res / sumw;
	double val = sumval / sumw;

	// normalize uncertainties (quadratic instead of linear addition)
	double err = sqrt(sumerr2 / sumw);
	double err_ref = sqrt(sumerr2_ref / sumw);
	double err_ref1 = sqrt(sumerr2_ref1 / sumw);
	double err_ref2 = sqrt(sumerr2_ref1 / sumw);
	double err_pt = sqrt(sumerr2_pt / sumw);
	double err_hcal = sqrt(sumerr2_hcal / sumw);
	double err_ecal = sqrt(sumerr2_ecal / sumw);
	double err_glu = sqrt(sumerr2_glu / sumw);
	double err_uds = sqrt(sumerr2_uds / sumw);
	double err_cha = sqrt(sumerr2_cha / sumw);
	double err_bot = sqrt(sumerr2_bot / sumw);
	double err_zjt = sqrt(sumerr2_zjt / sumw);
	double err_gjt = sqrt(sumerr2_gjt / sumw);
	double err_qcd = sqrt(sumerr2_qcd / sumw);
	double err_pu = sqrt(sumerr2_pu / sumw);
	double err_noflv = sqrt(sumerr2_noflv / sumw);

	// center uncertainties around JEC central value
	herr->SetBinContent(ipt, val);
	herr_l2l3res->SetBinContent(ipt, vall2l3res);
	herr_ref->SetBinContent(ipt, val);
	herr_spr->SetBinContent(ipt, val);
	herr_pu->SetBinContent(ipt, 1);
	herr_noflv->SetBinContent(ipt, val);
	herr_mpf->SetBinContent(ipt, val);
	herr_pt->SetBinContent(ipt, val);

	herr->SetBinError(ipt, val*err);
	herr_l2l3res->SetBinError(ipt, vall2l3res*err_ref);
	herr_ref->SetBinError(ipt, val*err_ref);
	herr_spr->SetBinError(ipt,val*sqrt(err_hcal*err_hcal
					   + err_ecal*err_ecal));
	herr_glu->SetBinError(ipt, val*err_glu);
	herr_uds->SetBinError(ipt, val*err_uds);
	herr_cha->SetBinError(ipt, val*err_cha);
	herr_bot->SetBinError(ipt, val*err_bot);
	herr_zjt->SetBinError(ipt, val*err_zjt);
	herr_gjt->SetBinError(ipt, val*err_gjt);
	herr_qcd->SetBinError(ipt, val*err_qcd);
	herr_pu->SetBinError(ipt, val*err_pu);
	herr_noflv->SetBinError(ipt, val*err_noflv);
	herr_mpf->SetBinError(ipt, val*err_pt);
	herr_pt->SetBinError(ipt, val*sqrt(err_pt*err_pt+err_pu*err_pu));

	double run1 = (sumrun1 / sumw);
	hrun1->SetBinContent(ipt, run1);
	hrun1->SetBinError(ipt, run1*err_ref1);

	double run2 = (sumrun2 / sumw);
	hrun2->SetBinContent(ipt, run2);
	hrun2->SetBinError(ipt, run2*err_ref2);

	double jes = (sumjes / sumw);
	hjes->SetBinContent(ipt, jes);
      } // ipt

      if(rp_debug) cout << "done creating reference JES bands" << endl;

      dout2->cd();

      herr->SetMarkerSize(0);
      herr->SetFillStyle(1001);
      herr->SetFillColorAlpha(kYellow+1,0.5);
      herr->Write();

      herr_l2l3res->SetMarkerSize(0);
      herr_l2l3res->SetFillStyle(1001);
      herr_l2l3res->SetFillColorAlpha(kYellow+1,0.5);
      herr_l2l3res->Write();

      herr_ref->SetMarkerSize(0);
      herr_ref->SetFillStyle(1001);
      herr_ref->SetFillColorAlpha(kYellow+1,0.5);
      herr_ref->Write();

      herr_spr->SetMarkerSize(0);
      herr_spr->SetFillStyle(1001);
      herr_spr->SetFillColorAlpha(kYellow,0.5);
      herr_spr->Write();

      herr_glu->SetMarkerSize(0);
      herr_glu->SetFillStyle(1001);
      herr_glu->SetFillColorAlpha(kYellow+1,0.5);
      herr_glu->Write();

      herr_uds->SetMarkerSize(0);
      herr_uds->SetFillStyle(1001);
      herr_uds->SetFillColorAlpha(kYellow+1,0.5);
      herr_uds->Write();

      herr_cha->SetMarkerSize(0);
      herr_cha->SetFillStyle(1001);
      herr_cha->SetFillColorAlpha(kYellow+1,0.5);
      herr_cha->Write();

      herr_bot->SetMarkerSize(0);
      herr_bot->SetFillStyle(1001);
      herr_bot->SetFillColorAlpha(kYellow+1,0.5);
      herr_bot->Write();

      herr_zjt->SetMarkerSize(0);
      herr_zjt->SetFillStyle(1001);
      herr_zjt->SetFillColorAlpha(kYellow+1,0.5);
      herr_zjt->Write();

      herr_gjt->SetMarkerSize(0);
      herr_gjt->SetFillStyle(1001);
      herr_gjt->SetFillColorAlpha(kYellow+1,0.5);
      herr_gjt->Write();

      herr_qcd->SetMarkerSize(0);
      herr_qcd->SetFillStyle(1001);
      herr_qcd->SetFillColorAlpha(kYellow+1,0.5);
      herr_qcd->Write();

      herr_pu->SetMarkerSize(0);
      herr_pu->SetFillStyle(1001);
      herr_pu->SetFillColorAlpha(kYellow+1,0.5);
      herr_pu->Write();

      herr_noflv->SetMarkerSize(0);
      herr_noflv->SetFillStyle(1001);
      herr_noflv->SetFillColorAlpha(kYellow+1,0.5);
      herr_noflv->Write();

      herr_mpf->SetMarkerSize(0);
      herr_mpf->SetFillStyle(1001);
      herr_mpf->SetFillColorAlpha(kYellow+1,0.5);
      herr_mpf->Write();

      herr_pt->SetMarkerSize(0);
      herr_pt->SetFillStyle(1001);
      herr_pt->SetFillColorAlpha(kYellow+1,0.5);
      herr_pt->Write();

      hrun1->SetMarkerSize(0);
      hrun1->SetFillStyle(1001);
      hrun1->SetFillColorAlpha(kCyan+1,0.5);
      hrun1->Write();

      hrun2->SetMarkerSize(0);
      hrun2->SetFillStyle(1001);
      hrun2->SetFillColorAlpha(kMagenta+1,0.5);
      hrun2->Write();

      hjes->Write();

      hruncl->Write();
    } // ieta
  } // JEC+sys

  fout->Close();

} // reprocess

// Helper function to retrieve FactorizedJetCorrector
FactorizedJetCorrector *getFJC(string l1, string l2, string res, string path) {

  // Set default jet algo
  if (l1!="" && !(TString(l1.c_str()).Contains("_AK")))
    l1 += "_AK4PFPuppi";
  if (l2!="" && !(TString(l2.c_str()).Contains("_AK")))
    l2 += "_AK4PFPuppi";
  if (res!="" && !(TString(res.c_str()).Contains("_AK")))
    res += "_AK4PFPuppi";

  if (path=="") path = "CondFormats/JetMETObjects/data";
  const char *cd = path.c_str();
  const char *cl1 = l1.c_str();
  const char *cl2 = l2.c_str();
  const char *cres = res.c_str();
  string s("");

  vector<JetCorrectorParameters> v;
  if (l1!=""){
    s = Form("%s/%s.txt",cd,cl1);
    cout << s << endl << flush;
    JetCorrectorParameters *pl1 = new JetCorrectorParameters(s);
    v.push_back(*pl1);
  }
  if (l2!="") {
    s = Form("%s/%s.txt",cd,cl2);
    cout << s << endl << flush;
    JetCorrectorParameters *pl2 = new JetCorrectorParameters(s);
    v.push_back(*pl2);
  }
  if (res!="") {
    s = Form("%s/%s.txt",cd,cres);
    cout << s << endl << flush;
    JetCorrectorParameters *pres = new JetCorrectorParameters(s);
    v.push_back(*pres);
  }
  FactorizedJetCorrector *jec = new FactorizedJetCorrector(v);

  return jec;
} // getJFC

// Invert JEC to get differences as a function of pTcorr
// Helper functions to find JEC for corrected pt
void setEtaPtRho(FactorizedJetCorrector *jec, double eta, double pt,
                 double rho){

  assert(jec);
  jec->setJetEta(eta);
  jec->setJetPt(pt);
  jec->setRho(rho);
  jec->setJetA(0.50265);

  return;
}

Double_t funcCorrPt(Double_t *x, Double_t *p) {
  
  double eta = p[0];
  double pt = x[0];
  double rho = p[1];
  setEtaPtRho(_thejec, eta, pt, rho);

  return (_thejec->getCorrection() * pt);
}

double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho) {

  // iterate to solve ptreco for given ptcorr
  _thejec = jec;
  if (!fCorrPt) fCorrPt = new TF1("fCorrPt",funcCorrPt,5,6500,2);
  fCorrPt->SetParameters(eta, rho);
  // Find ptreco that gives pTreco*JEC = pTcorr
  double ptreco = fCorrPt->GetX(ptcorr,5,6500);

  setEtaPtRho(jec, eta, ptreco, rho);

  return (jec->getCorrection());
} // getEtaPtE
