// Purpose: Produce jet veto maps from JMENANO analyzer on dijets
//          Updated version for 2024. Sum absolute pulls to avoid cancellation
// Author:  Mikko Voutilainen (at) cern (dot) ch
// Date:    2024-06-14
// NB: Update of minitools/doJetVetoV2.C to match L2Res.C, JERSF.C, L3Res.C
#include "TFile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TMath.h"

#include <string>
#include <vector>

#include "tdrstyle_mod22.C"

// Threshold for veto maps
double pullThreshold = 70; // default: 70

// Separate threshold for outer HF
double pullThresholdHF45 = 70; // default: 70

// Minimum of non-empty towers to consider eta strip
int nMinTowers = 70; // default: 70 (out of 72)
int nMinTowersJES = 70; // default: 70 (out of 72)

// Minimum events to consider bin for Pull and JES average calculation
int nMinPull = 4; // default: 4
int nMinJES = 5; // default: 5

// Use pulls instead of relative changes
bool doPull = true; // default: true

// plot results from each step
bool plotJetVeto = true; // default: true
bool plotJetVetoOnJES = true;
bool plotJetVetoOnJESnorm = false;

// If this is non-zero (not ""), will use just this one trigger
string oneTrig = ""; // default: ""
// Possible options to comment out
//string oneTrig = "HLT_PFJet500";
//string oneTrig = "HLT_PFJet450";
//string oneTrig = "HLT_PFJet320";
//string oneTrig = "HLT_PFJet40";
//string oneTrig = "HLT_PFJet60";
//string oneTrig = "HLT_PFJet80";
//string oneTrig = "HLT_PFJet140";
//string oneTrig = "HLT_PFJet260";

// If this is non-zero (not ""), will use just this one histogram
string oneHist = ""; // default: ""
// Possible options to comment out
//string oneHist = "p2asymm";
//string oneHist = "h2pt";
//string oneHist = "p2chf";
//string oneHist = "p2nhf";

// Helper script to clean out unreliable ranges
void cleanHist(TH2D *h2, string hist, string trg, string run);

void JetVetos(string run, string version);
void JetVeto(string run = "", string version = "vx") {

  if (run!="") { JetVeto(run,version); exit(0); }
  /*
  JetVetos("2022CD",version);
  JetVetos("2022EFG",version);
  JetVetos("2023BC",version);
  JetVetos("2023D",version);
  */
  
  //JetVetos("2024BCD",version);
  //JetVetos("2024E",version);

  //JetVetos("2024BCDE",version);
  //JetVetos("2024FG",version);
  
  //JetVetos("2024BCDEFG",version);
  //JetVetos("2024HI",version);

  //JetVetos("2024BCDEFGHI",version);
  //JetVetos("2024","V9M");

  JetVetos("2025C","V3M");
  JetVetos("2025D","V3M");
  JetVetos("2025E","V3M");
  JetVetos("2025F","V3M");
  JetVetos("2025CDEF","V3M");
}

void JetVetos(string run, string version) {

  #include "Config.C"
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  string sr = (run+"_"+version);
  const char *cr = sr.c_str();

  TFile *fout = new TFile(Form("rootfiles/jetveto%s.root",cr),
			  "RECREATE");
  fout->mkdir("trigs");

  gROOT->ProcessLine(".! touch pdf/JetVeto");
  
  TFile *f(0), *fg(0), *fpark(0);
  if (run=="2022CD") { // Summer22
    lumi_136TeV = "Run2022CD, 8.1 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230816/jmenano_data_out_2022CD_v1.root","READ");
  }
  if (run=="2022EFG") { // Summer22BPix
    lumi_136TeV = "Run2022EFG, 27.0 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230816/jmenano_data_out_2022EF_v1.root","READ");
  }  
  if (run=="2023BC") { // Summer23
    lumi_136TeV = "Run2023BC, 13.6 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230814/nano_data_out_2023BC_v1.root","READ");
  }
  if (run=="2023D") { // Summer23BPix
    lumi_136TeV = "Run2023D, 9.5 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230814/nano_data_out_2023D_v1.root","READ");
  }
  if (run=="2024BCD") {
    lumi_136TeV = "Run2024BCD, 15.3 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v50_2024/jmenano_data_cmb_2024BCD_JME_v50_2024.root","READ"); // May 16 golden, 12.3/fb => some problem with cmb?
    //f = new TFile("rootfiles/Prompt2024/v50_2024/jmenano_data_out_2024BCD_JME_v50_2024.root","READ"); // May 16 golden, 12.3/fb
    f = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_data_out_2024BCD_JME_v83_2024.root","READ"); // Aug 2 hybrid, 15.3/fb
    //pullThreshold = 250;//200;//150;//100;//85;//100;//70;//50;
    //pullThresholdHF45 = 300;
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024E") {
    lumi_136TeV = "Run2024E, 11.3 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_2024E_JME_v76_2024.root","READ"); // June 6 hybrid, 27.0/fb => some problem with cmb?
    //f = new TFile("rootfiles/Prompt2024/v76_2024/jmenano_data_out_2024E_JME_v76_2024.root","READ"); // June 6 hybrid, 27.0/fb
    f = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_data_out_2024E_JME_v83_2024.root","READ"); // Aug 2 hybrid, 11.3/fb
    //pullThreshold = 250;//200;//150;//100;//85;//100;//70;//50;
    //pullThresholdHF45 = 300;
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024BCDE") {
    lumi_136TeV = "Run2024BCDE, 26.6 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_data_out_2024BCDE_JME_v83_2024.root","READ"); // Aug 2 hybrid, 26.6/fb
    f = new TFile("rootfiles/Prompt2024/v109_2024/jmenano_data_out_2024BCDE_JME_v109_2024.root","READ"); // V6M
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024FG") {
    //lumi_136TeV = "Run2024F, 19.4 fb^{-1}";
    lumi_136TeV = "Run2024FG, 65.5 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v86_2024/jmenano_data_out_2024F_JME_v86_2024.root","READ"); // Aug 2 hybrid
    f = new TFile("rootfiles/Prompt2024/v109_2024/jmenano_data_out_2024FG_JME_v109_2024.root","READ"); // V6M
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024BCDEFG") {
    lumi_136TeV = "Run2024BCDEFG, 92.1 fb^{-1}";
    f = new TFile("rootfiles/Prompt2024/v110_2024/jmenano_data_out_2024BCDEFG_JME_v110_2024.root","READ"); // V6M closure
    fg = new TFile("rootfiles/Prompt2024/GamHistosFill_data_2024BCDEFG_w39.root","READ"); // V6M closure
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024HI") {
    lumi_136TeV = "Run2024HI, 17.0 fb^{-1}";
    f = new TFile("rootfiles/Prompt2024/v111_2024/jmenano_data_out_2024HI_JME_v111_2024.root","READ"); // V6M
    fg = new TFile("rootfiles/Prompt2024/GamHistosFill_data_2024HIskim_w40.root","READ"); // V6M
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024BCDEFGHI") {
    lumi_136TeV = "Run2024BCDEFGHI, 109.2 fb^{-1}";
    f = new TFile("rootfiles/Prompt2024/v111_2024/jmenano_data_out_2024BCDEFG_JME_v110_2024HI_JME_v111_2024.root","READ"); // V6M closure + V6M
    fg = new TFile("rootfiles/Prompt2024/GamHistosFill_data_2024BCDEFG_w39_2024HIskim_w40.root","READ"); // V6M closure + V6M
    fpark = new TFile("rootfiles/vbfparking/vbfparking_data_2024.root","READ");
    //nMinTowers = 50; // for BPix hole
  }
  //
  if (run=="2024") {
    lumi_136TeV = "Run2024 partial re-reco, 109.2 fb^{-1}";
    f = new TFile("rootfiles/Prompt2024/v121_v2_Jet/jmenano_data_out_2024CDEFGHI_nib_Rereco_JME_v121_v2.root","READ"); // V8M
    fg = new TFile("rootfiles/Prompt2024/w48_Gam/GamHistosFill_data_2024BCDEFGHI_w48.root","READ"); // V8M
    fpark = new TFile("rootfiles/vbfparking/vbfparking_data_2024.root","READ");
  }
  /*
  if (run=="2025C") {
    lumi_136TeV = "Run2025C, 20.8 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2025/Jet_v131/jmenano_data_out_2025C_JME_v131.root","READ"); // V1M
    //f = new TFile("rootfiles/Prompt2025/Jet_v138/jmenano_data_out_2025C_JME_v138.root","READ");
    f = new TFile("rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025C_JME_v141.root","READ");
    //fg = new TFile("rootfiles/Prompt2025/Gam_w54/GamHistosFill_data_2025C_w54.root","READ"); // V1M
    //fg = new TFile("rootfiles/Prompt2025/Gam_w58/GamHistosFill_data_2025D_w58.root","READ");
    fg = new TFile("rootfiles/Prompt2025/Gam_w60/GamHistosFill_data_2025C_w60.root","READ");
    pullThreshold = 70;
    //pullThresholdHF45 = 250;
  }
  if (run=="2025D") {
    lumi_136TeV = "Run2025D, 24.0 fb^{-1}";
    f = new TFile("rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025D_JME_v141.root","READ");
    fg = new TFile("rootfiles/Prompt2025/Gam_w58/GamHistosFill_data_2025D_w60.root","READ");
    pullThreshold = 70;
  }
  if (run=="2025E") {
    lumi_136TeV = "Run2025E, X.X fb^{-1}";
    f = new TFile("rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025E_JME_v141.root","READ");
    fg = new TFile("rootfiles/Prompt2025/Gam_w58/GamHistosFill_data_2025E_w60.root","READ");
    pullThreshold = 70;
  }
  if (run=="2025E") {
    lumi_136TeV = "Run2025E, X.X fb^{-1}";
    f = new TFile("rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025E_JME_v141.root","READ");
    fg = new TFile("rootfiles/Prompt2025/Gam_w58/GamHistosFill_data_2025E_w60.root","READ");
    pullThreshold = 70;
  }
  */
  
  // Read input files from Config.C if available there
  const char *ccr = run.c_str();
  if (mfile.find(Form("JET_%s_DATA_OUT",ccr))!=mfile.end()) {
    cout << "Reading JET_" << run << "_DATA_OUT from Config.C" << endl;
    f = new TFile(mfile[Form("JET_%s_DATA_OUT",ccr)].c_str(),"READ");
  }
  if (mfile.find(Form("GAM_%s_DATA",ccr))!=mfile.end()) {
    cout << "Reading GAM_" << run << "_DATA from Config.C" << endl;
    fg = new TFile(mfile[Form("GAM_%s_DATA",ccr)].c_str(),"READ");
  }
  //if (mlum[run]!="") {
  if (mlum.find(run)!=mlum.end()) {
    cout << "Reading LUM_" << run << " from Config.C" << endl;
    lumi_136TeV = Form("%s, %s",ccr,mlum[run].c_str());
  }

  // Don't use *cmb* files, trigger folder uncertainties messed up!!
  assert(!TString(f->GetName()).Contains("_cmb_"));

  assert(f && !f->IsZombie());

  vector<string> vtrg;
  //vtrg.push_back("HLT_ZeroBias"); // analyze correct PD

  vtrg.push_back("HLT_PFJet40");
  vtrg.push_back("HLT_PFJet60");
  vtrg.push_back("HLT_PFJet80");
  vtrg.push_back("HLT_PFJet140");
  vtrg.push_back("HLT_PFJet200");
  vtrg.push_back("HLT_PFJet260");
  vtrg.push_back("HLT_PFJet320");
  vtrg.push_back("HLT_PFJet450");
  vtrg.push_back("HLT_PFJet500");
  //vtrg.push_back("HLT_PFJet550");

  //if (run!="2024BC") {
  vtrg.push_back("HLT_PFJetFwd40");
  vtrg.push_back("HLT_PFJetFwd60");
  vtrg.push_back("HLT_PFJetFwd80");
  // NB: these had wrong |eta| threshold in v12+v13. To be fixed in v14
  vtrg.push_back("HLT_PFJetFwd140");
  vtrg.push_back("HLT_PFJetFwd200");
  vtrg.push_back("HLT_PFJetFwd260");
  vtrg.push_back("HLT_PFJetFwd320");
  vtrg.push_back("HLT_PFJetFwd400");
  vtrg.push_back("HLT_PFJetFwd450");
  vtrg.push_back("HLT_PFJetFwd500");
  //} 
  //else { // 2024BCD
  vtrg.push_back("HLT_ZeroBias");
  //}
  
  // Add dijet average triggers for better h2jes
  vtrg.push_back("HLT_DiPFJetAve40");
  vtrg.push_back("HLT_DiPFJetAve60");
  vtrg.push_back("HLT_DiPFJetAve80");
  vtrg.push_back("HLT_DiPFJetAve140");
  vtrg.push_back("HLT_DiPFJetAve200");
  vtrg.push_back("HLT_DiPFJetAve260");
  vtrg.push_back("HLT_DiPFJetAve320");
  vtrg.push_back("HLT_DiPFJetAve400");
  vtrg.push_back("HLT_DiPFJetAve500");
  
  vtrg.push_back("HLT_DiPFJetAve60_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve80_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve100_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve160_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve220_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve300_HFJEC");

  if (fg) vtrg.push_back("HLT_Photon50EB_TightID_TightIso");

  // VBF Parking data from Carlos Francisco Erice Cid
  if (fpark) vtrg.push_back("HLT_VBF_DiPFJet125_45_Mjj1050");
  
  if (oneTrig!="") {
    vtrg.clear();
    vtrg.push_back(oneTrig);
  }
  int ntrg = vtrg.size();

  vector<string> vh;
  vh.push_back("h2phieta");//h2pt");
  //vh.push_back("p2mpf"); // for GamJet instead of p2asymm?
  //vh.push_back("p2asymm");
  vh.push_back("p2asymm_noveto");
  vh.push_back("p2nef");
  vh.push_back("p2chf");
  vh.push_back("p2nhf");
  if (oneHist!="") {
    vh.clear();
    vh.push_back(oneHist);
  }
  int nh = vh.size();

  map<string,TH2D*> mh2phieta;
  TH2D *h2nomsums(0), *h2nomrefsum(0), *h2jes(0);
  TH2D *h2abssums(0), *h2absrefsum(0);
  for (int ih = 0; ih != nh; ++ih) {

    string hname = vh[ih];
    const char *ch = hname.c_str();
    TH2D *h2nomsum(0), *h2abssum(0);
  
    for (int itrg = 0; itrg != ntrg; ++itrg) {

      string trg = vtrg[itrg];
      //string trg = "HLT_PFJet500";
      //string trg = "HLT_ZeroBias";
      //string hname = "p2nef";
      //isProfile2D = (hname[0]='p' && hname[1]='2');

      // Enable photon+jet and VBF parking as well
      const char *ct = trg.c_str();
      TString tt(ct);
      TFile *fs = (tt.Contains("Photon") ? fg : tt.Contains("VBF") ? fpark : f);

      // Patch: until separate photon triggers enabled, use root path
      const char *cts = (tt.Contains("Photon") ? "" : ct);
      
      string objname = Form("%s/Jetveto/%s",cts,ch);
      TObject *obj = fs->Get(objname.c_str());
      if (!obj) {
	cout << "Missing " << objname << endl << flush;
      }
      assert(obj);
      assert(obj->InheritsFrom("TH2D"));
      bool isProf2D = obj->InheritsFrom("TProfile2D");

      TH2D *h2 = (isProf2D ? ((TProfile2D*)obj)->ProjectionXY(Form("p2raw_%s_%s_%s",cr,ct,ch)) : (TH2D*)obj);

      // Clone h2 to avoid changing original, then clean eta range
      h2 = (TH2D*)h2->Clone(Form("h2raw_%s_%s_%s",cr,ct,ch));

      // Clean out unreliable parts of histograms and/or triggers
      cleanHist(h2, ch, ct, cr);
      
      h2->UseCurrentStyle();
      TH2D *h2nom = (TH2D*)h2->Clone(Form("h2nom_%s_%s_%s",cr,ct,ch));
      TH2D *h2abs = (TH2D*)h2->Clone(Form("h2abs_%s_%s_%s",cr,ct,ch));

      // Store h2phieta for later filtering out low-statistics bins
      if (hname=="h2phieta") {
	mh2phieta[ct] = h2;
      }
      
      // Calculate average JES shift (don't use VBF for this)
      if ((hname=="p2asymm" || hname=="p2asymm_noveto") &&
	  !tt.Contains("VBF")) {
	
	if (!h2jes) {
	  h2jes = (TH2D*)h2->Clone(Form("h2jes_%s",cr)); h2jes->Reset();
	}
	//else {
	TH2D *h2phieta = mh2phieta[ct];
	for (int i = 1; i != h2->GetNbinsX()+1; ++i) {

	  // Check that eta strip well covered before using it
	  int ntower(0);
	  for (int j = 0; j != h2->GetNbinsY()+1; ++j) {
	    if (h2->GetBinError(i,j)!=0) {
	      if (!h2phieta || h2phieta->GetBinContent(i,j)>=nMinJES) {
		++ntower;
	      }
	    }
	  } // for j
	  if (ntower<=nMinTowersJES) continue;
	    
	  for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	    //if (h2->GetBinError(i, j)!=0) {
	    double val1 = h2jes->GetBinContent(i,j);
	    double err1 = h2jes->GetBinError(i,j);
	    double n1 = (err1!=0 ? 1./pow(err1,2) : 0);
	    double val2 = h2->GetBinContent(i,j);
	    double err2 = h2->GetBinError(i,j);
	    double n2 = (err2!=0 ? 1./pow(err2,2) : 0);
	    double val = (n1+n2!=0 ? (n1*val1 + n2*val2) / (n1+n2) : 0);
	    double err = (n1+n2!=0 ? (n1*err1 + n2*err2) / (n1+n2) : 0);

	    // Check that tower has enough entries to not mess with err
	    if (h2phieta && h2phieta->GetBinContent(i,j)<nMinJES)
	      continue;
	    
	    h2jes->SetBinContent(i,j,val);
	    h2jes->SetBinError(i,j,err);
	    //}
	  } // for j
	} // for i
	//}
	
	// JES with normalized eta strips
	// Useful for later phi symmetry studies and jet substructure
	TH2D *h2jesnorm = (TH2D*)h2->Clone(Form("h2jesnorm_%s_%s",cr,cts));
	for (int i = 1; i != h2jesnorm->GetNbinsX()+1; ++i) {
	  double norm = h2jesnorm->Integral(i,i,1,72) / 72.;
	  for (int j = 1; j != h2jesnorm->GetNbinsY()+1; ++j) {
	    h2jesnorm->SetBinContent(i,j,(1+h2jesnorm->GetBinContent(i,j))/(1+norm)-1);
	    h2jesnorm->SetBinError(i,j,h2jesnorm->GetBinError(i,j)/(1+norm));
	  }
	}

	// Save JES results to output file
	fout->cd("trigs");
	h2->Write(Form("jetasymmetrymap_%s",ct));
	h2jesnorm->Write(Form("jetasymmetrymap_norm_%s",ct));
	curdir->cd();
	
      } // p2asymm
      
      // Normalize eta strips vs phi
      for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
	
	// Recreate histogram with automatic binning for each eta bin
	TH1D *htmp = new TH1D("htmp","Distribution of values",100,-1,-1);
	htmp->SetBuffer(72.);

	TH2D *h2phieta = mh2phieta[ct];
	int ny = h2->GetNbinsY();
	int n(0);
	// Check that eta strip is well covered before using it
	for (int j = 1; j != ny+1; ++j) {
	  if (h2->GetBinError(i,j)!=0) {
	    if (!h2phieta || h2phieta->GetBinContent(i,j)>=nMinPull) {
	      htmp->Fill(h2->GetBinContent(i, j));
	      ++n;
	    }
	  }
	} // for j
	//htmp->Draw();
	
	// If not enough entries, clear column and continue
	//if (htmp->GetEntries()<70. || htmp->Integral()!=0) {
	if (n<nMinTowers) {
	  for (int j = 1; j != ny+1; ++j) {
	    h2->SetBinContent(i, j, 0.);
	    h2->SetBinError(i, j, 0.);
	  } // for j
	}
	// 
	else {
	  // Determine median and 68% quantile ([0.16,0.84] range divided by 2)
	  const Int_t nprobSum = 3;
	  Double_t probSum[nprobSum] = {0.16, 0.50, 0.84};
	  Double_t q[nprobSum];
	  htmp->GetQuantiles(nprobSum, q, probSum);
	  double median = q[1];
	  double q68 = 0.5*(q[2]-q[0]);
	  
	  // Determine pull width
	  for (int j = 1; j != ny+1; ++j) {
	    if (h2->GetBinError(i,j)!=0 && q68!=0 &&
		(!h2phieta || h2phieta->GetBinContent(i,j)>nMinPull)) {
	      if (doPull) {
		h2abs->SetBinContent(i, j, fabs(h2->GetBinContent(i, j)-median) / q68 - 1); // -1 to model the normal cancellation of +1 and -1 fluctuations
		h2abs->SetBinError(i, j, h2->GetBinError(i, j) / q68);
		h2nom->SetBinContent(i, j, (h2->GetBinContent(i, j)-median) / q68);
		h2nom->SetBinError(i, j, h2->GetBinError(i, j) / q68);
		
	      }
	      else {
		// TBD: code not checked for h2abs
		if (median<0.5) {
		  h2abs->SetBinContent(i,j,fabs(h2->GetBinContent(i,j)-median)/
				       (1+fabs(median)));
		  h2abs->SetBinError(i, j, h2->GetBinError(i, j) /
				     (1+fabs(median)));
		  h2nom->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/
				       (1+fabs(median)));
		  h2nom->SetBinError(i, j, h2->GetBinError(i, j) /
				     (1+fabs(median)));

		}
		else {
		  h2abs->SetBinContent(i,j,fabs(h2->GetBinContent(i,j)-median)/
				       median);
		  h2abs->SetBinError(i, j, h2->GetBinError(i, j) / median);
		  h2nom->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/median);
		  h2nom->SetBinError(i, j, h2->GetBinError(i, j) / median);
		}
	      }
	    }
	  } // for j
	}
	
	delete htmp;
      } // for i
      
      //h2->Draw("COLZ");
      if (!h2nomsum && !h2abssum) {
	h2nomsum = (TH2D*)h2nom->Clone(Form("h2nomsum_%s_%s",cr,ch));
	h2abssum = (TH2D*)h2abs->Clone(Form("h2abssum_%s_%s",cr,ch));
      }
      else {
	h2nomsum->Add(h2nom);
	h2abssum->Add(h2abs);
      }

      // Save results to output file
      fout->cd("trigs");
      h2nom->Write(Form("jetpullmap_nom_%s_%s",hname.c_str(),trg.c_str()));
      curdir->cd();
      //delete h2;
					
    } // for itrg

    double ntrg2 = (oneTrig!="" ? 1 : max(1, ntrg/2));

    const char *c1nomname = Form("c1nom_%s_%s",cr,ch);
    TH1D *h1nom = tdrHist(Form("h1nom_%s_%s_%s",cr,oneTrig.c_str(),ch),
			  "#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
    //if (doPull) h1nom->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    //else        h1nom->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c1nom = tdrCanvas(c1nomname,h1nom,8,11,kRectangular);
    c1nom->SetRightMargin(0.15);
    h2nomsum->Draw("COLZ SAME");
    if (doPull) h2nomsum->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h2nomsum->GetZaxis()->SetRangeUser(-0.5,0.5);

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.035);
    tex->DrawLatex(0.15,0.87,hname.c_str());
    tex->DrawLatex(0.15,0.83,oneTrig.c_str());

    gPad->RedrawAxis();
    gPad->Update();

    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_Run%s.pdf",
		       hname.c_str(), doPull ? "nompull" : "nomrel", run.c_str()));


    const char *c1absname = Form("c1abs_%s_%s",cr,ch);
    TH1D *h1abs = tdrHist(Form("h1abs_%s_%s_%s",cr,oneTrig.c_str(),ch),
			  "#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
    //if (doPull) h1abs->GetZaxis()->SetRangeUser(0,5*ntrg2);
    //else        h1abs->GetZaxis()->SetRangeUser(0.,0.50);

    TCanvas *c1abs = tdrCanvas(c1absname,h1abs,8,11,kRectangular);
    c1abs->SetRightMargin(0.15);
    h2abssum->Draw("COLZ SAME");
    if (doPull) h2abssum->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h2abssum->GetZaxis()->SetRangeUser(-0.5,0.5);
    //if (doPull) h2abssum->GetZaxis()->SetRangeUser(0,5*ntrg2);
    //else        h2abssum->GetZaxis()->SetRangeUser(0,0.5);

    tex->DrawLatex(0.15,0.87,hname.c_str());
    tex->DrawLatex(0.15,0.83,oneTrig.c_str());

    gPad->RedrawAxis();
    gPad->Update();

    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_Run%s.pdf",
		       hname.c_str(), doPull ? "abspull" : "absrel",
		       run.c_str()));


    
    if (!h2nomrefsum && !h2absrefsum &&
	(hname=="p2asymm" || hname=="p2asymm_noveto")) {
      h2nomrefsum = (TH2D*)h2nomsum->Clone("h2nomrefsum");
      h2absrefsum = (TH2D*)h2abssum->Clone("h2absrefsum");
    }

    if (!h2nomsums && !h2abssums) {
      h2nomsums = (TH2D*)h2nomsum->Clone("h2nomsums");
      h2abssums = (TH2D*)h2abssum->Clone("h2abssums");
    }
    else {
      h2nomsums->Add(h2nomsum);
      h2abssums->Add(h2abssum);
    }
  } // for ih

  TH1D *h1nom = tdrHist(Form("h1nom_%s",cr),
			"#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
  TH1D *h1abs = tdrHist(Form("h1abs_%s",cr),
			"#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
  /*
  if (doPull) {
    h1nom->GetZaxis()->SetRangeUser(-100,100);
    h1abs->GetZaxis()->SetRangeUser(0,100);
    //if (run=="2024BCD") {
    //h1nom->GetZaxis()->SetRangeUser(-150,150);
    //h1abs->GetZaxis()->SetRangeUser(0,150);
    //}
  }
  else {
    h1nom->GetZaxis()->SetRangeUser(-0.50,0.50);
    h1abs->GetZaxis()->SetRangeUser(0.,0.50);
  }
  */

  TCanvas *c1nom = tdrCanvas(Form("c1nom_%s",cr),h1nom,8,11,kRectangular);
  gPad->SetRightMargin(0.15);
  h2nomsums->Draw("COLZ SAME");

  TCanvas *c1abs = tdrCanvas(Form("c1abs_%s",cr),h1abs,8,11,kRectangular);
  gPad->SetRightMargin(0.15);
  h2abssums->Draw("COLZ SAME");

  if (doPull) {
    h2nomsums->GetZaxis()->SetRangeUser(-100,100);
    h2abssums->GetZaxis()->SetRangeUser(-100,100);
    //h2abssums->GetZaxis()->SetRangeUser(0,100);
    //if (run=="2024BCD") {
    //h2nomsums->GetZaxis()->SetRangeUser(-150,150);
    //h2abssums->GetZaxis()->SetRangeUser(-150,150);
    //}
  }
  else {
    h2nomsums->GetZaxis()->SetRangeUser(-0.50,0.50);
    h2abssums->GetZaxis()->SetRangeUser(-0.50,0.50);
    //h2abssums->GetZaxis()->SetRangeUser(0,0.50);
  }

  
  TH2D *h2veto = (TH2D*)h2nomsums->Clone("jetvetomap"); // default map
  TH2D *h2hot = (TH2D*)h2nomsums->Clone("jetvetomap_hot");
  TH2D *h2cold = (TH2D*)h2nomsums->Clone("jetvetomap_cold");
  TH2D *h2old = (TH2D*)h2nomsums->Clone("jetvetomap_old");
  TH2D *h2hotandcold = (TH2D*)h2nomsums->Clone("jetvetomap_hotandcold");
  TH2D *h2eep = (TH2D*)h2nomsums->Clone("jetvetomap_eep");
  TH2D *h2bpix = (TH2D*)h2nomsums->Clone("jetvetomap_bpix");
  TH2D *h2fpix = (TH2D*)h2nomsums->Clone("jetvetomap_fpix");
  TH2D *h2all = (TH2D*)h2nomsums->Clone("jetvetomap_all");
  h2veto->Reset();
  h2hot->Reset();
  h2cold->Reset();
  h2old->Reset();
  h2hotandcold->Reset();
  h2eep->Reset();
  h2bpix->Reset();
  h2fpix->Reset();
  h2all->Reset();

  h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold");
  h2hot ->SetTitle("Hot zones. Use for steep jet pT spectra and MET tails");
  h2cold->SetTitle("Cold zones. Mostly ECAL holes. Use for MET tails");
  h2old->SetTitle("Old zones. Already in previous map.");
  h2hotandcold->SetTitle("Union of hot and cold zone maps");
  h2eep->SetTitle("EE+ water leak region. Complete loss of ECAL energy in 2022EFG");
  h2bpix->SetTitle("Barrel pixel failure. Reduction of tracking efficiency in 2023D and later");
  h2fpix->SetTitle("Forward pixel failure. Reduction of tracking efficiency in 2024F and later");
  h2all->SetTitle("Union of hot and cold maps");
  if (run=="2022E" || run=="2022F" || run=="2022G" ||
      run=="2022EF" || run=="2022EFG") {
    h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold+EEP");
    h2all->SetTitle("Union of hot, cold and EEP maps");
  }
  if (run=="2023D") {
    h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold+BPIX");
    h2all->SetTitle("Union of hot, cold and BPIX maps");
  }
  if (run=="2024BCD" || run=="2024BCDE") {
    h2veto->SetTitle("JME recommended map. Hot+Cold(+not BPIX)");
    h2all->SetTitle("Union of all maps, used for JEC. Hot+cold+BPIX");
  }
  if (run=="2024F" || run=="2024G" || run=="2024FG" ||
      run=="2024BCDEFG" || run=="2024HI" || run=="2024BCDEFGHI" ||
      run=="2024") {
    h2veto->SetTitle("JME recommended map. Hot+Cold+FPIX(+not BPIX)");
    h2all->SetTitle("Union of all maps, used for JEC. Hot+cold+FPIX+BPIX");
  }
  if (run=="2025C" || run=="2025D" || run=="2025E" || run=="2025CDE" ||
      run=="2025F" || run=="2025CDEF") {
    h2veto->SetTitle("JME recommended map. Hot+Cold+FPIX(+not BPIX)");
    h2all->SetTitle("Union of all maps, used for JEC. Hot+cold+FPIX+BPIX");
  }
  h2nomsums->SetTitle("Raw nominal pull map. NomPull=(x_i-mu)/sigma. Summed over triggers");
  h2abssums->SetTitle("Raw absolute pull map. AbsPull=|(x_i-mu)/sigma|-1. Summed over triggers");
  h2nomrefsum->SetTitle("Nominal pull map from dijet asymmetry. Used for hot vs cold classification");

  for (int i = 1; i != h2nomsums->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2nomsums->GetNbinsY()+1; ++j) {
      double eta = h2nomsums->GetXaxis()->GetBinCenter(i);
      double phi = h2nomsums->GetYaxis()->GetBinCenter(j);

      // Match positive asymmetry with hot region
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (!h2nomrefsum || h2nomrefsum->GetBinContent(i,j)>0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2hot->SetBinContent(i,j,100);
      }

      // Match negative asymmetry with cold region
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (h2nomrefsum && h2nomrefsum->GetBinContent(i,j)<0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2cold->SetBinContent(i,j,-100);
      }

      // Match zero asymmetry with old region (probably had p2asymm masked)
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (h2nomrefsum && h2nomrefsum->GetBinContent(i,j)==0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2old->SetBinContent(i,j,-100);
      }

      // Extra manual margin for EE+ water leak region (2022E,F,G)
      if (eta>1.5 && phi>1.85 &&
	  eta<2.2 && phi<2.7 &&
	  phi<2.7-(2.7-2.1)/(2.1-1.6)*(eta-1.6)) {
	if (run=="2022E" || run=="2022F" || run=="2022G" ||
	    run=="2022EF" || run=="2022EFG") {
	  h2veto->SetBinContent(i,j,100);
	  h2eep->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
      } // EEP

      // BPIX reference region (2023D and later)
      if (eta>-1.5 && phi>-1.22 &&
	  eta<0.1 && phi<-0.87) {
	if (run=="2023D") { // Remove BPix
	  h2veto->SetBinContent(i,j,100);
	  h2bpix->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
	if (run=="2024BCD" || run=="2024BCDE" || run=="2024F" ||
	    run=="2024G" || run=="2024FG" ||
	    run=="2024BCDEFG" || run=="2024HI" || run=="2024BCDEFGHI" ||
	    run=="2024") {
	  // Keep BPix for recommended, drop for JEC
	  h2veto->SetBinContent(i,j,0);   // keep!
	  h2bpix->SetBinContent(i,j,100); // drop
	  h2all->SetBinContent(i,j,100);  // drop
	}
	if (run=="2025C" || run=="2025D" || run=="2025E" || run=="2025CDE" ||
	    run=="2025F" || run=="2025CDEF") {
	  // Keep BPix for recommended, drop for JEC
	  h2veto->SetBinContent(i,j,0);   // keep!
	  h2bpix->SetBinContent(i,j,100); // drop
	  h2all->SetBinContent(i,j,100);  // drop
	}
      } // BPIX

      // FPIX reference region (2024F and later)
      if (eta>-2.043 && phi>2.53 &&
	  eta<-1.653 && phi<2.71) {
	if (run=="2024F" || run=="2024G" || run=="2024FG" ||
	    run=="2024BCDEFG" || run=="2024HI" || run=="2024BCDEFGHI" ||
	    run=="2024") {
	  // Remove FPIX
	  h2veto->SetBinContent(i,j,100);
	  h2fpix->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
	if (run=="2025C" || run=="2025D" || run=="2025E" || run=="2025CDE" ||
	    run=="2025F" || run=="2025CDEF"){
	  h2veto->SetBinContent(i,j,100);
	  h2fpix->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
      } // FPIX

    } // for j
  } // for j

  if (plotJetVeto) {
    h2veto->SetLineColor(kRed);
    h2hotandcold->SetLineColor(kOrange+1);
    h2hot->SetLineColor(kRed);

    c1nom->cd();
    h2veto->Draw("SAME BOX");
    h2hotandcold->Draw("SAME BOX");
    h2hot->Draw("SAME BOX");

    c1abs->cd();
    h2veto->Draw("SAME BOX");
    h2hotandcold->Draw("SAME BOX");
    h2hot->Draw("SAME BOX");
  }
  gPad->RedrawAxis();
  gPad->Update();
  
  if (oneTrig!="" && oneHist!="") {
    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_%s_Run%s.pdf",
		       oneTrig.c_str(), oneHist.c_str(),
		       doPull ? "nompull" : "absrel",run.c_str()));
    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_%s_Run%s.pdf",
		       oneTrig.c_str(), oneHist.c_str(),
		       doPull ? "abspull" : "absrel",run.c_str()));
  }
  else {
    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_h2nommap_Run%s.pdf",run.c_str()));
    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_h2absmap_Run%s.pdf",run.c_str()));
  }

  TH2D *h2jesnorm(0);
  if (h2jes) {
    TH1D *h2 = tdrHist(Form("h2_%s",cr),
		       "#phi",-TMath::Pi(),+TMath::Pi(),"#eta",-5.2,5.2);
    //h2->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c2 = tdrCanvas(Form("c2_%s",cr),h2,8,11,kRectangular);
    gPad->SetRightMargin(0.15);
    h2jes->Draw("SAME COLZ");
    h2jes->GetZaxis()->SetRangeUser(-0.15,0.15);
    gPad->RedrawAxis();
    gPad->Update();

    if (plotJetVetoOnJES) {
      h2veto->SetLineColor(kRed);
      h2veto->Draw("SAME BOX");
    }
      
    c2->SaveAs(Form("pdf/JetVeto/JetVeto_x2jes_%s.pdf",run.c_str()));

    h2jes->SetTitle("Dijet asymmetry map, (pTprobe-pTtag)/pTave for |eta,tag|<1.3");


    // Normalize eta strips
    h2jesnorm = (TH2D*)h2jes->Clone(Form("h2jesnorm_%s",cr));
    for (int i = 1; i != h2jes->GetNbinsX()+1; ++i) {
      double norm = h2jes->Integral(i,i,1,72) / 72.;
      for (int j = 1; j != h2jes->GetNbinsY()+1; ++j) {
	h2jesnorm->SetBinContent(i,j,(1+h2jes->GetBinContent(i,j))/(1+norm)-1);
	h2jesnorm->SetBinError(i,j,h2jes->GetBinError(i,j)/(1+norm));
      }
    }
    
    TH1D *h2norm = tdrHist(Form("h2norm_%s",cr),
			   "#phi",-TMath::Pi(),+TMath::Pi(),"#eta",-5.2,5.2);
    //h2norm->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c2norm = tdrCanvas(Form("c2norm_%s",cr),h2norm,8,11,kRectangular);
    gPad->SetRightMargin(0.15);
    h2jesnorm->Draw("SAME COLZ");
    h2jesnorm->GetZaxis()->SetRangeUser(-0.15,0.15);
    gPad->RedrawAxis();
    gPad->Update();

    if (plotJetVetoOnJESnorm) {
      h2veto->SetLineColor(kRed);
      h2veto->Draw("SAME BOX");
    }
      
    c2norm->SaveAs(Form("pdf/JetVeto/JetVeto_x2jesnorm_%s.pdf",run.c_str()));

    h2jesnorm->SetTitle("Dijet asymmetry map, (pTprobe-pTtag)/pTave for |eta,tag|<1.3 normalized vs #phi_{probe}");
  } // if (h2jes)

  fout->cd();
  h2veto->Write("jetvetomap");
  h2hot->Write("jetvetomap_hot");
  h2cold->Write("jetvetomap_cold");
  if (h2old->Integral()!=0) h2old->Write("jetvetomap_old");
  h2hotandcold->Write("jetvetomap_hotandcold");
  if (h2eep->Integral()!=0) h2eep->Write("jetvetomap_eep");
  if (h2bpix->Integral()!=0) h2bpix->Write("jetvetomap_bpix");
  if (h2fpix->Integral()!=0) h2fpix->Write("jetvetomap_fpix");
  h2all->Write("jetvetomap_all");
  if (h2jes) h2jes->Write("jetasymmetrymap");
  if (h2jesnorm) h2jesnorm->Write("jetasymmetrymap_norm");
  h2nomsums->Write("jetpullsummap_nom");
  h2abssums->Write("jetpullsummap_abs");
  if (h2nomrefsum) h2nomrefsum->Write("jetpullsummap_ref");
  fout->Write();
  fout->Close();
} // JetVeto


// Clean out unreliable ranges of histograms, e.g. NHF in HF
void cleanHist(TH2D *h2, string hist, string trg, string run) {

  TString th(hist.c_str());
  TString tt(trg.c_str());
  TString tr(run.c_str());
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      double eta = h2->GetXaxis()->GetBinCenter(i);
      double abseta = fabs(eta);
      if ((th.Contains("p2nhf") && abseta>2.964) ||
	  (th.Contains("p2nef") && abseta>2.964) ||
	  (th.Contains("p2chf") && abseta>2.650) ||
	  
	  (tt.Contains("HFJEC") && abseta<2.964) ||
	  (tt.Contains("Fwd") && abseta<2.964) ||

	  // Trigger ranges assessed from jetpullmap_nom_h2phieta
	  (tt.Contains("PFJet40") && abseta>2.964) ||
	  (tt.Contains("PFJet60") && abseta>2.964) ||
	  (tt.Contains("PFJet80") && abseta>2.964) ||
	  (tt.Contains("PFJet140") && abseta>2.964) ||
	  (tt.Contains("PFJet200") && abseta>2.964) ||
	  (tt.Contains("PFJet260") && abseta>2.964) ||
	  (tt.Contains("PFJet320") && abseta>2.964) ||
	  (tt.Contains("PFJet450") && abseta>2.853) ||
	  (tt.Contains("PFJet500") && abseta>2.853) ||

	  (tt.Contains("PFJetFwd60") && abseta>4.716) ||
	  (tt.Contains("PFJetFwd80") && abseta>4.538) ||
	  (tt.Contains("PFJetFwd140") && abseta>4.191) ||
	  (tt.Contains("PFJetFwd200") && abseta>4.013) ||
	  (tt.Contains("PFJetFwd260") && abseta>3.664) ||
	  (tt.Contains("PFJetFwd320") && abseta>3.314) ||
	  (tt.Contains("PFJetFwd400") && abseta>3.139) ||
	  (tt.Contains("PFJetFwd450") && abseta>3.139) ||
	  (tt.Contains("PFJetFwd500") && abseta>2.964) ||

	  // ZeroBias ok
	  // DiPFJetAve 40 ok
	  (tt.Contains("DiPFJetAve60") && abseta>4.538) ||
	  (tt.Contains("DiPFJetAve80") && abseta>4.538) ||
	  (tt.Contains("DiPFJetAve140") && abseta>4.013) ||
	  (tt.Contains("DiPFJetAve200") && abseta>3.839) ||
	  (tt.Contains("DiPFJetAve260") && abseta>3.489) ||
	  (tt.Contains("DiPFJetAve320") && abseta>2.853) ||
	  (tt.Contains("DiPFJetAve400") && abseta>2.853) ||
	  (tt.Contains("DiPFJetAve500") && abseta>2.853) ||
	  
	  (tt.Contains("DiPFJetAve60_HFJEC") && abseta>4.716) ||
	  (tt.Contains("DiPFJetAve80_HFJEC") && abseta>4.538) ||
	  (tt.Contains("DiPFJetAve100_HFJEC") && abseta>4.538) ||
	  (tt.Contains("DiPFJetAve160_HFJEC") && abseta>4.191) ||
	  (tt.Contains("DiPFJetAve220_HFJEC") && abseta>3.839) ||
	  (tt.Contains("DiPFJetAve300_HFJEC") && abseta>3.664) ||

	  (tr.Contains("2025") && tt.Contains("Photon") && abseta>4.716) ||
	  (tr.Contains("2025") && tt.Contains("PFJet140") && abseta>2.853) ||
	  (tr.Contains("2025") && tt.Contains("PFJet200") && abseta>2.853) ||
	  (tr.Contains("2025") && tt.Contains("PFJet260") && abseta>2.853) ||
	  (tr.Contains("2025") && tt.Contains("PFJetFwd40") && abseta>4.716) ||
	  (tr.Contains("2025") && tt.Contains("PFJetFwd60") && abseta>4.191) ||
	  (tr.Contains("2025") && tt.Contains("PFJetFwd80") && abseta>4.191) ||
	  (tr.Contains("2025") && tt.Contains("Ave40") && abseta>4.716) ||
	  (tr.Contains("2025") && tt.Contains("Ave200") && abseta>2.853) ||
	  (tr.Contains("2025") && tt.Contains("Ave260") && abseta>2.853) ||
	  (tr.Contains("2025") && tt.Contains("Ave60_HFJEC") && abseta>4.191) ||
	  (tr.Contains("2025") && tt.Contains("Ave80_HFJEC") && abseta>4.191) ||
	  (tr.Contains("2025") && tt.Contains("Ave100_HFJEC") && abseta>4.191)
	  
	  ) {
	h2->SetBinContent(i,j,0);
	h2->SetBinError(i,j,0);
      }
    } // for j
  } // for i
}
