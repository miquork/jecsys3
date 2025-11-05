// Purpose: Draw Z+jet, gamma+jet and inclusive jet xsec and MPF vs time
//          Use breakup made by rootfiles/brilcalc/clusterRuns.C
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

#include "../tdrstyle_mod22.C"

bool debug = false;
bool addExtraText = false;
bool addExtraLine = false;

bool addNibLine = true;
bool addYearLine = true;

void rebinProfileCustom(TProfile* p, TProfile* p_new);
map<string,pair<int,int>> BuildRanges(const string& fname);
void addLines(map<string,pair<int,int>>& ml, TH1D *hcumlum2, double kf,
	      double xmin, double xmax, double ymin, double ymax,
	      TH1D *hbreaks = 0);

void drawTimeStability() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fout = new TFile("rootfiles/drawTimeStability.root","RECREATE");
  
  TFile *flum = new TFile("rootfiles/brilcalc/clusterRuns.root","READ");
  assert(flum && !flum->IsZombie());

  // Load luminosity information
  TH1D *hlum = (TH1D*)flum->Get("hlum"); assert(hlum); // per run
  TH1D *hlum2 = (TH1D*)flum->Get("hlum2"); assert(hlum2); // per range
  TH1D *hcumlum = (TH1D*)flum->Get("hcumlum"); assert(hcumlum); // per run
  TH1D *hcumlum2 = (TH1D*)flum->Get("hcumlum2"); assert(hcumlum2); // per range
  TH1D *hbins = (TH1D*)flum->Get("hcumlumbins2"); assert(hbins); // per range

  // Listing all files here. GamJet and ZmmJet switch filetype
  string vf[] = {

    // dijet
    "DiJet",
    "2022_v113","2023_v113",
    //"2024C_nib1","2024D_nib1","2024E_nib1",
    //"2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    //"2024H_nib1","2024I_nib1",
    "2024CDEFGHI_nib",
    //
    "2025C","2025D","2025E", "2025F","2025G",
    
    // gamma+jet
    "GamJet",
    "2022C","2022D","2022E","2022F","2022G",
    "2023Cv123","2023Cv4","2023D",
    "2024C_nib1-rereco","2024D_nib1-rereco","2024E_nib1-rereco",
    "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    "2024H_nib1","2024I_nib1",
    //
    "2025C","2025D","2025E","2025F","2025G",

    // Z+jet
    "ZmmJet",
    "2022CD","2022E","2022F","2022G",
    "2023C123","2023C4","2023D",
    "2024CDEReprocessing_v1_2024C_nib1",
    "2024CDEReprocessing_v1_2024D_nib1",
    "2024CDEReprocessing_v1_2024Ev1_nib1",
    "2024CDEReprocessing_v1_2024Ev2_nib1",
    "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    "2024H_nib1","2024I_nib1",
    //
    "2025C","2025D","2025E","2025F","2025G",

    "TTBar",
    "2024"
  };
  const int nf = sizeof(vf)/sizeof(vf[0]);

  // Listing all histograms here. GamJet and ZmmJet switch filetype
  string vh[] = {

    // Dijet
    "DiJet",
    "h1jetrate",


    // Gamma+jet
    "GamJet",
    //"pr50m_eta08hi","pr50m_eta08lo",
    //"pr50m_eta3to4","pr50m_eta4to5",
    //"pr110m_eta3to4","pr110m_eta4to5",

    "pr50n","pr110n","pr230n",
 
    "pr230m","pr110m","pr50m",
    "pr230b","pr110b","pr50b",
    "pr50nhf","pr50nef","pr50chf",
    "pr110nhf","pr110nef","pr110chf",
    "pr230nhf","pr230nef","pr230chf",


    // Zmm+jet
    "ZmmJet",
    "mz_run_zpt0","mz_run_zpt30","mz_run_zpt50","mz_run_zpt110",
    "nz_run_zpt0","nz_run_zpt30","nz_run_zpt50","nz_run_zpt110",

    "mpf_run_zpt30","mpf_run_zpt110","mpf_run_zpt50",
    "db_run_zpt30","db_run_zpt110","db_run_zpt50",

    //"chf_run_zpt30","nef_run_zpt30","nhf_run_zpt30",
    //"chf_run_zpt50","nef_run_zpt50","nhf_run_zpt50",
    //"chf_run_zpt110","nef_run_zpt110","nhf_run_zpt110",


    // ttbar
    //"TTBar",
    //"prof_W","prof_W_inWindow",
    //"prof_top","prof_top_inWindow","prof_top_improved"

  };
  const int nh = sizeof(vh)/sizeof(vh[0]);

  map<string, int> mcolor;
  mcolor["pr50n"] = kRed;
  mcolor["pr110n"] = kBlue;
  mcolor["pr230n"] = kGreen+2;
  mcolor["pr50m"] = kRed;
  mcolor["pr110m"] = kBlue;
  mcolor["pr230m"] = kGreen+2;
  mcolor["pr50chf"] = kRed;
  mcolor["pr50nhf"] = kGreen+2;
  mcolor["pr50nef"] = kBlue;
  mcolor["pr110chf"] = kRed;
  mcolor["pr110nhf"] = kGreen+2;
  mcolor["pr110nef"] = kBlue;
  mcolor["mpf_run_zpt30"] = kRed;
  mcolor["mpf_run_zpt50"] = kBlue;
  mcolor["mpf_run_zpt110"] = kGreen+2;
  mcolor["chf_run_zpt30"] = kRed;
  mcolor["nef_run_zpt30"] = kBlue;
  mcolor["nhf_run_zpt30"] = kGreen+2;

  map<string,const char*> mhead;
  mhead["pr50n"] = "#gamma+jet Xsec 50EB";
  mhead["pr110n"] = "#gamma+jet Xsec 110EB";
  mhead["pr230n"] = "#gamma+jet Xsec 200";
  mhead["pr50m"] = "#gamma+jet MPF 50EB";
  mhead["pr110m"] = "#gamma+jet MPF 110EB";
  mhead["pr230m"] = "#gamma+jet MPF 200";
  mhead["pr50b"] = "#gamma+jet DB 50EB";
  mhead["pr110b"] = "#gamma+jet DB 110EB";
  mhead["pr230b"] = "#gamma+jet DB 200";
  mhead["pr50chf"] = "#gamma+jet CHF 50EB";
  mhead["pr50nhf"] = "#gamma+jet NHF 50EB";
  mhead["pr50nef"] = "#gamma+jet NEF 50EB";
  mhead["pr110chf"] = "#gamma+jet CHF 110EB";
  mhead["pr110nhf"] = "#gamma+jet NHF 110EB";
  mhead["pr110nef"] = "#gamma+jet NEF 110EB";
  mhead["db_run_zpt30"] = "Z(#mu#mu)+jet DB 30";
  mhead["db_run_zpt50"] = "Z(#mu#mu)+jet DB 50";
  mhead["db_run_zpt110"] = "Z(#mu#mu)+jet DB 110";
  mhead["mpf_run_zpt30"] = "Z(#mu#mu)+jet MPF 30";
  mhead["mpf_run_zpt50"] = "Z(#mu#mu)+jet MPF 50";
  mhead["mpf_run_zpt110"] = "Z(#mu#mu)+jet MPF 110";
  mhead["chf_run_zpt30"] = "Z(#mu#mu)+jet CHF 30";
  mhead["nef_run_zpt30"] = "Z(#mu#mu)+jet NEF 30";
  mhead["nhf_run_zpt30"] = "Z(#mu#mu)+jet NHF 30";
  

  // Store normalizations to histogra
  TH1D *hscales = new TH1D("hscales",";Observable;Scale",2*nh,0,2*nh);

  // Load year, era, nib ranges from file
  map<string,pair<int,int> > ml = BuildRanges("rootfiles/brilcalc/fibs.txt");
  
  // Loop over files to retrieve stuff
  bool useGam(false), useZmm(false), useJet(false), useTT(false);
  for (int ih = 0; ih != nh; ++ih) {

    string sh = vh[ih];
    const char *ch = sh.c_str();
    TString th(ch);
    cout << "Analyzing " << ch << endl << flush;
    // Select the correct flag for following profiles/histograms
    // Then continue to move to next entries (profile/histogram names)
    if (sh=="GamJet") { useGam = true; useZmm = false; useJet = false; useTT = false; continue; }
    if (sh=="ZmmJet") { useGam = false; useZmm = true; useJet = false; useTT = false; continue; }
    if (sh=="DiJet")  { useGam = false; useZmm = false; useJet = true; useTT = false; continue; }
    if (sh=="TTBar")  { useGam = false; useZmm = false; useJet = false; useTT = true; continue; }
    bool isPF = (th.Contains("nhf") || th.Contains("nef") || th.Contains("chf"));
    bool isZmass = (th.Contains("mz_run"));
    bool isEta = (th.Contains("eta"));
    double kpf(1);
    if (th.Contains("nhf")) kpf = 20;  // 5%
    if (th.Contains("nef")) kpf = 4;   // 25%
    if (th.Contains("chf")) kpf = 1.5; // 65%
    
    TH1D *hsum = (TH1D*)hlum2->Clone(Form("hsum_%s",ch)); hsum->Reset();
    double vx[hbins->GetNbinsX()+1];
    for (int i = 1; i != hbins->GetNbinsX()+2; ++i) vx[i-1] = hbins->GetBinLowEdge(i);
    TProfile *psum = new TProfile(Form("psum_%s",ch),"",hbins->GetNbinsX(),&vx[0]);
    psum->Sumw2();
    TProfile *psumjes = new TProfile(Form("psumjes_%s",ch),"",hbins->GetNbinsX(),&vx[0]);
    psumjes->Sumw2();
    
    // Keep track of JES (1=L2L3Res) for gamma+jet MPF
    //TH1D *hjes = new TH1D(Form("hjes_%s",ch),";JES;Cum.Lum(/fb)",hbins->GetNbinsX(),&vx[0]);

    bool isGam(false), isZmm(false), isJet(false), isTT(false);
    for (int iff = 0; iff != nf; ++iff) {

      string sf = vf[iff];
      const char *cf = sf.c_str();
      TString tf(cf);
      
      if (sf=="GamJet") { isGam=true; isZmm=false; isJet=false; isTT=false; continue; }
      if (sf=="ZmmJet") { isGam=false; isZmm=true; isJet=false; isTT=false; continue; }
      if (sf=="DiJet")  { isGam=false; isZmm=false; isJet=true; isTT=false; continue; }
      if (sf=="TTBar")  { isGam=false; isZmm=false; isJet=false; isTT=true; continue; }
      if (useGam && !isGam) continue;
      if (useZmm && !isZmm) continue;
      if (useJet && !isJet) continue;
      if (useTT  && !isTT) continue;
      
      TFile *f(0);
      // Photon+jet files
      if (useGam) {
	
	if (tf.Contains("2023")) f = new TFile(Form("rootfiles/Summer23_L2L3Res/GamHistosFill_data_%s_w8.root",cf),"READ");
	if (tf.Contains("2022")) f = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s_v32.root",cf),"READ");
	if (tf.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/w48_Gam/GamHistosFill_data_%s_w48.root",cf),"READ"); // V9M prompt
	if (tf.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Gam_w65/GamHistosFill_data_%s_w65.root",cf),"READ"); // V2M prompt
      }
      
      // Zmm+jet files
      if (useZmm) {

	if (tf.Contains("2022F") || tf.Contains("2022G")) f = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v78.root",cf),"READ");
	else if (tf.Contains("2022")) f = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v76.root",cf),"READ");
	if (tf.Contains("2023")) f = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v84.root",cf),"READ");
	if (tf.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/v97_Zmm/jme_Zj_%s_Zmm_pileup_69200_V8M_v97.root",cf), "READ"); // V9M prompt
	if (tf.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Zmm_v102/jme_Zj_%s_Zmm_v102_nomu.root",cf), "READ"); // V2M prompt
      }

      // Inclusive jet files
      if (useJet) {

	char cv[256], ce[256];
	TString tf(cf); tf.ReplaceAll("_"," ");
	sscanf(tf.Data(),"%s %s",ce,cv);
	cout << "s="<<cf<<", ce="<<ce<<", cv="<<cv<<endl<<flush;

	if (tf.Contains("2022")) {
	  f = new TFile(Form("rootfiles/Prompt2024/%s_2022/jmenano_data_out_%s_%s.root",cv,ce,cv),"READ");
	}
	if (tf.Contains("2023")) {
	  f = new TFile(Form("rootfiles/Prompt2024/%s_2023/jmenano_data_out_%s_%s.root",cv,ce,cv),"READ");
	}
	if (tf.Contains("2024CDEFGHI_nib")) {
	  f = new TFile(Form("rootfiles/Prompt2024/v121_v2_Jet/jmenano_data_out_%s_Rereco_JME_v121_v2.root",cf),"READ"); // V9M re-reco
	}
	else if (tf.Contains("24C") || tf.Contains("24D") || tf.Contains("24E"))
	  f = new TFile(Form("rootfiles/Prompt2024/v121_v2_Jet/jmenano_data_out_%s_Rereco_JME_v121_v2.root",cf),"READ"); // V9M re-reco
	else if (tf.Contains("2024"))
	  f = new TFile(Form("rootfiles/Prompt2024/v121_v2_Jet/jmenano_data_out_%s_JME_v121_v2.root",cf),"READ"); // V9M prompt
	if (tf.Contains("2025F") || tf.Contains("2025G")) {
	  f = new TFile(Form("rootfiles/Prompt2025/Jet_v147/jmenano_data_out_%s_JME_v147.root",cf),"READ"); // V9M prompt
	}
	else if (tf.Contains("2025")) {
	  f = new TFile(Form("rootfiles/Prompt2025/Jet_v146/jmenano_data_out_%s_JME_v146.root",cf),"READ"); // V9M prompt
	}
      }

      // mW, mTop, ttbar xsec
      if (useTT) {
	f = new TFile("rootfiles/Prompt2024/Emilia_tt_2024.root","READ");
      }
      
      if (!f || f->IsZombie()) cout << "Missing " << cf << endl << flush;
      assert(f && !f->IsZombie());
      curdir->cd();

      TH1D *h(0);
      TProfile *p(0);
      TObject *o(0);
      if (useGam) o = f->Get(Form("runs%s/%s",isEta ? "_high-jeteta" : "",ch));
      if (useZmm) o = f->Get(Form("data/%s",ch));
      if (useJet) o = f->Get(Form("HLT_PFJet500/JetsperRuns/%s",ch));
      if (useTT)  o = f->Get(Form("%s",ch));
      // Detect object type automatically
      if (o) {
	if (o->InheritsFrom("TProfile"))  p = (TProfile*)o;
	else if (o->InheritsFrom("TH1D")) h = (TH1D*)o;
      }

      // Determine the average JES that was applied to each run
      TProfile2D *p2res(0);
      if (useGam) p2res = (TProfile2D*)f->Get("Gamjet2/p2res");
      if (useZmm) p2res = (TProfile2D*)f->Get("data/l2res/p2res");
      if (useJet) p2res = (TProfile2D*)f->Get("HLT_PFJet500/Dijet2/p2res");
      if (useTT)  { p2res = 0; }
      if (!p2res) cout << "Missing p2res in " << cf << endl << flush;
      double jes(1.);
      if (p2res) {
	// Slice out [110,230]x[-1.3,1.3] region to get mean JES here
	//int ipt1 = p2res->GetYaxis()->FindBin(110.);
	//int ipt2 = p2res->GetYaxis()->FindBin(230.);
	double pt1(110.), pt2(230.);
	TString th(ch);
	if (th.Contains("230"))      { pt1 = 230; pt2=500; }
	else if (th.Contains("110")) { pt1 = 110; pt2=230; }
	else if (th.Contains("50"))  { pt1 = 50;  pt2=110; }
	else if (th.Contains("30"))  { pt1 = 30;  pt2=50;  }
	int ipt1 = p2res->GetYaxis()->FindBin(pt1);
	int ipt2 = p2res->GetYaxis()->FindBin(pt2);
	
	TProfile *px = p2res->ProfileX("px",ipt1,ipt2);
	int ieta1 = p2res->GetXaxis()->FindBin(-1.3);
	int ieta2 = p2res->GetXaxis()->FindBin(+1.3);
	px->GetXaxis()->SetRangeUser(0,1.305);//ieta1,ieta2);
	jes = px->GetMean(2);
	delete px;
      }
      
      if (!h && !p) cout << "Missing " << ch << " in " << cf << endl << flush;
      if (h) {
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  int j = hsum->FindBin(h->GetBinCenter(i));
	  double k = 1;//(useJet && tf.Contains("2022") ? 12 : 1);
	  hsum->SetBinContent(j, hsum->GetBinContent(j) + k*h->GetBinContent(i));
	  hsum->SetBinError(j, sqrt(pow(hsum->GetBinError(j),2) + pow(k*h->GetBinError(i),2)));
	} // for i
      } // if h
	       
      if (p) {

	if (debug) cout << "Fill p for file " << cf << endl << flush;
	//rebinProfileCustom(psum, p);

	// Scale out L2L3Res
	TProfile *pjes = (TProfile*)p->Clone(Form("pjes_%s",ch));
	pjes->Scale(jes);
	
	for (int i = 1; i != p->GetNbinsX()+1; ++i) {
	  int run = p->GetBinCenter(i);
	  double cumlum = hcumlum2->GetBinContent(hcumlum2->FindBin(run));
	  //double y = p->GetBinContent(i);
	  //double w = p->GetBinEntries(i);
	  //psum->Fill(cumlum, y, w);

	  // Improved calculation to keep track of RMS
	  int j = psum->FindBin(cumlum);

	  // Add content and entries from the original bin to the new one
	  (*psum)[j] = (*p)[i] + (*psum)[j]; // sumwy
	  (*psum->GetSumw2())[j] = (*p->GetSumw2())[i] + (*psum->GetSumw2())[j]; // sumwy2
	  psum->SetBinEntries(j, p->GetBinEntries(i) + psum->GetBinEntries(j)); // sumw

	  (*psumjes)[j] = (*pjes)[i] + (*psumjes)[j]; // sumwy
	  (*psumjes->GetSumw2())[j] = (*pjes->GetSumw2())[i] + (*psumjes->GetSumw2())[j]; // sumwy2
	  psumjes->SetBinEntries(j, pjes->GetBinEntries(i) + psumjes->GetBinEntries(j)); // sumw
      
	  // Copy (if needed) bin sum of weight square
	  if (p->GetBinSumw2()->fN > i) {
	    //psum->Sumw2();
	    (*psum->GetBinSumw2())[j] = (*p->GetBinSumw2())[i] + (*psum->GetBinSumw2())[j]; // sum2
	  }
	  
	  if (pjes->GetBinSumw2()->fN > i) {
	    psumjes->Sumw2();
	    (*psumjes->GetBinSumw2())[j] = (*pjes->GetBinSumw2())[i] + (*psumjes->GetBinSumw2())[j]; // sum2
	  }
      
	  // Accumulate overall profile entries
	  psum->SetEntries(psum->GetEntries() + p->GetEntries());
	  psumjes->SetEntries(psumjes->GetEntries() + pjes->GetEntries());
	} // for i
      } // if p
    } // for iff
    hsum->Divide(hlum2);

    double xmin = hbins->GetXaxis()->GetXmin();
    double xmax = hbins->GetXaxis()->GetXmax();
    double kf = (psum->Integral()!=0 ? (isPF ? kpf : 1.) : 12.);
    //double ymin = (psum->Integral()!=0 ? (isPF ? -2.5*kpf : -2.5) : -2.5*kf);
    //double ymax = (psum->Integral()!=0 ? (isPF ? +5.5*kpf : +5.5) : +5.5*kf);
    //double ymin = (psum->Integral()!=0 ? (isPF ? -2.5*kpf : -7.0) : -2.5*kf);
    //double ymax = (psum->Integral()!=0 ? (isPF ? +5.5*kpf : +17.0) : +5.5*kf);
    double ymin = (psum->Integral()!=0 ? (isPF ? -2.5*kpf : -11.0) : -2.5*kf);
    double ymax = (psum->Integral()!=0 ? (isPF ? +5.5*kpf : +14.0) : +5.5*kf);
    if (isZmass) {
      ymin = -0.2; ymax = +0.3;
    }
    if (isEta) {
      ymin = -10; ymax = +22;
    }
    TH1D *h1 = tdrHist(Form("h1_%s",ch),"JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
    lumi_136TeV = "Run3, 2022-25";
    extraText = "Private";
    TCanvas *c1 = tdrCanvas(Form("c1_%s",ch),h1,8,11,kRectangular);

    if (addNibLine || addYearLine) {

      TH1D *hbreaks = new TH1D("hbreaks",";Break point;Cum.Lum. (/fb)",ml.size(),0,ml.size());
      
      addLines(ml, hcumlum2, kf, xmin, xmax, ymin, ymax, hbreaks);

      fout->cd();
      hbreaks->GetXaxis()->SetTitleOffset(3.0);
      hbreaks->Write("hbreaks",TObject::kOverwrite);
      curdir->cd();
      delete hbreaks;
    }
    else if (false) { // Draw extra lines etc.

      // Horizontal lines at 0 and +/-1% for reference
      TLine *l = new TLine();
      l->SetLineStyle(kDashed);
      l->SetLineColor(kGray+1);
      l->DrawLine(xmin, 0, xmax, 0);
      l->SetLineStyle(kDotted);
      l->DrawLine(xmin, +1*kf, xmax, +1*kf);
      l->DrawLine(xmin, -1*kf, xmax, -1*kf);

      if (isZmass) {
	l->DrawLine(xmin, +0.03*kf, xmax, +0.03*kf);
	l->DrawLine(xmin, -0.03*kf, xmax, -0.03*kf);
      }

      // Lines at era breaks
      // Define era boundaries as pairs of (start_run, end_run)
      std::vector<std::pair<int, int>> era_boundaries = {
        // 2022 (CDE) (FG)
        {355100, 355793}, {357734, 358219}, {366442, 367079}, 
        // 2024
        {378971, 379411}, {380948, 381383}, {381384, 381943}, {381944, 383779},
        {383780, 385813}, {385814, 386408}, {386409, 386951}
      };
      // List of eras as (start_run, era name)
      std::vector<std::pair<int, string>> eras = {
        // 2022
	{355374, "22"}, {359569, "E"}, {360390, "F"}, {362437, "G"}, // Updated actual first run
        // 2023
	{366727, "23"}, {367770, "Cv4"}, {369927, "D"}, // Updated actual first run
        // 2024
        //{378971, "24"}, /*{379412, "24C"}, {380253, "24D"},*/ {380948, "E"}, /*{381384, "24Ev2"},*/ {381944, "F"},
	{378985, "24"}, {380963, "E"}, {382229, "F"}, // Updated actual first run
	{383811, "G"}, {385836, "H"}, {386478, "I"} // Updated actual first run
      };

      // Additional HCAL breaks
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalRespCorrsTagsRun3
      // https://twiki.cern.ch/twiki/bin/view/CMS/HcalRespCorrsTags2011
      // => https://cms-talk.web.cern.ch/t/fast-track-validation-hlt-prompt-hcal-respcorrs-condition-update-from-hb-time-adjustment/25302/5
      eras.push_back(pair<int,string>(368822,"V1.0")); // (actual run, 23D-)
      eras.push_back(pair<int,string>(382298,"V2.0")); // 24_v2.0 (actual first run)
      eras.push_back(pair<int,string>(383247,"V2.1")); // 24_v2.1 (actual first run)
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
      eras.push_back(pair<int,string>(386401,"V1.0")); // Old tag?
      eras.push_back(pair<int,string>(386732,"VHP")); // HCAL pedestals
      eras.push_back(pair<int,string>(386936,"VHP")); // HCAL pedestals
      
      // Additional ECAL intercalibration updates
      // https://cms-talk.web.cern.ch/c/ppd/alca/108
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-368782-w24/25407
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-384719-w34/46354/10
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386025/57306
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386945/61017
      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379660/39889/5
      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379616/39588/4

      // Some more from:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
      eras.push_back(pair<int,string>(386489,"ICP")); // Pulseshapes
      eras.push_back(pair<int,string>(385728,"ICU")); // PulseShapes+TimeCalib
      eras.push_back(pair<int,string>(385286,"ICT")); // TimeCalib
      eras.push_back(pair<int,string>(384935,"ICN")); // ChannelStatus (noise)
      eras.push_back(pair<int,string>(383756,"ICU")); // PulseShapes+TimeCalib
      eras.push_back(pair<int,string>(382298,"ICF")); // 24_v2.0 (24F DD)
      eras.push_back(pair<int,string>(381379,"ICN")); // ChannelStatus EB+13
      eras.push_back(pair<int,string>(380000,"ICT")); // TimeCalib

      eras.push_back(pair<int,string>(379434,"ICT")); // TimeCalib
      eras.push_back(pair<int,string>(366727,"IC3")); // 2023 start
      eras.push_back(pair<int,string>(362475,"ICP")); // Pedestals
      eras.push_back(pair<int,string>(362012,"ICP")); // Pedestals
      eras.push_back(pair<int,string>(360390,"IC2F")); // 22F start
      eras.push_back(pair<int,string>(359569,"IC2E")); // 22E start
      // ...
      
      // Tracker voltage changes (Martin Delcourt, CMS-PPD-Muon..., 28 Nov 2024)
      eras.push_back(pair<int,string>(380481,"TRK")); // one down
      eras.push_back(pair<int,string>(382770,"TRK")); // one down
      eras.push_back(pair<int,string>(385012,"TRK")); // one down
      eras.push_back(pair<int,string>(385986,"TRK")); // next up => ok!

      // Tracker optical power changes (Martin Delcourt, CMS-PPD-Muon..., 3 Dec 2024)
      eras.push_back(pair<int,string>(380029,"TRP"));
      eras.push_back(pair<int,string>(382209,"TRP"));
      eras.push_back(pair<int,string>(383417,"TRP"));
      eras.push_back(pair<int,string>(383900,"TRP")); // big step in TIDplus
      eras.push_back(pair<int,string>(386278,"TRP")); // ? or 386478?

      // Follow up on BPix_L1 inefficiency (Martin Delcourt, CMS-PPD-Muon, 12 Dec)
      eras.push_back(pair<int,string>(385260,"TRB")); // BPix_L1 inefficiency
      

      // Keep record of breaks for drawTimeStabilityPairs
      TH1D *hbreaks = new TH1D("hbreaks",";Break point;Cum.Lum. (/fb)",eras.size(),0,eras.size());
      
      TLatex *tex = new TLatex();
      tex->SetTextSize(0.045);
      int ks(0), kh(0);
      for (unsigned int i = 0; i != eras.size(); ++i) {

	int run = eras[i].first;
	string s = eras[i].second;
	const char *cn = s.c_str();
	TString t = TString(cn);

	int j = hcumlum2->FindBin(run);
	if (hcumlum2->GetBinLowEdge(j)!=run) cout << "run="<<run<<", bin edges="<<hcumlum2->GetBinLowEdge(j)<<","<<hcumlum2->GetBinLowEdge(j+1)<<" (lum="<<hcumlum2->GetBinContent(j)<<")"<<endl<<flush;
	double cumlum = hcumlum2->GetBinContent(j);

	l->SetLineColor(kGray);
	if (t.Contains("V")) l->SetLineColor(kRed);
	if (t.Contains("IC")) l->SetLineColor(kBlue);
	if (t.Contains("EP")) l->SetLineColor(kBlue);
	if (t.Contains("XX")) l->SetLineColor(kMagenta+1);
	if (t.Contains("DJ")) l->SetLineColor(kOrange+1);
	if (t.Contains("TR")) l->SetLineColor(kOrange+1);
	l->SetLineStyle(kSolid);
	if (addExtraLine || l->GetLineColor()==kGray) {
	  if (l->GetLineColor()==kGray) l->SetLineColorAlpha(kGray,0.3);
	  l->DrawLine(cumlum,ymin,cumlum,ymax);
	}
	if (s=="22"||s=="23"||s=="24") ks = 0;
	if (t.Contains("V") || t.Contains("IC") || t.Contains("XX") || t.Contains("EP") || t.Contains("TR")) {
	  if (addExtraText) {
	    if (isZmass) tex->DrawLatex(cumlum+1,(-0.12-0.03*(kh++%3))*kf,cn);
	    else if (isEta) tex->DrawLatex(cumlum+1,(-3.2-0.8*(kh++%3))*kf,cn);
	    else tex->DrawLatex(cumlum+1,(-1.6-0.4*(kh++%3))*kf,cn);
	  }
	}
	else {
	  if (isZmass) tex->DrawLatex(cumlum+1,(+0.15-0.03*(ks++))*kf,cn);
	  else if (isEta) tex->DrawLatex(cumlum+1,(+7.0-0.8*(ks++))*kf,cn);
	  else tex->DrawLatex(cumlum+1,(+3.5-0.4*(ks++))*kf,cn);
	}
	hbreaks->SetBinContent(i+1,cumlum);
	hbreaks->GetXaxis()->SetBinLabel(i+1,cn);
      }
      
      fout->cd();
      hbreaks->Write("hbreaks",TObject::kOverwrite);
      curdir->cd();
      delete hbreaks;
    }
      

    TF1 *f1 = new TF1("f1","[0]",0,1e6); f1->SetParameter(0,0);
    TF1 *f2 = new TF1("f2","[0]",0,1e6); f2->SetParameter(0,0);
    if (hsum->Integral()!=0) {
      //f1->SetRange(382298,1e6); // V2.0 onwards
      f1->SetRange(379416, 392112+0.5); // 2024CDEFGH
      hsum->Fit(f1,"QRN");
      hsum->Scale(1./f1->GetParameter(0));
      
      // Turn Xsec into Xsec-1 (%) and map to cumlum2
      TH1D  *h = (TH1D*)hbins->Clone(Form("h_%s",ch)); h->Reset();
      for (int i = 1; i != hsum->GetNbinsX()+1; ++i) {
	int run = hsum->GetBinCenter(i);
	double cumlum = hcumlum2->GetBinContent(hcumlum2->FindBin(run));
	int j = hbins->FindBin(cumlum);
	if (hsum->GetBinContent(i)!=0) {

	  // Patch for jets
	  double k(1);
	  if (useJet) {
	    if      (run<=362760) k = 1.7*34.8/109.; // 2022
	    else if (run<=370790) k = 1.7*28.4/109.; // 2023
	    else if (run<=392112) k = 109./109.; // 2024
	    else if (run<=400000) k = 109./109.; // 2025
	  }
	  h->SetBinContent(j, (k*hsum->GetBinContent(i)-1)*100.);
	  h->SetBinError(j, k*hsum->GetBinError(i)*100.);
	}
      }
      
      h1->SetYTitle("Rel.Xsec-1 (%)");

      tdrDraw(h,"PE",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.6);

      fout->cd();
      h->Write(Form("%s",ch),TObject::kOverwrite);
      curdir->cd();
    } // hsum
    
    if (psum->Integral()!=0) {
      //f1->SetRange(82,1e6); // V2.0 onwards => wrong /fb
      //f1->SetRange(90.7,1e6); // V2.0 onwards
      //f2->SetRange(90.7,1e6); // V2.0 onwards

      int run1 = 379416, run2 = 392112; // 2024CDEFGH
      int j1 = hcumlum2->FindBin(run1);
      double cumlum1 = hcumlum2->GetBinContent(j1);
      int j2 = hcumlum2->FindBin(run2);
      double cumlum2 = hcumlum2->GetBinContent(j2);
      f1->SetRange(cumlum1,cumlum2);
      f2->SetRange(cumlum1,cumlum2);
      
      // Normalize average to unity to focus on time dependence
      psum->Fit(f1,"QRN");
      TH1D *h = psum->ProjectionX(Form("h_%s",ch));
      h->Scale(1./f1->GetParameter(0));
      
      psumjes->Fit(f2,"QRN");
      TH1D *hjes = psumjes->ProjectionX(Form("hjes_%s",ch));
      hjes->Scale(1./f2->GetParameter(0));

      // Turn JES into JES-1 (%)
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	if (h->GetBinContent(i)!=0) {
	  h->SetBinContent(i, (h->GetBinContent(i)-1)*100.);
	  h->SetBinError(i, h->GetBinError(i)*100.);
	  
	  hjes->SetBinContent(i, (hjes->GetBinContent(i)-1)*100.);
	  hjes->SetBinError(i, hjes->GetBinError(i)*100.);
	}
      }

      if (!isZmass)
	tdrDraw(hjes,"PE",kOpenSquare,kRed,kSolid,-1,kNone,0,0.6);
      tdrDraw(h,"PE",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.6);
      
      TLegend *leg = tdrLeg(0.65,0.85-0.05*3,0.90,0.85);
      //leg->SetHeader("#gamma+jet MPF 110EB");
      leg->SetHeader(mhead[ch]);
      if (isZmass) {
	leg->AddEntry(h,"m_{Z} / #LTm_{Z}#GT");
      }
      else {
	leg->AddEntry(hjes,"Before JEC");
	leg->AddEntry(h,"After JEC");
      }
	
      fout->cd();
      if (!isZmass)
	hjes->Write(Form("%s_jes",ch),TObject::kOverwrite);
      h->Write(Form("%s",ch),TObject::kOverwrite);
      curdir->cd();
    } // psum

    hscales->SetBinContent(2*ih+1, f1->GetParameter(0));
    hscales->SetBinError(2*ih+1, f1->GetParError(0));
    hscales->GetXaxis()->SetBinLabel(2*ih+1,ch);
    hscales->SetBinContent(2*ih+2, f2->GetParameter(0));
    hscales->SetBinError(2*ih+2, f2->GetParError(0));
    hscales->GetXaxis()->SetBinLabel(2*ih+2,Form("%s_jes",ch));
    
    c1->RedrawAxis();
    c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStability_%s.pdf",ch));
  } // for ih

  fout->cd();
  hscales->Write("hscales",TObject::kOverwrite);
  curdir->cd();
  
  cout << "Results stored in rootfiles/drawTimeStability.root" << endl << flush;
  fout->Close();
} // drawTimeStability



void rebinProfileCustom(TProfile* p, TProfile* p_new) {
  // Iterate through the bins of the original TProfile
  for (int ix = 1; ix <= p->GetNbinsX(); ix++) {
                
    // Original bin index
    int bin_orig = ix;
      
    // Get bin centers from original profile
    double x_center = p->GetXaxis()->GetBinCenter(ix);

    // Find corresponding bin in the new profile
    int new_ix = p_new->GetXaxis()->FindBin(x_center);
      
    // Get the new bin index
    int bin_new = new_ix;
      
    // profile keeps track of sumy, sumwy, sumwy2, sumw2
    // sumw=fArray, sumwy=fBinEntries.fArray, 
    // sumwy2 = fBinSumw2.fArray, sumw2 = fSum2.fArray
    // GetBinContent = sumwy/sumw
    // https://root-forum.cern.ch/t/copy-entries-of-tprofile/11828
      
    // Add content and entries from the original bin to the new one
    (*p_new)[bin_new] = (*p)[bin_orig] + (*p_new)[bin_new]; // sumwy
    (*p_new->GetSumw2())[bin_new] = (*p->GetSumw2())[bin_orig] +
      (*p_new->GetSumw2())[bin_new]; // sumwy2
    p_new->SetBinEntries(bin_new, p->GetBinEntries(bin_orig) + 
                          p_new->GetBinEntries(bin_new)); // sumw
      
    // Copy (if needed) bin sum of weight square
    if (p->GetBinSumw2()->fN > bin_orig) {
      p_new->Sumw2();
      (*p_new->GetBinSumw2())[bin_new] = (*p->GetBinSumw2())[bin_orig] +
        (*p_new->GetBinSumw2())[bin_new]; // sum2
    }
      
    // Accumulate overall profile entries
    p_new->SetEntries(p_new->GetEntries() + p->GetEntries());
  } // for ix
} // rebinProfileCustom


//#include <bits/stdc++.h>
using namespace std;

map<string,pair<int,int>> BuildRanges(const string& fname){

  ifstream in(fname); string L; map<string,pair<int,int>> M;
  auto upd = [&](const string& k, int a, int b){
    auto it=M.find(k);
    if(it==M.end()) M[k]={a,b};
    else { if(a<it->second.first) it->second.first=a; if(b>it->second.second) it->second.second=b; }
  };
  auto trim=[&](string &s){ size_t i=s.find_first_not_of(" \t"), j=s.find_last_not_of(" \t"); s=(i==string::npos)?"":s.substr(i,j-i+1); };

  getline(in,L); // remove header line
  while(getline(in,L)){
    size_t lb=L.find('['), rb=L.find(']'); if(lb==string::npos||rb==string::npos) continue;
    int a=0,b=0; sscanf(L.c_str()+lb,"[%d ,%d",&a,&b);
    size_t p1=L.find('|',rb), p2=(p1==string::npos)?string::npos:L.find('|',p1+1);
    if(p1==string::npos||p2==string::npos) continue;
    string tag=L.substr(p1+1,p2-p1-1); trim(tag);

    string nib=tag; size_t q=nib.rfind("-fib"); if(q!=string::npos) nib.resize(q);
    string era=nib; q=nib.rfind("-nib");     if(q!=string::npos) era.resize(q);
    string year=era; year.resize(4);

    upd(nib,a,b);    // e.g. "2024G-nib2"
    upd(era,a,b);   // e.g. "2024G"
    upd(year,a,b);   // e.g. "2024"
  }
  return M;
}


void addLines(map<string,pair<int,int>> &ml, TH1D *hcumlum2, double kf,
	      double xmin, double xmax, double ymin, double ymax,
	      TH1D *hbreaks) {

  //map<string,pair<int,int> > ml = BuildRanges("rootfiles/brilcalc/fibs.txt");
  typedef map<string,pair<int,int> >::const_iterator IT; int i(0);
  double prevlum(0);
  for (IT it = ml.begin(); it != ml.end(); ++it) {
    string tag = it->first;
    int run = it->second.first;
    bool isYear = (tag.size()==4);
    
    int j = hcumlum2->FindBin(run);
    double cumlum = hcumlum2->GetBinContent(j);
    
    // Horizontal lines at 0 and +/-1% for reference
    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->SetLineColor(kGray+1);
    l->DrawLine(xmin, 0, xmax, 0);
    l->SetLineStyle(kDotted);
    l->DrawLine(xmin, +1*kf, xmax, +1*kf);
    l->DrawLine(xmin, -1*kf, xmax, -1*kf);
    
    TLatex *tex = new TLatex();
    tex->SetTextSize(0.015);
    tex->SetTextAngle(270);
    tex->SetTextAlign(31);//33);//12);
    
    l->SetLineStyle(kSolid);
    if (hbreaks) {
      TString t(tag.c_str());
      hbreaks->SetBinContent(++i,cumlum);
      if (!t.Contains("-nib1")) hbreaks->GetXaxis()->SetBinLabel(i,tag.c_str());
    }
    // Avoid label overlaps
    if (cumlum-prevlum>3.0 || isYear) {
      //if (hbreaks) hbreaks->GetXaxis()->SetBinLabel(i,tag.c_str());
      tex->SetTextColorAlpha(kGray,isYear ? 1.0 : 0.3);
      l->SetLineColorAlpha(kGray,isYear ? 1.0 : 0.1);
      l->DrawLine(cumlum,ymin,cumlum,ymax);
      tex->DrawLatex(cumlum+0.01,ymin,tag.c_str());
      prevlum = cumlum;
    }
    else {
      l->SetLineColorAlpha(kGray,0.075);
      l->DrawLine(cumlum,ymin+3.0,cumlum,ymax);
    }
  }
} // void addLines
