// Purpose: Draw Z+jet, gamma+jet and inclusive jet xsec and MPF vs time
//          Use breakup made by rootfiles/brilcalc/clusterRuns.C
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include <iostream>

#include "../tdrstyle_mod22.C"

bool debug = false;
bool addExtraText = false;
bool addExtraLine = false;
void rebinProfileCustom(TProfile* p, TProfile* p_new);

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
    // //"2022C_v39","2022D_v39","2022E_v39","2022F_v39","2022G_v39",
    //"2022C_v113","2022D_v113","2022E_v113","2022F_v113","2022G_v113",
    //"2023Cv123_v39","2023Cv4_v39","2023D_v39",
    // //"2024B_v110","2024C_v110","2024D_v110","2024E_v110",
    // //"2024F_v110","2024G_v110","2024H_v111","2024I_v111",
    //"2024B_v112","2024C_v112","2024D_v112","2024E_v112",
    //"2024F_v112","2024G_v112","2024H_v112","2024I_v112",

    //"2022_v113","2023_v113","2024_v113",
    "2024B_nib1","2024C_nib1","2024D_nib1","2024Ev1_nib1","2024Ev2_nib1",
    "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    "2024H_nib1","2024I_nib1",

    // gamma+jet
    "GamJet",
    //"2022C","2022D","2022E","2022F","2022G",
    //"2023Cv123","2023Cv4","2023D",
    //"2024B","2024C","2024D","2024Ev1","2024Ev2","2024F","2024G","2024H","2024Iv1","2024Iv2",
    "2024B_nib1","2024C_nib1","2024D_nib1","2024Ev1_nib1","2024Ev2_nib1",
    "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    "2024H_nib1","2024I_nib1",

    //"2022CDE_v32","2022FG_v32",
    //"2023Cv123_w8","2023Cv4_w8","2023D_w8",
    //"2024BCD_w39","2024E_w39","2024F_w39","2024G_w39","2024H_w40","2024I_w40",


    // Z+jet
    "ZmmJet",
    //"2022CD","2022E",//"2022FG",
    //"2022F","2022G",
    //"2023C123","2023C4","2023D",
    //"2024BCD","2024E","2024F","2024G","2024H","2024I"
    "2024B_nib1","2024C_nib1","2024D_nib1","2024Ev1_nib1","2024Ev2_nib1",
    "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
    "2024H_nib1","2024I_nib1"

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
    //"mz_run_zpt0","mz_run_zpt30","mz_run_zpt50","mz_run_zpt110",
    "nz_run_zpt0","nz_run_zpt30","nz_run_zpt50","nz_run_zpt110",

    "mpf_run_zpt30","mpf_run_zpt110","mpf_run_zpt50",
    "db_run_zpt30","db_run_zpt110","db_run_zpt50",

    /*
    "chf_run_zpt30","nef_run_zpt30","nhf_run_zpt30",
    "chf_run_zpt50","nef_run_zpt50","nhf_run_zpt50",
    "chf_run_zpt110","nef_run_zpt110","nhf_run_zpt110",
    */
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
  mhead["pr50chf"] = "#gamma+jet CHF 50EB";
  mhead["pr50nhf"] = "#gamma+jet NHF 50EB";
  mhead["pr50nef"] = "#gamma+jet NEF 50EB";
  mhead["pr110chf"] = "#gamma+jet CHF 110EB";
  mhead["pr110nhf"] = "#gamma+jet NHF 110EB";
  mhead["pr110nef"] = "#gamma+jet NEF 110EB";
  mhead["mpf_run_zpt30"] = "Z(#mu#mu)+jet MPF 30";
  mhead["mpf_run_zpt50"] = "Z(#mu#mu)+jet MPF 50";
  mhead["mpf_run_zpt110"] = "Z(#mu#mu)+jet MPF 110";
  mhead["chf_run_zpt30"] = "Z(#mu#mu)+jet CHF 30";
  mhead["nef_run_zpt30"] = "Z(#mu#mu)+jet NEF 30";
  mhead["nhf_run_zpt30"] = "Z(#mu#mu)+jet NHF 30";
  

  // Store normalizations to histogra
  TH1D *hscales = new TH1D("hscales",";Observable;Scale",2*nh,0,2*nh);
  
  // Loop over files to retrieve stuff
  bool useGam(false), useZmm(false), useJet(false);
  for (int ih = 0; ih != nh; ++ih) {

    string sh = vh[ih];
    const char *ch = sh.c_str();
    TString th(ch);
    cout << "Analyzing " << ch << endl << flush;
    if (sh=="GamJet") { useGam = true; useZmm = false; useJet=false; continue; }
    if (sh=="ZmmJet") { useGam = false; useZmm = true; useJet=false; continue; }
    if (sh=="DiJet")  { useGam = false; useZmm = false; useJet=true; continue; }
    bool isPF = (th.Contains("nhf") || th.Contains("nef") || th.Contains("chf"));
    bool isZmass = (th.Contains("mz_run"));
    bool isEta = (th.Contains("eta"));
    double kpf(1);
    if (th.Contains("nhf")) kpf = 20;  // 5%
    if (th.Contains("nef")) kpf = 4;   // 25%
    if (th.Contains("chf")) kpf = 1.5; // 65%
    
    TH1D *hsum = (TH1D*)hlum2->Clone(Form("hsum_%s",ch)); hsum->Reset();
    //double vx[hsum->GetNbinsX()+1];
    //for (int i = 1; i != hsum->GetNbinsX()+2; ++i) vx[i-1] = hsum->GetBinLowEdge(i);
    //TProfile *psum = new TProfile(Form("psum_%s",ch),"",hsum->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());
    double vx[hbins->GetNbinsX()+1];
    for (int i = 1; i != hbins->GetNbinsX()+2; ++i) vx[i-1] = hbins->GetBinLowEdge(i);
    TProfile *psum = new TProfile(Form("psum_%s",ch),"",hbins->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());
    TProfile *psumjes = new TProfile(Form("psumjes_%s",ch),"",hbins->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());

    // Keep track of JES (1=L2L3Res) for gamma+jet MPF
    //TH1D *hjes = new TH1D(Form("hjes_%s",ch),";JES;Cum.Lum(/fb)",hbins->GetNbinsX(),&vx[0]);

    bool isGam(false), isZmm(false), isJet(false);
    for (int iff = 0; iff != nf; ++iff) {

      string sf = vf[iff];
      const char *cf = sf.c_str();
      TString tf(cf);
      
      if (sf=="GamJet") { isGam=true; isZmm=false; isJet=false; continue; }
      if (sf=="ZmmJet") { isGam=false; isZmm=true; isJet=false; continue; }
      if (sf=="DiJet")  { isGam=false; isZmm=false; isJet=true; continue; }
      if (useGam && !isGam) continue;
      if (useZmm && !isZmm) continue;
      if (useJet && !isJet) continue;
      
      TFile *f(0);
      // Photon+jet files
      if (useGam) {
	/*
	if (tf.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s.root",cf),"READ");
	if (tf.Contains("2023")) f = new TFile(Form("rootfiles/Summer23_L2L3Res/GamHistosFill_data_%s.root",cf),"READ");
	if (tf.Contains("2022")) f = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s.root",cf),"READ");
	*/
	//f = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w42.root",cf),"READ");
	f = new TFile(Form("rootfiles/Prompt2024/v45_Gam/GamHistosFill_data_%s_w45.root",cf),"READ");
      }
      // Zmm+jet files
      if (useZmm) {
	/*
	if (tf.Contains("2022FG")) f = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v78.root",cf),"READ");
	else if (tf.Contains("2022"))   f = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v76.root",cf),"READ");
	if (tf.Contains("2023"))   f = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v84.root",cf),"READ");
	if (tf.Contains("2024"))   f = new TFile(Form("rootfiles/Prompt2024/v88/jme_bplusZ_%s_Zmm_sync_v88.root",cf),"READ");
	*/
	
	//if (tf.Contains("2024F")) f = new TFile(Form("rootfiles/Prompt2024/v88/jme_bplusZ_%s_Zmm_sync_v88.root",cf),"READ");
	/*
	if (tf.Contains("2024BC") || tf.Contains("2024F") || tf.Contains("2024G") || tf.Contains("2024H")) f = new TFile(Form("rootfiles/Prompt2024/v91/jme_bplusZ_%s_Zmm_sync_v91d.root",cf),"READ");
	else
	  f = new TFile(Form("rootfiles/Prompt2024/v91/jme_bplusZ_%s_Zmm_sync_v91.root",cf),"READ");
	*/
	f = new TFile(Form("rootfiles/Prompt2024/v94_Zmm/jme_bplusZ_%s_Zmm_v94_Summer24.root",cf),"READ");
      }
      if (useJet) {
	char cv[256], ce[256];
	TString tf(cf); tf.ReplaceAll("_"," ");
	sscanf(tf.Data(),"%s %s",ce,cv);
	cout << "s="<<cf<<", ce="<<ce<<", cv="<<cv<<endl<<flush;
	/*
	if (tf.Contains("2024")) {
	  f = new TFile(Form("rootfiles/Prompt2024/%s_2024/jmenano_data_out_%s_JME_%s_2024.root",cv,ce,cv),"READ");
	}
	if (tf.Contains("2023")) {
	  f = new TFile(Form("rootfiles/Summer23_L2L3ResJERSF/%s_2023_etabin_SFv2/jmenano_data_out_%s_JME_%s_2023_etabin_SFv2.root",cv,ce,cv),"READ");
	}
	if (tf.Contains("2022")) {
	  //f = new TFile(Form("rootfiles/Summer22_L2L3ResJERSF/%s_2022_etabin_SFv/jmenano_data_out_%s_JME_%s_2022_etabin_SFv.root",cv,ce,cv),"READ");
	  f = new TFile(Form("rootfiles/Prompt2024/%s_2022/jmenano_data_out_%s_nib1_JME_%s_2022.root",cv,ce,cv),"READ");
	}
	*/
	//f = new TFile(Form("rootfiles/Prompt2024/%s_%s/jmenano_data_out_%s_%s.root",cv,ce,ce,cv));
	f = new TFile(Form("rootfiles/Prompt2024/v116_Jet/jmenano_data_out_%s_JME_v116.root",cf),"READ");
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
	px->GetXaxis()->SetRangeUser(ieta1,ieta2);
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
	    psum->Sumw2();
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
    //double kf = (psum->Integral()!=0 ? (isPF ? kpf : 1.) : 20.);
    double kf = (psum->Integral()!=0 ? (isPF ? kpf : 1.) : 12.);
    double ymin = (psum->Integral()!=0 ? (isPF ? -2.5*kpf : -2.5) : -2.5*kf);
    double ymax = (psum->Integral()!=0 ? (isPF ? +5.5*kpf : +5.5) : +5.5*kf);
    if (isZmass) {
      ymin = -0.2; ymax = +0.3;
    }
    //if (useJet) {
    //ymin = -200; ymax=+400;
    //}
    if (isEta) {
      //ymin = -5; ymax = +11;
      //ymin = -7.5; ymax = +16.5;
      ymin = -10; ymax = +22;
    }
    //TH1D *h1 = tdrHist(Form("h1_%s",ch),"JES",0.975,1.025,"Cumulative luminosity (fb^{-1})",0,hbins->GetXaxis()->GetXmax());
    TH1D *h1 = tdrHist(Form("h1_%s",ch),"JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
    lumi_136TeV = "Run3, 2022-24";
    extraText = "Private";
    TCanvas *c1 = tdrCanvas(Form("c1_%s",ch),h1,8,11,kRectangular);

    if (true) { // Draw extra lines etc.

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
        {355100, 355793}, /*{355794, 357486}, {357487, 357733},*/ {357734, 358219}, /*{358220, 359021}, {359022, 360331},*/
        //{360332, 362180}, {362181, 362349}, {362350, 362760},
        // 2023
        {366442, 367079}, //{367080, 367515}, {367516, 367620}, {367621, 367763}, {367765, 369802}, {369803, 370602}, {370603, 370790},
        // 2024
        {378971, 379411}, /*{379412, 380252}, {380253, 380947},*/ {380948, 381383}, {381384, 381943}, {381944, 383779},
        {383780, 385813}, {385814, 386408}, {386409, 386951}
      };
      // List of eras as (start_run, era name)
      std::vector<std::pair<int, string>> eras = {
        // 2022
        //{355100, "22"},/* {355794, "22C"}, {357487, "22D"}, {357734, "22Dv2"}, {358220, "22Dv3"},*/ {359022, "E"},
        //{360332, "F"}, /*{362181, "22HI"}, {362350, "G"},*/
	{355374, "22"}, {359569, "E"}, {360390, "F"}, /*{362181, "22HI"},*/ {362437, "G"}, // Updated actual first run
        // 2023
        //{366442, "23"}, /*{367080, "23C"}, {367516, "23Cv2"}, {367621, "23Cv3"},*/ {367765, "Cv4"}, {369803, "D"}, /*{370603, "23Dv2"},*/
	{366727, "23"}, {367770, "Cv4"}, {369927, "D"}, // Updated actual first run
        // 2024
        //{378971, "24"}, /*{379412, "24C"}, {380253, "24D"},*/ {380948, "E"}, /*{381384, "24Ev2"},*/ {381944, "F"},
	{378985, "24"}, {380963, "E"}, {382229, "F"}, // Updated actual first run
        //{383780, "G"}, {385814, "H"}, {386409, "I"}
	{383811, "G"}, {385836, "H"}, {386478, "I"} // Updated actual first run
      };

      // Additional HCAL breaks
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalRespCorrsTagsRun3
      // https://twiki.cern.ch/twiki/bin/view/CMS/HcalRespCorrsTags2011
      //eras.push_back(pair<int,string>(383195,"24_v2.0")); // HLT
      //eras.push_back(pair<int,string>(383219,"24_v2.1")); // HLT
      //eras.push_back(pair<int,string>(386401,"24_v3.0")); // HLT, eraH

      // => https://cms-talk.web.cern.ch/t/fast-track-validation-hlt-prompt-hcal-respcorrs-condition-update-from-hb-time-adjustment/25302/5
      //eras.push_back(pair<int,string>(368775,"V1.0")); // 368765, 368822 =>
      eras.push_back(pair<int,string>(368822,"V1.0")); // (actual run, 23D-)
      
      //eras.push_back(pair<int,string>(367765,"V1.0")); // 23_v1.0 => 23Cv4
      //eras.push_back(pair<int,string>(380637,"V1.0")); // Something in mid-24C?
      //eras.push_back(pair<int,string>(380852,"V1.0")); // 24E
      //eras.push_back(pair<int,string>(382287,"V2.0")); // 24_v2.0 =>
      eras.push_back(pair<int,string>(382298,"V2.0")); // 24_v2.0 (actual first run)
      //eras.push_back(pair<int,string>(383219,"V2.1")); // 24_v2.1 =>
      eras.push_back(pair<int,string>(383247,"V2.1")); // 24_v2.1 (actual first run)
      //eras.push_back(pair<int,string>(386401,"V3.0")); // 24_v3.0 => 24I =>
      //eras.push_back(pair<int,string>(386478,"V3.0")); // 24_v3.0 / 24I (actual first run)
      //eras.push_back(pair<int,string>(383811,"VG")); // 24G tes

      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
      eras.push_back(pair<int,string>(386401,"V1.0")); // Old tag?
      eras.push_back(pair<int,string>(386732,"VHP")); // HCAL pedestals
      eras.push_back(pair<int,string>(386936,"VHP")); // HCAL pedestals
      
      // Additional ECAL intercalibration updates
      // https://cms-talk.web.cern.ch/c/ppd/alca/108
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-368782-w24/25407
      //era_boundaries.push_back(pair<int,int>(368919,368919)); // 23D-, 368823, 369927 =>
      //eras.push_back(pair<int,string>(369927,"EP")); // =? 23D (actual run)
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-384719-w34/46354/10
      //eras.push_back(pair<int,string>(384719,"IC")); // mid-24G! =>
      //eras.push_back(pair<int,string>(384933,"IC1")); // mid-24G! (actual first run) 384644, 384933
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386025/57306
      //eras.push_back(pair<int,string>(386025,"IC")); // also actual run? or deployed later?
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386945/61017
      // eras.push_back(pair<int,string>(386945,"IC")); // post 24I

      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379660/39889/5
      //eras.push_back(pair<int,string>(380115,"IC")); // mid-24C (too late)
      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379616/39588/4
      //eras.push_back(pair<int,string>(379956,"IC")); // mid-24C 379866, 379984 =>
      //eras.push_back(pair<int,string>(379984,"IC1")); // (actual run)

      // Some more from:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/AlCaTSGConditionsUpdate
      //eras.push_back(pair<int,string>(387696,"IC"));
      //eras.push_back(pair<int,string>(387048,"IC"));
      //eras.push_back(pair<int,string>(386663,"ICT")); // TimeCalib
      //eras.push_back(pair<int,string>(386661,"ICS")); // PulseShapes
      //eras.push_back(pair<int,string>(386455,"IC"));
      eras.push_back(pair<int,string>(386489,"ICP")); // Pulseshapes
      //eras.push_back(pair<int,string>(386357,"ICL")); // TPGLinearization??
      //eras.push_back(pair<int,string>(386319,"IC"));
      //eras.push_back(pair<int,string>(386186,"IC"));
      eras.push_back(pair<int,string>(385728,"ICU")); // PulseShapes+TimeCalib
      //eras.push_back(pair<int,string>(385620,"ICC")); // per crystal
      eras.push_back(pair<int,string>(385286,"ICT")); // TimeCalib
      //eras.push_back(pair<int,string>(385154,"ICU")); // PulseShapes+TimeCalib
      //eras.push_back(pair<int,string>(384982,"IC")); // before
      eras.push_back(pair<int,string>(384935,"ICN")); // ChannelStatus (noise)
      //eras.push_back(pair<int,string>(384756,"ICT")); // TimeCalib FEdelay,aft
      //eras.push_back(pair<int,string>(384719,"ICP")); // pedestals
      // 384719 -> 384933 (mid-24G)
      //eras.push_back(pair<int,string>(384578,"IC"));
      //eras.push_back(pair<int,string>(384413,"ICU")); // PulseShapes+TimeCalib
      //eras.push_back(pair<int,string>(384332,"IC"));
      //eras.push_back(pair<int,string>(384285,"IC"));
      //eras.push_back(pair<int,string>(384237,"IC"));
      //eras.push_back(pair<int,string>(384052,"IC2"));
      eras.push_back(pair<int,string>(383756,"ICU")); // PulseShapes+TimeCalib
      //eras.push_back(pair<int,string>(383227,"IC"));
      //eras.push_back(pair<int,string>(383155,"IC2"));
      //eras.push_back(pair<int,string>(382921,"IC2"));
      //eras.push_back(pair<int,string>(382777,"IC"));
      //eras.push_back(pair<int,string>(382712,"IC"));
      eras.push_back(pair<int,string>(382298,"ICF")); // 24_v2.0 (24F DD)
      //eras.push_back(pair<int,string>(381640,"IC"));
      eras.push_back(pair<int,string>(381379,"ICN")); // ChannelStatus EB+13
      //eras.push_back(pair<int,string>(381208,"IC")); // BIG, later
      //eras.push_back(pair<int,string>(380992,"IC"));
      //eras.push_back(pair<int,string>(380049,"ICN")); // ChannelStatus EB-03
      eras.push_back(pair<int,string>(380000,"ICT")); // TimeCalib
      //eras.push_back(pair<int,string>(379956,"IC")); // BIG or MID?
      // 379956 -> 379984
      //eras.push_back(pair<int,string>(379693,"ICN")); // ChannelStatus EB+16
      //eras.push_back(pair<int,string>(379661,"ICN")); // ChannelStatus EB+16,bef
      eras.push_back(pair<int,string>(379434,"ICT")); // TimeCalib
      //eras.push_back(pair<int,string>(379367,"ICT")); // TimeCalib CC,aft
      eras.push_back(pair<int,string>(366727,"IC3")); // 2023 start
      //eras.push_back(pair<int,string>(362523,"ICT")); // TimeCalib
      //eras.push_back(pair<int,string>(362181,"ICHI")); // HI
      eras.push_back(pair<int,string>(362475,"ICP")); // Pedestals
      eras.push_back(pair<int,string>(362012,"ICP")); // Pedestals
      eras.push_back(pair<int,string>(360390,"IC2F")); // 22F start
      eras.push_back(pair<int,string>(359569,"IC2E")); // 22E start
      // ...
      
      //eras.push_back(pair<int,string>(368824,"XX")); // 23D- mystery run
      //eras.push_back(pair<int,string>(368765,"XX")); // 23D- mystery run
      //eras.push_back(pair<int,string>(380030,"XX")); // mid-24C mystery runx

      // Dijet team
      //Runs 382298 - 383246
      //Runs 383247 – 384932
      //Runs 384933 – End of Era G
      //eras.push_back(pair<int,string>(382298,"DJ")); // HCAL 24_v2.0
      //eras.push_back(pair<int,string>(383247,"DJ")); // HCAL 24_v2.1
      //eras.push_back(pair<int,string>(384933,"DJ")); // IC

      // Tracker voltage changes (Martin Delcourt, CMS-PPD-Muon..., 28 Nov 2024)
      //eras.push_back(pair<int,string>(380504,"TRK"));
      //eras.push_back(pair<int,string>(380517,"TRK")); // next up => nah
      eras.push_back(pair<int,string>(380481,"TRK")); // one down
      //eras.push_back(pair<int,string>(382772,"TRK"));
      //eras.push_back(pair<int,string>(382834,"TRK")); // next up => nah
      eras.push_back(pair<int,string>(382770,"TRK")); // one down
      //eras.push_back(pair<int,string>(385063,"TRK"));
      //eras.push_back(pair<int,string>(385127,"TRK")); // next up => not ok
      eras.push_back(pair<int,string>(385012,"TRK")); // one down
      //eras.push_back(pair<int,string>(385959,"TRK"));
      eras.push_back(pair<int,string>(385986,"TRK")); // next up => ok!
      //eras.push_back(pair<int,string>(386968,"TRK"));
      //eras.push_back(pair<int,string>(387017,"TRK")); // next up (after I)

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
	//if (hcumlum2->GetBinLowEdge(j)<run) ++j;
	if (hcumlum2->GetBinLowEdge(j)!=run) cout << "run="<<run<<", bin edges="<<hcumlum2->GetBinLowEdge(j)<<","<<hcumlum2->GetBinLowEdge(j+1)<<" (lum="<<hcumlum2->GetBinContent(j)<<")"<<endl<<flush;
	double cumlum = hcumlum2->GetBinContent(j);
	//int k = hcumlum->FindBin(run)+1;
	//if (hcumlum->GetBinLowEdge(k)<run) ++k;
	//double cumlum1 = hcumlum->GetBinContent(k);

	l->SetLineColor(kGray);
	if (t.Contains("V")) l->SetLineColor(kRed);
	if (t.Contains("IC")) l->SetLineColor(kBlue);
	if (t.Contains("EP")) l->SetLineColor(kBlue);
	if (t.Contains("XX")) l->SetLineColor(kMagenta+1);
	if (t.Contains("DJ")) l->SetLineColor(kOrange+1);
	if (t.Contains("TR")) l->SetLineColor(kOrange+1);
	l->SetLineStyle(kSolid);
	//if (t.Contains("V")) l->DrawLine(cumlum1,ymin,cumlum1,ymax);
	//else
	if (addExtraLine || l->GetLineColor()==kGray)
	  l->DrawLine(cumlum,ymin,cumlum,ymax);
	if (s=="22"||s=="23"||s=="24") ks = 0;
	//tex->DrawLatex(cumlum+1,-1.2-0.2*(ks++),cn);
	//if (t.Contains("V")) tex->DrawLatex(cumlum1+1,-1.6-0.4*(kh++),cn);
	if (t.Contains("V") || t.Contains("IC") || t.Contains("XX") || t.Contains("EP") || t.Contains("TR")) {
	  if (addExtraText) {
	    if (isZmass) tex->DrawLatex(cumlum+1,(-0.12-0.03*(kh++%3))*kf,cn);
	    else if (isEta) tex->DrawLatex(cumlum+1,(-3.2-0.8*(kh++%3))*kf,cn);
	  //else if (isEta) tex->DrawLatex(cumlum+1,(-5.8-1.6*(kh++%3))*kf,cn);
	    else tex->DrawLatex(cumlum+1,(-1.6-0.4*(kh++%3))*kf,cn);
	  }
	}
	else {
	  if (isZmass) tex->DrawLatex(cumlum+1,(+0.15-0.03*(ks++))*kf,cn);
	  else if (isEta) tex->DrawLatex(cumlum+1,(+7.0-0.8*(ks++))*kf,cn);
	  //else if (isEta) tex->DrawLatex(cumlum+1,(+14.0-1.6*(ks++))*kf,cn);
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
      f1->SetRange(382298,1e6); // V2.0 onwards
      hsum->Fit(f1,"QRN");
      hsum->Scale(1./f1->GetParameter(0));
      
      // Turn Xsec into Xsec-1 (%) and map to cumlum2
      TH1D  *h = (TH1D*)hbins->Clone(Form("h_%s",ch)); h->Reset();
      for (int i = 1; i != hsum->GetNbinsX()+1; ++i) {
	int run = hsum->GetBinCenter(i);
	double cumlum = hcumlum2->GetBinContent(hcumlum2->FindBin(run));
	int j = hbins->FindBin(cumlum);
	if (hsum->GetBinContent(i)!=0) {
	  h->SetBinContent(j, (hsum->GetBinContent(i)-1)*100.);
	  h->SetBinError(j, hsum->GetBinError(i)*100.);
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
      f1->SetRange(90.7,1e6); // V2.0 onwards
      f2->SetRange(90.7,1e6); // V2.0 onwards
      
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

