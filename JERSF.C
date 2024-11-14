// Purpose: derive JER SF with MPFX method
//          combine multiple channels, start from TProfile2D *p2m0
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFitResultPtr.h"

#include "tdrstyle_mod22.C"
#include "tools.C"


bool plotEtaBins = true; // lots of plots

bool skipZ = true;//false;
bool skipG = false;

// Patch [2.853,2.964] minimum statistical uncertainty (for trig bias?)
double patchErr29 = 0.005; // added in quadrature to stat unc

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
  f1->SetParameters(10,1,0.05);
  // Fix sensible range for each parameter
  if (!f1ref) { // should be called for MC
    f1->SetParLimits(0,  3.,  12.); // Noise (min was -1)
    f1->SetParLimits(1, 0.5, 1.5); // Stochastic
    f1->SetParLimits(2,0.03,0.07); // Constant
    // Special treatment for EC2
    if (eta>2.500 && eta<3.139) {
      f1->SetParameters(10,1.5,0.10);
      f1->SetParLimits(0, 3., 15); // Noise (min was -1)
      f1->SetParLimits(1, 0, 2.0); // Stochastic
      f1->SetParLimits(2, 0.03, 0.12); // Constant
    }
    // Special treatment for HF
    if (eta>3.139) f1->FixParameter(2,0.065);
    if (eta>4.013) f1->SetParLimits(0, 9., 12.);
  }
  if (f1ref) { // should be called for data
    double N = f1ref->GetParameter(0);
    double S = f1ref->GetParameter(1);
    double C = f1ref->GetParameter(2);
    f1->SetParameters(N, 1.1*S, 1.1*C);
    //f1->SetParLimits(0, N-0.2*fabs(N), N+0.2*fabs(N));
    //f1->SetParLimits(0, N-0.05*max(fabs(N),1.), N+0.05*max(fabs(N),1.));
    f1->SetParLimits(0, N, N+0.05*max(fabs(N),1.));
    //f1->SetParLimits(1, S, 1.5*S);
    //f1->SetParLimits(2, C, 5.0*C);
    //f1->SetParLimits(0, N,  10); // Noise
    f1->SetParLimits(1, S, 1.5); // Stochastic
    f1->SetParLimits(2, C,0.12); // Constant
    // Special treatment for EC2
    if (eta>2.500 && eta<3.139) {
      f1->SetParameters(N, S, min(2*C,0.24));
      f1->SetParLimits(0, N, N+0.10*max(fabs(N),1.));
      f1->SetParLimits(1, S, 2.5); // Stochastic (min 0.75?)
      f1->SetParLimits(2, C, 0.24); // Constant
    }
    // Special treatment for HF
    if (eta>3.489) f1->FixParameter(0, N);
  }
  if (f1ptr) (*f1ptr) = hjer->Fit(f1,"QRNS");
  else                  hjer->Fit(f1,"QRN");
  // If fit was empty, disregard f1ptr
  if (f1->GetNDF()<=0 && f1ptr) (*f1ptr) = 0;
  
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

  //gROOT->ProcessLine(".! mkdir pdf/JERSF");
  gROOT->ProcessLine(".! mkdir pdf/JERSF/vsEta");
  //gROOT->ProcessLine(".! mkdir pdf/JERSF/vsPt");
  
  // Set output directory;
  TFile *fout = new TFile("rootfiles/JERSF.root","RECREATE");
  fout->mkdir("Fits");
  fout->mkdir("Dijet");
  fout->mkdir("Extras");
  fout->mkdir("Raw");

  //string vrun[] = {"2023Cv123","2023Cv4","2023D"};//,"2023Cv4D"};
  //string vrun[] = {"2024B","2024C"};
  //string vrun[] = {"2024BC"};
  /*
  string vrun[] = {"2022C_Prompt2022","2022D_Prompt2022", //Pr22
                   "2022C_22Sep2023","2022D_22Sep2023","2022E_22Sep2023",//22Sep
		   "2022F_22Sep2023","2022G_22Sep2023", //22Sep
		   "2023Cv123_Prompt2023","2023Cv4_Prompt2023", //Pr23
		   "2023D_Prompt2023", //Pr23
		   "2022C_19Dec2023","2022D_19Dec2023","2022E_19Dec2023",//19Dec
		   "2022F_19Dec2023","2022G_19Dec2023", //19Dec
		   "2023Cv123_19Dec2023","2023Cv4_19Dec2023", //19Dec
		   "2023D_19Dec2023", //19Dec
		   "2024BC"}; //Pr24
  */
  //string vrun[] = {"2024C","2024D","2024CD"};
  //string vrun[] = {"2023Cv123","2023D","2024CD"};
  //string vrun[] = {"2024B","2024C","2024D","2024BCD"};
  //string vrun[] = {"2024BCD","2024BR"};
  //string vrun[] = {"2024C","2024BR","2024CR","2023D"};
  //string vrun[] = {"2024CS","2024CR","2024CP","2023D","2023Cv123"};
  //string vrun[] = {"2024BCD","2024C","2024E","2024Ev2","2024CP","2024CS","2024CR","2023D","2023Cv123"};
  //string vrun[] = {"2024FG","2024G","2024F","2024E","2024BCD"};
  string vrun[] = {"2024I","2024H","2024G","2024F","2024E","2024BCD"};
  //string vrun[] = {"2024CT","2024CS"};
  //string vrun[] = {"2024CT","2024CS","2024CP","2024F",
  //		   "2023D","2023Cv123"};
  //string vrun[] = {"2024C"};
  //string vrun[] = {"2024Ev2"};
  //string vrun[] = {"2024BCD"};
  //string vrun[] = {"2024BCD","2024CR"};
  //string vrun[] = {"2024E","2024C","2024D","2024BCD"};
  //string vrun[] = {"2024BCD","2024E","2024CR"};
  //string vrun[] = {"2023D"};
  //string vrun[] = {"2024Ev2"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  //string vmc[] = {"Summer23","Summer23","Summer23BPIX"};//,"Summer23BPIX"};
  //string vmc[] = {"Summer23MG","Summer23MG","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix"};
  /*
  string vmc[] = {"Summer23MGBPix","Summer23MGBPix", //Pr22
		  "Summer23MGBPix","Summer23MGBPix","Summer23MGBPix", //22Sep
		  "Summer23MGBPix","Summer23MGBPix", //22Sep
		  "Summer23MGBPix","Summer23MGBPix", //Pr23
		  "Summer23MGBPix", //Pr23
		  "Summer23MGBPix","Summer23MGBPix","Summer23MGBPix", //19Dec
		  "Summer23MGBPix","Summer23MGBPix", //19Dec
		  "Summer23MGBPix","Summer23MGBPix", //19Dec
		  "Summer23MGBPix", //19Dec
		  "Summer23MGBPix"}; //Pr24
  */
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix",
  //		  "Summer23MGBPix"};
  //"Winter24MCFlat"};
  //string vmc[] = {"Summer23MG_Cv123","Summer23MG_Cv4","Summer23MGBPix_D"};
  //string vmc[] = {"Summer23BPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix",
  //		  "Summer23MGBPix","Summer23MGBPix","Summer23MGBPix",
  //		  "Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  //string vmc[] = {"Winter24","Winter24","Winter24","Winter24","Winter24"};
  string vmc[] = {"Winter24","Winter24","Winter24","Winter24","Winter24","Winter24"};
  //string vmc[] = {"Winter24","Winter24"};
  //string vmc[] = {"Winter24","Winter24","Winter24","Winter24",
  //		  "Summer23MGBPix","Summer23MG"};
  //string vmc[] = {"Summer23MGBPix","Summer23MGBPix","Summer23MGBPix","Summer23MGBPix"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);
  assert(nmc==nrun);
  for (int irun = 0; irun != nrun; ++irun) {
    string run = vrun[irun];
    const char *cr = run.c_str();
    string mc = vmc[irun];
    const char *cm = mc.c_str();
    
    //gROOT->ProcessLine(Form(".! mkdir pdf/JERSF/%s/vsEta",cr));
  // (No indent here for the resf of the loop, maybe function call later)
    
  //string run = "2023D";
  //const char *cr = run.c_str();
  //string mc = "Summer23MGBPix";
  //const char *cm = mc.c_str();

    TFile *f(0), *fm(0), *fz(0), *fg(0), *fgm(0);
  if (TString(cr).Contains("Prompt2022")) {
    string run2 = TString(cr).ReplaceAll("_Prompt2022","").Data();
    const char *cr2 = run2.c_str();
    f = new TFile(Form("rootfiles/Prompt2022CD/v42_2022_Prompt/jmenano_data_cmb_%s_JME_v42_2022_Prompt.root",cr2),"READ");
    // Placeholder MC for all options
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
  }
  else if (TString(cr).Contains("Prompt2023")) {
    string run2 = TString(cr).ReplaceAll("_Prompt2023","").Data();
    const char *cr2 = run2.c_str();
    //f = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_data_cmb_%s_JME_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cr2),"READ");
    f = new TFile(Form("rootfiles/Summer23_L2L3ResJERSF/v39_2023_etabin_SFv2/jmenano_data_cmb_%s_JME_v39_2023_etabin_SFv2.root",cr2),"READ");
    // Placeholder MC for all options
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
  }
  else if (TString(cr).Contains("22Sep2023")) {
    string run2 = TString(cr).ReplaceAll("_22Sep2023","").Data();
    const char *cr2 = run2.c_str();
    f = new TFile(Form("rootfiles/jmenano_data_cmb_%s_JME_v39_2022_etabin_SFv.root",cr2),"READ");
    // Placeholder MC for all options
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
  }
  else if (TString(cr).Contains("19Dec2023")) {
    string run2 = TString(cr).ReplaceAll("_19Dec2023","").Data();
    const char *cr2 = run2.c_str();
    cout << "Reading in " << cr2 << endl << flush;
    f = new TFile(Form("rootfiles/19Dec2023/jmenano_data_cmb_%s_v39_eta_METcut_19Dec2023.root",cr2),"READ");
    // Placeholder MC for all options
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
  }
  else if (run=="2024CT") { // ECAL_CC_HCAL_DI
    f = new TFile(Form("rootfiles/Prompt2024/v90_2024/jmenano_data_cmb_%s_JME_v90_2024.root", cr),"READ"); // re-reco JEC
    fm = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_mc_cmb_Winter24MG_v83_2024.root","READ"); // no JER SF (v89 has JER SF)
    //
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v86.root","2024C"),"READ"); // regular 2024C, not re-reco
    //
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALCC-HCALDI_w35.root","2024C"),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w35.root","READ"); // as in regular 2024
  }
  else if (run=="2024BS" || run=="2024CS") { // ECALRATIO_HCALDI
    //f = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%srs_JME_v76_2024.root", run=="2024CS" ? "2024C" : "2024B"),"READ"); // prompt JEC
    f = new TFile(Form("rootfiles/Prompt2024/v77_2024/jmenano_data_cmb_%sS_JME_v77_2024.root", run=="2024CS" ? "2024C" : "2024B"),"READ"); // re-reco JEC
    //fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    fm = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_mc_cmb_Winter24MG_v83_2024.root","READ"); // no JER SF (v89 has JER SF)
    //
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ"); // May 16 golden, 12.3/fb
    //
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALR-HCALDI_w29.root",run=="2024CS" ? "2024C" : "2024B"),"READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w35.root","READ"); // as in regular 2024
  }
  else if (run=="2024BR" || run=="2024CR") { // ECALRATIO
    f = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_ECAL_JME_v76_2024.root", run=="2024CR" ? "2024C" : "2024B"),"READ");
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    //
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ"); // May 16 golden, 12.3/fb
    //
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w29.root",run=="2024CR" ? "2024C" : "2024B"),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
    assert(false); // not updated to Winter24 yet
  }
  else if (run=="2024BP" || run=="2024CP") { // Prompt (ECAL CC timing)
    f = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root", run=="2024CP" ? "2024C" : "2024B"),"READ");
    //fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    fm = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_mc_cmb_Winter24MG_v83_2024.root","READ"); // no JER SF (v89 has JER SF)
    //
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024BCD"),"READ"); // May 16 golden, 12.3/fb
    //
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w29.root",run=="2024CP" ? "2024C" : "2024B"),"READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w35.root","READ"); // as in regular 2024
  }
  /*
  else if (run=="2024BR" || run=="2024CR") {
    if (run=="2024BR")
      //f = new TFile(Form("rootfiles/Prompt2024/v62_2024/jmenano_data_cmb_%s_ECAL_v62_2024.root","2024B"),"READ");
      f = new TFile(Form("rootfiles/Prompt2024/v69_2024/jmenano_data_cmb_%s_ECAL_v69_2024.root","2024B"),"READ");
    else
      //f = new TFile(Form("rootfiles/Prompt2024/v66_2024/jmenano_data_cmb_%s_ECAL_v66_2024.root","2024C"),"READ"); // One half
      f = new TFile(Form("rootfiles/Prompt2024/v70_2024/jmenano_data_cmb_%s_ECAL_JME_v70_2024.root","2024C"),"READ");
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    //
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v81.root","2024BCD"),"READ"); // May 16 golden, 12.3/fb
    if (run=="2024BR") 
      //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w26.root","2024B"),"READ");
      fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w28.root","2024B"),"READ");
    else
      //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w26.root","2024C"),"READ");
      fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s-ECALRATIO_w28.root","2024C"),"READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ"); // as in regular 2024
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosMix_mc_2023-BPixP8QCD_w16.root","READ"); // QCD added as in 2023D below => missing p2gsm
    //fgm = new TFile(run=="2023D" ? "rootfiles/Prompt2024/GamHistosFill_mc_2023QCD-BPix_w12.root" : "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8_w6.root","READ"); // Summer23 (w3->w2)
  }
  */
  /*
  else if (run=="2024Ev2") {
    //f = new TFile(Form("rootfiles/Prompt2024/v71_2024/jmenano_data_cmb_%s_JME_v71_2024.root","2024E_v2"),"READ");
    f = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root","2024Ev2"),"READ");
    //fm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ"); // JER SF!!
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root","2024E"),"READ"); // v82
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w27.root",cr),"READ");
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w29.root",cr),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
  }
  else if (run=="2024E") {
    f = new TFile(Form("rootfiles/Prompt2024/v68_2024/jmenano_data_cmb_%s_JME_v68_2024.root",cr),"READ");
    //fm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ"); // JER SF!!
    fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root",cr),"READ"); // v82
    fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w27.root",cr),"READ");
    fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
  }
  */
  //else if (run=="2024F" || run=="2024EF") {
  //f = new TFile(Form("rootfiles/Prompt2024/v86_2024/jmenano_data_cmb_%s_JME_v86_2024.root",cr),"READ");
  //fm = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_mc_cmb_Winter24MG_v83_2024.root","READ"); // no JER SF?
  //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v85.root",cr),"READ");
  //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w33.root",cr),"READ");
  //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w33.root","READ"); // Winter24 photon+jet only
  //}
  else if (TString(cr).Contains("2024")) {
    //f = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_data_cmb_%s_v39_2024_Prompt_etabin_DCSOnly.root",cr),"READ");
    //f = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_data_cmb_%s_JME_v39_2024_Prompt_etabin_DCSOnly.root",cr),"READ"); // no ZB
    //f = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF/jmenano_data_cmb_%s_JME_v39_2024_Prompt_eta_SFD_DCSOnly_Filter_HLT_MPF.root",cr),"READ"); // MET/sumET<0.3
    //f = new TFile(Form("rootfiles/Prompt2024/jmenano_data_cmb_%s_JME_v39_2024_Prompt_Golden_29April.root",cr),"READ"); // 0.74/fb
    //f = new TFile(Form("rootfiles/Prompt2024/v41_2024_Golden/jmenano_data_cmb_%s_JME_v41_2024_Golden.root",cr),"READ"); // 3/fb
    //f = new TFile(Form("rootfiles/Prompt2024/v43_2024_Golden/jmenano_data_cmb_%s_JME_v43_2024_Golden.root",cr),"READ"); // 3/fb closure
    //f = new TFile(Form("rootfiles/Prompt2024/v50_1_2024_DCSOnly/jmenano_data_cmb_%s_JME_v50_1_2024_DCSOnly.root",cr),"READ"); // May 16 DCSOnly
    //f = new TFile(Form("rootfiles/Prompt2024/v50_2024/jmenano_data_cmb_%s_JME_v50_2024.root",cr),"READ"); // May 16 golden, 12.3/fb
    //f = new TFile(Form("rootfiles/Prompt2024/v67_2024/jmenano_data_cmb_%s_JME_v67_2024.root",cr),"READ");
    //f = new TFile(Form("rootfiles/Prompt2024/v83_2024/jmenano_data_cmb_%s_JME_v83_2024.root",cr),"READ");
    //f = new TFile(Form("rootfiles/Prompt2024/v83_2024/jmenano_data_cmb_%s_JME_v83_2024.root",cr),"READ");
    //f = new TFile(Form("rootfiles/Prompt2024/v89_2024/jmenano_data_cmb_%s_JME_v89_2024.root",cr),"READ"); // Aug 8 hybrid, V5M JEC
    //f = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_data_cmb_%s_JME_v109_2024.root",cr),"READ"); // V5M->V6M JEC
    if (run=="2024H" || run=="2024I")
      f = new TFile(Form("rootfiles/Prompt2024/v111_2024/jmenano_data_cmb_%s_JME_v111_2024.root",cr),"READ"); // V6M->V7M
    else
      f = new TFile(Form("rootfiles/Prompt2024/v110_2024/jmenano_data_cmb_%s_JME_v110_2024.root",cr),"READ"); // V6M closure
    //fm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ"); // with JER SF?
    //fm = new TFile("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_Summer23MGBPix_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root","READ"); // no JER SF
    //fm = new TFile("rootfiles/Prompt2024/v83_2024/jmenano_mc_cmb_Winter24MG_v83_2024.root","READ"); // no JER SF?
    //fm = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_mc_out_Winter24MGV14_v109_%s.root",cr),"READ"); // no JER SF?
    if (run=="2024H" || run=="2024I")
      fm = new TFile(Form("rootfiles/Prompt2024/v111_2024/jmenano_mc_out_Winter24MGV14_v111_%s.root",cr),"READ"); // Winter24
    else
      fm = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_mc_out_Winter24MGV14_v109_%s.root",cr),"READ"); // Winter24
    fm = (TFile*)fm->GetDirectory("HLT_MC");
    //fm = new TFile("rootfiles/Prompt2024/v89_2024/jmenano_mc_cmb_Winter24MG_v89_2024.root","READ"); // MCV14, JER SF
    //if (run=="2024BCD") {
    //fm = new TFile("rootfiles/Prompt2024//v51_2024_Golden/jmenano_mc_cmb_Winter24MCFlat_v51_2024_Golden.root","READ");
    //}
    //
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v78.root",cr),"READ"); // v77: Golden 2024B, DCSOnly 2024C
    //fz = new TFile("rootfiles/Prompt2024/jme_bplusZ_2024BC_Zmm_sync_v78golden.root","READ"); // 0.74/fb
    //fz = new TFile("rootfiles/Prompt2024/jme_bplusZ_2024BC_Zmm_sync_v79golden.root","READ"); // 3/fb
    //fz = new TFile("rootfiles/Prompt2024/jme_bplusZ_2024BC_Zmm_sync_v80DCSOnly.root","READ"); // 3/fb DCSOnly
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v81.root",cr),"READ"); // May 16 golden, 12.3/fb
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v82.root",cr),"READ");
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root",cr),"READ");
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v85.root",cr),"READ");
    //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v86.root",cr),"READ"); // Aug 8 hybrid, V5M
    //fz = new TFile(Form("rootfiles/Prompt2024/v87/jme_bplusZ_%s_Zmm_sync_v87.root",cr),"READ"); // V5M->V6M
    if (run=="2024H" || run=="2024I")
      fz = new TFile(Form("rootfiles/Prompt2024/v88/jme_bplusZ_%s_Zmm_sync_v88b.root",cr),"READ"); // V6M->V7M
    else
      fz = new TFile(Form("rootfiles/Prompt2024/v88/jme_bplusZ_%s_Zmm_sync_v88.root",cr),"READ"); // V6M closure
    //fz = new TFile("rootfiles/Prompt2024/jme_bplusZ_2024BC_Zmm_sync_v80golden.root","READ"); // 3/fb closure
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w12.root",cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w13.root",cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w14.root",cr),"READ"); // 0.74/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w16.root",cr),"READ"); // 3/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w17.root",cr),"READ");//string(cr)=="2024CD"? "2024C" : cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w21.root",cr),"READ");//string(cr)=="2024CD"? "2024C" : cr),"READ"); // DCSOnly
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w22.root",cr),"READ"); // May 16 golden, 12.3/fb
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w27.root",cr),"READ");
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w33.root",cr),"READ");
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w35.root",cr),"READ"); // Aug 8 hybrid, V5M
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w38.root",cr),"READ"); // V5M->V6M
    if (run=="2024H" || run=="2024I")
      fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%sskim_w40.root",cr),"READ"); // V6M->V7M
    else
      fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w39.root",cr),"READ"); // // V6M closure
    //fg = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s_w18.root",cr),"READ"); // 3/fb closure
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w12.root","READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023QCD-BPix_w12.root","READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root","READ");
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w33.root","READ"); // Winter24 photon+jet only
    //fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w35.root","READ"); // Winter24 photon+jet only
    fgm = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_%s-pu_w38.root",cr),"READ"); // Winter24 photon+jet only
  }
  else if (TString(cr).Contains("2023")) {
    f = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_data_cmb_%s_JME_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cr),"READ");
    fm = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_%s_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",cm),"READ");
    //
    fz = new TFile(Form("rootfiles/Summer23_L2ResOnly/jme_bplusZ_%s_Zmm_sync_v70.root",cr),"READ");
    fg = new TFile(Form("rootfiles/Summer23_L2ResOnly/GamHistosFill_data_%s_w6.root",cr),"READ");
    //fgm = new TFile(run=="2023D" ? "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8-BPix_w6.root" : "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8_w6.root","READ"); // Summer23 (w3->w2)
    //fgm = new TFile(run=="2023D" ? "rootfiles/Prompt2024/GamHistosFill_mc_2023QCD-BPix_w12.root" : "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8_w6.root","READ"); // Summer23 (w3->w2)
    fgm = new TFile(run=="2023D" ? "rootfiles/Prompt2024/GamHistosFill_mc_2023P8-BPix_w16.root" : "rootfiles/Summer23_L2ResOnly/GamHistosFill_mc_2023P8_w6.root","READ"); // Summer23 (w3->w2), without QCD and w16 as for 2024BR
  }
  else if (TString(cr).Contains("2022")) {
    assert(false);
  }
  assert(f && !f->IsZombie());
  assert(fm && !fm->IsZombie());
  assert(fz && !fz->IsZombie() || skipZ);
  assert(fg && !fg->IsZombie() || skipG);
  assert(fgm && !fgm->IsZombie() || skipG);

  curdir->cd();

  TProfile2D *p2s, *p2x, *p2sm, *p2xm;
  p2s = (TProfile2D*)f->Get("Dijet2/p2m0");  assert(p2s);
  p2x = (TProfile2D*)f->Get("Dijet2/p2m0x"); assert(p2x);
  p2sm = (TProfile2D*)fm->Get("Dijet2/p2m0");  assert(p2sm);
  p2xm = (TProfile2D*)fm->Get("Dijet2/p2m0x"); assert(p2xm);

  TProfile2D *p2zs(0), *p2zx(0), *p2zsm(0), *p2zxm(0);
  if (!skipZ) {
    p2zs = (TProfile2D*)fz->Get("data/l2res/p2m0");  assert(p2zs);
    p2zx = (TProfile2D*)fz->Get("data/l2res/p2m0x"); assert(p2zx);
    p2zsm = (TProfile2D*)fz->Get("mc/l2res/p2m0");  assert(p2zsm);
    p2zxm = (TProfile2D*)fz->Get("mc/l2res/p2m0x"); assert(p2zxm);
  }
  
  TProfile2D *p2gs(0), *p2gx(0), *p2gsm(0), *p2gxm(0);
  if (!skipG) {
    p2gs = (TProfile2D*)fg->Get("Gamjet2/p2m0");  assert(p2gs);
    p2gx = (TProfile2D*)fg->Get("Gamjet2/p2m0x"); assert(p2gx);
    p2gsm = (TProfile2D*)fgm->Get("Gamjet2/p2m0");  assert(p2gsm);
    p2gxm = (TProfile2D*)fgm->Get("Gamjet2/p2m0x"); assert(p2gxm);
  }
  
  TH1D *hjer13  = getJER(p2s, p2x, 0, 0,1.3,Form("hjer13_%s",cr));
  TH1D *hjer13m = getJER(p2sm,p2xm,0, 0,1.3,Form("hjer13m_%s",cr));
  TH1D *hsf13 = (TH1D*)hjer13->Clone(Form("hsf13_%s",cr));
  hsf13->Divide(hjer13m);

  TH1D *hjer13g(0), *hjer13gm(0), *hsf13g(0);
  if (!skipG) {
    hjer13g  = getJERZ(p2gs, p2gx,  0,1.3,Form("hjer13g_%s",cr));
    hjer13gm = getJERZ(p2gsm,p2gxm, 0,1.3,Form("hjer13gm_%s",cr));
    hsf13g = (TH1D*)hjer13g->Clone(Form("hsf13g_%s",cr));
    hsf13g->Divide(hjer13gm);
  }

  TH1D *hjer13z(0), *hjer13zm(0), *hsf13z(0);
  if (!skipZ) {
    hjer13z  = getJERZ(p2zs, p2zx,  0,1.3,Form("hjer13z_%s",cr));
    hjer13zm = getJERZ(p2zsm,p2zxm, 0,1.3,Form("hjer13zm_%s",cr));
    hsf13z = (TH1D*)hjer13z->Clone(Form("hsf13z_%s",cr));
    hsf13z->Divide(hjer13zm);
  }

  TF1 *f13m = fitJER(hjer13m,50,1500,0,Form("f13m_%s",cr),kBlue-9);
  TF1 *f13  = fitJER(hjer13, 50,1500,0,Form("f13_%s",cr),kBlue,f13m);
  TF1 *f13r = ratioJER(f13,f13m,Form("f13r_%s",cr),kBlue);

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  
  // Draw reference region JER
  double eps = 1e-4;
  TH1D *h0 = tdrHist("h0","JER",0.+eps,0.3-eps);
  TH1D *h0d = tdrHist("h0d","JER SF",0.5+eps,3.0-eps);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *c0 = tdrDiCanvas(Form("c13_%s",cr),h0,h0d,8,11);

  c0->cd(1);
  gPad->SetLogx();

  tex->DrawLatex(0.40,0.85,"|#eta| < 1.3");

  tdrDraw(hjer13m,"HIST",kNone,kBlack,kSolid,-1,kNone,0);
  f13m->SetLineColor(kBlack);
  f13m->SetLineStyle(kDashed);
  f13m->Draw("SAME");
  tdrDraw(hjer13,"Pz",kFullCircle,kBlack);
  f13->SetLineColor(kBlack);
  f13->Draw("SAME");

  if (!skipG) {
    tdrDraw(hjer13gm,"HISTE",kNone,kBlue,kSolid,-1,kNone,0);
    tdrDraw(hjer13g,"Pz",kOpenCircle,kBlue);

    TLegend *leg0 = tdrLeg(0.65,0.90-(skipZ ? 3 : 4)*0.05,0.90,0.90);

    if (!skipZ) {
      //tdrDraw(hjer13zm,"HISTE",kNone,kRed,kSolid,-1,kNone,0);
      //tdrDraw(hjer13z,"Pz",kOpenDiamond,kRed,kSolid,-1,kNone,0);
      //leg0->AddEntry(hjer13z,"Z(#mu#mu) + jet", "PLE");
      c0->cd(2);
      tdrDraw(hsf13z,"Pz",kOpenDiamond,kRed);
      leg0->AddEntry(hsf13z,"Z(#mu#mu) + jet", "PLE");
    }
    leg0->AddEntry(hjer13g,"#gamma + jet", "PLE");
    leg0->AddEntry(hjer13,"Dijet", "PLE");
    leg0->AddEntry(hjer13m,"MC", "L");

  }

  c0->cd(2);
  gPad->SetLogx();
  l->DrawLine(15,1,3500,1);
  
  tdrDraw(hsf13,"Pz",kFullCircle,kBlack);
  f13r->SetLineColor(kBlack);
  f13r->Draw("SAME");

  if (!skipG) {
    tdrDraw(hsf13g,"Pz",kOpenCircle,kBlue);

    //if (!skipZ) {
    //tdrDraw(hsf13z,"Pz",kOpenDiamond,kRed);
    //}
  }

  c0->SaveAs(Form("pdf/JERSF/JERSF_Eta13_%s.pdf",cr));
  if (run=="2024BCD") c0->SaveAs(Form("pdf/JERSF/JERSF_Eta13_%s.png",cr));
  
  //TCanvas *cx = new TCanvas("cx","cx",7*250,3*250);
  //cx->Divide(7,3,0,0);
  TCanvas *cx = new TCanvas(Form("cx_%s",cr),"cx",9*300,5*300);
  cx->Divide(9,5,0,0);
  TH2D *h2jersf = p2s->ProjectionXY(Form("h2jersf_%s",cr)); h2jersf->Reset();
  TH2D *h2jersf0 = p2s->ProjectionXY(Form("h2jersf0_%s",cr)); h2jersf0->Reset();
  //
  TH2D *h2jerdt = p2s->ProjectionXY(Form("h2jerdt_%s",cr)); h2jerdt->Reset();
  TH2D *h2jermc = p2s->ProjectionXY(Form("h2jermc_%s",cr)); h2jermc->Reset();
  //
  TH2D *h2zjer(0), *h2zjerdt(0), *h2zjermc(0);
  if (!skipZ) {
    h2zjer = p2zs->ProjectionXY(Form("h2zjer_%s",cr)); h2zjer->Reset();
    h2zjerdt =p2zs->ProjectionXY(Form("h2zjerdt_%s",cr)); h2zjerdt->Reset();
    h2zjermc =p2zs->ProjectionXY(Form("h2zjermc_%s",cr)); h2zjermc->Reset();
  }
  //
  TH2D *h2gjer(0), *h2gjerdt(0), *h2gjermc(0);
  if (!skipG) {
    h2gjer = p2gs->ProjectionXY(Form("h2gjer_%s",cr)); h2gjer->Reset();
    h2gjerdt =p2gs->ProjectionXY(Form("h2gjerdt_%s",cr)); h2gjerdt->Reset();
    h2gjermc =p2gs->ProjectionXY(Form("h2gjermc_%s",cr)); h2gjermc->Reset();
  }
  
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

  // Inflate [2.853,2.964] minimum uncertainty to avoid bad chi2 (trig bias?)
  if (patchErr29 && eta>2.853 && eta<2.964) {
    for (int i = 1; i!= hjer->GetNbinsX()+1; ++i) {
      double err = hjer->GetBinError(i);
      if (err!=0) {
	hjer->SetBinError(i, sqrt(pow(err,2)+pow(patchErr29,2)));
      }
    }
  } // patchErr29
   
  TH1D *hjerz(0), *hjerzm(0), *hsfz(0);
  if (!skipZ) {
    hjerz  = getJERZ(p2zs,p2zx,etamin,etamax,Form("hjerz_%d_%s",ieta,cr));
    hjerzm = getJERZ(p2zsm,p2zxm,etamin,etamax,Form("hjerzm_%d_%s",ieta,cr));
    hsfz = (TH1D*)hjerz->Clone(Form("hsfz_%d_%s",ieta,cr));
    hsfz->Divide(hjerzm);
  }

  TH1D *hjerg(0), *hjergm(0), *hsfg(0);
  if (!skipG) {
    hjerg  = getJERZ(p2gs,p2gx,etamin,etamax,Form("hjerg_%d_%s",ieta,cr));
    hjergm = getJERZ(p2gsm,p2gxm,etamin,etamax,Form("hjergm_%d_%s",ieta,cr));
    hsfg = (TH1D*)hjerg->Clone(Form("hsfg_%d_%s",ieta,cr));
    hsfg->Divide(hjergm);
  }
  
  //double ptmin = 50;
  double ptmin = 30;
  //double ptmax = 1500./cosh(eta1);
  double ptmax = min(2000.,5500./cosh(eta1));
  if (run=="2024BCD") ptmax = min(1500.,4890./cosh(eta1));
  // Clean up some high chi2 bins
  //if (eta>2.964 && eta<3.139) ptmax = 302;
  //if (eta>3.139 && eta<3.314) ptmax = 236;
  // Safeguard against bad trigger in EC2
  if (eta>2.500 && eta<2.964) ptmin = 50;
  // Get one more bin for pT dependence
  if (eta>4.889 && run=="2024F") { ptmin = 60; ptmax = 114.; }
  
  // Skip before empty bins give invalid TFitResultPtr
  //if (hjer->Integral()==0) continue;
  
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

  c1->cd(1);
  gPad->SetLogx();

  tex->DrawLatex(0.35,0.85,Form("[%1.3f,%1.3f]",eta1,eta2));
  tex->DrawLatex(0.35,0.80,Form("#chi^{2}/ndf=%1.1f/%d",
				f1->GetChisquare(),f1->GetNDF()));
  
  if (!skipZ) {
    tdrDraw(hjerzm,"HISTE",kNone,kRed,kSolid,-1,kNone);
    tdrDraw(hjerz,"Pz",kOpenCircle,kRed);
  }

  if (!skipG) {
    tdrDraw(hjergm,"HISTE",kNone,kBlue,kSolid,-1,kNone);
    tdrDraw(hjerg,"Pz",kOpenDiamond,kBlue);
  }

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
  //c1->SaveAs(Form("pdf/JERSF/%s/vsEta/JERSF_eta_%04d_%04d_%s.pdf",
  //		    cr,int(1000.*eta1),int(1000.*eta2),cr));

  // Also draw final results into a giant canvas
  cx->cd(ieta);
  TH1D *h6 = tdrHist(Form("h6_%d_%s",ieta,cr),"JER SF",0.5,2.5);
  double ptmaxe = 13600.*0.5/cosh(eta1);
  //if      (eta<1.653) h6->GetYaxis()->SetRangeUser(0.8,1.6);
  if      (eta<1.566) h6->GetYaxis()->SetRangeUser(0.8,1.6);
  //else if (eta<2.964) h6->GetYaxis()->SetRangeUser(0.5,3.0);
  //else if (eta<2.964) h6->GetYaxis()->SetRangeUser(0.5,3.0);
  //else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.5,2.5);
  else if (eta<5.191) h6->GetYaxis()->SetRangeUser(0.5,3.0);
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
  tex->DrawLatex(0.35,0.80,Form("#chi^{2}/ndf=%1.1f/%d",
				f1->GetChisquare(),f1->GetNDF()));
  
  // Store results with uncertainty to output histogram
  for (int ipt = 1; ipt != h2jersf->GetNbinsY()+1; ++ipt) {
    double pt = h2jersf->GetYaxis()->GetBinCenter(ipt);
    double pt1 = h2jersf->GetYaxis()->GetBinLowEdge(ipt);
    double emax1 = pt1*cosh(eta1);
    double jersf = f1r->Eval(pt);
    //double ejersf = 0.; // do properly later
    double d = f1->Eval(pt);
    double ed = (f1ptr ? tools::getFitErr(f1, f1ptr, pt) : 0);
    double m = f1m->Eval(pt);
    double em = (f1mptr ? tools::getFitErr(f1m, f1mptr, pt) : 0);
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
  if (!skipZ) {
    for (int ipt = 1; ipt != h2zjer->GetNbinsY()+1; ++ipt) {
      h2zjer->SetBinContent(ieta, ipt, hsfz->GetBinContent(ipt));
      h2zjer->SetBinError(ieta, ipt, hsfz->GetBinError(ipt));
      h2zjerdt->SetBinContent(ieta, ipt, hjerz->GetBinContent(ipt));
      h2zjerdt->SetBinError(ieta, ipt, hjerz->GetBinError(ipt));
      h2zjermc->SetBinContent(ieta, ipt, hjerzm->GetBinContent(ipt));
      h2zjermc->SetBinError(ieta, ipt, hjerzm->GetBinError(ipt)); 
    }
  }
  if (!skipG) {
    for (int ipt = 1; ipt != h2gjer->GetNbinsY()+1; ++ipt) {
      h2gjer->SetBinContent(ieta, ipt, hsfg->GetBinContent(ipt));
      h2gjer->SetBinError(ieta, ipt, hsfg->GetBinError(ipt));
      h2gjerdt->SetBinContent(ieta, ipt, hjerg->GetBinContent(ipt));
      h2gjerdt->SetBinError(ieta, ipt, hjerg->GetBinError(ipt));
      h2gjermc->SetBinContent(ieta, ipt, hjergm->GetBinContent(ipt));
      h2gjermc->SetBinError(ieta, ipt, hjergm->GetBinError(ipt)); 
    }
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
  //cy->SaveAs(Form("pdf/JERSF/JERSF_Summary_%s.png",cr));
  cy->SetName(Form("cy_%s",cr));

  
  // Draw chi2 to have quick control of fit quality
  TH1D *hz = tdrHist("hz","Fit #chi^{2} / NDF",0,10,"|#eta|",0,5.2);
  TCanvas *cz = tdrCanvas("cz",hz,8,33,kSquare);

  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);
  
  tdrDraw(hchi2,"HPz",kFullCircle,kBlack,kSolid,1,kNone);
  
  cz->SaveAs(Form("pdf/JERSF/JERSF_Chi2_%s.pdf",cr));


  // Produce output text file (FactorizedJetCorrecter style)
  //ofstream txt(Form("pdf/JERSF/Summer23_%s_JRV2_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV1M_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV2M_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV3M_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV4M_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV5M_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV5Mc_MC_SF_AK4PFPuppi.txt",cr));
  //ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV6M_MC_SF_AK4PFPuppi.txt",cr));
  ofstream txt(Form("textFiles/Prompt24/Prompt24_%s_JRV7M_MC_SF_AK4PFPuppi.txt",cr));
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
  if (!skipZ) {
    h2zjer->Write(Form("h2jersfZjet_%s_%s",cr,cm),TObject::kOverwrite);
    h2zjerdt->Write(Form("h2jerdtZjet_%s_%s",cr,cm),TObject::kOverwrite);
    h2zjermc->Write(Form("h2jermcZjet_%s_%s",cr,cm),TObject::kOverwrite);
  }
  //
  if (!skipG) {
    h2gjer->Write(Form("h2jersfGamjet_%s_%s",cr,cm),TObject::kOverwrite);
    h2gjerdt->Write(Form("h2jerdtGamjet_%s_%s",cr,cm),TObject::kOverwrite);
    h2gjermc->Write(Form("h2jermcGamjet_%s_%s",cr,cm),TObject::kOverwrite);
  }

  //if (!fout->FindObject("Raw")) fout->mkdir("Raw");
  fout->cd("Raw");
  p2s->Write(Form("p2m0dt_%s_%s",cr,cm),TObject::kOverwrite);
  p2sm->Write(Form("p2m0mc_%s_%s",cr,cm),TObject::kOverwrite);
  p2x->Write(Form("p2m0xdt_%s_%s",cr,cm),TObject::kOverwrite);
  p2xm->Write(Form("p2m0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  //
  if (!skipZ) {
    p2zs->Write(Form("p2zm0dt_%s_%s",cr,cm),TObject::kOverwrite);
    p2zsm->Write(Form("pzgm0mc_%s_%s",cr,cm),TObject::kOverwrite);
    p2zx->Write(Form("p2zm0xdt_%s_%s",cr,cm),TObject::kOverwrite);
    p2zxm->Write(Form("p2zm0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  }
  //
  if (!skipG) {
    p2gs->Write(Form("p2gm0dt_%s_%s",cr,cm),TObject::kOverwrite);
    p2gsm->Write(Form("p2gm0mc_%s_%s",cr,cm),TObject::kOverwrite);
    p2gx->Write(Form("p2gm0xdt_%s_%s",cr,cm),TObject::kOverwrite);
    p2gxm->Write(Form("p2gm0xmc_%s_%s",cr,cm),TObject::kOverwrite);
  }
  curdir->cd();

  } // for irunx

  fout->Close();
} // JERSF
