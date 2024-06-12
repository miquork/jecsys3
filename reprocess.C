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
double scaleEM = 1;//0.98; // scale EM+jet in data in absence of QCD bkg
bool scaleEMperEra = true;

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

// Clean L3Res: Z [30,400], photon [230,1200], multijet [800,1800], PF [15,2450]

// Override settings below to use full range (latest true applies)
const bool useFullRange(false);//true);
const bool usePrompt24Range(false);//true);
const bool usePrompt24RangeV2(true);//false);//true);

// jet+Z
double fjzptmin(20);//15);//20);//15);
double fjzptmax(70);//700);//70);//700.);
// Z+jet average
double fzjaptmin(20);
double fzjaptmax(400);//70);

// Z+jet
double fzptmin(70);//15);//70);//15);//30);//15.);//30);//(15.);//40.);   // Z+jet pTmin
double fzptmax(400);//700);//400);//700.);//1300.);  // Z+jet pTmax
//double fzmpfptmin(15); // Z+jet MPF pTmin
//double fzmpfptmax(700.);//1300);// Z+jet MPF pTmax
//double fzbalptmin(15); // Z+jet DB pTmin
//double fzbalptmax(700.);//1300);// Z+jet DB pTmax
//double fzbptmax(300.); // Z+b pTmax

// Photon+jet, pT bins: 700, 850, 1000, 1200, 1450, 1750, 2100, 2500, 3000
double fgptmin(175);//35);//175);//110);//60);//35);//60);//175);//110);//35);//60);//35);//110);//230);//30.);//110.);//50.);
double fgptmax(1000.);//2100.);//1200.);//1450.);//1750.); // extend 1750 when can (v22 up to 3000)


// Multijet and incjet pT bins (to-do: widen multijet high pT bins?)
// pT bins: 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
//    2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
//    4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000

// Multijet, pT bins (800-2500 GeV for V1)
bool rebinMultijet = true;
double fmjptmin(49);//30);//49);//97);//196.);//800);//100.);//700);//1000.);//500);//196);   // Multijet pTmin
double fmjptmax(2500);//2787);//2500.);//1800.);//2400.);//2785.);//3500.);//2000.);//2785.); // Multijet pTmax
bool doMultijetRecoil = false;//true; // recoil pT binning

// PF composition (incjet)
double fijptmin(15);//50);//97);
double fijptmax(3103);//2450);


// Helper functions to handle JEC
FactorizedJetCorrector *_thejec(0);
TF1 *fCorrPt(0);
FactorizedJetCorrector *getFJC(string l1, string l2="", string res="",
			       string path="");
void setEtaPtRho(FactorizedJetCorrector *jec, double eta,double pt,double rho);
Double_t funcCorrPt(Double_t *x, Double_t *p);
double getJEC(FactorizedJetCorrector *jec, double eta, double ptcorr,
	      double rho = 0);
void setEpoch(string epoch);

// put all the different methods in a single file for easy access for everybody
void reprocess(string epoch="") {

  setEpoch(epoch);
  TString tepoch = TString(epoch);
  
  // Set TDR style to have correct graphical settings when storing graphs
  setTDRStyle();

  // override range settings to check full range
  if (useFullRange) {
    fjzptmin = 15; fjzptmax = 700;
    fzptmin  = 15; fzptmax  = 1000;//700;
    fgptmin =  35; fgptmax = 1750;//2100;
    fmjptmin = 30; fmjptmax = 2787;//2000;//2787;
    fijptmin = 15; fijptmax = 2800;//3103;

    doMultijetRecoil = (epoch=="Run24CP" || epoch=="Run24BCD" || epoch=="Run24BCDE" || epoch=="Run24C" || epoch=="Run24E");
  }
  if (usePrompt24Range) {
    //fjzptmin = 15; //fjzptmax = 300;
    fjzptmin = 15; fjzptmax = 85;
    //fzptmin  = 40; fzptmax  = 700;
    fzptmin  = 85; fzptmax  = 400;
    //fgptmin =  35; fgptmax = 1000;
    fgptmin =  230; fgptmax = 1000;
    //fmjptmin = 30; fmjptmax = 2787;
    fmjptmin = 114; fmjptmax = 2000.;
    fijptmin = 15; fijptmax = 3103;
  }
  if (usePrompt24RangeV2) {
    fjzptmin = 15; fjzptmax = 70;//15;
    fzjaptmin = 25; fzjaptmax = 70;//15;
    fzptmin  = 15; fzptmax  = 1000;//700;
    fgptmin =  60;//40;//130;//40;
    fgptmax = 1750;//300;//1750;
    //doMultijetRecoil = true;
    //fmjptmin = 133;//600;//133;
    fmjptmax = 2500;//2787;//2000.;
    fijptmin = 15; fijptmax = 2787;//2500;//3103;
    //if (epoch=="Run24CP" || epoch=="Run24CR" || epoch=="Run24CS") {
    //fmjptmin = 133;
    //}
    //if (epoch=="Run24CR" || epoch=="Run24CS") { // back to jet+Z for low pT
      //fjzptmax = 70;
      //fzptmin = 70;
      //fgptmin = 60;
      //fmjptmin = 56;
      //fmjptmax = 2787;
    //}
    doMultijetRecoil = (epoch=="Run24CP" || epoch=="Run24BCD" || epoch=="Run24BCDE" || epoch=="Run24C" || epoch=="Run24E");
    fmjptmin = (doMultijetRecoil ? 600 : 133);
  }
  
  TDirectory *curdir = gDirectory;

  const char *cep = epoch.c_str();
  TFile *fout = new TFile(Form("rootfiles/jecdata%s.root", cep), "RECREATE");

  ////////////////////////////
  // Z+jet                  //
  ////////////////////////////

  map<string,const char*> mz;
  mz["Run24B"] = "2024B";
  mz["Run24C"] = "2024C";
  mz["Run24D"] = "2024D";
  mz["Run24E"] = "2024E";
  mz["Run24BC"] = "2024BC";
  mz["Run24BCD"] = "2024BCD";
  mz["Run24BCDE"] = "2024BCDE";
  mz["Run24CR"] = "2024BCD";
  mz["Run24CS"] = "2024BCD";
  mz["Run24CP"] = "2024BCD";
  
  TFile *fz(0), *fzjes(0);
  TH1D *hzjes(0);
  //const char *cdz = "../JERCProtoLab/Winter22Run3/L3Residual_Z";
  const char *cdzul = "rootfiles/Sami_20230630";
  if (tepoch.Contains("UL")) {
    fz = new TFile(Form("%s/jme_bplusZ_%s_Zmm_sync_v53.root",
			cdzul,epoch.c_str()),"READ");
  }
  
  //const char *cdz58p1 = "rootfiles/Sami_20230912"; // v58p1
  // v66 files used for 22Sep2023_V3 JECs
  if (epoch=="Run22CD") {
    fz = new TFile("rootfiles/jme_bplusZ_2022CD_Zmm_sync_v66.root","READ");
  }
  if (epoch=="Run22E") {
    fz = new TFile("rootfiles/jme_bplusZ_2022E_Zmm_sync_v66.root","READ");
  }
  if (epoch=="Run22FG") {
    //fz = new TFile("rootfiles/jme_bplusZ_2022FG_Zmm_sync_v66.root","READ");
    // 19Dec2023 update
    fz = new TFile("rootfiles/jme_bplusZ_2022F_Zmm_sync_v66_19Dec.root","READ"); // TMP until FG available
    //exit(0);
  }
  if (epoch=="Run23C123") {
    //fz = new TFile("rootfiles/jme_bplusZ_2023C123_Zmm_sync_v66.root","READ");
    //fz = new TFile("rootfiles/Summer23_noL2L3Res/jme_bplusZ_2023Cv123_Zmm_sync_v69.root","READ"); // Summer23
    //fz = new TFile("rootfiles/Summer23_L2ResOnly/jme_bplusZ_2023Cv123_Zmm_sync_v70.root","READ"); // Summer23 L2Res_V1
    fz = new TFile("rootfiles/Summer23_L2L3Res/jme_bplusZ_2023Cv123_Zmm_sync_v73smearoff.root","READ"); // Summer23 L2L3Res_V2
  }
  if (epoch=="Run23C4") {
    //fz = new TFile("rootfiles/jme_bplusZ_2023C4_Zmm_sync_v66.root","READ");
    //fz = new TFile("rootfiles/Summer23_noL2L3Res/jme_bplusZ_2023Cv4_Zmm_sync_v69.root","READ");
    //fz = new TFile("rootfiles/Summer23_L2ResOnly/jme_bplusZ_2023Cv4_Zmm_sync_v70.root","READ"); // Summer23 L2Res_V1
    fz = new TFile("rootfiles/Summer23_L2L3Res/jme_bplusZ_2023Cv4_Zmm_sync_v73smearoff.root","READ"); // Summer23 L2L3Res_V2
  }
  if (epoch=="Run23D") {
    //fz = new TFile("rootfiles/jme_bplusZ_2023D_Zmm_sync_v66.root","READ");
    //fz = new TFile("rootfiles/Summer23_noL2L3Res/jme_bplusZ_2023D_Zmm_sync_v69.root","READ");
    //fz = new TFile("rootfiles/Summer23_L2ResOnly/jme_bplusZ_2023D_Zmm_sync_v70.root","READ"); // Summer23 L2Res_V1
    fz = new TFile("rootfiles/Summer23_L2L3Res/jme_bplusZ_2023D_Zmm_sync_v73smearoff.root","READ"); // Summer23 L2L3Res_V2
  }
  if (epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" ||epoch=="Run24E" ||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP") {
  //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v78.root",
  //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v78golden.root",
  //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v79golden.root", // 3/fb golden
  //fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v80golden.root", // 3/fb golden (V2M) closure
    fz = new TFile(Form("rootfiles/Prompt2024/jme_bplusZ_%s_Zmm_sync_v83.root", // May16 golden JSON 12.3/fb golden (V3M + V2M closure), v81->v83
			mz[epoch]),"READ");
  }
  if (epoch=="Run3") {
    //fz = new TFile(Form("%s/jme_bplusZ_Run3_Zmm_sync_v59.root",cdz58p1),"READ"); // Sami's combo
    fz = new TFile("rootfiles/jecdataRun3Data.root","READ"); // manual combo
  }

  assert(fz && !fz->IsZombie());
  
  TFile *fmz = fz; 
  assert(fmz && !fmz->IsZombie());
  
  string sr = "eta_00_13";
  const char *cr = sr.c_str();
  const char *cl = CorLevel.c_str();

  TH1D *hcounts(0);
  TH1D *hmz_dt(0), *hmz_mc(0);
  if (tepoch.Contains("UL") ||
      epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G"||
      epoch=="Run22FG" || epoch=="Run22EFG" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C123" ||
      //epoch=="Run23C" ||
      //epoch=="Run23C1" ||epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
      epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" || epoch=="Run24E"||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
      (epoch=="Run3" && false)
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
  if (epoch=="Run3" && true) {
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


  ////////////////////////////
  // Photon+jet             //
  ////////////////////////////

  //TString sfp = epoch.c_str();
  //if (epoch!="Run3") sfp.ReplaceAll("Run","");
  map<string,const char*> mp;
  mp["Run22C"] = "2022C";
  mp["Run22D"] = "2022D";
  mp["Run22CD"] = "2022CD";
  mp["Run22CDE"] = "2022CDE";
  mp["Run22E"] = "2022E";
  mp["Run22F"] = "2022F";
  mp["Run22G"] = "2022G";
  mp["Run22FG"] = "2022FG";
  mp["Run23BC123"] = "2023Cv123"; // clone below
  mp["Run23C123"] = "2023Cv123";
  mp["Run23C4"] = "2023Cv4";
  mp["Run23D"] = "2023D";
  mp["Run23C4D"] = "2023Cv4D";
  mp["Run24B"] = "2024B";
  mp["Run24C"] = "2024C";
  mp["Run24D"] = "2024D";
  mp["Run24E"] = "2024Ev1"; //TMP
  mp["Run24BC"] = "2024BC";
  mp["Run24BCD"] = "2024BCD";
  mp["Run24BCDE"] = "2024BCDE";
  mp["Run24CR"] = "2024C";
  mp["Run24CS"] = "2024C";
  mp["Run24CP"] = "2024C";
  mp["Run3"] = "Run3";
  //TFile *fp = new TFile(Form("../gamjet/files/GamHistosRatio_%s_P8_v21.root",
  //TFile *fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8_v24.root",mp[epoch]),"READ");
  //TFile *fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8_v27.root",mp[epoch]),"READ");
  //TFile *fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_v27.root",mp[epoch]),"READ");
  //TFile *fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_v29.root",mp[epoch]),"READ"); // no L2L3Res

  // 22Sep2023 V3 used v32 files
  //TFile *fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_v32.root",mp[epoch]),"READ"); // L2L3Res_V3
  TFile *fp(0);
  if (epoch=="Run22FG") {
    //fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_19Dec_v33.root",mp[epoch]),"READ"); // 19Dec2023
    assert(false);
  }
  else if (epoch=="Run23C123" || epoch=="Run23C4" || epoch=="Run23D") {
    //fp = new TFile(Form("rootfiles/Summer23_noL2L3Res/GamHistosRatio_%s_P8_w2.root",mp[epoch]),"READ"); // Summer23 no QCD
    //fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_w4.root",mp[epoch]),"READ"); // Summer23 no QCD with L2Res
    //fp = new TFile(Form("rootfiles/Summer23_L2ResOnly/GamHistosRatio_%s_P8%s-noQCD_w4.root",mp[epoch],epoch=="Run23D" ? "BPix" : ""),"READ"); // Summer23 L2Res_V1 (noQCD)
    //fp = new TFile(Form("rootfiles/Summer23_L2ResOnly/GamHistosRatio_%s_P8%sQCD_w6.root",mp[epoch],epoch=="Run23D" ? "BPix" : ""),"READ"); // Summer23 L2Res_V1 (withQCD)
    fp = new TFile(Form("rootfiles/Summer23_L2L3Res/GamHistosRatio_%s_P8%sQCD_w8.root",mp[epoch],epoch=="Run23D" ? "BPix" : ""),"READ"); // Summer23 L2L3Res_V2 (withQCD)
  }
  else if (epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" ||
	   epoch=="Run24E" ||
	   epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE") {
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w12.root",mp[epoch]),"READ");
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w13.root",mp[epoch]),"READ"); // DCSOnly JSON
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w14.root",mp[epoch]),"READ"); // golden JSON
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPix-noQCD_w14.root",mp[epoch]),"READ"); // golden JSON 0.74/fb
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPix-noQCD_w16.root",mp[epoch]),"READ"); // golden JSON 3/fb
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPix-noQCD_w18.root",mp[epoch]),"READ"); // golden JSON 3/fb (V2M) closure
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPix-noQCD_w22.root",mp[epoch]),"READ"); // May 16 golden JSON 12.3/fb (V3M + V2M closure)
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w22.root",mp[epoch]),"READ"); // May 16 golden JSON 12.3/fb (V3M + V2M closure)
    fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w29.root",mp[epoch]),"READ"); // June 16 hybrid JSON (V3M closure)
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPix-noQCD_w12.root",mp[epoch]),"READ");
  }
  else if (epoch=="Run24CR") {
    //fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s-ECALRATIO_P8BPixQCD_w26.root","2024C"),"READ");
    fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s-ECALRATIO_P8BPixQCD_w29.root","2024C"),"READ");
  }
  else if (epoch=="Run24CS") {
    fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s-ECALR-HCALDI_P8BPixQCD_w29.root","2024C"),"READ");
  }
  else if (epoch=="Run24CP") {
    fp = new TFile(Form("rootfiles/Prompt2024/GamHistosRatio_%s_P8BPixQCD_w29.root","2024C"),"READ");
  }
  else {
    fp = new TFile(Form("../gamjet/rootfiles/GamHistosRatio_%s_P8QCD_v32.root",mp[epoch]),"READ"); // L2L3Res_V3
    assert(false);
  }
  assert(fp && !fp->IsZombie());

  ////////////////////////////
  // Multijet               //
  ////////////////////////////

  /*
  const char *cdmj = "rootfiles/Iita_20230824_jetveto";
  TFile *fmj(0);
  map<string,const char*> mmj, mmjm;
  mmj["Run22C"] = "jmenano_data_cmb_2022C_JME_v1.root";
  mmj["Run22D"] = "jmenano_data_cmb_2022D_JME_v1.root";
  mmj["Run22CD"] = "jmenano_data_cmb_2022CD_JME_v1.root";
  mmj["Run22E"] = "jmenano_data_cmb_2022E_JME_v1.root";
  mmj["Run22F"] = "jmenano_data_cmb_2022F_JME_v1.root";
  mmj["Run22G"] = "jmenano_data_cmb_2022G_JME_v1.root";
  mmj["Run22FG"] = "jmenano_data_cmb_2022FG_JME_v1.root";
  mmj["Run22EFG"] = "jmenano_data_cmb_2022EFG_JME_v1.root";
  mmj["Run23B"] = "nano_data_cmb_2023B_JME_v1.root";
  mmj["Run23BC123"] = "nano_data_cmb_2023BCv123_JME_v1.root";
  mmj["Run23C"] = "nano_data_cmb_2023Cv123_JME_v1.root";
  mmj["Run23C4"] = "nano_data_cmb_2023Cv4_JME_v1.root";
  mmj["Run23D"] = "nano_data_cmb_2023D_JME_v1.root";
  mmj["Run23C4D"] = "nano_data_cmb_2023Cv4D_JME_v1.root";
  mmj["Run3"] = "jmenano_data_cmb_2223_JME_v1.root";
  //
  mmjm["Run22C"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run22D"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run22CD"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run22E"] = "jmenano_mc_out_Summer22EE_v1.root";
  mmjm["Run22F"] = "jmenano_mc_out_Summer22EE_v1.root";
  mmjm["Run22G"] = "jmenano_mc_out_Summer22EE_v1.root";
  mmjm["Run22FG"] = "jmenano_mc_out_Summer22EE_v1.root";
  mmjm["Run22EFG"] = "jmenano_mc_out_Summer22EE_v1.root";
  mmjm["Run23B"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run23BC123"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run23C"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run23C4"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run23D"] = "jmenano_mc_out_Summer22_v1.root";
  mmjm["Run23C4D"] = "jmenano_mc_out_Summer22_v1.root";
  //
  mmjm["Run3"] = "jmenano_mc_out_Summer22Both_v1.root";
  TFile *fmjd = new TFile(Form("%s/%s",cdmj,mmj[cep]),"READ");
  assert(fmjd);
  TFile *fmjm = new TFile(Form("%s/%s",cdmj,mmjm[cep]),"READ");
  assert(fmjm);
  */

  map<string,const char*> mmjd;
  mmjd["Run22C"] = "2022C";
  mmjd["Run22D"] = "2022D";
  mmjd["Run22CD"] = "2022CD";
  mmjd["Run22CDE"] = "2022CDE";
  mmjd["Run22E"] = "2022E";
  mmjd["Run22F"] = "2022F";
  mmjd["Run22G"] = "2022G";
  mmjd["Run22FG"] = "2022FG";
  //mmjd["Run22FG"] = "2022G"; // 19Dec patch
  mmjd["Run23BC123"] = "2023BCv123";
  mmjd["Run23C123"] = "2023Cv123";//"2023BCv123";//"2023Cv123";
  mmjd["Run23C4"] = "2023Cv4";
  //mmjd["Run23C4"] = "2023BCv123"; //19Dec patch
  mmjd["Run23D"] = "2023D";
  mmjd["Run23C4D"] = "2023Cv4D";
  mmjd["Run24B"] = "2024B";
  mmjd["Run24C"] = "2024C";
  mmjd["Run24D"] = "2024D";
  mmjd["Run24E"] = "2024E";
  mmjd["Run24BC"] = "2024BC";
  mmjd["Run24BCD"] = "2024BCD";
  mmjd["Run24BCDE"] = "2024BCDE";
  mmjd["Run24CR"] = "2024C";
  mmjd["Run24CS"] = "2024C";
  mmjd["Run24CP"] = "2024C";
  mmjd["Run3"] = "Run3";
  map<string,const char*> mmjm;
  mmjm["Run22C"] = "Summer22MG";
  mmjm["Run22D"] = "Summer22MG";
  mmjm["Run22CD"] = "Summer22MG";
  mmjm["Run22CDE"] = "Summer22MG";
  mmjm["Run22E"] = "Summer22EEMG";
  mmjm["Run22F"] = "Summer22EEMG";
  mmjm["Run22G"] = "Summer22EEMG";
  mmjm["Run22FG"] = "Summer22EEMG";
  mmjm["Run23BC123"] = "Summer22MG";
  mmjm["Run23C123"] = "Summer23MG_Cv123";//GBPix";//"Summer22MG";
  mmjm["Run23C4"] = "Summer23MG_Cv4";//BPix";//"Summer22MG";
  mmjm["Run23D"] = "Summer23MGBPix_D";//"Summer22MG";
  mmjm["Run23C4D"] = "Summer23MGBPix";//22MG";
  mmjm["Run24B"] = "Summer23MGBPix";
  mmjm["Run24C"] = "Summer23MGBPix";
  mmjm["Run24D"] = "Summer23MGBPix";
  mmjm["Run24E"] = "Summer23MGBPix";
  mmjm["Run24BC"] = "Summer23MGBPix";
  mmjm["Run24BCD"] = "Summer23MGBPix";
  mmjm["Run24BCDE"] = "Summer23MGBPix";
  mmjm["Run24CR"] = "Summer23MGBPix";
  mmjm["Run24CS"] = "Summer23MGBPix";
  mmjm["Run24CP"] = "Summer23MGBPix";
  mmjm["Run3"] = "Summer22MG";
  //TFile *fmjd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v32.root",mmjd[epoch]),"READ"); // no L2L3Res (L2Res for 2022 only)
  //assert(fmjd && !fmjd->IsZombie());
  //TFile *fmjm = new TFile(Form("../dijet/rootfiles/jmenano_mc_out_%s_v32.root",mmjm[epoch]),"READ"); // no L2L3Res (L2Res for 2022 only)

  TFile *fmjd(0), *fmjm(0);
  if (epoch=="Run24CR") {
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v66_2024/jmenano_data_cmb_%s_ECAL_v66_2024.root","2024C"),"READ"); // One half
    fmjd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_ECAL_JME_v76_2024.root","2024C"),"READ");
    fmjm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (epoch=="Run24CS") {
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%srs_JME_v76_2024.root","2024C"),"READ"); // prompt JEC
    fmjd = new TFile(Form("rootfiles/Prompt2024/v77_2024/jmenano_data_cmb_%sS_JME_v77_2024.root","2024C"),"READ"); // re-reco JEC
    fmjm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (epoch=="Run24CP") {
    fmjd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root","2024C"),"READ");
    fmjm = new TFile("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_Summer23MGBPix_v39_2023_etabin_SFv2.root","READ");
  }
  else if (TString(epoch.c_str()).Contains("Run24")) {
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_data_cmb_%s_v39_2024_Prompt_etabin_DCSOnly.root",mmjd[epoch]),"READ");
    //fmjd = new TFile(Form("rootfiles/Prompt2024/jmenano_data_cmb_%s_JME_v39_2024_Prompt_Golden_29April.root",mmjd[epoch]),"READ"); // golde 0.74/fb
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v41_2024_Golden/jmenano_data_cmb_%s_JME_v41_2024_Golden.root",mmjd[epoch]),"READ"); // golden 3/fb
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v43_2024_Golden/jmenano_data_cmb_%s_JME_v43_2024_Golden.root",mmjd[epoch]),"READ"); // golden 3/fb closure
    //fmjd = new TFile(Form("rootfiles/Prompt2024/v50_2024/jmenano_data_cmb_%s_JME_v50_2024.root",mmjd[epoch]),"READ"); // May 16 golden, 12.3/fb (V3M + V2M closure)
    fmjd = new TFile(Form("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_%s_JME_v76_2024.root",mmjd[epoch]),"READ"); // June 6 hybrid (V3M closure)
    fmjm = new TFile(Form("rootfiles/Prompt2024/v39_2024_Prompt_etabin_DCSOnly/jmenano_mc_cmb_%s_v39_2023_etabin_SFv2.root",mmjm[epoch]),"READ");
  }
  else {
  // 22Sep2023_V3 used v35a files
  //TFile *fmjd = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_data_cmb_%s_JME_v36_Summer23.root",mmjd[epoch]),"READ"); // Summer23_V1
  //TFile *fmjm = new TFile(Form("rootfiles/Summer23_noL2L3Res/jmenano_mc_cmb_%s_v36_Summer23.root",mmjm[epoch]),"READ"); // Summer23_V1
  //TFile *fmjd = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_data_cmb_%s_JME_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",mmjd[epoch]),"READ"); // Summer23_V1
  //TFile *fmjm = new TFile(Form("rootfiles/Summer23_L2ResOnly/jmenano_mc_cmb_%s_v39_noRwPU_noSmearJets_25Feb2024_L2Res_v1.root",mmjm[epoch]),"READ"); // Summer23 L2Res_V1
  fmjd = new TFile(Form("rootfiles/Summer23_L2L3Res/jmenano_data_cmb_%s_JME_v39_L2Rel_L2L3Res_v2_SF.root",mmjd[epoch]),"READ"); // Summer23 L2L3Res_V2
  fmjm = new TFile(Form("rootfiles/Summer23_L2L3Res/jmenano_mc_cmb_%s_RwPU_v39_SmearJets_L2Res_v1_SF.root",mmjm[epoch]),"READ"); // Summer23 L2L3Res_V2
  //
  //TFile *fmjd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v35a.root",mmjd[epoch]),"READ"); // L2L3Res_V3
  //TFile *fmjm = new TFile(Form("../dijet/rootfiles/jmenano_mc_out_%s_v35a.root",mmjm[epoch]),"READ"); // L2L3Res_V3
  // 19Dec2023 special JetMET re-reco (v35 not working for multijet balance)
  //TFile *fmjd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v35_19Dec2023.root",mmjd[epoch]),"READ"); // 19Dec2023
  //TFile *fmjm = new TFile(Form("../dijet/rootfiles/jmenano_mc_cmb_%s_v35_19Dec2023.root",mmjm[epoch]),"READ"); // 19Dec2023 bugged
  //TFile *fmjm = new TFile(Form("../dijet/rootfiles/jmenano_mc_out_%s_v35a.root",mmjm[epoch]),"READ"); // 19Dec2023 using 22Sep2023
  }
  assert(fmjd && !fmjd->IsZombie());
  assert(fmjm && !fmjm->IsZombie());

  // Inclusive jets from same file as multijets
  TFile *fijd = fmjd;
  TFile *fijm = fmjm;
  // 19Dec2023 special JetMET re-reco first iteration works for incjet
  //TFile *fijd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v35_19Dec2023.root",mmjd[epoch]),"READ"); // 19Dec2023
  //if (epoch=="Run23C4")
  //fijd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v35_19Dec2023.root",mmjd["Run23C123"]),"READ"); // 19Dec2023
  //TFile *fijm = new TFile(Form("../dijet/rootfiles/jmenano_mc_out_%s_v35a.root",mmjm[epoch]),"READ"); // 19Dec2023 using 22Sep2023
  
  assert(fmjd && !fmjd->IsZombie());
  assert(fmjm && !fmjm->IsZombie());

  // Run 2 2D pT-eta distribution down to 5 GeV based on P8 dijet sample
  TFile *feta = new TFile("rootfiles/P8_dijet_45M_TH2D_correct_weighting.root",
			  "READ");
  assert(feta && !feta->IsZombie());
  TH2D *h2pteta = (TH2D*)feta->Get("pTjet_etajet"); assert(h2pteta);

  // Store pointers to all files in a map for easy handling later
  map<string, TFile*> files;
  files["jetz"] = fz;
  files["zjav"] = fz;
  files["zjet"] = fz;
  files["gamjet"] = fp;
  files["multijet_data"] = fmjd;
  files["multijet_mc"] = fmjm;
  files["multijet_ratio"] = fmjd;
  files["incjet_data"] = fijd;
  files["incjet_mc"] = fijm;
  files["incjet_ratio"] = fijm;

  // Map variable names used in different files to common scheme
  map<string, map<string, const char*> > rename;

  // Results from Sami's Z+b analysis
  if (tepoch.Contains("UL") ||
      epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run22FG" || epoch=="Run22EFG" ||
      epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C123" ||
      //epoch=="Run23C" ||
      //epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
      epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" || epoch=="Run24E"||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
      epoch=="Run3"
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

  rename["gamjet"]["ratio"] = "";
  rename["gamjet"]["data"] = "_DATA"; 
  rename["gamjet"]["mc"] = "_MC"; 
  rename["gamjet"]["mpfchs"] = "resp_MPFchs";
  rename["gamjet"]["mpfchs1"] = "resp_MPFchs"; 
  rename["gamjet"]["ptchs"] = "resp_DBchs";
  rename["gamjet"]["counts"] = "RawNEvents_data_vs_pt";
  //
  rename["gamjet"]["mpf1"] = "resp_MPFR1chs";
  rename["gamjet"]["mpfn"] = "resp_MPFRnchs";
  rename["gamjet"]["mpfu"] = "resp_MpfRuchs";
  rename["gamjet"]["rho"] = "resp_Rho_CHS";

  rename["gjet"]["ratio"] = "";
  rename["gjet"]["data"] = "_DATA";
  rename["gjet"]["mc"] = "_MC";
  rename["gjet"]["mpfchs"] = "respMPFchs";
  rename["gjet"]["mpfchs1"] = "respMPFchs";
  rename["gjet"]["ptchs"] = "resp_DBchs";
  rename["gjet"]["counts"] = "resp_";
  rename["gjet"]["chf"] = "chHEF";
  rename["gjet"]["nef"] = "neEmEF";
  rename["gjet"]["nhf"] = "neHEF";
  rename["gjet"]["cef"] = "chEmEF";
  rename["gjet"]["muf"] = "muEF";

  // Multijet average
  rename["multijet"]["mpfchs"] = "pm0m";
  rename["multijet"]["mpfchs1"] = "pm0m";
  rename["multijet"]["ptchs"] = "pm2m";
  rename["multijet"]["counts"] = "hptm_sel";
  rename["multijet"]["mpf1"] = "pm2m";
  rename["multijet"]["mpfn"] = "pmnm";
  rename["multijet"]["mpfu"] = "pmum";
  rename["multijet"]["crecoil"] = "pcrecoilm";
  if (doMultijetRecoil) {
    rename["multijet"]["mpfchs"] = "pm0r";
    rename["multijet"]["mpfchs1"] = "pm0r";
    rename["multijet"]["ptchs"] = "pm2r";
    rename["multijet"]["counts"] = "hptr_sel";
    rename["multijet"]["mpf1"] = "pm2r";
    rename["multijet"]["mpfn"] = "pmnr";
    rename["multijet"]["mpfu"] = "pmur";
    rename["multijet"]["crecoil"] = "pcrecoilr";
  }
  //rename["multijet"]["chf"] = "PFcomposition/pchf13";
  //rename["multijet"]["nef"] = "PFcomposition/pnef13";
  //rename["multijet"]["nhf"] = "PFcomposition/pnhf13";
  //rename["multijet"]["cef"] = "PFcomposition/pcef13";
  //rename["multijet"]["muf"] = "PFcomposition/pmuf13";
  //rename["multijet"]["rho"] = "PFcomposition/prho13";
  rename["multijet"]["chf"] = "../Dijet/PFcomposition/pchf13";
  rename["multijet"]["nef"] = "../Dijet/PFcomposition/pnef13";
  rename["multijet"]["nhf"] = "../Dijet/PFcomposition/pnhf13";
  rename["multijet"]["cef"] = "../Dijet/PFcomposition/pcef13";
  rename["multijet"]["muf"] = "../Dijet/PFcomposition/pmuf13";
  rename["multijet"]["rho"] = "../Dijet/PFcomposition/prho13";
  
  
  // color and style codes
  map<string, map<string, int> > style;

  style["jetz"]["mpfchs1"] = kFullDiamond;
  style["jetz"]["ptchs"] = kOpenDiamond;
  style["jetz"]["mpf1"] = kFullTriangleUp;
  style["jetz"]["mpfn"] = kFullTriangleUp;
  style["jetz"]["mpfu"] = kFullTriangleDown;
  style["jetz_mc"]["mpf1"] = kOpenTriangleUp;
  style["jetz_mc"]["mpfn"] = kOpenTriangleUp;
  style["jetz_mc"]["mpfu"] = kOpenTriangleDown;

  style["zjav"]["mpfchs1"] = kFullDiamond;
  style["zjav"]["ptchs"] = kOpenDiamond;
  style["zjav"]["mpf1"] = kFullTriangleUp;
  style["zjav"]["mpfn"] = kFullTriangleUp;
  style["zjav"]["mpfu"] = kFullTriangleDown;
  style["zjav_mc"]["mpf1"] = kOpenTriangleUp;
  style["zjav_mc"]["mpfn"] = kOpenTriangleUp;
  style["zjav_mc"]["mpfu"] = kOpenTriangleDown;
  
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

  style["gamjet"]["mpfchs1"] = kFullSquare;
  style["gamjet"]["ptchs"] = kOpenSquare;
  style["gamjet"]["mpf1"] = kFullTriangleUp;
  style["gamjet"]["mpfn"] = kFullTriangleUp;
  style["gamjet"]["mpfu"] = kFullTriangleDown;
  style["gamjet"]["rho"] = kFullTriangleDown;
  style["gamjet_mc"]["mpf1"] = kOpenTriangleUp;
  style["gamjet_mc"]["mpfn"] = kOpenTriangleUp;
  style["gamjet_mc"]["mpfu"] = kOpenTriangleDown;
  style["gamjet_mc"]["rho"] = kOpenTriangleDown;

  style["multijet"]["mpfchs1"] = kFullDiamond;
  style["multijet"]["ptchs"] = kOpenDiamond;
  style["multijet"]["chf"] = kFullCircle;
  style["multijet"]["nhf"] = kFullDiamond;
  style["multijet"]["nef"] = kFullSquare;
  style["multijet"]["cef"] = kFullDiamond;
  style["multijet"]["muf"] = kFullDiamond;
  style["multijet_mc"]["chf"] = kOpenCircle;
  style["multijet_mc"]["nhf"] = kOpenDiamond;
  style["multijet_mc"]["nef"] = kOpenSquare;
  style["multijet_mc"]["cef"] = kOpenDiamond;
  style["multijet_mc"]["muf"] = kOpenDiamond;
  style["multijet"]["mpf1"] = kFullTriangleUp;
  style["multijet"]["mpfn"] = kFullTriangleUp;
  style["multijet"]["mpfu"] = kFullTriangleDown;
  style["multijet"]["rho"] = kFullTriangleDown;
  style["multijet"]["npv"] = kFullTriangleDown;
  style["multijet_mc"]["mpf1"] = kOpenTriangleUp;
  style["multijet_mc"]["mpfn"] = kOpenTriangleUp;
  style["multijet_mc"]["mpfu"] = kOpenTriangleDown;
  style["multijet_mc"]["rho"] = kOpenTriangleDown;
  style["multijet_mc"]["npv"] = kOpenTriangleDown;

  map<string, int> color;
  color["jetz"] = kRed+1;
  color["jetz_mpf1"] = kRed;
  color["jetz_mpfn"] = kGreen+2;
  color["jetz_mpfu"] = kBlue;

  color["zjav"] = kRed+1;
  color["zjav_mpf1"] = kRed;
  color["zjav_mpfn"] = kGreen+2;
  color["zjav_mpfu"] = kBlue;
  
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

  style["gamjet"]["chf"] = kFullCircle;
  style["gamjet"]["nhf"] = kFullDiamond;
  style["gamjet"]["nef"] = kFullSquare;
  style["gamjet"]["cef"] = kFullDiamond;
  style["gamjet"]["muf"] = kFullDiamond;
  style["gamjet_mc"]["chf"] = kOpenCircle;
  style["gamjet_mc"]["nhf"] = kOpenDiamond;
  style["gamjet_mc"]["nef"] = kOpenSquare;
  style["gamjet_mc"]["cef"] = kOpenDiamond;
  style["gamjet_mc"]["muf"] = kOpenDiamond;

  color["gamjet"] = kBlue;
  color["gamjet_mpf1"] = kRed;
  color["gamjet_mpfn"] = kGreen+2;
  color["gamjet_mpfu"] = kBlue;
  color["gamjet_rho"] = kBlack;
  color["gamjet_chf"] = kRed;
  color["gamjet_nhf"] = kGreen+2;
  color["gamjet_nef"] = kBlue;
  color["gamjet_cef"] = kCyan+1;
  color["gamjet_muf"] = kMagenta+1;

  color["multijet"] = kRed+1;
  color["multijet_muf"] = kMagenta+1;
  color["multijet_mpf1"] = kRed;
  color["multijet_mpfn"] = kGreen+2;
  color["multijet_mpfu"] = kBlue;
  color["multijet_rho"] = kBlack;
  color["multijet_chf"] = kRed;
  color["multijet_nhf"] = kGreen+2;
  color["multijet_nef"] = kBlue;
  color["multijet_cef"] = kCyan+1;
  color["multijet_muf"] = kMagenta+1;

  style["incjet"]["chf"] = kFullCircle;
  style["incjet"]["nhf"] = kFullDiamond;
  style["incjet"]["nef"] = kFullSquare;
  style["incjet"]["cef"] = kFullDiamond;
  style["incjet"]["muf"] = kFullDiamond;

  color["incjet_chf"] = kRed;
  color["incjet_nhf"] = kGreen+2;
  color["incjet_nef"] = kBlue;
  color["incjet_cef"] = kCyan+1;
  color["incjet_muf"] = kMagenta+1;


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
  if (tepoch.Contains("UL") ||
      epoch=="RunCD" || 
      epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run22FG" || epoch=="Run22EFG" ||
      epoch=="Run23B" || epoch=="Run23BC123" ||  epoch=="Run23C123" ||
      //epoch=="Run23C" ||
      //epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
      epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" || epoch=="Run24E"||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
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
  types.push_back("crecoil");
  if (tepoch.Contains("UL") ||
      epoch=="RunCD" || 
      epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run22FG" || epoch=="Run22EFG" ||
      epoch=="Run23B" || epoch=="Run23BC123" ||  epoch=="Run23C123" ||
      //epoch=="Run23C" ||
      //epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
      epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" || epoch=="Run24E"||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
      epoch=="Run3") {
    types.push_back("rho");
  }
  //types.push_back("npv");

  // <pT,reco> and <pT,gen> vs ref pT (MC only)
  if (tepoch.Contains("UL") ||
      epoch=="RunCD" ||
      epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
      epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
      epoch=="Run22FG" || epoch=="Run22EFG" ||
      epoch=="Run23B" || epoch=="Run23BC123" ||  epoch=="Run23C123" ||
      //epoch=="Run23C" ||
      //epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
      epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
      epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" || epoch=="Run24E"||
      epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
      epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
      epoch=="Run3") {
    types.push_back("rjet");
    types.push_back("gjet");
  }

  vector<string> sets;
  sets.push_back("jetz");
  sets.push_back("zjav");
  sets.push_back("zjet");
  sets.push_back("gamjet");
  sets.push_back("multijet");
  sets.push_back("incjet");

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
	  if (!f) f = files[s+"_"+d];
	  if (!f) cout << "File not found for "<<s<<"_"<<d<<endl<<flush;
	  assert(f);

	  // C_recoil is only meant for multijets
	  if (t=="crecoil" && s!="multijet") continue;

	  bool ismpfc = (t=="mpf1" || t=="mpfn" || t=="mpfu");// || t=="rho");
	  bool isfrac = (t=="chf"||t=="nef"||t=="nhf"||
			 t=="cef"||t=="muf"||t=="puf");
	  if (s=="incjet" && !isfrac) continue;

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
	  //if (isflavormc && d=="multijet") continue;

	  // true responses only for gamjet for now => add to Z+jet
	  if ((t=="rjet" || t=="gjet") && d!="mc") continue;
	  if ((t=="rjet" || t=="gjet") && s=="multijet") continue;
	  if ((t=="rjet" || t=="gjet") && s=="gamjet") continue;

	  // jetz only for HDM inputs
	  if (s=="jetz" && !(t=="mpfchs1" || t=="ptchs" || ismpfc ||
			     t=="counts")) continue;
	  // zjav only for HDM inputs
	  if (s=="zjav" && !(t=="mpfchs1" || t=="ptchs" || ismpfc ||
			     t=="counts")) continue;

	  double alpha = 1.00;
	  eta1 = etas[ieta].first; // reset to avoid trouble with below
	  eta2 = etas[ieta].second; // reset to avoid trouble with below
	  
	  // Reconstruct naming scheme used in each of the files
	  // If non-conventional naming schemes, patch here
	  const char *c(0);
	  if (s=="jetz") {
	    c = Form("%s/eta_%02.0f_%02.0f/%s_jetpt_a%1.0f",
		     rename["zjet"][d],10*eta1,10*eta2,
		     rename["zjet"][t],100.*alpha);
	  } // "jetz"
	  if (s=="zjav") {
	    c = Form("%s/eta_%02.0f_%02.0f/%s_ptave_a%1.0f",
		     rename["zjet"][d],10*eta1,10*eta2,
		     rename["zjet"][t],100.*alpha);
	  } // "zjav"
	  if (s=="zjet") {
	    if (tepoch.Contains("UL") ||
		epoch=="RunCD" ||
		epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD" ||
		epoch=="Run22E" || epoch=="Run22F" || epoch=="Run22G" ||
		epoch=="Run22FG" || epoch=="Run22EFG" ||
		epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C123" ||
		//epoch=="Run23C" ||
		//epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
		epoch=="Run23C4" || epoch=="Run23D" || epoch=="Run23C4D" ||
		epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" ||
		epoch=="Run24E" ||
		epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE" ||
		epoch=="Run24CR" || epoch=="Run24CS" || epoch=="Run24CP" ||
		epoch=="Run3"
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
	    if (epoch=="Run3" && true) {
	      c = Form("%s/eta00-13/%s_%s_a100",d.c_str(),t.c_str(),s.c_str());
	    }
	  } // "zjet"
	  if (s=="gamjet" && t=="counts") {
	    c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f_%s",
		     rename[s]["mpfchs1"],
		     d=="ratio" ? rename[s]["data"] : rename[s][d],
		     100.*alpha, 10.*eta1, 10.*eta2, rename[s][t]);
	  }
	  //else if (s=="gamjet" && (ismpfc||t=="mpfchs1"||t=="ptchs") &&
	  //	   d=="ratio") { // patch missing => in v21 now
	  //  c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",
	  //           rename[s][t], rename[s]["data"], // ratio->data
	  //           100.*alpha, 10.*eta1, 10.*eta2);
	  //} // gamjet
	  else if (s=="gamjet" && isfrac) { // new PF composition
	    c = Form("pf/p%s_%s",tt,dd);
	  }
	  else if (s=="gamjet") {
	    c = Form("%s%s_a%1.0f_eta%02.0f_%02.0f",//%s",
		     rename[s][t], rename[s][d],
		     100.*alpha, 10.*eta1, 10.*eta2);//,
	    //(d!="data" ? sPSWgtG.c_str() : ""));
	  } // gamjet
	  if (s=="multijet") {
	    // Override file selection from files[s]
            //if (d=="data") f = fmjd;
            //if (d=="mc") f = fmjm;
            //if (d=="ratio") f = fmjd; // patch
	    c = Form("Multijet/%s",rename[s][t]);
	    //if (d=="mc") // patch 22Sep2023, not 19Dec2023
	    //c = Form("HLT_MC/Multijet/%s",rename[s][t]);
	  }
	  if (s=="incjet") {
	    // Override file selection from files[s]
            //if (d=="data") f = fijd;
            //if (d=="mc") f = fijm;
            //if (d=="ratio") f = fijd; // patch
	    c = Form("Incjet/PFcomposition/p%s13",tt);//rename[s][t]);
	    //if (d=="mc") // patch 22Sep2023, not 19Dec2023
	    //c = Form("HLT_MC/Incjet/PFcomposition/p%s13",tt);//rename[s][t]);
	  }
	  
	  assert(f);
	  TObject *obj = f->Get(c);
	  if (!obj) {
	    cout << "Graph " << c << " not found for "
		 << s << " " << t << " " << d << "!" <<endl << flush;
	    cout << "File: " << f->GetName() << endl << flush;
	    cout << "Eta " << eta1 << " - " << eta2 <<  endl << flush;
	  }
	  if (t=="counts" && !obj) obj = hcounts;
	  
	  // Calculate data/MC ratio
	  if ((s=="jetz" || s=="zjav" || s=="zjet" || sp=="zjet" ||
	       s=="multijet") &&
	      d=="ratio" && (t=="mpfchs1"||t=="ptchs")) {
	    TGraphErrors *gd = grs["data"][t][s][ieta];
	    TGraphErrors *gm = grs["mc"][t][s][ieta];
	    assert(gd);
	    assert(gm);
	    
	    TGraphErrors *g = tools::ratioGraphs(gd, gm);
	    obj = (TObject*)g;
	  }
	  if ((s=="jetz" || s=="zjav" || s=="zjet" ||
	       s=="multijet" || s=="incjet") &&
	      d=="ratio" &&
	      // bug: 20231114...
	      //(isfrac || t=="rmpf1" || t=="rmpfn" || t=="rpmfu")) {
	      (isfrac || ismpfc)) {
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

	  // Rebin multijets
	  if (s=="multijet" && (d=="data" || d=="mc") &&
	      (t=="mpfchs1" || t=="ptchs" || ismpfc) &&
	      rebinMultijet) {

	    Double_t vx[] = 
	      {1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49,
	       56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245,
	       272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686,
	       737, 790, 846, 905, 967, 1032,
	       // 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
	       // v1 binning
	       1248, 1497, 1784, 2116, 2500, 3273,  4252, 5492, 5777, 7000};
	    // Photon+jet
  	    // 1000, 1200, 1450, 1750, 2100, 2500, 3000};
	    
	    const int nx = sizeof(vx)/sizeof(vx[0])-1;
	    
	    if (obj->InheritsFrom("TProfile")) {
	      obj = ((TProfile*)obj)->Rebin(nx,Form("%s_%s_%s",ss,dd,tt),&vx[0]);
	    }
	  } // rebinMultijet
	    
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
	    //if (g->GetEY()[i]>0.2 && !isflavormc)  g->RemovePoint(i);
	    // Clean out point outside good ranges
	    //else
	    if (s=="jetz" &&
		(g->GetX()[i]<fjzptmin || g->GetX()[i]>fjzptmax))
	      g->RemovePoint(i);
	    else if (s=="zjav" &&
		     (g->GetX()[i]<fzjaptmin || g->GetX()[i]>fzjaptmax))
	      g->RemovePoint(i);
	    else if (s=="zjet" &&
		     (g->GetX()[i]<fzptmin || g->GetX()[i]>fzptmax))
	      g->RemovePoint(i);
	    else if (s=="gamjet" &&
		     (g->GetX()[i]<fgptmin || g->GetX()[i]>fgptmax))
	      g->RemovePoint(i);
	    else if (s=="multijet" &&
		     (g->GetX()[i]<fmjptmin || g->GetX()[i]>fmjptmax))
	      g->RemovePoint(i);
	    else if (s=="incjet" &&
		     (g->GetX()[i]<fijptmin || g->GetX()[i]>fijptmax))
	      g->RemovePoint(i);
	  } // for i
	  
	  // Mass corrections for Z+jet
	  if (correctZMass && (s=="jetz" || s=="zjav" || s=="zjet") &&
	      (d=="data" || d=="ratio") &&
	      (ismpfc || t=="ptchs" || t=="mpfchs1"))  { // 20240227 add mpfchs1
	    assert(false); // should not do this for mpu, mpfn?
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

	  // PATCH EM scale corrections for gamma+jet in absence of QCD bkg
	  if (s=="gamjet" && (d=="data" || d=="ratio") &&
	      (ismpfc || t=="ptchs" || t=="mpfchs1"))  { // 20240227 add mpfchs1
	    for (int i = 0; i != g->GetN(); ++i) {
	      double pt = g->GetX()[i];
	      if (scaleEMperEra) {
		scaleEM = 1;
		if (epoch=="Run23C123") scaleEM = 0.996;//0.998;
		if (epoch=="Run23C4")   scaleEM = 0.994;
		if (epoch=="Run23D")    scaleEM = 0.990;
	      }
	      g->SetPoint(i, g->GetX()[i], g->GetY()[i]*scaleEM);
	    }
	  }
	  // PATCH

	  // Temporary PATCH for multijet crecoil
	  //if (s=="multijet" && t=="crecoil") {
	  //for (int i = 0; i != g->GetN(); ++i) {
	  //  double c = g->GetY()[i];
	  //  if (fabs(c-1)<0.01) {
	  //	g->SetPoint(i, g->GetX()[i], 0.40);
	  //  }
	  //}
	  //}
	  // Temporary PATCH

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
    if (tepoch.Contains("UL")) {
      if (epoch=="2016ULAPV") {
	jec = getFJC("","","2016BCDEF_DATA_L2L3Residual_AK4PFchs","textFiles");
	mcjec = getFJC("","Summer19UL16APV_RunBCDEF_V7_DATA_L2Relative_AK4PFchs");
      }
      if (epoch=="2016UL") {
	jec = getFJC("","","Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs");
	//jec = getFJC("","","2016BCDEF_DATA_L2L3Residual_AK4PFchs","textFiles");
	jec = getFJC("","","2016GH_DATA_L2L3Residual_AK4PFchs","textFiles");
	mcjec = getFJC("","Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs");
      }
      if (epoch=="2017UL") {
	//jec = getFJC("","","Summer19UL17_RunE_V6_DATA_L2L3Residual_AK4PFchs");
	jec = getFJC("","","2017BCDEF_DATA_L2L3Residual_AK4PFchs","textFiles");
	mcjec = getFJC("","Summer19UL17_RunE_V6_DATA_L2Relative_AK4PFchs");
      }
      if (epoch=="2018UL") {
	//jec = getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs");
	jec = getFJC("","","2018ABCD_DATA_L2L3Residual_AK4PFchs","textFiles");
	mcjec = getFJC("","Summer19UL18_RunD_V5_DATA_L2Relative_AK4PFchs");
      }
    }
    if (epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD") {
      /*
      jec = getFJC("","","Winter22Run3_RunC_V2_DATA_L2L3Residual");
      mcjec = getFJC("","Winter22Run3_RunC_V2_DATA_L2Relative");
      */
      jec = getFJC("","","Summer22-22Sep2023_Run2022CD_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22Run3_V1_MC_L2Relative");
    }
    if (epoch=="Run22E") {
      jec = getFJC("","","Summer22EE-22Sep2023_Run2022E_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22EEVetoRun3_V1_MC_L2Relative");
    }
    if (epoch=="Run22F" || epoch=="Run22FG") {
      // because Sami used these for Zb (up to v59 at least?)
      /*
      jec = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
      mcjec = getFJC("","Winter23Prompt23_RunA_V1_DATA_L2Relative");
      */
      jec = getFJC("","","Summer22EEPrompt22_Run2022F_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22EEVetoRun3_V1_MC_L2Relative");
    }
    if (epoch=="Run22G") {
      jec = getFJC("","","Summer22EEPrompt22_Run2022G_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22EEVetoRun3_V1_MC_L2Relative");
    }
    if (epoch=="Run22EFG") { assert(false); exit(0); }
    if (epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C123") {
      /*
      jec = getFJC("","","Summer22Prompt23_Run2023Cv123_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22Run3_V1_MC_L2Relative");
      */
    //jec = getFJC("","","Summer22Prompt23_Run2023Cv123_V3_DATA_L2L3Residual");
      jec = getFJC("","","Summer23Prompt23_Run2023Cv123_V2_DATA_L2L3Residual");
      mcjec = getFJC("","Summer23Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    //epoch=="Run23C" ||
    //epoch=="Run23C1" || epoch=="Run23C2" || epoch=="Run23C3" ||
    if (epoch=="Run23C4") {
      /*
      jec = getFJC("","","Summer22Prompt23_Run2023Cv4_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22Run3_V1_MC_L2Relative");
      */
      //jec = getFJC("","","Summer22Prompt23_Run2023Cv4_V3_DATA_L2L3Residual");
      jec = getFJC("","","Summer23Prompt23_Run2023Cv4_V2_DATA_L2L3Residual");
      mcjec = getFJC("","Summer23Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run23D") {
      /*
      jec = getFJC("","","Summer22Prompt23_Run2023D_V3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22Run3_V1_MC_L2Relative");
      */
      //jec = getFJC("","","Summer22Prompt23_Run2023D_V3_DATA_L2L3Residual");
      // also 2024BC_V1M
      jec = getFJC("","","Summer23Prompt23_Run2023D_V2_DATA_L2L3Residual");
      mcjec = getFJC("","Summer23BPixRun3_V3_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" ||
	epoch=="Run24E" ||
	epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE") {
    //jec = getFJC("","","Prompt24_Run2024BC_V1M_DATA_L2L3Residual_AK4PFPuppi");
      jec = getFJC("","","Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi");
      mcjec = getFJC("","Winter24Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run24CR") {
      //jec = getFJC("","","Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi");
      jec = getFJC("","","Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi");
      mcjec = getFJC("","Winter24Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run24CS") {
      jec = getFJC("","","Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi");
      mcjec = getFJC("","Winter24Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run24CP") {
      jec = getFJC("","","Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi");
      mcjec = getFJC("","Winter24Run3_V1_MC_L2Relative_AK4PUPPI");
    }
    if (epoch=="Run23C4D") { assert(false); exit(0); }
    if (epoch=="Run3") {
      //jec = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
      //mcjec = getFJC("","Winter23Prompt23_RunA_V1_DATA_L2Relative");
      // Combo done with jecdataRun3Data.root and createL2L3ResTextFile.C
      jec = getFJC("","","Summer22Prompt23_Run3_reV3_DATA_L2L3Residual");
      mcjec = getFJC("","Summer22Run3_V1_MC_L2Relative");
    }
    assert(jec);
    assert(mcjec);

    //#define PAIR(a,b) (make_pair<double,FactorizedJetCorrector*>((a),getFJC("","",(b))))
    //add vjec here later to combine IOVs

    if (rp_debug) cout << "Loading reference JECs..." << endl << flush;
    
    FactorizedJetCorrector *jecrun1, *jecrun2, *jecold(0);

    // Reference Run I and Run II JEC for plotting on the back
    jecrun1 = getFJC("","","Winter14_V8_DATA_L2L3Residual_AK5PFchs"); 
    jecrun2 = getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs"); 
    
    // Store old JEC for undoing it in global fit (JEC from closure files)
    // NB: I think 'jec' is now used instaed of 'jecold' so this may be obsolete
    //jecold = getFJC("","",Form("Winter22Run3_Run%s_V1_DATA_L2L3Residual","A"));
    if (tepoch.Contains("UL")) {
      if (epoch=="2016ULAPV") {
	//jecold = getFJC("","","Summer19UL16_RunBCDEF_V7_DATA_L2L3Residual_AK4PFchs");
	jecold = getFJC("","","2016BCDEF_DATA_L2L3Residual_AK4PFchs","textFiles");
      }
      if (epoch=="2016UL") {
	//jecold = getFJC("","","Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs");
	jecold = getFJC("","","2016GH_DATA_L2L3Residual_AK4PFchs","textFiles");
      }
      if (epoch=="2017UL") {
	//jecold = getFJC("","","Summer19UL17_RunE_V6_DATA_L2L3Residual_AK4PFchs");
	jecold = getFJC("","","2017BCDEF_DATA_L2L3Residual_AK4PFchs","textFiles");
      }
      if (epoch=="2018UL") {
	//jecold = getFJC("","","Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs");
	jecold = getFJC("","","2018ABCD_DATA_L2L3Residual_AK4PFchs","textFiles");
      }
    }

    if (epoch=="Run22C" || epoch=="Run22D" || epoch=="Run22CD") {
      //jecold = getFJC("","","Winter22Run3_RunC_V2_DATA_L2L3Residual");
      jecold = getFJC("","","Summer22-22Sep2023_Run2022CD_V3_DATA_L2L3Residual");
    }
    if (epoch=="Run22E") {
      jecold = getFJC("","","Summer22EE-22Sep2023_Run2022E_V3_DATA_L2L3Residual");
    }
    if (epoch=="Run22F" || epoch=="Run22FG") {
      // because Sami used this (at least up to v59?)
      //jecold = getFJC("","","Winter23Prompt23_RunA_V1_DATA_L2L3Residual");
      jecold = getFJC("","","Summer22EEPrompt22_Run2022F_V3_DATA_L2L3Residual");
    }
    if (epoch=="Run22F" || epoch=="Run22FG") {
      jecold = getFJC("","","Summer22EEPrompt22_Run2022G_V3_DATA_L2L3Residual");
    }
    if (epoch=="Run23B" || epoch=="Run23BC123" || epoch=="Run23C123") {
      //jecold = getFJC("","","Summer22Prompt23_Run2023Cv123_V3_DATA_L2L3Residual");
      jecold = getFJC("","","Summer23Prompt23_Run2023Cv123_V2_DATA_L2L3Residual");
    }
    if (epoch=="Run23C4") {
      //jecold = getFJC("","","Summer22Prompt23_Run2023Cv4_V3_DATA_L2L3Residual");
      jecold = getFJC("","","Summer23Prompt23_Run2023Cv4_V2_DATA_L2L3Residual");
    }
    if (epoch=="Run23D") {
      // also 2023BC_V1M
      //jecold = getFJC("","","Summer22Prompt23_Run2023D_V3_DATA_L2L3Residual");
      jecold = getFJC("","","Summer23Prompt23_Run2023D_V2_DATA_L2L3Residual");
    }
    if (epoch=="Run24B" || epoch=="Run24C" || epoch=="Run24D" ||
	epoch=="Run24E" ||
	epoch=="Run24BC" || epoch=="Run24BCD" || epoch=="Run24BCDE") {
    //jecold = getFJC("","","Prompt24_Run2024BC_V1M_DATA_L2L3Residual_AK4PFPuppi");
      //jecold = getFJC("","","Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi");
      jecold = getFJC("","","Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi");
    }
    if (epoch=="Run24CR") {
      //jecold = getFJC("","","Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi");
      jecold = getFJC("","","Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi");
    }
    if (epoch=="Run24CS") {
      jecold = getFJC("","","Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi");
    }
    if (epoch=="Run24CP") {
      jecold = getFJC("","","Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi");
    }
    if (epoch=="Run3") {
      // Combo done with jecdataRun3Data.root and createL2L3ResTextFile.C
      jecold = getFJC("","","Summer22Prompt23_Run3_reV3_DATA_L2L3Residual");
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
  jec->setJetPhi(0.); // Summer23 addition

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

void setEpoch(string epoch) {

  if (TString(epoch).Contains("UL")) {
    fzptmin = 15;
  }
  //if (epoch=="Run23C4") {
  //fzptmin = 20.;
  //}
  //else
  return;
}
