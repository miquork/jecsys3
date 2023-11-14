// Purpose: Create L2L3Res text file with simple parameterization
//          Takes as input previous L2Res and complex 9p global JES fit
//          Outputs same L2Res and "simple" 7p->8p fit
//          ("simple" as in removing main parameter degeneracies)
//          For merging IOVs together at text file level, use
//          [minitools/mergeL2L3ResTextFiles.C]
#include "TString.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TLine.h"

#include "tdrstyle_mod22.C"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const bool debug = true;

void createL2L3ResTextFiles(string set);

TLegend *_leg(0);
void createL2L3ResTextFile() {

  setTDRStyle();

  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    0.90+1e-4,1.30-1e-4,"p_{T} (GeV)",ptmin,ptmax);
		    //0.88+1e-4,1.05-1e-4,"p_{T} (GeV)",ptmin,ptmax);
  lumi_136TeV = "Run3, 63 fb^{-1}"; // Not including 23B
 TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  c1->SetLeftMargin(0.17);
  c1->SetRightMargin(0.03);
  h->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  _leg = tdrLeg(0.45,0.90,0.75,0.90);

  //createL2L3ResTextFiles("RunC"); // crash
  //createL2L3ResTextFiles("RunCD");
  //createL2L3ResTextFiles("Ref");

  //createL2L3ResTextFiles("Run22C");
  //createL2L3ResTextFiles("Run22D");
  //createL2L3ResTextFiles("Run22E");
  //createL2L3ResTextFiles("Run22F");
  //createL2L3ResTextFiles("Run22G");


  //createL2L3ResTextFiles("Run23B"); // empty
  //createL2L3ResTextFiles("Run23C1");
  //createL2L3ResTextFiles("Run23C2");
  //createL2L3ResTextFiles("Run23C3");
  //createL2L3ResTextFiles("Run23C4");

  h->SetMinimum(0.70+1e-4);
  createL2L3ResTextFiles("Run22CD-22Sep2023");
  createL2L3ResTextFiles("Run22E-22Sep2023");
  //createL2L3ResTextFiles("Run22FG-Prompt");
  createL2L3ResTextFiles("Run22F-Prompt");
  createL2L3ResTextFiles("Run22G-Prompt");
  createL2L3ResTextFiles("Run23C123-Prompt");
  createL2L3ResTextFiles("Run23C4D-Prompt");
  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Run22CDE-22Sep2023_22FG_23BC-Prompt.pdf");
  
  //c1->SaveAs("pdf/createL2L3ResTextFile_2022C.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run22F.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run22FG_prompt.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run22FG_Run23BC_prompt.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run22CDEFG_23BC_prompt.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Run22CDEFG_23BC_prompt.pdf");


  /*
  // Produce Run2 year-averages for reprocess.C reference (jec, jecold)
  // Use files in jecsys2020 directory, not jecsys3
  h->GetYaxis()->SetRangeUser(0.94+1e-4,1.06-1e-4);
  createL2L3ResTextFiles("2016BCDEF");
  createL2L3ResTextFiles("2016GH");
  createL2L3ResTextFiles("2017BCDEF");
  createL2L3ResTextFiles("2018ABCD");
  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Run2UL_JECV7.pdf");
  */
}

void createL2L3ResTextFiles(string set) {

  //if (debug) 
  cout << "****************************************************************\n";
  cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";
  cout << "****************************************************************\n";

  cout << "** Processing " << set << " **" << endl << flush;
  cout << "******************************\n";

  // Simplify complex sum into an effective formula
  // Need good starting values and/or a few iterations to converge
  // Fit done to hjesfit from each IOV
  TDirectory *curdir = gDirectory;
  TFile *f(0);
  bool isRun3(false);
  //TFile *f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");
  if (set=="Ref") { // reference JEC
    f = new TFile("rootfiles/jecdataRunCD_v6.root","READ");
  }
  else if (set=="2016BCDEF") {
    f = new TFile("../jecsys2020/rootfiles/jecdata2016BCDEF.root","READ");
  }
  else if (set=="2016GH") {
    f = new TFile("../jecsys2020/rootfiles/jecdata2016GH.root","READ");
  }
  else if (set=="2017BCDEF") {
    f = new TFile("../jecsys2020/rootfiles/jecdata2017BCDEF.root","READ");
  }
  else if (set=="2018ABCD") {
    f = new TFile("../jecsys2020/rootfiles/jecdata2018ABCD.root","READ");
  }
  else if (set=="Run22CD-22Sep2023") {
    f = new TFile("rootfiles/jecdataRun22CD.root","READ"); isRun3=true;
  }
  else if (set=="Run22E-22Sep2023") {
    f = new TFile("rootfiles/jecdataRun22E.root","READ"); isRun3=true;
  }
  else if (set=="Run22F-Prompt" || set=="Run22G-Prompt" ||
	   set=="Run22FG-Prompt") {
    f = new TFile("rootfiles/jecdataRun22FG.root","READ"); isRun3=true;
  }
  else if (set=="Run23C123-Prompt") {
    f = new TFile("rootfiles/jecdataRun23C123.root","READ"); isRun3=true;
  }
  else if (set=="Run23C4D-Prompt") {
    f = new TFile("rootfiles/jecdataRun23C4D.root","READ"); isRun3=true;
  }
  else
    f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");

  /*
  if (set=="Run22F") {
    f = new TFile(Form("../JERCProtoLab/textFiles/Summer22EERun3_RunF_V2_DATA/jecdata%s.root",set.c_str()),"READ");
  }
  if (set=="Run22G") {
    f = new TFile(Form("../JERCProtoLab/textFiles/Summer22EERun3_RunG_V2_DATA/jecdata%s.root",set.c_str()),"READ");
  }
  if (set=="Run22G") {
    f = new TFile(Form("../JERCProtoLab/textFiles/Summer22EERun3_RunG_V2_DATA/jecdata%s.root",set.c_str()),"READ");
  }

  if (set=="Run23BC123") {
    f = new TFile(Form("../JERCProtoLab/textFiles/Winter23Prompt23_RunBC123_V2_DATA/jecdata%s.root",set.c_str()),"READ");
  }
  if (set=="Run23C4") {
    f = new TFile(Form("../JERCProtoLab/textFiles/Winter23Prompt23_RunC4_V2_DATA/jecdata%s.root",set.c_str()),"READ");
  }
  */  
  assert(f && !f->IsZombie());
  TH1D *h(0);
  //h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit");
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2"); // Run2 refit
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/herr_l2l3res"); // Run2 ref. JES
  assert(h);
  curdir->cd();

  TF1 *f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x",15,4500);
  f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001, 0);
  if (set=="2017H") 
    f1->SetParameters(0.99, 1.5,0.01, 0.01,1000.,1.3, 0.001, 0.);
  if (!isRun3)
    f1->FixParameter(7, 0.);
  
  // To avoid division by zero errors
  f1->SetParLimits(4,10.,6500.);
  f1->SetParLimits(5,0.,10.);
  
  map<string,int> color;
  color["2016BCDEF"] = kYellow+2;
  color["2016GH"] = kRed;
  color["2017BCDEF"] = kGreen+2;
  color["2018ABCD"] = kBlue;

  color["RunC"] = kYellow+2;
  color["RunCD"] = kOrange+1;
  color["Ref"] = kBlue;
  /*
  color["Run22F"] = kOrange+1;
  color["Run22G"] = kRed;
  color["Run23BC123"] = kGreen+2;
  color["Run23C4"] = kBlue;
  */

  color["Run22CD"] = kGreen+4;//+2;
  color["Run22C"] = kGreen+1;
  color["Run22D"] = kGreen+2;
  color["Run22E"] = kGreen+3;//kMagenta+1;
  color["Run22FG"] = kBlue;
  color["Run22F"] = kCyan+2;
  color["Run22G"] = kCyan+3;

  color["Run23BC123"] = kRed;
  color["Run23B"] = kRed-6;
  color["Run23C1"] = kOrange;
  color["Run23C2"] = kOrange+1;
  color["Run23C3"] = kOrange+2;
  color["Run23C4"] = kMagenta+2;//kGray+2;

  /*
  // Old colors
  color["Run22CD-22Sep2023"] = kGreen+2;
  color["Run22E-22Sep2023"] = kGreen+3;//kCyan+2;
  color["Run22FG-Prompt"] = kCyan+2;//kMagenta+3;
  color["Run22F-Prompt"] = kCyan+2;//kMagenta+2;
  color["Run22G-Prompt"] = kCyan+3;//kMagenta+1;
  color["Run23C123-Prompt"] = kRed;
  color["Run23C4D-Prompt"] = kMagenta+2;//kBlue;
  */
  // Newer colors
  color["Run22CD-22Sep2023"] = kGreen+2;//kBlue;
  color["Run22E-22Sep2023"] = kCyan+2;
  color["Run22FG-Prompt"] = kRed;
  color["Run22F-Prompt"] = kRed;
  color["Run22G-Prompt"] = kRed+2;
  color["Run23C123-Prompt"] = kOrange+1;
  color["Run23C4D-Prompt"] = kBlue;//kGreen+2;
  
  h->Fit(f1,"QRN");
  h->Fit(f1,"QRNM");
  h->Fit(f1,"QRNM");
  tdrDraw(h,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  h->SetFillColorAlpha(color[set]-9,0.7);
  f1->SetLineColor(color[set]);
  f1->Draw("SAME");

  _leg->SetY1(_leg->GetY1()-0.05);
  if (set=="Ref")
    _leg->AddEntry(h,"Old 22C (ref.)","FL");
  else
    _leg->AddEntry(h,set.c_str(),"FL");

  const int np = 8;//7;
  double p[np];
  for (int i = 0; i != np; ++i) {
    p[i] = f1->GetParameter(i);
  } // for i


  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  const char *run = set.c_str();

  string sin, sout;
  if (set=="2016BCDEF") {
    sin = "CondFormats/JetMETObjects/data/Summer19UL16APV_RunBCDEF_V7_DATA_L2L3Residual_AK4PFchs.txt";
  }
  if (set=="2016GH") {
    sin = "CondFormats/JetMETObjects/data/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs.txt";
  }
  if (set=="2017BCDEF") {
    sin = "CondFormats/JetMETObjects/data/Summer19UL17_RunE_V6_DATA_L2L3Residual_AK4PFchs.txt";
  }
  if (set=="2018ABCD") {
    sin = "CondFormats/JetMETObjects/data/Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs.txt";
  }
  if (set=="RunC" || set=="RunCD" || set=="Ref")
    sin = "../JERCProtoLab/Winter22Run3/L2Residual/Winter22Run3_V1/Winter22Run3_V1_MPF_LOGLIN_L2Residual_pythia8_AK4PFPuppi.txt"; // For V1 L2L3Res
  //sout = "../JERCProtoLab/Winter22Run3/global_fit/Winter22Run3_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V1 L2L3Res
  //sin = "../JERCProtoLab/Winter22Run3/L2Residual/WinterRun3_V2_PtJER/Winter22RunC_V2_MPF_LOGLIN_L2Residual_pythia8_AK4PFPuppi.txt"; // For V2 L2L3Res
  //sout = "../JERCProtoLab/Winter22Run3/global_fit/Winter22Run3_RunC_V2_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V2 L2l3Res (L2Res RunC, L3Res RunCD)
  //
  if (set=="Run22C" || set=="Run22D" || set=="Run22CD") {
    //sin = "CondFormats/JetMETObjects/data/Winter22Run3_V1_MC_L2Residual_AK4PFPuppi.txt";
    sin = "CondFormats/JetMETObjects/data/Winter22Run3_RunC_V2_DATA_L2Residual_AK4PFPuppi.txt";
  }
  if (set=="Run22E" || set=="Run22F" || set=="Run22G" || set=="Run22FG") {
    sin = "CondFormats/JetMETObjects/data/Winter23Prompt23_RunA_V1_DATA_L2Residual_AK4PFPuppi.txt";
  }
  if (set=="Run23B" || set=="Run23C1" || set=="Run23C2" || set=="Run23C3" || 
      set=="Run23BC123" || set=="Run23C4") {
    sin = "CondFormats/JetMETObjects/data/Winter23Prompt23_RunA_V1_DATA_L2Residual_AK4PFPuppi.txt";
  }
  if (set=="2016BCDEF" || set=="2016GH" || set=="2017BCDEF" || set=="2018ABCD")
   sout = Form("textfiles/%s_DATA_L2L3Residual_AK4PFchs.txt",run);
  else   
    //sout = Form("textfiles/Sami_20230630/%s_DATA_L2L3Residual_AK4PFPuppi.txt",run);
    sout = Form("textfiles/Run3partial/%s_DATA_L2L3Residual_AK4PFPuppi.txt",run);
  /*
  if (set=="Run22F") {
    sin = "../JERCProtoLab/textFiles/Summer22EERun3_RunF_V2_DATA/Summer22EERun3_RunF_V2_DATA_L2Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
    sout = "../JERCProtoLab/textFiles/Summer22EERun3_RunF_V2_DATA/Summer22EERun3_RunF_V2_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
  }
  if (set=="Run22G") {
    sin = "../JERCProtoLab/textFiles/Summer22EERun3_RunG_V2_DATA/Summer22EERun3_RunG_V2_DATA_L2Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
    sout = "../JERCProtoLab/textFiles/Summer22EERun3_RunG_V2_DATA/Summer22EERun3_RunG_V2_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
  }
  if (set=="Run23BC123") {
    sin = "../JERCProtoLab/textFiles/Winter23Prompt23_RunBC123_V2_DATA/Winter23Prompt23_RunA_V1_DATA_L2Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
    sout = "../JERCProtoLab/textFiles/Winter23Prompt23_RunBC123_V2_DATA/Winter23Prompt23_RunBC_V2_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
  }
  if (set=="Run23C4") {
    sin = "../JERCProtoLab/textFiles/Winter23Prompt23_RunC4_V2_DATA/Winter23Prompt23_RunA_V1_DATA_L2Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
    sout = "../JERCProtoLab/textFiles/Winter23Prompt23_RunC4_V2_DATA/Winter23Prompt23_RunC4_V2_DATA_L2L3Residual_AK4PFPuppi.txt"; // For V2 L2L3Res
  }
  */

  bool isL2Res(false);
  if (set=="Run22CD-22Sep2023") {
    sin = "CondFormats/JetMETObjects/data/Summer22_RunCD_V2_MPF_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
  }
  if (set=="Run22E-22Sep2023") {
    sin = "CondFormats/JetMETObjects/data/Summer22EE_RunE_V2_MPF_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
  }
  if (set=="Run22FG-Prompt") {
    sin = "CondFormats/JetMETObjects/data/Summer22EEPrompt22_RunFG_V2_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
    assert(false); // 22FG_L2Res not yet available
  }
  if (set=="Run22F-Prompt") {
    sin = "CondFormats/JetMETObjects/data/Summer22EEPrompt22_RunF_V2_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
  }
  if (set=="Run22G-Prompt") {
    sin = "CondFormats/JetMETObjects/data/Summer22EEPrompt22_RunG_V2_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
  }
  if (set=="Run23C123-Prompt" || set=="Run23C4D-Prompt") {
    sin = "CondFormats/JetMETObjects/data/Summer22EEPrompt22_RunG_V2_L2Residual_AK4PFPuppi.txt"; // Placeholder
    isL2Res = true;
  }

  cout << "Reading in L2Residual file:" << endl
       << "   " << sin << endl << flush;  
  ifstream fin(sin.c_str());
  assert(fin.is_open());
  
  cout << "Writing out L2L3Residual file:" << endl
       << "   " << sout << endl << flush;
  ofstream fout(sout.c_str());

  // 2018: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018

  string header;
  getline(fin, header);
  if (debug) cout << "Old L2L3Residual header:" << endl;
  if (debug) cout << header << endl;

  //header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}";
  if (isRun3) {
    //Old L2Residual header: [2]*([3]*([4]+TMath::Log(max([0],min([1],x)))*([5]+TMath::Log(max([0],min([1],x)))*[6])+[7]/x))
    header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+TMath::Log(max([0],min([1],x)))*([5]+TMath::Log(max([0],min([1],x)))*[6])+[7]/x))*1./([8]+[9]/x+[10]*log(x)/x+[11]*(pow(x/[12],[13])-1)/(pow(x/[12],[13])+1)+[14]*pow(x,-0.3051)+[15]*x) Correction L2Relative}";
  }
  else {
    header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}";
  } 
  //([6]+[7]*100./3.*(TMath::Max(0.,1.03091-0.051154*pow(x,-0.154227))-TMath::Max(0.,1.03091-0.051154*TMath::Power(208.,-0.154227)))+[8]*0.021*(-1.+1./(1.+exp(-(TMath::Log(x)-5.030)/0.395))))) Correction L2Relative}
  if (debug) cout << "New L2L3Residual header:" << endl;
  if (debug) cout << header << endl;
  fout << header << endl;
  
  string line;
  double etamin, etamax;
  int npar, xmin, xmax, ptmin0, ptmax1;
  double p2, p3, p4, p5;
  double p6, p7(0), p8(0);
  double p9(0), p10(0); // isRun3
  int cnt(0); int cntmax(0);
  while (getline(fin,line)) {
    if (cnt<cntmax && debug) cout << line << endl;
    if (isRun3 && isL2Res) {
      assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %d %d  %lf %lf %lf %lf"
		    "  %lf %lff",
		    &etamin, &etamax, &npar, &xmin, &xmax, &ptmin0, &ptmax1,
		    &p2, &p3, &p4, &p5,  &p6, &p7)==13);
    }
    else if (isRun3 && !isL2Res) {
      assert(false);
    }
    else {
      assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %d %d  %lf %lf %lf %lf"
		    "  %lf %lf %lf",
		    &etamin, &etamax, &npar, &xmin, &xmax, &ptmin0, &ptmax1,
		    &p2, &p3, &p4, &p5,  &p6, &p7, &p8)==14);
      assert(false); // takin L2Residual only for 22Sep2023
    }
    assert(!(etamin==0 && etamax==0));
    if (fabs(etamin)<0.01 && fabs(etamax)<0.01) {
      cout << "sscanf failed!" << endl << flush;
      exit(1);
    }

    if (isRun3) {
      int nparnew = 18;//17;
      if (cnt<cntmax && debug)
	cout << Form("  %9.6f %9.6f   %d   %d %d"
		     "   %d   %d"
		     "   %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"
		     "   %5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.5f",
		     etamin, etamax, nparnew, xmin, xmax,
		     ptmin0, ptmax1,
		     p2, p3, p4, p5, p6, p7,
		     p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])
	     << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d"
		   "   %d   %d"
		   "   %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"
		   "   %5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.5f",
		   etamin, etamax, nparnew, xmin, xmax,
		   ptmin0, ptmax1,
		   p2, p3, p4, p5, p6, p7,
		   p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])
	   << endl;
      ++cnt;
    }
    else {
      int nparnew = 15;
      if (cnt<cntmax && debug)
	cout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		     "   %8.6f %8.6f   "
		     "%5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f",
		     etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		     p2, p3, p4, p5,
		     p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   "
		   "%5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f",
		   etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,
		   p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
      ++cnt;
    }

  } // while getline

} // creataL2L3ResTextFile
