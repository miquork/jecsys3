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
#include "TGraphErrors.h"

#include "tdrstyle_mod22.C"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

const bool debug = true;

void createL2L3ResTextFiles(string set, bool leg2=false);

TCanvas *_c1(0), *_c3(0);
TLegend *_leg(0), *_leg2(0), *_leg3(0);
void createL2L3ResTextFile() {

  setTDRStyle();

  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    0.88+1e-4,1.21-1e-4,"p_{T,ref} (GeV)",ptmin,ptmax);
		    //0.89+1e-4,1.11-1e-4,"p_{T} (GeV)",ptmin,ptmax);
		    //0.87+1e-4,1.07-1e-4,"p_{T} (GeV)",ptmin,ptmax);
		    //0.75+1e-4,1.20-1e-4,"p_{T} (GeV)",ptmin,ptmax);
		    //0.90+1e-4,1.30-1e-4,"p_{T} (GeV)",ptmin,ptmax);
		    //0.88+1e-4,1.05-1e-4,"p_{T} (GeV)",ptmin,ptmax);
  //lumi_136TeV = "Run3, 63 fb^{-1}"; // Not including 23B
  //lumi_136TeV = "2024, 12.3 fb^{-1}"; // Not including 23B
  //lumi_136TeV = "2024, 27.0 fb^{-1}"; // June 6 hybrid
  //lumi_136TeV = "2024, 46.0 fb^{-1}"; // Aug 2 hybrid
  //lumi_136TeV = "2024, 92.1 fb^{-1}"; // Aug 2 hybrid
  //#include "Config.C"
  lumi_136TeV = Form("2024, %s","109 fb^{-1}");//mlum["2024"].c_str());
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  c1->SetLeftMargin(0.17);
  c1->SetRightMargin(0.03);
  h->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();
  _c1 = c1;

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  _leg = tdrLeg(0.45,0.90,0.75,0.90);
  _leg2 = tdrLeg(0.20,0.15,0.50,0.15);

  TH1D *h_3 = tdrHist("h_3","Absolute response at |#eta| < 1.3",
		      0.88+1e-4,1.21-1e-4,"p_{T,raw} (GeV)",5.,ptmax);
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);
  c3->SetLeftMargin(0.17);
  c3->SetRightMargin(0.03);
  h_3->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();
  _c3 = c3;

  l->SetLineStyle(kDashed);
  l->DrawLine(5.,1,ptmax,1);

  _leg3 = tdrLeg(0.45,0.90,0.75,0.90);
  

  
  //h->SetMinimum(0.70+1e-4); // Summer22
  /*
  createL2L3ResTextFiles("Run22CD-22Sep2023");
  createL2L3ResTextFiles("Run22E-22Sep2023");
  //createL2L3ResTextFiles("Run22FG-Prompt");
  createL2L3ResTextFiles("Run22F-Prompt");
  createL2L3ResTextFiles("Run22G-Prompt");
  createL2L3ResTextFiles("Run23C123-Prompt",true);
  createL2L3ResTextFiles("Run23C4-Prompt",true);
  createL2L3ResTextFiles("Run23D-Prompt",true);
  //createL2L3ResTextFiles("Run23C4D-Prompt",true);
  createL2L3ResTextFiles("Run3-Combo",true);
  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Summer22_reV3_Run22CDE-22Sep2023_22FG-Prompt_23CD-Prompt_Run3-Combo.pdf");
  */

  /*
  createL2L3ResTextFiles("Run23C123-Summer23",true);
  createL2L3ResTextFiles("Run23C4-Summer23",true);
  createL2L3ResTextFiles("Run23D-Summer23",true);
  c1->Update();
  //c1->SaveAs("pdf/createL2L3ResTextFile_Summer23_V1.pdf");
  c1->SaveAs("pdf/createL2L3ResTextFile_Summer23_V2.pdf");
  */

  /*
  //createL2L3ResTextFiles("Run23C123-Summer23",true);
  //createL2L3ResTextFiles("Run23D-Summer23",true);
  //createL2L3ResTextFiles("Run24BC-Prompt",true);
  createL2L3ResTextFiles("Run24BCD-Prompt",true);
  createL2L3ResTextFiles("Run24E-Prompt",true);
  //createL2L3ResTextFiles("Run24CR-ECALRATIO",true);
  //createL2L3ResTextFiles("Run24CS-HCALDI",true);
  createL2L3ResTextFiles("Run24F-Prompt",true);
  createL2L3ResTextFiles("Run24G-Prompt",true);
  createL2L3ResTextFiles("Run24H-Prompt",true);
  createL2L3ResTextFiles("Run24I-Prompt",true);
  gPad->RedrawAxis();
  c1->Update();
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V1M_Summer23_V2.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V2M.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V3M.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V4M.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V5M.pdf");
  //c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V6M.pdf");
  c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V7M.pdf");
  */

  createL2L3ResTextFiles("2024I_nib1",true);

  createL2L3ResTextFiles("2024H_nib1",true);
  createL2L3ResTextFiles("2024G_nib2",true);
  createL2L3ResTextFiles("2024G_nib1",true);
  createL2L3ResTextFiles("2024F_nib3",true);
  createL2L3ResTextFiles("2024F_nib2",true);
  createL2L3ResTextFiles("2024F_nib1",true);
  createL2L3ResTextFiles("2024Ev2_nib1",true);
  createL2L3ResTextFiles("2024Ev1_nib1",true);

  createL2L3ResTextFiles("2024D_nib1",true);
  
  createL2L3ResTextFiles("2024C_nib1",true);
  createL2L3ResTextFiles("2024B_nib1",true);

  c1->cd();
  gPad->RedrawAxis();
  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V8M_VsPtRef.pdf");
  
  c3->cd();
  gPad->RedrawAxis();
  c3->Update();
  c3->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V8M_VsPtRaw.pdf");

  
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

void createL2L3ResTextFiles(string set, bool leg2) {

  //if (debug) 
  cout << "****************************************************************\n";
  cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";
  cout << "****************************************************************\n";

  cout << "** Processing " << set << " **" << endl << flush;
  cout << "******************************\n";

  const char *cs = set.c_str();
  TString ts(cs);
  
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
  //else if (set=="Run23C123-Prompt") {
  else if (set=="Run23C123-Summer23") {
    f = new TFile("rootfiles/jecdataRun23C123.root","READ"); isRun3=true;
  }
  //else if (set=="Run23C4-Prompt") {
  else if (set=="Run23C4-Summer23") {
    f = new TFile("rootfiles/jecdataRun23C4.root","READ"); isRun3=true;
  }
  //else if (set=="Run23D-Prompt") {
  else if (set=="Run23D-Summer23") {
    f = new TFile("rootfiles/jecdataRun23D.root","READ"); isRun3=true;
  }
  //else if (set=="Run23C4D-Prompt") {
  else if (set=="Run23C4D-Summer23") {
    f = new TFile("rootfiles/jecdataRun23C4D.root","READ"); isRun3=true;
  }
  else if (set=="Run24BC-Prompt") {
    f = new TFile("rootfiles/jecdataRun24BC.root","READ"); isRun3=true;
  }
  else if (set=="Run24BCD-Prompt") {
    //f = new TFile("rootfiles/jecdataRun24BCD.root","READ"); isRun3=true;
    //f = new TFile("rootfiles/jecdataRun24BCDE.root","READ"); isRun3=true;
    f = new TFile("rootfiles/jecdataRun24BCD.root","READ"); isRun3=true;
  }
  else if (set=="Run24E-Prompt") {
    //f = new TFile("rootfiles/jecdataRun24BCDE.root","READ"); isRun3=true;
    f = new TFile("rootfiles/jecdataRun24E.root","READ"); isRun3=true;
  }
  else if (set=="Run24F-Prompt") {
    f = new TFile("rootfiles/jecdataRun24F.root","READ"); isRun3=true;
  }
  else if (set=="Run24G-Prompt") {
    f = new TFile("rootfiles/jecdataRun24G.root","READ"); isRun3=true;
  }
  else if (set=="Run24H-Prompt") {
    f = new TFile("rootfiles/jecdataRun24H.root","READ"); isRun3=true;
  }
  else if (set=="Run24I-Prompt") {
    f = new TFile("rootfiles/jecdataRun24I.root","READ"); isRun3=true;
  }
  else if (set=="Run24CR-ECALRATIO") {
    f = new TFile("rootfiles/jecdataRun24CR.root","READ"); isRun3=true;
  }
  else if (set=="Run24CS-HCALDI") {
    f = new TFile("rootfiles/jecdataRun24CS.root","READ"); isRun3=true;
  }
  else if (set=="Run3-Combo") {
    // use input to the global fit for pre-fit average JEC
    f = new TFile("rootfiles/jecdataRun3Data.root","READ"); isRun3=true;
  }
  else if (ts.Contains("nib")) {
    f = new TFile(Form("rootfiles/jecdata%s.root",cs),"READ"); isRun3=true;
  }
  else
    f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");

  assert(f && !f->IsZombie());
  TH1D *h(0);
  //h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit");
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2"); // Run2 refit
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/herr_l2l3res"); // Run2 ref. JES
  assert(h);
  curdir->cd();

  TF1 *f1(0), *f1raw(0);
  bool isECALCC = (ts.Contains("nib") &&
		   (ts.Contains("B") || ts.Contains("C") || ts.Contains("D") ||
		    ts.Contains("E")));
  bool isHBoff = (ts.Contains("nib") &&
		  (ts.Contains("F_nib2") || ts.Contains("F_nib3") ||
		   ts.Contains("G") || ts.Contains("H") || ts.Contains("I")));
  if (ts.Contains("nib")) {
    f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/(0.1*x)+[2]*log10(0.1*x)/(0.1*x)+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow(x/15.,2.)/(1+0.5*pow(x/15.,4.))",15,4500);
    f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001, 0, 0.10);
    if (!isECALCC) {
      f1->FixParameter(3,0.);
      f1->FixParameter(4,1500.);
      f1->FixParameter(5,0.);
      f1->FixParameter(7,0.);
    }
    if (isECALCC) {
      //f1->SetParameters(0.9182, 0.5000, 0.80155, -0.13124, 1680.79, 2.2711, -0.45400, 0.00002991); // BCD V5M
      //f1->SetParameters(0.938231, 0.047432, 0.063449, -0.0616438, 1565.16, 2.97389, -0.118774, 0.0126795, 0.0432517); // 2024D_nib1 chi2/NDF=4.6/55
      //f1->SetParameters(0.98,0.,0., -0.06,1600.,3., -0.12, 0.04,0.);
      f1->SetParameters(0.938231, 0.047432, 0.063449, -0.0616438, 1565.16, 2.97389, -0.118774, 0.0126795, 0.0432517); // 2024D_nib1 chi2/NDF=4.6/55
    }
    
    f1raw = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/(0.1*x)+[2]*log10(0.1*x)/(0.1*x)+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow((x+[9])/15.,2.)/(1+0.5*pow((x+[9])/15.,4.))",15,4500);
    //f1raw->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001, 0, 0.10, 0.);
    /*
    if (!isECALCC) {
      f1raw->FixParameter(3,0.);
      f1raw->FixParameter(4,1500.);
      f1raw->FixParameter(5,0.);
      f1raw->FixParameter(7,0.);
    }
    */
  }
  else {
    f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*x",15,4500);
    f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001, 0);
    if (set=="2017H") 
      f1->SetParameters(0.99, 1.5,0.01, 0.01,1000.,1.3, 0.001, 0.);
    if (!isRun3)
      f1->FixParameter(7, 0.);
  }
  /*
  // Avoid divergence at pT=15 GeV
  if (set=="Run24BCD-Prompt" || set=="Run24E-Prompt" ||
      ts.Contains("nib") && isECALCC) {
      //ts.Contains("B") || ts.Contains("C") || ts.Contains("D") ||
      //ts.Contains("E") || ts.Contains("nib")) {
    //f1->SetParameters(0.918151, 0.5, 0.801532, -0.131235, 1680.79, 2.27111, -0.453995, 2.99145e-05); // chi2/NDF=20.9/56
    //f1->SetParameters(1.02315, -0.0420837, 0.41949, 0.0577525, 11.1207, 0.0618524, -0.148952, -6.17298e-07); // Run24F-Prompt chi2/NDF=0.0/56
    //f1->SetParLimits(1,0.0,+0.5);

    f1->SetParameters(0.9182, 0.5000, 0.80155, -0.13124, 1680.79, 2.2711, -0.45400, 0.00002991); // BCD V5M
  }
  */
  /*
  // Extra adjustment for 2024H_nib1 V8M low pT, plus other FGHI
  if (ts.Contains("nib") && ts.Contains("H") ||
      ts.Contains("nib") && ts.Contains("G") ||
      ts.Contains("nib") && ts.Contains("I") ||
      ts.Contains("nib") && ts.Contains("F")) {
    f1->SetParameters(2.27559, -12.8966, 16.782, -0.658209, 315.177, 0.760886, -8.9014, -2.37033e-05); // 2024G_nib2 chi2/NDF=21.6/56
    //f1->FixParameter(3, 0.); // breaks?
    // => replace with off_nhf shape for better fit
  }

  if (ts.Contains("nib")) {
    if (ts.Contains("B") || ts.Contains("C") || ts.Contains("D") ||
	ts.Contains("E")) {
      f1->SetParameters(0.9182, 0.5000, 0.80155, -0.13124, 1680.79, 2.2711, -0.45400, 0.00002991, 0.05, 0.);
    }
    else if (ts.Contains("F") || ts.Contains("G") || ts.Contains("H") ||
	     ts.Contains("I")) {
      f1->SetParameters(0.9182, 0.5000, 0.80155, -0.13124, 1680.79, 2.2711, -0.45400, 0.00002991, 0.20, 0.);
    }   
  }
  */
  
  // To avoid very weird fits
  f1->SetParLimits(0,0.5,2.5);
  
  // To avoid division by zero errors
  f1->SetParLimits(4,10.,6500.);
  f1->SetParLimits(5,0.,10.);
  
  map<string,int> color;
  color["2016BCDEF"] = kYellow+2;
  color["2016GH"] = kRed;
  color["2017BCDEF"] = kGreen+2;
  color["2018ABCD"] = kBlue;

  color["Run22CD-22Sep2023"] = kGreen+2;
  color["Run22E-22Sep2023"] = kCyan+2;
  //color["Run22FG-Prompt"] = kRed;
  color["Run22F-Prompt"] = kRed;
  color["Run22G-Prompt"] = kRed+2;
  color["Run23C123-Prompt"] = kOrange+1;
  //color["Run23C4D-Prompt"] = kMagenta+2;
  color["Run23C4-Prompt"] = kBlue;
  color["Run23D-Prompt"] = kMagenta;
  color["Run24BC-Prompt"] = kGreen+2;
  color["Run24BCD-Prompt"] = kRed+2;
  color["Run24E-Prompt"] = kRed+1;
  color["Run24F-Prompt"] = kGreen+2;
  color["Run24G-Prompt"] = kCyan+2;
  color["Run24H-Prompt"] = kBlue+1;
  color["Run24I-Prompt"] = kMagenta+1;
  color["Run24CR-ECALRATIO"] = kGreen+2;
  color["Run24CS-HCALDI"] = kBlue+2;
  color["Run3-Combo"] = kYellow+2;

  color["Run23C123-Summer23"] = kOrange+1;
  color["Run23C4-Summer23"] = kBlue;
  color["Run23D-Summer23"] = kMagenta;
  
  color["2024B_nib1"] = kRed;
  color["2024C_nib1"] = kRed+2;
  color["2024D_nib1"] = kRed+3;
  color["2024Ev1_nib1"] = kOrange;
  color["2024Ev2_nib1"] = kOrange+1;
  color["2024F_nib1"] = kGreen+2;
  color["2024F_nib2"] = kBlue;
  color["2024F_nib3"] = kCyan+2;
  color["2024G_nib1"] = kCyan+4;
  color["2024G_nib2"] = kMagenta+1;
  color["2024H_nib1"] = kMagenta+2;
  color["2024I_nib1"] = kPink;
  
  h->Fit(f1,"QRN");
  h->Fit(f1,"QRNM");
  h->Fit(f1,"QRNM");

  _c1->cd();
  
  tdrDraw(h,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  h->SetFillColorAlpha(color[set]-9,0.7);
  f1->SetLineColor(color[set]);
  f1->Draw("SAME");

  if (ts.Contains("nib")) {
    //if (_leg->GetNRows()<6) {
    _leg->SetTextSize(0.03);
    _leg->SetY1(_leg->GetY1()-0.03);
    _leg->AddEntry(h,set.c_str(),"FL");
    //}
    //else {
    //_leg2->SetTextSize(0.03);
    //_leg2->SetY2(_leg2->GetY2()+0.03);
    //_leg2->AddEntry(h,set.c_str(),"FL");
    //}
  }
  else {
    if (leg2)
      _leg2->SetY2(_leg2->GetY2()+0.05);
    else
      _leg->SetY1(_leg->GetY1()-0.05);
    
    if (set=="Ref")
      _leg->AddEntry(h,"Old 22C (ref.)","FL");
    else if (leg2)
      _leg2->AddEntry(h,set.c_str(),"FL");
    else
      _leg->AddEntry(h,set.c_str(),"FL");
  }
  
  cout << "f1->SetParameters(";
  const int np = 9;//8;//7;
  double p[np];
  for (int i = 0; i != min(f1->GetNpar(),np); ++i) {
    p[i] = f1->GetParameter(i);
    cout << (i==0 ? "" : ", ") << p[i];
  } // for i
  if (f1raw) {
    //"[0]+[1]/(0.1*x)+[2]*log10(0.1*x)/(0.1*x)+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow((x+[9])/15.,2.)/(1+0.5*pow((x+[9])/15.,4.))",15,4500);
    double k_10 = f1->Eval(10.);
    double k_100 = f1->Eval(100.);
    double k_1500 = f1->Eval(1500.);
    f1raw->SetParameter(0, f1->GetParameter(0)); // 1
    f1raw->SetParameter(1, f1->GetParameter(1)/k_10); // 1/x
    f1raw->SetParameter(2, f1->GetParameter(2)/k_10); // log(x)/x
    f1raw->SetParameter(3, f1->GetParameter(3));
    f1raw->SetParameter(4, f1->GetParameter(4)/k_1500);
    f1raw->SetParameter(5, f1->GetParameter(5));
    f1raw->SetParameter(6, f1->GetParameter(6)*pow(k_100,-0.3051));
    f1raw->SetParameter(7, f1->GetParameter(7)*k_1500);
    f1raw->SetParameter(8, f1->GetParameter(8));
    f1raw->SetParameter(9, 0.);
    
    if (!isECALCC) {
      f1raw->FixParameter(3,0.);
      f1raw->FixParameter(4,1500.);
      f1raw->FixParameter(5,0.);
      f1raw->FixParameter(7,0.);
    }
  }
  cout << Form("); // %s chi2/NDF=%1.1f/%d", set.c_str(),
	       f1->GetChisquare(), f1->GetNDF()) << endl << flush;


  _c3->cd();
  
  // Refit vs pT,raw
  TGraphErrors *graw3 = new TGraphErrors();
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double pt = h->GetBinCenter(i);
    double jes = h->GetBinContent(i);
    double jese = h->GetBinError(i);
    graw3->SetPoint(i-1, pt*jes, jes);
    graw3->SetPointError(i-1, pt*jese, jese);
  } // for i

  tdrDraw(graw3,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  graw3->SetFillColorAlpha(color[set]-9,0.7);

  if (f1raw) {
    graw3->Fit(f1raw,"QRN");
    graw3->Fit(f1raw,"QRNM");
    graw3->Fit(f1raw,"QRNM");    
    f1raw->SetLineColor(color[set]);
    f1raw->Draw("SAME");
  }
  
  
  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  const char *run = set.c_str();

  string sin(""), sout(""), sout2("");
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
  //
  if (set=="2016BCDEF" || set=="2016GH" || set=="2017BCDEF" || set=="2018ABCD")
   sout = Form("textfiles/%s_DATA_L2L3Residual_AK4PFchs.txt",run);

  bool isL2Res(false); bool isNewL2Res(false);
  if (set=="Run22CD-22Sep2023") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22-22Sep2023_Run2022CD_V3_DATA_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22-22Sep2023_Run2022CD_reV3_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run22E-22Sep2023") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22EE-22Sep2023_Run2022E_V3_DATA_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22EE-22Sep2023_Run2022E_reV3_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run22F-Prompt") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22EEPrompt22_Run2022F_V3_DATA_L2Residual_AK4PFPuppi.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22EEPrompt22_Run2022F_reV3_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run22G-Prompt") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22EEPrompt22_Run2022G_V3_DATA_L2Residual_AK4PFPuppi.txt";  isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22EEPrompt22_Run2022G_reV3_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  /*
  if (set=="Run23C123-Prompt") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22Prompt23_Run2023Cv123_V3_DATA_L2Residual_AK4PFPUPPI.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22Prompt23_Run2023Cv123_reV3_DATA_L2L3Residual_AK4PFPUPPI.txt";
  }
  if (set=="Run23C4-Prompt") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22Prompt23_Run2023Cv4_V3_DATA_L2Residual_AK4PFPUPPI.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22Prompt23_Run2023Cv4_reV3_DATA_L2L3Residual_AK4PFPUPPI.txt";
  }
  if (set=="Run23D-Prompt") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22Prompt23_Run2023D_V3_DATA_L2Residual_AK4PFPUPPI.txt"; isL2Res = true;
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22Prompt23_Run2023D_reV3_DATA_L2L3Residual_AK4PFPUPPI.txt";
  }
  */
  if (set=="Run23C123-Summer23") {
    sin = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023Cv123_V1_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023Cv123_V2_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run23C4-Summer23") {
    sin = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023Cv4_V1_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023Cv4_V2_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run23D-Summer23") {
    sin = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023D_V1_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Summer23_L2ResOnly/Summer23Prompt23_Run2023D_V2_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24BC-Prompt") {
    //sin = "textFiles/Prompt24/Prompt24_Run2024BC_V1M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024BC_V1M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024BC_V2M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024BC_V2M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24BCD-Prompt") {
    //sin = "textFiles/Prompt24/Prompt24_Run2024BCD_V3M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024BCD_V3M_DATA_L2L3Residual_AK4PFPuppi.txt";
    //sin = "textFiles/Prompt24/Prompt24_Run2024BCD_V4M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sin = "textFiles/Prompt24_V5M/Prompt24_Run2024BCD_V5M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024BCD_V5M_DATA_L2L3Residual_AK4PFPuppi.txt";
    //sin = "textFiles/Prompt24_V6M/Prompt24_Run2024BCD_V6M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024BCD_V6M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024BCD_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024BCD_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24E-Prompt") {
    //sin = "textFiles/Prompt24/Prompt24_Run2024E_V4M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sin = "textFiles/Prompt24_V5M/Prompt24_Run2024E_V5M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024E_V5M_DATA_L2L3Residual_AK4PFPuppi.txt";
    //sin = "textFiles/Prompt24_V6M/Prompt24_Run2024E_V6M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024E_V6M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024E_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024E_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24F-Prompt") {
    //sin = "textFiles/Prompt24_V5M/Prompt24_Run2024F_V5M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024F_V5M_DATA_L2L3Residual_AK4PFPuppi.txt";
    //sin = "textFiles/Prompt24_V6M/Prompt24_Run2024F_V6M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024F_V6M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024F_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024F_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24G-Prompt") {
    //sin = "textFiles/Prompt24_V6M/Prompt24_Run2024G_V6M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024G_V6M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024G_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024G_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24H-Prompt") {
    sin = "textFiles/Prompt24/Prompt24_Run2024H_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024H_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24I-Prompt") {
    sin = "textFiles/Prompt24/Prompt24_Run2024I_V7M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024I_V7M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24CR-ECALRATIO") {
    //sin = "textFiles/Prompt24/Prompt24_Run2024CR_V3M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    //sout = "textFiles/Prompt24/Prompt24_Run2024CR_V3M_DATA_L2L3Residual_AK4PFPuppi.txt";
    sin = "textFiles/Prompt24/Prompt24_Run2024CR_V4M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024CR_V4M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (set=="Run24CS-HCALDI") {
    sin = "textFiles/Prompt24/Prompt24_Run2024CS_V4M_DATA_L2Residual_AK4PFPuppi.txt"; isNewL2Res = true;
    sout = "textFiles/Prompt24/Prompt24_Run2024CS_V4M_DATA_L2L3Residual_AK4PFPuppi.txt";
  }
  if (ts.Contains("nib")) {
    sin = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2ResidualVsPtRef_AK4PFPuppi.txt",cs); isNewL2Res = true;
    sout = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2L3ResidualVsPtRef_AK4PFPuppi.txt",cs);
    sout2 = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2L3Residual_AK4PFPuppi.txt",cs);
  }
  
  if (set=="Run3-Combo") {
    sin = "textFiles/Run3_22Sep2023_v3/Summer22Prompt23_Run2023D_V3_DATA_L2Residual_AK4PFPUPPI.txt"; isL2Res = true; // placeholder L2Res
    sout = "textFiles/Run3_22Sep2023_reV3/Summer22Prompt23_Run3_reV3_DATA_L2L3Residual_AK4PFPUPPI.txt";
  }

  assert(sin!="");
  assert(sout!="");
  assert(sout!=sin);
						       
  cout << "Reading in L2Residual file:" << endl
       << "   " << sin << endl << flush;  
  ifstream fin(sin.c_str());
  assert(fin.is_open());
  
  cout << "Writing out L2(L3)Residual file:" << endl
       << "   " << sout << endl << flush;
  ofstream fout(sout.c_str());

  // 2018: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018

  string header;
  getline(fin, header);
  if (debug) cout << "Old L2L3Residual header:" << endl;
  if (debug) cout << header << endl;

  //header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}";
  if (isRun3 && isNewL2Res && ts.Contains("nib")) {
    header =
      "{ 1 JetEta 1 JetPt 1./([0]+[1]*log10(0.01*x)+[2]/(0.1*x)+"
      "[3]*log10(0.1*x)/(0.1*x)+[4]*(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1)+"
      "[7]*pow(x,-0.3051)+[8]*x+"
      "[9]*pow(x/15.,2.)/(1+0.5*pow(x/15.,4.)))"
      " Correction L2Relative}";
  }
  else if (isRun3 && isNewL2Res) {
    header =
      "{ 1 JetEta 1 JetPt 1./([0]+[1]*log10(0.01*x)+[2]/(x/10.))"
      "*1./([3]+[4]/x+[5]*log(x)/x+[6]*(pow(x/[7],[8])-1)/(pow(x/[7],[8])+1)+[9]*pow(x,-0.3051)+[10]*x)"
      " Correction L2Relative}";
  }
  else if (isRun3) {
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
  double etamin(0), etamax(0);
  int npar(0), xmin(0), xmax(0), ptmin0(0), ptmax1(0);
  double p0(0), p1(0); // isNewL2Res
  double p2(0), p3(0), p4(0), p5(0);
  double p6(0), p7(0), p8(0);
  double p9(0), p10(0); // isRun3
  double p11(0); // nib
  int cnt(0); int cntmax(0);
  while (getline(fin,line)) {
    if (cnt<cntmax && debug) cout << line << endl;
    if (isRun3 && isNewL2Res) {
      assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %lf %lf %lf",
		    &etamin, &etamax, &npar, &xmin, &xmax,
		    &p0, &p1, &p2)==8);
    }
    else if (isRun3 && isL2Res) {
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

    if (isRun3 && isNewL2Res && ts.Contains("nib")) {
      int nparnew = 2 + 3 + 9; // 14
      string s = Form("  %6.3f %6.3f %2d  %2d %4d   "
		      "%7.4f %7.4f %7.4f   "
		      "%5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.8f  "
		      "%5.4f",
		      etamin, etamax, nparnew, xmin, xmax,
		      p0, p1, p2,
		      p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
      if (cnt<cntmax && debug)
	cout << s << endl;
      fout << s << endl;
      ++cnt;
    }
    else if (isRun3 && isNewL2Res) {

      int nparnew = 2 + 3 + 8; // 13
      string s = Form("  %6.3f %6.3f %2d  %2d %4d   "
		      "%7.4f %7.4f %7.4f   "
		      "%5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.8f",
		      etamin, etamax, nparnew, xmin, xmax,
		      p0, p1, p2,
		      p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
      if (cnt<cntmax && debug)
	cout << s << endl;
      fout << s << endl;
      ++cnt;
    }
    else if (isRun3) {
      if (isL2Res) {
	// patch L2L3Res files
	if (etamin > etamax) {
	  assert(etamin<=0);
	  assert(etamax<0);
	  cout << "Swap sign of " << etamin << " and " << etamax << endl;
	  etamin *= -1;
	  etamax *= -1;
	}
      }
      int nparnew = 18;//17;
      if (cnt<cntmax && debug)
	cout << Form("  %9.6f %9.6f   %d   %d %d"
		     "   %d   %d"
		     "   %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"
		     "   %5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.8f",
		     etamin, etamax, nparnew, xmin, xmax,
		     ptmin0, ptmax1,
		     p2, p3, p4, p5, p6, p7,
		     p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])
	     << endl;
      fout << Form("  %9.6f %9.6f   %d   %d %d"
		   "   %d   %d"
		   "   %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f"
		   "   %5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.8f",
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


  // New code to merge L2Res and L3Res vs pTref, then re-map combination
  // to pT,raw. This needs to be done separately for each |eta| bin
  if (sout2!="") {
    assert(isRun3);
    assert(isNewL2Res);

    cout << "Running new code for merging L2Res+L3Res vs pTref," << endl
	 << "then remapping to vs pT,raw for each eta bin" << endl;
    cout << "Reading in L2Residual file:" << endl
	 << "   " << sin << endl << flush;  
    ifstream fin(sin.c_str());
    assert(fin.is_open());

    // Input L2Residual header
    string header;
    getline(fin, header);
    if (debug) cout << "Input L2Residual header:" << endl;
    if (debug) cout << header << endl;
    const int nparold = 2 + 3 + 3; // 8
    if (debug) cout << "Input L3Residual function:" << endl;
    //if (debug) cout << "1./("<<f1->GetExpFormula().Data()<<")" << endl;
    if (debug) cout << "1./("<<f1raw->GetExpFormula().Data()<<")" << endl;

    // New combined function and header
    string func = "[0]+[1]*log10(0.01*x)+[2]/(0.1*x)+[3]*log10(0.1*x)/(0.1*x)+[4]*(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1)+[7]*pow(x,-0.3051)+[8]*(0.001*x)+[9]*pow((x+[10])/15.,2.)/(1+0.5*pow((x+[10])/15.,4.))";
    /*
    if (ts.Contains("B") || ts.Contains("C") || ts.Contains("D") ||
	ts.Contains("E")) func +=  "+[8]*(x/1000.)";
    else
      //func += "+[8]*pow(x/15.6,1.927)/(1+0.5*pow(x/15.6,3.853))";
      //func += Form("+[8]*pow((x%+1.2f)/16.,2.)/(1+0.5*pow((x%+1.2f)/16.,4.))",
      //		   f1raw->GetParameter(9),f1raw->GetParameter(9));
      func += "+[8]*pow((x+[9])/16.,2.)/(1+0.5*pow((x+[9])/16.,4.))";
    */
    header = "{1 JetEta 1 JetPt 1./("+func+") Correction L2Relative";
    
    cout << "Writing out L2L3Residual file:" << endl
	 << "   " << sout2 << endl << flush;
    if (debug) cout << "Output L2L3Residual header:" << endl;
    if (debug) cout << header << endl;
    ofstream fout(sout2.c_str());
    const int nparnew = 2 + 3 + 11;//8;//9; // 14
    fout << header << endl;
    
    string line;
    double etamin(0), etamax(0);
    int npar(0), xmin(0), xmax(0), ptmin0(0), ptmax1(0);
    double p0(0), p1(0), p2(0);
    double p3(0), p4(0), p5(0), p6(0), p7(0), p8(0), p9(0), p10(0), p11(0);
    int cnt(0), ieta(0), cntmax(0);

    const int nxy = 41;
    TCanvas *cx = new TCanvas(Form("cx_%s",cs),"cx",9*300,5*300);
    cx->Divide(7,6,0,0);
      
    while (getline(fin,line)) {

      // Read in L2Res
      if (cnt<cntmax && debug) cout << line << endl;
      assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %lf %lf %lf",
		    &etamin, &etamax, &npar, &xmin, &xmax,
		    &p0, &p1, &p2)==nparold);
      assert(npar==2+3);
      
      // TF1 for relative JES (1./L2Residual)
      TF1 *f2 = new TF1(Form("f2_%s_%d",cs,cnt),
			"[0]+[1]*log10(0.01*x)+[2]/(x/10.)",xmin,xmax);
      f2->SetParameters(p0, p1, p2);

      // Generate graph with 100 points logarithmically distributed between
      // pTref=15 GeV / 1.40 (JER 40%*1) and E=sqrt(s)/2. * 1.2 (JER 10%*2)
      // (limited to pTmax<4500 GeV)
      const int npt = 100;
      const double ptrefmin = 15;
      const double ptrefmin2 = ptrefmin / 1.4;
      const double erefmax = (13600./2.);
      const double erefmax2 = erefmax * 1.2;
      double eta = min(fabs(etamin),fabs(etamax));
      double ptrefmax = erefmax / cosh(eta);
      double ptrefmax2 = min(4500., erefmax2 / cosh(eta));

      TGraph *gref = new TGraph();
      TGraph *graw = new TGraph();
      //TGraph *graw = new TGraphErrors();
      TGraph *gr = new TGraph();
      double c = pow(ptrefmax2/ptrefmin2,1./(npt-1));
      for (int i = 0; i != npt; ++i) {
	double ptref = ptrefmin2 * pow(c, i);
	//double jes = f1->Eval(ptref) * f2->Eval(ptref);
	//double jesref = h->Interpolate(ptref);
	double jesref = f1->Eval(ptref);
	double jes = jesref * f2->Eval(ptref);
	double ptraw = ptref*jes;
	gref->SetPoint(i, ptref, jes);
	graw->SetPoint(i, ptraw, jes);
	//graw->SetPoint(i, 0., jes);
	// Bias in current corrections when f1,f2 evaluated at wrong pt
	//double jes2 = f1->Eval(ptraw) * f2->Eval(ptraw);
	//double jesraw = h->Interpolate(ptraw);
	double jesraw = f1->Eval(ptraw);
	double jes2 = jesraw * f2->Eval(ptraw);
	gr->SetPoint(i, ptref, jes2/jes);
      }

      // Reference JES refit vs <pT,raw>
      // "[0]+[1]*log10(x/100.)+[2]/(x/10.)+[3]*log10(x/10.)/(x/10.)+"
      // "[4]*(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1)+[7]*pow(x,-0.3051)+"
      // "[8]*(x/1000.)+[9]*pow((x+[10])/15.,2)/(1+0.5*pow((x+[10])/15.,4))"
      // f1raw:
      // "[0]+[1]/(x/10.)+[2]*log10(x/10.)/(x/10.)+"
      // "[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)"
      // "+[7]*(x/1000.)+[8]*pow((x-[9])/15.,2)/(1+0.5*pow((x+[9])/15.,4))"
      // f2 (p1,p2,p3):
      // "[0]+[1]*log10(0.01*x)+[2]/(x/10.)"
      //double ptminraw = ptrefmin2*f1->Eval(ptrefmin2)*f2->Eval(ptrefmin2);
      //double ptmaxraw = ptrefmax2*f1->Eval(ptrefmax2)*f2->Eval(ptrefmax2);
      //double jesrefmin2  = h->Interpolate(ptrefmin2);
      //double jesrefmax2 = h->Interpolate(ptrefmax2);
      double jesrefmin2  = f1->Eval(ptrefmin2);
      double jesrefmax2 = f1->Eval(ptrefmax2);
      double ptminraw = ptrefmin2*jesrefmin2*f2->Eval(ptrefmin2);
      double ptmaxraw = ptrefmax2*jesrefmax2*f2->Eval(ptrefmax2);
      const int nparnew = 11;//9;
      //double ptrefmid2 = sqrt(ptrefmin2*ptrefmax2);
      double k2_1000 = f2->Eval(1000.);
      double k1_1000 = f1->Eval(1000.);
      double kjes_1000 = k2_1000*k1_1000;
      //double k1_1000 = h->Interpolate(1000.);
      double k2_100 = f2->Eval(100.);
      double k1_100 = f1->Eval(100.);
      double kjes_100 = k2_100*k1_100;
      //double k1_100 = h->Interpolate(100.);
      double k2_10 = f2->Eval(10.);
      double k1_10 = f1->Eval(10.);
      double kjes_10 = k2_10*k1_10;
      //double k1_10 = h->Interpolate(10.);
      TF1 *f23 = new TF1(Form("f23_%s_%d",cs,cnt),func.c_str(),
			 floor(ptminraw), ceil(ptmaxraw));
      f23->SetParameter(0, k2_100*f1raw->GetParameter(0)); // 1
      //f23->SetParameter(0, k1_100*p0); // 1
      //f23->SetParameter(1, k2_100*p1); // log10(x/100.)
      f23->SetParameter(1, k1_100*p1); // log10(x/100.)
      //f23->SetParameter(2, k2_10*0.1*f1raw->GetParameter(1)*k1_10 +
      f23->SetParameter(2, k2_10*f1->GetParameter(1)*k1_10 +
			k1_10*p2*k2_10); // 1/(x/10)
      //f23->SetParameter(2, k2_10*0.1*f1raw->GetParameter(1) +
      //		k1_10*p2); // 1/(x/10)
      f23->SetParameter(3, 0.); // log10(x/10)/(x/10)

      /*
      f23->SetParameter(4, k2_100*f1raw->GetParameter(3));
      //f23->SetParameter(5, k2_100*k1_100*f1raw->GetParameter(4));
      f23->SetParameter(5, k2_100*f1raw->GetParameter(4));
      f23->SetParameter(6, f1raw->GetParameter(5));
      f23->SetParameter(7, k2_100*f1raw->GetParameter(6));
      */
      f23->SetParameter(4, k2_100*f1->GetParameter(3));
      f23->SetParameter(5, f1->GetParameter(4)/kjes_100);
      f23->SetParameter(6, f1->GetParameter(5));
      f23->SetParameter(7, k2_100*f1raw->GetParameter(6)*pow(kjes_100,-0.3051));
      
      //if (ts.Contains("B") || ts.Contains("C") || ts.Contains("D") ||
      //  ts.Contains("E")) {
      //if (isECALCC) {
      //f23->SetParameter(8, k2_1000*1000.*f1raw->GetParameter(7)/k1_1000);
      //f23->SetParameter(8, k2_1000*f1raw->GetParameter(7));///k1_1000);
      f23->SetParameter(8, k2_1000*f1raw->GetParameter(7)*kjes_100);
      //f23->SetParameter(8, k2_1000*1000.*f1raw->GetParameter(7));
	//else {
      //f23->SetParameter(9, k2_10*f1raw->GetParameter(8));
      //f23->SetParameter(10, f1raw->GetParameter(9));
      f23->SetParameter(9, k2_10*f1->GetParameter(8)*pow(kjes_10,2));
      f23->SetParameter(10, 0.);//f1->GetParameter(9));
      
      //}
	
      // To avoid very weird fits
      f23->SetParameter(0, max(0.3,min(0.25,f23->GetParameter(0))));
      f23->SetParLimits(0,0.3,2.5);
	
      // To avoid division by zero errors
      f23->SetParameter(4, max(10.,min(6500.,f23->GetParameter(4))));
      f23->SetParLimits(4,10.,6500.);
      f23->SetParameter(5, max(0.,min(10.,f23->GetParameter(5))));
      f23->SetParLimits(5,0.,10.);
      //if (f23->GetNpar()>8)
      f23->SetParLimits(10,-10,+10.);

      if (!isECALCC) {
	f23->FixParameter(4,0.);
	f23->FixParameter(5,1500.);
	f23->FixParameter(6,0.);
	f23->FixParameter(8,0.);
      }	
      
      
      graw->Fit(f23,"QRN");
      graw->Fit(f23,"QRNM");
      graw->Fit(f23,"QRNM");
      
      // Output L2L3Res
      double p[nparnew];
      for (int ip = 0; ip != nparnew; ++ip) {
	p[ip] = f23->GetParameter(ip);
      }

      //string s = Form("  %6.3f %6.3f %2d  %2d %4d   "
      //	      "%5.4f %5.4f %5.5f %5.5f %5.2f %5.4f %5.5f %5.8f %5.8f",
      string s = Form("  %6.3f %6.3f %2d   %2d %4d   "
		      "%6.4f %+7.4f   %+7.4f %+7.4f   "
		      "%7.4f %5.1f %6.4f   %+7.3f   %+7.4f   "
		      "%+7.3f   %+7.4f",
		      etamin, etamax, nparnew, int(ptminraw), int(ptmaxraw),
		      p[0],p[1], p[2],p[3], p[4],p[5],p[6],  p[7], p[8],
		      p[9],p[10]);
      if (cnt<cntmax && debug)
	cout << s << endl;
      fout << s << endl;
      ++cnt;

      // Plotting only for positive side (is anyway symmetric)
      if (etamin<0) continue;
      ++ieta;
      
      // Calculate ratio of fit and graph to check quality
      TGraph *gr2 = new TGraph(graw->GetN());
      for (int i = 0; i != npt; ++i) {
	double ptref = gref->GetX()[i];
	double ptraw = graw->GetX()[i];
	double jes1 = graw->GetY()[i];
	double jes2 = f23->Eval(ptraw);
	//gr2->SetPoint(i, ptref, jes2/jes1);
	gr2->SetPoint(i, ptref, 1+(jes2/jes1-1)*10.);
      }

      cx->cd(ieta);
      double eps = 1e-4;
      TH1D *hx = tdrHist(Form("hx_%s_%d",cs,ieta),"Rel. JES Data/MC refit",
			 //0.65+eps,1.35-eps,"p_{T} (GeV)",5.,3500.);
			 0.65+eps,1.35-eps,"p_{T} (GeV)",5.,4500.);
      // 8x6
      if      (eta<1.218) hx->GetYaxis()->SetRangeUser(0.85+eps,1.25-eps);
      else if (eta<1.830) hx->GetYaxis()->SetRangeUser(0.55+eps,1.20-eps);
      else if (eta<2.853) hx->GetYaxis()->SetRangeUser(0.30+eps,1.20-eps);
      else if (eta<4.013) hx->GetYaxis()->SetRangeUser(0.70+eps,1.25-eps);
      else if (eta<5.191) hx->GetYaxis()->SetRangeUser(0.60+eps,1.35-eps);
      // Adjust HF for 2024BCD (and earlier)
      if (ts.Contains("B") || ts.Contains("C") || ts.Contains("D")) {
	if (eta>4.013 && eta<5.191)
	  hx->GetYaxis()->SetRangeUser(0.30+eps,1.25-eps);
      }
      if (ts.Contains("E") || ts.Contains("F_nib1")) {
	if (eta>2.853 && eta<4.013)
	  hx->GetYaxis()->SetRangeUser(0.30+eps,1.30-eps);
      }
      
      hx->Draw();
      gPad->SetLogx();

      TLine *l = new TLine();
      l->SetLineColor(kGray+1);
      l->SetLineStyle(kDotted);
      l->DrawLine(ptrefmax,0.30,ptrefmax,1.35);
      l->DrawLine(ptrefmin,0.30,ptrefmin,1.35);
      l->SetLineStyle(kDashed);
      l->DrawLine(5.,1,3500.,1);

      f1->SetRange(ptrefmin2,ptrefmax2);
      f1->SetLineColor(kBlue-9);
      f1->DrawClone("SAME");
      TH1D *hc = (TH1D*)h->Clone(Form("%s_clone",h->GetName()));
      hc->GetXaxis()->SetRangeUser(ptrefmin2,ptrefmax2);
      tdrDraw(hc,"L",kNone,kBlue-9,kSolid,-1,kNone);

      f2->SetRange(ptrefmin2,ptrefmax2);
      f2->SetLineColor(kCyan+1);
      f2->Draw("SAME");

      tdrDraw(gref, "L", kNone, kBlue, kSolid, -1, kNone, 0, 0.5);
      tdrDraw(graw, "Pz", kOpenCircle, kBlack, kSolid, -1, kNone, 0, 0.5);
      tdrDraw(gr, "L", kNone, kGreen+2, kSolid, -1, kNone, 0, 0.5);
      tdrDraw(gr2, "L", kNone, kMagenta+1, kSolid, -1, kNone, 0, 0.5);

      f23->SetLineColor(kRed);
      f23->Draw("SAME");

      TLatex *tex = new TLatex();
      tex->SetNDC(); tex->SetTextSize(0.045);
      tex->DrawLatex(0.50,0.85,Form("[%1.3f,%1.3f]",etamin,etamax));
      if (ieta==nxy) {
	cx->cd(ieta+1);
	//TLegend *legx = tdrLeg(0.05,0.90-0.06*2.5*5,0.30,0.90);
	TLegend *legx = tdrLeg(0.05,0.90-0.06*2.0*7,0.30,0.90);
	legx->SetTextSize(0.045*2.0);//2.5);
	legx->AddEntry(f1,"JES-L3Res fit vs p_{T,ref}","L");
	legx->AddEntry(f2,"JES-L2Res fit vs p_{T,ref}","L");
	legx->AddEntry(gref,"JES vs p_{T,ref}","P");
	legx->AddEntry(graw,"JES vs p_{T,raw}","P");
	legx->AddEntry(f23,"Fit-L2L3Res vs p_{T,raw}","L");
	legx->AddEntry(gr,"Ratio of JES vs p_{T,ref}","L");
	legx->AddEntry(gr2,"Ratio of fit/JES vs p_{T,raw} #times 10","L");

	cx->cd(ieta);

        #include "Config.C"
	/*
	std::map<std::string, std::string> mlum;
	// NIB-level:
	mlum["2022B_nib1"] = "0.097 fb^{-1}";
	mlum["2022C_nib1"] = "5.0 fb^{-1}";
	mlum["2022D_nib1"] = "2.97 fb^{-1}";
	mlum["2022E_nib1"] = "5.8 fb^{-1}";
	mlum["2022F_nib1"] = "17.8 fb^{-1}";
	mlum["2022G_nib1"] = "3.08 fb^{-1}";
	mlum["2023B_nib1"] = "0.64 fb^{-1}";
	mlum["2023C_nib1"] = "7.2 fb^{-1}";
	mlum["2023Cv4_nib1"] = "10.4 fb^{-1}";
	mlum["2023Cv4_nib2"] = "0.407 fb^{-1}";
	mlum["2023D_nib1"] = "9.7 fb^{-1}";
	mlum["2024B_nib1"] = "0.130 fb^{-1}";
	mlum["2024C_nib1"] = "7.2 fb^{-1}";
	mlum["2024D_nib1"] = "8.0 fb^{-1}";
	mlum["2024Ev1_nib1"] = "6.3 fb^{-1}";
	mlum["2024Ev2_nib1"] = "5.0 fb^{-1}";
	mlum["2024F_nib1"] = "0.88 fb^{-1}";
	mlum["2024F_nib2"] = "14.4 fb^{-1}";
	mlum["2024F_nib3"] = "12.5 fb^{-1}";
	mlum["2024G_nib1"] = "16.9 fb^{-1}";
	mlum["2024G_nib2"] = "20.9 fb^{-1}";
	mlum["2024H_nib1"] = "5.4 fb^{-1}";
	mlum["2024I_nib1"] = "11.5 fb^{-1}";
	*/
	double siz = tex->GetTextSize();
	tex->SetTextSize(siz*2.0);//2.5);
	tex->DrawLatex(0.50,0.70,Form("%s vs",cs));
	tex->DrawLatex(0.50,0.55,"Winter 24");
	tex->DrawLatex(0.50,0.40,mlum[cs].c_str());
	tex->SetTextSize(siz);
      }
    } // for ieta

    cx->SaveAs(Form("pdf/createL2L3ResTextFile/createL2L3ResTextFile_PtRawReFit_%s.pdf",cs));
  } // sout2!=""
  if (_c1) _c1->cd();
  f1->SetLineColor(color[set]);
  f1->SetRange(15,4500.);
} // createL2L3ResTextFile
