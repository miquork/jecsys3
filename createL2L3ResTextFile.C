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

void createL2L3ResTextFiles(string set);

TCanvas *_c1(0), *_c3(0);
TLegend *_leg(0), *_leg1(0), *_leg2(0), *_leg3(0);
void createL2L3ResTextFile() {

  //if (debug)
  cout << endl;
  cout << "****************************************************************\n";
  cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";
  cout << "****************************************************************\n";
  cout << endl;
    
  setTDRStyle();

  // Plots vs pT,ref
  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    0.88+1e-4,1.21-1e-4,"p_{T,ref} (GeV)",ptmin,ptmax);
  lumi_136TeV = Form("2024, %s","109 fb^{-1}");
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
  _leg1 = tdrLeg(0.45,0.90,0.75,0.90);
  _leg2 = tdrLeg(0.20,0.15,0.50,0.15);

  // Plots vs pT,raw
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
  

  // Listing of individual IOVs
  createL2L3ResTextFiles("2024I_nib1");

  createL2L3ResTextFiles("2024H_nib1");
  createL2L3ResTextFiles("2024G_nib2");
  createL2L3ResTextFiles("2024G_nib1");
  createL2L3ResTextFiles("2024F_nib3");
  createL2L3ResTextFiles("2024F_nib2");
  createL2L3ResTextFiles("2024F_nib1");
  createL2L3ResTextFiles("2024Ev2_nib1");
  createL2L3ResTextFiles("2024Ev1_nib1");

  createL2L3ResTextFiles("2024D_nib1");

  createL2L3ResTextFiles("2024C_nib1");
  createL2L3ResTextFiles("2024B_nib1");

  c1->cd();
  gPad->RedrawAxis();
  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V8M_VsPtRef.pdf");
  
  c3->cd();
  gPad->RedrawAxis();
  c3->Update();
  c3->SaveAs("pdf/createL2L3ResTextFile_Prompt24_V8M_VsPtRaw.pdf");

} // createL2L3ResTextFile


void createL2L3ResTextFiles(string set) {

  cout << endl;
  cout << "******************************\n";
  cout << "** Processing " << set << " **" << endl << flush;
  cout << "******************************\n";

  const char *cs = set.c_str();
  TString ts(cs);
  
  // Simplify complex sum into an effective formula
  // Need good starting values and/or a few iterations to converge
  // Fit done to hjesfit from each IOV
  TDirectory *curdir = gDirectory;
  TFile *f(0);
  f = new TFile(Form("rootfiles/jecdata%s.root",cs),"READ");

  assert(f && !f->IsZombie());
  TH1D *h(0);
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2"); // Run2 refit
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/herr_l2l3res"); // Run2 ref. JES
  assert(h);
  curdir->cd();

  bool isECALCC = (ts.Contains("24B") || ts.Contains("24C") ||
		   ts.Contains("24D") || ts.Contains("24E"));
  bool isHBoff = (ts.Contains("24F_nib2") || ts.Contains("24F_nib3") ||
		  ts.Contains("24G") || ts.Contains("24H") ||
		  ts.Contains("24I"));
  
  TF1 *f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/(0.1*x)+[2]*log10(x)/(0.1*x)+[3]*(1+(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1))+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow(x*[9]/15.,2.)/(1+0.5*pow(x*[9]/15.,4.))",15,4500);
  TF1 *f1raw = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/(0.1*x)+[2]*log10(x)/(0.1*x)+[3]*(1+(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1))+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow((x*[9])/15.,2.)/(1+0.5*pow((x*[9])/15.,4.))",15,4500);

  string n1[10] =
  {"[0]","[1]/x","[2]*l/x","[3]*x^5","x/[4]^5",
   "x^[5]","[6]/x^.3","[7]*x","[8]*x^2","(x+[9])^2"};

  // Set initial guesses and limit or switch off unnecessary parts of fit
  if (!isECALCC) {
    f1->SetParameters(0.96, 0.,0., 0,1600,0, -0.12, 0., 0.16, 1.);

    // Remove "erf" component for "ecalcc"
    f1->FixParameter(3, 0.);
    f1->FixParameter(4, 1500.);
    f1->FixParameter(5, 0.);
    f1->FixParameter(7, 0.);
  }
  else {
    f1->SetParameters(0.94, 0.,0., -0.06,1600,3, -0.12, 0.01, 0.04, 1.);

    // Limit "erf" component for "ecalcc"
    //f1->SetParLimits(3, -0.3, 0.3);
    // To avoid division by zero errors
    f1->SetParLimits(3, -0.2, 0.2);
    f1->SetParLimits(4, 1000., 2000.);
    f1->SetParLimits(5, 1., 5.);
    f1->SetParLimits(7, -0.2, 0.2);
  }

  // To avoid very weird fits
  f1->SetParLimits(0, 0.5, 2.0);   // 1
  f1->SetParLimits(1, -0.5, +0.5); // 1/x
  f1->SetParLimits(2, -0.5, +0.5); // log(x)/x
  f1->SetParLimits(6, -0.5, +0.5); // x^-0.3
  f1->SetParLimits(8, 0., 0.5); // nhf_off
  f1->SetParLimits(9, 0.5, 2.0); // nhf_off peak
  
  // Other reasonable limitations
  //f1->SetParLimits(0,  0.5, 1.5);
  //f1->SetParLimits(1, -0.3, 0.3);
  //f1->SetParLimits(2, -0.3, 0.3);
  //f1->SetParLimits(6, -0.3, 0.3);
  //f1->SetParLimits(8, -0.3, 0.3);
    
  map<string,int> color;
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

  _leg1->SetTextSize(0.03);
  _leg1->SetY1NDC(_leg1->GetY1NDC()-0.03);
  _leg1->AddEntry(h,set.c_str(),"FL");

  cout << Form("Writing out textfiles/createL2L3ResTextFile_%s.txt",
	       set.c_str()) << endl;
  ofstream ftxt(Form("textfiles/createL2L3ResTextFile_%s.txt",set.c_str()));
  
  ftxt << f1->GetExpFormula() << endl;
  ftxt << "   f1->SetParNames  (";
  const int np = 11;//10;
  double p1ref[np];
  for (int i = 0; i != min(f1->GetNpar(),np); ++i) {
    if (i==1) ftxt << ",        ";
    ftxt << (i==0 ? "" : ", ") << Form("%8s",n1[i].c_str());
  } // for i
  ftxt << endl;
  
  cout << "f1->SetParameters(";
  ftxt << "   f1->SetParameters(";
  for (int i = 0; i != min(f1->GetNpar(),np); ++i) {
    p1ref[i] = f1->GetParameter(i);
    cout << (i==0 ? "" : ", ") << p1ref[i];
    if (i==1) ftxt << ",        ";
    ftxt << (i==0 ? "" : ", ") << Form("%+8.3f",p1ref[i]);
  } // for i
  cout << Form("); // %s chi2/NDF=%1.1f/%d", set.c_str(),
	       f1->GetChisquare(), f1->GetNDF()) << endl << flush;
  ftxt << Form("); // %s chi2/NDF=%1.1f/%d", set.c_str(),
	       f1->GetChisquare(), f1->GetNDF()) << endl << flush;

  // Functional form as reminder:
  // "[0]+[1]/(0.1*x)+[2]*log10(x)/(0.1*x)+[3]*(1+(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1))+[6]*pow(x,-0.3051)+[7]*(0.001*x)+[8]*pow((x+[9])/15.,2.)/(1+0.5*pow((x+[9])/15.,4.))",15,4500);
  // Set initial parameters carefully to ensure good convergence of the fit
  double k_10 = f1->Eval(10.);
  double k_20 = f1->Eval(20.); // => k10->k_20 better?
  double k_100 = f1->Eval(100.);
  double k_1500 = f1->Eval(1500.);
  double k_3000 = f1->Eval(3000.);
  f1raw->SetParameter(0, f1->GetParameter(0)); // 1
  f1raw->SetParameter(1, f1->GetParameter(1)/k_20 +
		      f1->GetParameter(2)/k_20*log10(k_20)); // 1/x
  f1raw->SetParameter(2, f1->GetParameter(2)/k_20); // log(x)/x
  f1raw->SetParameter(3, f1->GetParameter(3)); // "erf" heigh
  f1raw->SetParameter(4, f1->GetParameter(4)/k_3000); // "erf" mean
  f1raw->SetParameter(5, f1->GetParameter(5)); // "erf" steepness
  f1raw->SetParameter(6, f1->GetParameter(6)*pow(k_100,-0.3051)); // "hcal"
  f1raw->SetParameter(7, f1->GetParameter(7)*k_3000); // "trk" x
  // Last part has essentially fixed heigh vs JES, but peak position shifts
  f1raw->SetParameter(8, f1->GetParameter(8)); // "nhf_off" peak height
  //f1raw->SetParameter(9, 18.*(k_10-1)); // "nhf_off" peak position
  f1raw->SetParameter(9, k_20*f1->GetParameter(9)); // "nhf_off" peak position

  // Remove "erf" when no bias from ECALCC
  if (!isECALCC) {
    f1raw->FixParameter(3,0.);
    f1raw->FixParameter(4,1500.);
    f1raw->FixParameter(5,0.);
    f1raw->FixParameter(7,0.);
  }
  else {
    f1raw->SetParLimits(3, -0.2, 0.);
    f1raw->SetParLimits(4, 1000., 6500.);
    f1raw->SetParLimits(5, 1., 10.);
    f1raw->SetParLimits(7, -0.2, 0.2);
  }
  
  // To avoid very weird fits
  f1raw->SetParLimits(0, 0.5, 2.0);   // 1
  f1raw->SetParLimits(1, -0.5, +0.5); // 1/x
  f1raw->SetParLimits(2, -0.5, +0.5); // log(x)/x
  f1raw->SetParLimits(6, -0.5, +0.5); // x^-0.3
  f1raw->SetParLimits(8, 0., 0.5); // nhf_off
  f1raw->SetParLimits(9, 0.5, 2.0); // nhf_off peak

  
  _c3->cd();
  
  // Refit L3Res vs pT,raw for reference
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

  _leg3->SetTextSize(0.03);
  _leg3->SetY1(_leg3->GetY1()-0.03);
  _leg3->AddEntry(graw3,set.c_str(),"FL");
  
  cout << "f1raw->SetParameters(";
  ftxt << "f1raw->SetParameters(";
  const int npraw = 10;
  double p1raw[np];
  for (int i = 0; i != min(f1raw->GetNpar(),npraw); ++i) {
    p1raw[i] = f1raw->GetParameter(i);
    cout << (i==0 ? "" : ", ") << p1raw[i];
    if (i==1) ftxt << ",        ";
    ftxt << (i==0 ? "" : ", ") << Form("%+8.3f",p1raw[i]);
  } // for i
  cout << "); // Pre-set values" << endl;
  ftxt << "); // Pre-set values" << endl;
    
  graw3->Fit(f1raw,"QRN");
  graw3->Fit(f1raw,"QRNM");
  graw3->Fit(f1raw,"QRNM");
  
  f1raw->SetLineColor(color[set]);
  f1raw->Draw("SAME");

  cout << "f1raw->SetParameters(";
  ftxt << "f1raw->SetParameters(";
  for (int i = 0; i != min(f1raw->GetNpar(),npraw); ++i) {
    p1raw[i] = f1raw->GetParameter(i);
    cout << (i==0 ? "" : ", ") << p1raw[i];
    if (i==1) ftxt << ",        ";
    ftxt << (i==0 ? "" : ", ") << Form("%+8.3f",p1raw[i]);
  } // for i
  cout << Form("); // %s chi2/NDF=%1.1f/%d", set.c_str(),
	       f1raw->GetChisquare(), f1raw->GetNDF()) << endl << flush;
  ftxt << Form("); // %s chi2/NDF=%1.1f/%d", set.c_str(),
	       f1raw->GetChisquare(), f1raw->GetNDF()) << endl << flush;
  
  
  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  const char *run = set.c_str();

  string sin(""), sout2(""), sout3("");
  sin = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2ResidualVsPtRef_AK4PFPuppi.txt",cs);
  sout2 = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2L3ResidualVsPtRef_AK4PFPuppi.txt",cs);
  sout3 = Form("textFiles/Prompt24/Prompt24_Run%s_V8M_DATA_L2L3Residual_AK4PFPuppi.txt",cs);
  
  assert(sin!="");
  assert(sout2!="");
  assert(sout3!="");
  assert(sout2!=sin);
  assert(sout2!=sout3);
  assert(sout3!=sin);
						       
  // New code to merge L2Res and L3Res vs pTref, then re-map combination
  // to pT,raw. This needs to be done separately for each |eta| bin
  cout << "Running new code for merging L2Res+L3Res vs pTref," << endl
       << "then remapping to vs pT,raw for each eta bin" << endl;

  cout << "Reading in L2Residual file:" << endl
       << "   " << sin << endl << flush;  
  ifstream fin(sin.c_str());
  assert(fin.is_open());
  
  // Input L2Residual header
  string header, header2, header3;
  getline(fin, header);
  if (debug) cout << "Input L2Residual header:" << endl;
  if (debug) cout << header << endl;
  const int nparold = 2 + 3 + 3; // 8

  if (debug) cout << "Input L3Residual function(s):" << endl;
  if (debug) cout << "1./("<<f1->GetExpFormula().Data()<<")" << endl;
  if (debug) cout << "1./("<<f1raw->GetExpFormula().Data()<<")" << endl;

  // New combined function(s) and header(s)
  //string func = "[0]+[1]*log10(0.01*x)+[2]/(0.1*x)+[3]*log10(x)/(0.1*x)+[4]*(1+(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1))+[7]*pow(x,-0.3051)+[8]*(0.001*x)+[9]*pow(x/15.,2.)/(1+0.5*pow(x/15.,4.))";
  //string func3 = "[0]+[1]*log10(0.01*x)+[2]/(0.1*x)+[3]*log10(x)/(0.1*x)+[4]*(1+(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1))+[7]*pow(x,-0.3051)+[8]*(0.001*x)+[9]*pow((x+[10])/15.,2.)/(1+0.5*pow((x+[10])/15.,4.))";
  string func3 = "[0]+[1]*log10(0.01*x)+[2]/(0.1*x)+[3]*log10(x)/(0.1*x)+[4]*(1+(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1))+[7]*pow(x,-0.3051)+[8]*(0.001*x)+[9]*pow(x*[10]/15.,2.)/(1+0.5*pow(x*[10]/15.,4.))";

  header2 = "{1 JetEta 1 JetPt 1./("+func3+") Correction L2Relative";
  header3 = "{1 JetEta 1 JetPt 1./("+func3+") Correction L2Relative";
    
  cout << "Writing out L2L3Residual files:" << endl
       << "   " << sout2 << endl
       << "   " << sout3 << endl << flush;
  if (debug) cout << "Output L2L3Residual headers:" << endl;
  if (debug) cout << header2 << endl;
  if (debug) cout << header3 << endl;

  ofstream fout22ref(sout2.c_str());
  //const int nparnew = 2 + 3 + 10; // 15
  //fout << header << endl;
  fout22ref << header2 << endl;
  
  ofstream fout23raw(sout3.c_str());
  const int nparnew3 = 2 + 3 + 11; // 16
  fout23raw << header3 << endl;
    
  string line;
  double etamin(0), etamax(0);
  int npar(0), xmin(0), xmax(0), ptmin0(0), ptmax1(0);
  double p0(0), p1(0), p2(0);
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
    
    // Generate graphs with 100 points logarithmically distributed between
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
    
    //TGraph *gref = new TGraph();
    //TGraph *graw = new TGraph();
    TGraphErrors *gref = new TGraphErrors();
    TGraphErrors *graw = new TGraphErrors();
    TGraph *gr = new TGraph();
    double c = pow(ptrefmax2/ptrefmin2,1./(npt-1));
    for (int i = 0; i != npt; ++i) {

      double ptref = ptrefmin2 * pow(c, i);
      //double jesref = h->Interpolate(ptref);
      double jesref = f1->Eval(ptref);
      double jes = jesref * f2->Eval(ptref);
      double ptraw = ptref*jes;
      gref->SetPoint(i, ptref, jes);
      gref->SetPointError(i, ptref*0.001, jes*0.001);
      graw->SetPoint(i, ptraw, jes);
      graw->SetPointError(i, ptraw*0.001, jes*0.001);
      // Bias in current corrections when f1,f2 evaluated at wrong pt
      //double jesraw = h->Interpolate(ptraw);
      double jesraw = f1->Eval(ptraw);
      double jes2 = jesraw * f2->Eval(ptraw);
      gr->SetPoint(i, ptref, jes2/jes);
    }

    // Reference JES refit vs <pT,ref>:
    // "[0]+[1]*log10(x/100.)+[2]/(x/10.)+[3]*log10(x)/(x/10.)+"
    // "[4]*(1+(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1))+[7]*pow(x,-0.3051)+"
    // "[8]*(x/1000.)+[9]*pow((x*[10])/15.,2)/(1+0.5*pow((x*[10])/15.,4))"
    // // "[8]*(x/1000.)+[9]*pow((x*[10])/15.,2)/(1+0.5*pow((x*[10])/15.,4))"
    // f1/f1raw:
    // "[0]+[1]/(x/10.)+[2]*log10(x)/(x/10.)+"
    // "[3]*(1+(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1))+[6]*pow(x,-0.3051)"
    // "+[7]*(x/1000.)+[8]*pow((x*[9])/15.,2)/(1+0.5*pow((x*[9])/15.,4))"
    // // "+[7]*(x/1000.)+[8]*pow((x*[9])/15.,2)/(1+0.5*pow((x*[9])/15.,4))"
    // f2 (p1,p2,p3):
    // "[0]+[1]*log10(0.01*x)+[2]/(x/10.)"

    // First combined fit on L2Res*L3Res vs pTref so we factor out
    // issues caused by mapping on the x-axis and set better initial values
    // Three ranges:
    // 1) x=20 GeV for offset (1/x, log(x)/x, "nhf_off")
    // 2) x~300 GeV for absolute scale and "hcal" (-"erf")
    // 3) x~4000 GeV for log, "erf" and "trk" x
    // Consider product of F(x) = A + Bl + C/x and f(x) = a + b/x + cl/x + d(x)
    // => F(x)*f(x) = (A + Bl + C/x)*(a + b/x + cl/x + d(x))
    // Multiplying and regrouping can get
    // Aa + Ba*l + (Ab+Ca)*1/x + (Ac+Bb+B*c*l)*l/x + F(x)*d(x) + O(1/x^2)
    
    TF1 *f22 = new TF1(Form("f22_%s_%d",cs,cnt),func3.c_str(),
		       floor(ptrefmin2), ceil(ptrefmax2));
    /*
    // Works well for 2023D, works for 2024I
    f22->SetParameter(0, f2->GetParameter(0)*f1->GetParameter(0)); // 1
    f22->SetParameter(1, f2->GetParameter(1)*f1->GetParameter(0)); // log
    f22->SetParameter(2, f2->GetParameter(0)*f1->GetParameter(1) +
    		      f2->GetParameter(2)*f1->GetParameter(0)); // 1/x
    f22->SetParameter(3, f2->GetParameter(0)*f1->GetParameter(2) +
		      f2->GetParameter(1)*f1->GetParameter(1) +
		      f2->GetParameter(1)*f1->GetParameter(2)*log10(20.));// l/x
    */
    double k1_300 = f1->Eval(300.);
    double k2_300 = f2->Eval(300.);
    double k1_20 = f1->Eval(20.);
    double k2_20 = f2->Eval(20.);

    // This not working so well for 2024D
    f22->SetParameter(0, k2_300*f1->GetParameter(0)); // 1
    f22->SetParameter(1, k1_300*p1); // log(x)

    f22->SetParameter(2, k2_20*f1->GetParameter(1) + k1_20*p2); // 1/x
    f22->SetParameter(3, k2_20*f1->GetParameter(2) + k1_20*p1); // log(x)/x
    
    double k2_4000 = f2->Eval(4000.);
    // "erf" for "ecalcc"
    f22->SetParameter(4, k2_4000*f1->GetParameter(3));
    f22->SetParameter(5, f1->GetParameter(4));
    f22->SetParameter(6, f1->GetParameter(5));

    // x^-0.3 + x parts for "hcal" + "trk"
    //double k2_300 = f2->Eval(300.);
    f22->SetParameter(7, k2_300*f1->GetParameter(6));
    f22->SetParameter(8, k2_4000*f1->GetParameter(7));

    //double k2_20 = f2->Eval(20.);
    f22->FixParameter(9, k2_20*f1->GetParameter(8));
    f22->FixParameter(10, k2_20*f1->GetParameter(9));

    // Limitations on ecalcc shapes
    if (!isECALCC) {
      f22->FixParameter(4, 0.);
      f22->FixParameter(5, 1500.);
      f22->FixParameter(6, 0.);
      f22->FixParameter(8, 0.);
    }
    else {
      f22->SetParameter(4, max(-0.2,min(0.2,f22->GetParameter(4))));
      f22->SetParLimits(4, -0.2, 0.2);
      f22->SetParameter(5, max(1000.,min(2000.,f22->GetParameter(5))));
      f22->SetParLimits(5, 1000., 2000.);
      f22->SetParameter(6, max(1.,min(5.,f22->GetParameter(6))));
      f22->SetParLimits(6, 1., 5.);
      f22->SetParameter(7, max(-0.2,min(0.2,f22->GetParameter(7))));
      f22->SetParLimits(7, -0.2, 0.2);
      f22->SetParameter(8, max(-0.2,min(0.2,f22->GetParameter(8))));
      f22->SetParLimits(8, -0.2, 0.2);
    }

    // To avoid very weird fits
    f22->SetParameter(0, max(0.5,min(2.0,f22->GetParameter(0))));
    f22->SetParLimits(0, 0.5, 2.0);   // 1
    f22->SetParameter(1, max(-0.25,min(0.25,f22->GetParameter(1))));
    f22->SetParLimits(1, -0.25, +0.25);   // log(x)
    f22->SetParameter(2, max(-0.5,min(0.5,f22->GetParameter(2))));
    f22->SetParLimits(2, -0.5, +0.5); // 1/x
    f22->SetParameter(3, max(-0.5,min(0.5,f22->GetParameter(3))));
    f22->SetParLimits(3, -0.5, +0.5); // log(x)/x
    f22->SetParameter(7, max(-0.5,min(0.5,f22->GetParameter(7))));
    f22->SetParLimits(7, -0.5, +0.5); // x^-0.3
    f22->SetParameter(9, max(0.,min(0.5,f22->GetParameter(9))));
    f22->SetParLimits(9, 0., 0.5); // nhf_off
    f22->SetParameter(10, max(0.5,min(2.0,f22->GetParameter(10))));
    f22->SetParLimits(10, 0.5, 2.0); // nhf_off peak
      
    // Input parameters
    if (etamin>=0) {
      ftxt << "  f22->SetParameters(";
      for (int ip = 0; ip != f22->GetNpar(); ++ip) {
	ftxt << (ip==0 ? "" : ", ") << Form("%+8.3f",f22->GetParameter(ip));
      }
      ftxt << "); // Pre-set values" << endl;
    }
    
    gref->Fit(f22,"QRN");
    gref->Fit(f22,"QRNM");
    gref->Fit(f22,"QRNM");
    
    // Output L2L3Res
    double p22ref[np];
    assert(np>=f22->GetNpar());
    for (int ip = 0; ip != min(f22->GetNpar(),np); ++ip) {
      p22ref[ip] = f22->GetParameter(ip);
    }
    if (etamin>=0) {
      ftxt << "  f22->SetParameters(";
      for (int ip = 0; ip != f22->GetNpar(); ++ip) {
	ftxt << (ip==0 ? "" : ", ") << Form("%+8.3f",f22->GetParameter(ip));
      }
      ftxt << Form("); // %s chi2/NDF=%1.1f/%d [%1.3f,%1.3f]\n", set.c_str(),
		   f22->GetChisquare(), f22->GetNDF(),etamin,etamax) << flush;
    }

    // Calculate ratio of fit and graph to check quality
    TGraph *gr22 = new TGraph(gref->GetN());
    for (int i = 0; i != npt; ++i) {
      double ptref = gref->GetX()[i];
      double jes1 = gref->GetY()[i];
      double jes2 = f22->Eval(ptref);
      gr22->SetPoint(i, ptref, 1+(jes2/jes1-1)*10.);
    }

    // Reference JES refit vs <pT,raw>:
    // "[0]+[1]*log10(x/100.)+[2]/(x/10.)+[3]*log10(x)/(x/10.)+"
    // "[4]*(1+(pow(x/[5],[6])-1)/(pow(x/[5],[6])+1))+[7]*pow(x,-0.3051)+"
    // "[8]*(x/1000.)+[9]*pow((x+[10])/15.,2)/(1+0.5*pow((x+[10])/15.,4))"
    
    //double jesrefmin2  = h->Interpolate(ptrefmin2);
    //double jesrefmax2 = h->Interpolate(ptrefmax2);
    double jesrefmin2  = f1->Eval(ptrefmin2);
    double jesrefmax2 = f1->Eval(ptrefmax2);
    double jesmin2  = jesrefmin2*f2->Eval(ptrefmin2);
    double jesmax2 = jesrefmax2*f2->Eval(ptrefmax2);
    double ptminraw = ptrefmin2*jesmin2;
    double ptmaxraw = ptrefmax2*jesmax2;
    const int nparnew = 11;//9;

    //double k2_10 = f2->Eval(10.);
    //double k1_10 = f1->Eval(10.);
    //double kjes_10 = k2_10*k1_10;

    TF1 *f23 = new TF1(Form("f23_%s_%d",cs,cnt),func3.c_str(),
		       floor(ptminraw), ceil(ptmaxraw));

    //double k1_100 = h->Interpolate(100.);
    //double k1_300 = f1->Eval(300.);
    //double k2_300 = f2->Eval(300.);
    double kjes_300 = k1_300*k2_300;
    f23->SetParameter(0, k2_300*f1->GetParameter(0)); // 1
    f23->SetParameter(1, k1_300*p1); // log10(x/100.)

    //double k1_20 = h->Interpolate(20.);
    //double k1_20 = f1->Eval(20.);
    //double k2_20 = f2->Eval(20.);
    double kjes_20 = k2_20*k1_20;
    f23->SetParameter(2, k2_20*(f1->GetParameter(1)/kjes_20 +
				f1->GetParameter(2)/kjes_20*log10(kjes_20)) +
		      k1_20*(p2/kjes_20 +
			     p2/kjes_20*log10(kjes_20))); // 1/(x/10)
    f23->SetParameter(3, k2_20*f1->GetParameter(2)/kjes_20 +
		      k1_20*p2/kjes_20); // log10(x)/(x/10)

    //double k1_1500 = h->Interpolate(1500.);
    //double k1_1500 = f1->Eval(1500.);
    //double k2_1500 = f2->Eval(1500.);
    //double kjes_1500 = k2_1500*k1_1500;
    // "erf" shape for "ecalcc" 
    double k1_4000 = f1->Eval(4000.);
    //double k2_4000 = f2->Eval(4000.);
    double kjes_4000 = k2_4000*k1_4000;
    f23->SetParameter(4, k2_4000*f1->GetParameter(3)); // [-1,1] "erf"
    f23->SetParameter(5, f1->GetParameter(4)/kjes_4000); // "erf" position
    f23->SetParameter(6, f1->GetParameter(5)); // "erf" steepness

    // "hcal" x^=0.3 and "trk" x shapes
    f23->SetParameter(7, k2_300*f1->GetParameter(6)*pow(kjes_300,-0.3051));
    f23->SetParameter(8, k2_4000*f1->GetParameter(7)*kjes_4000); // x

    // Last part has essentially fixed heigh vs JES, but peak position shifts
    f23->SetParameter(9, k2_20*f1->GetParameter(8));
    //f23->SetParameter(10, 20*(kjes_20-1));
    f23->SetParameter(10, kjes_20);
    
    // To avoid very weird fits
    //f23->SetParameter(0, max(0.3,min(0.25,f23->GetParameter(0))));
    //f23->SetParLimits(0,0.3,2.5);
    
    // Some sanity check limits
    //f23->SetParameter(9, max(0., min(0.25, f23->GetParameter(9))));
    //f23->SetParLimits(9, 0., 0.25);
    //f23->SetParameter(10, max(-9., min(+18., f23->GetParameter(10))));
    //f23->SetParLimits(10, -9, +18);
    
    // To avoid division by zero errors
    //f23->SetParameter(5, max(10.,min(6500.,f23->GetParameter(5))));
    //f23->SetParLimits(5, 10., 6500.);
    //f23->SetParameter(6, max(0.,min(10.,f23->GetParameter(6))));
    //f23->SetParLimits(6, 0., 10.);
    //if (f23->GetNpar()>8)
    //f23->SetParLimits(10,-10,+10.);

    //if (isECALCC) {
    //f23->SetParameter(5, 1500);
    //f23->SetParLimits(5, 1000., 2000.);
    //}
    if (!isECALCC) {
      f23->FixParameter(4, 0.);
      f23->FixParameter(5, 1500.);
      f23->FixParameter(6, 0.);
      f23->FixParameter(8, 0.);
    }
    else {
      f23->SetParameter(4, max(-0.2,min(0.2,f23->GetParameter(4))));
      f23->SetParLimits(4, -0.2, 0.2);
      f23->SetParameter(5, max(1000.,min(2000.,f23->GetParameter(5))));
      f23->SetParLimits(5, 1000., 2000.);
      f23->SetParameter(6, max(1.,min(5.,f23->GetParameter(6))));
      f23->SetParLimits(6, 1., 5.);
      f23->SetParameter(8, max(-0.2,min(0.2,f23->GetParameter(8))));
      f23->SetParLimits(8, -0.2, 0.2);
    }

    // To avoid very weird fits
    f23->SetParameter(0, max(0.5,min(2.0,f23->GetParameter(0))));
    f23->SetParLimits(0, 0.5, 2.0);   // 1
    f23->SetParameter(1, max(-0.25,min(0.25,f23->GetParameter(1))));
    f23->SetParLimits(1, -0.25, +0.25); // log(x)
    f23->SetParameter(2, max(-0.5,min(0.5,f23->GetParameter(2))));
    f23->SetParLimits(2, -0.5, +0.5); // 1/x
    f23->SetParameter(3, max(-0.5,min(0.5,f23->GetParameter(3))));
    f23->SetParLimits(3, -0.5, +0.5); // log(x)/x
    f23->SetParameter(7, max(-0.5,min(0.5,f23->GetParameter(7))));
    f23->SetParLimits(7, -0.5, +0.5); // x^-0.3
    f23->SetParameter(9, max(0.,min(0.5,f23->GetParameter(9))));
    f23->SetParLimits(9, 0., 0.5); // nhf_off
    f23->SetParameter(10, max(0.5,min(2.0,f23->GetParameter(10))));
    f23->SetParLimits(10, 0.5, 2.0); // nhf_off peak
    
    string n23[nparnew] =
      {"[0]","[1]*lnx","[2]/x","[3]*l/x","[4]*x^6",
       "x/[5]^6","x^[6]","[7]/x^.3","[8]*x","[9]*x^2",
       "(x+[10])^2"};
    
    // Input parameters
    if (etamin==0) {
      ftxt << f23->GetExpFormula() << endl;
      ftxt << "  f23->SetParNames  (";
      for (int ip = 0; ip != f23->GetNpar(); ++ip) {
	ftxt << (ip==0 ? "" : ", ") << Form("%8s",n23[ip].c_str());
      }
      ftxt << endl;
    }
    
    double p23raw[nparnew];
    assert(nparnew>=f23->GetNpar());
    for (int ip = 0; ip != f23->GetNpar(); ++ip) {
      p23raw[ip] = f23->GetParameter(ip);
    }
    
    // Input parameters
    if (etamin>=0) {
      ftxt << "  f23->SetParameters(";
      for (int ip = 0; ip != f23->GetNpar(); ++ip) {
	ftxt << (ip==0 ? "" : ", ") << Form("%+8.3f",p23raw[ip]);
      }
      ftxt << "); // Pre-set values" << endl;
    }
    
    graw->Fit(f23,"QRN");
    graw->Fit(f23,"QRNM");
    graw->Fit(f23,"QRNM");
    
    // Output L2L3Res
    if (etamin>=0) ftxt << "  f23->SetParameters(";
    for (int ip = 0; ip != f23->GetNpar(); ++ip) {
      p23raw[ip] = f23->GetParameter(ip);
      if (etamin>=0)
	ftxt << (ip==0 ? "" : ", ") << Form("%+8.3f",p23raw[ip]);
    }
    if (etamin>=0)
      ftxt << Form("); // %s chi2/NDF=%1.1f/%d [%1.3f,%1.3f]\n", set.c_str(),
		   f23->GetChisquare(), f23->GetNDF(),etamin,etamax) << flush;
    
    string s23raw = Form("  %6.3f %6.3f %2d   %2d %4d   "
			 "%6.4f %+7.4f   %+7.4f %+7.4f   "
			 "%7.4f %5.1f %6.4f   %+7.3f   %+7.4f   "
			 "%+7.3f   %+7.4f",
			 etamin, etamax, nparnew, int(ptminraw), int(ptmaxraw),
			 p23raw[0],p23raw[1], p23raw[2],p23raw[3],
			 p23raw[4],p23raw[5],p23raw[6],  p23raw[7], p23raw[8],
			 p23raw[9],p23raw[10]);
    if (cnt<cntmax && debug)
      cout << s23raw << endl;
    fout23raw << s23raw << endl;
    ++cnt;

    string s22ref = Form("  %6.3f %6.3f %2d   %2d %4d   "
			 "%6.4f %+7.4f   %+7.4f %+7.4f   "
			 "%7.4f %5.1f %6.4f   %+7.3f   %+7.4f   "
			 "%+7.3f   %+7.4f",
			 etamin, etamax, nparnew, int(ptminraw), int(ptmaxraw),
			 p22ref[0],p22ref[1], p22ref[2],p22ref[3],
			 p22ref[4],p22ref[5],p22ref[6],  p22ref[7], p22ref[8],
			 p22ref[9],p22ref[10]);
    fout22ref << s22ref << endl;
    
    
    // Plotting only for positive side (is anyway symmetric)
    if (etamin<0) continue;
    ++ieta;
    
    // Calculate ratio of fit and graph to check quality
    TGraph *gr23 = new TGraph(graw->GetN());
    for (int i = 0; i != npt; ++i) {
      double ptref = gref->GetX()[i];
      double ptraw = graw->GetX()[i];
      double jes1 = graw->GetY()[i];
      double jes2 = f23->Eval(ptraw);
      gr23->SetPoint(i, ptref, 1+(jes2/jes1-1)*10.);
    }
    
    cx->cd(ieta);
    double eps = 1e-4;
    TH1D *hx = tdrHist(Form("hx_%s_%d",cs,ieta),"Rel. JES Data/MC refit",
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
    
    tdrDraw(gref, "Pz", kOpenCircle, kGray+1, kSolid, -1, kNone, 0, 0.5);
    tdrDraw(graw, "Pz", kOpenCircle, kBlack, kSolid, -1, kNone, 0, 0.5);
    tdrDraw(gr, "L", kNone, kGreen+2, kSolid, -1, kNone, 0, 0.5);
    tdrDraw(gr22, "L", kNone, kMagenta-9, kSolid, -1, kNone, 0, 0.5);
    tdrDraw(gr23, "L", kNone, kMagenta+1, kSolid, -1, kNone, 0, 0.5);

    f22->SetLineColor(kOrange+1);
    f22->Draw("SAME");
    
    f23->SetLineColor(kRed);
    f23->Draw("SAME");
    
    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045*2.0);
    tex->DrawLatex(0.50,0.90,Form("[%1.3f,%1.3f]",etamin,etamax));
    if (ieta==nxy) {
      cx->cd(ieta+1);
      TLegend *legx = tdrLeg(0.05,0.95-0.05*2.0*9,0.30,0.95);
      legx->SetTextSize(0.045*2.0);//2.5);
      legx->AddEntry(f1,"JES-L3Res fit vs p_{T,ref}","L");
      legx->AddEntry(f2,"JES-L2Res fit vs p_{T,ref}","L");
      legx->AddEntry(gref,"JES vs p_{T,ref}","P");
      legx->AddEntry(f22,"Fit-L2L3Res vs p_{T,ref}","L");
      legx->AddEntry(graw,"JES vs p_{T,raw}","P");
      legx->AddEntry(f23,"Fit-L2L3Res vs p_{T,raw}","L");
      legx->AddEntry(gr,"Ratio of JES vs p_{T,ref}","L");
      legx->AddEntry(gr22,"Ratio of fit/JES vs p_{T,ref} #times 10","L");
      legx->AddEntry(gr23,"Ratio of fit/JES vs p_{T,raw} #times 10","L");
      
      cx->cd(ieta);
      
       #include "Config.C"

      tex->DrawLatex(0.50,0.65,Form("%s vs",cs));
      tex->DrawLatex(0.50,0.50,"Winter 24");
      tex->DrawLatex(0.50,0.35,mlum[cs].c_str());
      //tex->SetTextSize(siz);
    }
  } // for ieta
  
  cx->SaveAs(Form("pdf/createL2L3ResTextFile/createL2L3ResTextFile_PtRawReFit_%s.pdf",cs));

  f1->SetLineColor(color[set]);
  f1->SetRange(15,4500.);
} // createL2L3ResTextFile
