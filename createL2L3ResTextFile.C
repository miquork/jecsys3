// Purpose: Create L2L3Res text file with simple parameterization
//          Takes as input previous L2Res and complex 9p global JES fit
//          Outputs same L2Res and "simple" 7p fit
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

void createL2L3ResTextFiles(string set="RunC");

TLegend *_leg(0);
void createL2L3ResTextFile() {

  setTDRStyle();

  double ptmin = 15;
  double ptmax = 4500;
  TH1D *h = tdrHist("h","Absolute response at |#eta| < 1.3",
		    0.88+1e-4,1.05-1e-4,"p_{T} (GeV)",ptmin,ptmax);
  lumi_136TeV = "X.X fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  c1->SetLeftMargin(0.17);
  c1->SetRightMargin(0.03);
  h->SetTitleOffset(1.5,"Y");
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(ptmin,1,ptmax,1);

  _leg = tdrLeg(0.45,0.90,0.75,0.90);

  createL2L3ResTextFiles("RunC");

  c1->Update();
  c1->SaveAs("pdf/createL2L3ResTextFile_2022C.pdf");
}

void createL2L3ResTextFiles(string set) {

  //if (debug) 
  cout << "Warning: sscanf only works correctly when code is compiled (.C+)\n";

  cout << "Processing " << set << endl << flush;

  // Simplify complex sum into an effective formula
  // Need good starting values and/or a few iterations to converge
  // Fit done to hjesfit from each IOV
  TDirectory *curdir = gDirectory;
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",set.c_str()),"READ");

  assert(f && !f->IsZombie());
  TH1D *h(0);
  //h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit2");
  //if (!h) h = (TH1D*)f->Get("ratio/eta00-13/sys/hjesfit");
  if (!h) h = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet");
  assert(h);
  curdir->cd();

  TF1 *f1 = new TF1(Form("f1_%s",set.c_str()),"[0]+[1]/x+[2]*log(x)/x+[3]*(pow(x/[4],[5])-1)/(pow(x/[4],[5])+1)+[6]*pow(x,-0.3051)",15,4500);
  f1->SetParameters(0.98, 0.1,0.01, 0.01,500.,1.3, 0.001);
  if (set=="2017H") 
    f1->SetParameters(0.99, 1.5,0.01, 0.01,1000.,1.3, 0.001);

  map<string,int> color;
  color["RunC"] = kYellow+2;

  h->Fit(f1,"QRN");
  h->Fit(f1,"QRNM");
  h->Fit(f1,"QRNM");
  tdrDraw(h,"LE3",kNone,color[set],kSolid,-1,1001,color[set]-9);
  h->SetFillColorAlpha(color[set]-9,0.7);
  f1->SetLineColor(color[set]);
  f1->Draw("SAME");

  _leg->SetY1(_leg->GetY1()-0.05);
  _leg->AddEntry(h,set.c_str(),"FL");

  const int np = 7;
  double p[np];
  for (int i = 0; i != np; ++i) {
    p[i] = f1->GetParameter(i);
  } // for i


  /////////////////////////////////////////////////////////////
  // Generate input and output file names semi-automatically  //
  //////////////////////////////////////////////////////////////
  const char *run = set.c_str();

  string sin, sout;
  sin = "../JERCProtoLab/Winter22Run3/L2Residual/Winter22Run3_V1/Winter22Run3_V1_MPF_LOGLIN_L2Residual_pythia8_AK4PFPuppi.txt";
  sout = "../JERCProtoLab/Winter22Run3/global_fit/Winter22Run3_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.txt";

  ifstream fin(sin.c_str());
  assert(fin.is_open());
  ofstream fout(sout.c_str());

  // 2018: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsOctober2018

  cout << "Reading in L2Residual file:" << endl
       << "   " << sin << endl;
  cout << "Writing out L2L3Residual file:" << endl
       << "   " << sout << endl;

  string header;
  getline(fin, header);
  if (debug) cout << "Old L2L3Residual header:" << endl;
  if (debug) cout << header << endl;

  header = "{ 1 JetEta 1 JetPt [2]*([3]*([4]+[5]*TMath::Log(max([0],min([1],x))))*1./([6]+[7]/x+[8]*log(x)/x+[9]*(pow(x/[10],[11])-1)/(pow(x/[10],[11])+1)+[12]*pow(x,-0.3051))) Correction L2Relative}";
  if (debug) cout << "New L2L3Residual header:" << endl;
  if (debug) cout << header << endl;
  fout << header << endl;
  
  string line;
  double etamin, etamax;
  int npar, xmin, xmax, ptmin0, ptmax1;
  double p2, p3, p4, p5;
  double p6, p7, p8;
  int cnt(0); int cntmax(0);
  while (getline(fin,line)) {
    if (cnt<cntmax && debug) cout << line << endl;
    assert(sscanf(line.c_str(),"%lf %lf  %d  %d %d  %d %d  %lf %lf %lf %lf"
		  "  %lf %lf %lf",
		  &etamin, &etamax, &npar, &xmin, &xmax, &ptmin0, &ptmax1,
		  &p2, &p3, &p4, &p5,  &p6, &p7, &p8)==14);
    assert(!(etamin==0 && etamax==0));
    if (fabs(etamin)<0.01 && fabs(etamax)<0.01) {
      cout << "sscanf failed!" << endl << flush;
      exit(1);
    }
    int nparnew = 15;
    if (cnt<cntmax && debug)
      cout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		   "   %8.6f %8.6f   "
		   "%5.4f %5.4f %5.5f %5.5f %5.1f %5.4f %5.5f",
		   etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		   p2, p3, p4, p5,
		   p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
    fout << Form("  %9.6f %9.6f   %d   %d %d   %d   %d   %8.6f %8.6f"
		 "   %8.6f %8.6f   "
		 "%5.4f %5.4f %5.5f %5.5f %5.1f %5.4f %5.5f",
		 etamin, etamax, nparnew, xmin, xmax, ptmin0, ptmax1,
		 p2, p3, p4, p5,
		 p[0], p[1], p[2], p[3], p[4], p[5], p[6]) << endl;
    ++cnt;
  }

} // creataL2L3ResTextFile
