// Purpose: Merge L2Res(pTref;eta) and L3Res(pTref) into L2L3Res(pTL2;eta),
//          where <pTL2*L2L3Res(pTL2;eta)> = <pTgen> = DBref*pTref and
//          DBref = <pTgen>/<pTref>
//
//          Step 1. Map L2Res and L3Res vs pTref to L2L3Res vs pTL2 numerically
//          Step 2. Re-fit function with effective physical JES parameterization
//          Step 3. Monitor remapping-fit difference to 0.1% level
//
//          Some considerations on effective physical JES:
//          1. Foundation: offset (1/x), scale (const), slope (x^-0.3)
//            1b. Offset mods: log(x)/x, 1/x^2
//            2b. Slope mods: log(x), log(x)^2
//          2. Extras: PU efficiency/noise (erf(x)), high-pT tracking (x)
//          3. Iterative fitting to layer more complex shapes on top of baseline
//
//          Also account for non-linear L2L3Res producing average correction
//          different from single value at <pTL2>?
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>
#include <fstream>

#include "tdrstyle_mod22.C"

bool debug = false;

void createL2L3ResTextFileV2s(string run = "2025CDEFG");

void createL2L3ResTextFileV2() {

  createL2L3ResTextFileV2s("2024CDE_nib");
  createL2L3ResTextFileV2s("2024FGHI_nib");
  createL2L3ResTextFileV2s("2025CDEFG");

  createL2L3ResTextFileV2s("2024C_nib1");
  createL2L3ResTextFileV2s("2024D_nib1");
  createL2L3ResTextFileV2s("2024E_nib1");
  createL2L3ResTextFileV2s("2024F_nib1");
  createL2L3ResTextFileV2s("2024F_nib2");
  createL2L3ResTextFileV2s("2024F_nib3");
  createL2L3ResTextFileV2s("2024G_nib1");
  createL2L3ResTextFileV2s("2024G_nib2");
  createL2L3ResTextFileV2s("2024H_nib1");
  createL2L3ResTextFileV2s("2024I_nib1");

  createL2L3ResTextFileV2s("2025C");
  createL2L3ResTextFileV2s("2025E");
  createL2L3ResTextFileV2s("2025D");
  createL2L3ResTextFileV2s("2025F");
  createL2L3ResTextFileV2s("2025G");
  
  createL2L3ResTextFileV2s("2026B");
  createL2L3ResTextFileV2s("2026C");
  createL2L3ResTextFileV2s("2026D");
  
} // createL2L3ResTextFileV2

void createL2L3ResTextFileV2s(string run) {

  cout << endl;
  cout << "******************************\n";
  cout << "** Processing " << run << " **\n";
  cout << "******************************\n";
  cout << flush;
  
  TString tr(run.c_str());

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/createL2L3ResTextFileV2");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/createL2L3ResTextFileV2");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  ///////////////////
  // Read in L2Res //
  ///////////////////

  // Inputs stored in these three vectors
  vector<double> vetamin, vetamax;
  vector<TF1*> vfl2res;

  // Reading in L2Res from text file
  // Alternative is using rootfiles/L2Res.root and TH2D "h2jes1_const_eta_%s"
  // but this loses detailed information on original functional form
  
  string sin(""), sout("");
  string path("textFiles/Prompt");
  string name2("DATA_L2ResidualVsPtRefAsymm_AK4PFPuppi");
  string name3("DATA_L2L3Residual_AK4PFPuppi");
  const char *cp = path.c_str();
  const char *cr = run.c_str();
  const char *cn2 = name2.c_str();
  const char *cn3 = name3.c_str();

  // Input L2Residual files
  if (tr.Contains("24")) sin = Form("%s/Prompt24_Run%s_V10M_%s.txt",cp,cr,cn2);
  if (tr.Contains("25")) sin = Form("%s/Prompt25_Run%s_V4M_%s.txt",cp,cr,cn2);
  if (tr.Contains("26")) sin = Form("%s/Prompt26_Run%s_V1M_%s.txt",cp,cr,cn2);

  // Output L2L3Residual files
  if (tr.Contains("24")) sout = Form("%s/Prompt24_Run%s_V10M_%s.txt",cp,cr,cn3);
  if (tr.Contains("25")) sout = Form("%s/Prompt25_Run%s_V4M_%s.txt",cp,cr,cn3);
  if (tr.Contains("26")) sout = Form("%s/Prompt26_Run%s_V1M_%s.txt",cp,cr,cn3);
  
  cout << "Reading in L2Residual file:" << endl
       << "   " << sin << endl << flush;  
  ifstream fin(sin.c_str());
  assert(fin.is_open());

  string header, line;
  getline(fin, header);
  if (debug) cout << "Input L2Residual header:" << endl;
  if (debug) cout << header << endl;

  int nbinvar(0), nparvar(0);
  string func2;
  char cfunc2[512];
  int nread(0);
  nread=sscanf(header.c_str(),"{ %d JetEta %d JetPt %s Correction L2Relative}",
	       &nbinvar, &nparvar, cfunc2);
  assert(nread==3);
  if (debug) cout << Form("nbinvar=%d, nparvar=%d, string=%s",
			  nbinvar, nparvar, cfunc2) << endl;

  const int npar2(5);
  int npar, xmin, xmax;
  double etamin, etamax;
  double p0, p1, p2, p3, p4;
  while (getline(fin,line)) {
    nread=sscanf(line.c_str(),"%lf %lf  %d  %d %d  %lf %lf %lf %lf %lf",
		 &etamin, &etamax, &npar, &xmin, &xmax,
		 &p0, &p1, &p2, &p3, &p4);
    assert(nread==2+3+5);   // sscanf succeeded
    assert(npar==npar2+2);  // npar matches expectation
    //sscanf();

    TF1 *f2 = new TF1(Form("f2_%+04dto%+04d_%s",
			   int(etamin*1000),int(etamax*1000),cr),
		      cfunc2,10,4500);
    f2->SetParameters(p0,p1,p2,p3,p4);

    vetamin.push_back(etamin);
    vetamax.push_back(etamax);
    vfl2res.push_back(f2);
  }
  
  //////////////////////////
  // Read in L3Res and DB //
  //////////////////////////

  string sin3 = Form("rootfiles/jecdata%s.root",cr);
  cout << "Reading in L3Residual file:" << endl
       << "   " << sin3 << endl << flush;  

  TFile *f = new TFile(sin3.c_str(),"READ");
  assert(f && !f->IsZombie());

  curdir->cd();
  
  TH1D *hl3jes = (TH1D*)f->Get("ratio/eta00-13/run3/hFit_Rjet"); assert(hl3jes);



  ////////////////////////////////////////////////////////////
  // Merge L2Res and L2L3Res into a single function vs pTL2 //
  ////////////////////////////////////////////////////////////

  // Inclusive jets pT binning
  double vx[] =
    {5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
     97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
     507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
     1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
     2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
     4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const int nx = sizeof(vx)/sizeof(vx[0]);

  // Listing of increasingly complex functional forms for iterative fitting
  const int npar3(9); // maximum number of fit parameters
  vector<string> vfit;
  vfit.push_back("[0]+[1]/x+[2]*pow(x,-0.3)");
  vfit.push_back("[0]+[1]/x+[2]*pow(x,-0.3)+[3]*log(0.01*x)");
  vfit.push_back("[0]+[1]/x+[2]*pow(x,-0.3)+[3]*log(0.01*x)+[4]*pow(log(0.01*x),2)");
  vfit.push_back("[0]+[1]/x+[2]*pow(x,-0.3)+[3]*log(0.01*x)+[4]*pow(log(0.01*x),2)+[5]/pow(0.1*x,2)");
  vfit.push_back("[0]+[1]/x+[2]*pow(x,-0.3)+[3]*log(0.01*x)+[4]*pow(log(0.01*x),2)+[5]/pow(0.1*x,2)+[6]*(1-erf((log(x)-log([7]*45.))/log([8]*1.67)))");
  //const char *cfunc = func.c_str();

  // Set initial values and fit ranges for each parameter
  map<int, double> mparmin;
  map<int, double> mparmax;
  map<int, double> mparset;
  mparmin[0] =  0.8;  mparmax[0] =  1.2; mparset[0] = 1.00; // x^0
  mparmin[1] = -4.0;  mparmax[1] = +8.0; mparset[1] = +0.1; // x^-1
  mparmin[2] = -0.7;  mparmax[2] = +0.7; mparset[2] = -0.1; // x^-0.3
  mparmin[3] = -0.2;  mparmax[3] = +0.2; mparset[3] = 0.01; // log(x)
  mparmin[4] = -0.2;  mparmax[4] = +0.2; mparset[4] = 0.01; // log(x)^2
  mparmin[5] = -0.7;  mparmax[5] = +1.4; mparset[5] = 0.01; // x^-2
  mparmin[6] =  0.0;  mparmax[6] =  0.2; mparset[6] = 0.10; // c*erf
  mparmin[7] =  0.3;  mparmax[7] =  1.5; mparset[7] = 1.00; // erf(mu)
  mparmin[8] =  0.5;  mparmax[8] =  2.0; mparset[8] = 1.00; // erf(sigma)
		 
  vector<TGraphErrors*> vgl2l3res(vetamin.size());
  vector<TGraphErrors*> vgdiff(vetamin.size());
  vector<TF1*> vfl2l3res(vetamin.size());
  TGraphErrors *gdmax = new TGraphErrors(vetamin.size());
  TGraphErrors *gd30max = new TGraphErrors(vetamin.size());
  for (int ieta = 0; ieta != int(vetamin.size()); ++ieta) {
    
    double eta = 0.5*(vetamin[ieta]+vetamax[ieta]);
    double etamin = vetamin[ieta];
    double etamax = vetamax[ieta];
    TF1 *f2 = vfl2res[ieta];
    
    TGraphErrors *g = new TGraphErrors();
    g->SetName(Form("g_%+04dto%+04d_%s",int(etamin*1000),int(etamax*1000),cr));
    double ptmin(7000), ptmax(0.);
    for (int i = 0; i != nx-1; ++i) {
      
      double ptref = 0.5*(vx[i]+vx[i+1]);
      double dptref = 0.5*(vx[i+1]-vx[i+1]);
      double eref = ptref*cosh(eta);
      if (ptref>5000. || ptref<12. || eref>0.5*13600.) continue;
      
      double l2res = f2->Eval(ptref);
      double l2jes = (l2res!=0 ? 1./l2res : 0);
      double l3jes = hl3jes->Interpolate(ptref);
      double jes = l2jes*l3jes;

      const double db = 0.9;
      double ptgen = db * ptref;
      double ptl2 = jes * ptgen;
      double ptl2min = ptl2 * vx[i]/ptref;
      double ptl2max = ptl2 * vx[i+1]/ptref;

      double eptl2 = sqrt(pow(0.5*fabs(1-db)*ptl2,2)+pow(dptref,2));
      double ejes = sqrt(pow(0.1/ptref,2) + pow(0.001,2) + pow(0.001*ptref/1000.,2));

      int n = g->GetN();
      g->SetPoint(n, ptl2, jes);
      double eptfactor(0.25); // tune to constrain fit vs pt better
      double efactor(0.5); // tune to force fit follow points better
      g->SetPointError(n, eptfactor*eptl2, efactor*ejes);
      if (ptl2min<ptmin) ptmin = ptl2min;
      if (ptl2max>ptmax) ptmax = ptl2max;
    } // for i
    vgl2l3res[ieta] = g;

    TF1 *f23prev(0), *f23best(0);
    double chi2min(9999.);
    for (int i = 0; i != int(vfit.size()); ++i) {

      // Try both fit iterative fit and a fresh fit; keep better one
      const char *cfunc = vfit[i].c_str();
      TF1 *f23a = new TF1(Form("f23a_%+04dto%+04d_%d_%s",
			       int(etamin*1000),int(etamax*1000),i,cr),
			  cfunc,ptmin,ptmax);
      TF1 *f23b = new TF1(Form("f23b_%+04dto%+04d_%d_%s",
			       int(etamin*1000),int(etamax*1000),i,cr),
			  cfunc,ptmin,ptmax);
      int npar = f23a->GetNpar();
      for (int j = 0; j != npar; ++j) {
	f23a->SetParameter(j, mparset[j]);
	f23a->SetParLimits(j, mparmin[j], mparmax[j]);
	f23b->SetParameter(j, mparset[j]);
	f23b->SetParLimits(j, mparmin[j], mparmax[j]);
      } // for j
      if (f23prev!=0) {
	for (int j = 0; j != f23prev->GetNpar(); ++j) {
	  f23b->SetParameter(j, f23prev->GetParameter(j));
	} // for j
      }
      g->Fit(f23a,"QRN");
      g->Fit(f23b,"QRN");

      double chi2a = f23a->GetChisquare() / max(f23a->GetNDF(),1);
      double chi2b = f23b->GetChisquare() / max(f23a->GetNDF(),1);
      if (f23a->GetNDF()<1) chi2a = 9999.;
      if (f23b->GetNDF()<1) chi2b = 9999.;
      if (f23best==0 || chi2a<chi2min) { f23best = f23a; chi2min = chi2a; }
      if (f23best==0 || chi2b<chi2min) { f23best = f23b; chi2min = chi2b; }
      f23prev = (chi2a<=chi2b ? f23a : f23b);
    } // for i
    
    //if (vetamin[ieta]==0) {
    //g->Draw();
    //gPad->SetLogx();
    //f23best->Draw("SAME");
    //}
    vfl2l3res[ieta] = f23best;

    TGraphErrors *gd = new TGraphErrors(g->GetN());
    gd->SetName(Form("gd_%+04dto%+04d_%s",
		     int(etamin*1000),int(etamax*1000),cr));

    double diffmax(0.), diff30max(0);
    for (int i = 0; i != g->GetN(); ++i) {
      double pt = g->GetX()[i];
      double fit = f23best->Eval(pt);
      double ptref = pt/(0.9*fit);
      double diff = fit - g->GetY()[i];
      gd->SetPoint(i, pt, 1+10.*diff);
      gd->SetPointError(i, g->GetEY()[i], 10.*g->GetEY()[i]);
      if (fabs(diff)>fabs(diffmax)) diffmax = diff;
      if (fabs(diff)>fabs(diff30max) && ptref>30) diff30max = diff;
    }
    vgdiff[ieta] = gd;

    double deta = 0.5*(etamax-etamin);
    gdmax->SetPoint(ieta, eta, min(2.95,max(-1.95,100.*diffmax)));
    gdmax->SetPointError(ieta, deta, 0.001);
    gd30max->SetPoint(ieta, eta, min(2.95,max(-1.95,100.*diff30max)));
    gd30max->SetPointError(ieta, deta, 0.001);
  } // for ieta

  
  //////////////////////
  // Setup DB mapping //
  //////////////////////

  // DB=0.90+/0.01 is generally a good approximation for <pTgen>/<pTref>,
  // but at pT<40 GeV the ratio rises rapidly to 1.1 at 15 GeV
  // Is this a real effect we should account for, or merely phase space bias?
  // Accounting for this would compress the low pT 1/x rise and make the
  // hard-to-fit part even harder to fit, so assume the latter for now
  // Proper treatment would also have to include eta-dependence also

  // Fit from jecdata*.root:/mc/eta00-13/ptchs_zjet_a100
  // TF1 *f1db = new TF1(Form("f1db_%s",cr),"[0]+[1]*log(0.01*x)+[2]*pow(0.1*x,[3])",15);
  // f1db->SetParameters(0.8968, 0.0112, 1.085, -2.562);

  
  //////////////////////////////////////////////////
  // Draw control plots of fits and fit residuals //
  //////////////////////////////////////////////////

  #include "Config.C"
  
  TH1D *h0 = tdrHist(Form("h0_%s",cr),"Re-fit max(#DeltaJES) (%)",-2,+3,
		     "#eta",-5.2,5.2);
  lumi_136TeV = mlum["24to26C"];
  extraText = "Work-in-progress";
  TCanvas *c0 = tdrCanvas("c0",h0,8,11,kSquare);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(-5.2,0,+5.2,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(-5.2,+0.5,+5.2,+0.5);
  l->DrawLine(-5.2,-0.5,+5.2,-0.5);
  
  tdrDraw(gdmax,"Pz",kOpenSquare,kBlue,kSolid,-1,kNone,0);
  tdrDraw(gd30max,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);

  TLegend *leg = tdrLeg(0.55,0.85-0.05*2,0.80,0.85);
  leg->AddEntry(gdmax,"Any p_{T} > 15 GeV","PLE");
  leg->AddEntry(gd30max,"Only p_{T} > 30 GeV","PLE");

  gPad->RedrawAxis();
  
  c0->SaveAs(Form("pdf/createL2L3ResTextFileV2/createL2L3ResTextFileV2_c0_diffmax_%s.pdf",cr));

  
  const int nbins = vfl2l3res.size()/2;
  assert(int(vfl2l3res.size())==2*nbins);
  if (debug) cout << "Plotting 2*"<<nbins<<" eta bins" << endl;
  int ncx(7), ncy(6), nc(ncx*ncy);
  TCanvas *c1 = new TCanvas("c1","c1",ncx*300,ncy*300);
  c1->Divide(ncx,ncy,0,0);

  for (int i = 0; i != min(nc,nbins); ++i) {

    c1->cd(i+1);

    int j1 = nbins+i;
    int j2 = nbins-1-i;
    double etamin = vetamin[j1]; assert(etamin>=0);
    if (etamin!=-vetamax[j2]) cout<<etamin<<" vs "<<-vetamax[j2]<<endl<<flush;
    assert(etamin==-vetamax[j2]);
    double etamax = vetamax[j1]; assert(etamax>=0);
    if (etamax!=-vetamin[j2]) cout<<etamax<<" vs "<<-vetamax[j2]<<endl<<flush;
    assert(etamax==-vetamin[j2]);

    double ymin(0.7), ymax(1.4), xmin(10), xmax(4500);
    TH1D *h = tdrHist(Form("f23_%+04dto%+04d_%s",
			   int(1000*etamin),int(1000*etamax),cr),
		      "1/L2L3Res",ymin,ymax,"p_{T,L2} (GeV)",xmin,xmax);
    h->Draw("AXIS");

    l->SetLineStyle(kDashed);
    l->DrawLine(xmin,1,xmax,1);
    l->SetLineStyle(kDotted);
    l->DrawLine(xmin,1.05,xmax,1.05);
    l->DrawLine(xmin,0.95,xmax,0.95);
    l->DrawLine(30,ymin,30,ymax);
    
    TGraphErrors *gm = vgl2l3res[j1];
    TGraphErrors *gp = vgl2l3res[j2];

    TGraphErrors *gdm = vgdiff[j1];
    TGraphErrors *gdp = vgdiff[j2];
    
    TF1 *fm = vfl2l3res[j1];
    TF1 *fp = vfl2l3res[j2];

    tdrDraw(gm,"LPz",kFullCircle,kBlue-9,kSolid,-1,1001,kBlue-9,0.3);
    tdrDraw(gp,"LPz",kFullCircle,kRed-9,kSolid,-1,1001,kRed-9,0.3);

    fm->SetLineColor(kBlue);
    fm->Draw("SAME");
    fp->SetLineColor(kRed);
    fp->Draw("SAME");

    tdrDraw(gdm,"LPz",kNone,kCyan+1,kSolid,-1,1001,kCyan+1);
    tdrDraw(gdp,"LPz",kNone,kMagenta+1,kSolid,-1,1001,kMagenta+1);

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.045*1.5);
    tex->DrawLatex(0.50,0.80,Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    
    gPad->SetLogx();
    gPad->RedrawAxis();

    if (i==0) {
      c1->cd(nc);
      TLegend *leg = tdrLeg(0.05,0.90-0.05*2.0*7,0.55,0.90);
      leg->SetTextSize(0.045*2.0);
      leg->SetHeader(cr);
      leg->AddEntry(gm,"-|#eta| data","PLE");
      leg->AddEntry(gp,"+|#eta| data","PLE");
      leg->AddEntry(fm,"-|#eta| fit","L");
      leg->AddEntry(fp,"+|#eta| fit","L");
      leg->AddEntry(gdm,"-|#eta| 10#times(fit-data)+1","LE");
      leg->AddEntry(gdp,"+|#eta| 10#times(fit-data)+1","LE");
    }		
  }

  c1->SaveAs(Form("pdf/createL2L3ResTextFileV2/createL2L3ResTextFileV2_c1_%s.pdf",cr));
  
			    
  //////////////////////////////////////////////////
  // Draw control plots of fits and fit residuals //
  //////////////////////////////////////////////////

  cout << "Writing results to:\n"
       << "   " << sout << endl;
  
  string func = vfit[vfit.size()-1];
  ofstream fout(sout.c_str());
  string header3 = "{ 1 JetEta 1 JetPt 1./("+func+") Correction L2Relative }";

  fout << header3 << endl;
  for (int i = 0; i != int(vetamin.size()); ++i) {

    double etamin = vetamin[i];
    double etamax = vetamax[i];
    TF1 *f3 = vfl2l3res[i];

    double p[npar3];
    assert(f3->GetNpar()<=npar3);
    for (int j = 0; j != npar3; ++j) {
      if (j<f3->GetNpar()) p[j] = f3->GetParameter(j);
      else p[j] = 0;
    }
    double ptminraw = f3->GetXmin();
    double ptmaxraw = f3->GetXmax();

    string s23 = Form("  %6.3f %6.3f %2d   %4.1f %4.0f   "
		      "%6.4f  %6.3f  %7.4f   %7.4f  %7.4f  %7.4f   "
		      "%7.4f  %6.3f  %6.3f",
		      etamin, etamax, npar3 + 2,
		      ptminraw, ptmaxraw,
		      p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8]);
    fout << s23 << endl;
  }

} // createL2L3ResTextFileV2s
