// Purpose: Run 3 evolution of minitools/fullSimShapes.C for global fit shapes
//          Inputs use jet-by-jet matching for improved stability and precision
//          PF variations are expressed as deltaE vs pTgen, stabilizing them
//          All fits use uniform logpol6 shape for smoothing
//
//          Future options:
//          - data vs data reconstruction
//          - phi wedges vs others for tracking inefficiencies
//          - CHF vs NHF vs NEF from JetVeto maps
#include "TFile.h"
#include "TProfile.h"

#include <map>
#include <string>

#include "tdrstyle_mod22.C"

bool plotLinSum = false;
bool plotQuadSum = true;//false;

void shiftAndScale(TH1D *h, double shift = 0, double scale = 1);
void addInQuadrature(TH1D *hsum, TH1D *h);
void addToErrorInQuadrature(TH1D *h, double err);


void pfShapes() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/pfShapes");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/pfShapes");

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *frun2mc = new TFile("rootfiles/pfShapes/JME-RunIISummer20UL18_2M_variations.root","READ");
  assert(frun2mc && !frun2mc->IsZombie());

  TFile *frun3mc = new TFile("rootfiles/pfShapes/JME-Run3Summer22_1M_variations-v3.root","READ");
  assert(frun3mc && !frun3mc->IsZombie());

  TFile *ftrkmc = new TFile("rootfiles/pfShapes/JME-RunIISummer20UL18_2M_variations_newtracking.root","READ");

  TFile* vf[] = {frun3mc, frun2mc, ftrkmc};
  const int nf = sizeof(vf)/sizeof(vf[0]);
  
  curdir->cd();

  // Observables
  string vobs[] = {"Rjet","chErawdivEgen","nhErawdivEgen","nemErawdivEgen"};
  const int nobs = sizeof(vobs)/sizeof(vobs[0]);

  // PF shape variations
  string vvar[] = {"hadHcalm3","HBdepth1div1p2","HBdepth2div1p2",
		   "HBdepth3div1p2","HBdepth4div1p2","HOoff",
		   "ecalm3","ecalGain12p3","ecalGain6p3","ecalGain1p3",
		   "trkEffNtrk1m3","trkEffNtrk2ToInfm3","trkEff0998Nm1",
		   "hcalPFCutsUp"};
  /*
  string vvar[] = {"hadHcalm3",
		   "ecalm3","ecalGain12p3","ecalGain6p3","ecalGain1p3",
		   "trkEffNtrk1m3","trkEffNtrk2ToInfm3","trkEff0998Nm1"};

  string vvar[] = {"trkEffNtrk1m3","trkEffNtrk2ToInfm3","trkEff0998Nm1"};
  */
  const int nvar = sizeof(vvar)/sizeof(vvar[0]);

  map<string,int> mcolor;
  mcolor["Rjet"] = kBlack;
  mcolor["chErawdivEgen"] = kRed;
  mcolor["nhErawdivEgen"] = kGreen+2;
  mcolor["nemErawdivEgen"] = kBlue;
  //
  mcolor["hadHcalm3"] = kGreen+2;
  mcolor["hadHcalp3"] = kGreen+2;
  mcolor["HBdepth1div1p2"] = kGreen;
  mcolor["HBdepth2div1p2"] = kGreen+1;
  mcolor["HBdepth3div1p2"] = kGreen+2;
  mcolor["HBdepth4div1p2"] = kGreen+3;
  mcolor["HOoff"] = kGreen+1;
  //
  mcolor["ecalm3"] = kBlue;
  mcolor["ecalGain12p3"] = kCyan;
  mcolor["ecalGain6p3"] = kCyan+1;
  mcolor["ecalGain1p3"] = kCyan+2;
  //
  mcolor["trkEffNtrk1m3"] = kRed;
  mcolor["trkEffNtrk2ToInfm3"] = kRed+1;
  mcolor["trkEff0998Nm1"] = kRed+2;
  //
  mcolor["hcalPFCutsUp"] = kMagenta+1;
  
  map<string,int> mmarker;
  mmarker["Rjet"] = kFullCircle;
  mmarker["chErawdivEgen"] = kFullDiamond;
  mmarker["nhErawdivEgen"] = kFullDiamond;
  mmarker["nemErawdivEgen"] = kFullDiamond;
  //
  mmarker["hadHcalm3"] = kFullSquare;
  mmarker["hadHcalp3"] = kFullSquare;
  mmarker["HBdepth1div1p2"] = kOpenSquare;
  mmarker["HBdepth2div1p2"] = kOpenSquare;
  mmarker["HBdepth3div1p2"] = kOpenSquare;
  mmarker["HBdepth4div1p2"] = kOpenSquare;
  mmarker["HOoff"] = kFullSquare;
  //
  mmarker["ecalm3"] = kFullCircle;
  mmarker["ecalGain12p3"] = kOpenCircle;
  mmarker["ecalGain6p3"] = kOpenCircle;
  mmarker["ecalGain1p3"] = kOpenCircle;
  //
  mmarker["trkEffNtrk1m3"] = kOpenDiamond;
  mmarker["trkEffNtrk2ToInfm3"] = kFullDiamond;
  mmarker["trkEff0998Nm1"] = kOpenDiamond;
  //
  mmarker["hcalPFCutsUp"] = kFullTriangleUp;
  
  map<string,string> mlabel;
  mlabel["Rjet"] = "Jet";
  mlabel["chErawdivEgen"] = "Charged";
  mlabel["nhErawdivEgen"] = "Hadron";
  mlabel["nemErawdivEgen"] = "EM";
  //
  mlabel["hadHcalm3"] = "HCAL -3%";
  mlabel["hadHcalp3"] = "HCAL -3%"; // modified by x-1
  mlabel["HBdepth1div1p2"] = "HB depth1 -20%"; // really 1/1.2
  mlabel["HBdepth2div1p2"] = "HB depth2 -10%"; // really 1/1.2 x 0.5
  mlabel["HBdepth3div1p2"] = "HB depth3 -20%"; // really 1/1.2
  mlabel["HBdepth4div1p2"] = "HB depth4 -30%"; // really 1/1.2 x 1.5
  mlabel["HOoff"] = "HO -100%";
  //
  mlabel["ecalm3"] = "ECAL -1%";              // modified by x1/3
  mlabel["ecalGain12p3"] = "ECAL gain12 -1%"; // modified by x-1
  mlabel["ecalGain6p3"] = "ECAL gain6 -1%";   // modified by x-1
  mlabel["ecalGain1p3"] = "ECAL gain1 -1%";   // modified by x-1
  mlabel["trkEffNtrk1m3"] = "Ntrack=1 -3%";
  mlabel["trkEffNtrk2ToInfm3"] = "Ntrack>1 -10%";  // modified by x10/3
  mlabel["trkEff0998Nm1"] = "Ntrack #times -0.2%";
  //
  mlabel["hcalPFCutsUp"] = "HCalPFCuts up";
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045);
  tex->SetNDC();

  double eps = 1e-4;
  double xmin(15), xmax(4500);
  TH1D *h2 = tdrHist("h2","Response (%)",-4+eps,4-eps,
		     //TH1D *h2 = tdrHist("h2","Response (%)",-2.5+eps,1.0-eps, // Trk only
		     "p_{T,jet} (GeV)",xmin,xmax);
  TH1D *hrjetquad(0);
  TH1D *hrjetlin(0);
  lumi_136TeV = "pfshapes.C";
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);
  c2->SetLogx();
  l->DrawLine(xmin,0,xmax,0);

  int nadd = (plotLinSum && plotQuadSum ? 2 :
	      (plotLinSum || plotQuadSum) ? 1 : 0);
  double leg2ymin = 0.55;
  double leg2y2 = max(0.90-0.035*(nvar+nadd),leg2ymin);
  TLegend *leg2 = tdrLeg(0.40,leg2y2,0.65,0.90);
  leg2->SetTextSize(leg2y2==leg2ymin ? (0.90-leg2y2)/(nvar+nadd) : 0.035);
  
  // Four-point variation plots
  for (int ivar = 0; ivar != nvar; ++ivar) {

    string sv = vvar[ivar];
    const char *cv = sv.c_str();
    TString tv(cv);
    TH1D *h = tdrHist(Form("h_%s",cv),"Response (%)",-3+eps,3-eps,
		      "p_{T,jet} (GeV)",xmin,xmax);
    if (sv=="trkEffNtrk2ToInfm3") h->GetYaxis()->SetRangeUser(-7+eps,7-eps);
    lumi_136TeV = "pfShapes.C";
    TCanvas *c1 = tdrCanvas(Form("c1_%s",cv),h,8,11,kSquare);
    gPad->SetLogx();

    l->SetLineStyle(kDashed);
    l->DrawLine(xmin,0,xmax,0);
    l->SetLineStyle(kDotted);
    l->DrawLine(xmin,+1,xmax,+1);
    l->DrawLine(xmin,-1,xmax,-1);

    TLegend *leg = tdrLeg(0.60,0.90-4*0.05,0.85,0.90);

    //tex->DrawLatex(0.20,0.75,cv);
    
    for (int iobs = 0; iobs != nobs; ++iobs) {

      string so = vobs[iobs].c_str();
      const char *co = so.c_str();
      TProfile *p(0);
      for (int i = 0; i != nf && p==0; ++i) {
	TFile *f = vf[i]; assert(f);
	p = (TProfile*)f->Get(Form("%s_%s_eta1p3",co,cv));
      }
      if (p==0) cout << Form("Graph %s %s not found!\n",co,cv) << flush;
      assert(p);
      TH1D *ho = p->ProjectionX(Form("h_%s_%s",co,cv)); 

      double k(1);
      if (tv.Contains("ecalGain")) k = -1;
      if (tv.Contains("ecal")) k *= 1./3.;
      if (tv.Contains("hadHcalp3")) k *= -1;
      if (tv.Contains("HBdepth2div1p2")) k *= 0.5;
      if (tv.Contains("HBdepth4div1p2")) k *= 1.5;
      //if (tv.Contains("trkEff") && sv!="trkEffNtrk1m3") k *= 10./3.;
      if (sv=="trkEffNtrk2ToInfm3") k *= 10./3.;

      c1->cd();

      if (iobs==0) {
	tex->SetTextSize(0.035);
	tex->DrawLatex(0.20,0.70,Form("(%s%s)",cv,k==1 ? "" : Form(" #times %1.2g",k)));
	tex->SetTextSize(0.045);
	tex->DrawLatex(0.20,0.75,mlabel[cv].c_str());
      }
      
      shiftAndScale(ho, so=="Rjet" ? -1 : 0, 100.*k);
      tdrDraw(ho,"Pz",mmarker[so],mcolor[so],kSolid,-1,kNone,0);

      leg->AddEntry(ho,mlabel[so].c_str(),"PLE");

      if (so=="Rjet") {
	c2->cd();

	// Linear addition of biases
	if (!hrjetlin) {
	  hrjetlin = (TH1D*)ho->Clone("hrjetlin");
	  if (plotLinSum)
	    tdrDraw(hrjetlin,"LE3",kNone,kYellow+2,kSolid,-1,1001,kYellow+1);//-9);
	}
	else {
	  //addLinearly(hrjetlin, ho);
	  hrjetlin->Add(ho);
	}

	// Quadratic addition of uncertainties
	if (!hrjetquad) {
	  hrjetquad = (TH1D*)ho->Clone("hrjetquad");
	  if (plotQuadSum)
	    tdrDraw(hrjetquad,"E3",kNone,kGray+1,kSolid,-1,1001,kGray);
	}
	else {
	  addInQuadrature(hrjetquad, ho);
	}
	
	// Make sure error band is not too thin to see
	//if (hrjetquad && iobs==nvar-1)
	//addToErrorInQuadrature(hrjetquad, 1.5);

	// Clone to avoid style changes in variation plots
	TH1D *ho_copy = (TH1D*)ho->Clone(Form("ho_%s",cv));
	  
	
	tdrDraw(ho_copy,"Pz",mmarker[sv],mcolor[sv],kSolid,-1,kNone,0);
	leg2->AddEntry(ho_copy,mlabel[sv].c_str(),"PLE");
      }
    } // for iobs

    gPad->RedrawAxis();
    c1->Update();
    c1->SaveAs(Form("pdf/pfShapes/pfShapes_%s.pdf",cv));
  } // for ivar

  if (plotQuadSum) leg2->AddEntry(hrjetquad,"Quadratic sum","FL");
  if (plotLinSum)  leg2->AddEntry(hrjetlin,"Linear sum","FL");
  c2->Update();
  gPad->RedrawAxis();
  c2->SaveAs("pdf/pfShapes/pfShapes_Rjet.pdf");
  
} // void pfShapes


void shiftAndScale(TH1D *h, double shift, double scale) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (h->GetBinContent(i)!=0 || h->GetBinError(i)!=0) {
      h->SetBinContent(i, (h->GetBinContent(i)+shift)*scale);
      h->SetBinError(i, h->GetBinError(i)*scale);
    } // non-empty bin
  } // for i
} // shiftAndScale


void addInQuadrature(TH1D *hsum, TH1D *h) {

  for (int i = 1; i != hsum->GetNbinsX()+1; ++i) {
    double old = hsum->GetBinContent(i);
    double val = h->GetBinContent(i);
    double old_err = hsum->GetBinError(i);
    double err = h->GetBinError(i);
    if (old*val>=0) {
      double sum = TMath::Sign(sqrt(old*old + val*val), old);
      double err_sum = sqrt(old_err*old_err + err*err);
      hsum->SetBinContent(i, sum);
      hsum->SetBinError(i, err_sum);
    }
  } // for i
} // addInQuadrature


void addToErrorInQuadrature(TH1D *h, double err) {
  
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double val = h->GetBinContent(i);
    double err_old = h->GetBinError(i);
    double err_new = sqrt(err_old*err_old+err*err);
    if (err_old*val>=0) {
      h->SetBinError(i, err_new);
    }
  } // for i
} // addToErrorInQuadrature
