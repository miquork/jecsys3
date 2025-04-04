// Purpose: Derive MC truth JES + JER from dijet MC
//          1) Primary fit with low pT bias corrected JES+JER vs gen pT
//            1a) Arithmetic mean/RMS + low pT bias corrections
//            1b) Gaussian mean+sigma + total integral constraint (low pT bias)
//            1c) Median+CI68 + low pT bias corrections
//          2) Invert response vs reco pT for JEC
//            2a) Implement effective physical parameterization
//          3) Additional (small) JER+non-linearity correction to JEC 
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include "tdrstyle_mod22.C"

bool debug = false;
bool correctLowPtBias = true;//false;

void evalJES(TH1D *h);
void evalJER(TH1D *h, TProfile *p);
void evalJET(TH1D *h, const TProfile *p);

void fixJESandJER(TH1D *hjes, TH1D *hjer, TH1D *heff, const TF1 *f1jer);

void MCTruth() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/Prompt2024/v121_Jet/jmenano_mc_out_Summer24MG_v121_2024CDEFGHI.root","READ");
  assert(f && !f->IsZombie());

  gDirectory->cd("HLT_MC");
  gDirectory->cd("MCtruth");
  TDirectory *d = gDirectory;

  // <pTreco/pTgen> vs (eta_jet, pTgen); 41 bins in X=eta_jet
  TProfile2D *p2e = (TProfile2D*)d->Get("p2eff");   assert(p2e);
  TProfile2D *p2r = (TProfile2D*)d->Get("p2r_raw"); assert(p2r);
  TProfile2D *p2c = (TProfile2D*)d->Get("p2r");     assert(p2c);

  // 7*6=42>41 pads; last one for legend
  int size = 600;
  TCanvas *c1eff = new TCanvas("c1eff","c1eff",7*size,6*size);
  c1eff->Divide(7,6,0,0);
  TCanvas *c1jes = new TCanvas("c1jes","c1jes",7*size,6*size);
  c1jes->Divide(7,6,0,0);
  TCanvas *c1jer = new TCanvas("c1jer","c1jer",7*size,6*size);
  c1jer->Divide(7,6,0,0);
  TCanvas *c1jet = new TCanvas("c1jet","c1jet",7*size,6*size);
  c1jet->Divide(7,6,0,0);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045*1.5);

  c1eff->cd(42);
  int nentry = (correctLowPtBias ? 5 : 3);
  TLegend *legeff = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legeff->SetTextSize(0.045*1.5);
  
  c1jes->cd(42);
  TLegend *legjes = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legjes->SetTextSize(0.045*1.5);

  c1jer->cd(42);
  TLegend *legjer = tdrLeg(0.05,0.80-0.05*1.5*nentry,0.50,0.80);
  legjer->SetTextSize(0.045*1.5);

  double eps = 1e-4;
  
  for (int i = 1; i != p2r->GetNbinsX()+1; ++i) {
    //for (int j = 1; j != p2r->GetNbinsY()+1; ++j) {

    //////////////////////////////////////////////////////////////////////
    // Step 1. First round fitting low pT biased results at high pT end //
    //////////////////////////////////////////////////////////////////////
    
    // Initial EFF, JES, JER and JET counts
    TProfile *pe = p2e->ProfileY(Form("pe_%d",i),i,i);
    TH1D *he = pe->ProjectionX(Form("he_%d",i));
    TProfile *pr = p2r->ProfileY(Form("pr_%d",i),i,i);
    TH1D *hr = pr->ProjectionX(Form("hr_%d",i));
    evalJES(hr);
    TProfile *pc = p2c->ProfileY(Form("pc_%d",i),i,i);
    pc->SetErrorOption("S");
    TH1D *hs = pc->ProjectionX(Form("hs_%d",i));
    evalJER(hs, pc);
    TH1D *hn = pc->ProjectionX(Form("hn_%d",i));
    evalJET(hn, pc);

    
    ////////////////////////////////////////////
    // Step 1.1 Jet reconstruction EFFiciency //
    ////////////////////////////////////////////
    
    c1eff->cd(i);
    gPad->SetLogx();

    double effmin(0.), effmax(1.2);
    TH1D *heff = tdrHist(Form("heff_%d",i),"Efficiency",effmin+eps,effmax-eps,
		      "p_{T,gen} (GeV)",5,5000);
    heff->GetXaxis()->SetMoreLogLabels(kFALSE);
    heff->Draw("AXIS");
    l->DrawLine(5,1,5000,1);

    double etamin = p2r->GetXaxis()->GetBinLowEdge(i);
    double etamax = p2r->GetXaxis()->GetBinLowEdge(i+1);
    double ptmin_eff = 15.;
    double ptmin = 30.;
    double ptmin_fwd = 30.;
    double ptmax = min(4000.,0.7*6800./cosh(etamax));
    TF1 *f1eff = new TF1(Form("f1eff_%d",i),
			 "0.5*(1+erf((x-[0])/([1]*[0])))",
			 5.,ptmax);
    f1eff->SetParameters(15.,0.30);
    he->Fit(f1eff,"QRN");
    
    tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(he,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1eff->Draw("SAME");

    if (i==1) {
      legeff->SetHeader("Summer24, MG QCD MC");
      legeff->AddEntry(hr,"Data","PLE");
      legeff->AddEntry(f1eff,"Fit","L");
    }

    
    /////////////////////////////////////
    // Step 1.2 Jet Energy Scale (JES) //
    /////////////////////////////////////

    c1jes->cd(i);
    gPad->SetLogx();

    double jesmin(0.5), jesmax(1.5);
    // Range after low pT bias correction
    if (i<=7)       { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=14) { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=21) { jesmin = 0.75; jesmax = 1.05; }
    else if (i<=28) { jesmin = 0.50; jesmax = 1.10; }
    else if (i<=35) { jesmin = 0.35; jesmax = 1.20; }
    else if (i<=42) { jesmin = 0.35; jesmax = 1.15; }
    TH1D *hjes = tdrHist(Form("hjes_%d",i),"JES",jesmin+eps,jesmax-eps,
		      "p_{T,gen} (GeV)",5,5000);
    hjes->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjes->Draw("AXIS");
    l->DrawLine(5,1,5000,1);
    
    TF1 *f1jes = new TF1(Form("f1jes_%d",i),
			 "[0]+[1]*pow(x,[2])+[3]/x+[4]*1e-5*x",
			 ptmin,ptmax);
    // Enforce physical ranges
    f1jes->SetParameters(1,-1,-0.3,0,1);
    f1jes->SetParLimits(0,0.8,1.2);
    f1jes->SetParLimits(1,-10,0.);
    f1jes->SetParLimits(2,-0.5,-0.2);
    f1jes->SetParLimits(3,-10.,+10.);
    f1jes->SetParLimits(4,-10,+10.);
    // Switch off trk parameters outside tracker
    if (etamin>2.65) {
      f1jes->FixParameter(4, 0.);
      f1jes->SetRange(ptmin_fwd,ptmax);
    }
    hr->Fit(f1jes,"QRN");

    tex->DrawLatex(i%7==1 ? 0.30 : 0.20, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jes->Draw("SAME");

    if (i==1) {
      legjes->SetHeader("Summer24, MG QCD MC");
      legjes->AddEntry(hr,"Data","PLE");
      legjes->AddEntry(f1jes,"Fit","L");
    }

    
    //////////////////////////////////////////
    // Step 1.3 Jet Energy Resolution (JER) //
    //////////////////////////////////////////
    
    c1jer->cd(i);
    gPad->SetLogx();

    double jermin(0.), jermax(0.5);
    // Range after low pT bias correction
    if (i<=7)       { jermin = 0.; jermax = 0.50; }
    else if (i<=14) { jermin = 0.; jermax = 0.50; }
    else if (i<=21) { jermin = 0.; jermax = 0.55; }
    else if (i<=28) { jermin = 0.; jermax = 1.00; }
    else if (i<=35) { jermin = 0.; jermax = 1.10; }
    else if (i<=42) { jermin = 0.; jermax = 0.75; }
    TH1D *hjer = tdrHist(Form("hjer_%d",i),"JER",jermin+eps,jermax-eps,
			 "p_{T,gen} (GeV)",5,5000);
    hjer->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjer->Draw("AXIS");

    TF1 *f1jer = new TF1(Form("f1jer_%d",i),
			 "sqrt([0]*[0]/(x*x)+[1]*[1]/x+[2]*[2])",
			 ptmin,ptmax);
    f1jer->SetParameters(1,1,0.05);
    f1jer->SetParLimits(0,0,10.);
    f1jer->SetParLimits(1,0.5,1.5);
    f1jer->SetParLimits(2,0.03,0.12);
    if (etamin>3.139) {
      f1jer->FixParameter(2,0.08);
    }
    hs->Fit(f1jer,"QRN");

    tex->DrawLatex(i%7==1 ? 0.50 : 0.40, 0.90,
		   Form("%1.3f<|#eta|<%1.3f",etamin,etamax));
    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jer->Draw("SAME");

    if (i==1) {
      legjer->SetHeader("Summer24, MG QCD MC");
      legjer->AddEntry(hr,"Data","PLE");
      legjer->AddEntry(f1jes,"Fit","L");
    }


    /////////////////////////
    // Step 1.4 JET counts //
    /////////////////////////
    
    c1jet->cd(i);
    gPad->SetLogx();
    gPad->SetLogy();
    
    TH1D *hjet = tdrHist(Form("hjet_%d",i),
			 "Weighted N(jet)/d#eta/dp_{T} (1/GeV)",
			 1e-10,1e10,"p_{T,gen} (GeV)",5,5000);
    hjet->GetXaxis()->SetMoreLogLabels(kFALSE);
    hjet->GetYaxis()->SetMoreLogLabels(kFALSE);
    hjet->Draw("AXIS");

    TF1 *f1jet = new TF1("f1jet","[0]*pow(x,[1])*pow(1-2*x/[3],[2])",
			 ptmin,ptmax);
    f1jet->SetParameters(1e11,-5,5,13800./cosh(etamin));
    f1jet->FixParameter(3,13800./cosh(etamin));
    hn->Fit(f1jet,"QRN");
    
    tdrDraw(hn,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0);
    f1jet->Draw("SAME");


    /////////////////////////////////////////////////////////////////////////
    // Step 2. Second round with low pT bias correction and extending down //
    /////////////////////////////////////////////////////////////////////////
    if (!correctLowPtBias) continue;

    TF1 *f1eff_raw = (TF1*)f1eff->Clone(Form("f1eff_raw_%d",i));
    TH1D *he_raw = (TH1D*)he->Clone(Form("he_raw_%d",i));
    TF1 *f1jes_raw = (TF1*)f1jes->Clone(Form("f1jes_raw_%d",i));
    TH1D *hr_raw = (TH1D*)hr->Clone(Form("hr_raw_%d",i));
    TF1 *f1jer_raw = (TF1*)f1jer->Clone(Form("f1jer_raw_%d",i));
    TH1D *hs_raw = (TH1D*)hs->Clone(Form("hs_raw_%d",i));
    fixJESandJER(hr, hs, he, f1jer);


    c1eff->cd(i);

    f1eff->SetRange(15.,ptmax);
    he->Fit(f1eff,"QRN");

    tdrDraw(he_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1eff_raw->Draw("SAME");
    tdrDraw(he,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1eff->SetLineColor(kGreen+2);
    f1eff->Draw("SAME");

    if (i==1) {
      legeff->AddEntry(he_raw,"Orig. data","PLE");
      legeff->AddEntry(f1eff_raw,"Orig. fit","L");
    }

    
    c1jes->cd(i);

    f1jes->SetRange(15.,ptmax);
    hr->Fit(f1jes,"QRN");

    tdrDraw(hr_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1jes_raw->Draw("SAME");
    tdrDraw(hr,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1jes->SetLineColor(kGreen+2);
    f1jes->Draw("SAME");

    if (i==1) {
      legjes->AddEntry(hr_raw,"Orig. data","PLE");
      legjes->AddEntry(f1jes_raw,"Orig. fit","L");
    }
    

    c1jer->cd(i);

    f1jer->SetRange(15.,ptmax);
    hs->Fit(f1jer,"QRN");

    tdrDraw(hs_raw,"Pz",kOpenCircle,kRed+1,kSolid,-1,kNone,0,1.5,1);
    f1jer_raw->Draw("SAME");
    tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,2.0,1);
    f1jer->SetLineColor(kGreen+2);
    f1jer->Draw("SAME");

    if (i==1) {
      legjer->AddEntry(hr_raw,"Orig. data","PLE");
      legjer->AddEntry(f1jes_raw,"Orig. fit","L");
    }
    
    //} // for j
  } // for i


  c1eff->SaveAs("pdf/MCTruth/MCTruth_EFF.pdf");
  c1jes->SaveAs("pdf/MCTruth/MCTruth_JES.pdf");
  c1jer->SaveAs("pdf/MCTruth/MCTruth_JER.pdf");
  c1jet->SaveAs("pdf/MCTruth/MCTruth_JET.pdf");
}


// Evaluate JES
// - add minimum uncertainty to avoid overfitting high statistics
void evalJES(TH1D *h) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double jes = h->GetBinContent(i);
    if (jes!=0) {
      double errmin = 0.0020;
      double ejes = sqrt(pow(h->GetBinError(i),2)+errmin*errmin);
      h->SetBinContent(i, jes);
      h->SetBinError(i, ejes);
    }
  }
  return;
}

// Evaluate JER
// - assume input histogram projected from JES profile after JEC applied
// - put bin content to RMS extracted from profile width divided by JES~1
// - update error to error on the mean over JER
void evalJER(TH1D *h, TProfile *p) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double jer = h->GetBinError(i); // RMS with option "S"
    p->SetErrorOption("");
    double ejer = (p ? p->GetBinError(i) : 0); // RMS/sqrt(N) with option ""
    double jes = (p ? p->GetBinContent(i) : 1);
    if (jer!=0 && jes!=0) {
      h->SetBinContent(i, jer/jes);
      h->SetBinError(i, ejer/(jer/jes));
    }
  }
  return;
}

// Evaluate JET counts
// - get jet counts from profile entries per bin
// - normalize by bin width in eta and pT
void evalJET(TH1D *h, const TProfile *p) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double n = p->GetBinEntries(i);
    if (n!=0) {
      double effn = p->GetBinEffectiveEntries(i);
      double en = (effn>0 ? n/sqrt(effn) : 0);
      double dpt = h->GetBinWidth(i);
      double deta = 2*0.087;
      h->SetBinContent(i, n / dpt / deta);
      h->SetBinError(i, en / dpt / deta);
    }
  }
  return;
}

// Low pT bias correction for JES
// - fill missing part with Gaussian extension
// - efficiency as profile to get accurate value per bin
// - JER as fit to extrapolate robustly to biased region
//   (also to have eventual JES consistent with JER)
void fixJESandJER(TH1D *hjes, TH1D *hjer, TH1D *heff, const TF1 *f1jer) {

  for (int i = 1; i != hjes->GetNbinsX()+1; ++i) {

    double pt = hjes->GetBinCenter(i);
    double ptmin = hjes->GetBinLowEdge(i);
    double ptmax = hjes->GetBinLowEdge(i+1);
    double jes_raw = hjes->GetBinContent(i);
    double ejes = hjes->GetBinError(i);
    double jer_raw = hjer->GetBinContent(i);
    double ejer = hjer->GetBinError(i);
    double eff = heff->GetBinContent(i);
    double eeff = heff->GetBinError(i);
    // Cap jer at 100%
    double jer_fit = min(1.0, f1jer->Eval(pt));

    // Only run the extra correction if efficiency is less than 99%
    double jes(jes_raw), jer(jer_raw), eff_cut(eff);
    //if (eff<0.16) { // 1-sided 1-sigma cut
    //if (eff<0.5) { // cut at peak
    if (eff<0.25) { // 
      jes = ejes = jer = ejer = eff_cut = eeff = 0;
    }
    else if (eff<0.999) {
      
      // Starting points:
      // 1) we know the relative width sigma/mu of the Gaussian from JER (f1jer)
      // 2) we know what fraction of the Gaussian is lost from EFF (heff)
      // 3) we don't know unbiased mean from JES (hjes), but can start with mu=1
      // => estimate zcut that integrates to eff, 
      //    then solve for the mean and sigma of the truncated Gaussian
      double mu(1), sigma(1);
      double zcut = TMath::ErfInverse(1. - 2.*eff) * sqrt(2.);
      double phi = TMath::Gaus(zcut, 0, 1, true); // normalized
      double mu_bias = 1 + jer_fit * phi / eff;
      double sigma_bias = sqrt(1 + zcut * phi / eff - pow(phi / eff, 2))
	/ mu_bias;
      
      // Correct jes_raw based on mu/mu_raw
      jes = (mu_bias>0 ? jes_raw * mu / mu_bias : jes_raw);
      jer = (sigma_bias>0  && mu_bias>0 ?
	     (jer_raw * sigma / sigma_bias) : jer_raw);

      // Rough place-holder errors for now
      ejes = sqrt(pow(ejes,2)+pow(jes_raw*(mu/mu_bias-1)*0.5,2));
      ejer = sqrt(pow(ejer,2)+pow(jer_raw*(sigma/sigma_bias-1)*0.5,2));
      
      if (debug)
	cout << Form("pt=%1.0f [%1.0f,%1.0f] GeV, eff=%1.3f, jer=%1.3f, "
		     "zcut=%1.3f, phi=%1.3f, mu_bias=%1.3f\n"
		     "* jes_raw=%1.3f, jes=%1.3f, "
		     "* jer_raw=%1.3f, jer=%1.3f\n",
		     pt, ptmin, ptmax, eff, jer, zcut, phi, mu_bias,
		     jes_raw, jes, jer_raw, jer);
    }
    heff->SetBinContent(i, eff_cut);
    heff->SetBinError(i, eeff);
    hjes->SetBinContent(i, jes);
    hjes->SetBinError(i, ejes);
    hjer->SetBinContent(i, jer);
    hjer->SetBinError(i, ejer);
  } // for i
  return;
}
