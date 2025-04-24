// Purpose: Draw two PRL plots for inclusive jet xsec at 13.6 TeV (late 2024)
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

// Normalize histogram for pT and lumi
void normalizeHist(TH1D *h, double absy, string trigger="");
void unfoldHist(TH1D *h, TH1D *hu);
void cleanHist(TH1D *h, double xmin, double xmax, bool filterBins = false);

// Draw trigger efficiency in addition to data/theory ratio
const bool drawTrigEff = true;

// Fit function with some granularity vs eta
Double_t fitFunc(Double_t *xx, Double_t *p) {

  // Reference fit:
  // Form("[0]*pow(x-[3],[1])*pow(1-x*%1.4g,[2])"
  // 2.*cosh(y1)/13600.),28,6800./cosh(y1+0.125));
  double x = xx[0]; // pt
  double y1 = p[4];
  double y2 = p[5];
  assert(y2>y1);

  const int ny = 20;
  double dy = (y2-y1)/ny;
  double sumf(0);
  for (int i = 0; i != ny; ++i) {

    //double y = y1 + dy*i+0.5*dy;
    double y = y1 + dy*i;
    double erg = x*cosh(y);
    //double f = (erg<6800 ? p[0]*pow(x-p[3],p[1])*pow(1-erg/6800.,p[2]) : 0.);
    //double f = p[0]*pow(x-p[3],p[1])*pow(1-x*2.*cosh(y1)/13600.,p[2]);
    //double f = (erg<6800. ? p[0]*pow(x-p[3],p[1])*pow(1-x*2.*cosh(y)/13600.,p[2]) : 0);
    double f = (erg<6800. ? p[0]*pow(x-p[3],p[1])*pow(1-erg/6800.,p[2]) : 0);
    sumf += f * dy;
  }

  return (sumf / (y2-y1));

  //return (p[0]*pow(x-p[3],p[1])*pow(1-x*2.*cosh(y1)/13600.,p[2]));
} // fitFunc
       
void drawIncjetPRL() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //////////////////////////////////////////////
  // Step 0: get bin-by-bin unfolding from MC //
  //////////////////////////////////////////////

  TFile *fm = new TFile("rootfiles/Prompt2024/v114_2024/jmenano_mc_out_Winter24MGV14_v114_2024.root","READ");

  // Data / theory or Data / fit plots
  TH1D *h_0 = tdrHist("h_0","Gen / Reco",0,2,"Jet p_{T} (GeV)",24,5000);
  TCanvas *c0 = new TCanvas("c0","c0",5*600,2*600);
  c0->Divide(5,2,0,0);

  extraText = "Private";
  lumi_136TeV = "Winter24MGV14";
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  tex->DrawLatex(0.045,0.95,"CMS Private");
  tex->DrawLatex(0.77,0.95,lumi_136TeV + " (13.6 TeV)");

  const int nybins(10);
  int marker[nybins] =
    {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullSquare,
     kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullSquare};
  int color[nybins] =
    {kBlack, kBlack, kBlack, kBlack, kBlack,
     kRed, kBlue, kBlue, kBlue, kBlue};
  double size[nybins] =
      {0.5, 0.5, 0.5, 0.5, 0.5,
     0.5, 0.5, 0.5, 0.5, 0.5};

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);

  fm->cd("HLT_MC");
  gDirectory->cd("Incjet");
  gDirectory->cd("Unfolding");
  TDirectory *dm = gDirectory;

  map<int,TH1D*> mmc;
  map<int,TH1D*> munf;
  for (int ybin = 0; ybin != nybins; ++ybin) {

    int ieta = 5*(ybin+1);
    double y1 = 0+ybin*0.5;
    double y2 = 0.5+ybin*0.5;
    double y = 0.5*(y1+y2);

    // Original spectrum
    TH1D *hptrec = (TH1D*)dm->Get(Form("hpt%02d_rec",ieta)); assert(hptrec);
    TH1D *hptgen = (TH1D*)dm->Get(Form("hpt%02d_gen",ieta)); assert(hptgen);
    TH1D *hptr = (TH1D*)hptgen->Clone(Form("hptr_%d",ybin));
    hptr->Divide(hptrec);
    //normalizeHist(hpt, y);
    //cleanHist(hpt, 28., 6800./cosh(y1));

    TH1D *hmc = (TH1D*)hptgen->Clone(Form("hmc_%d",ybin));
    hmc->Scale(1,"width");
    mmc[ieta] = hmc;
    
    TH1D *hptr2 = (TH1D*)hptr->Clone(Form("hptr2_%d",ybin));
    unfoldHist(hptr2, hptr);
    hptr2->Divide(hptr);
    
    c0->cd(ybin+1);
    gPad->SetLogx();

    if (ybin<=4) {
      gPad->SetTopMargin(0.10);
      h_0->GetYaxis()->SetTitleSize(0.045*1.7);
      h_0->GetYaxis()->SetLabelSize(0.045*1.7);
    }
    if (ybin>=5) {
      gPad->SetBottomMargin(0.17);
      h_0->GetYaxis()->SetTitleSize(0.045*1.5);
      h_0->GetYaxis()->SetLabelSize(0.045*1.5);
      h_0->GetXaxis()->SetTitleSize(ybin>5 ? 0.045*1.7 : 0.045*1.5);
      h_0->GetXaxis()->SetLabelSize(ybin>5 ? 0.045*1.7 : 0.045*1.5);
      h_0->GetXaxis()->SetTitleOffset(ybin>5 ? 1.1 : 1.2);
    }

    h_0->DrawClone("AXIS");
    l->DrawLine(28,1.,(ybin==9 ? 200. : 5000.),1.);
    if (true) {
      tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.85 : 0.95,
		     Form("[%1.1f, %1.1f]",y1,y2));
    }
      
    //tdrDraw(hptr,"H][",kNone,color[ybin],kSolid,-1,kNone,0,1.5);
    //tdrDraw(hrf,"Pz",marker[ybin],color[ybin],kSolid,-1,kNone,0,1.5);

    hptr2->GetXaxis()->SetRangeUser(24,5000);
    hptr->GetXaxis()->SetRangeUser(24,5000);
    tdrDraw(hptr2,"H][",kNone,color[ybin],kSolid,-1,kNone,0,1.5);
    tdrDraw(hptr,"Pz",marker[ybin],color[ybin],kSolid,-1,kNone,0,1.5);
    
    //if (y1==0) leg2->AddEntry(hrf,Form("|y| < %1.1f",y2),"PLE");
    //else leg2->AddEntry(hrf,Form("%1.1f #leq |y| < %1.1f",y1,y2),"PLE");
    gPad->RedrawAxis();

    munf[ieta] = hptr;
  } // for ybin

  c0->cd();
  gPad->RedrawAxis();
  c0->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_unfolding.pdf");

  
  //////////////////////////////////////////
  // Step 1: draw (unfolded) data spectra //
  //         (vs theory or fit, or MC)    //
  //////////////////////////////////////////

  // Rerun combination with minitools/DijetHistosCombine.C if needed:
  //TFile *f = new TFile("rootfiles/Prompt2024/v114_2024/jmenano_data_cmb_2024FGHI_JME_v114_2024.root","READ");
  TFile *f = new TFile("rootfiles/Prompt2024/v116_Jet/jmenano_data_cmb_2024FGHI_JME_v116.root","READ"); // V8M+PtRaw
  assert(f && !f->IsZombie());
  curdir->cd();
  
  TH1D *h = tdrHist("h","Cross section d#sigma/dp_{T}dy (pb/GeV)",
  		    //1.001e-8,0.999e9,
		    1.001e-8,0.999e10,
  		    //"Jet p_{T} (GeV)",28,5000);
		    "Jet p_{T} (GeV)",25,5000);
  extraText = "Private";
  lumi_136TeV = "Late 2024, 82.4 fb^{-1}";

  // Spectra all on the same plot
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  h->GetYaxis()->SetTitleOffset(1.19);
  c1->SetLogx();
  c1->SetLogy();

  
  //////////////////////////////////////////
  // Step 2: draw (unfolded) data spectra //
  //         over theory (or fit, or MC)  //
  //////////////////////////////////////////
  
  // Data / theory or Data / fit plots
  TH1D *h2 = tdrHist("h2","Data / fit",0,2,"Jet p_{T} (GeV)",24,5000);
  TCanvas *c2 = new TCanvas("c2","c2",5*600,2*600);
  c2->Divide(5,2,0,0);

  //TLatex *tex = new TLatex();
  //tex->SetNDC(); tex->SetTextSize(0.045);

  tex->DrawLatex(0.045,0.95,"CMS Private");
  tex->DrawLatex(0.77,0.95,lumi_136TeV + " (13.6 TeV)");

  /* const int nybins(10); */
  /* int marker[nybins] = */
  /*   {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullSquare, */
  /*    kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullSquare}; */
  /* int color[nybins] = */
  /*   {kBlack, kBlack, kBlack, kBlack, kBlack, */
  /*    kRed, kBlue, kBlue, kBlue, kBlue}; */
  /* double size[nybins] = */
  /*     {0.5, 0.5, 0.5, 0.5, 0.5, */
  /*    0.5, 0.5, 0.5, 0.5, 0.5}; */

  //TLine *l = new TLine();
  //l->SetLineStyle(kDashed);
  //l->SetLineColor(kGray+2);
  
  c1->cd();
  TLegend *leg = tdrLeg(0.85-0.25,0.90-0.035*nybins,0.85,0.90);
  leg->SetTextSize(0.040);
  c2->cd(10);
  TLegend *leg2 = tdrLeg(0.85-0.25*1.5,0.90-0.035*nybins*1.5,0.85,0.90);
  leg2->SetTextSize(0.040*1.5);

  vector<TF1*> vf1(nybins);
  for (int ybin = 0; ybin != nybins; ++ybin) {

    int ieta = 5*(ybin+1);
    double y1 = 0+ybin*0.5;
    double y2 = 0.5+ybin*0.5;
    double y = 0.5*(y1+y2);

    // Original spectrum
    TH1D *hpt = (TH1D*)f->Get(Form("Incjet/hpt%02d",ieta)); assert(hpt);
    hpt = (TH1D*)hpt->Clone(Form("hpt_%d",ybin));
    normalizeHist(hpt, y);
    cleanHist(hpt, 28., 6800./cosh(y1));

    // Bin-by-bin unfolding
    unfoldHist(hpt, munf[ieta]);
    
    TH1D *hptf = (TH1D*)hpt->Clone(Form("hptf_%d",ybin));
    cleanHist(hptf, 28, 13600., true);

    TH1D *hmc = mmc[ieta];
    cleanHist(hmc, 28, 13600., false);
    
    // Fits for ratio
    TF1 *f1 = new TF1(Form("f1_%d",ybin),
		      Form("[0]*pow(x-[3],[1])*pow(1-x*%1.4g,[2])",
			   2.*cosh(y1)/13600.),28,6800./cosh(y1+0.125));
    f1->SetParameters(1e11,-5,10,0.);
    if (ybin>4) {
      f1->FixParameter(3,0.);
    }
    hptf->Fit(f1,"QRN");
    hptf->Fit(f1,"QRNM");

    // More advanced fits for ratio
    TF1 *f2 = new TF1(Form("f2_%d",ybin), fitFunc, //28, 6800, 6);
		      28, 6800./cosh(y1+0.125), 6);
    f2->SetParameters(f1->GetParameter(0),f1->GetParameter(1),
		      f1->GetParameter(2), f1->GetParameter(3), y1,y2);
    f2->FixParameter(4, y1);
    f2->FixParameter(5, y2);
    if (ybin>4) {
      f2->FixParameter(3,0.);
    }
    //else
    //f2->SetParLimits(3, -4,4);
    //if (ybin==9) {
    //f2->FixParameter(1,-3.);
    //}
    hptf->Fit(f2,"QRN");
    //f2->FixParameter(3,0.);
    //hptf->Fit(f2,"QRN");
    //f2->FixParameter(1,-5.);
    hptf->Fit(f2,"QRN");
    // Integral: slow, but more accurate for later
    hptf->Fit(f2,"QRNI");
    vf1[ybin] = f2;

    c1->cd();
    tdrDraw(hmc,"H][",kNone,color[ybin],kSolid,-1,kNone,0,size[ybin]);
    //hmc->SetLineWidth(2);
    tdrDraw(hpt,"H][",kNone,color[ybin],kSolid,-1,kNone,0,size[ybin]);
    tdrDraw(hptf,"Pz",marker[ybin],color[ybin],kSolid,-1,kNone,0,size[ybin]);

    if (y1==0) leg->AddEntry(hptf,Form("|y| < %1.1f",y2),"PLE");
    else leg->AddEntry(hptf,Form("%1.1f #leq |y| < %1.1f",y1,y2),"PLE");

    f1->SetLineColor(hpt->GetLineColor());
    f1->SetLineStyle(kDotted);
    f1->Draw("SAME");

    f2->SetLineColor(hpt->GetLineColor());
    //f2->SetLineWidth(2);
    f2->Draw("SAME");
    

    c2->cd(ybin+1);
    gPad->SetLogx();

    TH1D *hr = (TH1D*)hpt->Clone(Form("hr_%d",ybin));
    TH1D *hrf = (TH1D*)hptf->Clone(Form("hrf_%d",ybin));
    TH1D *hrmc = (TH1D*)hmc->Clone(Form("hrmc_%d",ybin));
    for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
      double x1 = min(hr->GetBinLowEdge(i),6800./cosh(y1));
      double x2 = min(hr->GetBinLowEdge(i+1),6800./cosh(y1));
      if (x2==x1) {
	hr->SetBinContent(i, 0.);
	hr->SetBinError(i, 0.);
	hrf->SetBinContent(i, 0.);
	hrf->SetBinError(i, 0.);
	hrmc->SetBinContent(i, 0.);
	hrmc->SetBinError(i, 0.);
	continue;
      }
      //double xsec = f1->Integral(x1, x2) / (x2 - x1);
      double xsec = f2->Integral(x1, x2) / (x2 - x1);
      if (xsec>0) {
	hr->SetBinContent(i, hpt->GetBinContent(i) / xsec);
	hr->SetBinError(i, hpt->GetBinError(i) / xsec);
	hrf->SetBinContent(i, hptf->GetBinContent(i) / xsec);
	hrf->SetBinError(i, hptf->GetBinError(i) / xsec);
	hrmc->SetBinContent(i, hmc->GetBinContent(i) / xsec);
	hrmc->SetBinError(i, hmc->GetBinError(i) / xsec);
      }
    } // for i

    if (ybin<=4) {
      gPad->SetTopMargin(0.10);
      h2->GetYaxis()->SetRangeUser(0+1e-4,2);
      h2->GetYaxis()->SetTitleSize(0.045*1.7);
      h2->GetYaxis()->SetLabelSize(0.045*1.7);
    }
    if (ybin>=5) {
      gPad->SetBottomMargin(0.17);
      h2->GetYaxis()->SetRangeUser(0,4-1e-4);
      h2->GetYaxis()->SetTitleSize(0.045*1.5);
      h2->GetYaxis()->SetLabelSize(0.045*1.5);
      h2->GetXaxis()->SetTitleSize(ybin>5 ? 0.045*1.7 : 0.045*1.5);
      h2->GetXaxis()->SetLabelSize(ybin>5 ? 0.045*1.7 : 0.045*1.5);
      h2->GetXaxis()->SetTitleOffset(ybin>5 ? 1.1 : 1.2);
    }

    h2->DrawClone("AXIS");
    l->DrawLine(28,1.,(ybin==9 ? 200. : 5000.),1.);
    if (true) {
      tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.85 : 0.95,
		     Form("[%1.1f, %1.1f]",y1,y2));
      tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.80 : 0.90,
		     Form("N_{0} = %1.3g #pm %1.3g",
			  //f1->GetParameter(0), f1->GetParError(0)));
			  f2->GetParameter(0), f2->GetParError(0)));
      tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.75 : 0.85,
		     Form("#alpha = %1.3f #pm %1.3f",
			  //f1->GetParameter(1), f1->GetParError(1)));
			  f2->GetParameter(1), f2->GetParError(1)));
      tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.70 : 0.80,
		     Form("#beta = %1.2f #pm %1.2f",
			  //f1->GetParameter(2), f1->GetParError(2)));
			  f2->GetParameter(2), f2->GetParError(2)));
      //if (f1->GetParameter(3)!=0 || f1->GetParError(3)!=0)
      if (f2->GetParameter(3)!=0 || f2->GetParError(3)!=0)
	tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.65 : 0.75,
		       Form("#Deltap_{T} = %1.2f #pm %1.2f GeV",
			    //f1->GetParameter(3), f1->GetParError(3)));
			    f2->GetParameter(3), f2->GetParError(3)));
    }

    tdrDraw(hrmc,"H][",kNone,color[ybin],kSolid,-1,kNone,0,1.5);
    tdrDraw(hr,"H][",kNone,color[ybin],kSolid,-1,kNone,0,1.5);
    tdrDraw(hrf,"Pz",marker[ybin],color[ybin],kSolid,-1,kNone,0,1.5);
    if (y1==0) leg2->AddEntry(hrf,Form("|y| < %1.1f",y2),"PLE");
    else leg2->AddEntry(hrf,Form("%1.1f #leq |y| < %1.1f",y1,y2),"PLE");
    gPad->RedrawAxis();
  } // for ybin

  c1->cd();
  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_spectra.pdf");

  c2->cd(10);
  leg2->Draw();
  c2->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_ratio.pdf");

  // Bonus: calculate probability of another pileup interaction being harder
  //        with jet at |eta|<2.5, which causes inefficiency with PUPPI
  TH1D *hmiss = (TH1D*)munf[5]->Clone("hmiss");
  for (int i = 1; i != hmiss->GetNbinsX()+1; ++i) {
    double ptmin = hmiss->GetBinLowEdge(i);
    double sumf(0);
    for (int ybin = 0; ybin != 5; ++ybin) {
      TF1 *f1 = vf1[ybin]; assert(tf1);
      sumf += f1->Integral(ptmin,6800.);
    } // for ybin
    hmiss->SetBinContent(i, 60.*sumf*1e-12 / 80e-3);
  } // for i

  c2->cd(1);
  tdrDraw(hmiss,"H][",kNone,kRed,kSolid,-1,kNone,0,1.5);
  c2->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_ratio_hmiss.pdf");

  
  ////////////////////////////////////////////////////////////
  // Step 3: Check trigger spectra vs combined spectrum     //
  //         This is basic sanity check of trigger turn-ons //
  ////////////////////////////////////////////////////////////
  
  if (drawTrigEff) {

    c2->cd(10);
  
    cout << "Do drawTrigEff" << endl << flush;
    //TFile *ft = new TFile("rootfiles/Prompt2024/v114_2024//jmenano_data_out_2024FGHI_JME_v114_2024.root","READ");
    TFile *ft = new TFile("rootfiles/Prompt2024/v116_Jet/jmenano_data_out_2024FGHI_JME_v116.root","READ"); // V8M+PtRaw
    assert(ft && !ft->IsZombie());

    string vt[] =
    //{"HLT_PFJet40"};
      {"HLT_ZeroBias", "HLT_PFJet40", "HLT_PFJet60", "HLT_PFJet80",
       "HLT_PFJet140","HLT_PFJet200", "HLT_PFJet260",
       "HLT_PFJetFwd40", "HLT_PFJetFwd60" , "HLT_PFJetFwd80",
       "HLT_PFJetFwd140", "HLT_PFJetFwd200", "HLT_PFJetFwd260"};
    const int nt = sizeof(vt)/sizeof(vt[0]);
    int color[nt] =
      {kGray+2, kMagenta+1, kBlue, kRed,
       kGreen+1, kCyan+1, kOrange+1,
       kMagenta-9, kBlue-9, kRed-9,
       kGreen-9, kRed-9, kOrange-9};

    TLegend *leg2t = tdrLeg(0.85-0.25*1.5,0.90-0.035*nt*1.5,0.85,0.90);
    leg2t->SetTextSize(0.040*1.5);
    
    for (int it = 0; it != nt; ++it) {
      const char *trg = vt[it].c_str();
      
      for (int ybin = 0; ybin != nybins; ++ybin) {

	c2->cd(ybin+1);
	int ieta = 5*(ybin+1);
	double y1 = 0+ybin*0.5;
	double y2 = 0.5+ybin*0.5;
	double y = 0.5*(y1+y2);
	    
	// Original spectrum
	TH1D *hpt = (TH1D*)f->Get(Form("Incjet/hpt%02d",ieta)); assert(hpt);
	hpt = (TH1D*)hpt->Clone(Form("hpt2_%d",ybin));
	normalizeHist(hpt, y);
	cleanHist(hpt, 28., 13600.);
	
	TH1D *ht = (TH1D*)ft->Get(Form("%s/Incjet/hpt%02d",trg,ieta));
	assert(ht);
	TH1D *hrt = (TH1D*)ht->Clone(Form("hrt_%d",ybin));
	normalizeHist(hrt, y, trg);
	cleanHist(hrt, 28., 13600.);

	hrt->Divide(hpt);

	if (it==0) {
	  h2->GetYaxis()->SetRangeUser(0.+1e-3,2-1e-3);
	  h2->DrawClone("AXIS");
	  l->DrawLine(28,1.,(ybin==9 ? 200. : 5000.),1.);
	  tex->DrawLatex(ybin%5==0 ? 0.25 : 0.10, ybin<5 ? 0.85 : 0.95,
			 Form("[%1.1f, %1.1f]",y1,y2));
	}
	if (ybin==0) {
	  leg2t->AddEntry(hrt,TString(trg).ReplaceAll("HLT_",""),"FL");
	}
	
	tdrDraw(hrt,"HIST][",kNone,color[it],kSolid,-1,kNone,0,1,1);

	if (it==nt-1) gPad->RedrawAxis();
	if (ybin==nybins-1) leg2t->Draw();
      } // for ybin
    } // for it
    
    c2->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_trigEff.pdf");

  }

} // drawIncjetPRL

#include "../minitools/DijetHistosCombine.C"
void normalizeHist(TH1D *h, double absy, string trg) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    double x = h->GetBinCenter(i);
    double ex = h->GetBinWidth(i);
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);

    /*
    mi["HLT_PFJet40"]  = range{49,  84,  0, fwdeta0};
    mi["HLT_PFJet60"]  = range{84,  114, 0, fwdeta};
    mi["HLT_PFJet80"]  = range{114, 196, 0, fwdeta};
    mi["HLT_PFJet140"] = range{196, 272, 0, fwdeta};
    mi["HLT_PFJet200"] = range{272, 330, 0, fwdeta0};
    mi["HLT_PFJet260"] = range{330, 395, 0, fwdeta0};
    mi["HLT_PFJet320"] = range{395, 468, 0, fwdeta0};
    mi["HLT_PFJet400"] = range{468, 548, 0, fwdeta0};
    mi["HLT_PFJet450"] = range{548, 686, 0, fwdeta0};
    mi["HLT_PFJet500"] = range{686,6500, 0, fwdeta0};

    double fwdeta = 3.139;  // was 2.853. 80% (100%) on negative (positive) side
    double fwdeta0 = 2.964; // 2.853; // 40 and 260 up
    mi["HLT_PFJetFwd40"]  = range{49,  84,  fwdeta0, 5.2};
    mi["HLT_PFJetFwd60"]  = range{84,  114, fwdeta, 5.2};
    mi["HLT_PFJetFwd80"]  = range{114, 196, fwdeta, 5.2};
    mi["HLT_PFJetFwd140"] = range{196, 272, fwdeta, 5.2};
    mi["HLT_PFJetFwd200"] = range{272, 330, fwdeta0, 5.2};
    mi["HLT_PFJetFwd260"] = range{330, 395, fwdeta0, 5.2};
    mi["HLT_PFJetFwd320"] = range{395, 468, fwdeta0, 5.2};
    mi["HLT_PFJetFwd400"] = range{468, 548, fwdeta0, 5.2};
    mi["HLT_PFJetFwd450"] = range{548, 686, fwdeta0, 5.2};
    mi["HLT_PFJetFwd500"] = range{686,6500, fwdeta0, 5.2};
    */

    // HLT_ZeroBias_v* fudge factor for luminosity
    //double kzb = 0.5;
    double kzb = 0.5*0.82;//*0.9;
    double kfwd80 = 0.5*1.3;//1.2;//1.1;//*0.9;
    double kfwd140 = 0.5*1.3;//1.1;
    double kfwd200 = 0.5*1.3;
    double kfwd260 = 0.5*1.3;
    double kfwd320 = 0.5*1.3;
    double kfwd400 = 0.5*1.3;

    map<string, double> ml;
    ml["HLT_PFJet500"] = 82419.755946457;
    ml["HLT_PFJet450"] = 10302.469493307;
    ml["HLT_PFJet400"] = 5151.234746654;
    ml["HLT_PFJet320"] = 2575.617373327;
    ml["HLT_PFJet260"] = 643.904343332;
    ml["HLT_PFJet200"] = 235.485016990;
    ml["HLT_PFJet140"] = 54.946503964;
    ml["HLT_PFJet80"] = 4.488592109;
    ml["HLT_PFJet60"] = 1.180534644;
    ml["HLT_PFJet40"] = 0.157404619;
    ml["HLT_ZeroBias"] =  0.157404619*kzb;

    ml["HLT_PFJetFwd500"] = 82419.755946457;
    ml["HLT_PFJetFwd450"] = 82419.755946457;
    ml["HLT_PFJetFwd400"] = 82419.755946457 * kfwd400;
    ml["HLT_PFJetFwd320"] = 20604.938986614 * kfwd320;
    ml["HLT_PFJetFwd260"] = 10302.469493307 * kfwd260;
    ml["HLT_PFJetFwd200"] = 1831.550132143 * kfwd200;
    ml["HLT_PFJetFwd140"] = 312.196045252 * kfwd140;
    ml["HLT_PFJetFwd80"] = 15.095193397 * kfwd80;
    ml["HLT_PFJetFwd60"] = 2.361069288;
    ml["HLT_PFJetFwd40"] = 0.190857498;

    // Copied from minitools/DijetHistosCombine.C
    // Need tables from there as well
    //#include "../minitools/DijetHistosCombine.C"

    double lum(0);
    if (trg=="") {

      string trigger("none");
      bool covered = false;
      for (auto& box : incjetBoxes) {
	// Check if (pt, eta) is in the region
	if (inRange(x, box.ptMin, box.ptMax) &&
	    inRange(absy, box.absyMin, box.absyMax)) {
	  
	  // Copy bin content for this trigger
	  trigger = box.triggerName;
	  
	  // Sanity check for overlapping boxes
	  if (covered) {
	    std::cerr << "Warning: Trigger " << trigger
		      << " bin (pt=" << x << ", absy=" << absy << ")"
		      << " already covered by another trigger!\n";
	  }
	  
	  // Mark as covered and break if you expect only one coverage
	  covered = true;
	  //break;
	}
      }
      
      // Sanity check for missing phase space corners
      if (!covered) {
	std::cerr << "Warning: Bin (pt=" << x
		  << ", absy=" << absy << ") not covered by any trigger!\n";
      }
      
      lum = ml[trigger];
      assert(lum!=0);
    }
    else
      lum = ml[trg];
    
    /*
    double lum(82438.188401497); // 2024FGHI_all.txt online
    // Prescale trigger luminosities calculated with brilcalc, e.g.:
    // brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json --datatag online -u /pb --begin 382229 -i Cert_Collisions2024_378981_386951_Golden.json --hltpath "HLT_PFJet500_v*" > 2024FGHI_500.txt
    if (absy<3.0) {
      if (x>686) { // HLT_PFJet500_v*
	lum = 82419.755946457;
      } else if (x>548) { // HLT_PFJet450_v*
	lum = 10302.469493307;
      } else if (x>468) { // HLT_PFJet400_v*
	lum = 5151.234746654;
      } else if (x>395) { // HLT_PFJet320_v*
	lum = 2575.617373327;
      } else if (x>330) { // HLT_PFJet260_v*
	lum = 643.904343332;
      } else if (x>272) { // HLT_PFJet200_v*
	lum = 235.485016990;
      } else if (x>196) { // HLT_PFJet140_v*
	lum = 54.946503964;
      } else if (x>114) {  // HLT_PFJet80_v*
	lum = 4.488592109;
      } else if (x>84) { // HLT_PFJet60_v*
	lum = 1.180534644;
      } else if (x>49) { // HLT_PFJet40_v*
	lum = 0.157404619;
      } else if (x>=0) { // HLT_ZeroBias_v*
	lum = 0.157432157*kzb;
      }
    }
    else {
      if (x>686) { // HLT_PFJetFwe500_v*
	lum = 82419.755946457;//82419.755946457;
      } else if (x>548) { // HLT_PFJetFwd450_v*
	lum = 82419.755946457;//10302.469493307;
      } else if (x>468) { // HLT_PFJetFwd400_v*
	lum = 82419.755946457;//5151.234746654;
      } else if (x>395) { // HLT_PFJetFwd320_v*
	lum = 20604.938986614;//2575.617373327;
      } else if (x>330) { // HLT_PFJetFwd260_v*
	lum = 10302.469493307;//643.904343332;
      } else if (x>272) { // HLT_PFJetFwd200_v*
	lum = 1831.550132143;//235.485016990;
      } else if (x>196) { // HLT_PFJetFwd140_v*
	lum = 312.196045252;//54.946503964;
      } else if (x>114) {  // HLT_PFJetFwd80_v*
	lum = 15.095193397;//4.488592109;
      } else if (x>84) { // HLT_PFJetFwd60_v*
	lum = 2.361069288;//1.180534644;
      } else if (x>49) { // HLT_PFJetFwd40_v*
	lum = 0.190857498;//0.157404619;
      } else if (x>=0) { // HLT_ZeroBias_v*
	lum = 0.157432157*kzb;
      }
    }
    */
    
    h->SetBinContent(i, y / ex / lum);
    h->SetBinError(i, ey / ex / lum);
  } // for i
  
} // normalizeHist

void unfoldHist(TH1D *h, TH1D *hu) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double x = h->GetBinCenter(i);
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);
    double r(1), er(0);
    if (x>84) {
      r = hu->GetBinContent(i);
      er = hu->GetBinError(i);
    }
    else {
      int j = hu->GetXaxis()->FindBin(84.+1);
      double r1 = hu->GetBinContent(j);
      double r2 = hu->GetBinContent(j+1);
      //r = max(0.05,r1 - (j-i)*(r2-r1));
      r = r1*pow(r1/r2,j-i);
      er = sqrt(pow(er,2)+pow(r-r1,2));
    }
    //if (y>0 && r>0) {
    //if (y>0 && r>0.1) {
    if (y>0 && r>0.25) {
    //if (y>0 && r>1./3.) {
      h->SetBinContent(i, y * r);
      h->SetBinError(i, y * r * sqrt(pow(ey/y,2)+pow(er/r,2)));
    }
    else {
      h->SetBinContent(i, 0.);
      h->SetBinError(i, 0.);
    }
  } // for i in h
} // unfoldHist

void cleanHist(TH1D *h, double xmin, double xmax, bool filterBins) {

  int ybin(0);
  if (filterBins) {
    sscanf(h->GetName(),"hptf_%d",&ybin);
    cout << "Filter ybin " << ybin << endl << flush;
  }
  
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double x = h->GetBinCenter(i);
    if (x<xmin || x>xmax) {
      h->SetBinContent(i, 0.);
      h->SetBinError(i, 0.);
    }

    if (filterBins) {
      double y = ybin*0.5;
      if ( (ybin==5 && (x<74)) || //|| (x>110 && x<200))) ||
	   //(ybin>5 && (x<84 || (x>110 && x<160) || x>0.8*6800./cosh(y))) ) {
	   //(ybin>5 && (x<28 || (x>110 && x<160) || x>0.6*6800./cosh(y))) ||
	   //(ybin<5 && (x>49 && x<57)) ||
	   //(ybin>5 && (x>49 && x<57)) ||
	   //(ybin==9 && x<49) ) {
	   //(ybin==9 && x<56) ) {
	   (ybin==9 && x<37) ) {
	h->SetBinContent(i, 0.);
	h->SetBinError(i, 0.);
      }
    }
    
  } // for i
} // cleanHist
