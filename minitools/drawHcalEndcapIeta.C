// Purpose: Draw variations of HCAL per |ieta| tower in HE (ieta=16-29)
//          Compare to gamma+jet and dijet data
//          Do fit of HE scale per ieta if feasible
#include "TFile.h"
#include "TProfile2D.h"
#include "TF2.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

// Switch off to fit only E_HCAL/E_ref component, not full response
// Helps to reduce effects from tracking and ECAL scales
bool useHcalOnly = false;//true;
double ptmin = 115;//50;//115; // 50

// Fit gamma+jet data to sum of |ieta| bin variations
map<int, TH2D*> *_mh2mc;
Double_t fit2D(Double_t *x, Double_t *p) {

  double pt = x[0];
  double eta = x[1];
  double zsum(0);//, djessum(0);
  for (int ieta = 16; ieta != 29; ++ieta) {

    TH2D *h2 = (*_mh2mc)[ieta]; assert(h2);
    int i = h2->GetXaxis()->FindBin(pt);
    int j = h2->GetYaxis()->FindBin(eta);
    int ipar = ieta-16;
    double z = h2->GetBinContent(i, j); // JES change in % for HE scale 1/1.5
    zsum += p[ipar] * z;

    // Turn parameters into scales direcly
    //double djes_dc = (0.01*z) / (1./1.5-1.);
    //double c = p[ipar]; // HE scale per ieta
    //double djes = djes_dc * (c-1.);
    //djessum += djes;
  } // for ieta

  return zsum;
  //return (1+djessum)*100.;
} // fit2D

// Zero out unnecessary regions, turn ratio to (ratio-1)*100 or to diff*100
void cleanAndScale(TH2D *h2, double x1, double x2, double y1, double y2,
		   bool useHcalOnly) {

  for (int ix = 1; ix != h2->GetNbinsX()+1; ++ix) {
    for (int iy = 1; iy != h2->GetNbinsY()+1; ++iy) {
      double x = h2->GetXaxis()->GetBinCenter(ix);
      double y = h2->GetYaxis()->GetBinCenter(iy);
      if (x>=x1 && x<x2 && y>=y1 && y<y2) {
	double z = h2->GetBinContent(ix, iy);
	double ez = h2->GetBinError(ix, iy);
	if (z!=0 && ez!=0) {
	  if (useHcalOnly)
	    h2->SetBinContent(ix, iy, z*100.);
	  else
	    h2->SetBinContent(ix, iy, (z-1)*100.);
	  h2->SetBinError(ix, iy, ez*100.);
	}
      }
      else {
	h2->SetBinContent(ix, iy, 0.);
	h2->SetBinError(ix, iy, 0.);
      }
    } // for jeta
  } // for ipt
} // cleanAndScale

void drawHcalEndcapIeta() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  string st = Form("%s_pt%1.0f",(useHcalOnly ? "HcalOnly" : "AllCalo"),ptmin);
  const char *ct = st.c_str();
  
  TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_variations-ieta16to29.root","READ");
  assert(f && !f->IsZombie());

  TFile *fg = new TFile("rootfiles/Prompt2024/GamHistosFill_data_2024G_w39.root","READ");
  assert(fg && !fg->IsZombie());
  TFile *fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_2024G-pu_w38.root","READ"); // PU reweighing => less MC stats?
  //TFile *fgm = new TFile("rootfiles/Prompt2024/GamHistosFill_mc_winter2024P8_w39.root","READ"); // No PU reweighing => more MC stats?
  assert(fg && !fg->IsZombie());
  
  TFile *fd = new TFile("rootfiles/Prompt2024/v112_2024/jmenano_data_cmb_2024G_JME_v112_2024.root","READ");
  assert(fd && !fd->IsZombie());

  TFile *fz = new TFile("rootfiles/Prompt2024/v88/jme_bplusZ_2024F_Zmm_sync_v88.root","READ");
  assert(fz && !fz->IsZombie());
  
  curdir->cd();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045*1.5);
  tex->SetNDC();
  
  double eps = 1e-4;
  //double ymin(0.80+eps), ymax(1.00);
  double ymin(0.70+eps), ymax(1.00-eps);
  double ymin2(-30+eps), ymax2(0-eps);
  //double ymin(0.50+eps), ymax(1.00);
  if (useHcalOnly) {
    ymin = -0.30+eps; ymax = 0-eps;
  }
  TCanvas *c1 = new TCanvas("c1","c1",300*5,300*3);
  c1->Divide(5,3,0,0);

  TH2D *h2ref(0);
  //double y1(1.479), y2(2.964);
  double y1(1.305), y2(3.139);
  double x1(ptmin), x2(450);
  //double x1(60), x2(450);
  map<int, TH2D*> mh2mc;
  _mh2mc = &mh2mc;
  for (int ieta = 16; ieta != 30; ++ieta) {

    TProfile2D *p2;
    if (useHcalOnly)
      p2 = (TProfile2D*)f->Get(Form("nhErawdivEgen_HEieta%ddiv1p5",ieta));
    else
      p2 = (TProfile2D*)f->Get(Form("Rjet_HEieta%ddiv1p5",ieta));
    assert(p2);
    TH2D *h2 = p2->ProjectionXY(Form("h2_%d",ieta));
    
    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
    h2->SetXTitle("p_{T} (GeV)");
    h2->SetYTitle("#eta");

    int ipad = ieta-16+1;
    c1->cd(ipad);
    gPad->SetLogx();
    
    if (ipad%5==0) gPad->SetRightMargin(0.15);
    h2->Draw((ipad%5==0) ? "COLZ" : "COL");
    h2->GetZaxis()->SetRangeUser(ymin,ymax);

    //cleanAndScale(h2,x1,x2,y1,y2);
    
    l->DrawLine(x1,y1,x2,y1);
    l->DrawLine(x1,y2,x2,y2);
    l->DrawLine(x1,y1,x1,y2);
    l->DrawLine(x2,y1,x2,y2);

    tex->DrawLatex(0.50,0.90,Form("|i#eta|=%d",ieta));

    TH2D *h2mc = (TH2D*)h2->Clone(Form("h2mc_%d",ieta));
    cleanAndScale(h2mc,x1,x2,y1,y2,useHcalOnly);
    mh2mc[ieta] = h2mc;
    
    // Reference histogram for data that has eta,pT swapped
    if (ieta==28) {
      h2ref = (TH2D*)h2->Clone("h2ref"); h2ref->Reset();
    }
  } // for ietax

  // Draw  data in last pad
  TH2D *h2data(0);
  if (true) {
    c1->cd(15);
    gPad->SetLogx();
    gPad->SetRightMargin(0.15);
    
    TProfile2D *p2(0), *p2m(0), *p2res(0), *p2jes(0);
    TH2D *h2(0), *h2m(0), *h2res(0);
    p2 = (TProfile2D*)fg->Get("Gamjet2/p2m0");
    //p2 = (TProfile2D*)fd->Get("Dijet2/p2m0tc");
    assert(p2);
    h2 = p2->ProjectionXY("h2_data");
    
    p2res = (TProfile2D*)fg->Get("Gamjet2/p2res");
    //p2res = (TProfile2D*)fd->Get("Dijet2/p2res");
    assert(p2res);
    h2res = p2res->ProjectionXY("h2res_data");
    h2->Multiply(h2res);

    p2m = (TProfile2D*)fgm->Get("Gamjet2/p2m0");
    assert(p2m);
    h2m = p2m->ProjectionXY("h2_mc");
    
    if (useHcalOnly) {
      TProfile2D *p2nhf = (TProfile2D*)fg->Get("Gamjet2/PFcomposition/p2nhf");
      assert(p2nhf);
      TH2D *h2nhf = p2nhf->ProjectionXY("h2nhf_data");
      h2->Multiply(h2nhf);
      // (MPF = E_raw/E_ref) * (NHF = Eraw,nh / Eraw) = E_raw,nh / E_ref

      TProfile2D *p2nhfm = (TProfile2D*)fgm->Get("Gamjet2/PFcomposition/p2nhf");
      assert(p2nhfm);
      TH2D *h2nhfm = p2nhfm->ProjectionXY("h2nhfm_mc");
      h2m->Multiply(h2nhfm);

      // Get denominator to Egen by factoring out MC JEC corrected Eraw
      TProfile2D *p2jes = (TProfile2D*)fz->Get("mc/l2res/p2jes");
      assert(p2jes);
      TH2D *h2jes = p2jes->ProjectionXY("h2jes_mc");
      h2->Multiply(h2jes);
      h2m->Multiply(h2jes);
      
      h2->Add(h2m,-1);
      // Caveat: only calibrated out L2L3Res, not MC truth L2L3
      // so (E_raw,L2L3 / E_raw) factor still left
      // This is large and pT-dependent in HE

    }
    else {
      h2->Divide(h2m);
    }
    
    //h2->GetXaxis()->SetRangeUser(15.,7000.);
    //h2->Draw("COLZ");

    // Swap data axes from eta,pT to pT,eta like variations
    for (int ieta = 1; ieta != h2->GetNbinsX()+1; ++ieta) {
      for (int jpt = 1; jpt != h2->GetNbinsY()+1; ++jpt) {

	double eta = h2->GetXaxis()->GetBinCenter(ieta);
	double pt = h2->GetYaxis()->GetBinCenter(jpt);
	int ipt = h2ref->GetXaxis()->FindBin(pt);
	int jeta = h2ref->GetYaxis()->FindBin(eta);
	//if (pt>=53) {
	  h2ref->SetBinContent(ipt, jeta, h2->GetBinContent(ieta, jpt));
	  h2ref->SetBinError(ipt, jeta, h2->GetBinError(ieta, jpt));
	  //}
      }
    } // for ipt
    h2ref->Draw("COLZ");
    h2ref->GetZaxis()->SetRangeUser(ymin,ymax);

    l->DrawLine(x1,y1,x2,y1);
    l->DrawLine(x1,y2,x2,y2);
    l->DrawLine(x1,y1,x1,y2);
    l->DrawLine(x2,y1,x2,y2);

    tex->DrawLatex(0.50,0.90,"#gamma+jet data");

    //cleanAndScale(h2ref,x1,x2,y1,y2);

    h2data = (TH2D*)h2ref->Clone("h2ref_data");
    cleanAndScale(h2data,x1,x2,y1,y2,useHcalOnly);
    
  } // draw data

  //if (useHcalOnly)
  //c1->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_raw_HcalOnly.pdf");
  //else
  //c1->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_raw.pdf");
  c1->SaveAs(Form("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_raw_%s.pdf",ct));


  TCanvas *c2 = new TCanvas("c2","c2",300*5,300*3);
  c2->Divide(5,3,0,0);

  for (int ieta = 17; ieta!=30; ++ieta) {
    int ipad = ieta-17+1;
    c2->cd(ipad);
    gPad->SetLogx();
    TH2D *h2 = mh2mc[ieta]; assert(h2);
    h2->GetXaxis()->SetRangeUser(x1,x2);
    h2->GetYaxis()->SetRangeUser(y1,y2);
    h2->GetZaxis()->SetRangeUser(ymin2,ymax2);
    if (ipad%5==0) gPad->SetRightMargin(0.15);
    h2->Draw(ipad%5==0 ? "COLZ" : "COL");

    double ytex(0.10);
    if (ipad>10) ytex = 0.20;
    tex->DrawLatex(0.50,ytex,Form("|i#eta|=%d",ieta));
  }

  TH2D *h2fit(0);
  TF2 *f2(0);
  if (h2data) {
    c2->cd(15);
    gPad->SetLogx();
    h2data->GetXaxis()->SetRangeUser(x1,x2);
    h2data->GetYaxis()->SetRangeUser(y1,y2);
    h2data->GetZaxis()->SetRangeUser(ymin2,ymax2);
    gPad->SetRightMargin(0.15);
    h2data->Draw("COLZ");

    tex->DrawLatex(0.50,0.20,"#gamma+jet data");

    // Fit MC to data
    c2->cd(14);
    gPad->SetLogx();
    
    f2 = new TF2("fit2D",fit2D,x1,x2,y1,y2,14);
    for (int i = 0; i != f2->GetNpar(); ++i) {
      int ieta = 16+i;
      if (ieta==28) {
	//f2->SetParameter(i,0.9);
	f2->SetParameter(i,1);//0.1);
	//f2->SetParError(i,0.1);
      }
      else {
	f2->SetParameter(i,0);//0.9);
	//f2->SetParameter(i,0.1);
	//f2->SetParError(i,0.1);
      }
      f2->SetParLimits(i,-3.,3.);
      f2->SetParName(i,Form("ieta=%d",ieta));
    }
    h2data->Fit(f2,"QRN");
    h2data->Fit(f2,"RNM");

    h2fit = (TH2D*)h2data->Clone("h2fit"); h2fit->Reset();
    for (int i = 1; i != h2fit->GetNbinsX()+1; ++i) {
      for (int j = 1; j != h2fit->GetNbinsY()+1; ++j) {
	double pt = h2fit->GetXaxis()->GetBinCenter(i);
	double eta = h2fit->GetYaxis()->GetBinCenter(j);
	h2fit->SetBinContent(i, j, f2->Eval(pt, eta));
      } // for j
    } // for i
    h2fit->Draw("COLZ");
    h2fit->GetXaxis()->SetRangeUser(x1,x2);
    h2fit->GetYaxis()->SetRangeUser(y1,y2);
    h2fit->GetZaxis()->SetRangeUser(ymin2,ymax2);

    tex->DrawLatex(0.50,0.20,"#gamma+jet fit");
  } // h2data

  //if (useHcalOnly)
  //c2->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_box_HcalOnly.pdf");
  //else
  //c2->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_box.pdf");
  c2->SaveAs(Form("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_box_%s.pdf",ct));

  TH1D *hcorr(0);
  if (h2fit && f2) {
    TH1D *h4 = tdrHist("h4","HE scale",0,2,"|i#eta|",15.5,29.5);
    lumi_136TeV = "2024G";
    extraText = "Private";
    TCanvas *c4 = tdrCanvas("c4",h4,8,11,kSquare);

    double k = max(1.,sqrt(f2->GetChisquare()/f2->GetNDF()));
    TH1D *hfit = new TH1D("hfit",";|i#eta|;HE scale",14,15.5,29.5);
    for (int ipar = 0; ipar != f2->GetNpar(); ++ipar) {
      hfit->SetBinContent(ipar+1, 1+(1/1.5-1)*f2->GetParameter(ipar));
      hfit->SetBinError(ipar+1, k*(1/1.5-1)*f2->GetParError(ipar));
      // After turning fit parameters directly to HE scales
      //hfit->SetBinContent(ipar+1, f2->GetParameter(ipar));
      //hfit->SetBinError(ipar+1, k*f2->GetParError(ipar));
    } // if ipar

    l->DrawLine(15.5,1,29.5,1);
    tdrDraw(hfit,"HP",kFullCircle,kBlack,kSolid,-1,kNone,0);

    gPad->RedrawAxis();
    //if (useHcalOnly)
    //c4->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEscale_HcalOnly.pdf");
    //else
    //c4->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEscale.pdf");
    c4->SaveAs(Form("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEscale_%s.pdf",ct));

    TH1D *h4b = tdrHist("h4b","HcalRespCorr from #gamma+jet",0,3,//8,//5,
			"|i#eta|",15.5,29.5);
    //if (useHcalOnly) h4b->GetYaxis()->SetRangeUser(0,10);
    TCanvas *c4b = tdrCanvas("c4b",h4b,8,11,kSquare);

    hcorr = (TH1D*)hfit->Clone("hcorr");
    for (int i = 1; i != hcorr->GetNbinsX()+1; ++i) {
      if (hfit->GetBinContent(i)!=0) {
	hcorr->SetBinContent(i, 1./hfit->GetBinContent(i));
	hcorr->SetBinError(i, 1./hfit->GetBinContent(i) *
			   hfit->GetBinError(i)/hfit->GetBinContent(i));
      }
    }

    l->DrawLine(15.5,1,29.5,1);
    tdrDraw(hcorr,"HP",kFullCircle,kBlack,kSolid,-1,kNone,0);

    TF1 *f1 = new TF1("f1","[0]",19.5,25.5);
    hcorr->Fit(f1,"QRNW");
    f1->DrawClone("SAME");
    f1->SetLineStyle(kDashed);
    f1->SetRange(19.5,29.5);
    f1->Draw("SAME");
    
    gPad->RedrawAxis();
    //if (useHcalOnly)
    //c4b->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEcorr_HcalOnly.pdf");
    //else
    //c4b->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEcorr.pdf");
    c4b->SaveAs(Form("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_HEcorr_%s.pdf",ct));
  }

  // Save correction for later comparisons
  if (hcorr) {
    TFile *fout = new TFile("rootfiles/drawHcalEndcapIeta.root","UPDATE");
    assert(fout && !fout->IsZombie());
    hcorr->SetName(Form("HcalRespCorr_%s",ct));
    hcorr->Write(Form("HcalRespCorr_%s",ct),TObject::kOverwrite);
    fout->Close();
  }

  // Compare corrections
  if (true) {
    
    TH1D *h5 = tdrHist("h5","HcalRespCorr from #gamma+jet",0.5,2.5,
		       "|i#eta|",15.5,29.5);
    lumi_136TeV = "2024G";
    extraText = "Private";
    TCanvas *c5 = tdrCanvas("c5",h5,8,11,kSquare);

    l->DrawLine(15.5,1,29.5,1);
    
    TFile *fin = new TFile("rootfiles/drawHcalEndcapIeta.root","READ");
    assert(fin && !fin->IsZombie());
    
    TIter next(fin->GetListOfKeys());
    while (TKey *key = (TKey*)next()) {
      // Recurse directory structure
      if (string(key->GetClassName())=="TDirectoryFile") {
	// Do nothing here
      }
      else {
	TObject *obj = key->ReadObj();
	
	if (obj->InheritsFrom("TH1D")) {
	  cout << obj->GetName() << endl << flush;
	  TH1D *hcorr = (TH1D*)obj;
	  tdrDraw(hcorr,"HP",kFullCircle,kBlack,kSolid,-1,kNone,0);

	  TF1 *f1c = new TF1("f1c","[0]",19.5,25.5);
	  hcorr->Fit(f1c,"QRNW");
	  f1c->DrawClone("SAME");
	  f1c->SetLineStyle(kDashed);
	  f1c->SetRange(19.5,29.5);
	  f1c->Draw("SAME");
	}
      } // else	
    } // while key
    c5->SaveAs("pdf/drawHcalEndcapIeta/drawHcalEndcapIeta_xcomparisons.pdf");
  } // compare corrections
  
} // Drawhcalendcapieta
