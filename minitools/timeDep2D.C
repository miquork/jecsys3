// L2L3Res time dependence with inclusive jets
// Use Incjet/h2pteta_lumi
// Projection sanity check with Incjet/p2pteta vs hpt13
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"

#include <fstream>
#include <cstdio>

#include "../tdrstyle_mod22.C"

string triggerString(double pt, double eta) {

  // Quick printout from files:
  // for (int i = 1; i != hpt13->GetNbinsX()+1; ++i) cout << hpt13->GetBinLowEdge(i) << (hpt13->GetBinContent(i)>hpt13->GetBinContent(i-1) ? "*":"") << ", "; cout << endl;
  
  // Inclusive jet thresholds, mi[]
  if (pt>=686. && fabs(eta)>2.964) return "HLT_PFJetFwd500";
  if (pt>=548. && pt<686 && fabs(eta)>2.964) return "HLT_PFJetFwd450";
  if (pt>=468. && pt<548 && fabs(eta)>2.964) return "HLT_PFJetFwd400";
  if (pt>=395. && pt<468 && fabs(eta)>2.964) return "HLT_PFJetFwd320";
  if (pt>=330. && pt<395 && fabs(eta)>2.964) return "HLT_PFJetFwd260";
  if (pt>=272. && pt<330 && fabs(eta)>2.964) return "HLT_PFJetFwd200";
  if (pt>=196. && pt<272 && fabs(eta)>3.139) return "HLT_PFJetFwd140"; // x
  if (pt>=114. && pt<196 && fabs(eta)>3.139) return "HLT_PFJetFwd80"; // x
  if (pt>=84.  && pt<114 && fabs(eta)>3.239) return "HLT_PFJetFwd60"; // x
  if (pt>=49.  && pt<84  && fabs(eta)>2.965)  return "HLT_PFJetFwd40";
    
  //if (pt>686.) return "HLT_PFJet550";
  if (pt>=686.) return "HLT_PFJet500";
  if (pt>=548.) return "HLT_PFJet450";
  if (pt>=468.) return "HLT_PFJet400";
  if (pt>=395.) return "HLT_PFJet320";
  if (pt>=330.) return "HLT_PFJet260";
  if (pt>=272.) return "HLT_PFJet200";
  if (pt>=196.) return "HLT_PFJet140";
  //if (pt>196.) return "HLT_PFJet110";
  if (pt>=114.) return "HLT_PFJet80";
  if (pt>=84.)  return "HLT_PFJet60";
  if (pt>=49.)  return "HLT_PFJet40"; // 15 now?
  return "HLT_PFJet40"; // "HLT_ZeroBias"

}

map<string, map<string, double> > _trgLumCache;
double triggerLumi(string run, string trg) {
  if (_trgLumCache[run][trg]!=0) return _trgLumCache[run][trg];

  const char *cr = run.c_str();
  //string sf = Form("rootfiles/Prompt2024/lumiText/lumi_HLT_CombinedJSONAug01_%s_PFJet.csv",cr);
  string sf = Form("rootfiles/Prompt2024/lumiText/lumi_HLT_CombinedJSONSep15_%s_PFJet.csv",cr);
  cout << "Opening " << sf << endl << flush;
  ifstream f(sf.c_str());
  assert(f.is_open());

  double lumi(0);
  float rec(0), del(0); // in /fb
  int ver,nfill,nrun,ncms,nruntot(0);
  string s;
  string sRegExp = Form("#%s_v%%d,%%d,%%d,%%d,%%f,%%f",trg.c_str());
  const char *csr = sRegExp.c_str();
  cout << "Reading regExp " << csr << endl << flush;
  while (f >> s) {
    const char *cs = s.c_str();
    if (sscanf(cs,csr,&ver,&nfill,&nrun,&ncms,&del,&rec)==6) {
      cout << "Parsing string " << cs << endl << flush;
      cout << trg << "_v" << ver<<": " << rec << endl << flush;
      lumi += rec;
      nruntot += nrun;
    }
  } // while f

  cout << trg << ": " << lumi << " ("<<nruntot<<" runs)" << endl << flush;
  if (lumi==0) lumi = -1;

  // Ad-hoc factors to match eras (brilcalc by hand over .csv)
  if (run=="2024BCD") lumi *= 15.3/14.7;
  if (run=="2024E") lumi *= 11.3/11.3;
  if (run=="2024F") lumi *= 27.8/27.8;
  if (run=="2024G") lumi *= 37.8/33.1;// * 0.97;
  
  _trgLumCache[run][trg] = lumi * 1e3; // fb-1 to pb-1

  return lumi;
}

void normalizeHistByLumi(TH1D *h, string run) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double eta = 0.;
    double pt = h->GetBinCenter(i);
    string trg = triggerString(pt, eta);
    double lumi = triggerLumi(run, trg);
    h->SetBinContent(i, h->GetBinContent(i)/lumi);
    h->SetBinError(i, h->GetBinError(i)/lumi);
  }
} // normalizeHistByLumi (TH1D)

void normalizeHistByLumiAndBin(TH2D *h2, string run) {
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      double eta = h2->GetXaxis()->GetBinCenter(i);
      double deta = h2->GetXaxis()->GetBinWidth(i);
      double pt = h2->GetYaxis()->GetBinCenter(j);
      double dpt = h2->GetYaxis()->GetBinWidth(j);
      string trg = triggerString(pt, eta);
      double lumi = triggerLumi(run, trg);
      h2->SetBinContent(i, j, h2->GetBinContent(i,j)/lumi/deta/dpt);
      h2->SetBinError(i, j, h2->GetBinError(i,j)/lumi/deta/dpt);
    }
  }
} // normalizeHistByLumi (TH2D)

TH1D* timeDep1Ds(string run="2024BCD") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TFile *f = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_data_cmb_%s_JME_v109_2024.root",cr),"READ");
  assert(f && !f->IsZombie());

  TH1D *hpt13 = (TH1D*)f->Get("Incjet/hpt13"); assert(hpt13);
  TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta"); assert(h2pteta);
  TH2D *h2pteta_lum = (TH2D*)f->Get("Incjet/h2pteta_lumi"); assert(h2pteta_lum);

  hpt13 = (TH1D*)hpt13->Clone(Form("hpt13_%s",cr));
  hpt13->Scale(1.,"width");

  int i1 = h2pteta->GetXaxis()->FindBin(-1.305+1e-2);
  int i2 = h2pteta->GetXaxis()->FindBin(+1.305-1e-2);
  TH1D *hpt13p = h2pteta->ProjectionY(Form("hpt13p_%s",cr),i1,i2);
  hpt13p->Scale(1.,"width");

  TH1D *hpt13l = h2pteta_lum->ProjectionY(Form("hpt13l_%s",cr),i1,i2);
  double lum500 = triggerLumi(run,"HLT_PFJet500");
  hpt13l->Scale(22.5/lum500*1e-5,"width"); // ??
  
  TH1D *h = tdrHist(Form("h_%s",cr),"Cross section (pb / GeV)",0.2e-7,2e11,
		    "p_{T,jet} (GeV)",15,4500.);
  TH1D *hd = tdrHist(Form("hd_%s",cr),"Ratio",0.95,1.05,
		     "p_{T,jet} (GeV)",15,4500.);
  lumi_136TeV = Form("%s, %1.1f fb^{-1}",cr,lum500);
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cr),h,hd,8,11);

  normalizeHistByLumi(hpt13,cr);
  normalizeHistByLumi(hpt13p,cr);
  
  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(hpt13,"Pz",kOpenCircle);
  tdrDraw(hpt13p,"Pz",kFullDiamond);
  //tdrDraw(hpt13l,"Pz",kOpenDiamond);

  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();
  TH1D *hpt13r = (TH1D*)hpt13p->Clone(Form("hpt13r_%s",cr));
  hpt13r->Divide(hpt13);
  tdrDraw(hpt13r,"Pz",kFullDiamond);

  TH1D *hpt13lr = (TH1D*)hpt13l->Clone(Form("hpt13lr_%s",cr));
  hpt13lr->Divide(hpt13p);
  hpt13lr->Scale(1.04);//0.22);
  //tdrDraw(hpt13lr,"Pz",kOpenDiamond);
  
  gPad->RedrawAxis();

  c1->SaveAs(Form("pdf/timeDep2D/timeDep2D_1DBB_%s.pdf",cr));
  c1->Close();

  return hpt13p;
} // TH1D* timeDep1Ds

void timeDep2Ds(string run="2024BCD", string ref="2024F") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TFile *f = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_data_cmb_%s_JME_v109_2024.root",cr),"READ");
  assert(f && !f->IsZombie());

  const char *cref = ref.c_str();
  TFile *fref = new TFile(Form("rootfiles/Prompt2024/v109_2024/jmenano_data_cmb_%s_JME_v109_2024.root",cref),"READ");
  assert(fref && !fref->IsZombie());

  TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta"); assert(h2pteta);
  TH2D *h2pteta_ref = (TH2D*)fref->Get("Incjet/h2pteta"); assert(h2pteta_ref);
  curdir->cd();
  
  double lum500 = triggerLumi(run,"HLT_PFJet500") * 1e-3;
  TH1D *h = tdrHist(Form("h_%s",cr),"p_{T,jet} (GeV)",15,4500,
		    "#eta_{jet} (rad)",-5.2,5.2);
  lumi_136TeV = Form("%s, %1.1f fb^{-1}",cr,lum500);
  extraText = "Private";
  TCanvas *c2 = tdrCanvas(Form("c2_%s",cr),h,8,11,kSquare);
  gPad->SetRightMargin(0.20);
  
  TH2D *h2 = (TH2D*)h2pteta->Clone(Form("h2pteta_%s",cr));
  normalizeHistByLumiAndBin(h2,cr);

  TH2D *h2ref = (TH2D*)h2pteta_ref->Clone(Form("h2pteta_%s_%s",cr,cref));
  normalizeHistByLumiAndBin(h2ref,cref);

  gPad->SetLogy();
  gPad->SetLogz();
  h2->Draw("SAME COLZ");
  h2->GetZaxis()->SetRangeUser(0.2e-6,2.5e8);
  h2->GetZaxis()->SetTitle("Cross section (pb / GeV / rad)");
  h2->GetZaxis()->SetTitleOffset(1.6);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetLabelSize(0.05);

  gPad->RedrawAxis();

  gPad->Update();
  c2->SaveAs(Form("pdf/timeDep2D/timeDep2D_xsec_%s.pdf",cr));

  
  TH1D *h_3 = tdrHist(Form("h_3_%s",cr),"p_{T} (GeV)",15,4500,
		      "#eta (rad)",-5.2,5.2);
  h_3->GetYaxis()->SetMoreLogLabels();
  h_3->GetYaxis()->SetNoExponent();
  TCanvas *c3 = tdrCanvas(Form("c3_%s",cr),h_3,8,11,kSquare);
  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.20);
  h_3->GetYaxis()->SetTitleOffset(1.5);

  TH2D *h2r = (TH2D*)h2->Clone(Form("h2r_%s",cr));
  h2r->Divide(h2ref);

  gPad->SetLogy();
  h2r->Draw("SAME COLZ");
  double eps = 1e-4;
  //h2r->GetZaxis()->SetRangeUser(0.5+eps,1.5-eps);
  h2r->GetZaxis()->SetRangeUser(0.75+eps,1.25-eps);
  h2r->GetZaxis()->SetTitle(Form("%s / %s",cr,cref));
  h2r->GetZaxis()->SetTitleOffset(1.6);
  h2r->GetZaxis()->SetTitleSize(0.05);
  h2r->GetZaxis()->SetLabelSize(0.05);

  gPad->RedrawAxis();

  gPad->Update();
  c3->SaveAs(Form("pdf/timeDep2D/timeDep2D_ratio_%s.pdf",cr));

  c2->Close();
  c3->Close();

} // timeDep2Ds

void timeDep2D() {

  ///////////////////////////////////////
  // Start with 1D spectra as a warmup //
  ///////////////////////////////////////
  
  // Unprescaled luminosities calculated by hand
  double lum24bcd = 15.3; // fb-1
  double lum24e = 11.3; // fb-1
  double lum24f = 27.8; // fb-1
  double lum24g = 37.8; // fb-1
  double lumsum = lum24bcd + lum24e + lum24f + lum24g;

  // Normalized spectra for each eta
  TH1D *h24bcd = timeDep1Ds("2024BCD");
  TH1D *h24e   = timeDep1Ds("2024E");
  TH1D *h24f   = timeDep1Ds("2024F");
  TH1D *h24g   = timeDep1Ds("2024G");

  TH1D *h24 = (TH1D*)h24bcd->Clone("h24"); h24->Reset();
  h24->Add(h24bcd, lum24bcd);
  h24->Add(h24e,   lum24e);
  h24->Add(h24f,   lum24f);
  h24->Add(h24g,   lum24g);
  h24->Scale(1./lumsum);

  //TH1D *h = tdrHist("h","Cross section (fb / TeV)",0.2,2e17,
  TH1D *h = tdrHist("h","Cross section (pb / GeV)",0.2e-7,2e11,
		    "p_{T,jet} (GeV)",15,4500.);
  TH1D *hd = tdrHist("hd","Ratio to 2024",0.75,1.25,
		     "p_{T,jet} (GeV)",15,4500.);
  extraText = "Private";
  lumi_136TeV = Form("2024, %1.1f fb^{-1}",lumsum);
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(h24,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(h24bcd,"Pz",kFullDiamond,kRed);
  tdrDraw(h24e,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(h24f,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(h24g,"Pz",kOpenDiamond,kBlack);

  TLegend *leg = tdrLeg(0.60,0.90-5*0.05,0.85,0.90);
  leg->AddEntry(h24,Form("2024, %1.1f fb^{-1}",lumsum),"L");
  leg->AddEntry(h24bcd,Form("BCD, %1.1f fb^{-1}",lum24bcd),"PLE");
  leg->AddEntry(h24e,Form("E, %1.1f fb^{-1}",lum24e),"PLE");
  leg->AddEntry(h24f,Form("F, %1.1f fb^{-1}",lum24f),"PLE");
  leg->AddEntry(h24g,Form("G, %1.1f fb^{-1}",lum24g),"PLE");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.35,0.85,"|#eta| < 1.305");
  
  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();

  TH1D *hr24bcd = (TH1D*)h24bcd->Clone("hr24bcd");
  TH1D *hr24e = (TH1D*)h24e->Clone("hr24e");
  TH1D *hr24f = (TH1D*)h24f->Clone("hr24f");
  TH1D *hr24g = (TH1D*)h24g->Clone("hr24g");
  TH1D *hr24 = (TH1D*)h24->Clone("hr24");
  hr24bcd->Divide(h24);
  hr24e->Divide(h24);
  hr24f->Divide(h24);
  hr24g->Divide(h24);
  hr24->Divide(h24);

  tdrDraw(hr24,"HIST",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hr24bcd,"Pz",kFullDiamond,kRed);
  tdrDraw(hr24e,"Pz",kOpenDiamond,kOrange+2);
  tdrDraw(hr24f,"Pz",kFullDiamond,kGreen+2);
  tdrDraw(hr24g,"Pz",kOpenDiamond,kBlack);
    
  gPad->RedrawAxis();

  c1->SaveAs("pdf/timeDep2D/timeDep2D_1D.pdf");


  ///////////////////
  // On to full 2D //
  ///////////////////
  timeDep2Ds("2024BCD");
  timeDep2Ds("2024E");
  timeDep2Ds("2024F");
  timeDep2Ds("2024G");
  
}
