// Purpose: derive JEC shape from NHF-NEF/CHF difference
#include "TFile.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

void drawLowPtPFs(string run);

void drawLowPtPF() {
  drawLowPtPFs("2024B_nib1");
  drawLowPtPFs("2024C_nib1");
  drawLowPtPFs("2024D_nib1");
  drawLowPtPFs("2024Ev1_nib1");
  drawLowPtPFs("2024Ev2_nib1");
  drawLowPtPFs("2024F_nib1");
  drawLowPtPFs("2024F_nib2");
  drawLowPtPFs("2024F_nib3");
  drawLowPtPFs("2024G_nib1");
  drawLowPtPFs("2024G_nib2");
  drawLowPtPFs("2024H_nib1");
  drawLowPtPFs("2024I_nib1");
}

void drawLowPtPFs(string run) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",cr),"READ");
  assert(f && !f->IsZombie());

  f->cd("ratio");
  gDirectory->cd("eta00-13");
  TDirectory *d = gDirectory;
  f->cd("data");
  gDirectory->cd("eta00-13");
  TDirectory *dd = gDirectory;
  curdir->cd();
  
  TGraphErrors *gc, *gn, *ge;
  gc = (TGraphErrors*)d->Get("chf_incjet_a100"); assert(gc);
  gn = (TGraphErrors*)d->Get("nhf_incjet_a100"); assert(gn);
  ge = (TGraphErrors*)d->Get("nef_incjet_a100"); assert(ge);

  TGraphErrors *gcd, *gnd, *ged;
  gcd = (TGraphErrors*)dd->Get("chf_incjet_a100"); assert(gcd);
  gnd = (TGraphErrors*)dd->Get("nhf_incjet_a100"); assert(gnd);
  ged = (TGraphErrors*)dd->Get("nef_incjet_a100"); assert(ged);

  TH1D *hjes = (TH1D*)d->Get("run3/hFit_Rjet"); assert(hjes);
  
  assert(gcd->GetN()==gn->GetN());
  assert(gnd->GetN()==gn->GetN());
  assert(ged->GetN()==gn->GetN());
  
  // [Assume that NEF is absolutely right, and scale denominator to fix it]
  // Assume CHF is right, and scale denominator to fix it
  TGraphErrors *gcn = new TGraphErrors(gn->GetN());
  TGraphErrors *gnn = new TGraphErrors(gn->GetN());
  TGraphErrors *gen = new TGraphErrors(gn->GetN());
  for (int i = 0; i != gn->GetN(); ++i) {
    double pt = ge->GetX()[i];
    double pte = ge->GetEX()[i];
    /*
    double dfe = ge->GetY()[i];
    double fed = ged->GetY()[i];
    double rjet = (fed-dfe) / fed;
    */
    /*
    double dfc = gc->GetY()[i];
    double fcd = gcd->GetY()[i];
    double rjet = (fcd-dfc) / fcd;
    */
    double rjet = hjes->GetBinContent(hjes->FindBin(pt));
    // Limit response to 1.12 as original jets maybe overshot
    //rjet = min(rjet,1.12);
    
    gcn->SetPoint(i, pt, gc->GetY()[i] + gcd->GetY()[i]*(1 - 1./rjet));
    gcn->SetPointError(i, pte, gc->GetEY()[i]);
    gnn->SetPoint(i, pt, gn->GetY()[i] + gnd->GetY()[i]*(1 - 1./rjet));
    gnn->SetPointError(i, pte, gn->GetEY()[i]);
    gen->SetPoint(i, pt, ge->GetY()[i] + ged->GetY()[i]*(1 - 1./rjet));
    gen->SetPointError(i, pte, ge->GetEY()[i]);
  }

  // Then assume lost charged hadrons are recovered in NHF
  // (except at lowest pT end => model with Rjet,MC?)
  TGraphErrors *gnn2 = new TGraphErrors(gn->GetN());
  for (int i = 0; i != gn->GetN(); ++i) {
    double pt = gn->GetX()[i];
    double pte = gn->GetEX()[i];
    gnn2->SetPoint(i, pt, gnn->GetY()[i]+min(0.,gcn->GetY()[i]));
    gnn2->SetPointError(i, pte, sqrt(pow(gnn->GetEY()[i],2) +
				     pow(gcn->GetEY()[i],2)));
  }
    
  #include "../Config.C"
  double eps = 1e-4;
  TH1D *h = tdrHist("h","PF fraction (Data-MC)",-0.08+eps,0.18-eps);
  if (run=="2024B_nib1")
    lumi_136TeV = Form("%s + %s, %s","Winter24",cr,mlum[run].c_str());
  else
    lumi_136TeV = Form("%s + %s, %s","Summer24",cr,mlum[run].c_str());
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,0,3500,0);

  tdrDraw(gnn,"LPz",kFullDiamond,kGreen+2-9,kSolid,-1);
  
  tdrDraw(ge,"LPz",kOpenSquare,kBlue,kSolid,-1);
  tdrDraw(gc,"LPz",kOpenCircle,kRed,kSolid,-1);
  tdrDraw(gn,"LPz",kOpenDiamond,kGreen+2,kSolid,-1);

  tdrDraw(gen,"LPz",kFullSquare,kBlue,kSolid,-1);
  tdrDraw(gcn,"LPz",kFullCircle,kRed,kSolid,-1);
  //tdrDraw(gnn,"LPz",kFullDiamond,kGreen+2-9,kSolid,-1);
  tdrDraw(gnn2,"LPz",kFullDiamond,kGreen+2,kSolid,-1);

  //TF1 *f1 = new TF1("f1","[0]*pow(x/[4],[1])/(1+[2]*pow(x/[4],[3]))",15,3500);
  //f1->SetParameters(0.18,2,0.5,4,15.);
  //TF1 *f1 = new TF1("f1","[0]*pow(x/15.,[1])/(1+[2]*pow(x/15.,2*[1]))",15,3500);
  //f1->SetParameters(0.15,1,1);
  //f1->FixParameter(1,2);
  //f1->FixParameter(2,0.5);

  //TF1 *f1 = new TF1("f1","[0]*pow(x/15.,[1])/(1+0.5*pow(x/15.,2*[1]))",15,3500);
  TF1 *f1 = new TF1("f1","[0]*pow(x/[2],[1])/(1+0.5*pow(x/[2],2*[1]))",15,3500);
  f1->SetParameters(0.18,2,15);
  gnn2->Fit(f1,"RNW");
  f1->SetLineColor(kGreen+2);
  f1->SetLineWidth(3);
  f1->Draw("SAME");

  TLegend *leg = tdrLeg(0.55,0.90-0.04*10,0.80,0.90);
  leg->SetTextSize(0.040);
  leg->AddEntry(gc,"CHF raw","PLE");
  leg->AddEntry(gn,"NHF raw","PLE");
  leg->AddEntry(ge,"NEF raw","PLE");
  leg->AddEntry(gcn,"CHF corr","PLE");
  leg->AddEntry(gnn,"NHF corr","PLE");
  leg->AddEntry(gen,"NEF corr","PLE");
  leg->AddEntry(gnn2,"NHF - |#DeltaCHF|","PLE");
  leg->AddEntry(f1,"Fit","L");

  
  TString t(f1->GetExpFormula());
  t.ReplaceAll("[p0]",Form("%1.2f",100.*f1->GetParameter(0)));
  t.ReplaceAll("2*[p1]",Form("%1.3f",2.*f1->GetParameter(1)));
  t.ReplaceAll("[p1]",Form("%1.3f",f1->GetParameter(1)));
  t.ReplaceAll("[p2]",Form("%1.1f",f1->GetParameter(2)));

  cout << "Function to copy-paste to globalFitSettings.h" << endl;
  cout << Form("    {\"off_nhf\",\"Rjet\",\"%s\"}, =>",
	       f1->GetExpFormula().Data()) << endl;
  cout << Form("    {\"off_nhf\",\"Rjet\",\"%s\"},",
	       t.Data()) << endl;
  
  //c1->SaveAs("pdf/drawLowPtPF.pdf");
  c1->SaveAs(Form("pdf/drawLowPtPF/drawLowPtPF_%s.pdf",cr));
  
} // drawLowPtPF
