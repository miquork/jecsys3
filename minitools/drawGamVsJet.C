// Purpose: Draw photon scale vs jet in barrel (ultimately calibrated by Zmm)
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "../tdrstyle_mod22.C"

void drawGamVsJet() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf/drawGamVsJet");
  gROOT->ProcessLine(".! touch pdf/drawGamVsJet");
  
  TFile *f = new TFile("rootfiles/Prompt2024/w43/GamHistosFill_data_2024_w43.root","READ");
  //TFile *f = new TFile("rootfiles/Prompt2024/w43/GamHistosFill_data_2024_w43.root","READ");
  assert(f && !f->IsZombie());

  TProfile *pg50 = (TProfile*)f->Get("control/pr50mpfvseta"); assert(pg50);
  TProfile *pg110 = (TProfile*)f->Get("control/pr110mpfvseta"); assert(pg110);
  TProfile *pg200 = (TProfile*)f->Get("control/pr200mpfvseta"); assert(pg200);

  TH1D *hr110 = pg110->ProjectionX("hr110");
  hr110->Divide(pg50);
  TH1D *hr200 = pg200->ProjectionX("hr200");
  hr200->Divide(pg50);
  double k200 = 1;//0.992;
  for (int i = 1; i != hr110->GetNbinsX()+1; ++i) {
    hr110->SetBinContent(i, (hr110->GetBinContent(i)-1)*100);
    hr110->SetBinError(i, hr110->GetBinError(i)*100);
    hr200->SetBinContent(i, (hr200->GetBinContent(i)*k200-1)*100);
    hr200->SetBinError(i, hr200->GetBinError(i)*k200*100*k200);
  }
  
  extraText = "Private";
  lumi_136TeV = "2024, 110 fb^{-1}";
  TH1D *h = tdrHist("h","#LTp_{T,jet}#GT / p_{T,#gamma} (MPF)^{ }",0.975,1.040,
		    "#eta_{#gamma}",-1.4795,+1.4795);
  TH1D *hd = tdrHist("hd","X/50-1 (%)",-1,+2.5,"#eta_{#gamma}",-1.4795,+1.4795);
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(-1.4795,1,+1.4795,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(-1.4795,1.005,+1.4795,1.005);
  l->DrawLine(-1.4795,0.995,+1.4795,0.995);

  tdrDraw(pg200,"Pz",kNone,kRed);
  tdrDraw(pg110,"Pz",kNone,kBlue);
  tdrDraw(pg50,"Pz",kNone,kGreen+2);

  TLegend *leg = tdrLeg(0.60,0.89-3*0.06,0.85,0.89);
  leg->AddEntry(pg200,"Photon200","PLE");
  leg->AddEntry(pg110,"Photon110EB","PLE");
  leg->AddEntry(pg50,"Photon50EB","PLE");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.35,0.85,"#gamma + jet");
  tex->DrawLatex(0.35,0.79,"|#eta_{jet}| < 1.3");
  
  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->DrawLine(-1.4795,0,+1.4795,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(-1.4795,+0.5,+1.4795,+0.5);
  l->DrawLine(-1.4795,-0.5,+1.4795,-0.5);

  tdrDraw(hr200,"Pz",kNone,kRed);
  tdrDraw(hr110,"Pz",kNone,kBlue);

  c1->SaveAs("pdf/drawGamVsJet/drawGamVsJet_eta.pdf");


  TProfile *pgpt = (TProfile*)f->Get("resp_MPFchs_DATA_a100_eta00_13");
  assert(pgpt);

  TFile *fz = new TFile("rootfiles/Prompt2024/v93/jme_bplusZ_2024_Zmm_v93.root","READ");
  assert(fz && !fz->IsZombie());
  TProfile2D *p2z = (TProfile2D*)fz->Get("data/l2res/p2m0"); assert(p2z);
  int i13 = p2z->GetXaxis()->FindBin(1.305-0.06);
  TProfile *pzpt = p2z->ProfileY("pzpt",1,i13);

  TH1D *hr = pgpt->ProjectionX("hr");
  //hr->Divide(pzpt);
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    //hr->SetBinContent(i, (hr->GetBinContent(i)-1)*100);
    //hr->SetBinError(i, hr->GetBinError(i)*100);
    int j = pzpt->GetXaxis()->FindBin(hr->GetBinCenter(i));
    if (pzpt->GetBinContent(j)!=0) {
      hr->SetBinContent(i, (pgpt->GetBinContent(i)/pzpt->GetBinContent(j)-1)*100);
      hr->SetBinError(i, sqrt(pow(pgpt->GetBinError(i),2)+pow(pzpt->GetBinError(j),2))*100);
    }
  }

  
  TH1D *h_2 = tdrHist("h_2","#LTp_{T,jet}#GT / p_{T,Z/#gamma} (MPF)^{ }",
		      0.975, 1.025+1e-3, "p_{T,Z/#gamma}", 40, 300);
  TH1D *h_2d = tdrHist("h_2d","(Z/#gamma)-1 (%)", -1, +2.5,
		       "p_{T,Z/#gamma}", 40, 300);
  TCanvas *c2 = tdrDiCanvas("c2",h_2,h_2d,8,11);

  c2->cd(1);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->DrawLine(40,1,300,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(40,1.005,300,1.005);
  l->DrawLine(40,0.995,300,0.995);
  
  tdrDraw(pzpt,"Pz",kNone,kRed);
  tdrDraw(pgpt,"Pz",kNone,kBlue);
  
  TLegend *leg2 = tdrLeg(0.60,0.89-2*0.06,0.85,0.89);
  leg2->AddEntry(pgpt,"#gamma + jet","PLE");
  leg2->AddEntry(pzpt,"Z(#mu#mu) + jet","PLE");

  tex->DrawLatex(0.35,0.84,"|#eta_{jet}| < 1.3");
  
  c2->cd(2);
  gPad->SetLogx();
  
  l->SetLineStyle(kDashed);
  l->DrawLine(40,0,300,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(40,+0.5,300,+0.5);
  l->DrawLine(40,-0.5,300,-0.5);
  
  tdrDraw(hr,"Pz",kNone,kBlue);
  
  c2->SaveAs("pdf/drawGamVsJet/drawGamVsJet_pt.pdf");
} // drawGamVsJet
