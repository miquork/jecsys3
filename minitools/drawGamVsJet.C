// Purpose: Draw photon scale vs jet in barrel (ultimately calibrated by Zmm)
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "../tdrstyle_mod22.C"

// Scales for 2025 (adjust 2024 in code)
double k200 = 0.989;//0.992;
double k110 = 0.997;
double k50 = 1;

double etamin = -2.5;//-1.5;
double etamax = +2.5;//+1.5;
double rmax = 1.10;//1.035-1e-4;//1.05;
double rmin = 0.87;//0.955;
double rrmin = -13;//-1.6;
double rrmax = +6;//+1.2;
double retamin = -3.5;//-2.2;//-1.5;
double retamax = +6;//+3.2;//+2.0;



void drawGamVsJet(string year = "2025") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf/drawGamVsJet");
  gROOT->ProcessLine(".! touch pdf/drawGamVsJet");

  TFile *f(0);
  if (year=="2024") f = new TFile("rootfiles/Prompt2024/w43/GamHistosFill_data_2024_w43.root","READ");
  if (year=="2025") f = new TFile("rootfiles/Prompt2025/Gam_w64/GamHistosFill_data_2025CDEFG_w64.root","READ");
  assert(f && !f->IsZombie());

  TProfile *pg50 = (TProfile*)f->Get("control/pr50mpfvseta"); assert(pg50);
  TProfile *pg110 = (TProfile*)f->Get("control/pr110mpfvseta"); assert(pg110);
  TProfile *pg200 = (TProfile*)f->Get("control/pr200mpfvseta"); assert(pg200);

  curdir->cd();

  // Adjust scales for 2024
  if (year=="2024") {
    k200 = 0.984;//0.989
    k110 = 0.992;//0.997;
  }
  
  TH1D *hg50 = pg50->ProjectionX("hg50");     hg50->Scale(k50);
  TH1D *hg110 = pg110->ProjectionX("hg110");  hg110->Scale(k110);
  TH1D *hg200 = pg200->ProjectionX("hg200");  hg200->Scale(k200);
  
  //TH1D *hr110 = pg110->ProjectionX("hr110");
  TH1D *hr110 = (TH1D*)hg110->Clone("hr110");
  hr110->Divide(hg50);
  //TH1D *hr200 = pg200->ProjectionX("hr200");
  TH1D *hr200 = (TH1D*)hg200->Clone("hr200");
  hr200->Divide(hg50);

  for (int i = 1; i != hr110->GetNbinsX()+1; ++i) {
    hr110->SetBinContent(i, (hr110->GetBinContent(i)-1)*100);
    hr110->SetBinError(i, hr110->GetBinError(i)*100);
    hr200->SetBinContent(i, (hr200->GetBinContent(i)-1)*100);
    hr200->SetBinError(i, hr200->GetBinError(i)*100);
    // Patch for hr200 beyond EB
    double eta = hr110->GetXaxis()->GetBinCenter(i);
    if (fabs(eta)>1.479) {
      hg50->SetBinContent(i, 0);
      hg50->SetBinError(i, 0);
      hg110->SetBinContent(i, 0);
      hg110->SetBinError(i, 0);
      hr110->SetBinContent(i, 0);
      hr110->SetBinError(i, 0);
      int j = hg50->GetXaxis()->FindBin(eta<0 ? -1.479+1e-3 : +1.479-1e-3);
      double keta50 = hg50->GetBinContent(j);
      double keta200 = hg200->GetBinContent(j);
      double k = k200/k50;
      hr200->SetBinContent(i, (pg200->GetBinContent(i)*k/keta50-1)*100);
      hr200->SetBinError(i, pg200->GetBinError(i)*k/keta50*100);
    }
  }
  
  extraText = "Private";
  if (year=="2024") lumi_136TeV = "2024, 110 fb^{-1}";
  if (year=="2025") lumi_136TeV = "2025, ~100 fb^{-1}";
  //TH1D *h = tdrHist("h","#LTp_{T,jet}#GT / p_{T,#gamma} (MPF)^{ }",0.975,1.040,
  TH1D *h = tdrHist("h","#LTp_{T,jet}#GT / p_{T,#gamma} (MPF)^{ }", rmin, rmax,
		    "Photon #eta_{#gamma}",etamin,etamax);//-1.5,1.5);//-1.4795,+1.4795);
  //TH1D *hd = tdrHist("hd","X/50-1 (%)",-1,+2.5,"#eta_{#gamma}",-1.4795,+1.4795);
  TH1D *hd = tdrHist("hd","X/50-1 (%)",rrmin,rrmax,"Photon #eta_{#gamma}",
		     etamin,etamax);//-1.5,1.5);//-1.4795,+1.4795);
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);
  h->GetYaxis()->SetTitleOffset(1.04);

  c1->cd(1);
  
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(etamin,1,etamax,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(etamin,1.005,etamax,1.005);
  l->DrawLine(etamin,0.995,etamax,0.995);
  l->DrawLine(-0.8,rmin,-0.8,rmax-(rmax-rmin)*0.23);
  l->DrawLine(+0.8,rmin,+0.8,rmax-(rmax-rmin)*0.23);
  l->DrawLine(-1.4795,rmin,-1.4795,rmax-(rmax-rmin)*0.23);
  l->DrawLine(+1.4795,rmin,+1.4795,rmax-(rmax-rmin)*0.23);

  tdrDraw(hg50,"HIST][",kNone,kGreen+2,kSolid,-1,kNone,0);
  tdrDraw(hg200,"Pz",kNone,kRed);
  tdrDraw(hg110,"Pz",kNone,kBlue);
  tdrDraw(hg50,"Pz",kNone,kGreen+2,kSolid,-1,kNone,0);

  TLegend *leg = tdrLeg(0.50,0.89-3*0.06,0.75,0.89);
  leg->AddEntry(hg200,Form("Photon200    #times %1.3f",k200),"PLE");
  leg->AddEntry(hg110,Form("Photon110EB #times %1.3f",k110),"PLE");
  leg->AddEntry(hg50,Form("Photon50EB  #times %1.3f",k50),"PLE");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.33,0.84,"#gamma + jet");
  tex->DrawLatex(0.33,0.78,"|#eta_{jet}| < 1.3");
  gPad->RedrawAxis();
  
  c1->cd(2);

  l->SetLineStyle(kDashed);
  l->DrawLine(etamin,0,etamax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(etamin,+0.5,etamax,+0.5);
  l->DrawLine(etamin,-0.5,etamax,-0.5);
  l->DrawLine(-0.8,rrmin,-0.8,rrmax-(rrmax-rrmin)*0.);
  l->DrawLine(+0.8,rrmin,+0.8,rrmax-(rrmax-rrmin)*0.);
  l->DrawLine(-1.4795,rrmin,-1.4795,rrmax-(rrmax-rrmin)*0.);
  l->DrawLine(+1.4795,rrmin,+1.4795,rrmax-(rrmax-rrmin)*0.);

  tdrDraw(hr200,"Pz",kNone,kRed);
  tdrDraw(hr110,"Pz",kNone,kBlue);

  gPad->RedrawAxis();
  
  c1->SaveAs(Form("pdf/drawGamVsJet/drawGamVsJet_eta_%s.pdf",year.c_str()));


  TProfile *pgpt = (TProfile*)f->Get("resp_MPFchs_DATA_a100_eta00_13");
  assert(pgpt);

  TFile *fz(0);
  if (year=="2024") fz = new TFile("rootfiles/Prompt2024/v93/jme_bplusZ_2024_Zmm_v93.root","READ");
  if (year=="2025") fz = new TFile("rootfiles/Prompt2025/Zmm_v102/jme_Zj_2025_Zmm_v102_nomu.root","READ");
  assert(fz && !fz->IsZombie());
  TProfile2D *p2z = (TProfile2D*)fz->Get("data/l2res/p2m0"); assert(p2z);
  int i13 = p2z->GetXaxis()->FindBin(1.305-0.06);
  TProfile *pzpt = p2z->ProfileY("pzpt",1,i13);

  TH1D *hgpt = pgpt->ProjectionX("hgpt");
  TH1D *hzpt = pzpt->ProjectionX("hzpt");
  //TH1D *hr = pgpt->ProjectionX("hr");
  TH1D *hr = (TH1D*)hgpt->Clone("hr");
  //hr->Divide(pzpt);
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    //hr->SetBinContent(i, (hr->GetBinContent(i)-1)*100);
    //hr->SetBinError(i, hr->GetBinError(i)*100);
    int j = hzpt->GetXaxis()->FindBin(hr->GetBinCenter(i));
    if (hzpt->GetBinContent(j)!=0) {
      hr->SetBinContent(i, (hgpt->GetBinContent(i)/hzpt->GetBinContent(j)-1)*100);
      hr->SetBinError(i, sqrt(pow(hgpt->GetBinError(i),2)+pow(hzpt->GetBinError(j),2))*100);
    }
  }

  
  TH1D *h_2 = tdrHist("h_2","#LTp_{T,jet}#GT / p_{T,Z/#gamma} (MPF)^{ }",
		      0.975, 1.025+1e-3, "p_{T,Z/#gamma}", 40, 300);
  //TH1D *h_2d = tdrHist("h_2d","(#gamma/Z+jet)-1 (%)", -1, +2.5,
  TH1D *h_2d = tdrHist("h_2d","(#gamma/Z+jet)-1 (%)", -1.0, +2.0,
		       "p_{T,Z/#gamma}", 40, 300);
  TCanvas *c2 = tdrDiCanvas("c2",h_2,h_2d,8,11);
  h_2->GetYaxis()->SetTitleOffset(1.04);
  
  c2->cd(1);
  gPad->SetLogx();

  l->SetLineStyle(kDashed);
  l->DrawLine(40,1,300,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(40,1.005,300,1.005);
  l->DrawLine(40,0.995,300,0.995);

  hgpt->GetXaxis()->SetRangeUser(40,300);
  hzpt->GetXaxis()->SetRangeUser(40,300);
  tdrDraw(hzpt,"Pz",kNone,kRed);
  tdrDraw(hgpt,"Pz",kNone,kBlue);
  
  TLegend *leg2 = tdrLeg(0.60,0.89-2*0.06,0.85,0.89);
  leg2->AddEntry(hgpt,"#gamma + jet","PLE");
  leg2->AddEntry(hzpt,"Z(#mu#mu) + jet","PLE");

  tex->DrawLatex(0.35,0.84,"|#eta_{jet}| < 1.3");
  
  c2->cd(2);
  gPad->SetLogx();
  
  l->SetLineStyle(kDashed);
  l->DrawLine(40,0,300,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(40,+0.5,300,+0.5);
  l->DrawLine(40,-0.5,300,-0.5);

  hr->GetXaxis()->SetRangeUser(40,300);
  tdrDraw(hr,"Pz",kNone,kBlue);
  
  c2->SaveAs(Form("pdf/drawGamVsJet/drawGamVsJet_pt_%s.pdf",year.c_str()));


  // Look at eta-asymmetry
  TH1D *h_3 = tdrHist("h_3","(#eta+) / (#eta-) - 1 (%)",retamin,retamax,
		      "|#eta_{#gamma}|",0,etamax);//1.5);
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0,0,etamax,0);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,+0.5,etamax,+0.5);
  l->DrawLine(0,-0.5,etamax,-0.5);
  //l->DrawLine(-0.8,retamin,-0.8,retamax-(retamax-retamin)*0.0);
  l->DrawLine(+0.8,retamin,+0.8,retamax-(retamax-retamin)*0.0);
  l->DrawLine(+1.4795,retamin,+1.4795,retamax-(retamax-retamin)*0.0);

  
  TH1D *hreta50 = (TH1D*)hg50->Clone("hreta50");
  TH1D *hreta110 = (TH1D*)hg110->Clone("hreta110");
  TH1D *hreta200 = (TH1D*)hg200->Clone("hreta200");
  for (int i = 1; i != hg50->GetNbinsX()+1; ++i) {
    double eta = hg50->GetXaxis()->GetBinCenter(i); // had a bug :(
    if (eta>0) {
      int j = hg50->GetXaxis()->FindBin(-eta);
      double reta50 = hg50->GetBinContent(j);
      if (reta50>0) {
	hreta50->SetBinContent(i, (hg50->GetBinContent(i)/hg50->GetBinContent(j)-1)*100);
	hreta50->SetBinError(i, sqrt(pow(hg50->GetBinError(i)/hg50->GetBinContent(i),2) + pow(hg50->GetBinError(j)/hg50->GetBinContent(j),2))*(1+0.01*hreta50->GetBinContent(i))*100);
      }
      double reta110 = hg110->GetBinContent(j);
      if (reta110>0) {
	hreta110->SetBinContent(i, (hg110->GetBinContent(i)/hg110->GetBinContent(j)-1)*100);
	hreta110->SetBinError(i, sqrt(pow(hg110->GetBinError(i)/hg110->GetBinContent(i),2) + pow(hg110->GetBinError(j)/hg110->GetBinContent(j),2))*(1+0.01*hreta110->GetBinContent(i))*100);
      }
      double reta200 = hg200->GetBinContent(j);
      if (reta200>0) {
	hreta200->SetBinContent(i, (hg200->GetBinContent(i)/hg200->GetBinContent(j)-1)*100);
	hreta200->SetBinError(i, sqrt(pow(hg200->GetBinError(i)/hg200->GetBinContent(i),2) + pow(hg200->GetBinError(j)/hg200->GetBinContent(j),2))*(1+0.01*hreta200->GetBinContent(i))*100);
      }
    } // eta>0
  } // for i


  tdrDraw(hreta50,"HIST",kNone,kGreen+2,kSolid,-1,kNone,0);
  tdrDraw(hreta200,"Pz",kNone,kRed);
  tdrDraw(hreta110,"Pz",kNone,kBlue);
  tdrDraw(hreta50,"Pz",kNone,kGreen+2,kSolid,-1,kNone,0);

  TLegend *leg3 = tdrLeg(0.60,0.90-0.05*3,0.85,0.90);
  leg3->AddEntry(hreta200,"Photon200EB","PLE");
  leg3->AddEntry(hreta110,"Photon110EB","PLE");
  leg3->AddEntry(hreta50,"Photon50EB","PLE");

  tex->DrawLatex(0.35,0.87,"#gamma + jet");
  tex->DrawLatex(0.35,0.81,"|#eta_{jet}| < 1.3");
  
  gPad->RedrawAxis();
  
  c3->SaveAs(Form("pdf/drawGamVsJet/drawGamVsJet_etaasymm_%s.pdf",year.c_str()));
  
} // drawGamVsJet
