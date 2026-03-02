// Purpose: Draw jet asymmetry vs phi (from gamma+jet, dijet) for HCAL DPG
//          Goal is to show large variation vs iphi
//          This could be potential source of large JER SF
//          Could also add jet rate
// Acshually: NEF also has big holes correlating to low response, so
//            ECAL intercalibration is a mess also?
//            Need to replace p2asymm by energy component
#include "TProfile2D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

TH1D *getPEF(TFile *f, double eta) {
  assert(f);
  TProfile2D *p2 = (TProfile2D*)f->Get("Jetveto/p2asymm_noveto");
  assert(p2);
  //TProfile2D *p2nhf = (TProfile2D*)f->Get("Jetveto/p2nhf");
  TProfile2D *p2nhf = (TProfile2D*)f->Get("Jetveto/p2nef");
  assert(p2nhf);

  int ieta = p2->GetXaxis()->FindBin(eta);
  TProfile *p = p2->ProfileY("p",ieta,ieta);
  TProfile *p_nhf = p2nhf->ProfileY("p_nhf",ieta,ieta);

  TF1 *f1 = new TF1("f1","[0]",-TMath::Pi(),TMath::Pi());
  p->Fit(f1,"QRN");
  double r0 = (1+f1->GetParameter(0));
  p_nhf->Fit(f1,"QRN");
  double f0 = f1->GetParameter(0);
  
  TH1D *h = p->ProjectionX(Form("h_%s%04d",(eta>0 ? "p" : ","),int(1000.*eta)));
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    // p is asymmetry, i.e. (pTjet-pTgam)/pTgam = pTjet/pTgam-1
    double r = 1+p->GetBinContent(i);
    double er = p->GetBinError(i);
    double f = p_nhf->GetBinContent(i);
    double ef = p_nhf->GetBinError(i);
    double f0 = f1->GetParameter(0);

    double y = (r*f - r0*f0);// / (r0*f0);
    double ey = sqrt(pow(er*f,2) + pow(r*ef,2));// / (r0*f0);
    h->SetBinContent(i, y);
    h->SetBinError(i, ey);
  }
  
  delete p;
  delete p_nhf;
  delete f1;
  
  return h;
}

void drawJetRespVsPhi() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  //TFile *fg = new TFile("rootfiles/Prompt2024/w48_Gam/GamHistosFill_data_2024Inib1_w48.root","READ");
  //TFile *fg = new TFile("rootfiles/Prompt2024/w48_Gam/GamHistosFill_data_2024BCDEFGHI_w48.root","READ");
  TFile *fg = new TFile("rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root","READ");
  //TFile *fg = new TFile("rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_2022P8-PTG_no-pu_w48.root","READ");
  assert(fg && !fg->IsZombie());

  /*
  TProfile2D *p2 = (TProfile2D*)fg->Get("Jetveto/p2asymm_noveto"); assert(p2);
  TProfile2D *p2nhf = (TProfile2D*)fg->Get("Jetveto/p2nhf"); assert(p2nhf);

  // ieta=28,29, ref ieta=25
  // see pdf/Gains2025_preX2.pdf
  // 2.85-3.0
  int ieta29p = p2->GetXaxis()->FindBin(2.9);
  TProfile *p29p = p2->ProfileY("p29p",ieta29p,ieta29p);
  TProfile *p29p_nhf = p2nhf->ProfileY("p29p_nhf",ieta29p,ieta29p);
  int ieta29m = p2->GetXaxis()->FindBin(-2.9);
  TProfile *p29m = p2->ProfileY("p29m",ieta29m,ieta29m);
  TProfile *p29m_nhf = p2nhf->ProfileY("p29m_nhf",ieta29m,ieta29m);
  // 2.65-2.85
  int ieta28p = p2->GetXaxis()->FindBin(2.7);
  TProfile *p28p = p2->ProfileY("p28p",ieta28p,ieta28p);
  TProfile *p28p_nhf = p2nhf->ProfileY("p28p_nhf",ieta28p,ieta28p);
  int ieta28m = p2->GetXaxis()->FindBin(-2.7);
  TProfile *p28m = p2->ProfileY("p28m",ieta28m,ieta28m);
  TProfile *p28m_nhf = p2nhf->ProfileY("p28m_nhf",ieta28m,ieta28m);
  // 2.5-2.65: 27, 2.322-2.5: 26
  int ieta27p = p2->GetXaxis()->FindBin(2.55);
  TProfile *p27p = p2->ProfileY("p27p",ieta27p,ieta27p);
  TProfile *p27p_nhf = p2nhf->ProfileY("p27p_nhf",ieta27p,ieta27p);
  int ieta27m = p2->GetXaxis()->FindBin(-2.4);
  TProfile *p27m = p2->ProfileY("p27m",ieta27m,ieta27m);
  TProfile *p27m_nhf = p2nhf->ProfileY("p27m_nhf",ieta27m,ieta27m);
  */
  TH1D *p29p = getPEF(fg, +2.9);
  TH1D *p29m = getPEF(fg, -2.9);
  TH1D *p28p = getPEF(fg, +2.7);
  TH1D *p28m = getPEF(fg, -2.7);
  TH1D *p27p = getPEF(fg, +2.55);
  TH1D *p27m = getPEF(fg, -2.55);
  
  //TH1D *h = tdrHist("h","Asymmetry #LTp_{T,jet}/p_{T,#gamma}#GT-1",-0.30,0.25,
  //TH1D *h = tdrHist("h","HCAL asymmetry #LTp_{T,HCAL}/p_{T,#gamma}#GT-1",
  TH1D *h = tdrHist("h","ECAL asymmetry #LTp_{T,HCAL}/p_{T,#gamma}#GT-1",
		    -0.40,0.40,"#phi_{jet}",-TMath::Pi(),+TMath::Pi());
  extraText = "Private";
  lumi_136TeV = "#gamma+jet EB50, 2024";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(-TMath::Pi(),0,+TMath::Pi(),0);
  l->DrawLine(-TMath::Pi()/2.,-0.2,-TMath::Pi()/2.,0.2);
  l->DrawLine(0.,-0.2,0.,0.2);
  l->DrawLine(+TMath::Pi()/2.,-0.2,+TMath::Pi()/2.,0.2);
  
  TLegend *leg = tdrLeg(0.40,0.90-0.05*3,0.90,0.90);
  leg->SetNColumns(2);

  tdrDraw(p29m,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0);
  tdrDraw(p29p,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0);
  tdrDraw(p28m,"Pz",kOpenCircle,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p28p,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0);
  tdrDraw(p27m,"Pz",kOpenDiamond,kGreen+2,kSolid,-1,kNone,0);
  tdrDraw(p27p,"Pz",kFullDiamond,kGreen+2,kSolid,-1,kNone,0);

  leg->AddEntry(p29m,"i#eta=-29","PLE");
  leg->AddEntry(p29p,"i#eta=+29","PLE");
  leg->AddEntry(p28m,"i#eta=-28","PLE");
  leg->AddEntry(p28p,"i#eta=+28","PLE");
  leg->AddEntry(p27m,"i#eta=-27","PLE");
  leg->AddEntry(p27p,"i#eta=+27","PLE");

  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawJetRespVsPhi/drawJetRespVsPhi.pdf");


  ////////////////////////////////////////////////////////
  // Step 2. Populate RMS distribution and expected RMS //
  ////////////////////////////////////////////////////////

  if (true) {
    
    //TH1D *h_2 = tdrHist("h_2","HCAL towers",0,55,"HCAL asymmetry",-0.30,0.40);
    TH1D *h_2 = tdrHist("h_2","ECAL towers",0,55,"ECAL asymmetry",-0.30,0.40);
  extraText = "Private";
  lumi_136TeV = "#gamma+jet EB50, 2024";
  TCanvas *c2 = tdrCanvas("c2",h_2,8,11,kSquare);

  TLegend *leg2 = tdrLeg(0.57,0.90-0.05*6,0.82,0.90);


  TH1D *h29 = new TH1D("h29",";Asymmetry;HCAL towers",70,-0.35,0.35);
  TH1D *h28 = new TH1D("h28",";Asymmetry;HCAL towers",70,-0.35,0.35);
  TH1D *h27 = new TH1D("h27",";Asymmetry;HCAL towers",70,-0.35,0.35);

  TH1D *h29s = new TH1D("h29s",";Stat.;HCAL towers",70,-0.35,0.35);
  TH1D *h28s = new TH1D("h28s",";Stat.;HCAL towers",70,-0.35,0.35);
  TH1D *h27s = new TH1D("h27s",";Stat.;HCAL towers",70,-0.35,0.35);

  TRandom3 rnd;
  for (int i = 1; i != p29m->GetNbinsX()+1; ++i) {
    h29->Fill(p29m->GetBinContent(i));
    h29->Fill(p29p->GetBinContent(i));
    h28->Fill(p28m->GetBinContent(i));
    h28->Fill(p28p->GetBinContent(i));
    h27->Fill(p27m->GetBinContent(i));
    h27->Fill(p27p->GetBinContent(i));

    h29s->Fill(rnd.Gaus(0.,p29m->GetBinError(i)));
    h29s->Fill(rnd.Gaus(0.,p29p->GetBinError(i)));
    h28s->Fill(rnd.Gaus(0.,p28m->GetBinError(i)));
    h28s->Fill(rnd.Gaus(0.,p28p->GetBinError(i)));
    h27s->Fill(rnd.Gaus(0.,p27m->GetBinError(i)));
    h27s->Fill(rnd.Gaus(0.,p27p->GetBinError(i)));
  }

  tdrDraw(h27s,"HIST",kNone,kGreen+2,kSolid,-1,3005,kGreen+2-9);
  h27s->SetFillColorAlpha(kGreen+2-9,0.70);
  tdrDraw(h28s,"HIST",kNone,kBlue,kSolid,-1,3005,kBlue-9);
  h28s->SetFillColorAlpha(kBlue-9,0.30);
  tdrDraw(h29s,"HIST",kNone,kRed,kSolid,-1,3005,kRed-9);
  h29s->SetFillColorAlpha(kRed-9,0.30);
  
  tdrDraw(h27,"HIST",kNone,kGreen+2,kSolid,-1,1001,kGreen+2-9);
  h27->SetFillColorAlpha(kGreen+2-9,0.70);
  tdrDraw(h28,"HIST",kNone,kBlue,kSolid,-1,1001,kBlue-9);
  h28->SetFillColorAlpha(kBlue-9,0.30);
  tdrDraw(h29,"HIST",kNone,kRed,kSolid,-1,1001,kRed-9);
  h29->SetFillColorAlpha(kRed-9,0.30);

  leg2->AddEntry(h27,Form("|i#eta|=27, #sigma=%1.1f%%",h27->GetRMS()*100.),"F");
  leg2->AddEntry(h28,Form("|i#eta|=28, #sigma=%1.1f%%",h28->GetRMS()*100.),"F");
  leg2->AddEntry(h29,Form("|i#eta|=29, #sigma=%1.1f%%",h29->GetRMS()*100.),"F");
  
  leg2->AddEntry(h27s,Form("#sigma_{stat}(27)=%1.2f%%",h27s->GetRMS()*100),"F");
  leg2->AddEntry(h28s,Form("#sigma_{stat}(28)=%1.2f%%",h28s->GetRMS()*100),"F");
  leg2->AddEntry(h29s,Form("#sigma_{stat}(29)=%1.2f%%",h29s->GetRMS()*100),"F");
  
  gPad->RedrawAxis();
  c2->SaveAs("pdf/drawJetRespVsPhi/drawJetRespVsPhi_Asymm.pdf");
  }

  ////////////////////////////////////////////////////////////////////
  // Step 3. Correlate high responses with high gains in depths 2,3 //
  ////////////////////////////////////////////////////////////////////

  TFile *fgain = new TFile("rootfiles/HE_gains_from_John.root","READ");
  assert(fgain && !fgain->IsZombie());

  TH2D *h2g2 = (TH2D*)fgain->Get("gh_2");
  TH2D *h2g3 = (TH2D*)fgain->Get("gh_3");

  int i29m = h2g3->GetXaxis()->FindBin(-28);//-29.);
  int i29p = h2g3->GetXaxis()->FindBin(+28);//+29.);

  TH1D *hg2p = h2g2->ProjectionY("hg2p",i29p,i29p);
  TF1 *f1gain2p = new TF1("f1gain2p","[0]",0,74);
  hg2p->Fit(f1gain2p,"QRNW");

  TH1D *hg2m = h2g2->ProjectionY("hg2m",i29m,i29m);
  TF1 *f1gain2m = new TF1("f1gain2m","[0]",0,74);
  hg2m->Fit(f1gain2m,"QRNW");

  TH1D *hg3p = h2g3->ProjectionY("hg3p",i29p,i29p);
  TF1 *f1gain3p = new TF1("f1gain3p","[0]",0,74);
  hg3p->Fit(f1gain3p,"QRNW");

  TH1D *hg3m = h2g3->ProjectionY("hg3p",i29m,i29m);
  TF1 *f1gain3m = new TF1("f1gain3m","[0]",0,74);
  hg3m->Fit(f1gain3m,"QRNW");

  /*
  TH1D *h29pg2 = p29p->ProjectionX("h29pg2"); h29pg2->Reset();
  TH1D *h29pg3 = p29p->ProjectionX("h29pg3"); h29pg3->Reset();
  TH1D *h29mg2 = p29m->ProjectionX("h29mg2"); h29mg2->Reset();
  TH1D *h29mg3 = p29m->ProjectionX("h29mg3"); h29mg3->Reset();
  */
  
  TH1D *h29pg2 = (TH1D*)p29p->Clone("h29pg2"); h29pg2->Reset();
  TH1D *h29pg3 = (TH1D*)p29p->Clone("h29pg3"); h29pg3->Reset();
  TH1D *h29mg2 = (TH1D*)p29m->Clone("h29mg2"); h29mg2->Reset();
  TH1D *h29mg3 = (TH1D*)p29m->Clone("h29mg3"); h29mg3->Reset();
  
  for (int i = 1; i != p29p->GetNbinsX()+1; ++i) {
    
    double phi = p29p->GetBinCenter(i)+TMath::Pi(); // [-pi,pi]->[0,2pi]
    int iphi = int(phi/(TMath::TwoPi()/72.))+1; // phi[0,2pi]->iphi[1,72]
    int biny = 2*((iphi-1)/2)+2; // iphi[1,72]->biny[2,4,..,70,72]

    double gain2p = hg2p->GetBinContent(biny);
    double normgain2p = gain2p / f1gain2p->GetParameter(0);
    double gain3p = hg3p->GetBinContent(biny);
    double normgain3p = gain3p / f1gain3p->GetParameter(0);
    
    h29pg2->SetBinContent(i, (normgain2p-1)*0.35);
    h29pg3->SetBinContent(i, (normgain3p-1)*0.15);

    double gain2m = hg2m->GetBinContent(biny);
    double normgain2m = gain2m / f1gain2m->GetParameter(0);
    double gain3m = hg3m->GetBinContent(biny);
    double normgain3m = gain3m / f1gain3m->GetParameter(0);
    
    h29mg2->SetBinContent(i, (normgain2m-1)*0.35);
    h29mg3->SetBinContent(i, (normgain3m-1)*0.15);
  }

  TGraphErrors *g29v3 = new TGraphErrors(72*2);
  assert(p29p->GetNbinsX()==72);
  for (int i = 1; i != 72+1; ++i) {

    double phi = p29p->GetBinCenter(i)+TMath::Pi(); // [-pi,pi]->[0,2pi]
    int iphi = int(phi/(TMath::TwoPi()/72.))+1; // phi[0,2pi]->iphi[1,72]
    int biny = 2*((iphi-1)/2)+2; // iphi[1,72]->biny[2,4,..,70,72]

    //double gain2m = h2g2->GetBinContent(i29m, biny);
    double gain2m = h29mg2->GetBinContent(biny);
    double egain2m =0;// h2g2->GetBinError(i29, biny);
    //double gain3m = h2g3->GetBinContent(i29m, biny);
    double gain3m = h29mg3->GetBinContent(biny);
    double egain3m =0;// h2g3->GetBinError(i29, biny);
    double gainm = (0.35*gain2m+0.15*gain3m)/(0.35+0.15);
    double egainm = 0;
    
    g29v3->SetPoint(i-1, gain2m, p29m->GetBinContent(i));
    g29v3->SetPointError(i-1, egain2m, p29m->GetBinError(i));

    //double gain2p = h2g2->GetBinContent(i29p, biny);
    double gain2p = h29pg2->GetBinContent(biny);
    double egain2p =0;// h2g2->GetBinError(i29, biny);
    //double gain3p = h2g3->GetBinContent(i29p, biny);
    double gain3p = h29pg3->GetBinContent(biny);
    double egain3p =0;// h2g3->GetBinError(i29, biny);
    double gainp = (0.35*gain2p+0.15*gain3p)/(0.35+0.15);
    double egainp = 0;

    /*
    g29v3->SetPoint(i-1, gain2p, p29p->GetBinContent(i));
    g29v3->SetPointError(i-1, egain2p, p29p->GetBinError(i));
    */
    g29v3->SetPoint(i-1, gain2p, p28p->GetBinContent(i));
    g29v3->SetPointError(i-1, egain2p, p28p->GetBinError(i));
  }

  //TH1D *h_3 = tdrHist("h_3","HCAL asymmetry",-0.30,0.40,"Gain",0.000,0.010);
  TH1D *h_3 = tdrHist("h_3","ECAL asymmetry",-0.30,0.40,"Gain",0.000,0.010);
  extraText = "Private";
  lumi_136TeV = "#gamma+jet EB50, 2024";
  TCanvas *c3 = tdrCanvas("c3",h_3,8,11,kSquare);

  TLegend *leg3 = tdrLeg(0.57,0.90-0.05*6,0.82,0.90);

  //tdrDraw(g29v3,"Pz",kFullCircle,kBlack);
  g29v3->Draw("APz");
  /*
  //hg2p->Draw("HIST");
  h29g2->Draw("HIST");
  h29g3->Draw("HIST SAME");
  f1gain2p->Draw("SAME");
  f1gain3p->Draw("SAME");
  */
  
  c1->cd();
  h29pg2->SetLineColor(kRed);
  h29pg2->Draw("HIST SAME");
  h29pg3->SetLineColor(kRed-9);
  h29pg3->Draw("HIST SAME");

  TF1 *f14p = new TF1("f14p","[0]+[1]*cos(4*x)*fabs(cos(4*x))",
		      -TMath::Pi(),TMath::Pi());
  p29p->Fit(f14p,"QRN");
  f14p->Draw("SAME");
  TF1 *f14m = new TF1("f14m","[0]+[1]*cos(4*x)*fabs(cos(4*x))",
		      -TMath::Pi(),TMath::Pi());
  p29m->Fit(f14m,"QRN");
  f14m->SetLineColor(kRed-9);
  f14m->Draw("SAME");

} // drawJetRespVsPhi


