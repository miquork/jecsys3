// Purpose: Use compareHcalDepths.root as a guide to extrapolate
//          results from CorrFact.root
#include "TFile.h"
#include "../tdrstyle_mod22.C"

void extrapolateHCalDepths() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/extrapolateHCalDepths");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/extrapolateHCalDepths");
  
  // B for barrel
  //TFile *fb = new TFile("rootfiles/CorrFact.root","READ");
  //TFile *fb = new TFile("rootfiles/CorrFact_25Oct29_EG25DEF.root","READ");
  TFile *fb = new TFile("rootfiles/CorrFact_25Sep17_EG25DEFC.root","READ");
  assert(fb && !fb->IsZombie());

  // E for endcap
  //TFile *fe = new TFile("rootfiles/compareHcalDepths.root","READ");
  TFile *fe = new TFile("rootfiles/compareHcalDepths_25Sep17_For_Yildiray.root","READ");
  assert(fe && !fe->IsZombie());

  curdir->cd();

  const int ndepth = 10;
  TCanvas *c1a = new TCanvas("c1a","c1a",5*400,2*400);
  c1a->Divide(5,2,0,0);

  TH1D* vhb[ndepth+1];
  TH1D* vhe[ndepth+1];
  TH1D* vhe2[ndepth+1];
  for (int depth = 1; depth != ndepth+1; ++depth) {
  
    TH1D *hb = (TH1D*)fb->Get(Form("hf_dd_%d",depth)); assert(hb);
    TH1D *he = (TH1D*)fe->Get(Form("h_depth_%d",depth)); assert(he);

    TH1D *h = tdrHist(Form("h_%d",depth),"Correction",0.5,3.5,"i#eta",-42,42);
    lumi_136TeV = "2025DEF, EGamma";
    extraText = "Private";
    //TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
    c1a->cd(depth);
    h->Draw("AXIS");
    
    tdrDraw(hb,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(he,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.5,1);

    vhb[depth] = hb;
    vhe[depth] = he;
  }

  TF1* vf1[ndepth+1];
  
  TCanvas *c2a = new TCanvas("c2a","c2a",5*400,2*400);
  c2a->Divide(4,2,0,0);
  const int nhedepth = 7;
  const double depthscale[nhedepth+1] =
    {0,  1.40, 1.07, 1.07, 1.05,  1.00, 1.00, 1.00};
  for (int depth = 1; depth != nhedepth+1; ++depth) {
    
    TH1D *hb = vhb[depth]; assert(hb);
    TH1D *hb_sub1 = vhb[max(depth-1,1)];
    TH1D *he = vhe[depth]; assert(he);
    TH1D *he_sub1 = vhe[max(depth-1,1)];
    
    TH1D *hb_ref = (TH1D*)hb->Clone(Form("hb_ref_%d",depth));
    TH1D *hb_inv = (TH1D*)hb->Clone(Form("hb_inv_%d",depth));
    TH1D *hb_avg = (TH1D*)hb->Clone(Form("hb_avg_%d",depth));
    TH1D *hb2 = (TH1D*)hb->Clone(Form("hb2_%d",depth));
    
    TH1D *he_ref = (TH1D*)he->Clone(Form("he_ref_%d",depth));
    TH1D *he_inv = (TH1D*)he->Clone(Form("he_inv_%d",depth));
    TH1D *he_avg = (TH1D*)he->Clone(Form("he_avg_%d",depth));
    TH1D *he2 = (TH1D*)he->Clone(Form("he2_%d",depth));

    int i1 = he->GetXaxis()->FindBin(20);
    int i2 = he->GetXaxis()->FindBin(30);    
    for (int i = i1; i != i2+1; ++i) {
      double c(0), cb(0);
      int ieta = int(he->GetXaxis()->GetBinCenter(i)+0.5);
      int j = he->GetXaxis()->FindBin(-ieta);
      int ib = hb->GetXaxis()->FindBin(+ieta);
      int jb = hb->GetXaxis()->FindBin(-ieta);
      if (depth==1 || depth==2 || depth==7 || ieta<=25) {
	c = 0.5 * (he->GetBinContent(i) + he->GetBinContent(j));
	cb = 0.5 * (hb->GetBinContent(ib) + hb->GetBinContent(jb));
      }
      if (depth==3 && ieta>=26) {
	c = 0.5 * 0.5 * (he->GetBinContent(i) + he->GetBinContent(j)) +
	  0.5 * 0.5 * (he_sub1->GetBinContent(i) + he_sub1->GetBinContent(j));
	cb = 0.5 * 0.5 * (hb->GetBinContent(ib) + hb->GetBinContent(jb)) +
	  0.5 * 0.5 * (hb_sub1->GetBinContent(ib) + hb_sub1->GetBinContent(jb));
      }
      if ((depth==4 || depth==5 || depth==6 || depth==7) && ieta>=26) {
	c = 0.5 * (he_sub1->GetBinContent(i) + he_sub1->GetBinContent(j));
	cb = 0.5 * (hb_sub1->GetBinContent(ib) + hb_sub1->GetBinContent(jb));
      }

      double k = depthscale[depth];
      if (depth==1 && ieta==28) k *= 1./1.5;//0.75;
      if (depth==1 && ieta==29) k *= 1./1.5;//0.75;
      if (ieta>=26) he2->SetBinContent(i, k*c);
      else          he2->SetBinContent(i, 0);
      he_ref->SetBinContent(i, k*he->GetBinContent(i));
      he_inv->SetBinContent(i, k*he->GetBinContent(j));
      he_avg->SetBinContent(i, k*0.5*(he->GetBinContent(i)+he->GetBinContent(j)));
      if (ieta>=26) hb2->SetBinContent(ib, cb);
      else          hb2->SetBinContent(ib, 0);
      hb_ref->SetBinContent(ib, hb->GetBinContent(ib));
      hb_inv->SetBinContent(ib, hb->GetBinContent(jb));
      hb_avg->SetBinContent(i, 0.5*(hb->GetBinContent(ib)+hb->GetBinContent(jb)));

      double ke = 0.2;
      hb_ref->SetBinError(ib, ke*hb->GetBinError(ib));
      hb_inv->SetBinError(ib, ke*hb->GetBinError(jb));
      hb_avg->SetBinError(ib, 0.5*ke*(hb->GetBinError(ib)+hb->GetBinError(jb)));
      hb2->SetBinError(ib, 0.5*ke*(hb->GetBinError(ib)+hb->GetBinError(jb)));
    }
    
    TH1D *h = tdrHist(Form("h2_%d",depth),"Correction",
		      //0.5+1e-4,(depth<=4 ? 3.5 : 2.0)-1e-4,"i#eta",19.5,30);
		      0.5+1e-4,(depth<=4 ? 2.5 : 2.0)-1e-4,"i#eta",19.5,29.5);
    lumi_136TeV = "2025DEF, EGamma";
    extraText = "Private";
    //TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
    c2a->cd(depth);
    h->Draw("AXIS");

    TLine *l = new TLine();
    l->SetLineStyle(kDashed);
    l->SetLineColor(kGray+1);
    l->DrawLine(19.5,1,29.5,1);
    if (depth==1) l->DrawLine(27.5,0.5,27.5,2.5);

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.05*1.5);
    double x = (depth%4==1 ? 0.25 : 0.17);
    tex->DrawLatex(x,0.88,Form("Depth %d",depth));
    tex->DrawLatex(x,0.80,Form("#gamma+jet scale %1.2f",depthscale[depth]));
    if (depth==1) tex->DrawLatex(x,0.72,"|i#eta|=28,29 #times 1/1.5");
    
    tdrDraw(hb_ref,"Pz",kOpenTriangleUp,kRed,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(hb_inv,"Pz",kOpenTriangleDown,kRed,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(hb_avg,"Pz",kFullSquare,kRed,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(hb2,"Pz",kOpenSquare,kMagenta+1,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(he_ref,"Pz",kOpenTriangleUp,kBlue,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(he_inv,"Pz",kOpenTriangleDown,kBlue,kSolid,-1,kNone,0,0.7,1);
    tdrDraw(he_avg,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,1.0,1);
    tdrDraw(he2,"Pz",kFullDiamond,kGreen+2,kSolid,-1,kNone,0,1.0,1);

    TF1 *f1 = new TF1(Form("f1_%d",depth),"[0]",20.5,25.5);
    hb_avg->Fit(f1,"QRN");
    f1->SetLineColor(kBlack);
    f1->SetLineWidth(3);
    f1->Draw("SAME");
    vf1[depth] = f1;
    
    TF1 *f2 = new TF1(Form("f2_%d",depth),"[0]",25.5,29.5);
    if (depth==1) f2->SetParameter(0, f1->GetParameter(0));
    if (depth==2) f2->SetParameter(0, f1->GetParameter(0));
    if (depth==3) f2->SetParameter(0, 0.5*(f1->GetParameter(0)+
					   vf1[depth-1]->GetParameter(0)));
    if (depth==4) f2->SetParameter(0, vf1[depth-1]->GetParameter(0));
    if (depth==5) f2->SetParameter(0, vf1[depth-1]->GetParameter(0));
    if (depth==6) f2->SetParameter(0, vf1[depth-1]->GetParameter(0));
    if (depth==7) f2->SetParameter(0, vf1[depth-1]->GetParameter(0));
    f2->SetLineColor(kBlack);
    f2->SetLineStyle(kDotted);
    f2->SetLineWidth(3);
    f2->Draw("SAME");

    if (depth==7) {
      c2a->cd(8);
      tex->DrawLatex(0.15, 0.87, "Prompt EGamma 2025DEF");
      TLegend *leg = tdrLeg(0.25,0.82-0.05*1.5*10,0.60,0.82);
      leg->SetTextSize(0.045*1.5);
      leg->AddEntry(hb_avg,"IsoTrack |#eta|","PLE");
      leg->AddEntry(hb_ref,"IsoTrack #eta+","PLE");
      leg->AddEntry(hb_inv,"IsoTrack #eta-","PLE");
      leg->AddEntry(he_avg,"#gamma+jet |#eta|","PLE");
      leg->AddEntry(he_ref,"#gamma+jet #eta+","PLE");
      leg->AddEntry(he_inv,"#gamma+jet #eta-","PLE");
      leg->AddEntry(hb2,"IsoTrack Scheme","PLE");
      leg->AddEntry(he2,"#gamma+jet Scheme","PLE");
      leg->AddEntry(f1,"21 #leq i#eta #leq 25 fit","L");
      leg->AddEntry(f2,"i#eta #geq 26 extrap.","L");
    }
  } // for depth

  c2a->SaveAs("pdf/extrapolateHCalDepths/extrapolateHCalDepths.pdf");
} // extrapolateHCalDepths
