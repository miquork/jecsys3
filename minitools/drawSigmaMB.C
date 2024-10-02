// Purpose draw estimate of sigmaMB scale factor based on gamjet+dijet+Zjet data
#include "TH1D.h"

#include "../tdrstyle_mod22.C"

#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <iostream>

bool useModObs = true;

int drawSigmaMB() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  // Read in results
  //ifstream fd("rootfiles/dijet_pileup_mean-values_2024DEFG_nestor.txt");
  ifstream fd("rootfiles/dijet_pileup_mean-values_2024DEFG_nestor_v2.txt");
  if (!fd.is_open()) return -2;

  ifstream fg("rootfiles/gamjet_pileup_mean-values_2024DEF_bettina.txt");
  if (!fd.is_open()) return -2;

  ifstream fz("rootfiles/zmmjet_pileup_mean-values_2024DEFG_sami.txt");
  if (!fz.is_open()) return -2;
  
  char buf[512];
  fd.getline(buf,512);
  //cout << buf << endl << flush;
  fg.getline(buf,512);
  fz.getline(buf,512);
  cout << buf << endl << flush;
  
  char cera[256], cmod[256], cobs[256];
  float mean, err;
  map<string, map<string, map<string, pair<float,float> > > > md;
  while (fd.getline(buf,512)) {
    //cout << buf << flush;
    int n = sscanf(buf,"%s %s %s %f %f",cera,cmod,cobs,&mean,&err);
    //cout << " => n="<<n << endl << flush;
    if (n!=5) return -1;
    md[cobs][cera][cmod] = make_pair(mean, err);
  }
  map<string, map<string, map<string, pair<float,float> > > > mg;
  while (fg.getline(buf,512)) {
    //cout << buf << flush;
    int n = sscanf(buf,"%s %s %s %f %f",cera,cmod,cobs,&mean,&err);
    //cout << " => n="<<n << endl << flush;
    if (n!=5) return -1;
    mg[cobs][cera][cmod] = make_pair(mean, err);
  }
  map<string, map<string, map<string, pair<float,float> > > > mz;
  while (fz.getline(buf,512)) {
    //cout << buf << flush;
    int n = sscanf(buf,"%s %s %s %f %f",cera,cmod,cobs,&mean,&err);
    //cout << " => n="<<n << endl << flush;
    if (n!=5) return -1;
    mz[cobs][cera][cmod] = make_pair(mean, err);
  }

  
  int nera = md.size();
  cout << "Found " << nera << " eras in dijet" << endl << flush;
  int nerag = mg.size();
  cout << "Found " << nerag << " eras in gamjet" << endl << flush;
  int neraz = mz.size();
  cout << "Found " << neraz << " eras in zmmjet" << endl << flush;

  TH1D *hnpv = new TH1D("hnpv",";#sigma_{MB} (mb);Observations",20,65,85);
  TH1D *hrho = new TH1D("hrho",";#sigma_{MB} (mb);Observations",20,65,85);
  
  TH1D *h = tdrHist("h","Data / MC",0.87,1.42,"Era",0.5,nera+0.5);
  //TH1D *h = tdrHist("h","Data / MC",0.65,1.55,"Era",0.5,nera+0.5);
  lumi_136TeV = "Prompt2024";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0.5,1,0.5+nera,1);
  double sf = 1.06;
  l->SetLineStyle(kDotted);
  l->DrawLine(0.5,sf,0.5+nera,sf);

  // TWiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Recommended_cross_section
  // https://arxiv.org/pdf/2011.14909 (Fig.2, Eq.6)
  TF1 *f1 = new TF1("f1","[0]+[1]*pow(log(x*x),[2])",0.01,100);
  //f1->SetParameters(28.84,0.05,2.37); // Fig.2
  f1->SetParameters(28.84,0.0458,2.374); // Eq.(6)
  double sigmamb8 = f1->Eval(8000.);
  double sigmamb13 = f1->Eval(13000.);
  double sigmamb136 = f1->Eval(13600.);
  cout << "sigmaMB(8 TeV)="<<sigmamb8 << endl;
  cout << "sigmaMB(13 TeV)="<<sigmamb13 << endl;
  cout << "sigmaMB(13.6 TeV)="<<sigmamb136 << endl;
  double kt = sigmamb136/sigmamb8;

  l->SetLineStyle(kSolid);
  l->DrawLine(0.5,kt,0.5+nera,kt);
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  tex->DrawLatex(0.19,0.75,"Default #sigma_{MB}=69.2 mb");
  tex->DrawLatex(0.19,0.71,Form("Better #sigma_{MB}=%1.1f mb",69.2*sf));
  tex->DrawLatex(0.19,0.67,Form("Best #sigma_{MB}=%1.1f mb",69.2*kt));
  //tex->DrawLatex(0.19,0.71,"HIST before reweighing");
  //tex->DrawLatex(0.19,0.67,"MARK after reweighing");
 
  TLegend *leg = tdrLeg(0.68,0.90-(nera+1)*0.05,0.93,0.90);
  leg->SetHeader("J");
  TLegend *legg = tdrLeg(0.61,0.90-(nera+1)*0.05,0.86,0.90);
  legg->SetHeader("#gamma");
  //TLegend *legz = tdrLeg(0.54,0.90-(nera+1)*0.05,0.79,0.90);
  //legz->SetHeader("#gamma");
  
  typedef map<string, map<string, map<string, pair<float,float> > > >::iterator IT;
  typedef map<string, map<string, pair<float,float> > >::iterator JT;

  int n(0);
  for (IT it = md.begin(); it != md.end(); ++it) {

    const char *cobs = it->first.c_str();
    TH1D *hobs = new TH1D(Form("h_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);
    TH1D *hobs2 = new TH1D(Form("h2_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);

    int j(0);
    for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {

      const char *cera = jt->first.c_str();
      h->GetXaxis()->SetBinLabel(++j, jt->first.c_str());
      h->GetXaxis()->SetLabelSize(0.07);//0.05);

      double mc = jt->second["mc"].first;
      double mcpuoff = jt->second["mcpuoff"].first;
      double dt = jt->second["data"].first;
      double edt = jt->second["data"].second;

      // Subtract hard scatter
      string sobs = string(cobs);
      if (useModObs && (sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	mc -= 1;
	mcpuoff -= 1;
	dt -= 1;
      }
      if (useModObs && sobs=="rho" && mc!=0) {
	mc -= 2.0;
	mcpuoff -= 2.0;
	dt -= 2.0;
      }
      
      hobs->SetBinContent(j, (mc!=0 ? dt/mc : 0));
      //hobs->SetBinError(j, (mc!=0 ? edt/mc : 0));
      hobs2->SetBinContent(j, (mcpuoff!=0 ? dt/mcpuoff : 0));

      string sera = string(cera);
      if ((sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	hnpv->Fill(dt/mc*69.2);
      }
      if ((sobs=="rho" && sera!="D") && mc!=0) {
	hrho->Fill(dt/mc*69.2);
      }
    } // for JT

    //tdrDraw(hobs2,"HIST][",kFullCircle+n,kBlack+n,kSolid,-1,kNone,0,0.8);
    tdrDraw(hobs,"Pz",kFullCircle+n,kBlack+n,kSolid,-1,kNone,0,0.8); ++n;
    leg->AddEntry(hobs,cobs,"PLE");
  } // for IT

  int ng(0);
  for (IT it = mg.begin(); it != mg.end(); ++it) {

    const char *cobs = it->first.c_str();
    TH1D *hobs = new TH1D(Form("hg_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);
    TH1D *hobs2 = new TH1D(Form("h2g_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);

    int j(0);
    for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {

      const char *cera = jt->first.c_str();
      h->GetXaxis()->SetBinLabel(++j, jt->first.c_str());

      double mc = jt->second["mc"].first;
      if (string(cobs)=="mu") mc = jt->second["mcpuoff"].first;
      double mcpuoff = jt->second["mcpuoff"].first;
      if (string(cobs)=="mu") mcpuoff = jt->second["mc"].first;
      double dt = jt->second["data"].first;
      double edt = jt->second["data"].second;

      // Subtract hard scatter
      string sobs = string(cobs);
      if (useModObs && (sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	mc -= 1;
	mcpuoff -= 1;
	dt -= 1;
      }
      if (useModObs && sobs=="rho" && mc!=0) {
	mc -= 2.0;
	mcpuoff -= 2.0;
	dt -= 2.0;
      }
      
      hobs->SetBinContent(j, (mc!=0 ? dt/mc : 0));
      //hobs->SetBinError(j, (mc!=0 ? edt/mc : 0));
      hobs2->SetBinContent(j, (mcpuoff!=0 ? dt/mcpuoff : 0));

      string sera = string(cera);
      if ((sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	hnpv->Fill(dt/mc*69.2);
      }
      if ((sobs=="rho" && sera!="D") && mc!=0) {
	hrho->Fill(dt/mc*69.2);
      }
    } // for JT

    //tdrDraw(hobs2,"HIST][",kOpenCircle,kBlack+ng,kDashed,-1,kNone,0,1.2);
    tdrDraw(hobs,"Pz",kOpenCircle+ng,kBlack+ng,kSolid,-1,kNone,0,1.2); ++ng;
    legg->AddEntry(hobs," ","PLE");
  } // for IT

  int nz(0);
  for (IT it = mz.begin(); it != mz.end(); ++it) {

    const char *cobs = it->first.c_str();
    TH1D *hobs = new TH1D(Form("hz_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);
    TH1D *hobs2 = new TH1D(Form("h2z_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);

    int j(0);
    for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {

      const char *cera = jt->first.c_str();
      h->GetXaxis()->SetBinLabel(++j, jt->first.c_str());

      double mcpuoff = jt->second["mcpuoff"].first;
      double pumc = md["mu"][cera]["mcpuoff"].first;
      double purw = md["mu"][cera]["mc"].first;
      double mc = (pumc!=0 ? mcpuoff * purw / pumc : mcpuoff);
      double dt = jt->second["data"].first;
      double edt = jt->second["data"].second;

      // Subtract hard scatter
      string sobs = string(cobs);
      if (useModObs && (sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	mc -= 1;
	mcpuoff -= 1;
	dt -= 1;
      }
      if (useModObs && sobs=="rho" && mc!=0) {
	mc -= 2.0;
	mcpuoff -= 2.0;
	dt -= 2.0;
      }
      
      hobs->SetBinContent(j, (mc!=0 ? dt/mc : 0));
      hobs2->SetBinContent(j, (mcpuoff!=0 ? dt/mcpuoff : 0));

      string sera = string(cera);
      if ((sobs=="npvall" || sobs=="npvgood") && mc!=0) {
	hnpv->Fill(dt/mc*69.2);
      }
      if ((sobs=="rho" && sera!="D") && mc!=0) {
	hrho->Fill(dt/mc*69.2);
      }
    } // for JT

    //tdrDraw(hobs2,"HIST][",kOpenCircle,kBlack+nz,kDashed,-1,kNone,0,1.2);
    //tdrDraw(hobs,"Pz",kOpenDiamond+1+nz,kBlack+nz,kSolid,-1,kNone,0,1.2); ++nz;
    //legz->AddEntry(hobs," ","PLE");
  } // for IT

  tex->DrawLatex(0.36,0.16,"E: HF scale up (rho up)");
  tex->DrawLatex(0.56,0.21,"F: HCALDI on (rho up)");


  c1->SaveAs("pdf/drawSigmaMB/drawSigmaMB.pdf");


  TH1D *h2 = tdrHist("h2","Observations",0,5.5,"#sigma_{MB} (mb)",65,85);
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);

  l->SetLineStyle(kDashed);
  l->DrawLine(69.2,0,69.2,4.);
  l->DrawLine(75.3,0,75.3,1.8);
  l->DrawLine(75.3,2.45,75.3,3.7);
  
  tdrDraw(hnpv,"HIST",kFullSquare,kRed,kSolid,kRed+1,1001,kRed);
  hnpv->SetFillColorAlpha(kRed,0.70);
  tdrDraw(hrho,"HIST",kFullTriangleDown,kBlue,kSolid,kBlue+1,1001,kBlue);
  hrho->SetFillColorAlpha(kBlue,0.70);

  TGraphErrors *gnpv = new TGraphErrors(1);
  gnpv->SetPoint(0,hnpv->GetMean(),2.5);
  gnpv->SetPointError(0,hnpv->GetRMS(),0);
  TGraphErrors *grho = new TGraphErrors(1);
  grho->SetPoint(0,hrho->GetMean(),2.5);
  grho->SetPointError(0,hrho->GetRMS(),0);
  TGraphErrors *gsmb = new TGraphErrors(1);
  gsmb->SetPoint(0,0.5*(hrho->GetMean()+hnpv->GetMean()),2.25);
  gsmb->SetPointError(0,0.5*(hrho->GetRMS()+hnpv->GetRMS()),0);

  tdrDraw(gnpv,"Pz",kFullSquare,kRed+1); gnpv->SetLineWidth(2);
  tdrDraw(grho,"Pz",kFullTriangleDown,kBlue+1); grho->SetLineWidth(2);
  tdrDraw(gsmb,"Pz",kFullCircle,kBlack); gsmb->SetLineWidth(2);
  tex->DrawLatex(0.55,0.41,Form("#sigma_{MB} = %1.1f#pm%1.1f mb",
			       gsmb->GetX()[0],gsmb->GetEX()[0]));
  tex->DrawLatex(0.30,0.75,"#sigma_{MB} = 69.2 mb (default)");
  tex->DrawLatex(0.54,0.70,"#sigma_{MB} = 75.3 mb (theory)");

  TLegend *leg2 = tdrLeg(0.53,0.90-2*0.05,0.78,0.90);
  leg2->AddEntry(hnpv,"Vertices (N_{PV})","FPL");
  leg2->AddEntry(hrho,"Offset density (#rho)","FPL");
  
  gPad->RedrawAxis();
  c2->SaveAs("pdf/drawSigmaMB/drawSigmaMB_dist.pdf");
  
  return 0;
} // drawSigmaMB
