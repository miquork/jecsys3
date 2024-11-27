// Purpose draw estimate of sigmaMB scale factor based on gamjet+dijet+Zjet data
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"

#include "../tdrstyle_mod22.C"

#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <iostream>

bool useModObs = true;
bool drawPuoff = false;//true;//false;

int drawSigmaMBV2() {

  TDirectory *curdir = gDirectory;
  setTDRStyle();

  string vera[] = {
    //"24B","24C","24D","24E",
    //"24BCD","24E",
    //"24F","24G","24H","24I"
    "2022-2023","2024"
  };
  const int nera0 = sizeof(vera)/sizeof(vera[0]);

  // Read in results
  char cera2[256], cmod[256], cobs[256];
  float mean, err;
  map<string, map<string, map<string, pair<float,float> > > > md;
  map<string, map<string, map<string, pair<float,float> > > > mg;
  map<string, map<string, map<string, pair<float,float> > > > mz;

  for (int i = 0; i != nera0; ++i) {

    const char *cera = vera[i].c_str();
    //string name = Form("rootfiles/PU2024_Xsec75mb/20%s/histogram_means.csv",cera);
    //string name = Form("rootfiles/PURW_2024/20%s/histogram_means.csv",cera);
    string name = Form("rootfiles/histogram_means_%s.csv",cera);
    ifstream fd(name.c_str());
    if (!fd.is_open()) {
      cout << "\nCould not open " << name << endl << flush;
      return -2*(i+1);
    }

    char buf[512];
    // Remove header line
    fd.getline(buf,512);
    cout << buf << endl << flush;

    // Read in data
    while (fd.getline(buf,512)) {
      //cout << buf << endl << flush;
      TString ts(buf);
      ts.ReplaceAll(","," ");
      ts.ReplaceAll(";"," ");
      cout << ts << endl << flush;
      int n = sscanf(ts.Data(),"%s %s %s %f %f",cera2,cmod,cobs,&mean,&err);
      //cout << " => n="<<n << endl << flush;
      if (n!=5) return -10;
      string sera = cera2;
      if (sera=="2022F"||sera=="2022G") continue;
      //md[cobs][cera][cmod] = make_pair(mean, err);
      string sobs(cobs);
      if (sobs=="Mu") continue;
      md[cobs][sera][cmod] = make_pair(mean, err);
    } // while getline
  } // for i!=nera
		

  // Graphical settings
  map<string, map<string, int> > mm; // markers
  mm["j"]["mu"] = kFullCircle;
  mm["j"]["npvall"] = kFullDiamond;
  mm["j"]["npvgood"] = kFullSquare;
  mm["j"]["rho"] = kFullCircle;
  //
  mm["j"]["Mu"] = kFullCircle;
  mm["j"]["NPV"] = kOpenSquare;//kFullDiamond;
  mm["j"]["NPVGood"] = kFullSquare;
  mm["j"]["RhoAll"] = kFullCircle;
  mm["j"]["Rho_C"] = kOpenCircle;//kFullStar;
  mm["j"]["Rho_CCPU"] = kOpenStar;
  
  map<string, int> mc; // colors
  mc["mu"] = kBlack;
  mc["npvall"] = kOrange+1;
  mc["npvgood"] = kRed;
  mc["rho"] = kBlue;
  //
  mc["Mu"] = kBlack;
  mc["NPV"] = kOrange+2;
  mc["NPVGood"] = kRed;
  mc["RhoAll"] = kBlue;
  mc["Rho_C"] = kMagenta+1;
  mc["Rho_CCPU"] = kCyan+2;//kOrange+1;//kYellow+2;

  //map<string, double> ms; // size-1
  
  
  int nobs = md.size();
  int nera = md["RhoAll"].size();
  cout << "Found " << nobs << " observables in dijet" << endl << flush;
  cout << "Found " << nera << " eras vs " << nera << " input" << endl << flush;

  TH1D *hnpv = new TH1D("hnpv",";#sigma_{MB} (mb);Observations",30,60,90);
  TH1D *hrho = new TH1D("hrho",";#sigma_{MB} (mb);Observations",30,60,90);
  TH1D *hrhoc = new TH1D("hrhoc",";#sigma_{MB} (mb);Observations",30,60,90);
  TH1D *hnpv_sel = new TH1D("hnpv_sel",";#sigma_{MB} (mb);N",30,60,90);
  TH1D *hrho_sel = new TH1D("hrho_sel",";#sigma_{MB} (mb);N",30,60,90);
  
  //TH1D *h = tdrHist("h","Data / MC",0.87,1.42,"Era",0.5,nera+0.5);
  TH1D *h = tdrHist("h","Data / MC",0.86,1.16,"Era",0.5,nera+0.5);
  //TH1D *h = tdrHist("h","Data / MC",0.65,1.55,"Era",0.5,nera+0.5);
  //if (drawPuoff) h->GetYaxis()->SetRangeUser(0.65,1.55);
  if (drawPuoff) h->GetYaxis()->SetRangeUser(0.65,1.30);
  // 2022: 34.748619392/fb 
  // 2023: 28.397263466/fb
  // 2024: 109.083612776/fb
  // Run3: 172.22950/fb
  lumi_136TeV = "Run 3, 172 fb^{-1}";//Prompt2024";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kRectangular);//kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0.5,1,0.5+nera,1);
  double sf = 1.06;
  l->SetLineStyle(kDotted);
  l->DrawLine(0.5,1.0+0.046,0.5+nera,1.0+0.046);
  l->DrawLine(0.5,1.0-0.046,0.5+nera,1.0-0.046);
  //l->SetLineColor(kRed);
  //double ksmb = 69.2/75.3;
  //l->DrawLine(0.5,ksmb,0.5+nera,ksmb);
  //l->DrawLine(0.5,sf,0.5+nera,sf);

  // TWiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Recommended_cross_section
  // https://arxiv.org/pdf/2011.14909 (Fig.2, Eq.6)
  TF1 *f1 = new TF1("f1","[0]+[1]*pow(log(x*x),[2])",0.01,100);
  //f1->SetParameters(28.84,0.05,2.37); // Fig.2
  f1->SetParameters(28.84,0.0458,2.374); // Eq.(6)
  double sigmamb8 = f1->Eval(8000.);
  double sigmamb13 = f1->Eval(13000.);
  double sigmamb136 = f1->Eval(13600.);
  cout << "sigmaMB(8 TeV)    = "<<sigmamb8 << endl;
  cout << "sigmaMB(13 TeV)   = "<<sigmamb13 << endl;
  cout << "sigmaMB(13.6 TeV) = "<<sigmamb136 << endl;
  double kt = sigmamb136/sigmamb8;

  //l->SetLineStyle(kSolid);
  //l->DrawLine(0.5,kt,0.5+nera,kt);
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035);
  //tex->DrawLatex(0.19,0.75,"Default #sigma_{MB}=69.2 mb");
  if (drawPuoff) {
    tex->DrawLatex(0.19,0.765,"#sigma_{MB} = 75.3 mb");
    tex->DrawLatex(0.19,0.72,"MARK: with weights");
    tex->DrawLatex(0.19,0.68,"HIST: no PU weights");
  }
  else {
    tex->DrawLatex(0.19,0.73,"#sigma_{MB} = 75.3#pm3.5 mb (4.6%)");
    tex->DrawLatex(0.19,0.68,"vs #sigma_{MB} = 69.2 mb (old default)");
    tex->DrawLatex(0.19,0.63,"vs #sigma_{MB} = 78.8 mb (theory)");
    //tex->DrawLatex(0.19,0.71,Form("Better #sigma_{MB}=%1.1f mb",69.2*sf));
    //tex->DrawLatex(0.19,0.67,Form("Best #sigma_{MB}=%1.1f mb",69.2*kt));
  }
 
  //TLegend *leg = tdrLeg(0.68,0.90-(nera+1)*0.05,0.93,0.90);
  TLegend *leg = tdrLeg(0.65,0.90-(nobs)*0.05,0.90,0.90);
  //leg->SetHeader("J");
  //TLegend *legg = tdrLeg(0.61,0.90-(nera+1)*0.05,0.86,0.90);
  //legg->SetHeader("#gamma");
  //TLegend *legz = tdrLeg(0.54,0.90-(nera+1)*0.05,0.79,0.90);
  //legz->SetHeader("Z");
  
  typedef map<string, map<string, map<string, pair<float,float> > > >::iterator IT;
  typedef map<string, map<string, pair<float,float> > >::iterator JT;

  int n(0), i24e(0), i24f(0), i24g(0);
  for (IT it = md.begin(); it != md.end(); ++it) {

    const char *cobs = it->first.c_str();
    cout << "Observable " << cobs << ": eras ";
    
    TH1D *hobs = new TH1D(Form("h_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);
    TH1D *hobs2 = new TH1D(Form("h2_%s",cobs),";Era;Data/MC",nera,0.5,nera+0.5);

    int j(0);
    for (JT jt = it->second.begin(); jt != it->second.end(); ++jt) {

      const char *cera = jt->first.c_str();
      cout << cera << ", ";

      //h->GetXaxis()->SetBinLabel(++j, jt->first.c_str());
      TString label(jt->first.c_str());
      label.ReplaceAll("20","");
      label.ReplaceAll("_nib1","");
      label.ReplaceAll("_nib2","");
      label.ReplaceAll("_nib3","");
      label.ReplaceAll("v","");

      if (label.Contains("24E")) i24e = j;
      if (label.Contains("24F")) i24f = j;
      if (label.Contains("24G")) i24g = j;
      h->GetXaxis()->SetBinLabel(++j, label.Data());
      h->GetXaxis()->SetLabelSize(0.05);//0.07);//0.05);

      //double mc = jt->second["mc"].first;
      //double mcpuoff = jt->second["mcpuoff"].first;
      double mc_withrw = jt->second["mcpuoff"].first;
      double mc_norw = jt->second["mc"].first;
      double dt = jt->second["data"].first;
      double edt = jt->second["data"].second;

      // Subtract hard scatter
      string sobs = string(cobs);
      if (useModObs &&
	  (sobs=="npvall" || sobs=="npvgood" ||
	   sobs=="NPV" || sobs=="NPVGood") &&
	  //mc!=0) {
	  (mc_withrw!=0 || mc_norw!=0)) {
	//mc -= 1;
	//mcpuoff -= 1;
	mc_withrw -= 1;
	mc_norw -= 1;
	dt -= 1;
      }
      if (useModObs &&
	  (sobs=="rho" ||
	   sobs=="RhoAll" || sobs=="Rho_C" || sobs=="Rho_CCPU") &&
	  //mc!=0) {
	  (mc_withrw!=0 || mc_norw!=0)) {
	//mc -= 2.0;
	//mcpuoff -= 2.0;
	mc_norw -= 2.0;
	mc_withrw -= 2.0;
	dt -= 2.0;
      }
      
      //hobs->SetBinContent(j, (mc!=0 ? dt/mc : 0));
      hobs->SetBinContent(j, (mc_withrw!=0 ? dt/mc_withrw : 0));
      //hobs->SetBinError(j, (mc!=0 ? edt/mc : 0));
      //hobs2->SetBinContent(j, (mcpuoff!=0 ? dt/mcpuoff : 0));
      hobs2->SetBinContent(j, (mc_norw!=0 ? dt/mc_norw : 0));


      //string sera = string(cera);
      string sera = string(label.Data());
      if ((sobs=="npvall" || sobs=="npvgood" ||
	   sobs=="NPV" || sobs=="NPVGood") &&
	  // mc!=0) {
	  mc_withrw!=0) {
	//hnpv->Fill(dt/mc*69.2);
	hnpv->Fill(dt/mc_withrw*75.3);
	if (sobs=="NPVGood" && (sera=="24E" || sera=="24F" || sera=="24G"))
	  //|| sera=="24H" || sera=="24I"))
	  hnpv_sel->Fill(dt/mc_withrw*75.3);
      }
      //if ((sobs=="rho" && sera!="D") && mc!=0) {
      if ((sobs=="rho" ||
	   sobs=="RhoAll") && // || sobs=="Rho_C" || sobs=="Rho_CCPU") &&
	  //mc!=0) {
	  mc_withrw!=0) {
	//hrho->Fill(dt/mc*69.2);
	hrho->Fill(dt/mc_withrw*75.3);
	if (sobs=="RhoAll" && (sera=="24E" || sera=="24F" || sera=="24G"))
	  //|| sera=="24H" || sera=="24I"))
	  hrho_sel->Fill(dt/mc_withrw*75.3);
      }
      if ((sobs=="Rho_C" || sobs=="Rho_CCPU") && mc_withrw!=0) {
	hrhoc->Fill(dt/mc_withrw*75.3);
      }

    } // for JT
    cout << endl;

    if (drawPuoff)
      //tdrDraw(hobs2,"HIST][",kFullCircle+n,kBlack+n,kSolid,-1,kNone,0,0.8);
      tdrDraw(hobs2,"HIST][",mm["j"][cobs],mc[cobs],kSolid,-1,kNone,0,0.8);
    //tdrDraw(hobs,"Pz",kFullCircle+n,kBlack+n,kSolid,-1,kNone,0,0.8); ++n;
    tdrDraw(hobs,"Pz",mm["j"][cobs],mc[cobs],kSolid,-1,kNone,0,0.8); ++n;
    leg->AddEntry(hobs,cobs,"PLE");
  } // for IT
  /*
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

    if (drawPuoff)
      //tdrDraw(hobs2,"HIST][",kOpenCircle,kBlack+ng,kDashed,-1,kNone,0,1.2);
      tdrDraw(hobs2,"HIST][",mm["g"][cobs],mc[cobs],kDashed,-1,kNone,0,0.8);
    //tdrDraw(hobs,"Pz",kOpenCircle+ng,kBlack+ng,kSolid,-1,kNone,0,1.2); ++ng;
    tdrDraw(hobs,"Pz",mm["g"][cobs],mc[cobs],kSolid,-1,kNone,0,1.2); ++ng;
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
      //double pumc = md["mu"][cera]["mcpuoff"].first;
      //double purw = md["mu"][cera]["mc"].first;
      double pumc = mg["mu"][cera]["mc"].first;
      double purw = mg["mu"][cera]["mcpuoff"].first;
      if (string(cera)=="G") {
	pumc = md["mu"][cera]["mcpuoff"].first;
	purw = md["mu"][cera]["mc"].first;
      }
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
      
      if (sera=="D") {
	cout << sera << "," << sobs << ": " << mcpuoff << " " << mc << " " << dt << " " << endl;
      }
    } // for JT

    if (drawPuoff)
      //tdrDraw(hobs2,"HIST][",kOpenCircle,kBlack+nz,kDotted,-1,kNone,0,1.2);
      tdrDraw(hobs2,"HIST][",mm["z"][cobs],mc[cobs],kDotted,-1,kNone,0,1.2);
    //tdrDraw(hobs,"Pz",kOpenDiamond+1+nz,kBlack+nz,kSolid,-1,kNone,0,1.2); ++nz;
    tdrDraw(hobs,"Pz",mm["z"][cobs],mc[cobs],kSolid,-1,kNone,0,1.2); ++nz;
    legz->AddEntry(hobs," ","PLE");
  } // for IT
  */
  //tex->DrawLatex(0.36,0.16,"E: HF scale up (rho up)");
  //tex->DrawLatex(0.56,0.21,"F: HCALDI on (rho up)");
  //tex->DrawLatex(0.56,0.19,"24E: HF scale fixed");
  tex->DrawLatex(0.72,0.15,"HF fixed");
  if (drawPuoff) {
    l->DrawLine(0.5+i24e,0.65,0.5+i24e,1.08);
    l->DrawLine(1.5+i24g,0.65,1.5+i24g,1.08);
  }
  else {
    l->DrawLine(0.5+i24e,0.86,0.5+i24e,1.055);
    l->DrawLine(1.5+i24g,0.86,1.5+i24g,1.055);
  }

  if (drawPuoff)
    c1->SaveAs("pdf/drawSigmaMB/drawSigmaMBV2_drawPUoff.pdf");
  else
    c1->SaveAs("pdf/drawSigmaMB/drawSigmaMBV2.pdf");


  //TH1D *h2 = tdrHist("h2","Observations",0,5.5,"#sigma_{MB} (mb)",65,85);
  //TH1D *h2 = tdrHist("h2","Observations",0,9.5,"#sigma_{MB} (mb)",65,85);
  //TH1D *h2 = tdrHist("h2","Observations",0,8.0,"#sigma_{MB} (mb)",55.5,89.5);
  //TH1D *h2 = tdrHist("h2","Observations",0,8.0,"#sigma_{MB} (mb)",61.5,83.5);
  //TH1D *h2 = tdrHist("h2","Observations",0,8.0,"#sigma_{MB} (mb)",64.5,80.5);
  TH1D *h2 = tdrHist("h2","Observations",0,13.5,"#sigma_{MB} (mb)",63.5,83.5);
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);

  l->SetLineStyle(kDashed);
  /*
  l->DrawLine(69.2,0,69.2,5.4);
  l->DrawLine(75.3,0,75.3,2.0);
  l->DrawLine(75.3,2.5,75.3,3.1);
  //l->DrawLine(75.3,2.45,75.3,5.0);
  l->DrawLine(78.8,0,78.8,5.0);
  */
  
  tdrDraw(hrhoc,"HIST",kFullTriangleUp,kCyan+2,kSolid,kCyan+2,1001,kCyan+1);
  hrhoc->SetFillColorAlpha(kCyan+1,0.50);
  
  tdrDraw(hnpv,"HIST",kFullSquare,kRed,kSolid,kRed+1,1001,kRed);
  hnpv->SetFillColorAlpha(kRed,0.50);
  tdrDraw(hnpv_sel,"HIST",kFullSquare,kRed+1,kSolid,kRed+1,1001,kRed);
  hnpv_sel->SetFillColorAlpha(kRed,0.50);

  tdrDraw(hrho,"HIST",kFullTriangleDown,kBlue,kSolid,kBlue+1,1001,kBlue);
  hrho->SetFillColorAlpha(kBlue,0.50);
  tdrDraw(hrho_sel,"HIST",kFullTriangleDown,kBlue,kSolid,kBlue+1,1001,kBlue);
  hrho_sel->SetFillColorAlpha(kBlue,0.50);

  TGraphErrors *gnpv = new TGraphErrors(1);
  //gnpv->SetPoint(0,hnpv->GetMean(),2.5);
  //gnpv->SetPointError(0,hnpv->GetRMS(),0);
  gnpv->SetPoint(0,hnpv_sel->GetMean(),2.5);
  gnpv->SetPointError(0,max(hnpv_sel->GetRMS(),0.5),0);
  TGraphErrors *grho = new TGraphErrors(1);
  //grho->SetPoint(0,hrho->GetMean(),2.5);
  //grho->SetPointError(0,hrho->GetRMS(),0);
  grho->SetPoint(0,hrho_sel->GetMean(),2.5);
  grho->SetPointError(0,max(hrho_sel->GetRMS(),0.5),0);
  TGraphErrors *gsmb = new TGraphErrors(1);
  //gsmb->SetPoint(0,0.5*(hrho->GetMean()+hnpv->GetMean()),2.25);
  //gsmb->SetPointError(0,0.5*(hrho->GetRMS()+hnpv->GetRMS()),0);
  gsmb->SetPoint(0,0.5*(hrho_sel->GetMean()+hnpv_sel->GetMean()),2.25);
  //gsmb->SetPointError(0,max(0.5*(hrho_sel->GetRMS()+hnpv_sel->GetRMS()),
  //			    0.5*fabs(hrho_sel->GetMean()-hnpv_sel->GetMean())),0);
  gsmb->SetPointError(0,0.5*fabs((grho->GetX()[0]+grho->GetEX()[0]) -
				 (gnpv->GetX()[0]-gnpv->GetEX()[0])),0);

  tdrDraw(gnpv,"Pz",kFullSquare,kRed+2); gnpv->SetLineWidth(2);
  tdrDraw(grho,"Pz",kFullTriangleDown,kBlue+1); grho->SetLineWidth(2);
  tdrDraw(gsmb,"Pz",kFullCircle,kBlack); gsmb->SetLineWidth(2);

  //tex->DrawLatex(0.55,0.41,Form("#sigma_{MB} = %1.1f#pm%1.1f mb",
  double relerr = gsmb->GetEX()[0] / gsmb->GetX()[0] * 100.;
  tex->DrawLatex(0.55,0.58,Form("#sigma_{MB} = %1.1f#pm%1.1f mb (%1.1f%%)",
				gsmb->GetX()[0],gsmb->GetEX()[0],relerr));
  //tex->DrawLatex(0.30,0.70,"#sigma_{MB} = 69.2 mb (old default)");
  tex->DrawLatex(0.185,0.65,"#sigma_{MB} = 69.2 mb");
  tex->DrawLatex(0.215,0.60,"(old default)");
  //tex->DrawLatex(0.54,0.65,"#sigma_{MB} = 75.3 mb (theory)");
  //tex->DrawLatex(0.59,0.66,"#sigma_{MB} = 78.8 mb (theory)");
  tex->DrawLatex(0.67,0.50,"#sigma_{MB} = 78.8 mb");
  tex->DrawLatex(0.71,0.45,"(theory)");

  l->DrawLine(69.2,0,69.2,7.5);
  l->DrawLine(75.3,0,75.3,1.8);
  l->DrawLine(75.3,2.7,75.3,7.2);
  l->DrawLine(78.8,0,78.8,5.0);

  TLegend *leg2 = tdrLeg(0.53,0.90-3*0.05,0.78,0.90);
  leg2->AddEntry(hnpv,"Vertices (N_{PV})","FPL");
  leg2->AddEntry(hrhoc,"Central offset (#rho_{C})","FPL");
  leg2->AddEntry(hrho,"Full offset (#rho)","FPL");
  
  gPad->RedrawAxis();
  c2->SaveAs("pdf/drawSigmaMB/drawSigmaMBV2_dist.pdf");
  
  return 0;
} // drawSigmaMB
