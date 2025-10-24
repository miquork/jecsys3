// Purpose: compare HcalDepths from PhotonJet, IsoTrack and IsoTrackOrig
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TF1.h"
#include "TMultiGraph.h"

#include <vector>
#include <string>

#include "../tdrstyle_mod22.C"

#include "hcal_depth.C"

// From Long Wang, Re: ROOT file with HCalRespCorrs for comparison, 2024/2/1
double etaVal(int ieta) {
  //if (ieta>=0) ++ieta; // patch iEta!=0 (Mikko's addition)
  double etavl(0.);
  if (ieta <= -24)
    etavl = .1695*ieta + 1.9931;
  else if (ieta <= -1)
    etavl = .0875*ieta + .0489;
  else if  (ieta < 24)
    etavl = .0875*ieta - .0489;
  else
    etavl = .1695*ieta - 1.9931;
  return etavl;
} // etaVal

void compareHcalDepths() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const int nfile(3);
  vector<TFile*> vf(nfile);
  // TBD: run EG25CC for more stats
  //vf[0] = new TFile("../isotrack/rootfiles/HcalDepthsFromIsoTrackOrig_25Sep04_EG25C2C.root","READ"); assert(vf[0]);
  vf[0] = new TFile("../isotrack/rootfiles/HcalDepthsFromIsoTrackOrig_25Sep04_EG25CC.root","READ"); assert(vf[0]);
  vf[1] = new TFile("rootfiles/HcalDepthFromIsoTrack.root","READ");  assert(vf[1]);
  vf[2] = new TFile("rootfiles/HcalDepthFromPhotonJet.root","READ");  assert(vf[2]);

  TFile *fout = new TFile("rootfiles/compareHcalDepths.root","RECREATE");
  curdir->cd();
  
  // Data/MC difference fit in minitools/compareHcalDepths.C for IsoTrack
  // vs Lambda=((d-1)+0.5)*cosh(eta). For PhotonJet, Lambda+=1, scale 0.84
  TF1 *f3 = new TF1("f3","[p0]+[p1]*x+[p2]*x*x",0,8);
  //f3->SetParameters(1.238, -0.1903, 0.01734);
  f3->SetParameters(1.24, -0.2152, 0.02356);
  TF1 *f4 = new TF1("f4","[p0]+[p1]*x+[p2]*x*x",0,8);
  f4->SetParameters(1.102, -0.0709, 0.0009969);

  const int nd(11);
  std::string depth_str[nd] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "hfLong", "hfShort"};
  const char* depth_label[nd] = {"ECAL","Depth 1","Depth 2","Depth 3","Depth 4",
				 "Depth 5 (HE)","Depth 6 (HE)","Depth 7 (HE)","HO","HF Long","HF Short"};
  int marker[nfile+2] = {kFullSquare, kOpenCircle, kFullDiamond, kOpenDiamond, kFullStar};
  int color[nfile+2] = {kBlue, kRed, kGreen+2, kCyan+1, kMagenta+1};
  string label[nfile+2] = {"IsoTrack (EGamma)","EGamma (IsoTrack)","EGamma (PhotonJet)","PhotonJet #times IsoTrack",
			   "PhotonJet #times f(#Lambda)"};

  // For c1f energy fraction plots
  int markerf[nfile+2] = {kNone, kOpenSquare, kOpenCircle, kFullCircle, kNone};
  int colorf[nfile+2] = {kNone, kRed, kBlue, kBlue, kNone};
  
  // Output directories
  std::string plot_dir = "pdf/compareHcalDepths/";
  gSystem->mkdir(plot_dir.c_str(), kTRUE);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed); l->SetLineColor(kGray+1);
  TLine *l2 = new TLine();
  l2->SetLineStyle(kDotted); l->SetLineColor(kGray);
  
  TCanvas *c1 = new TCanvas("c1","c1",4*300,2*300);
  c1->Divide(4,2,0,0);

  TCanvas *c1c = new TCanvas("c1c","c1c",4*300,2*300);
  c1c->Divide(4,2,0,0);
  
  TCanvas *c1cr = new TCanvas("c1cr","c1cr",4*300,2*300);
  c1cr->Divide(4,2,0,0);

  TCanvas *c1f = new TCanvas("c1f","c1f",4*300,2*300);
  c1f->Divide(4,2,0,0);

  TCanvas *c2 = new TCanvas("c2","c2",4*300,2*300);
  c2->Divide(4,2,0,0);
  
  TCanvas *c2c = new TCanvas("c2c","c2c",4*300,2*300);
  c2c->Divide(4,2,0,0);
  

  extraText = "Private";
  lumi_136TeV = "2025C";
  TH1D *h3 = tdrHist("h3","Fraction ratio (HB)",0.5,1.5,"Depth (Layer or #Lambda)",0,8);
  TCanvas *c3 = tdrCanvas("c3",h3,8,11,kSquare);
  TGraphErrors *g3_ehb = new TGraphErrors();
  TGraphErrors *g3_hb = new TGraphErrors();
  TGraphErrors *g3_eh = new TGraphErrors();
  TGraphErrors *g3_h = new TGraphErrors();

  TH1D *h4 = tdrHist("h4","Fraction ratio (HE)",0.2,2.0,"Depth (Layer or #Lambda)",0,8);
  TCanvas *c4 = tdrCanvas("c4",h4,8,11,kSquare);
  TGraphErrors *g4_ehb = new TGraphErrors();
  TGraphErrors *g4_hb = new TGraphErrors();
  TGraphErrors *g4_eh = new TGraphErrors();
  TGraphErrors *g4_h = new TGraphErrors();

  TH1D *h5 = tdrHist("h5","Fraction (HB)",0.0,0.7,"Depth (Layer or #Lambda)",0,8);
  TCanvas *c5 = tdrCanvas("c5",h5,8,11,kSquare);
  TGraphErrors *g5_ehb = new TGraphErrors();
  TGraphErrors *g5_hb = new TGraphErrors();
  TGraphErrors *g5_eh = new TGraphErrors();
  TGraphErrors *g5_h = new TGraphErrors();

  vector<TH1D*> vcorr(nd);
  
  const int ndepth = 7+1;
  //int ietamax(29);
  double ietamax(30-1e-4);
  double hb_sf(0.90);//0.84);
  double he_sf(1.00);//1.05);
  for (int d = 1; d != ndepth+1; ++d) {

    c1->cd(d);
    
    TH1D *h = tdrHist(Form("h_%d",d),"Correction factor",0.3,2-1e-4,"i#eta",-ietamax,ietamax);
    //TH1D *h = new TH1D(Form("h_%d",d),";i#eta;Correction factor",100,-ietamax,ietamax);
    h->GetYaxis()->SetRangeUser(0.0+1e-4,2.5-1e-4);
    //h->GetYaxis()->SetRangeUser(0.3,2-1e-4);
    //if (d<=4) h->SetMinimum(0.5);
    h->Draw("AXIS");
    l->DrawLine(-ietamax,1,ietamax,1);
    l2->DrawLine(-ietamax,1.2,ietamax,1.2);
    l2->DrawLine(-ietamax,0.8,ietamax,0.8);

    tex->SetTextSize(d<=4 ? 0.045*1.2 : 0.045);
    tex->DrawLatex(d%4==1 ? 0.45 : 0.35, 0.85, depth_label[d]);
    
    TH1D* vhd[nfile+2];
    for (unsigned int ifile = 0; ifile != nfile; ++ifile) {
      TH1D *hd = (TH1D*)vf[ifile]->Get(Form("h_depth_%d",d)); assert(hd);
      tdrDraw(hd,"Pz",marker[ifile],color[ifile],kSolid,-1,kNone,0,0.7);
      vhd[ifile] = hd;
    }

    // Solve for Corrected EGamma (PhotonJet)
    TH1D *hd = (TH1D*)vhd[2]->Clone(Form("h_depth_%d_corr",d));
    TH1D *hd3 = (TH1D*)vhd[2]->Clone(Form("h_depth_%d_corr3",d));
    for (int i = 1; i != hd->GetNbinsX()+1; ++i) {
      double c1 = vhd[0]->GetBinContent(i); // nominal IsoTrack correction
      double c1f1 = vhd[1]->GetBinContent(i); // correction times depth fraction for IsoTrack
      double f1 = (c1>0 ? c1f1 / c1 : 0); // depth fraction from IsoTrack

      // Patch missing mumbers outside tracver coverage
      int jeta = int(hd->GetBinCenter(i)+0.5);
      if (f1==0 && jeta<25) {
	int j0 = vhd[0]->GetXaxis()->FindBin(-25);
	int j1 = vhd[1]->GetXaxis()->FindBin(-25);
	double c1 = vhd[0]->GetBinContent(j0);
	double c1f1 = vhd[1]->GetBinContent(j1);
	f1 = (c1>0 ? c1f1 / c1 : 0); // depth fraction from IsoTrack
      }
      if (f1==0 && jeta>25) {
	int j0 = vhd[0]->GetXaxis()->FindBin(+25);
	int j1 = vhd[1]->GetXaxis()->FindBin(+25);
	double c1 = vhd[0]->GetBinContent(j0);
	double c1f1 = vhd[1]->GetBinContent(j1);
	f1 = (c1>0 ? c1f1 / c1 : 0); // depth fraction from IsoTrack
      }
      
      double c2f2 = vhd[2]->GetBinContent(i); // correction times depth fraction for PhotonJet
      double ec2f2 = vhd[2]->GetBinError(i);
      double c2 = (f1>0 ? c2f2 / f1 : c2f2);
      double ec2 = (f1>0 ? ec2f2 / f1 : ec2f2);
      hd->SetBinContent(i, c2);
      hd->SetBinError(i, c2, ec2);

      // Simplified models
      int ieta = int(hd->GetBinCenter(i)+0.5);
      double eta = etaVal(ieta);
      double lambda_isotrack = ((d-1)+0.5)*cosh(eta); // HB
      double lambda_photonjet = ((d-1)+0.5)*cosh(eta)+1; // HB + ECAL
      double scale_photonjet = hb_sf;//0.84;
      double f2 = f3->Eval(lambda_photonjet) / scale_photonjet;
      if (fabs(eta)>1.305) {
	// Need to add proper mapping of depths
	// ...
	lambda_isotrack = ((d-1)+0.5)*cosh(eta)/fabs(sinh(eta))+1; // HE + ECAL
	lambda_photonjet = ((d-1)+0.5)*cosh(eta)/fabs(sinh(eta))+1; // HE + ECAL
	scale_photonjet = he_sf;
	//f2 = f4->Eval(lambda_photonjet) / scale_photonjet;
	f2 = f1 / f4->Eval(lambda_isotrack) * (f4->Eval(lambda_photonjet) / scale_photonjet);
      }
      
      double c3 = (f2>0 ? c2f2 / f2 : c2f2);
      double ec3 = (f1>0 ? ec2f2 / f2 : ec2f2);
      hd3->SetBinContent(i, c3);
      hd3->SetBinError(i, c2, ec3);
    }
    tdrDraw(hd,"Pz",marker[nfile],color[nfile],kSolid,-1,kNone,0,0.7);
    vhd[nfile] = hd;

    tdrDraw(hd3,"Pz",marker[nfile+1],color[nfile+1],kSolid,-1,kNone,0,0.7);
    vhd[nfile+1] = hd3;
    
    if (d==ndepth) {
      TLegend *leg = tdrLeg(0.35,0.80-(nfile+2)*0.05,0.60,0.80);
      for (int ifile = 0; ifile != nfile+2; ++ifile) {
	leg->AddEntry(vhd[ifile],label[ifile].c_str(),"PLE");
      }
    }

    c1c->cd(d);

    TH1D *h1c = tdrHist(Form("h1c_%d",d),"Correction factor",0.0,3-1e-4,"i#eta",-ietamax,ietamax);
    h1c->Draw("AXIS");
    
    l->DrawLine(-ietamax,1,ietamax,1);
    l2->DrawLine(-ietamax,1.2,ietamax,1.2);
    l2->DrawLine(-ietamax,0.8,ietamax,0.8);
    //l->DrawLine(-16.5,0.3,-16.5,1.7);
    //l->DrawLine(+16.5,0.3,+16.5,1.7);
    if (d<5 || d>7) {
      l->DrawLine(-15.5,0.,-15.5,3.);
      l->DrawLine(+15.5,0.,+15.5,3.);
    }
    if (d<7) {
      l->DrawLine(-18.5,0.,-18.5,3.);
      l->DrawLine(+18.5,0.,+18.5,3.);
    }
    l2->DrawLine(-25.5,0.,-25.5,3.0);
    l2->DrawLine(+25.5,0.,+25.5,3.0);
    tex->DrawLatex(d%4==1 ? 0.45 : 0.35, 0.85, depth_label[d]);

    TH1D *h_isotrack = (TH1D*)vhd[0]->Clone(Form("h_isotrack_%d",d));
    TH1D *h_photonjet = (TH1D*)hd3->Clone(Form("h_photonjet_%d",d));
    //tdrDraw(h_isotrack,"Pz",kOpenSquare,kGray+1,kSolid,-1,kNone,0,0.7);
    //tdrDraw(h_photonjet,"Pz",kFullCircle,kBlack,kSolid,-1,kNone,0,0.7);
    tdrDraw(h_isotrack,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,0.7);
    tdrDraw(h_photonjet,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7);

    //if (d==ndepth) {
    double dxc1 = (d%4==1 ? 0.10 : 0);
    TLegend *legc1 = tdrLeg(0.35+dxc1,0.80-2*0.05,0.60+dxc1,0.80);
    legc1->AddEntry(h_isotrack,"IsoTrack","PLE");
    legc1->AddEntry(h_photonjet,"PhotonJet","PLE");
    //}

    vcorr[d] = h_photonjet;


    c1cr->cd(d);

    double ycrmin(0.45), ycrmax(1.55);
    TH1D *h1cr = tdrHist(Form("h1cr_%d",d),"Correction factor ratio",ycrmin+1e-4,ycrmax-1e-4,"i#eta",-ietamax,ietamax);
    h1cr->Draw("AXIS");
    
    l->DrawLine(-ietamax,1,ietamax,1);
    l2->DrawLine(-ietamax,1.05,ietamax,1.05);
    l2->DrawLine(-ietamax,0.95,ietamax,0.95);
    //l->DrawLine(-16.5,0.3,-16.5,1.7);
    //l->DrawLine(+16.5,0.3,+16.5,1.7);
    if (d<5 || d>7) {
      l->DrawLine(-15.5,ycrmin,-15.5,ycrmax);
      l->DrawLine(+15.5,ycrmin,+15.5,ycrmax);
    }
    if (d<7) {
      l->DrawLine(-18.5,ycrmin,-18.5,ycrmax);
      l->DrawLine(+18.5,ycrmin,+18.5,ycrmax);
    }
    l2->DrawLine(-25.5,ycrmin,-25.5,ycrmax);
    l2->DrawLine(+25.5,ycrmin,+25.5,ycrmax);
    tex->DrawLatex(d%4==1 ? 0.45 : 0.35, 0.85, depth_label[d]);

    TH1D *hr_isotrack = (TH1D*)h_isotrack->Clone(Form("hr_isotrack_%d",d));
    hr_isotrack->Divide(h_photonjet);
    tdrDraw(hr_isotrack,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7);

    TLegend *legc1r = tdrLeg(0.35+dxc1,0.80-2*0.05,0.60+dxc1,0.80);
    legc1r->AddEntry(hr_isotrack,"IsoTrack / PhotonJet","PLE");
    
    c1f->cd(d);
    
    TH1D *h1f = tdrHist(Form("h1f_%d",d),"Energy fraction",0.+1e-4,0.6-1e-4,"i#eta",-ietamax,ietamax);
    h1f->Draw("AXIS");

    if (d<5 || d>7) {
      l->DrawLine(-15.5,0.0,-15.5,0.6);
      l->DrawLine(+15.5,0.0,+15.5,0.6);
    }
    if (d!=7) {
      l->DrawLine(-18.5,0.0,-18.5,0.6);
      l->DrawLine(+18.5,0.0,+18.5,0.6);
    }
    l2->DrawLine(-25.5,0.0,-25.5,0.6);
    l2->DrawLine(+25.5,0.0,+25.5,0.6);
    
    int d2 = (d==8 ? 0 : d);
    tex->SetTextSize(d<=4 ? 0.045*1.2 : 0.045);
    tex->DrawLatex(d%4==1 ? 0.45 : 0.35, 0.85, depth_label[d2]);
    
    TH1D* vhf[nfile];
    TH1D* vhf2[nfile];
    for (unsigned int ifile = 0; ifile != nfile; ++ifile) {
      TH1D *hf = (TH1D*)vf[ifile]->Get(Form("h_frac_%d",d2)); //assert(hf);
      vhf[ifile] = vhf2[ifile] = 0;
      if (!hf) continue;
      tdrDraw(hf,"Pz",markerf[ifile],colorf[ifile],kSolid,-1,kNone,0,0.7);
      vhf[ifile] = hf;

      TH1D *hfe = (TH1D*)vf[ifile]->Get(Form("h_frac_%d",0)); assert(hfe);
      TH1D *hf2 = (TH1D*)hf->Clone(Form("hf2_%d_%d",ifile,d));
      for (int i = 0; i != hf2->GetNbinsX()+1; ++i) {
	double fh = hf->GetBinContent(i);
	double fe = hfe->GetBinContent(i);
	double fe2 = (ifile==2 ? 0.5 : fe);//0.25 + max(fe-0.25,0.)*2;
	//hf2->SetBinContent(i, fh / (1 - fe));
	hf2->SetBinContent(i, fh / (1 - fe2));
      }
      if (d2!=0 && ifile!=1)
	tdrDraw(hf2,"Pz",markerf[ifile+1],colorf[ifile+1],kSolid,-1,kNone,0,0.9);
      vhf2[ifile] = hf2;
    }

    if (d==4) {
      TLegend *leg = tdrLeg(0.35,0.80-3*0.05,0.60,0.80);
      leg->AddEntry(vhf[1],"IsoTrack","PLE");
      leg->AddEntry(vhf[2],"PhotonJet","PLE");
      leg->AddEntry(vhf2[2],"PhotonJet/(1-f_{E}=0.5)","PLE");
    }
    
    c2->cd(d);
    
    TH1D *h2 = tdrHist(Form("h2_%d",d),"Ratio to IsoTrack",
		       0.5,1.6-1e-4,"i#eta",-ietamax,ietamax);
    if (d>=4) h2->GetYaxis()->SetRangeUser(0.3,1.6);
    h2->Draw("AXIS");
    l->DrawLine(-ietamax,1,ietamax,1);
    l2->DrawLine(-ietamax,1.2,ietamax,1.2);
    l2->DrawLine(-ietamax,0.8,ietamax,0.8);

    tex->SetTextSize(d<=4 ? 0.045*1.2 : 0.045);
    tex->DrawLatex(d%4==1 ? 0.45 : 0.35, 0.85, depth_label[d]);

    TH1D* vhdr[nfile+1];
    for (unsigned int ifile = 0; ifile != nfile+2; ++ifile) {
      TH1D *hdr = (TH1D*)vhd[ifile]->Clone(Form("hdr_%d_%d",ifile,d));
      hdr->Divide(vhd[0]);
      tdrDraw(hdr,ifile==0 ? "HIST][" : "Pz",marker[ifile],color[ifile],kSolid,-1,kNone,0,0.7);
      vhdr[ifile] = hdr;

      TH1D *hdf = (ifile<nfile ? vhf2[ifile] : 0);

      // Approximate depth profile difference    
      for (int i = 1; i != hdr->GetNbinsX()+1; ++i) {
	//double ieta = hdr->GetBinCenter(i);
	int ieta = int(hdr->GetBinCenter(i)+0.5);
	double eta = etaVal(ieta);
	double df = hdr->GetBinContent(i);
	double ff = (hdf ? hdf->GetBinContent(i) : 0);

	//double width0 = (d==1 ? 1.5 : d==2 ? 4 : d==3 ? 5 : d==4 ? 7 : 1);
	double width0 = (d==1 ? 1.5 : d==2 ? 4 : d==3 ? 5 : d==4 ? 7 : 1);
	//double width0 = (d==1 ? 1.25 : d==2 ? 3.5 : d==3 ? 4.5 : d==4 ? 7.5 : 1);
	double width = width0*cosh(eta)/4; // depth==2 equivalent

	double thickness;
	double lambda = hcal_depth(ieta, d, thickness)-HCAL_ECAL_LAMBDA;
	//lambda -= thickness*0.5-0.5*cosh(eta); // center->front+0.5*lambda
	if (fabs(eta)<1.305) lambda += HB_D1_EXTRA_RADIAL*cosh(eta);
	double hcal_frac = (ifile==2 ? 2.3*(0.42-0.25)/(1-0.25) : 1.);
	
	if (fabs(eta)<1.305 && df!=0) {
	  double lambda0 = d; // HB
	  //double lambda = ((d-1)+0.5)*cosh(eta); // HB
	  //if (ifile>=2) lambda += 1; // HB+ECAL
	  //if (ifile==1) lambda += 1; // shift, but why?
	  if (ifile==1) {
	    g3_h->SetPoint(g3_h->GetN(), lambda0, df);
	    g3_hb->SetPoint(g3_hb->GetN(), lambda, df);
	    g3_hb->SetPointError(g3_hb->GetN()-1, 0.05, 0.02);

	    if (hdf) {
	      g5_h->SetPoint(g5_h->GetN(), lambda0, ff);
	      g5_hb->SetPoint(g5_hb->GetN(), lambda, ff / width * hcal_frac);
	      g5_hb->SetPointError(g5_hb->GetN()-1, 0.05, 0.02);
	    }
	  }
	  if (ifile==2) {
	    g3_eh->SetPoint(g3_eh->GetN(), lambda0, df);
	    g3_ehb->SetPoint(g3_eh->GetN(), lambda+HCAL_ECAL_LAMBDA, df*hb_sf);
	    g3_ehb->SetPointError(g3_h->GetN(), 0.05, 0.02);

	    if (hdf) {
	      g5_eh->SetPoint(g5_eh->GetN(), lambda0, ff);
	      g5_ehb->SetPoint(g5_ehb->GetN(), lambda+HCAL_ECAL_LAMBDA, ff / width * hcal_frac);//*hb_sf);
	      g5_ehb->SetPointError(g5_hb->GetN()-1, 0.05, 0.02);
	    }
	  }
	}
	if (fabs(eta)>=1.305 && df!=0) {
	  double lambda0 = d; // HE
	  //double lambda = ((d-1)+0.5)*cosh(eta)/fabs(sinh(eta)); // HE
	  if (ifile==1) {
	    g4_h->SetPoint(g4_h->GetN(), lambda0, df);
	    g4_hb->SetPoint(g4_hb->GetN(), lambda+HCAL_ECAL_LAMBDA, df);
	    g4_hb->SetPointError(g4_hb->GetN()-1, 0.05, 0.02);
	  }
	  if (ifile==2) {
	    g4_eh->SetPoint(g4_eh->GetN(), lambda0, df);
	    g4_ehb->SetPoint(g4_ehb->GetN(), lambda+HCAL_ECAL_LAMBDA, df*he_sf);
	    g4_ehb->SetPointError(g4_hb->GetN()-1, 0.05, 0.02);
	  }
	}
      } // for i
    } // for ifile

    if (d==ndepth) {
      TLegend *leg = tdrLeg(0.35,0.80-(nfile+2)*0.05,0.60,0.80);
      for (int ifile = 0; ifile != nfile+2; ++ifile) {
	leg->AddEntry(vhd[ifile],label[ifile].c_str(),"PLE");
      }
    }

  } // for d

  c3->cd();

  tdrDraw(g3_eh, "Pz", kFullSquare, kBlue, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g3_h, "Pz", kFullCircle, kRed, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g3_hb, "Pz", kOpenCircle, kRed, kSolid, -1);
  tdrDraw(g3_ehb, "Pz", kOpenSquare, kBlue, kSolid, -1);

  TMultiGraph *g3m = new TMultiGraph();
  g3m->Add(g3_hb);
  g3m->Add(g3_ehb);
  //TF1 *f3 = new TF1("f3","[0]+[1]*x+[2]*x*x+[3]/x+[4]/(x*x)",0,8);
  //TF1 *f3 = new TF1("f3","[0]+[1]*x+[2]*x*x",0,8);
  //f3->SetParameters(1.09,-0.11,+0.01,0,0);
  //f3->SetParameters(1.238, -0.1903, 0.01734);
  //f3->SetParameters(1.097, -0.1306, 0.01314);
  f3->SetParameters(1.24, -0.2152, 0.02356);

  g3m->Fit(f3,"RN");

  f3->SetLineColor(kGreen+2);
  f3->SetLineWidth(3);
  f3->Draw("SAME");

  cout << "  // Data/MC difference fit in minitools/compareHcalDepths.C for IsoTrack" << endl;
  cout << "  // vs Lambda=((d-1)+0.5)*cosh(eta). For PhotonJet, Lambda+=1, scale 0.84" << endl;
  cout << "  TF1 *f3 = new TF1(\"f3\",\"" << f3->GetExpFormula() << "\",0,8);" << endl;
  cout << Form("  f3->SetParameters(%1.4g, %1.4g, %1.4g);",
	       f3->GetParameter(0), f3->GetParameter(1), f3->GetParameter(2)) << endl;

  TGraph *g3 = new TGraph();
  g3->SetPoint(0, 0.5, 1/hb_sf);//1.16);
  tdrDraw(g3,"Pz",kFullStar,kGreen+3,kSolid,-1,kNone,0,2);
  
  TLegend *leg3 = tdrLeg(0.32,0.90-6*0.05,0.57,0.90);
  //leg3->SetTextSize(0.045);
  leg3->AddEntry(g3_h,"IsoTrack H-hadrons vs layer","PLE");
  leg3->AddEntry(g3_eh,"PhotonJet EH-had. vs layer","PLE");
  leg3->AddEntry(g3_hb,"IsoTrack H vs #Lambda","PLE");
  leg3->AddEntry(g3_ehb,Form("PhotonJet EH #times %1.2f vs #Lambda",hb_sf),"PLE");
  leg3->AddEntry(f3,"Quadratic fit to vs #Lambda","L");
  leg3->AddEntry(g3,Form("ECAL center: 0.5 #Lambda, f=1/%1.2f",hb_sf),"P");


  c4->cd();
  
  tdrDraw(g4_eh, "Pz", kFullSquare, kBlue, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g4_h, "Pz", kFullCircle, kRed, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g4_hb, "Pz", kOpenCircle, kRed, kSolid, -1);
  tdrDraw(g4_ehb, "Pz", kOpenSquare, kBlue, kSolid, -1);

  TMultiGraph *g4m = new TMultiGraph();
  g4m->Add(g4_hb);
  g4m->Add(g4_ehb);
  //TF1 *f4 = new TF1("f4","[0]+[1]*x+[2]*x*x",0,8);
  f4->SetParameters(1.24, -0.2152, 0.02356);

  g4m->Fit(f4,"RN");

  f4->SetLineColor(kGreen+2);
  f4->SetLineWidth(3);
  f4->Draw("SAME");

  cout << "  // Data/MC difference fit in minitools/compareHcalDepths.C for IsoTrack" << endl;
  cout << "  // vs Lambda=((d-1)+0.5)*cosh(eta). For PhotonJet, Lambda+=1, scale 0.84" << endl;
  cout << "  TF1 *f4 = new TF1(\"f4\",\"" << f4->GetExpFormula() << "\",0,8);" << endl;
  cout << Form("  f4->SetParameters(%1.4g, %1.4g, %1.4g);",
	       f4->GetParameter(0), f4->GetParameter(1), f4->GetParameter(2)) << endl;

  TGraph *g4 = new TGraph();
  g4->SetPoint(0, 0.5, 1/he_sf);//1.16);
  tdrDraw(g4,"Pz",kFullStar,kGreen+3,kSolid,-1,kNone,0,2);
  
  TLegend *leg4 = tdrLeg(0.32,0.90-6*0.05,0.57,0.90);
  leg4->AddEntry(g4_h,"IsoTrack H-hadrons vs layer","PLE");
  leg4->AddEntry(g4_eh,"PhotonJet EH-had. vs layer","PLE");
  leg4->AddEntry(g4_hb,"IsoTrack H vs #Lambda","PLE");
  leg4->AddEntry(g4_ehb,Form("PhotonJet EH #times %1.2f vs #Lambda",he_sf),"PLE");
  leg4->AddEntry(f4,"Quadratic fit to vs #Lambda","L");
  leg4->AddEntry(g4,Form("ECAL center: 0.5 #Lambda, f=1/%1.2f",he_sf),"P");


  c5->cd();
  
  tdrDraw(g5_eh, "Pz", kFullSquare, kBlue, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g5_h, "Pz", kFullCircle, kRed, kSolid, -1, kNone, 0, 0.7);
  tdrDraw(g5_hb, "Pz", kOpenCircle, kRed, kSolid, -1);
  tdrDraw(g5_ehb, "Pz", kOpenSquare, kBlue, kSolid, -1);

  TLegend *leg5 = tdrLeg(0.32,0.90-4*0.05,0.57,0.90);
  leg5->AddEntry(g5_h,"IsoTrack H-hadrons vs layer","PLE");
  leg5->AddEntry(g5_eh,"PhotonJet EH-had. vs layer","PLE");
  leg5->AddEntry(g5_hb,"IsoTrack H vs #Lambda","PLE");
  leg5->AddEntry(g5_ehb,Form("PhotonJet EH #times %1.2f vs #Lambda",he_sf),"PLE");
  
  c1->SaveAs(Form("%s/compareHcalDepths.pdf",plot_dir.c_str()));
  c1c->SaveAs(Form("%s/compareHcalDepths_Final.pdf",plot_dir.c_str()));
  c1cr->SaveAs(Form("%s/compareHcalDepths_FinalRatio.pdf",plot_dir.c_str()));
  c1f->SaveAs(Form("%s/compareHcalDepthsFractions.pdf",plot_dir.c_str()));
  c2->SaveAs(Form("%s/compareHcalDepthsRatios.pdf",plot_dir.c_str()));
  c3->SaveAs(Form("%s/compareHcalDepthsLambda_HB.pdf",plot_dir.c_str()));
  c4->SaveAs(Form("%s/compareHcalDepthsLambda_HE.pdf",plot_dir.c_str()));
  c5->SaveAs(Form("%s/compareHcalDepthsFractionVsLambda_HB.pdf",plot_dir.c_str()));

  // Load missing HF depths from PhotonJet
  TH1D *hecal = (TH1D*)vf[2]->Get(Form("h_depth_%d",0)); assert(hecal);
  vcorr[0] = hecal;
  for (int d = ndepth+1; d != nd; ++d) {
    TH1D *h = (TH1D*)vf[2]->Get(Form("h_depth_%d",d)); assert(h);
    vcorr[d] = h;
  }

  fout->cd();
  for (int d = 0; d != nd; ++d) {
    TH1D *h = vcorr[d]; assert(h);
    h->SetTitle(depth_label[d]);
    h->Write(Form("h_depth_%d",d));
  }
  fout->Write();
  fout->Close();
    
} // compareHcalDepths
