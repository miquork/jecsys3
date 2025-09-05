// Purpose: compare HcalDepths from PhotonJet, IsoTrack and IsoTrackOrig
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"

#include <vector>
#include <string>

#include "../tdrstyle_mod22.C" 

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
  vf[0] = new TFile("../isotrack/rootfiles/HcalDepthsFromIsoTrackOrig_25Sep04_EG25C2C.root","READ"); assert(vf[0]);
  vf[1] = new TFile("rootfiles/HcalDepthFromIsoTrack.root","READ");  assert(vf[1]);
  vf[2] = new TFile("rootfiles/HcalDepthFromPhotonJet.root","READ");  assert(vf[2]);

  // Data/MC difference fit in minitools/compareHcalDepths.C for IsoTrack
  // vs Lambda=((d-1)+0.5)*cosh(eta). For PhotonJet, Lambda+=1, scale 0.84
  TF1 *f3 = new TF1("f3","[p0]+[p1]*x+[p2]*x*x",0,8);
  //f3->SetParameters(1.238, -0.1903, 0.01734);
  f3->SetParameters(1.24, -0.2152, 0.02356);

  
  std::string depth_str[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "hfLong", "hfShort"};
  const char* depth_label[11] = {"ECAL","Depth 1","Depth 2","Depth 3","Depth 4",
				 "Depth 5 (HE)","Depth 6 (HE)","Depth 7 (HE)","HO","HF Long","HF Short"};
  int marker[nfile+2] = {kFullSquare, kOpenCircle, kFullDiamond, kOpenDiamond, kFullStar};
  int color[nfile+2] = {kBlue, kRed, kGreen+2, kCyan+1, kMagenta+1};
  string label[nfile+2] = {"IsoTrack (EGamma)","EGamma (IsoTrack)","EGamma (PhotonJet)","PhotonJet #times IsoTrack",
			   "PhotonJet #times f(#Lambda)"};
  
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

  TCanvas *c2 = new TCanvas("c2","c2",4*300,2*300);
  c2->Divide(4,2,0,0);

  extraText = "Private";
  lumi_136TeV = "2025C";
  TH1D *h3 = tdrHist("h3","Fraction ratio",0.5,1.5,"Depth (Layer or #Lambda)",0,8);
  TCanvas *c3 = tdrCanvas("c3",h3,8,11,kSquare);
  TGraphErrors *g3_ehb = new TGraphErrors();
  TGraphErrors *g3_hb = new TGraphErrors();
  TGraphErrors *g3_eh = new TGraphErrors();
  TGraphErrors *g3_h = new TGraphErrors();
  
  const int ndepth = 7+1;
  int ietamax(29);
  double eh_sf(0.84);
  for (int d = 1; d != ndepth+1; ++d) {

    c1->cd(d);
    
    TH1D *h = tdrHist(Form("h_%d",d),"Correction factor",0.3,2-1e-4,"i#eta",-ietamax,ietamax);
    if (d<=4) h->SetMinimum(0.5);
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
      double scale_photonjet = 0.84;
      double f2 = f3->Eval(lambda_photonjet) / scale_photonjet;

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


    c2->cd(d);
    
    TH1D *h2 = tdrHist(Form("h2_%d",d),"Ratio to IsoTrack",
		       0.5,1.6-1e-4,"i#eta",-ietamax,ietamax);
    if (d>=4) h->GetYaxis()->SetRangeUser(0.3,1.6);
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

      // Approximate depth profile difference    
      for (int i = 1; i != hdr->GetNbinsX()+1; ++i) {
	//double ieta = hdr->GetBinCenter(i);
	int ieta = int(hdr->GetBinCenter(i)+0.5);
	double eta = etaVal(ieta);
	double df = hdr->GetBinContent(i);
	if (fabs(eta)<1.305 && df!=0) {
	  double lambda0 = d; // HB
	  double lambda = ((d-1)+0.5)*cosh(eta); // HB
	  //if (ifile>=2) lambda += 1; // HB+ECAL
	  //if (ifile==1) lambda += 1; // shift, but why?
	  if (ifile==1) {
	    g3_h->SetPoint(g3_h->GetN(), lambda0, df);
	    g3_hb->SetPoint(g3_h->GetN(), lambda, df);
	    g3_hb->SetPointError(g3_h->GetN(), 0.05, 0.02);
	  }
	  if (ifile==2) {
	    g3_eh->SetPoint(g3_eh->GetN(), lambda0, df);
	    g3_ehb->SetPoint(g3_eh->GetN(), lambda+1, df*eh_sf);
	    g3_ehb->SetPointError(g3_h->GetN(), 0.05, 0.02);
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

  TMultiGraph *gm = new TMultiGraph();
  gm->Add(g3_hb);
  gm->Add(g3_ehb);
  //TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x+[3]/x+[4]/(x*x)",0,8);
  TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x",0,8);
  //f1->SetParameters(1.09,-0.11,+0.01,0,0);
  //f1->SetParameters(1.238, -0.1903, 0.01734);
  //f1->SetParameters(1.097, -0.1306, 0.01314);
  f1->SetParameters(1.24, -0.2152, 0.02356);

  gm->Fit(f1,"RN");

  f1->SetLineColor(kGreen+2);
  f1->SetLineWidth(3);
  f1->Draw("SAME");

  cout << "  // Data/MC difference fit in minitools/compareHcalDepths.C for IsoTrack" << endl;
  cout << "  // vs Lambda=((d-1)+0.5)*cosh(eta). For PhotonJet, Lambda+=1, scale 0.84" << endl;
  cout << "  TF1 *f3 = new TF1(\"f3\",\"" << f1->GetExpFormula() << "\",0,8);" << endl;
  cout << Form("  f3->SetParameters(%1.4g, %1.4g, %1.4g);",
	       f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2)) << endl;

  TGraph *g3 = new TGraph();
  g3->SetPoint(0, 0.5, 1/0.84);//1.16);
  tdrDraw(g3,"Pz",kFullStar,kGreen+3,kSolid,-1,kNone,0,2);
  
  TLegend *leg3 = tdrLeg(0.32,0.90-6*0.05,0.57,0.90);
  //leg3->SetTextSize(0.045);
  leg3->AddEntry(g3_h,"IsoTrack H-hadrons vs layer","PLE");
  leg3->AddEntry(g3_eh,"PhotonJet EH-had. vs layer","PLE");
  leg3->AddEntry(g3_hb,"IsoTrack H vs #Lambda","PLE");
  leg3->AddEntry(g3_ehb,Form("PhotonJet EH #times %1.2f vs #Lambda",eh_sf),"PLE");
  leg3->AddEntry(f1,"Quadratic fit to vs #Lambda","L");
  leg3->AddEntry(g3,Form("ECAL center: 0.5 #Lambda, f=1/%1.2f",eh_sf),"P");
  
  c1->SaveAs(Form("%s/compareHcalDepths.pdf",plot_dir.c_str()));
  c2->SaveAs(Form("%s/compareHcalDepthsRatios.pdf",plot_dir.c_str()));
  c3->SaveAs(Form("%s/compareHcalDepthsLambda.pdf",plot_dir.c_str()));
    
} // compareHcalDepths
