// Purpose: Draw HCAL response corrections for 19Dec vs 22Sep to compare them
#include "TFile.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"

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

// the bins are set through this way: 
// hist = new TH2F(Form("hnew_depth%d",d), Form("hnew_depth%d; iEta; iPhi",d), 86, -43, 43, 72, 0, 72);
// hist.SetBinContent( (ieta_<0?ieta_+44:ieta_+43), iphi_, respcorr_);

void drawHcalRespCorrs(string run);

int _iera(0);
TH1D *_hera(0);
TFile *_fout(0);
void drawHcalRespCorr() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  _fout = new TFile("rootfiles/hcalRespCorr_1D.root","UPDATE");
  
  _iera = 1;
  _hera = new TH1D("_hera",";Era;Avg. HcalRespCorr",7,0.5,7.5);
  _hera->GetYaxis()->SetRangeUser(0.85,1.20);
  /*
  drawHcalRespCorrs("2022CD");
  //drawHcalRespCorrs("2022CD_vs_22Sep");
  //drawHcalRespCorrs("2022CD_vs_prompt");
  drawHcalRespCorrs("2022E");
  //drawHcalRespCorrs("2022E_vs_22Sep");
  //drawHcalRespCorrs("2022E_vs_prompt");
  drawHcalRespCorrs("2022F");
  drawHcalRespCorrs("2022G");
  drawHcalRespCorrs("2023Cv123");
  drawHcalRespCorrs("2023Cv4");
  drawHcalRespCorrs("2023D");
  */
  drawHcalRespCorrs("ECAL_R_HCAL_DI");
  
  extraText = "Private";
  lumi_136TeV = "Run3";
  TCanvas *c1 = tdrCanvas("c1",_hera,8,11,kSquare);
  TLine *l = new TLine();

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0.5,1,7.5,1);

  l->SetLineStyle(kDotted);
  l->SetLineColor(kGray+2);
  l->DrawLine(0.5,1.015,3.5,1.015);
  l->DrawLine(3.5,0.940,7.5,0.940);
  
  tdrDraw(_hera,"HE",kFullCircle,kBlack,kSolid,-1,kNone);

  gPad->RedrawAxis();
  
  c1->SaveAs("pdf/drawHcalRespCorr/drawHcalRespCorr_eras.pdf");
} // drawHcalRespCorr

void scale(TH2D *h2, int id) {

  //double weight[4] = {0.2,0.3,0.3,0.2}; // must sum to 1.000
  // Hui Wang, PdmV General, 7 Nov 2023: https://indico.cern.ch/event/1337844/contributions/5631869/attachments/2747658/4781510/HCAL%20timing%20in%20Run%203.pdf
  // Page 7, data (HO=depth8)
  double weight_hb[8] =    {0.17, 0.37, 0.27, 0.16,  0, 0, 0,  0.03}; // 1.000

  // Guessing other regions
  // depth5: 17-27
  // depth6: 18-27
  // depth7: 25-27
  double weight_ec5[8] =   {0.15, 0.35, 0.24, 0.14,  0.12, 0.00,  0, 0.};
  double weight_ec56[8] =  {0.14, 0.32, 0.22, 0.13,  0.12, 0.07,  0, 0.};
  double weight_ec567[8] = {0.13, 0.31, 0.21, 0.12,  0.11, 0.07,  0.05, 0.};
  double weight_hf[8] =    {0.17, 0.38, 0.28, 0.17,  0, 0, 0,  0.00};
  
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {

    int ieta = (i<=43 ? i-44 : i-43);
    //double eta = etaVal(ieta);
    double *w(0);      
    if (abs(ieta)<=16)      w = weight_hb;
    else if (abs(ieta)<=17) w = weight_ec5;
    else if (abs(ieta)<=24) w = weight_ec56;
    else if (abs(ieta)<=30) w = weight_ec567;
    else if (abs(ieta)<=43) w = weight_hf;

    // count empty bins
    int nbins(0);
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)!=0) ++nbins;
    } // for j
    int ntot = h2->GetNbinsY();
    
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      h2->SetBinContent(i, j, h2->GetBinContent(i, j) * w[id]
			* (nbins!=0 ? 72./double(nbins) : 0));
    } // for j
  } // for i

} // scale

void drawHcalRespCorrs(string run) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TFile *f(0);
  // ECALRATIO: HcalRespCorrs_2023_v3.0_data = old
  // ECAL_R_HCAL_DI: HcalRespCorrs_2024_v2.0_data = new
  if (run=="ECAL_R_HCAL_DI")
    f = new TFile("rootfiles/hcalrespcorrs_2023-v3_2024-v2.root","READ");
  else
    f = new TFile(Form("rootfiles/hcalrespcorr_%s.root",cr),"READ");
  assert(f && !f->IsZombie());

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);

  // create new histogram for mapping iEta->detector eta
  // Regular L2Relative eta binning
  double vx[] =
      {-5.191,
       -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
       -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
       -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
       -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435,
       -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
       0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
       1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
       2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
       4.363, 4.538, 4.716, 4.889, 5.191};
  const int nx = sizeof(vx)/sizeof(vx[0]) - 1;
  //cout << "Nbins vx = " << nx << endl << flush;

  TH1D *hetan = new TH1D("hetan",";#eta;HcalRespCorr avg.",nx,vx);
  TProfile *petan = new TProfile("petan",";#eta;HcalRespCorr avg.",nx,vx);
  TH1D *hetao = new TH1D("hetao",";#eta;HcalRespCorr avg.",nx,vx);
  TProfile *petao = new TProfile("petao",";#eta;HcalRespCorr avg.",nx,vx);
  
  
  //TH1D *h = tdrHist("h","HcalRespCorr",34,38,"iEta",-42,42);
  //TH1D *hd = tdrHist("hd","19Dec/22Sep",0.8,1.2,"iEta",-42,42);
  //TH1D *h = tdrHist("h","HcalRespCorr",0.+1e-4,3.-1e-4,"iEta",-43,43);
  //TH1D *hd = tdrHist("hd","19Dec/22Sep",0.50,2.-1e-4,"iEta",-43,43);
  // ECAL_R_HCAL_DI
  TH1D *h = tdrHist("h","HcalRespCorr (DD)",0.+1e-4,3.-1e-4,"iEta",-43,43);
  TH1D *hd = tdrHist("hd","DI / DD",0.50,2.-1e-4,"iEta",-43,43);
  lumi_136TeV = cr;//"2022D";
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cr),h,hd,8,11);
  
  int color[8] = {kRed, kOrange+2, kMagenta+2, kBlue, // barrel 
		  kRed-9, kBlue-9, kGreen-9, // endcap
		  kGreen+2}; // HO

  c1->cd(1);
  l->DrawLine(-43,1,43,1);

  TLegend *leg = tdrLeg(0.40,0.88,0.65,0.88);
  TLegend *leg2 = tdrLeg(0.73,0.88,0.98,0.88);
  
  c1->cd(2);
  l->DrawLine(-43,1,43,1);

  TH2D *h2sumn(0), *h2sumo(0);

  const int nd = 8; // number of depths to analyze
  TH1D *vhd[nd];
  for (int id = 0; id != nd; ++id) {

    TH2D *h2n = (TH2D*)f->Get(Form("hnew_depth%d",id+1)); //assert(h2n);
    TH2D *h2o = (TH2D*)f->Get(Form("hold_depth%d",id+1)); //assert(h2o);
    if (!h2n && !h2o) {
    h2n = (TH2D*)f->Get(Form("hnew_capid0_depth%d",id+1));
    h2o = (TH2D*)f->Get(Form("hold_capid0_depth%d",id+1));
    }
    assert(h2n);
    assert(h2o);
    
    // Calculate average fractional contribution to jet energy
    TH2D *h2nw = (TH2D*)h2n->Clone(Form("h2nw_depth%d",id+1));
    TH2D *h2ow = (TH2D*)h2o->Clone(Form("h2ow_depth%d",id+1));
    scale(h2nw,id);
    scale(h2ow,id);

    // Calculate weighted average
    if (!h2sumn && !h2sumo) {
      h2sumn = (TH2D*)h2nw->Clone("h2sumn"); //h2sumn->Scale(weight[id]);
      h2sumo = (TH2D*)h2ow->Clone("h2sumo"); //h2sumo->Scale(weight[id]);
    }
    else {
      h2sumn->Add(h2nw);//, weight[id]);
      h2sumo->Add(h2ow);//, weight[id]);
    }
    
    //TProfile *pn = h2n->ProfileX(Form("pn_%d",id+1));
    //TProfile *po = h2o->ProfileX(Form("po_%d",id+1));
    TH1D *hn = h2n->ProjectionX(Form("hn_%d",id+1));
    TH1D *ho = h2o->ProjectionX(Form("ho_%d",id+1));
    //hn->Scale(1./72.);
    //ho->Scale(1./72.);
    
    // Divide by number of non-empty bins for each layer
    for (int i = 1; i != h2n->GetNbinsX()+1; ++i) {
      int nbins(0);
      for (int j = 1; j != h2n->GetNbinsY()+1; ++j) {
	if (h2n->GetBinContent(i,j)!=0) ++nbins;
      }
      if (nbins!=0) {
	hn->SetBinContent(i, hn->GetBinContent(i)/nbins);
	ho->SetBinContent(i, ho->GetBinContent(i)/nbins);
      }
    }
    
    
    c1->cd(1);
    //tdrDraw(pn,"HIST",kFullCircle,color[id],kSolid,-1,kNone);
    //tdrDraw(po,"HIST",kOpenCircle,color[id],kSolid,-1,kNone);
    //if (id+1<=4 || id+1==8) {
    //tdrDraw(ho,"HIST",kOpenCircle,color[id],kDashed,-1,kNone);
    tdrDraw(ho,"HIST",kFullDiamond,color[id],kSolid,-1,kNone); //ECAL_R_HCAL_DI
    if (run!="ECAL_R_HCAL_DI") {
    if (id+1==4) 
      tdrDraw(hn,"HP",kFullCircle,color[id],kSolid,-1,kNone);
    else {
      tdrDraw(hn,"HIST",kFullCircle,color[id],kSolid,-1,kNone);
      tdrDraw(hn,"HISTP",kFullCircle,color[id],kSolid,-1,kNone);
    }
    hn->SetMarkerSize(0.6); // ECAL_R_HCAL_DI
    //hn->SetMarkerSize(0.6);
    //ho->SetLineWidth(2); // ECAL_R_HCAL_DI

    if (id+1<=4 || id+1==8) {
      leg->SetY1(leg->GetY1()-0.045);
      leg->AddEntry(hn,Form("Depth %d",id+1),"PLE");
    }
    else {
      leg2->SetY1(leg2->GetY1()-0.045);
      leg2->AddEntry(hn,Form("Depth %d",id+1),"PLE");
    }
    }
    else {
      if (id+1<=4 || id+1==8) {
      leg->SetY1(leg->GetY1()-0.045);
      leg->AddEntry(ho,Form("Depth %d",id+1),"PLE");
      }
      else {
	leg2->SetY1(leg2->GetY1()-0.045);
	leg2->AddEntry(ho,Form("Depth %d",id+1),"PLE");
      }
    }
      
      
    c1->cd(2);
    //TH1D *hr = pn->ProjectionX(Form("hr_%d",id+1));
    //hr->Divide(po);
    TH1D *hr = (TH1D*)hn->Clone(Form("hr_%d",id+1));
    hr->Divide(ho);

    if (id+1<=4 || id+1==8) {
      tdrDraw(hr,"HIST",kFullCircle,color[id],kSolid,-1,kNone);
      
      TF1 *f1 = new TF1(Form("f1_%d",id+1),"[0]",-1.305/0.087,+1.305/0.087);
      f1->SetLineColor(hr->GetLineColor());
      hr->Fit(f1,"QRN");
      f1->Draw("SAME");
      if (id==0)
	cout << endl << "Run" << cr << endl;
      cout << Form("%s: %1.3f +/- %1.3f\n",f1->GetName(),f1->GetParameter(0),
		   f1->GetParError(0));
    }

    vhd[id] = ho;
  } // for id

  //cout << "Nbins hn = " << h2sumn->GetNbinsX() << endl << flush;
  
  // Draw weighted average
  c1->cd(1);

  TH1D *hn = h2sumn->ProjectionX("hn");
  TH1D *ho = h2sumo->ProjectionX("ho");
  hn->Scale(1./72.);
  ho->Scale(1./72.);
    
  //tdrDraw(ho,"HIST",kOpenDiamond,kBlack,kDashed,-1,kNone);
  tdrDraw(ho,"HISTP",kFullDiamond,kBlack,kSolid,-1,kNone); // ECAL_R_HCAL_DI
  //tdrDraw(hn,"HP",kFullDiamond,kBlack,kSolid,-1,kNone);
  //tdrDraw(hn,"HISTP",kFullDiamond,kBlack,kSolid,-1,kNone); // ECAL_R_HCAL_DI
  //hn->SetMarkerSize(0.6); // ECAL_R_HCAL_DI
  //ho->SetLineWidth(2); // ECAL_R_HCAL_DI

  leg->SetY1(leg->GetY1()-0.045);
  //leg->AddEntry(hn,"Weighted avg.","PLE");
  leg->AddEntry(ho,"Weighted avg.","PLE");

  gPad->Update();
  gPad->RedrawAxis();

  c1->cd(2);

  TH1D *hr = (TH1D*)hn->Clone("hr");
  hr->Divide(ho);

  //tdrDraw(hr,"HIST",kFullDiamond,kBlack,kSolid,-1,kNone);
  tdrDraw(hr,"HISTP",kFullDiamond,kBlack,kSolid,-1,kNone);

  TF1 *f1 = new TF1("f1_avg","[0]",-1.305/0.087,+1.305/0.087);
  f1->SetLineColor(hr->GetLineColor());
  hr->Fit(f1,"QRN");
  f1->Draw("SAME");
  cout << Form("%s: %1.3f +/- %1.3f\n",f1->GetName(),f1->GetParameter(0),
	       f1->GetParError(0));
  
  gPad->Update();
  gPad->RedrawAxis();

  c1->SaveAs(Form("pdf/drawHcalRespCorr/drawHcalRespCorr_%s.pdf",cr));


  TH1D *h2 = tdrHist("h2","HcalRespCorr",0.3,1.6,"#eta",-5.2,5.2);
  //TH1D *h2d = tdrHist("h2d","19Dec/22Sep",0.7,1.3,"#eta",-5.2,5.2);
  //TH1D *h2d = tdrHist("h2d","19Dec/22Sep",0.55,1.30,"#eta",-5.2,5.2);
  // ECAL_R_HCAL_DI
  TH1D *h2d = tdrHist("h2d","DI / DD",0.55,1.30,"#eta",-5.2,5.2);
  lumi_136TeV = cr;//"2022D";
  extraText = "Private";
  TCanvas *c2 = tdrDiCanvas(Form("c2_%s",cr),h2,h2d,8,11);

  for (int i = 1; i != hn->GetNbinsX()+1; ++i) {
    //int ieta = floor(hn->GetBinLowEdge(i)+0.5); // incl. ieta=0, fix in etaVal
    //hist.SetBinContent( (ieta_<0?ieta_+44:ieta_+43), iphi_, respcorr_);
    int ieta = (i<=43 ? i-44 : i-43);
    double eta = etaVal(ieta);
    hetan->Fill(eta, hn->GetBinContent(i));
    petan->Fill(eta, hn->GetBinContent(i));
    hetao->Fill(eta, ho->GetBinContent(i));
    petao->Fill(eta, ho->GetBinContent(i));
  }
  
  c2->cd(1);
  l->DrawLine(-5.2,1,5.2,1);

  tdrDraw(petao,"HIST",kOpenCircle,kBlack,kDashed,-1,kNone);
  //tdrDraw(petao,"HISTP",kOpenCircle,kBlack,kSolid,-1,kNone);
  //tdrDraw(hetao,"HISTP",kOpenCircle,kBlack,kSolid,-1,kNone);
  //petao->SetMarkerSize(0.6);

  tdrDraw(petan,"HIST",kFullCircle,kBlack,kSolid,-1,kNone);
  tdrDraw(petan,"HISTP",kFullCircle,kBlack,kSolid,-1,kNone);
  //tdrDraw(hetan,"HISTP",kOpenCircle,kBlack,kSolid,-1,kNone);
  petan->SetMarkerSize(0.6);


  c2->cd(2);
  l->DrawLine(-5.2,1,5.2,1);

  TH1D *heta = petan->ProjectionX("heta");
  heta->Divide(petao);

  tdrDraw(heta,"HIST",kFullCircle,kBlack,kSolid,-1,kNone);
  tdrDraw(heta,"HISTP",kFullCircle,kBlack,kSolid,-1,kNone);
  heta->SetMarkerSize(0.6);

  TF1 *f1cmb = new TF1(Form("f1_cmb_%s",cr),"[0]",-1.305,+1.305);
  f1cmb->SetLineColor(heta->GetLineColor());
  heta->Fit(f1cmb,"QRNW");
  f1cmb->Draw("SAME");
  cout << Form("%s: %1.3f +/- %1.3f\n",f1cmb->GetName(),f1cmb->GetParameter(0),
	       f1cmb->GetParError(0));

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(2*0.045);
  tex->DrawLatex(0.44,0.32,Form("%1.3f #pm %1.3f",f1->GetParameter(0),
				f1->GetParError(0)));
  
  c2->SaveAs(Form("pdf/drawHcalRespCorr/drawHcalRespCorr_avg_%s.pdf",cr));

  if (_hera && _iera>0 && _iera<_hera->GetNbinsX()+1) {
    _hera->SetBinContent(_iera, f1cmb->GetParameter(0));
    if (!(run=="ECAL_R_HCAL_DI"))
      _hera->SetBinError(_iera, f1cmb->GetParError(0));
    _hera->GetXaxis()->SetBinLabel(_iera, cr);
    ++_iera;
  }

  if (_fout) {
    _fout->cd();
    for (int i = 0; i != nd; ++i) {
      if (vhd[i]) {
	vhd[i]->SetLineStyle(kSolid);
	vhd[i]->Write(Form("depth%d_%s",i+1,cr), TObject::kOverwrite);
      }
    }
    ho->Write(Form("avg_%s",cr), TObject::kOverwrite);
    curdir->cd();
  }
  
} // drawHcalRespCorr§
