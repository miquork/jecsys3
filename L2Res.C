// Purpose: derive L2Res from Z+jet, gamma+jet, dijet based on new TProfile2D
#include "TFile.h"
#include "TProfile2D.h"
#include "TF1.h"

#include "tdrstyle_mod22.C"

TLegend *_leg(0);

// Step 1. Slice 1D profile out of 2D in given range and draw it
TProfile* drawEta(TProfile2D *p2, double ptmin, double ptmax,
		  string draw, int marker, int color, string label="") {
  
  assert(p2);
  int iy1 = p2->GetYaxis()->FindBin(ptmin);
  int iy2 = p2->GetYaxis()->FindBin(ptmax);
  string id = Form("p%s_%1.0f_%1.0f_%s_%d_%d", p2->GetName(), ptmin, ptmax,
		   draw.c_str(), marker, color);
  TProfile *p = p2->ProfileX(Form("p%s",id.c_str()),iy1,iy2);
  tdrDraw(p,draw,marker,color,kSolid,-1,kNone);

  if (_leg && label!="") {
    int pt1 = int(p2->GetYaxis()->GetBinLowEdge(iy1));
    int pt2 = int(p2->GetYaxis()->GetBinLowEdge(iy2+1));
    if (label == "Z")
      _leg->AddEntry(p, Form("%s+jet [%d,%d]",label.c_str(),pt1,pt2),"PLE");
    else
      _leg->AddEntry(p, Form("%s+jet",label.c_str()), "PLE");
    _leg->SetY1(_leg->GetY1()-0.05);
  }

  return p;
} // drawEta

// Step 2. Project profile to histogram, normalize by |eta|<1.3 and draw it
TH1D *drawNormEta(TProfile *p, string draw, int marker, int color) {

  assert(p);
  string id = Form("%s_%s_%d_%d",p->GetName(),draw.c_str(),marker,color);
  TH1D *h = p->ProjectionX(Form("h%s",id.c_str()));
  TF1 *f1 = new TF1(Form("f1%s_",id.c_str()),"[0]",0,1.305);
  h->Fit(f1,"QRN");
  f1->SetRange(0,5.2);
  h->Divide(f1);
  tdrDraw(h,draw,marker,color,kSolid,-1,kNone);
  
  return h;
}

// Step 3,4. Take a ratio and draw it
TH1D *drawRatio(TH1D *h, TH1D *hm, string draw, int marker, int color) {
  assert(h);
  assert(hm);
  string id = Form("%s_%s_%s_%d_%d", h->GetName(), hm->GetName(),
		   draw.c_str(), marker, color);
  TH1D *hr = (TH1D*)h->Clone(Form("hr_%s",id.c_str()));
  hr->Divide(hm);

  tdrDraw(hr,draw,marker,color,kSolid,-1,kNone);
  
  return hr;
}

void L2Res() {

  // Set graphical styles
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // Set output directory;
  TFile *fout = new TFile("rootfiles/L2Res.root","RECREATE");
  
  string vrun[] = {"2023Cv123","2023Cv4","2023D"};
  const int nrun = sizeof(vrun)/sizeof(vrun[0]);
  //string vmc[] = {"Summer22","Summer22","Summer22"};
  string vmc[] = {"Summer23","Summer23","Summer23BPIX"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);
  assert(nmc==nrun);
  for (int irun = 0; irun != nrun; ++irun) {
    string run = vrun[irun];//"2023D";
    const char *cr = run.c_str();
    string mc = vmc[irun];//"Summer22";
    const char *cm = mc.c_str();
  // (No indent here for the resf of the loop, maybe function call later)

  /*
  // Mikko's code to create L2Res.root file
  // Load Z+jet
  //TFile *fz = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_v66.root",cr),"READ"); // Summer22
    TFile *fz = new TFile(Form("rootfiles/Winter23_noL2L3Res/jme_bplusZ_%s_Zmm_sync_v67.root",cr),"READ"); // Winter23
  assert(fz && !fz->IsZombie());
  fz->cd("data/l2res");
  TDirectory *dz = gDirectory;
  fz->cd("mc/l2res");
  TDirectory *dzm = gDirectory;

  // Load G+jet
  //TFile *fg = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s_v32.root",cr),"READ"); // Summer22
  TFile *fg = new TFile(Form("rootfiles/Winter23_noL2L3Res/GamHistosFill_data_%s_w1.root",cr),"READ"); // Winter23
  assert(fg && !fg->IsZombie());
  fg->cd("Gamjet2");
  TDirectory *dg = gDirectory;
  //
  //TFile *fgm = new TFile("../gamjet/rootfiles/GamHistosFill_mc_2022P8_v32.root","READ");
  TFile *fgm = new TFile(run=="2023D" ? "rootfiles/Winter23_noL2L3Res/GamHistosFill_mc_2023P8-BPix_w1.root" : "rootfiles/Winter23_noL2L3Res/GamHistosFill_mc_2023P8_w1.root","READ");
  assert(fgm && !fgm->IsZombie());
  fgm->cd("Gamjet2");
  TDirectory *dgm = gDirectory;

  // Load dijet
  //TFile *fd = new TFile(Form("../dijet/rootfiles/jmenano_data_cmb_%s_JME_v35a.root",cr),"READ"); // Summer22
  TFile *fd = new TFile(Form("rootfiles/Winter23_noL2L3Res/jmenano_data_cmb_%s_JME_v36_Summer23DT_NoL2L3Res.root",cr),"READ"); // Winter23
  assert(fd && !fd->IsZombie());
  fd->cd("Dijet2");
  TDirectory *dd = gDirectory;
  //
  //TFile *fdm = new TFile("../dijet/rootfiles/jmenano_mc_cmb_Summer22MG_v35a.root","READ"); // Summer22
  TFile *fdm = new TFile(run=="2023D" ? "rootfiles/Winter23_noL2L3Res/jmenano_mc_cmb_Summer23MGBPix_v36_Summer23DT_NoL2L3Res.root" : "rootfiles/Winter23_noL2L3Res/jmenano_mc_cmb_Summer23MG_v36_Summer23DT_NoL2L3Res.root","READ"); // Winter23
  assert(fdm && !fdm->IsZombie());
  fdm->cd("Dijet2");
  TDirectory *ddm = gDirectory;
  
  // Get the TProfile2D inputs. Z+jet, gam+jet, 3x dijet
  
  // Z+jet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2z = (TProfile2D*)dz->Get("p2m0tc"); assert(p2z);
  TProfile2D *p2zm = (TProfile2D*)dzm->Get("p2m0tc"); assert(p2zm);

  // G+jet: x:eta, y:pT,avp (no other variants yet)
  TProfile2D *p2g = (TProfile2D*)dg->Get("p2m0"); assert(p2g);
  TProfile2D *p2gm = (TProfile2D*)dgm->Get("p2m0"); assert(p2gm);

  // Dijet: x:eta, y:pT,avp (others p2m0tc for pT,tag, p2m0pf for pT,probe )
  TProfile2D *p2j = (TProfile2D*)dd->Get("p2m0tc"); assert(p2j);
  TProfile2D *p2jm = (TProfile2D*)ddm->Get("p2m0tc"); assert(p2jm);
  TProfile2D *p2p = (TProfile2D*)dd->Get("p2m0pf"); assert(p2p);
  TProfile2D *p2pm = (TProfile2D*)ddm->Get("p2m0pf"); assert(p2pm);
  TProfile2D *p2d = (TProfile2D*)dd->Get("p2m0"); assert(p2d);
  TProfile2D *p2dm = (TProfile2D*)ddm->Get("p2m0"); assert(p2dm);
  */

  // Use common file instead
  TFile *f = new TFile("rootfiles/L2Res_2023_Winter23.root","READ");
  assert(f && !f->IsZombie());

  TProfile2D *p2z, *p2zm, *p2g, *p2gm;
  p2z = (TProfile2D*)f->Get(Form("p2m0tc_zjet_data_%s",cr)); assert(p2z);
  p2zm = (TProfile2D*)f->Get(Form("p2m0tc_zjet_mc_%s",cm)); assert(p2zm);
  p2g = (TProfile2D*)f->Get(Form("p2m0tc_gjet_data_%s",cr)); assert(p2g);
  p2gm = (TProfile2D*)f->Get(Form("p2m0tc_gjet_mc_%s",cm)); assert(p2gm);

  TProfile2D *p2j, *p2jm, *p2p, *p2pm, *p2d, *p2dm;
  p2j = (TProfile2D*)f->Get(Form("p2m0tc_dijet_data_%s",cr)); assert(p2j);
  p2jm = (TProfile2D*)f->Get(Form("p2m0tc_dijet_mc_%s",cm)); assert(p2jm);
  p2p = (TProfile2D*)f->Get(Form("p2m0pf_dijet_data_%s",cr)); assert(p2p);
  p2pm = (TProfile2D*)f->Get(Form("p2m0pf_dijet_mc_%s",cm)); assert(p2pm);
  p2d = (TProfile2D*)f->Get(Form("p2m0ab_dijet_data_%s",cr)); assert(p2d);
  p2dm = (TProfile2D*)f->Get(Form("p2m0ab_dijet_mc_%s",cm)); assert(p2dm);  

    
  // Loop over gamjet pT bins for plotting
  for (int ipt = 1; ipt != p2g->GetNbinsY()+1; ++ipt) {
    double ptmin = p2g->GetYaxis()->GetBinCenter(ipt);
    double ptmax = ptmin;
    double pt1 = p2g->GetYaxis()->GetBinLowEdge(ipt);
    double pt2 = p2g->GetYaxis()->GetBinLowEdge(ipt+1);
  // (No indent here for the resf of the loop, maybe function call later)

    
  // Step 1. Slice pT, draw response vs eta. No other manipulation yet
  TH1D *h1 = tdrHist("h1","JES",0.65,1.45,"|#eta|",0,5.2);
  lumi_136TeV = Form("%s - %s",cr,cm);
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h1,8,33,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(0,1,5.2,1);

  TLegend *leg1 = tdrLeg(0.20,0.90,0.45,0.90);
  _leg = leg1;
  
  //double ptmin(60), ptmax(80);
  //double ptmin(80), ptmax(120);
  //double ptmin(100), ptmax(100);
  //double ptmin(120), ptmax(120); // photon max

  TProfile *pzm, *pgm, *pdm, *pjm, *ppm;
  pzm = drawEta(p2zm,ptmin,ptmax,"HISTE",kNone,kRed);
  pgm = drawEta(p2gm,ptmin,ptmax,"HISTE",kNone,kBlue);
  pjm = drawEta(p2jm,ptmin,ptmax,"HISTE",kNone,kGreen+2);
  ppm = drawEta(p2pm,ptmin,ptmax,"HISTE",kNone,kOrange+2);
  pdm = drawEta(p2dm,ptmin,ptmax,"HISTE",kNone,kBlack);

  TProfile *pz, *pg, *pd, *pj, *pp;
  pz = drawEta(p2z,ptmin,ptmax,"Pz",kFullSquare,kRed,"Z");
  pg = drawEta(p2g,ptmin,ptmax,"Pz",kFullCircle,kBlue,"#gamma");
  pj = drawEta(p2j,ptmin,ptmax,"Pz",kFullDiamond,kGreen+2,"Tag");
  pp = drawEta(p2p,ptmin,ptmax,"Pz",kFullDiamond,kOrange+2,"Probe");
  pd = drawEta(p2d,ptmin,ptmax,"Pz",kOpenDiamond,kBlack,"Dijet");

  c1->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c1"));
    

  // Step 2. Project profile to histogram, normalize by |eta|<1.3
  TH1D *h2 = tdrHist("h2","Rel. JES",0.75,1.35,"|#eta|",0,5.2);
  TCanvas *c2 = tdrCanvas("c2",h2,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
  
  TH1D *hzm, *hgm, *hdm, *hjm, *hpm;
  hzm = drawNormEta(pzm,"HISTE",kNone,kRed);
  hgm = drawNormEta(pgm,"HISTE",kNone,kBlue);
  hjm = drawNormEta(pjm,"HISTE",kNone,kGreen+2);
  hpm = drawNormEta(ppm,"HISTE",kNone,kOrange+2);
  hdm = drawNormEta(pdm,"HISTE",kNone,kBlack);

  TH1D *hz, *hg, *hd, *hj, *hp;
  hz = drawNormEta(pz,"Pz",kFullSquare,kRed);
  hg = drawNormEta(pg,"Pz",kFullCircle,kBlue);
  hj = drawNormEta(pj,"Pz",kFullDiamond,kGreen+2);
  hp = drawNormEta(pp,"Pz",kFullDiamond,kOrange+2);
  hd = drawNormEta(pd,"Pz",kOpenDiamond,kBlack);

  c2->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c2"));
    
  
  // Step 3. Draw data/MC ratio before normalization
  TH1D *h3 = tdrHist("h3","JES Data/MC",0.80,1.15,"|#eta|",0,5.2);
  TCanvas *c3 = tdrCanvas("c3",h3,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
  
  TH1D *hzr, *hgr, *hdr, *hjr, *hpr;
  hzr = drawRatio(pz->ProjectionX(),pzm,"Pz",kFullSquare,kRed);
  hgr = drawRatio(pg->ProjectionX(),pgm,"Pz",kFullCircle,kBlue);
  hjr = drawRatio(pj->ProjectionX(),pjm,"Pz",kFullDiamond,kGreen+2);
  hpr = drawRatio(pp->ProjectionX(),ppm,"Pz",kFullDiamond,kOrange+2);
  hdr = drawRatio(pd->ProjectionX(),pdm,"Pz",kOpenDiamond,kBlack);

  c3->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c3"));
  
  
  // Step 4. Draw data/MC ratio of normalized JES
  TH1D *h4 = tdrHist("h4","Rel. JES Data/MC",0.80,1.15,"|#eta|",0,5.2);
  TCanvas *c4 = tdrCanvas("c4",h4,8,33,kSquare);

  l->DrawLine(0,1,5.2,1);

  leg1->Draw("SAME");
	     
  TH1D *hzrn, *hgrn, *hdrn, *hjrn, *hprn;
  hzrn = drawRatio(hz,hzm,"Pz",kFullSquare,kRed);
  hgrn = drawRatio(hg,hgm,"Pz",kFullCircle,kBlue);
  hjrn = drawRatio(hj,hjm,"Pz",kFullDiamond,kGreen+2);
  hprn = drawRatio(hp,hpm,"Pz",kFullDiamond,kOrange+2);
  hdrn = drawRatio(hd,hdm,"Pz",kOpenDiamond,kBlack);

  c4->SaveAs(Form("pdf/L2Res/vsEta/L2Res_vsEta_%04d_%04d_%s_%s.pdf",
		  int(pt1),int(pt2),cr,"c4"));

  // Rename to avoid loop leakage and errors
  h1->SetName(Form("h1_%s_%d",cr,ipt));
  h2->SetName(Form("h2_%s_%d",cr,ipt));
  h3->SetName(Form("h3_%s_%d",cr,ipt));
  h4->SetName(Form("h4_%s_%d",cr,ipt));
  
  c1->SetName(Form("c1_%s_%d",cr,ipt));
  c2->SetName(Form("c2_%s_%d",cr,ipt));
  c3->SetName(Form("c3_%s_%d",cr,ipt));
  c4->SetName(Form("c4_%s_%d",cr,ipt));

  } // for ipt

  // Save results for easier book-keeping
  if (fout) {
    fout->cd();
    p2z->Write(Form("p2m0tc_zjet_data_%s",cr),TObject::kOverwrite);
    p2zm->Write(Form("p2m0tc_zjet_mc_%s",cm),TObject::kOverwrite);
    p2g->Write(Form("p2m0tc_gjet_data_%s",cr),TObject::kOverwrite);
    p2gm->Write(Form("p2m0tc_gjet_mc_%s",cm),TObject::kOverwrite);
    p2j->Write(Form("p2m0tc_dijet_data_%s",cr),TObject::kOverwrite);
    p2jm->Write(Form("p2m0tc_dijet_mc_%s",cm),TObject::kOverwrite);
    p2p->Write(Form("p2m0pf_dijet_data_%s",cr),TObject::kOverwrite);
    p2pm->Write(Form("p2m0pf_dijet_mc_%s",cm),TObject::kOverwrite);
    p2d->Write(Form("p2m0ab_dijet_data_%s",cr),TObject::kOverwrite);
    p2dm->Write(Form("p2m0ab_dijet_mc_%s",cm),TObject::kOverwrite);
    curdir->cd();
  }

  } // for irun

  if (fout) {
    fout->Write();
    fout->Close();
  }
} // L2Res
