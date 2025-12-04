// Purpose: estimate JER SF vs eta from Prompt and ReReco
//          using Fikri's photon+jet samples
#include "TFile.h"
#include "TProfile2D.h"

#include "../tdrstyle_mod22.C"

TH1D *getRMS(TProfile *p, string name) {
  TH1D *h = p->ProjectionX(name.c_str());
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    string s = p->GetErrorOption();
    p->SetErrorOption("");
    double mean = p->GetBinContent(i);
    double emean = p->GetBinError(i);
    p->SetErrorOption("S");
    double rms = p->GetBinError(i);
    p->SetErrorOption(s.c_str());
    h->SetBinContent(i, mean>0 ? rms/mean : 0);
    h->SetBinError(i, rms>0 ? emean/rms : 0);
  } // for i
  return h;
}

TH1D *subRMS(TH1D *hb, TH1D *hx, string name) {
  TH1D *h = (TH1D*)hb->Clone(name.c_str());
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double rmsb = hb->GetBinContent(i);
    double rmsx = hx->GetBinContent(i);
    double eb = hb->GetBinError(i);
    double ex = hx->GetBinError(i);
    double rms = sqrt(max(rmsb*rmsb - rmsx*rmsx,0.));
    //double rms = sqrt(fabs(rmsb*rmsb - rmsx*rmsx));
    double err = sqrt(pow(rms>0 ? eb/rms : 0., 2) +
		      pow(rms>0 ? ex/rms : 0., 2)) * rms;
    h->SetBinContent(i, rms);
    h->SetBinError(i, err);
  } // for i
  return h;
}

void scaleRMSerr(TH1D *h, double c) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    h->SetBinError(i, c * h->GetBinError(i));
  }
}

void FikriJERSF() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! mkdir pdf/FikriJERSF");
  gROOT->ProcessLine(".! touch pdf/FikriJERSF");

  TFile *fs = new TFile("rootfiles/Fikri_JetDepths_2025E_ReReco_v6p0/Histo_RAWToPFNano_Data25E_EGamma_FTV12Nov25.root","READ");
  //TFile *fs = new TFile("rootfiles/Fikri_JetDepths_2025E_ReReco_v6p0/Histo_RAWToPFNano_Data25E_EGamma_FTV31Oct25.root","READ"); // check old one
  assert(fs && !fs->IsZombie());
  //TFile *fr = new TFile("rootfiles/Fikri_JetDepths_2025E_ReRecoV2/Histo_RAWToPFNano_Data25E_EGamma_FTV31Oct25.root","READ"); //v1
  TFile *fr = new TFile("rootfiles/Fikri_JetDepths_2025E_ReReco_v6p0/Histo_RAWToPFNano_Data25E_EGamma_FTV31Oct25v2.root","READ");
  assert(fr && !fr->IsZombie());
  //TFile *fb = new TFile("rootfiles/Fikri_JetDepths_2025E_ReRecoV2/Histo_RAWToPFNano_Data25E_EGamma_Baseline.root","READ");
  TFile *fb = new TFile("rootfiles/Fikri_JetDepths_2025E_ReReco_v6p0/Histo_RAWToPFNano_Data25E_EGamma_Baseline.root","READ");
  assert(fb && !fb->IsZombie());
  //TFile *fm = new TFile("rootfiles/Fikri_JetDepths_2025E_ReRecoV2/Histo_MC25Winter_GJ_LO_MG.root","READ");
  TFile *fm = new TFile("rootfiles/Fikri_JetDepths_2025E_ReReco_v6p0/Histo_MC25Winter_GJ_LO_MG.root","READ");
  assert(fm && !fm->IsZombie());

  curdir->cd();
  
  TProfile2D *p2s(0), *p2xs(0);
  p2s = (TProfile2D*)fs->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpf"); assert(p2s);
  p2xs = (TProfile2D*)fs->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpfx"); assert(p2xs);
  TProfile2D *p2r(0), *p2xr(0), *p2b(0), *p2xb(0), *p2m(0), *p2xm;
  p2r = (TProfile2D*)fr->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpf"); assert(p2r);
  p2xr = (TProfile2D*)fr->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpfx"); assert(p2xr);
  p2b = (TProfile2D*)fb->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpf"); assert(p2b);
  p2xb = (TProfile2D*)fb->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpfx"); assert(p2xb);
  p2m = (TProfile2D*)fm->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpf"); assert(p2m);
  p2xm = (TProfile2D*)fm->Get("hp2_ProbeJet0_abseta_TagPhotonCand_pt_vs_resp_mpfx"); assert(p2xm);

  TH1D *h = tdrHist("h","JES proxy",0.80+1e-4,1.15-1e-4,"Jet |#eta|",0,5.2);
  //TH1D *h = tdrHist("h","JES proxy",0.75+1e-4,1.20-1e-4,"Jet |#eta|",0,5.2);
  lumi_136TeV = "2025E EGamma";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,1.02,5.2,1.02);
  l->DrawLine(0,0.98,5.2,0.98);
  
  int j1 = p2r->GetYaxis()->FindBin(130.);
  int j2 = p2r->GetYaxis()->FindBin(200.);
  TProfile *pb = p2b->ProfileX("pb",j1,j2,"");
  TProfile *pr = p2r->ProfileX("pr",j1,j2,"");
  TProfile *ps = p2s->ProfileX("ps",j1,j2,"");
  TProfile *pm = p2m->ProfileX("pm",j1,j2,"");
  TProfile *pxb = p2xb->ProfileX("pxb",j1,j2,"");
  TProfile *pxr = p2xr->ProfileX("pxr",j1,j2,"");
  TProfile *pxs = p2xs->ProfileX("pxs",j1,j2,"");
  TProfile *pxm = p2xm->ProfileX("pxm",j1,j2,"");

  tdrDraw(pm,"HIST",kNone,kRed+2,kSolid,-1,kNone,0);
  tdrDraw(pb,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
  tdrDraw(pr,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
  tdrDraw(ps,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.18,0.20,"Photon+jet 130<p_{T,#gamma}<200 GeV");
  tex->DrawLatex(0.18,0.15,"Photon110EB");
  
  TLegend *leg = tdrLeg(0.35,0.90-0.045*4,0.60,0.90);
  leg->AddEntry(pr,"FTV31Oct25v2 re-reco","PLE");
  leg->AddEntry(ps,"FTV25Nov25 re-reco","PLE");
  leg->AddEntry(pb,"Prompt data","PLE");
  leg->AddEntry(pm,"Winter25 MC (w/o QCD bkg)","LF");
  
  c1->SaveAs("pdf/FikriJERSF/FikriJERSF_Mean.pdf");

  
  TH1D *h2 = tdrHist("h2","RMS",0,0.7,"Jet |#eta|",0,5.2);
  TCanvas *c2 = tdrCanvas("c2",h2,8,11,kSquare);

  TH1D *hm = getRMS(pm,"hm");
  TH1D *hb = getRMS(pb,"hb");
  TH1D *hr = getRMS(pr,"hr");
  TH1D *hs = getRMS(ps,"hs");
  TH1D *hxm = getRMS(pxm,"hxm");
  TH1D *hxb = getRMS(pxb,"hxb");
  TH1D *hxr = getRMS(pxr,"hxr");
  TH1D *hxs = getRMS(pxs,"hxs");

  tdrDraw(hxm,"HE",kNone,kRed+2,kSolid,-1,kNone,0,1,2);
  tdrDraw(hxb,"HE",kNone,kRed+1,kSolid,-1,kNone,0,1,2);
  tdrDraw(hxr,"HE",kNone,kGray+2,kSolid,-1,kNone);
  tdrDraw(hxs,"HE",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(hm,"Pz",kOpenDiamond,kRed+1,kSolid,-1,kNone,0);
  tdrDraw(hb,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
  tdrDraw(hr,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
  tdrDraw(hs,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);

  TH1D *hm2 = subRMS(hm,hxm,"hm2");
  TH1D *hb2 = subRMS(hb,hxb,"hb2");
  TH1D *hr2 = subRMS(hr,hxr,"hr2");
  TH1D *hs2 = subRMS(hs,hxs,"hs2");

  TH1D *h3 = tdrHist("h3","JER proxy",0.1+1e-4,0.5-1e-4,"Jet |#eta|",0,5.2);
  TH1D *h3d = tdrHist("h3d","Ratio",0.8+1e-4,1.1-1e-4,"Jet |#eta|",0,5.2);
  TCanvas *c3 = tdrDiCanvas("c3",h3,h3d,8,11);

  c3->cd(1);
  tdrDraw(hm2,"HIST",kNone,kRed+1,kSolid,-1,kNone,0);
  tdrDraw(hb2,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
  tdrDraw(hr2,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone);
  tdrDraw(hs2,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);

  tex->DrawLatex(0.35,0.85,"Photon+jet 130<p_{T,#gamma}<200 GeV");
  tex->DrawLatex(0.35,0.80,"Photon110EB");
  
  TLegend *leg3 = tdrLeg(0.22,0.72-0.045*4,0.47,0.72);
  leg3->AddEntry(hr2,"FTV31Oct25v2 re-reco","PLE");
  leg3->AddEntry(hs2,"FTV12Nov25 re-reco","PLE");
  leg3->AddEntry(hb2,"Prompt data","PLE");
  leg3->AddEntry(hm2,"Winter25  MC (w/o QCD bkg)","LF");
  
  c3->cd(2);
  TH1D *hrb = (TH1D*)hr2->Clone("hrb");
  hrb->Divide(hb2);
  TH1D *hsb = (TH1D*)hs2->Clone("hsb");
  hsb->Divide(hb2);

  l->SetLineStyle(kDashed);
  l->DrawLine(0,1,5.2,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(0,0.9,5.2,0.9);
  
  scaleRMSerr(hrb,0.3);
  tdrDraw(hrb,"Pz",kOpenCircle,kGray+1,kSolid,-1,kNone);
  scaleRMSerr(hsb,0.3);
  tdrDraw(hsb,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);

  c3->SaveAs("pdf/FikriJERSF/FikriJERSF_RMS.pdf");



  /////////////////////////////////////////
  // Plots vs pT (at high eta)           //
  /////////////////////////////////////////

  double vy[][4] =
    {{0.000, 1.305, 130, 1200},
     {1.305, 2.043, 130, 1000},
     {2.043, 2.500, 130, 600},
     {2.650, 2.964, 130, 600},
     {2.964, 3.489, 130, 600},
     {3.489, 5.191, 130, 600}
    };
  const int ny = sizeof(vy)/sizeof(vy[0]);
  for (int iy = 0; iy != ny; ++iy) {
  
    double etamin(vy[iy][0]), etamax(vy[iy][1]);
    double ptmin(vy[iy][2]), ptmax(vy[iy][3]);
    
    TH1D *h4 = tdrHist(Form("h4_%d",iy),"JES proxy",0.80+1e-4,1.15-1e-4,
		       "Jet p_{T} (GeV)",ptmin,ptmax);
    lumi_136TeV = "2025E EGamma";
    extraText = "Private";
    TCanvas *c4 = tdrCanvas(Form("c4_%d",iy),h4,8,11,kSquare);
    gPad->SetLogx();
    
    l->SetLineColor(kGray+1);
    l->SetLineStyle(kDashed);
    l->DrawLine(ptmin,1,ptmax,1);
    l->SetLineStyle(kDotted);
    l->DrawLine(ptmin,1.02,ptmax,1.02);
    l->DrawLine(ptmin,0.98,ptmax,0.98);
    
    int k1 = p2r->GetXaxis()->FindBin(etamin+0.01);
    int k2 = p2r->GetXaxis()->FindBin(etamax-0.01);
    double ymin = p2r->GetXaxis()->GetBinLowEdge(k1);
    double ymax = p2r->GetXaxis()->GetBinLowEdge(k2+1);
    
    TProfile *pptb = p2b->ProfileY(Form("pptb_%d",iy),k1,k2,"");
    TProfile *pptr = p2r->ProfileY(Form("pptr_%d",iy),k1,k2,"");
    TProfile *ppts = p2s->ProfileY(Form("ppts_%d",iy),k1,k2,"");
    TProfile *pptm = p2m->ProfileY(Form("pptm_%d",iy),k1,k2,"");
    TProfile *pptxb = p2xb->ProfileY(Form("pptxb_%d",iy),k1,k2,"");
    TProfile *pptxr = p2xr->ProfileY(Form("pptxr_%d",iy),k1,k2,"");
    TProfile *pptxs = p2xs->ProfileY(Form("pptxs_%d",iy),k1,k2,"");
    TProfile *pptxm = p2xm->ProfileY(Form("pptxm_%d",iy),k1,k2,"");

    tdrDraw(pptm,"HIST",kNone,kRed+2,kSolid,-1,kNone,0);
    tdrDraw(pptb,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
    tdrDraw(pptr,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
    tdrDraw(ppts,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);
    
    tex->SetNDC(); tex->SetTextSize(0.045);
    tex->DrawLatex(0.18,0.20,Form("Photon+jet %1.3f<#eta_{jet}<%1.3f",ymin,ymax));
    tex->DrawLatex(0.18,0.15,"Photon110EB");
    
    TLegend *leg4 = tdrLeg(0.35,0.90-0.045*4,0.60,0.90);
    leg4->AddEntry(ppts,"FTV12Nov25 re-reco","PLE");
    leg4->AddEntry(pptr,"FTV31Oct25 re-reco","PLE");
    leg4->AddEntry(pptb,"Prompt data","PLE");
    leg4->AddEntry(pptm,"Winter25 MC (w/o QCD bkg)","LF");
  
    c4->SaveAs(Form("pdf/FikriJERSF/FikriJERSF_VsPt_eta%d_%d_Mean.pdf",
		    int(1000.*ymin),int(1000*ymax)));
    
  
    TH1D *h5 = tdrHist(Form("h5_%d",iy),"RMS",0,0.7,"Jet p_{T} (GeV)",ptmin,ptmax);
    TCanvas *c5 = tdrCanvas(Form("c5_%d",iy),h5,8,11,kSquare);
    gPad->SetLogx();
    
    TH1D *hptm = getRMS(pptm,Form("hptm_%d",iy));
    TH1D *hptb = getRMS(pptb,Form("hptb_%d",iy));
    TH1D *hptr = getRMS(pptr,Form("hptr_%d",iy));
    TH1D *hpts = getRMS(ppts,Form("hpts_%d",iy));
    TH1D *hptxm = getRMS(pptxm,Form("hptxm_%d",iy));
    TH1D *hptxb = getRMS(pptxb,Form("hptxb_%d",iy));
    TH1D *hptxr = getRMS(pptxr,Form("hptxr_%d",iy));
    TH1D *hptxs = getRMS(pptxs,Form("hptxs_%d",iy));

    tdrDraw(hptxm,"HE",kNone,kRed+2,kSolid,-1,kNone,0,1,2);
    tdrDraw(hptxb,"HE",kNone,kRed+1,kSolid,-1,kNone,0,1,2);
    tdrDraw(hptxr,"HE",kNone,kGray+2,kSolid,-1,kNone);
    tdrDraw(hptxs,"HE",kNone,kBlack,kSolid,-1,kNone);
    tdrDraw(hptm,"Pz",kOpenDiamond,kRed+1,kSolid,-1,kNone,0);
    tdrDraw(hptb,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
    tdrDraw(hptr,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
    tdrDraw(hpts,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);
    
    TH1D *hptm2 = subRMS(hptm,hptxm,Form("hm2_%d",iy));
    TH1D *hptb2 = subRMS(hptb,hptxb,Form("hb2_%d",iy));
    TH1D *hptr2 = subRMS(hptr,hptxr,Form("hr2_%d",iy));
    TH1D *hpts2 = subRMS(hpts,hptxs,Form("hs2_%d",iy));

    TH1D *h6 = tdrHist(Form("h6_%d",iy),"JER proxy",0.,0.60-1e-4,
		       "Jet p_{T} (GeV)",ptmin,ptmax);
    TH1D *h6d = tdrHist(Form("h6d_%d",iy),"Ratio",0.7+1e-4,1.3-1e-4,
			"Jet p_{T} (GeV)",ptmin,ptmax);
    TCanvas *c6 = tdrDiCanvas(Form("c6_%d",iy),h6,h6d,8,11);
    
    c6->cd(1);
    gPad->SetLogx();  
    
    tdrDraw(hptm2,"HIST",kNone,kRed+1,kSolid,-1,kNone,0);
    tdrDraw(hptb2,"Pz",kOpenSquare,kRed,kSolid,-1,kNone,0,1,2);
    tdrDraw(hptr2,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
    tdrDraw(hpts2,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);
    
    tex->DrawLatex(0.35,0.85,Form("Photon+jet %1.3f<#eta_{jet}<%1.3f",ymin,ymax));
    tex->DrawLatex(0.35,0.80,"Photon110EB");
    
    TLegend *leg6 = tdrLeg(0.22,0.72-0.045*4,0.47,0.72);
    leg6->AddEntry(hpts2,"FTV12Nov25 re-reco","PLE");
    leg6->AddEntry(hptr2,"FTV31Oct25 re-reco","PLE");
    leg6->AddEntry(hptb2,"Prompt data","PLE");
    leg6->AddEntry(hptm2,"Winter25  MC (w/o QCD bkg)","LF");
    
    gPad->RedrawAxis();
    
    c6->cd(2);
    gPad->SetLogx();
    
    TH1D *hptrb = (TH1D*)hptr2->Clone(Form("hptrb_%d",iy));
    hptrb->Divide(hptb2);
    TH1D *hptsb = (TH1D*)hpts2->Clone(Form("hptsb_%d",iy));
    hptsb->Divide(hptb2);
    
    l->SetLineStyle(kDashed);
    l->DrawLine(ptmin,1,ptmax,1);
    l->SetLineStyle(kDotted);
    l->DrawLine(ptmin,0.9,ptmax,0.9);
    
    scaleRMSerr(hptrb,0.3);
    tdrDraw(hptrb,"Pz",kOpenCircle,kGray+2,kSolid,-1,kNone);
    scaleRMSerr(hptsb,0.3);
    tdrDraw(hptsb,"Pz",kFullCircle,kBlack,kSolid,-1,kNone);

    c6->SaveAs(Form("pdf/FikriJERSF/FikriJERSF_VsPt_eta%d_%d_RMS.pdf",
		    int(1000*ymin),int(1000*ymax)));
  } // for iy
} // FikriJERSF
