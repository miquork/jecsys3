// Purpose: Plots from Fikri's Z/gamma+jet HCAL depth data
//          1) Plot HB depth 1 data/MC vs pT for 24(CD)E and 24F(GHI)
//            1b) HO data/MC
//          2) Plot Data/MC for all depths for 24(CD)E (and 24F(GHI)?)
#include "TFile.h"
#include "TProfile2D.h"

#include <map>

#include "../tdrstyle_mod22.C"

const bool scaleFSR = true;
const double kfsr_data = 0.90;//0.91;//0.89;
const double kfsr_mc = 0.90;

void FikriDepths() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  gROOT->ProcessLine(".! touch pdf/FikriDepths");

  //string se("24E");
  //string se("25C");
  string se("25Cv2");
  //string sf("24F");
  //string sf("25C");
  //string sf("24E");
  string sf("24Iv2");
  string sm("S24");
  //string sm("W25");

  //string shfd("24Iv2");
  string shfd("25Cv2");
  string shfm("S24");
  
  map<string,const char*> mapsm;
  mapsm["S24"] = "MC24_GJ_Sherpa";
  mapsm["W25"] = "MC25Winter_GJ_LO_MG";

  const char *ce = se.c_str();
  const char *cf = sf.c_str();
  const char *cm = sm.c_str();
  const char *cm2 = mapsm[sm];

  const char *chfd = shfd.c_str();
  const char *chfm = shfm.c_str();
  const char *chfm2 = mapsm[shfm];
  
  //TFile *fe = new TFile("rootfiles/Fikri_JetDepths/Histo_Data24E_EGamma.root",
  //TFile *fe = new TFile(Form("rootfiles/Fikri_JetDepths_v2/Histo_Data%s_EGamma.root",ce),
  TFile *fe = new TFile(Form("rootfiles/Fikri_JetDepths_2025C/Histo_Data%s_EGamma.root",ce),
			"READ");
  assert(fe && !fe->IsZombie());

  //TFile *ff = new TFile("rootfiles/Fikri_JetDepths/Histo_Data24F_EGamma.root",
  TFile *ff = new TFile(Form("rootfiles/Fikri_JetDepths_v2/Histo_Data%s_EGamma.root",cf),
			"READ");
  assert(ff && !ff->IsZombie());

  //TFile *fm = new TFile("rootfiles/Fikri_JetDepths/Histo_MC24_GJ_Sherpa.root",
  TFile *fm = new TFile(Form("rootfiles/Fikri_JetDepths_v2/Histo_%s.root",cm2),
			"READ");
  assert(fm && !fm->IsZombie());

  //TFile *fm2 = new TFile("rootfiles/Fikri_JetDepths_v2/Histo_MC25Winter_GJ_LO_MG.root","READ");
  //assert(fm2 && !fm2->IsZombie());

  //TFile *fhfd = new TFile(Form("rootfiles/Fikri_JetDepths_v2/Histo_Data%s_EGamma.root",chfd),"READ");
  TFile *fhfd = new TFile(Form("rootfiles/Fikri_JetDepths_2025C/Histo_Data%s_EGamma.root",chfd),"READ");
  assert(fhfd && !fhfd->IsZombie());
  TFile *fhfm = new TFile(Form("rootfiles/Fikri_JetDepths_v2/Histo_%s.root",chfm2),"READ");
  assert(fhfm && !fhfm->IsZombie());
  
  curdir->cd();

  /////////////////////////////////////////////////
  // Step 1. Check HB depth 1 time/pT dependence //
  /////////////////////////////////////////////////

  TProfile2D *p2e(0), *p2f(0), *p2c(0), *p2m(0);
  //p2e = (TProfile2D*)fe->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth1"); assert(p2e);
  p2e = (TProfile2D*)fe->Get("hp2_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth1"); assert(p2e);
  p2f = (TProfile2D*)ff->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth1"); assert(p2f);
  p2f = (TProfile2D*)ff->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth1"); assert(p2f);
  p2m = (TProfile2D*)fm->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth1"); assert(p2m);

  // 1b. Same for HO
  TProfile2D *p2e_ho(0), *p2f_ho(0), *p2m_ho(0);
  //p2e_ho = (TProfile2D*)fe->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth8"); assert(p2e_ho);
  p2e_ho = (TProfile2D*)fe->Get("hp2_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth8"); assert(p2e_ho);
  p2f_ho = (TProfile2D*)ff->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth8"); assert(p2f_ho);
  p2m_ho = (TProfile2D*)fm->Get("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth8"); assert(p2m_ho);

  // 1b. Same for HF
  TProfile2D *p2e_hf(0), *p2f_hf(0), *p2m_hf(0);
  p2e_hf = (TProfile2D*)fe->Get("hp2_tnpHFLongShort_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2e_hf);
  p2f_hf = (TProfile2D*)ff->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2f_hf);
  p2m_hf = (TProfile2D*)fm->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2m_hf);
  //TProfile2D *p2e_hfs(0), *p2f_hfs(0), *p2m_hfs(0);

  
  TH1D *hu = tdrHist("hu","Jet energy fraction vs tag",0,0.3,
		     "Jet #eta",-3.15,3.15);
  TH1D *hd = tdrHist("hd","Data / MC",0.1,1.7,
		     "Jet #eta",-3.15,3.15);
  extraText = "Private";
  //lumi_136TeV = "2024E (Re-reco) vs 2024F (Prompt)";
  lumi_136TeV = Form("20%s",ce);
  TCanvas *c1 = tdrDiCanvas("c1",hu,hd,8,11);

  TH1D *hu_ho = tdrHist("hu_ho","Jet energy fraction vs tag (%)"
			,1e-5,3.5-1e-4,"Jet #eta",-1.5,1.5);//-3.15,3.15);
  TH1D *hd_ho = tdrHist("hd_ho","Data / MC",0.,20.,
			"Jet #eta",-1.5,1.5);//-3.15,3.15);
  TCanvas *c1_ho = tdrDiCanvas("c1_ho",hu_ho,hd_ho,8,11);

  TH1D *hu_hf = tdrHist("hu_hf","Jet energy fraction vs tag (%)"
			,25+1e-5,75.-1e-4,"Jet #eta",-5.19,5.19);
  TH1D *hd_hf = tdrHist("hd_hf","Data / MC",0.75,1.25,
			"Jet #eta",-5.19,5.19);
  TCanvas *c1_hf = tdrDiCanvas("c1_hf",hu_hf,hd_hf,8,11);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);

  c1->cd(1);
  tex->DrawLatex(0.35,0.845,"Depth 1");
  tex->DrawLatex(0.35,0.795,"#gamma + jet DB");
  tex->DrawLatex(0.35,0.745,"Photon50EB");

  TLegend *leg55 = tdrLeg(0.74,0.88-0.05*4,0.99,0.88);
  TLegend *leg110 = tdrLeg(0.67,0.88-0.05*4,0.92,0.88);
  TLegend *leg230 = tdrLeg(0.60,0.88-0.05*4,0.85,0.88);

  c1_ho->cd(1);
  tex->DrawLatex(0.35,0.845,"HO");
  tex->DrawLatex(0.35,0.795,"#gamma + jet DB");
  tex->DrawLatex(0.35,0.745,"Photon50EB");

  c1_hf->cd(1);
  tex->DrawLatex(0.35,0.845,"HF");
  tex->DrawLatex(0.35,0.795,"#gamma + jet DB");
  tex->DrawLatex(0.35,0.745,"Photon50EB");

  //TLegend *leg55_ho = tdrLeg(0.74,0.88-0.05*4,0.99,0.88);
  //TLegend *leg110_ho = tdrLeg(0.67,0.88-0.05*4,0.92,0.88);
  //TLegend *leg230_ho = tdrLeg(0.60,0.88-0.05*4,0.85,0.88);
  
  const double vpt[] = {230, 110, 55};
  //const double vpt[] = {400, 110, 55};
  const int npt = sizeof(vpt)/sizeof(vpt[0]);
  int colorm[] = {kGray+1,kGray+2,kBlack};
  double sizem[] = {1, 2, 3};
  int markere[] = {kOpenSquare,kOpenSquare,kFullSquare};
  int colore[] = {kGreen+1-9,kGreen+1,kGreen+2};
  double sizee[] = {0.3,0.4,0.5};
  int markerf[] = {kOpenCircle,kOpenCircle,kFullCircle};
  int colorf[] = {kRed-9,kRed,kRed};
  double sizef[] = {0.3,0.4,0.5};
  for (int ipt = 0; ipt != npt; ++ipt) {
    
    c1->cd(1);
    
    //double ptref = 55;
    double ptref = vpt[ipt];
    int k = p2e->GetXaxis()->FindBin(ptref);
    TProfile *p1e = p2e->ProfileY(Form("p1e_%d",k),k,k);
    TProfile *p1f = p2f->ProfileY(Form("p1f_%d",k),k,k);
    TProfile *p1m = p2m->ProfileY(Form("p1m_%d",k),k,k);
  
    tdrDraw(p1m,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone,0,0.5,sizem[ipt]);
    tdrDraw(p1e,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(p1f,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);

    if (ipt==npt-1) {
      leg55->SetHeader(" 55 GeV");
      leg55->AddEntry(p1f,Form("20%s",cf),"PLE");
      leg55->AddEntry(p1e,Form("20%s",ce),"PLE");
      leg55->AddEntry(p1m,Form("%s MC",cm),"FL");
    }
    else if (ipt==1) {
      leg110->SetHeader("110, ");
      leg110->AddEntry(p1f," ","PLE");
      leg110->AddEntry(p1e," ","PLE");
      leg110->AddEntry(p1m," ","FL");
    }
    else if (ipt==0) {
      leg230->SetHeader("230, ");
      leg230->AddEntry(p1f," ","PLE");
      leg230->AddEntry(p1e," ","PLE");
      leg230->AddEntry(p1m," ","FL");
    }

    c1_ho->cd(1);
    int k_ho = p2e_ho->GetXaxis()->FindBin(ptref);
    TProfile *p1e_ho = p2e_ho->ProfileY(Form("p1e_ho_%d",k_ho),k_ho,k_ho);
    TProfile *p1f_ho = p2f_ho->ProfileY(Form("p1f_ho_%d",k_ho),k_ho,k_ho);
    TProfile *p1m_ho = p2m_ho->ProfileY(Form("p1m_ho_%d",k_ho),k_ho,k_ho);

    TH1D *h1m_ho = p1m_ho->ProjectionX(Form("h1m_ho_%d",k_ho),"e");
    TH1D *h1e_ho = p1e_ho->ProjectionX(Form("h1e_ho_%d",k_ho),"e");
    TH1D *h1f_ho = p1f_ho->ProjectionX(Form("h1f_ho_%d",k_ho),"e");
    h1m_ho->Scale(100.);
    h1e_ho->Scale(100.);
    h1f_ho->Scale(100.);
    
    tdrDraw(h1m_ho,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone,0,0.5,sizem[ipt]);
    tdrDraw(h1e_ho,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(h1f_ho,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);


    c1_hf->cd(1);
    int k_hf = p2e_hf->GetXaxis()->FindBin(ptref);
    TProfile *p1e_hf = p2e_hf->ProfileY(Form("p1e_hf_%d",k_hf),k_hf,k_hf);
    TProfile *p1f_hf = p2f_hf->ProfileY(Form("p1f_hf_%d",k_hf),k_hf,k_hf);
    TProfile *p1m_hf = p2m_hf->ProfileY(Form("p1m_hf_%d",k_hf),k_hf,k_hf);

    TH1D *h1m_hf = p1m_hf->ProjectionX(Form("h1m_hf_%d",k_hf),"e");
    TH1D *h1e_hf = p1e_hf->ProjectionX(Form("h1e_hf_%d",k_hf),"e");
    TH1D *h1f_hf = p1f_hf->ProjectionX(Form("h1f_hf_%d",k_hf),"e");
    h1m_hf->Scale(100.);
    h1e_hf->Scale(100.);
    h1f_hf->Scale(100.);
    
    tdrDraw(h1m_hf,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone,0,0.5,sizem[ipt]);
    tdrDraw(h1e_hf,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(h1f_hf,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);
    
    
    c1->cd(2);

    TH1D *h1m = p1m->ProjectionX(Form("h1m_%d",k_ho));
    h1m->Divide(p1m);
    TH1D *h1e = p1e->ProjectionX(Form("h1e_%d",k_ho));
    h1e->Divide(p1m);
    TH1D *h1f = p1f->ProjectionX(Form("h1f_%d",k_ho));
    h1f->Divide(p1m);
    
    tdrDraw(h1m,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone, 0,0.5,sizem[ipt]);
    tdrDraw(h1e,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(h1f,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);

    c1_ho->cd(2);

    TH1D *h1e_rho = (TH1D*)h1e_ho->Clone(Form("h1e_rho_%d",k_ho));
    TH1D *h1f_rho = (TH1D*)h1f_ho->Clone(Form("h1f_rho_%d",k_ho));
    TH1D *h1m_rho = (TH1D*)h1m_ho->Clone(Form("h1m_rho_%d",k_ho));
    h1e_rho->Divide(h1m_ho);
    h1f_rho->Divide(h1m_ho);
    h1m_rho->Divide(h1m_ho);
    
    tdrDraw(h1m_rho,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone, 0,0.5,sizem[ipt]);
    tdrDraw(h1e_rho,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(h1f_rho,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);

    c1_hf->cd(2);

    TH1D *h1e_rhf = (TH1D*)h1e_hf->Clone(Form("h1e_rhf_%d",k_hf));
    TH1D *h1f_rhf = (TH1D*)h1f_hf->Clone(Form("h1f_rhf_%d",k_hf));
    TH1D *h1m_rhf = (TH1D*)h1m_hf->Clone(Form("h1m_rhf_%d",k_hf));
    h1e_rhf->Divide(h1m_hf);
    h1f_rhf->Divide(h1m_hf);
    h1m_rhf->Divide(h1m_hf);
    
    tdrDraw(h1m_rhf,"HIST][",kNone,colorm[ipt],kSolid,-1,kNone, 0,0.5,sizem[ipt]);
    tdrDraw(h1e_rhf,"Pz",markere[ipt],colore[ipt],kSolid,-1,kNone,0,sizee[ipt]);
    tdrDraw(h1f_rhf,"Pz",markerf[ipt],colorf[ipt],kSolid,-1,kNone,0,sizef[ipt]);
  } // for ipt

  c1_ho->cd(1);
  leg55->DrawClone("SAME");
  leg110->DrawClone("SAME");
  leg230->DrawClone("SAME");
  c1_ho->SaveAs(Form("pdf/FikriDepths/FikriDepths_HO_20%s_20%s.pdf",ce,cf));

  c1_hf->cd(1);
  leg55->DrawClone("SAME");
  leg110->DrawClone("SAME");
  leg230->DrawClone("SAME");
  c1_hf->SaveAs(Form("pdf/FikriDepths/FikriDepths_HF_20%s_20%s.pdf",ce,cf));
  
  c1->SaveAs(Form("pdf/FikriDepths/FikriDepths_HBdepth1_20%s_20%s.pdf",ce,cf));
  //c1->SaveAs("pdf/FikriDepths/FikriDepths_HBdepth1.gif");
  

  
  /////////////////////////////////////////////////////
  // Step 2. HCAL depth scale comparison to IsoTrack //
  /////////////////////////////////////////////////////

  lumi_136TeV = Form("20%s",ce);//"2024E (re-reco)";
  TH1D *hu_2 = tdrHist("hu_2","Jet energy fraction vs tag",0,0.62,
		       "Jet #eta",-3.15,3.15);
  //TH1D *hd_2 = tdrHist("hd_2","Data / MC",0,3.2,
  TH1D *hd_2 = tdrHist("hd_2","Data / MC",0.4,1.8,
		       "Jet #eta",-3.15,3.15);
  TCanvas *c2 = tdrDiCanvas("c2",hu_2,hd_2,8,11);

  c2->cd(1);
  //tex->DrawLatex(0.35,0.845,"Depth 1");

  double ptref = 55;
  //double ptref = 42;
  //double ptref = 55;
  //double ptref = 110;
  //double ptref = 230;
  //double ptref = 400;
  //double ptref = 600;
  //double ptref2(0);
  //double ptref2 = 420;
  double ptref2 = 250;
  int k1 = p2e->GetXaxis()->FindBin(ptref);
  int k2 = (ptref2!=0 ? p2e->GetXaxis()->FindBin(ptref2) : k1);
  double ptmin = p2e->GetXaxis()->GetBinLowEdge(k1);
  double ptmax = p2e->GetXaxis()->GetBinLowEdge(k2+1);
  tex->DrawLatex(0.35,0.855,Form("%1.0f < p_{T,#gamma} < %1.0f GeV",
				 ptmin,ptmax));
  tex->DrawLatex(0.35,0.790,"#gamma + jet DB");
  tex->DrawLatex(0.35,0.740,"Photon50EB");  

  tex->SetTextSize(0.035);
  tex->DrawLatex(0.43,0.52,Form("k_{FSR,data} = %1.3f",kfsr_data));
  tex->DrawLatex(0.43,0.48,Form("k_{FSR,MC} = %1.3f",kfsr_mc)); 
  
  TLegend *leg2 = tdrLeg(0.75,0.88-0.035*9,0.95,0.88);
  leg2->SetTextSize(0.035);

  c2->cd(2);

  TLatex *tex2 = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.035*1.5);
  tex->SetTextColor(kGray+2);
  tex->DrawLatex(0.385,0.66,"HE depths #5-7 in HB from loopers");
  
  const int ndepth = 9;
  const char* label2e[] =
    {"ECAL",
     "Depth 1", "Depth 2", "Depth 3", "Depth 4",
     "#5 (HE)", "#6 (HE)", "#7 (HE)",
     "HO (HB)"};
  int color2e[ndepth] =
    {kCyan+2, //kBlue+2,  // ECAL
     kBlue, kOrange+2, kGreen+2, kRed, // HB
     kYellow+1, kOrange+3, kGray+1, // HE
     kBlack};//kGray+2}; // HO
  int marker2e[ndepth] =
    {kFullDiamond, // ECAL
     /*kOpen*/kFullTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp, //HB
     kOpenStar, kOpenDiamond, kOpenCross, // HE
     kFullStar}; // HO

  map<int,TProfile*> mp1e;
  map<int,TProfile*> mp1m;
  for (int id = 0; id != ndepth; ++id) {

    TProfile2D *p2e(0), *p2m(0);
    //p2e = (TProfile2D*)fe->Get(Form("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth%d",id)); assert(p2e);
    p2e = (TProfile2D*)fe->Get(Form("hp2_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth%d",id)); assert(p2e);
    p2m = (TProfile2D*)fm->Get(Form("hp2_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBDepth%d",id)); assert(p2m);

    TProfile *p1e = p2e->ProfileY(Form("p1e_depth%d",id),k1,k2);
    TProfile *p1m = p2m->ProfileY(Form("p1m_depth%d",id),k1,k2);

    mp1e[id] = p1e;
    mp1m[id] = p1m;

    TH1D *h1m = (TH1D*)p1m->ProjectionX(Form("h1m_depth%d",id));
    if (scaleFSR) h1m->Scale(1./kfsr_mc);
    TH1D *h1e = (TH1D*)p1e->ProjectionX(Form("h1e_depth%d",id));
    if (scaleFSR) h1e->Scale(1./kfsr_data);
    
    c2->cd(1);
    tdrDraw(h1m,"HIST][",kNone,color2e[id],kSolid,-1,kNone,0);
    tdrDraw(h1e,"Pz",marker2e[id],color2e[id],kSolid,-1,kNone,0,0.5);

    leg2->AddEntry(h1e,label2e[id],"PLE");
    
    c2->cd(2);
    TH1D *h1mr = (TH1D*)p1m->ProjectionX(Form("h1mr_depth%d",id));
    if (scaleFSR) h1mr->Scale(1./kfsr_mc);
    h1mr->Divide(h1m);
    TH1D *h1er = (TH1D*)p1e->ProjectionX(Form("h1er_depth%d",id));
    if (scaleFSR) h1er->Scale(1./kfsr_data);
    h1er->Divide(h1m);

    tdrDraw(h1mr,"HIST][",kNone,color2e[id],kSolid,-1,kNone,0);
    tdrDraw(h1er,"Pz",marker2e[id],color2e[id],kSolid,-1,kNone,0,0.5);
  } // for depth

  //c2->SaveAs(Form("pdf/FikriDepths/FikriDepths_20%s.pdf",ce));
  c2->SaveAs(Form("pdf/FikriDepths/FikriDepths_20%s_%03.0f-%03.0f.pdf",ce,
		  ptmin,ptmax));

  
  // Make a stack plot also just to be fancy
  TH1D *h_2 = tdrHist("h_2","Jet energy fraction vs tag",0,1.25,
		      "Jet #eta",-3.15,3.15);
  TCanvas *c2s = tdrCanvas("c2s",h_2,8,11,kSquare);

  TLegend *leg2s = tdrLeg(0.35,0.90-0.035*3,0.95,0.90);
  leg2s->SetTextSize(0.035);
  leg2s->SetNColumns(3);

  const char* label2se[] =
    {"ECAL",
     "#1", "#2", "#3", "#4",
     "#5 (HE)", "#6 (HE)", "#7 (HE)",
     "HO (HB)"};

  for (int id = ndepth-1; id != -1; --id) {
    TH1D *hm = mp1m[id]->ProjectionX(Form("hstackm_depth%d",id));
    hm->Scale(1./kfsr_mc);
    for (int jd = 0; jd != id; ++jd) {
      assert(mp1m[id]);
      TH1D *h2m = mp1m[jd]->ProjectionX(Form("hstackm_depth%d_depth%d",id,jd));
      h2m->Scale(1./kfsr_mc);
      hm->Add(h2m);
    }
    tdrDraw(hm,"HIST][",kNone,color2e[id],kSolid,-1,1001,color2e[id]);
    //leg2s->AddEntry(hm,label2se[id],"F");
  } // for jd

  for (int id = ndepth-1; id != -1; --id) {
    TH1D *he = mp1e[id]->ProjectionX(Form("hstacke_depth%d",id));
    he->Scale(1./kfsr_data);
    for (int jd = 0; jd != id; ++jd) {
      assert(mp1e[id]);
      TH1D *h2e = mp1e[jd]->ProjectionX(Form("hstacke_depth%d_depth%d",id,jd));
      h2e->Scale(1./kfsr_data);
      he->Add(h2e);
    }
    tdrDraw(he,"Pz",marker2e[id],color2e[id]+1,kSolid,-1,1001,color2e[id],0.5);
    leg2s->AddEntry(he,label2se[id],"FPLE");
  } // for jd

  l->DrawLine(-3.15,1,3.15,1);
  
  tex->SetTextSize(0.035);
  tex->SetTextColor(kCyan+3);
  tex->DrawLatex(0.43,0.30,"#gamma + jet");
  tex->DrawLatex(0.43,0.26,"Photon50EB");
  tex->DrawLatex(0.43,0.22,Form("k_{FSR,data} = %1.3f",kfsr_data));
  tex->DrawLatex(0.43,0.18,Form("k_{FSR,MC} = %1.3f",kfsr_mc)); 
  
  gPad->RedrawAxis();
  
  c2s->SaveAs(Form("pdf/FikriDepths/FikriDepthsStack_20%s.pdf",ce));


  //////////////////////////////////////
  // Step 3. HF long/short comparison //
  //////////////////////////////////////
  TProfile2D *p2d_hfe(0), *p2d_hfh(0), *p2d_hfl(0), *p2d_hfs(0);
  TProfile2D *p2m_hfe(0), *p2m_hfh(0), *p2m_hfl(0), *p2m_hfs(0);

  //p2d_hfl = (TProfile2D*)fhfd->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2d_hfl);
  p2d_hfl = (TProfile2D*)fhfd->Get("hp2_tnpHFLongShort_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2d_hfl);
  p2m_hfl = (TProfile2D*)fhfm->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfLong"); assert(p2m_hfl);

  //p2d_hfs = (TProfile2D*)fhfd->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfShort"); assert(p2d_hfs);
  p2d_hfs = (TProfile2D*)fhfd->Get("hp2_tnpHFLongShort_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfShort"); assert(p2d_hfs);
  p2m_hfs = (TProfile2D*)fhfm->Get("hp2_tnpHFLongShort_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfShort"); assert(p2m_hfs);

  
  //p2d_hfe = (TProfile2D*)fhfd->Get("hp2_tnpHF_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfem"); assert(p2d_hfe);
  //p2d_hfe = (TProfile2D*)fhfd->Get("hp2_tnpHF_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfem"); assert(p2d_hfe);
  //p2m_hfe = (TProfile2D*)fhfm->Get("hp2_tnpHF_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfem"); assert(p2m_hfe);

  //p2d_hfh = (TProfile2D*)fhfd->Get("hp2_tnpHF_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfhad"); assert(p2d_hfh);
  //p2d_hfh = (TProfile2D*)fhfd->Get("hp2_tnpHF_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfhad"); assert(p2d_hfh);
  //p2m_hfh = (TProfile2D*)fhfm->Get("hp2_tnpHF_ProbePhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDBhfhad"); assert(p2m_hfh);
  

  TH1D *hu_hf3 = tdrHist("hu_hf3","Jet energy fraction vs tag",0,1.0,
			"Jet #eta",-5.19,5.19);
  TH1D *hd_hf3 = tdrHist("hd_hf3","Data / MC",0.5,1.5,
			"Jet #eta",-5.19,5.19);
  extraText = "Private";
  lumi_136TeV = Form("20%s vs %s",chfd,chfm);
  TCanvas *c3 = tdrDiCanvas("c3",hu_hf3,hd_hf3,8,11);

  double pt_hf(60.);
  int khf = p2d_hfl->GetXaxis()->FindBin(pt_hf);
  int ptmin_hf = p2d_hfl->GetXaxis()->GetBinLowEdge(khf);
  int ptmax_hf = p2d_hfl->GetXaxis()->GetBinLowEdge(khf+1);
  TH1D *hd_hfl = p2d_hfl->ProjectionY("hd_hfl",khf,khf);
  TH1D *hm_hfl = p2m_hfl->ProjectionY("hm_hfl",khf,khf);
  TH1D *hd_hfs = p2d_hfs->ProjectionY("hd_hfs",khf,khf);
  TH1D *hm_hfs = p2m_hfs->ProjectionY("hm_hfs",khf,khf);
  //
  /*
  TH1D *hd_hfe = p2d_hfe->ProjectionY("hd_hfe",khf,khf);
  TH1D *hm_hfe = p2m_hfe->ProjectionY("hm_hfe",khf,khf);
  TH1D *hd_hfh = p2d_hfh->ProjectionY("hd_hfh",khf,khf);
  TH1D *hm_hfh = p2m_hfh->ProjectionY("hm_hfh",khf,khf);
  */
  //
  TH1D *hd_hf2  = (TH1D*)hd_hfl->Clone("hd_hf2");   hd_hf2->Add(hd_hfs);
  //TH1D *hd_hf2b = (TH1D*)hd_hfe->Clone("hd_hf2b");  hd_hf2b->Add(hd_hfh);
  TH1D *hm_hf2  = (TH1D*)hm_hfl->Clone("hm_hf2");   hm_hf2->Add(hm_hfs);
  //TH1D *hm_hf2b = (TH1D*)hm_hfe->Clone("hm_hf2b");  hm_hf2b->Add(hm_hfh);
  
  c3->cd(1);

  tdrDraw(hm_hfl,"HIST][",kNone,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hd_hfl,"Pz",kFullSquare,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hm_hfs,"HIST][",kNone,kRed,kSolid,-1,kNone,0);
  tdrDraw(hd_hfs,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);

  /*
  tdrDraw(hm_hfe,"HIST][",kNone,kOrange+2,kSolid,-1,kNone,0);
  tdrDraw(hd_hfe,"Pz",kOpenSquare,kOrange+2,kSolid,-1,kNone,0);
  tdrDraw(hm_hfh,"HIST][",kNone,kMagenta+1,kSolid,-1,kNone,0);
  tdrDraw(hd_hfh,"Pz",kOpenCircle,kMagenta+1,kSolid,-1,kNone,0);
  */

  tdrDraw(hm_hf2,"HIST][",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hd_hf2,"Pz",kFullDiamond,kBlack,kSolid,-1,kNone,0);
  //tdrDraw(hm_hf2b,"HIST][",kNone,kGray+1,kDashed,-1,kNone,0);
  //tdrDraw(hd_hf2b,"Pz",kOpenDiamond,kGray+1,kDashed,-1,kNone,0);

  c3->cd(2);

  l->DrawLine(-5.19,1,5.19,1);
  
  TH1D *hr_hfl = (TH1D*)hd_hfl->Clone("hr_hfl"); hr_hfl->Divide(hm_hfl);
  TH1D *hr_hfs = (TH1D*)hd_hfs->Clone("hr_hfs"); hr_hfs->Divide(hm_hfs);
  TH1D *hr_hf2 = (TH1D*)hd_hf2->Clone("hr_hf2"); hr_hf2->Divide(hm_hf2);

  tdrDraw(hr_hfl,"Pz",kFullSquare,kBlue,kSolid,-1,kNone,0);
  tdrDraw(hr_hfs,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0);
  tdrDraw(hr_hf2,"Pz",kFullDiamond,kBlack,kSolid,-1,kNone,0);
  
  c3->SaveAs(Form("pdf/FikriDepths/FikriDepthsHF_20%s_vs_%s_pt%03d_%03d.pdf",chfd,chfm,ptmin_hf,ptmax_hf));
} // FikriDepths
