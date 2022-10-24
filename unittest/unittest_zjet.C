// Purpose: Unit test for Z+jet inputs
//          Compare main results to a stable reference:
//          [MPF, DB] x [MC, data, ratio]
// Author:  Mikko Voutilainen at cern dot ch
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include "../tdrstyle_mod22.C"
#include "../tools.C"

void unittest_zjet() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // ID string for test and reference
  string id, idref;
  const char *cdir = "data";

  // Set file to be tested against reference
  const char *cd = "../JERCProtoLab/Winter22Run3/L3Residual_Z";
  TFile *fnew = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9c.root",cd),"READ"); id = "RunE";//"v9ce"; // RunE (no corrections)
  //TFile *fnew = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9b.root",cd),"READ"); id = "v9be"; // RunE + CD_V2_L2L3Res input (only L2L3Res, no L2L3! bigger DY)
  //TFile *fnew = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9.root",cd),"READ"); id = "v9e"; // RunE + CD_V2_L2L3Res input (only L2L3Res, no L2L3! smaller DY)
  //TFile *fnew = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p0_Run2022CD_v7.root",cd),"READ"); id = "v7cd"; // RunC_V2_L2L3Res closure (problematic?)
  assert(fnew  && !fnew->IsZombie());
  const char *cid = id.c_str();
  
  // Set reference file
  //TFile *fref = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9b.root",cd),"READ"); idref = "v9be"; // RunE + CD_V2_L2L3Res (only L2L3Res, no L2L3! bigger DY)
  //TFile *fref = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p1_Run2022E_v9.root",cd),"READ"); idref = "v9e"; // RunE + CD_V2_L2L3Res input (only L2L3Res, no L2L3! smaller DY)
  //TFile *fref = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p0_Run2022CD_v7.root",cd),"READ"); idref = "v7cd"; // RunC_V2_L2L3Res closure (problematic?)
  TFile *fref = new TFile(Form("%s/jme_bplusZ_merged_2022muon_EOS_v2p0_Run2022CD_v6.root",cd),"READ"); idref = "RunCD";//"v6cd"; // RunC_V2_L2L3Res input (stable, no corrections applied)
  assert(fref  && !fref->IsZombie());
  const char *cidref = idref.c_str();
  
  curdir->cd();
  
  // Load graphs to be tested
  const char *cy = "eta_00_13";
  TGraphErrors *gm(0), *gd(0), *gmref(0), *gdref(0);
  gm = (TGraphErrors*)fnew->Get(Form("%s/%s/rmpf_zmmjet_a100",cdir,cy));
  gd = (TGraphErrors*)fnew->Get(Form("%s/%s/rmpfjet1_zmmjet_a100",cdir,cy));
  gmref = (TGraphErrors*)fref->Get(Form("%s/%s/rmpf_zmmjet_a100",cdir,cy));
  gdref = (TGraphErrors*)fref->Get(Form("%s/%s/rmpfjet1_zmmjet_a100",cdir,cy));

  assert(gm);
  assert(gd);
  assert(gmref);
  assert(gdref);

  TGraphErrors *gmr(0), *gdr(0);
  gmr = tools::ratioGraphs(gm,gmref);
  gdr = tools::ratioGraphs(gd,gdref);

  TLine *l = new TLine(); l->SetLineStyle(kDashed);
  TH1D *hup = tdrHist("hup","MPF or DB",0.75,1.35);
  TH1D *hdw = tdrHist("hdw",Form("%s / %s",cid,cidref),0.9,1.1);
  lumi_136TeV = "Run3";
  TCanvas *c1 = tdrDiCanvas("c1",hup,hdw,8,11);

  c1->cd(1);
  gPad->SetLogx();
  l->DrawLine(15,1,3500,1);

  tdrDraw(gmref,"PLz",kOpenDiamond,kGreen+2);
  tdrDraw(gdref,"PLz",kOpenCircle,kGreen+2);
  tdrDraw(gm,"PLz",kFullDiamond,kRed+1);
  tdrDraw(gd,"PLz",kFullCircle,kRed+1);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.045); tex->SetNDC();
  tex->DrawLatex(0.50,0.84,cdir);

  TLegend *legu = tdrLeg(0.63,0.88-4*0.05,0.88,0.88);
  legu->AddEntry(gm,Form("MPF %s",cid),"PE");
  legu->AddEntry(gd,Form("DB %s",cid),"PE");
  legu->AddEntry(gmref,Form("MPF %s",cidref),"PE");
  legu->AddEntry(gdref,Form("DB %s",cidref),"PE");

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(15,1,3500,1);

  tdrDraw(gmr,"PLz",kFullDiamond,kRed+1);
  tdrDraw(gdr,"PLz",kFullCircle,kRed+1);

  TLegend *legd = tdrLeg(0.88-2*0.12,0.88-2*0.05,0.88,0.88);
  legd->SetNColumns(2);
  legd->SetTextSize(2*0.045);
  legd->AddEntry(gmr,"MPF","PE");
  legd->AddEntry(gdr,"DB","PE");

  c1->SaveAs(Form("pdf/unittest/unittest_zjet_%s_%s_vs_%s.pdf",
		  cdir,cid,cidref));
  
} // unittest_zjet
