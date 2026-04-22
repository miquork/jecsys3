// Purpose: Draw jet eta spikes, and their high pT resolution, with
//          - 2026C low PU data vs 2026B high PU data
//          - Run 2 simulation with PFHC off vs PFHC on
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"

#include "../tdrstyle_mod22.C"

void drawEtaSpikesJER();

void drawEtaSpikes() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/drawEtaSpikes");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/drawEtaSpikes");

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  drawEtaSpikesJER(); exit(0);
  
  //TFile *flowpu = new TFile("rootfiles/Prompt2026/Jet_v160/jmenano_data_out_2026C_JME_v160.root","READ");
  //TFile *flowpu = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer24MC_Flat22_NoPFH_v161.root","READ");
  TFile *flowpu = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer22MC_hadCalibOff_1M_v161.root","READ");
  assert(flowpu && !flowpu->IsZombie());

  //TFile *fhighpu = new TFile("rootfiles/Prompt2026/Jet_v158/jmenano_data_out_2026B_JME_v158.root","READ");
  //TFile *fhighpu = new TFile("rootfiles/Prompt2026/Jet_v158/jmenano_mc_out_Summer24MC_Flat_JMENANO_v158.root","READ");
  TFile *fhighpu = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer22MC_base_1M_v161.root","READ");
  assert(fhighpu && !fhighpu->IsZombie());

  curdir->cd();

  bool isMC = true;
  double pt = 15;//30;
  TH1D *hlowpu(0), *hhighpu(0);
  if (pt==30) {
    cout << "Reading 30 GeV special" << endl;
    hlowpu = (TH1D*)flowpu->Get("HLT_ZeroBias/Incjet/hpteta30");//20");
    hhighpu = (TH1D*)fhighpu->Get("HLT_ZeroBias/Incjet/hpteta30");//20");
  }
  else {
    cout << "Reading 2D histograms" << endl;
    TH2D *h2lowpu = (TH2D*)flowpu->Get("HLT_ZeroBias/Incjet/h2pteta");
    if (!h2lowpu) h2lowpu =  (TH2D*)flowpu->Get("HLT_MC/Incjet/h2pteta");
    assert(h2lowpu);
    int i1 = h2lowpu->GetYaxis()->FindBin(15.);
    int i2 = h2lowpu->GetYaxis()->GetNbins();
    hlowpu = h2lowpu->ProjectionX("hlowpu",i1,i2);

    TH2D *h2highpu = (TH2D*)fhighpu->Get("HLT_ZeroBias/Incjet/h2pteta");
    if (!h2highpu) h2highpu = (TH2D*)fhighpu->Get("HLT_MC/Incjet/h2pteta");
    assert(h2highpu);
    hhighpu = h2highpu->ProjectionX("hhighpu",i1,i2);
  }
  assert(hlowpu);
  assert(hhighpu);

  TH1D *hmu_highpu = (TH1D*)fhighpu->Get("HLT_ZeroBias/Pileup/h_PUProfile");
  if (!hmu_highpu) hmu_highpu=(TH1D*)fhighpu->Get("HLT_MC/Pileup/h_PUProfile");
  assert(hmu_highpu);
  TH1D *hmu_lowpu = (TH1D*)flowpu->Get("HLT_ZeroBias/Pileup/h_PUProfile");
  if (!hmu_lowpu) hmu_lowpu=(TH1D*)flowpu->Get("HLT_MC/Pileup/h_PUProfile");
  assert(hmu_lowpu);
  
  double khb(1.0), khf(10.);
  double lumlowpu  = 0.000003509*1e15*1e-3; // fb-1 -> mb-1
  double lumhighpu = 0.000037337*1e15*1e-3; // fb-1
  //               = 0.000014077*1e15*1e-3;
  hlowpu = (TH1D*)hlowpu->Clone("hwlopu");
  if (isMC) cout << "hlowpu->Integral()="<<hlowpu->Integral()<<endl;
  //hlowpu->Scale(isMC ? 1./3.*1e-4 : 1./(lumlowpu), "width");
  hlowpu->Scale(isMC ? 2.5e-4 : 1./(lumlowpu), "width"); // Summer22MC
  //if (isMC) hlowpu->Scale(10./hlowpu->Integral());
  hhighpu = (TH1D*)hhighpu->Clone("hhighpu");
  if (isMC) cout << "hhighpu->Integral()="<<hhighpu->Integral()<<endl;
  //hhighpu->Scale(isMC ? 0.14*1./3.*1e-4 : 1./(lumhighpu), "width"); // /fb->mb
  hhighpu->Scale(isMC ? 2.5e-4 : 1./(lumhighpu), "width"); // Summer22 MC
  //if (isMC) hhighpu->Scale(10./hhighpu->Integral());
  
  TH1D *hlowpu_hb = (TH1D*)hlowpu->Clone("hwlopu_hb");
  hlowpu_hb->Scale(khb);//2.6); // -10% JES -> -50% xsec?
  TH1D *hlowpu_hf = (TH1D*)hlowpu->Clone("hwlopu_hf");
  hlowpu_hf->Scale(khf);

  
  TH1D *h = tdrHist("h","Jet cross section d#sigma/d#eta (mb)",0,1.5,//3.5,//3*1.5,
		    "Jet #eta",-5.2,5.2);
  lumi_136TeV = "2026B+C, HLT_ZeroBias";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  //tex->DrawLatex(0.40,0.85,"p_{T,jet} > 20 GeV");
  tex->DrawLatex(0.40,0.85,"p_{T,jet} > 15 GeV");

  tdrDraw(hhighpu,"HIST",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(hlowpu,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hlowpu,"Pz",kOpenCircle,kBlack,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hlowpu_hb,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7);
  //tdrDraw(hlowpu_hf,"Pz",kFullCircle,kBlue,kSolid,-1,kNone,0,0.7);

  TF1 *f1 = new TF1("f1","[0]*pow(1-2*15.*cosh(x)/13600.,[1])",-5.2,5.2);
  f1->SetParameters(0.20,10);
  f1->Draw("SAME");

  TLegend *leg = tdrLeg(0.35,0.80-0.05*3,0.60,0.80);
  if (isMC) {
    leg->AddEntry(hlowpu,Form("No PFHC,#LT#mu#GT=%1.1f",hmu_lowpu->GetMean()),"PLE");
    leg->AddEntry(hhighpu,Form("Baseline, #LT#mu#GT=%1.1f",hmu_highpu->GetMean()),"F");

  }
  else {
    leg->AddEntry(hlowpu,Form("Low PU, %1.2g nb^{-1}, #LT#mu#GT=%1.1f",lumlowpu/1e6,hmu_lowpu->GetMean()),"PLE");
    leg->AddEntry(hhighpu,Form("High PU, %1.2g nb^{-1}, #LT#mu#GT=%1.1f",lumhighpu/1e6,hmu_highpu->GetMean()),"F");
  }
  leg->AddEntry(f1,"Naive empirical model","L");
  
  c1->SaveAs("pdf/drawEtaSpikes/drawEtaSpikes.pdf");
  
} // drawEtaSpikes



void drawEtaSpikesJER() {

  TDirectory *curdir = gDirectory;

  TFile *foff = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer24MC_Flat22_NoPFH_v161.root","READ");
  //TFile *foff = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer22MC_hadCalibOff_1M_v161.root","READ");
  assert(foff && !foff->IsZombie());

  TFile *fon = new TFile("rootfiles/Prompt2026/Jet_v158/jmenano_mc_out_Summer24MC_Flat_JMENANO_v158.root","READ");
  //TFile *fon = new TFile("rootfiles/Prompt2026/Jet_v161/jmenano_mc_out_Summer22MC_base_1M_v161.root","READ");
  assert(fon && !fon->IsZombie());

  lumi_136TeV = "Summer24MC PFHC on/off";
  TH1D *h = tdrHist("h","Area normalized probability",0,0.13,
		    "p_{T,reco}/p_{T,gen}",0.4,1.6);
  TCanvas *c1 = tdrCanvas("c1jer",h,8,11,kSquare);

  // Also "", " Match" and _raw variants (what did they do?)
  TH3D *h3on = (TH3D*)fon->Get("HLT_MC/MCtruth/Response3D");
  assert(h3on);
  TH3D *h3off = (TH3D*)foff->Get("HLT_MC/MCtruth/Response3D");
  assert(h3off);

  int ieta = h3on->GetXaxis()->FindBin(2.7);
  double eta1 = h3on->GetXaxis()->GetBinLowEdge(ieta);
  double eta2 = h3on->GetXaxis()->GetBinLowEdge(ieta+1);
  int ipt = h3on->GetYaxis()->FindBin(350.);
  double pt1 = h3on->GetYaxis()->GetBinLowEdge(ipt);
  double pt2 = h3on->GetYaxis()->GetBinLowEdge(ipt+1);

  h3on->GetXaxis()->SetRange(ieta,ieta);
  h3on->GetYaxis()->SetRange(ipt,ipt);
  TH1D *h1on = h3on->ProjectionZ("h1on",ieta,ieta,ipt,ipt);
  h3off->GetXaxis()->SetRange(ieta,ieta);
  h3off->GetYaxis()->SetRange(ipt,ipt);
  TH1D *h1off = h3off->ProjectionZ("h1off",ieta,ieta,ipt,ipt);

  h1on->Scale(1./h1on->Integral());
  h1off->Scale(1./h1off->Integral());

  tdrDraw(h1on,"HIST",kNone,kBlack,kSolid,-1,kNone,0);
  tdrDraw(h1off,"Pz",kFullCircle,kRed,kSolid,-1,kNone,0,0.7);

  TLegend *leg = tdrLeg(0.63,0.90-0.05*4,0.88,0.90);
  leg->AddEntry(h1off,"PFHC off","PLE");
  leg->AddEntry(h1off,Form("RMS=%0.3f",h1off->GetRMS()),"");
  leg->AddEntry(h1on,"PFHC on","F");
  leg->AddEntry(h1on,Form("RMS=%0.3f",h1on->GetRMS()),"");

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.18,0.75,Form("%1.3f<eta<%1.3f",eta1,eta2));
  tex->DrawLatex(0.18,0.70,Form("%1.0f<p_{T}<%1.0f GeV",pt1,pt2));
  
  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawEtaSpikes/drawEtaSpikesJER.pdf");
} // void drawEtaSpikesJER
