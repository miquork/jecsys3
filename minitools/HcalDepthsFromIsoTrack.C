// HcalDepthsFromIsoTrack.C
//
// Purpose: This macro derives HCAL depth intercalibrations using IsoTrack samples
// from provided ROOT files containing TProfile2D histograms of fractional responses
// per depth, eta, and pt.
//
// Input files:
// - ../isotrack/rootfiles/IsoTrack_25Aug20_EG25C2C_orig.root (data)
// - rootfiles/IsoTrack_25Sep04_25MCNOPU_orig.root (MC)
//
// Key histograms:
// - p2c for N=0..8
// - missing for HF
// - pc for total response
//
// Computation:
// - For each depth (0=ECAL,1-7=HCAL,8=HO,9=HF long,10=HF short)
// - Project to depth (ProfileY) for mean fractions in data and MC
// - Correction = (MC fraction / data fraction) * k_FSR_ratio (default 0.99)
// - Stat unc from propagation
// - Syst unc from pt dependence (low/high pt split) + k_FSR_syst (0.01)
// - Skip invalid bins (fraction <=1e-5 or error=0)
//
// Outputs:
// - Summary plots vs eta: Frac_Data_All, Frac_MC_All, Ratio_All, Corr_All.pdf
// - Plots vs pt per depth: FracRatioVsPt_Depth*.pdf
// - RespDB vs pt: RespDB_Data_VsPt, RespDB_MC_VsPt, RespDB_Ratio_VsPt.pdf (all ieta)
// - Text: HcalRespCorrs_IsoTrack_2025.txt (conditions format, skip depth0)
// - Summary: HcalDepthsFromIsoTrack_Summary.txt (ieta,depth,corr,stat,syst)
// - ROOT: rootfiles/HcalDepthFromIsoTrack.root with TH1D h_depth_<N> (corr vs ieta, error=stat)
//
// Usage: root -l -b -q minitools/HcalDepthsFromIsoTrack.C
// To continue: Load output ROOT, access h_depth_* for corrections per depth/ieta
// Compare to IsoTrack method, upload text to DB, etc.

#include <TFile.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TH1D.h>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>

#include "../tdrstyle_mod22.C"

/*
bool isValidEta(int ieta, int depth) {
  int absi = abs(ieta);
  if (depth==0 && absi>29) return false;
  if (depth==1 && (absi==17 || absi==18 || absi>29)) return false;
  if (depth==2 && absi>29) return false;
  if (depth==3 && absi>29) return false;
  if (depth==4 && (absi==17 || absi>28)) return false;
  if (depth==5 && (absi<18 || absi>28)) return false;
  if (depth==6 && (absi<19 || absi>28)) return false;
  if (depth==7 && (absi<26 || absi>28)) return false;
  if (depth==8 && (absi>16)) return false;
  if (depth==9 && (absi<29)) return false;
  if (depth==10 && (absi<29)) return false;
  return true;
}
*/
bool isValidEta(int ieta, int depth) {
  int absi = abs(ieta);
  if (depth==0 && absi<=29) return true;
  if (depth==1 && (absi<=16 || (absi>=19 && absi<=29))) return true;
  if (depth==2 && absi<=29) return true;
  if (depth==3 && absi<=29) return true;
  if (depth==4 && (absi<=16 || (absi>=18 && absi<=28))) return true;
  if (depth==5 && (absi>=18 && absi<=28)) return true;
  if (depth==6 && (absi>=19 && absi<=28)) return true;
  if (depth==7 && (absi>=26 && absi<=28)) return true;
  if (depth==8 && (absi<=15)) return true;
  if (depth==9 && (absi>=29)) return true;
  if (depth==10 && (absi>=29)) return true;
  return false;
}

void HcalDepthsFromIsoTrack() {
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  writeExtraText = true;
  extraText = "Private";
  lumi_136TeV = "IsoTrack 40-60 GeV, EGamma 2025C";
  //lumi_136TeV = "Summer24 + EGamma 2025Cv2";

  // Configurable parameters
  double ptmin = 55;
  double ptmid = 110;//78;
  double ptmax = 230;//110;
  double k_FSR_ratio = 0.99;  // k_FSR_data / k_FSR_MC
  double k_FSR_syst = 0.01;   // Fractional syst uncertainty on k_FSR_ratio

  // Depths (0=ECAL, 1-7=HCAL, 8=HO, 9=hfLong, 10=hfShort)
  int num_depths = 9;
  std::string depth_str[11] = {"0", "1", "2", "3", "4", "5", "6", "7", "8"};// "hfLong", "hfShort"};
  const char* label[11] = {"ECAL","Depth 1","Depth 2","Depth 3","Depth 4","5 (HE)","6 (HE)","7 (HE)","HO"};//,"HF Long","HF Short"};
  int colors[11] = {kCyan+1, kBlue, kOrange+1, kGreen+1, kRed, kYellow+1, kOrange+3, kGray+2, kGray+1};//, kBlue, kRed};
  int markers[11] = {20, 21, 22, 23, 24, 25, 26, 27, 28};//, 29, 30};

  // Hist prefixes
  std::string hist_prefix = "p2c";
  std::string hf_prefix = "";//hp2_tnpHFLongShort_TagPhotonCand_pt_vs_ProbeJet0_eta_vs_fracRespDB";
  std::string resp_prefix = "pc";

  // Output directories
  std::string plot_dir = "pdf/HcalDepthsFromIsoTrack/";
  std::string text_dir = "txt/HcalDepthsFromIsoTrack/";
  gSystem->mkdir(plot_dir.c_str(), kTRUE);
  gSystem->mkdir(text_dir.c_str(), kTRUE);

  // Open files
  //TFile *f_data = new TFile("../isotrack/rootfiles/IsoTrack_25Aug20_EG25C2C_orig.root","READ");
  //TFile *f_data = new TFile("rootfiles/IsoTrack_25Sep04_EG25C_orig.root","READ");
  //TFile *f_data = new TFile("rootfiles/IsoTrack_25Sep04_EG25C2C_orig.root","READ");
  TFile *f_data = new TFile("rootfiles/IsoTrack_25Sep04_EG25CC_orig.root","READ");
  //TFile *f_mc = new TFile("rootfiles/IsoTrack_25Sep04_25MCNOPU_orig.root","READ");
  TFile *f_mc = new TFile("rootfiles/IsoTrack_25Sep04_24MCPU_orig.root","READ");
  if (!f_data || f_data->IsZombie() || !f_mc || f_mc->IsZombie()) {
    return;
  }
  curdir->cd();
  
  // Define common ieta/eta binning for HB+HE and HF
  double veta[] = { // eta edges
    -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191
  };
  const int n_veta = sizeof(veta)/sizeof(veta[0])-1;
  
  double etas[] = { // eta centers
    -5.04, -4.8025, -4.627, -4.4505, -4.277, -4.102, -3.926, -3.7515, -3.5765, -3.4015, -3.2265, -3.0515, -2.9085, -2.7515, -2.575, -2.411, -2.247, -2.1075, -1.9865, -1.88, -1.785, -1.6965, -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, -1.1745, -1.0875, -1.0005, -0.918, -0.831, -0.7395, -0.6525, -0.5655, -0.4785, -0.3915, -0.3045, -0.2175, -0.1305, -0.0435, 0.0435, 0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 0.5655, 0.6525, 0.7395, 0.831, 0.918, 1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 1.4355, 1.5225, 1.6095, 1.6965, 1.785, 1.88, 1.9865, 2.1075, 2.247, 2.411, 2.575, 2.7515, 2.9085, 3.0515, 3.2265, 3.4015, 3.5765, 3.7515, 3.926, 4.102, 4.277, 4.4505, 4.627, 4.8025, 5.04
  };
  const int n_eta = sizeof(etas)/sizeof(etas[0]);
  assert(n_eta==n_veta);

  int ietas[] = { // ieta
    -41, -40, -39, -38, -37, -36, -35, -34, -33, -32, -31, -30, -29, -28, -27, -26, -25, -24, -23, -22, -21, -20, -19, -18, -17, -16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41
  };
  const int n_ieta = sizeof(ietas)/sizeof(ietas[0]);
  assert(n_eta==n_ieta);

  double vieta[n_ieta+1];
  //for (int i = 0; i != n_ieta+2; ++i) { vieta[i] = ietas[i]; }
  for (int i = 0; i != n_ieta/2; ++i) { vieta[i] = ietas[i]-0.5; }
  vieta[n_ieta/2] = 0.;
  for (int i = n_ieta/2; i != n_ieta+1; ++i) { vieta[i+1] = ietas[i]+0.5; }
  
  // Storage for corrs and uncs (depth x eta)
  std::vector<std::vector<double>> corrs(num_depths, std::vector<double>(n_eta, 1.0));
  std::vector<std::vector<double>> stat_unc(num_depths, std::vector<double>(n_eta, 0.0));
  std::vector<std::vector<double>> syst_unc(num_depths, std::vector<double>(n_eta, 0.0));

  // For summary plots
  std::vector<TGraphErrors*> g_datas;
  std::vector<TGraphErrors*> g_mcs;
  std::vector<TGraphErrors*> g_ratios;
  std::vector<TH1D*> h_corrs;
  std::vector<TH1D*> h_fracs;

  // Loop over depths
  for (int d = 0; d < num_depths; ++d) {
    //std::string dstr = depth_str[d];
    //std::string hname = (d <= 8) ? hist_prefix + "Depth" + dstr : hf_prefix + dstr;
    std::string hname = hist_prefix;
    
    TProfile2D *hp_data = (TProfile2D*)f_data->Get(hname.c_str());
    TProfile2D *hp_mc = (TProfile2D*)f_mc->Get(hname.c_str());
    if (!hp_data || !hp_mc) continue;

    int ix1 = hp_data->GetYaxis()->FindBin(d);//ptmin);
    int ix2 = hp_data->GetYaxis()->FindBin(d);//ptmax-0.1);

    TProfile *proj_data = hp_data->ProfileX(Form("proj_data_%d",d),ix1,ix2);
    TProfile *proj_mc = hp_mc->ProfileX(Form("proj_mc_%d",d),ix1,ix2);

    // Syst from pt dependence: split pt into low/high halves
    int ix_mid = hp_data->GetXaxis()->FindBin(d);//ptmid);
    TProfile *proj_data_low = hp_data->ProfileY(Form("proj_data_low_%d",d), ix1, ix_mid);
    TProfile *proj_mc_low = hp_mc->ProfileY(Form("proj_mc_low_%d",d), ix1, ix_mid);
    TProfile *proj_data_high = hp_data->ProfileY(Form("proj_data_high_%d",d), ix_mid + 1, ix2);
    TProfile *proj_mc_high = hp_mc->ProfileY(Form("proj_mc_high_%d",d), ix_mid + 1, ix2);

    // Graphs for plots
    TGraphErrors *g_data = new TGraphErrors();
    TGraphErrors *g_mc = new TGraphErrors();
    TGraphErrors *g_ratio = new TGraphErrors();
    TH1D *h_corr = new TH1D(Form("h_depth_%d",d),";Correction Factor;i#eta",n_eta,vieta);
    TH1D *h_frac = new TH1D(Form("h_frac_%d",d),";Energy Fraction;i#eta",n_eta,vieta);
    int np = 0;

    // Loop over eta bins
    for (int iy = 1; iy <= n_eta; ++iy) {
      int iy_hist = proj_data->GetXaxis()->FindBin(ietas[iy-1]);
      double f_data = proj_data->GetBinContent(iy_hist);
      double e_data = proj_data->GetBinError(iy_hist);
      double f_mc = proj_mc->GetBinContent(iy_hist);
      double e_mc = proj_mc->GetBinError(iy_hist);

      if (f_data <= 1e-5 || f_mc <= 1e-5 || e_data <= 0 || e_mc <= 0) continue;  // Invalid or no data

      // Update FSR bias correction to use RespDB
      double ratio = f_data / f_mc;
      double corr = f_mc / f_data * k_FSR_ratio;
      double stat = corr * sqrt( pow(e_mc/f_mc, 2) + pow(e_data/f_data, 2) );

      double f_data_low = proj_data_low->GetBinContent(iy_hist);
      double f_mc_low = proj_mc_low->GetBinContent(iy_hist);
      double corr_low = (f_data_low > 1e-5 && f_mc_low > 1e-5) ? (f_mc_low / f_data_low) * k_FSR_ratio : corr;

      double f_data_high = proj_data_high->GetBinContent(iy_hist);
      double f_mc_high = proj_mc_high->GetBinContent(iy_hist);
      double corr_high = (f_data_high > 1e-5 && f_mc_high > 1e-5) ? (f_mc_high / f_data_high) * k_FSR_ratio : corr;

      double syst_pt = fabs(corr_low - corr_high) / 2.0;
      double syst = sqrt( pow(syst_pt, 2) + pow(k_FSR_syst * corr, 2) );

      // Syst from excluding ECAL
      
      // Store
      int idx = iy - 1;
      corrs[d][idx] = corr;
      stat_unc[d][idx] = stat;
      syst_unc[d][idx] = syst;

      // Fill graphs
      int absi = abs(ietas[idx]);
      if ((d<=4 && absi<=29) || (d>=5 && d<=6 && absi>=18 && absi<=29) ||
	  (d==7 && absi>=26 && absi<=29) ||
	  (d==8 && absi<=16) || (d>=9 && absi>=29)) {
	g_data->SetPoint(np, etas[idx], f_data);
	g_data->SetPointError(np, 0, e_data);
	g_mc->SetPoint(np, etas[idx], f_mc);
	g_mc->SetPointError(np, 0, e_mc);
	g_ratio->SetPoint(np, etas[idx], ratio);
	g_ratio->SetPointError(np, 0, ratio * sqrt( pow(e_mc/f_mc, 2) + pow(e_data/f_data, 2) ));
	np++;
      }
      if (isValidEta(ietas[idx],d)) {
	h_corr->SetBinContent(iy, corr);
	h_corr->SetBinError(iy, stat);
	h_frac->SetBinContent(iy, f_data);
	h_frac->SetBinError(iy, e_data);
      }
    }

    // Save graphs for summary
    g_datas.push_back(g_data);
    g_mcs.push_back(g_mc);
    g_ratios.push_back(g_ratio);
    h_corrs.push_back(h_corr);
    h_fracs.push_back(h_frac);
  }

  // Common plotting tools
  TLine *l = new TLine();
  l->SetLineStyle(kDashed); l->SetLineColor(kGray+1);
  TLine *l2 = new TLine();
  l2->SetLineStyle(kDotted); l2->SetLineColor(kGray);
  
  // Summary Frac_Data_All
  double eta_min = veta[0];
  double eta_max = veta[n_eta];
  TH1D *h_eta_frac = tdrHist("h_eta_frac","Mean fraction",0,0.5,"#eta_{jet}",eta_min,eta_max);
  TCanvas *c_frac_data = tdrCanvas("frac_data_all",h_eta_frac,8,11,kRectangular);
  l->DrawLine(eta_min,1,eta_max,1);
  l2->DrawLine(eta_min,1.2,eta_max,1.2);
  l2->DrawLine(eta_min,0.8,eta_max,0.8);
  
  TLegend *leg_frac_data = tdrLeg(0.7,0.7,0.9,0.9);
  for (size_t i = 0; i < g_datas.size(); ++i) {
    auto g = g_datas[i];
    tdrDraw(g, "Pz", markers[i], colors[i]);
    //leg_frac_data->AddEntry(g, ("Depth " + depth_str[i]).c_str(), "lep");
    leg_frac_data->AddEntry(g, label[i], "lep");
  }
  leg_frac_data->Draw();
  c_frac_data->SaveAs((plot_dir + "HcalDepthsFromIsoTrack_Frac_Data_All.pdf").c_str());
  //delete c_frac_data;

  // Summary Frac_MC_All
  TCanvas *c_frac_mc = tdrCanvas("frac_mc_all",h_eta_frac,8,11,kRectangular);
  l->DrawLine(eta_min,1,eta_max,1);
  l2->DrawLine(eta_min,1.2,eta_max,1.2);
  l2->DrawLine(eta_min,0.8,eta_max,0.8);
  TLegend *leg_frac_mc = tdrLeg(0.7,0.7,0.9,0.9);
  for (size_t i = 0; i < g_mcs.size(); ++i) {
    auto g = g_mcs[i];
    tdrDraw(g, "Pz", markers[i], colors[i]);
    //leg_frac_mc->AddEntry(g, ("Depth " + depth_str[i]).c_str(), "lep");
    leg_frac_mc->AddEntry(g, label[i], "lep");
  }
  leg_frac_mc->Draw();
  c_frac_mc->SaveAs((plot_dir + "HcalDepthsFromIsoTrack_Frac_MC_All.pdf").c_str());
  //delete c_frac_mc;

  // Summary Ratio_All
  TH1D *h_eta_ratio = tdrHist("h_eta_ratio","Data/MC fraction",0.45,2.25,"#eta_{jet}",eta_min,eta_max);
  TCanvas *c_ratio_all = tdrCanvas("ratio_all",h_eta_ratio,8,11,kRectangular);
  l->DrawLine(eta_min,1,eta_max,1);
  l2->DrawLine(eta_min,1.2,eta_max,1.2);
  l2->DrawLine(eta_min,0.8,eta_max,0.8);
  TLegend *leg_ratio = tdrLeg(0.30,0.70,0.90,0.90);
  leg_ratio->SetNColumns(3);
  for (size_t i = 0; i < g_ratios.size(); ++i) {
    auto g = g_ratios[i];
    tdrDraw(g, "Pz", markers[i], colors[i]);
    //leg_ratio->AddEntry(g, ("Depth " + depth_str[i]).c_str(), "lep");
    leg_ratio->AddEntry(g, label[i], "lep");
  }
  leg_ratio->Draw();
  c_ratio_all->SaveAs((plot_dir + "HcalDepthsFromIsoTrack_Ratio_All.pdf").c_str());
  //delete c_ratio_all;

  // Summary Corr_All
  double ieta_min = vieta[0];
  double ieta_max = vieta[n_eta];
  double c_min(0.05), c_max(2.8), c_maxb(2.0);
  TH1D *h_eta_corr = tdrHist("h_eta_corr","Correction factor",c_min,c_max,"i#eta",ieta_min,ieta_max);
  TCanvas *c_all = tdrCanvas("corr_all",h_eta_corr,8,11,kRectangular);

  l->DrawLine(ieta_min,1,ieta_max,1);
  l2->DrawLine(ieta_min,1.2,ieta_max,1.2);
  l2->DrawLine(ieta_min,0.8,ieta_max,0.8);
  l2->DrawLine(-29.5,c_min,-29.5,c_maxb);
  l2->DrawLine(-25.5,c_min,-25.5,c_maxb);
  l2->DrawLine(-18.5,c_min,-18.5,c_maxb);
  l2->DrawLine(-16.5,c_min,-16.5,c_maxb);
  l2->DrawLine(+16.5,c_min,+16.5,c_maxb);
  l2->DrawLine(+18.5,c_min,+18.5,c_maxb);
  l2->DrawLine(+29.5,c_min,+29.5,c_maxb);
  l2->DrawLine(+25.5,c_min,+25.5,c_maxb);
    
  TLegend *leg_all = tdrLeg(0.30,0.70,0.90,0.90);
  leg_all->SetNColumns(3);

  for (size_t i = 0; i < h_corrs.size(); ++i) {
    auto h = h_corrs[i]; assert(h);
    tdrDraw(h, "Pz", markers[i], colors[i]);
    //leg_all->AddEntry(h, ("Depth " + depth_str[i]).c_str(), "lep");
    leg_all->AddEntry(h, label[i], "lep");
  }
  leg_all->Draw();
  c_all->SaveAs((plot_dir + "HcalDepthsFromIsoTrack_Corr_All.pdf").c_str());
  //delete c_all;


  // Output text file for HcalRespCorrs
  std::ofstream out(text_dir + "HcalRespCorrs_IsoTrack_2025C.txt");
  out << "HcalRespCorrs_IsoTrack_v1.0" << std::endl;  // Tag
  out << "1 infinite" << std::endl;  // IOV example

  for (int d = 1; d < num_depths; ++d) {  // Skip depth 0
    std::string dstr = depth_str[d];
    std::string subdet;
    int depth_out;
    int min_absieta = 1, max_absieta = 41;
    if (d <= 7) {
      depth_out = d;
    } else if (d == 8) {
      subdet = "HO";
      depth_out = 1;
      max_absieta = 15;
    } else {
      subdet = "HF";
      depth_out = (d == 9) ? 1 : 2;
      min_absieta = 29;
      max_absieta = 41;
    }

    for (int idx = 0; idx < n_eta; ++idx) {
      int ieta = ietas[idx];
      int absieta = TMath::Abs(ieta);
      if (absieta < min_absieta || absieta > max_absieta) continue;

      double corr = corrs[d][idx];
      if (corr == 1.0) continue;  // Not computed (invalid)

      if (d <= 7) {
        subdet = (absieta <= 16) ? "HB" : "HE";
      }

      for (int iphi = 1; iphi <= 72; ++iphi) {
        out << subdet << " " << ieta << " " << iphi << " " << depth_out << " " << corr << std::endl;
      }
    }
  }
  out.close();

  // Optional: Print summary table of corrs and uncs to console or separate file
  std::ofstream summary(text_dir + "HcalDepthsFromIsoTrack_Summary.txt");
  summary << "ieta\tdepth\tcorr\tstat_unc\tsyst_unc" << std::endl;
  for (int d = 0; d < num_depths; ++d) {  // Include depth 0
    for (int idx = 0; idx < n_eta; ++idx) {
      if (corrs[d][idx] != 1.0 || d == 0) {  // Include even if default for d=0
        summary << ietas[idx] << "\t" << depth_str[d] << "\t" << corrs[d][idx] << "\t" << stat_unc[d][idx] << "\t" << syst_unc[d][idx] << std::endl;
      }
    }
  }
  summary.close();

  // Store corrections as TH1D
  TFile *fout = new TFile("rootfiles/HcalDepthFromIsoTrack.root","RECREATE");
  for (int d = 0; d < num_depths; ++d) {
    //std::string hname = "h_depth_" + depth_str[d];
    //TH1D *h = new TH1D(hname.c_str(), ";i#eta;Correction", n_ieta, vieta);
    //for (int idx = 0; idx < n_eta; ++idx) {
    //int bin = h->FindBin(ietas[idx]);
    //h->SetBinContent(bin, corrs[d][idx]);
    //h->SetBinError(bin, stat_unc[d][idx]);
    //}
    TH1D *h = h_corrs[d]; assert(h);
    h->Write();
  }
  for (int d = 0; d < num_depths; ++d) {
    TH1D *h = h_fracs[d]; assert(h);
    h->Write();
  }
  fout->Close();

  // Cleanup
  f_data->Close();
  f_mc->Close();
}
