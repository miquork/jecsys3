// shower_profile.C
// Run with: root -l -b -q shower_profile.C
// Produces shower_profile.png

#include "../tdrstyle_mod22.C"
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TMath.h>
#include <vector>

void shower_profile() {
  setTDRStyle();

  // Parameters
  double r_ecal = 0.04;  // X0 / lambda for PbWO4 ~0.04
  double r_dead = 0.1;   // ~ for brass/steel
  double r_hcal = 0.1;   // ~ for brass/steel
  int ieta = 1;          // Start with ieta=1 (~eta=0)
  //double dead_lambda = 0.608;  // Adjustable dead material in lambda (default for ieta=1)
  double dead_lambda = 0.608-0.333;  // Adjustable dead material in lambda (default for ieta=1)
  double ecal_lambda = 0;//1.10;    // ECAL depth in lambda

  // HCAL incremental groups for ieta=1 (multiplied for multi-layer groups)
  //std::vector<double> inc_groups = {0.471, 0.319*8, 0.356*6, 0.571, 2.091, 1.150};
  std::vector<double> inc_groups = {0.608-dead_lambda, 0.471+0.319*3, 0.319*5, 0.356*6+0.571, 2.091+1.150};
  //std::vector<int> group_colors = {kBlue, kOrange, kGreen, kRed, kBrown, kYellow};
  std::vector<int> group_colors = {
    /* 1 */ TColor::GetColor("#7fb3d5"),
    /* 2 */ TColor::GetColor("#f5b041"),
    /* 3 */ TColor::GetColor("#82e0aa"),
    /* 4 */ TColor::GetColor("#f1948a"),
    /* HO */ kBrown,
    /* beyond */ kYellow
  };

  // Can add other ieta here, e.g.:
  // if (ieta == 2) { dead_lambda = 0.732; inc_groups = {0.472, 0.322*8, 0.359*6, 0.582, 2.008, 1.151}; }
  // Similarly for others from Table 1, adjusting for ieta>4 without last layer.

  // Compute total lambda
  double total_lambda = ecal_lambda + dead_lambda;
  for (auto dl : inc_groups) total_lambda += dl;

  // Shower profile: dE/dlambda ~ pow(x, alpha-1) * exp(-x/beta) normalized to 1
  //double alpha = 2.0;  // Example parameters (adjustable)
  //double beta = 1.5;
  double alpha = 2.2;  // Optimized for 50 GeV hadron
  double beta = 1.2; // Optimized for 50 GeV hadron
  double norm = 1.0 / (TMath::Gamma(alpha) * TMath::Power(beta, alpha));
  TF1 *profile = new TF1("profile", "[0]*pow(x,[1]-1)*exp(-x/[2])", 0, total_lambda);
  profile->SetParameters(norm, alpha, beta);

  // Regions
  std::vector<std::pair<double, double>> regions;
  std::vector<int> region_colors;
  // ECAL
  regions.push_back({0.0, ecal_lambda});
  region_colors.push_back(kBlue);
  // Dead
  double dead_start = ecal_lambda;
  double dead_end = dead_start + dead_lambda;
  regions.push_back({dead_start, dead_end});
  region_colors.push_back(kGray);
  // HCAL groups
  double cum_hcal = dead_end;
  for (size_t i = 0; i < inc_groups.size(); ++i) {
    double group_start = cum_hcal;
    double group_end = group_start + inc_groups[i];
    regions.push_back({group_start, group_end});
    region_colors.push_back(group_colors[i]);
    cum_hcal = group_end;
  }

  // Background hist and canvas in tdr style
  extraText = "Private";
  lumi_136TeV = "Toy simulation";
  //TH1D *h = tdrHist("h","dE / d#lambda_{INT} (arb. units, norm. to 1)",0,profile->GetMaximum()*1.2,"Depth [#lambda_{INT}]",0,total_lambda);
  TH1D *h = tdrHist("h","dE / d#lambda_{INT} (norm. to 1)",0,0.35,"Depth [#lambda_{INT}]",0,total_lambda);
  TCanvas *c = tdrCanvas("c",h,8,11,kRectangular);//kSquare);
  c->cd();

  // Fill colored areas under curve
  int npoints = 200;  // For smooth curve
  double max_y = profile->GetMaximum();
  double min_width = 0.5;  // For adding text
  for (size_t i = 0; i < regions.size(); ++i) {
    double x1 = regions[i].first;
    double x2 = regions[i].second;
    double dx = (x2 - x1) / npoints;
    TGraph *g = new TGraph(2 * npoints + 2);
    for (int j = 0; j <= npoints; ++j) {
      double x = x1 + j * dx;
      double y = profile->Eval(x);
      g->SetPoint(j, x, y);
      g->SetPoint(2 * npoints + 1 - j, x, 0.0);
    }
    g->SetFillColor(region_colors[i]);
    g->SetFillStyle(1001);
    g->Draw("f");

    // Vertical boundary lines (except at 0)
    if (x1 > 0) {
      TLine *l = new TLine(x1, 0, x1, 0.35);//max_y * 1.2);
      l->SetLineStyle(kDotted);
      l->SetLineColor(kGray+1);//kBlack);
      l->Draw();
    }
  }
  for (size_t i = 0; i < regions.size(); ++i) {
    double x1 = regions[i].first;
    double x2 = regions[i].second;
    // Fraction text if wide enough (always for dead)
    double fraction = profile->Integral(x1, x2);
    //if ((x2 - x1 > min_width) || (i == 1)) {  // i==1 is dead
    if (fraction>0) {
      //TLatex *tex = new TLatex(max(0.5, (x1 + x2) / 2.0), (i==1 ? 0.03 : 0) + max_y * 0.5 - 0.015*i, (i == 1) ? Form("%.1f%% lost", fraction * 100) : Form("%.1f%%", fraction * 100));
      TLatex *tex = new TLatex(max(0.5, (x1 + x2) / 2.0), max_y * 0.5 - 0.015*i, (i==0 ? Form("%.1f%% ECAL", fraction*100) : (i == 1) ? Form("%.1f%% lost", fraction * 100) : Form("%.1f%%", fraction * 100)));
      tex->SetTextAlign(22);
      tex->SetTextSize(0.03);
      tex->SetTextColor(kBlack);//kWhite);
      tex->Draw();
    }
  }

  // Draw profile line on top
  // tdrDraw(profile, "L", -1, kBlack, kSolid, kBlack);
  profile->Draw("SAME");
  
  // Second axis: X0 on top (piecewise due to different ratios)
  double x0_ecal = ecal_lambda / r_ecal;
  double x0_dead = dead_lambda / r_dead;
  double x0_hcal = (total_lambda - ecal_lambda - dead_lambda) / r_hcal;
  double total_x0 = x0_ecal + x0_dead + x0_hcal;
  double post_ecal_start_lambda = ecal_lambda;
  double post_ecal_start_x0 = x0_ecal;

  // Axis for ECAL part
  //TGaxis *ax_ecal = new TGaxis(0, h->GetMaximum(), ecal_lambda, h->GetMaximum(), 0, x0_ecal, 505, "-");
  TGaxis *ax_ecal = new TGaxis(0, h->GetMaximum(), ecal_lambda, h->GetMaximum(), 0, x0_ecal, 501, "-");
  ax_ecal->SetLineColor(kBlue);
  ax_ecal->SetLabelColor(kBlue);
  ax_ecal->SetTitleColor(kBlue);
  ax_ecal->SetLabelSize(0.03);
  ax_ecal->SetTitle("Depth [X_{0}]");
  ax_ecal->SetTitleOffset(0.8);
  ax_ecal->Draw();

  // Axis for post-ECAL (dead + HCAL, same r)
  TGaxis *ax_post = new TGaxis(post_ecal_start_lambda, h->GetMaximum(), total_lambda, h->GetMaximum(), post_ecal_start_x0, total_x0, 510, "+");
  ax_post->SetLineColor(kBlue);
  ax_post->SetLabelColor(kBlue);
  ax_post->SetTitleColor(kBlue);
  ax_post->SetLabelSize(0.03);
  ax_post->Draw();

  gPad->RedrawAxis();

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  if (ecal_lambda==0) tex->DrawLatex(0.44,0.60,"50 GeV H-hadron");
  if (ecal_lambda!=0) tex->DrawLatex(0.44,0.60,"50 GeV EH-hadron");
  
  if (ecal_lambda==0) c->SaveAs("pdf/shower_profile/shower_profile_H.pdf");
  if (ecal_lambda!=0) c->SaveAs("pdf/shower_profile/shower_profile_EH.pdf");
}
