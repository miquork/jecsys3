// hcal_depth_Run3.C — Run 3 HB+HE depth mapper from DN-2020/037 (Sunanda)
// Load:  .L hcal_depth_Run3.C+g
// Use :  test_hcal_depth();            // geometry stack
//        integrate_hcal_depth(50.0);   // per-depth fractions vs |ieta| (hist steps)
//
// Fixes vs v1:
//  • Clamp HE layer ranges to 0..17 (no index 18!).
//  • Validate segments (finite, y1>y0) before drawing / using in integrals.
//  • |ieta|=16 HE-slice now sits above ECAL (correct y0/y1).
//  • Fraction plot uses TH1D “hist” step lines (no jagged spikes).
//
// Numbers come from the DN tables (HB Table 1; HE Table 2). ECAL is a fixed λI offset.
// Tracker ignored as requested.

#ifndef HCAL_DEPTH_RUN3_C
#define HCAL_DEPTH_RUN3_C

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <cassert>

// ---------------- configuration ----------------
static double ECAL_LAMBDA = 1.10;     // ECAL contribution, fixed offset (λI)
static bool   DISABLE_HE_D1_AT_18 = true; // for plotting convenience

// ======== Raw absorber data (λI) from DN-2020/037 ========
// HB: Table 1 — per |ieta| row: {bL0, bL1, b2_9, b10_15, bL16} (HO columns dropped)
struct HBRow { double b0, b1, b2_9, b10_15, b16; };
static const std::array<HBRow,15> HB_TAB = {{
  /*1*/  {0.608,0.471,0.319,0.356,0.571},
  /*2*/  {0.732,0.472,0.322,0.359,0.582},
  /*3*/  {0.712,0.474,0.326,0.364,0.575},
  /*4*/  {0.751,0.490,0.334,0.372,0.612},
  /*5*/  {0.678,0.505,0.343,0.383,0.679},
  /*6*/  {0.824,0.454,0.356,0.397,0.722},
  /*7*/  {0.826,0.549,0.371,0.414,0.754},
  /*8*/  {0.849,0.543,0.389,0.434,0.717},
  /*9*/  {0.856,0.532,0.410,0.457,0.822},
  /*10*/ {0.952,0.606,0.434,0.484,0.823},
  /*11*/ {0.939,0.676,0.461,0.514,0.905},
  /*12*/ {1.023,0.700,0.492,0.549,0.929},
  /*13*/ {1.024,0.763,0.527,0.588,1.258},
  /*14*/ {1.093,0.795,0.565,0.630,1.647},
  /*15*/ {1.140,0.848,0.608,0.639,2.411}
}}; // DN-2020/037, Table 1

// HE: Table 2 — per |ieta| row: {bL0, bL1, b2_17}, for |ieta|=18..29.
struct HERow { double b0, b1, b2_17; };
static const std::array<HERow,12> HE_TAB = {{
  /*18*/ {0.363,0.491,0.533},
  /*19*/ {0.853,0.526,0.525},
  /*20*/ {0.859,0.525,0.519},
  /*21*/ {0.904,0.521,0.513},
  /*22*/ {0.894,0.516,0.508},
  /*23*/ {0.857,0.511,0.504},
  /*24*/ {0.901,0.506,0.499},
  /*25*/ {0.897,0.502,0.496},
  /*26*/ {0.882,0.496,0.493},
  /*27*/ {0.846,0.495,0.491},
  /*28*/ {0.640,0.492,0.488},
  /*29*/ {0.657,0.494,0.488}
}}; // DN-2020/037, Table 2

// ======== Build “front of layer” arrays (λI from ECAL backface) ========
// HB layers 0..16 (17 layers). Reuse row 15 for ai=16 (HB part).
inline std::vector<double> HB_fronts_from_row(int ai){
  const HBRow row = HB_TAB[ std::min(ai,15)-1 ];
  std::vector<double> F(17,0.0);
  double x = ECAL_LAMBDA;
  x += row.b0; F[0]=x;
  x += row.b1; F[1]=x;
  for(int l=2; l<=9;  ++l){ x += row.b2_9;   F[l]=x; }
  for(int l=10;l<=15; ++l){ x += row.b10_15; F[l]=x; }
  x += row.b16; F[16]=x;
  return F;
}

// HE layers 0..17 (18 layers).
inline std::vector<double> HE_fronts_from_row(int ai){
  const HERow row = HE_TAB[ai-18];
  std::vector<double> F(18,0.0);
  double x = ECAL_LAMBDA;
  x += row.b0; F[0]=x;
  x += row.b1; F[1]=x;
  for(int l=2; l<=17; ++l){ x += row.b2_17; F[l]=x; }
  return F;
}

// ======== Depth layout ========
inline int nDepths_at_ieta(int ai){
  if (ai>=1 && ai<=16)  return 4;
  if (ai==17)           return 2;   // only D2,D3 in HE
  if (ai==18)           return 5;
  if (ai>=19 && ai<=25) return 6;
  if (ai>=26 && ai<=28) return 7;
  if (ai==29)           return 3;
  return 0;
}

// HB layer lists per depth
inline std::vector<int> HB_layers_for_depth(int ai, int depth){
  if (ai<=15){
    if (depth==1) return {0};
    if (depth==2) return {1,2,3,4};
    if (depth==3) return {5,6,7,8,9};
    if (depth==4){
      if (ai==14 || ai==15) return {10,11,12,13,14,15}; // miss 16
      return {10,11,12,13,14,15,16};
    }
  } else if (ai==16){
    if (depth==1) return {1};
    if (depth==2) return {2,3,4};
    if (depth==3) return {5,6,7,8,9};
    if (depth==4) return {}; // provided by HE slice
  }
  return {};
}

// HE borders per your b5/b6/b7/b29; |iη|=17 uses only D2,D3 with borders {0,7,12}.
inline const std::vector<int>& HE_borders_for_ai(int ai){
  static const std::vector<int> b17 = {0,7,12};                 // → D2:[0..7), D3:[7..12)
  static const std::vector<int> b5  = {0,1,4,8,11,14};          // 18
  static const std::vector<int> b6  = {0,1,4,7,10,14,19};       // 19..25
  static const std::vector<int> b7  = {0,1,3,5,7,10,14,19};     // 26..28
  static const std::vector<int> b29 = {0,1,3,5};                // 29
  if (ai==17) return b17;
  if (ai==18) return b5;
  if (ai>=19 && ai<=25) return b6;
  if (ai>=26 && ai<=28) return b7;
  return b29; // 29
}

// Clamp borders to valid HE layer index range [0,18) when building lists
inline std::vector<int> HE_layers_for_depth(int ai, int dIndex){
  const auto& b = HE_borders_for_ai(ai);
  const int nD = (int)b.size()-1;
  const int nLay = 18; // 0..17
  if (dIndex<1 || dIndex>nD) return {};
  int L0 = std::clamp(b[dIndex-1], 0, nLay);
  int L1 = std::clamp(b[dIndex],   0, nLay);
  if (L1<=L0) return {};
  std::vector<int> v; v.reserve(L1-L0);
  for (int l=L0; l<L1; ++l) v.push_back(l);
  return v;
}
inline int HE_depth_label(int ai, int dIndex){ return (ai==17 ? dIndex+1 : dIndex); }

// ======== Public API ========
struct HcalSeg {
  double front{NAN}, back{NAN}, center{NAN}, thick{NAN}; // λI from ECAL backface
  int    region{0};   // 0=HB, 1=HE
  int    label{0};    // 1..7 (at |iη|=17 → 2,3)
  bool   valid{false};
};

inline HcalSeg hcal_segment_Run3(int ieta, int depth)
{
  const int ai = std::abs(ieta);
  HcalSeg s;

  // ---------- HB ----------
  if (ai>=1 && ai<=16){
    if (depth<1 || depth>4) return s;
    const auto F = HB_fronts_from_row(std::min(ai,15));

    // |ieta|=16, D4 from HE (front slice of |ieta|=18)
    if (ai==16 && depth==4){
      const auto Fhe = HE_fronts_from_row(18);
      const double y0 = Fhe[0];              // layer 0 front
      const double y1 = Fhe[3];              // up to front of layer 3
      if (std::isfinite(y0) && std::isfinite(y1) && y1>y0){
        s.front=y0; s.back=y1; s.center=0.5*(y0+y1); s.thick=y1-y0;
        s.region=1; s.label=4; s.valid=true;
      }
      return s;
    }

    const auto layers = HB_layers_for_depth(ai,depth);
    if (layers.empty()) return s;

    const int l0 = layers.front();
    double y0 = F[l0], y1 = 0.0;

    if      (depth==1) y1 = F[1];
    else if (depth==2) y1 = F[5];
    else if (depth==3) y1 = F[10];
    else { // D4 → stop at front of last layer in this group
      const int last = layers.back();
      y1 = F[last];
    }

    if (std::isfinite(y0) && std::isfinite(y1) && y1>y0){
      s.front=y0; s.back=y1; s.center=0.5*(y0+y1); s.thick=y1-y0;
      s.region=0; s.label=depth; s.valid=true;
    }
    return s;
  }

  // ---------- HE ----------
  if (ai>=17 && ai<=29){
    const int nD = nDepths_at_ieta(ai);
    if (depth<1 || depth>nD) return s;
    if (ai==18 && DISABLE_HE_D1_AT_18 && depth==1) return s;

    const auto F = HE_fronts_from_row( ai==17 ? 18 : ai );
    const auto layers = HE_layers_for_depth(ai, depth);
    if (layers.empty()) return s;

    const int l0 = layers.front();
    const int last = layers.back();
    const double y0 = F[l0];

    // y1 = front of first layer of the *next* depth; if none, front of last layer
    const auto nextLayers = HE_layers_for_depth(ai, depth+1);
    double y1 = nextLayers.empty() ? F[last] : F[nextLayers.front()];

    if (std::isfinite(y0) && std::isfinite(y1) && y1>y0){
      s.front=y0; s.back=y1; s.center=0.5*(y0+y1); s.thick=y1-y0;
      s.region=1; s.label=HE_depth_label(ai, depth); s.valid=true;
    }
    return s;
  }

  return s;
}

// Convenience wrapper
inline double hcal_depth_Run3(int ieta, int depth, double &thick){
  HcalSeg s = hcal_segment_Run3(ieta,depth);
  thick = s.thick; return s.center;
}
inline int hcal_nDepths_Run3(int ieta){
  const int ai = std::abs(ieta);
  if (ai>=1 && ai<=16) return 4;
  if (ai>=17 && ai<=29) return nDepths_at_ieta(ai);
  return 0;
}

// ---------------- plotting ----------------
#include "TCanvas.h"
#include "TBox.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TH1.h"
#include "TGraph.h"

inline int colorHB(int L){ static int c[]={ TColor::GetColor("#7fb3d5"),
                                           TColor::GetColor("#f5b041"),
                                           TColor::GetColor("#82e0aa"),
                                           TColor::GetColor("#f1948a") };
  return c[std::max(1,std::min(4,L))-1];
}
inline int colorHE(int L){ static int c[]={ TColor::GetColor("#7fb3d5"),
                                           TColor::GetColor("#f5b041"),
                                           TColor::GetColor("#82e0aa"),
                                           TColor::GetColor("#f1948a"),
                                           TColor::GetColor("#f7dc6f"),
                                           TColor::GetColor("#d2b48c"),
                                           TColor::GetColor("#b2babb") };
  return c[std::max(1,std::min(7,L))-1];
}

void test_hcal_depth(){
  const int ieta_min=1, ieta_max=29;

  double yMax = ECAL_LAMBDA;
  for (int i=ieta_min;i<=ieta_max;++i){
    const int nD = hcal_nDepths_Run3(i);
    for (int d=1; d<=nD; ++d){
      HcalSeg s = hcal_segment_Run3(i,d);
      if (s.valid && std::isfinite(s.back)) yMax = std::max(yMax, s.back);
    }
  }
  if (!std::isfinite(yMax) || yMax<ECAL_LAMBDA) yMax = ECAL_LAMBDA+5.0;
  yMax *= 1.05;

  TCanvas* c=new TCanvas("c_depth","HCAL stack (Run3 DN tables)",1200,620);
  TH1F* frame=(TH1F*)gPad->DrawFrame(ieta_min-0.5,0.0,ieta_max+0.5,yMax);
  frame->SetTitle(";|#eta_{i}|;Cumulative depth from ECAL backface  (#lambda_{I})");
  frame->GetYaxis()->SetTitleOffset(1.10);

  // ECAL band
  for (int i=ieta_min;i<=ieta_max;++i){
    TBox* b=new TBox(i-0.5,0.0,i+0.5,ECAL_LAMBDA);
    b->SetFillColor(TColor::GetColor("#e5e7e9")); b->SetLineColor(0); b->Draw("same");
  }

  // stacks (skip any invalid/degenerate segments)
  for (int i=ieta_min;i<=ieta_max;++i){
    const int nD = hcal_nDepths_Run3(i);
    for (int d=1; d<=nD; ++d){
      HcalSeg s = hcal_segment_Run3(i,d);
      if (!s.valid) continue;
      if (!(std::isfinite(s.front) && std::isfinite(s.back))) continue;
      if (s.back <= s.front) continue;
      TBox* b=new TBox(i-0.5, s.front, i+0.5, s.back);
      b->SetFillColor( (s.region==0) ? colorHB(s.label) : colorHE(s.label) );
      b->SetLineColor(kBlack); b->Draw("same");
    }
  }

  // separators
  TLine* hbhe=new TLine(16.5,0.0,16.5,yMax); hbhe->SetLineStyle(2); hbhe->SetLineColor(kGray+2); hbhe->Draw("same");
  TLine* hef = new TLine(29.5,0.0,29.5,yMax); hef->SetLineStyle(2); hef->SetLineColor(kGray+2); hef->Draw("same");

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
  tx.DrawLatex(0.19,0.92,"HB"); tx.DrawLatex(0.70,0.92,"HE");
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.14,0.85,Form("ECAL (constant): %.2f #lambda_{I}", ECAL_LAMBDA));
  tx.DrawLatex(0.14,0.80,"Run3 absorber spacings from DN-2020/037");

  c->SetGridx(); c->SetTicks(1,1); c->Update();
}

// --------------- toy longitudinal model (histogram style; normalized per |ieta|) ---------------
namespace prof {
  inline double gamma_pdf(double u,double k,double th){
    if(u<=0) return 0.0;
    return std::pow(u,k-1.0)*std::exp(-u/th)/( std::tgamma(k)*std::pow(th,k) );
  }
  struct Model { double fem=0.30; double kh=3.0, th=0.80; double ke=2.0, te=0.30; };
  inline double pdf(double u,const Model&m){
    return (1.0-m.fem)*gamma_pdf(u,m.kh,m.th) + m.fem*gamma_pdf(u,m.ke,m.te);
  }
  // Gauss-Legendre integrate on [a,b]
  inline double gauss16(double a,double b,const Model&m){
    static const double x[8]={0.0950125098376374,0.2816035507792589,0.4580167776572274,0.6178762444026438,
                              0.7554044083550030,0.8656312023878318,0.9445750230732326,0.9894009349916499};
    static const double w[8]={0.1894506104550685,0.1826034150449236,0.1691565193950025,0.1495959888165767,
                              0.1246289712555339,0.0951585116824928,0.0622535239386479,0.0271524594117541};
    if(b<=a) return 0.0; const double c=0.5*(b+a), h=0.5*(b-a); double s=0.0;
    for(int i=0;i<8;++i){ const double dx=h*x[i]; s+=w[i]*(pdf(c+dx,m)+pdf(c-dx,m)); }
    return s*h;
  }
  inline double integrate(double a,double b,const Model&m){
    a=std::max(0.0,a); b=std::max(b,a);
    const int N=std::max(1,int((b-a)/2.0)); double sum=0.0, x=a, step=(b-a)/N;
    for(int i=0;i<N;++i){ sum+=gauss16(x,x+step,m); x+=step; } return sum;
  }
}

void integrate_hcal_depth(double energy = 50.){
  using namespace prof; (void)energy; Model M;

  const int ieta_min=1, ieta_max=29, maxLbl=7;
  // one histogram per depth *label* (so HE@17 uses labels 2 and 3)
  std::vector<TH1D*> H(maxLbl+1,nullptr);
  for(int L=1; L<=maxLbl; ++L){
    H[L]=new TH1D(Form("hD%d",L), "", ieta_max, 0.5, ieta_max+0.5);
    H[L]->SetLineWidth(2);
    H[L]->SetLineColor( (L<=4)? colorHB(L) : colorHE(L) );
    H[L]->SetMarkerStyle(0);
  }
  TH1D* Hsum = new TH1D("hSum","", ieta_max, 0.5, ieta_max+0.5);
  Hsum->SetLineStyle(2); Hsum->SetLineColor(kBlack);

  for (int ieta=ieta_min; ieta<=ieta_max; ++ieta){
    const int nD = hcal_nDepths_Run3(ieta);
    std::vector<double> num(nD+1,0.0);
    std::vector<int>    lab(nD+1,0);
    double tot=0.0;

    for(int d=1; d<=nD; ++d){
      HcalSeg s = hcal_segment_Run3(ieta,d);
      if (!s.valid) continue;
      if (!(std::isfinite(s.front) && std::isfinite(s.back))) continue;
      if (s.back <= s.front) continue;
      const double u0 = s.front - ECAL_LAMBDA;
      const double u1 = s.back  - ECAL_LAMBDA;
      const double v  = integrate(u0,u1,M);
      num[d]=v; lab[d]=s.label; tot+=v;
    }
    if (tot<=0) continue;
    const double inv = 1.0/tot;
    double sum=0.0;
    for(int d=1; d<=nD; ++d){
      if (num[d]<=0) continue;
      const double f = num[d]*inv; sum+=f;
      H[ lab[d] ]->SetBinContent(ieta, f);
    }
    Hsum->SetBinContent(ieta, 1.0);
  }

  TCanvas* c=new TCanvas("c_profile","Depth fractions (normalized per |ieta|, hist style)",1100,650);
  TH1F* frame=(TH1F*)gPad->DrawFrame(0.5,0.0,29.5,1.05);
  frame->SetTitle(Form("HCAL depth fractions (Run3 DN tables; E=%.0f GeV);|#eta_{i}|;Fraction",energy));
  frame->GetYaxis()->SetTitleOffset(1.12);

  for(int L=1; L<=maxLbl; ++L) if (H[L]) H[L]->Draw("HIST SAME");
  Hsum->Draw("HIST SAME");
  c->SetGridx(); c->SetTicks(1,1); c->Update();
}

#endif // HCAL_DEPTH_RUN3_C
