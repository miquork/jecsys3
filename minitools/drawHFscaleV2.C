// Purpose: 1) Draw HF scale from James Natoli, V2 update on March 27:
//          https://indico.cern.ch/event/1398985/contributions/5882918/attachments/2828010/4940971/HF_ZeeCalibrationRun3_Full.pdf?#page=7
//          2) Extrapolate to higher ieta (Zee bias from HF electron pT cut)
#include "TF1.h"
#include "TFile.h"

#include "../tdrstyle_mod22.C"

const bool debug = false;

void drawHFscaleV2() {

  setTDRStyle();
  
  // From James Vincent Natoli, 21 March 2024 (same in 27 March 2024 talk)
  // 2022 PreEE (Era BCD)
  //
  // Notes: HF timing shifted -3.5ns at end of 22B, leading to -20% HF scale?
  // "355657 HF timing was moved by 3.5ns towards earlier position" from
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalDPGRun3Updates#Time_re_alignments_AN1
  string v2022[] =
    {
      "-39: 1.34513",
      "-38: 1.68973",
      "-37: 1.65798",
      "-36: 1.64308",
      "-35: 1.56459",
      "-34: 1.48143",
      "-33: 1.32341",
      "-32: 1.23729",
      "-31: 1.26103",
      "-30: 1.34506",
      "30: 1.37317",
      "31: 1.25929",
      "32: 1.27125",
      "33: 1.35888",
      "34: 1.47288",
      "35: 1.57894",
      "36: 1.59043",
      "37: 1.65269",
      "38: 1.50635",
      "39: 1.15942"
    };
  const int n2022 = sizeof(v2022)/sizeof(v2022[0]);
  const int neta = n2022;
  
  // From James Vincent Natoli's talk on 27 March 2024
  // 2023[sic] PostEE* (Era EF), *exclude Era E[sic]
  // (probably means 2022EE eras EF and excluding era G)
  //
  // Notes: HF timing fixed +1.5ns middle of 22E, fixing -20% scale?
  // "360144 HF timing shifted (around +1.5ns) according to DPG corrections https://indico.cern.ch/event/1198215/#2-hf-time-alignment" from
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalDPGRun3Updates#Time_re_alignments_AN1
  string v2022ee[] =
    {
      "-39: 1.09325",
      "-38: 1.30552",
      "-37: 1.41048",
      "-36: 1.3743",
      "-35: 1.27036",
      "-34: 1.21606",
      "-33: 1.06164",
      "-32: 0.990691",
      "-31: 1.02687",
      "-30: 1.13321",
      "30: 1.16627",
      "31: 1.02887",
      "32: 1.02586",
      "33: 1.09824",
      "34: 1.21672",
      "35: 1.31477",
      "36: 1.32714",
      "37: 1.39904",
      "38: 1.2932",
      "39: 1.13874"
    };
  const int n2022ee = sizeof(v2022ee)/sizeof(v2022ee[0]);
  assert(n2022ee==neta);
  
  // From James Vincent Natoli's talk on 27 March 2024
  string v2023[] =
    {
      "-39: 1.02711",
      "-38: 1.4323",
      "-37: 1.55933",
      "-36: 1.49696",
      "-35: 1.37296",
      "-34: 1.28859",
      "-33: 1.12613",
      "-32: 1.05715",
      "-31: 1.07569",
      "-30: 1.03853",
      "30: 1.07464",
      "31: 1.07371",
      "32: 1.09449",
      "33: 1.1675",
      "34: 1.29517",
      "35: 1.40232",
      "36: 1.44118",
      "37: 1.47066",
      "38: 1.43343",
      "39: 1.02197"
    };
  const int n2023 = sizeof(v2023)/sizeof(v2023[0]);
  assert(n2023==neta);
  
  // From James Vincent Natoli, 21 March 2024 (same as in 27 March 2024 talk)
  // 2023 PostBPix (Era D)
  // [renamed to 2023bpix compared to drawHFscale.C V1)
  string v2023bpix[] =
    {
      "-39: 1.01025",
      "-38: 1.51389",
      "-37: 1.54939",
      "-36: 1.48511",
      "-35: 1.37398",
      "-34: 1.29287",
      "-33: 1.1278",
      "-32: 1.05079",
      "-31: 1.06414",
      "-30: 1.02493",
      "30: 1.07568",
      "31: 1.07161",
      "32: 1.09498",
      "33: 1.16741",
      "34: 1.2854",
      "35: 1.40581",
      "36: 1.43495",
      "37: 1.50745",
      "38: 1.44211",
      "39: 1.49512"
    };
  const int n2023bpix = sizeof(v2023bpix)/sizeof(v2023bpix[0]);
  assert(n2023bpix==neta);

  string eras[] =
    {"2022","2022ee","2023","2023bpix"};
  const int nera = sizeof(eras)/sizeof(eras[0]);

  map<string, string*> mera;
  mera["2022"] = v2022;
  mera["2022ee"] = v2022ee;
  mera["2023"] = v2023;
  mera["2023bpix"] = v2023bpix;
  
  // Regular L2Relative eta binning => HF only
  double vx[] =
    {
      //0, 0.087, 0.174, 0.261, 0.348, 0.435,
      //0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
      //1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
      //2.65, 2.853,
      2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
      4.363, 4.538, 4.716, 4.889, 5.191
    };
  const int nx = sizeof(vx)/sizeof(vx[0])-1;

  TH1D *hf3 = new TH1D("hf3",";|#eta|;HF scale",nx,vx);

  TH1D *hf = new TH1D("hf",";|#eta|;HF scale",nx,vx);
  map<string, TH1D*> mh, mhm, mhp, mhe;
  map<string, TF1*> mf1, mf2;

  double k22 = 1.18;
  double k22ee = 0.95;
  for (int iera = 0; iera != nera; ++iera) {

    string era = eras[iera];
    string* vera = mera[era];

    const char *ce = era.c_str();
    TH1D *h = new TH1D(Form("hf%s",ce),";|#eta|;HF scale",nx,vx);
    TH1D *hm = new TH1D(Form("hf%sm",ce),";|#eta|;HF scale",nx,vx);
    TH1D *hp = new TH1D(Form("hf%sp",ce),";|#eta|;HF scale",nx,vx);
    TH1D *he = new TH1D(Form("hf%se",ce),";|#eta|;HF scale",nx,vx);
    mh[era] = h;
    mhm[era] = hm;
    mhp[era] = hp;
    mhe[era] = he;
    
    cout << "Processing " << era << endl << flush;
    for (int i = 0; i != neta; ++i) {
      int ieta;
      float corr;
      if (sscanf(vera[i].c_str(),"%d: %f",&ieta,&corr)!=2) {
	cout << "Could not parse '"<<vera[i]<<"', exit" << endl << flush;
	exit(0);
      }
      int ix = (abs(ieta)-30) + 1;
      if (debug) cout << "ieta="<<ieta<<", corr="<<corr<<", ix="<<ix<<endl;
      //<< ", i="<<i<<"/"<<neta<<endl;
      double k(1);
      if (era=="2022") k = k22;
      if (era=="2022ee") k = k22ee;
      double scale = k/corr;
      if (ieta<0) hm->SetBinContent(ix, scale);
      if (ieta>0) hp->SetBinContent(ix, scale);
      h->SetBinContent(ix, h->GetBinContent(ix)+0.5*scale); // average
    } // for i

    TF1 *f1 = new TF1("f1","[0]*(1-[1]*log(cosh(x)))",3.139,4.013);
    TF1 *f2 = new TF1("f2","[0]*(1-[1]*log(cosh(x)))",2.964,5.191);

    h->Fit(f1,"QRNW");
    mf1[era] = f1;

    f2->SetParameters(f1->GetParameter(0), f1->GetParameter(1));
    mf2[era] = f2;

    // Do extrapolated HF scales
    for (int i = 1; i != h->GetNbinsX()+1; ++i) {
      double scale(1);
      double x = h->GetBinCenter(i);
      double xmin = h->GetBinLowEdge(i);
      if (i<=6) scale = h->GetBinContent(i);
      if (i>6) scale = f1->Eval(x);
      he->SetBinContent(i, scale);
      cout << Form("|ieta| %d: %1.4f",29+i,1./scale) << endl;

      hf->SetBinContent(i, hf->GetBinContent(i)+0.5/nera*scale);
    } // for i
  } // for iera


  TH1D *h = tdrHist("h","HF scale",0.3,1.4,"|#eta|",2.964,5.191);
  lumi_136TeV = "Run3, 2022+2023";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(2.964,1,5.191,1);

  tdrDraw(mh["2022"],"Pz",kFullSquare,kBlue);
  tdrDraw(mh["2022ee"],"Pz",kOpenSquare,kGreen+2);
  tdrDraw(mh["2023"],"Pz",kOpenCircle,kRed);
  tdrDraw(mh["2023bpix"],"Pz",kFullCircle,kRed);

  mf1["2022"]->SetLineColor(kBlue);
  mf1["2022"]->Draw("SAME");
  mf2["2022"]->SetLineColor(kBlue);
  mf2["2022"]->SetLineStyle(kDashed);
  mf2["2022"]->Draw("SAME");

  mf1["2022ee"]->SetLineColor(kGreen+2);
  mf1["2022ee"]->Draw("SAME");
  mf2["2022ee"]->SetLineColor(kGreen+2);
  mf2["2022ee"]->SetLineStyle(kDashed);
  mf2["2022ee"]->Draw("SAME");

  mf1["2023"]->SetLineColor(kRed-9);
  mf1["2023"]->Draw("SAME");
  mf2["2023"]->SetLineColor(kRed-9);
  mf2["2023"]->SetLineStyle(kDashed);
  mf2["2023"]->Draw("SAME");

  mf1["2023bpix"]->SetLineColor(kRed);
  mf1["2023bpix"]->Draw("SAME");
  mf2["2023bpix"]->SetLineColor(kRed);
  mf2["2023bpix"]->SetLineStyle(kDashed);
  mf2["2023bpix"]->Draw("SAME");
  
  TLegend *leg = tdrLeg(0.37,0.90-5*0.045,0.57,0.90);
  leg->AddEntry(mh["2023bpix"],"2023BPix","PLE");
  leg->AddEntry(mh["2023"],"2023","PLE");
  leg->AddEntry(mf1["2023bpix"],"2023: p_{0}(1-p_{1}log(cosh(x)))","L");
  leg->AddEntry(mf1["2022ee"],"22EE: p_{0}(1-p_{1}log(cosh(x)))","L");
  leg->AddEntry(mf1["2022"],"2022: p_{0}(1-p_{1}log(cosh(x)))","L");

  TLegend *leg2 = tdrLeg(0.61,0.90-2*0.045,0.81,0.90);
  leg2->AddEntry(mh["2022ee"],Form("2022EE  #times %1.2f",k22ee),"PLE");
  leg2->AddEntry(mh["2022"],Form("2022  #times %1.2f",k22),"PLE");

  /*
  // Do combined HF scales
  cout << "2022+2023" << endl;
  for (int i = 1; i != hf->GetNbinsX()+1; ++i) {
    double scale(1);
    double x = hf->GetBinCenter(i);
    double xmin = hf->GetBinLowEdge(i);
    //if (i==1) scale = 1;
    //if (i==2) scale = f1->Eval(x);
    //if (i>=3 && i<=6) scale = hf3->GetBinContent(i);
    if (i<=6) scale = hf3->GetBinContent(i);
    if (i>6) scale = f1->Eval(x);
    hf->SetBinContent(i, scale);

    cout << Form("|ieta| %d: %1.4f",29+i,1./scale) << endl;
  } // for i
  cout << "2023e:" << endl;
  for (int i = 1; i != hf23e->GetNbinsX()+1; ++i) {
    double scale(1);
    double x = hf23e->GetBinCenter(i);
    double xmin = hf23e->GetBinLowEdge(i);
    if (i<=6) scale = hf23->GetBinContent(i);
    if (i>6) scale = f23->Eval(x);
    hf23e->SetBinContent(i, scale);
    cout << Form("|ieta| %d: %1.4f",29+i,1./scale) << endl;
  } // for i
  cout << "2022e:" << endl;
  for (int i = 1; i != hf22e->GetNbinsX()+1; ++i) {
    double scale(1);
    double x = hf22e->GetBinCenter(i);
    double xmin = hf22e->GetBinLowEdge(i);
    if (i<=6) scale = hf22->GetBinContent(i);
    if (i>6) scale = f22->Eval(x);
    hf22e->SetBinContent(i, scale);
    cout << Form("|ieta| %d: %1.4f",29+i,1./scale) << endl;
  } // for i
  */
  //tdrDraw(hf,"Pz",kFullStar,kBlack);
  tdrDraw(mhe["2022"],"Pz",kFullStar,kBlue+2);
  tdrDraw(mhe["2022ee"],"Pz",kFullStar,kGreen+3);
  tdrDraw(mhe["2023bpix"],"Pz",kFullStar,kRed+2);
  
  TLegend *leg3 = tdrLeg(0.18,0.17,0.43,0.17+1*0.045);
  hf->SetMarkerStyle(kFullStar);
  //leg3->AddEntry(hf,"Combined estimates","PLE");
  leg3->AddEntry(hf,"Extrapolated estimates","PLE");

  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawHFscaleV2.pdf");

  TFile *fout = new TFile("rootfiles/drawHFscaleV2.root","RECREATE");
  mh["2022"]->Write("hf2022");
  mhe["2022"]->Write("hf2022e");
  mh["2022ee"]->Write("hf2022EE");
  mhe["2022ee"]->Write("hf2022EEe");
  mh["2023"]->Write("hf2023NonBPIX");
  mhe["2023"]->Write("hf2023NonBPIXe");
  mh["2023bpix"]->Write("hf2023");
  mhe["2023bpix"]->Write("hf2023e");
  fout->Close();
} // drawHFscale
