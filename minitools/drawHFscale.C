// Purpose: draw HF scale from James Natoli
//          extrapolate to higher ieta
#include "TF1.h"

#include "../tdrstyle_mod22.C"

void drawHFscale() {

  setTDRStyle();
  
  // From James Vincent Natoli, 21 March 2024
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

  // From James Vincent Natoli, 21 March 2024
  string v2023[] =
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
  const int n2023 = sizeof(v2023)/sizeof(v2023[0]);
  

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
  TH1D *hf23e = new TH1D("hf23e",";|#eta|;HF scale",nx,vx);
  TH1D *hf22e = new TH1D("hf22e",";|#eta|;HF scale",nx,vx);
  
  TH1D *hf22m = new TH1D("hf22m",";|#eta|;HF scale",nx,vx);
  TH1D *hf22p = new TH1D("hf22p",";|#eta|;HF scale",nx,vx);
  TH1D *hf22 = new TH1D("hf22",";|#eta|;HF scale",nx,vx);

  double k22 = 1.18;//1;//1.17;
  for (int i = 0; i != n2022; ++i) {
    int ieta;
    float corr;
    if (sscanf(v2022[i].c_str(),"%d: %f",&ieta,&corr)!=2) exit(0);
    int ix = (abs(ieta)-30) + 1;
    cout << "ieta="<<ieta<<", corr="<<corr<<", ix="<<ix<<endl;
    double scale = k22/corr;
    if (ieta<0) hf22m->SetBinContent(ix, scale);
    if (ieta>0) hf22p->SetBinContent(ix, scale);
    hf22->SetBinContent(ix, hf22->GetBinContent(ix)+0.5*scale); // average
    hf3->SetBinContent(ix, hf3->GetBinContent(ix)+0.25*scale); // average
  }

  TH1D *hf23m = new TH1D("hf23m",";|#eta|;HF scale",nx,vx);
  TH1D *hf23p = new TH1D("hf23p",";|#eta|;HF scale",nx,vx);
  TH1D *hf23 = new TH1D("hf23",";|#eta|;HF scale",nx,vx);

  for (int i = 0; i != n2023; ++i) {
    int ieta;
    float corr;
    if (sscanf(v2023[i].c_str(),"%d: %f",&ieta,&corr)!=2) exit(0);
    int ix = (abs(ieta)-30) + 1;
    cout << "ieta="<<ieta<<", corr="<<corr<<", ix="<<ix<<endl;
    double scale = 1./corr;
    if (ieta<0) hf23m->SetBinContent(ix, scale);
    if (ieta>0) hf23p->SetBinContent(ix, scale);
    hf23->SetBinContent(ix, hf23->GetBinContent(ix)+0.5*scale); // average
    hf3->SetBinContent(ix, hf3->GetBinContent(ix)+0.25*scale); // average
  }

  TH1D *h = tdrHist("h","HF scale",0.3,1.4,"|#eta|",2.964,5.191);
  lumi_136TeV = "Run3, 2022+2023";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(2.964,1,5.191,1);

  /*
  tdrDraw(hf22m,"Pz",kOpenSquare,kBlue);
  tdrDraw(hf22p,"Pz",kOpenSquare,kRed);
  tdrDraw(hf22,"Pz",kFullSquare,kBlue);//kBlack);
  hf22->SetMarkerSize(0.8);

  tdrDraw(hf23m,"Pz",kOpenCircle,kBlue);
  tdrDraw(hf23p,"Pz",kOpenCircle,kRed);
  tdrDraw(hf23,"Pz",kFullCircle,kRed);//kBlack);
  hf23->SetMarkerSize(0.8);
  */

  tdrDraw(hf22m,"Pz",kOpenTriangleDown,kBlue);
  tdrDraw(hf22p,"Pz",kOpenTriangleUp,kBlue);
  tdrDraw(hf22,"Pz",kFullSquare,kBlue);
  hf22->SetMarkerSize(0.8);

  tdrDraw(hf23m,"Pz",kOpenTriangleDown,kRed);
  tdrDraw(hf23p,"Pz",kOpenTriangleUp,kRed);
  tdrDraw(hf23,"Pz",kFullCircle,kRed);
  hf23->SetMarkerSize(0.8);

  
  //tdrDraw(hf3,"Pz",kFullDiamond,kGreen+2);

  //TF1 *f22 = new TF1("f22","[0]*(1-[1]*cosh(x))",3.314,4.013);
  //TF1 *f22 = new TF1("f22","[0]*(1-[1]*log(cosh(x)))",3.314,4.013);
  TF1 *f22 = new TF1("f22","[0]*(1-[1]*log(cosh(x)))",3.139,4.013);
  //TF1 *f22 = new TF1("f22","[0]*(1-[1]*cosh(x)-[2]*log(cosh(x)))",3.314,4.013);
  f22->SetLineColor(kBlue);
  hf22->Fit(f22,"QRNW");
  f22->DrawClone("SAME");
  f22->SetRange(2.964,5.191);
  f22->SetLineStyle(kDashed);
  f22->DrawClone("SAME");
  f22->SetLineStyle(kSolid);

  cout << Form("2022: p0=%6.3f+/-%5.3f, p1=%6.3f+/-%5.3f",
	       f22->GetParameter(0),f22->GetParError(0),
  	       f22->GetParameter(1),f22->GetParError(1)) << endl;
	       
  
  //TF1 *f23 = new TF1("f23","[0]*(1-[1]*cosh(x))",2.964,4.013);
  //TF1 *f23 = new TF1("f23","[0]*(1-[1]*log(cosh(x)))",2.964,4.013);
  TF1 *f23 = new TF1("f23","[0]*(1-[1]*log(cosh(x)))",3.139,4.013);
  //TF1 *f23 = new TF1("f23","[0]*(1-[1]*cosh(x)-[2]*log(cosh(x)))",2.964,4.013);
  f23->SetLineColor(kRed);
  hf23->Fit(f23,"QRNW");
  f23->DrawClone("SAME");
  f23->SetRange(2.964,5.191);
  f23->SetLineStyle(kDashed);
  f23->DrawClone("SAME");
  f23->SetLineStyle(kSolid);

  cout << Form("2023: p0=%6.3f+/-%5.3f, p1=%6.3f+/-%5.3f",
	       f23->GetParameter(0),f23->GetParError(0),
	       f23->GetParameter(1),f23->GetParError(1)) << endl;
  TF1 *f1 = new TF1("f1","[0]*(1-[1]*log(cosh(x)))",3.139,4.013);

  f1->SetLineColor(kGreen+2);
  hf3->Fit(f1,"QRNW");
  //f1->DrawClone("SAME");
  f1->SetRange(2.964,5.191);
  f1->SetLineStyle(kDashed);
  //f1->DrawClone("SAME");
  f1->SetLineStyle(kSolid);

  cout << Form("Run3: p0=%6.3f+/-%5.3f, p1=%6.3f+/-%5.3f",
	       f1->GetParameter(0),f1->GetParError(0),
  	       f1->GetParameter(1),f1->GetParError(1)) << endl;

  TLegend *leg = tdrLeg(0.40,0.90-5*0.045,0.65,0.90);
  leg->AddEntry(hf23,"2023","PLE");
  leg->AddEntry(hf23m,"2023-","P");
  leg->AddEntry(hf23p,"2023+","P");
  //leg->AddEntry(hf3,Form("2023 + 2022#times%1.2f",k22),"P");
  //leg->AddEntry(f1,"p_{0}(1-p_{1}log(cosh(x)))","L");
  leg->AddEntry(f23,"2023: p_{0}(1-p_{1}log(cosh(x)))","L");
  leg->AddEntry(f22,"2022: p_{0}(1-p_{1}log(cosh(x)))","L");

  TLegend *leg2 = tdrLeg(0.60,0.90-3*0.045,0.85,0.90);
  if (k22!=1) {
    leg2->AddEntry(hf22,Form("2022  #times %1.2f",k22),"PLE");
    leg2->AddEntry(hf22m,Form("2022- #times %1.2f",k22),"P");
    leg2->AddEntry(hf22p,Form("2022+ #times %1.2f",k22),"P");
  }
  else {
    leg2->AddEntry(hf22,"2022","PLE");
    leg2->AddEntry(hf22m,"2022-","P");
    leg2->AddEntry(hf22p,"2022+","P");
  }
    
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

  //tdrDraw(hf,"Pz",kFullStar,kBlack);
  tdrDraw(hf22e,"Pz",kFullStar,kBlue+2);
  tdrDraw(hf23e,"Pz",kFullStar,kRed+2);
  
  TLegend *leg3 = tdrLeg(0.18,0.17,0.43,0.17+1*0.045);
  hf->SetMarkerStyle(kFullStar);
  //leg3->AddEntry(hf,"Combined estimates","PLE");
  leg3->AddEntry(hf,"Extrapolated estimates","PLE");
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawHFscale.pdf");

  TFile *fout = new TFile("rootfiles/drawHFscale.root","RECREATE");
  hf22->Write("hf2022");
  hf22e->Write("hf2022e");
  hf23->Write("hf2023");
  hf23e->Write("hf2023e");
  fout->Close();
} // drawHFscale
