// Purpose: draw comparison of JER SF for |eta|<0.783 for all eras
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include "../tdrstyle_mod22.C"
#include "../tools.C"

TH1D *profileY(const char *name, TH2D *h2, int ix1, int ix2, bool w1=false) {
  TH1D *h1 = h2->ProjectionY(name,ix1,ix2); h1->Scale(1./3.);
  // Do better merging
  for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
    double sumw(0), sum(0), esum2(0);
    for (int i = ix1; i != ix2; ++i) {
      double y = h2->GetBinContent(i,j);
      double ey = h2->GetBinError(i,j);
      double w = (ey!=0 ? 1./pow(ey,2) : 0);
      if (w1) w = 1;
      if (w1) ey = max(ey,0.025);
      sumw += w;
      sum += w*y;
      esum2 += pow(w*ey,2);
    } // for i
    double ynew = (sumw!=0 ? sum/sumw : 0);
    double eynew = (sumw!=0 ? sqrt(esum2/sumw)/sqrt(sumw) : 0);
    h1->SetBinContent(j, ynew);
    h1->SetBinError(j, eynew);
  } //for j

  return h1;
} // profileY

void drawJERSFcomparisonV2() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  //TFile *f = new TFile("rootfiles/JERSF_Run3_vs_Summer23BPix.root","READ");
  TFile *f = new TFile("rootfiles/JERSF_Run3_vs_Summer23BPix_upto_2024I.root","READ");
  assert(f && !f->IsZombie());

  string vera[] = {
    "2022C_Prompt2022","2022D_Prompt2022",
    "2022C_22Sep2023","2022D_22Sep2023","2022E_22Sep2023",
       "2022E_22Sep2023",
    "2022F_22Sep2023","2022G_22Sep2023",
    "2023Cv123_Prompt2023","2023Cv4_Prompt2023","2023D_Prompt2023",
       "2023D_Prompt2023",
    "2022C_19Dec2023","2022D_19Dec2023","2022E_19Dec2023",
      "2022E_19Dec2023",
    "2022F_19Dec2023","2022G_19Dec2023",
    "2023Cv123_19Dec2023","2023Cv4_19Dec2023","2023D_19Dec2023",
      "2023D_19Dec2023",
      //"2024BC"
      "2024BCD","2024E","2024F","2024G","2024H","2024I"
  };
  const int nera = sizeof(vera)/sizeof(vera[0]);

  map<string,int> color;
  color["2022C_Prompt2022"] = kGreen+1;
  color["2022D_Prompt2022"] = kGreen+2;
  color["2022C_22Sep2023"] = kGreen+1;
  color["2022D_22Sep2023"] = kGreen+2;
  color["2022E_22Sep2023"] = kCyan+2;
  color["2022F_22Sep2023"] = kRed;
  color["2022G_22Sep2023"] = kRed+1;
  color["2023Cv123_Prompt2023"] = kBlue;
  color["2023Cv4_Prompt2023"] = kMagenta+1;
  color["2023D_Prompt2023"] = kMagenta+2;
  color["2022C_19Dec2023"] = kGreen+1;
  color["2022D_19Dec2023"] = kGreen+2;
  color["2022E_19Dec2023"] = kCyan+2;
  color["2022F_19Dec2023"] = kRed;
  color["2022G_19Dec2023"] = kRed+2;
  color["2023Cv123_19Dec2023"] = kBlue;
  color["2023Cv4_19Dec2023"] = kMagenta+1;
  color["2023D_19Dec2023"] = kMagenta+2;
  color["2024BC"] = kBlack;
  color["2024BCD"] = kYellow+1;
  color["2024E"] = kYellow+2;
  color["2024F"] = kOrange;
  color["2024G"] = kOrange+1;
  color["2024H"] = kYellow+2;
  color["2024I"] = kYellow+3;

  map<string,int> marker;
  marker["2022C_Prompt2022"] = kOpenSquare;
  marker["2022D_Prompt2022"] = kOpenSquare;
  marker["2022C_22Sep2023"] = kFullSquare;
  marker["2022D_22Sep2023"] = kFullSquare;
  marker["2022E_22Sep2023"] = kFullSquare;
  marker["2022F_22Sep2023"] = kFullSquare;
  marker["2022G_22Sep2023"] = kFullSquare;
  marker["2023Cv123_Prompt2023"] = kOpenSquare;
  marker["2023Cv4_Prompt2023"] = kOpenSquare;
  marker["2023D_Prompt2023"] = kOpenSquare;
  marker["2022C_19Dec2023"] = kFullCircle;
  marker["2022D_19Dec2023"] = kFullCircle;
  marker["2022E_19Dec2023"] = kFullCircle;
  marker["2022F_19Dec2023"] = kFullCircle;
  marker["2022G_19Dec2023"] = kFullCircle;
  marker["2023Cv123_19Dec2023"] = kFullCircle;
  marker["2023Cv4_19Dec2023"] = kFullCircle;
  marker["2023D_19Dec2023"] = kFullCircle;
  marker["2024BC"] = kFullCircle;//kOpenSquare;
  marker["2024BCD"] = kOpenStar;
  marker["2024E"] = kOpenStar;
  marker["2024F"] = kFullStar;
  marker["2024G"] = kFullStar;
  marker["2024H"] = kFullStar;
  marker["2024I"] = kFullStar;

  map<string,const char*> label;
  label["2022C_Prompt2022"] = "22C_Pr22";
  label["2022D_Prompt2022"] = "22D_Pr22";
  label["2022C_22Sep2023"] = "22C_22Sep";
  label["2022D_22Sep2023"] = "22D_22Sep";
  label["2022E_22Sep2023"] = "22E_22Sep";
  label["2022F_22Sep2023"] = "22F_22Sep";
  label["2022G_22Sep2023"] = "22G_22Sep";
  label["2023Cv123_Prompt2023"] = "23C3_Pr23";
  label["2023Cv4_Prompt2023"] = "23C4_Pr23";
  label["2023D_Prompt2023"] = "23D_Pr23";
  label["2022C_19Dec2023"] = "22C_19Dec";
  label["2022D_19Dec2023"] = "22D_19Dec";
  label["2022E_19Dec2023"] = "22E_19Dec";
  label["2022F_19Dec2023"] = "22F_19Dec";
  label["2022G_19Dec2023"] = "22G_19Dec";
  label["2023Cv123_19Dec2023"] = "23C3_19Dec";
  label["2023Cv4_19Dec2023"] = "23C4_19Dec";
  label["2023D_19Dec2023"] = "23D_19Dec";
  label["2024BC"] = "24BC_Pr24";
  label["2024BCD"] = "24BCD_Pr24";
  label["2024E"] = "24E_Pr24";
  label["2024F"] = "24F_Pr24";
  label["2024G"] = "24G_Pr24";
  label["2024H"] = "24H_Pr24";
  label["2024I"] = "24I_Pr24";

  map<string,int> bin;
  bin["2022C_Prompt2022"] = 1;
  bin["2022D_Prompt2022"] = 2;
  bin["2022C_22Sep2023"] = 4;
  bin["2022D_22Sep2023"] = 5;
  bin["2022E_22Sep2023"] = 6;
  bin["2022F_22Sep2023"] = 10;
  bin["2022G_22Sep2023"] = 11;
  bin["2023Cv123_Prompt2023"] = 14;
  bin["2023Cv4_Prompt2023"] = 15;
  bin["2023D_Prompt2023"] = 16;
  bin["2022C_19Dec2023"] = 7;
  bin["2022D_19Dec2023"] = 8;
  bin["2022E_19Dec2023"] = 9;
  bin["2022F_19Dec2023"] = 12;
  bin["2022G_19Dec2023"] = 13;
  bin["2023Cv123_19Dec2023"] = 17;
  bin["2023Cv4_19Dec2023"] = 18;
  bin["2023D_19Dec2023"] = 19;
  bin["2024BC"] = 20;
  bin["2024BCD"] = 20;
  bin["2024E"] = 21;
  bin["2024F"] = 22;
  bin["2024G"] = 23;
  bin["2024H"] = 24;
  bin["2024I"] = 25;
  
  TH1D *h = tdrHist("h","JER SF",0.5,3.0);
  lumi_136TeV = "Run3 vs Summer23BPix MC";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,1,3500,1);

  const double etamax = 0.783;
  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.35,0.87,Form("|#eta| < %1.3f",etamax));
  
  double lsize = 0.025;
  TLegend *leg = tdrLeg(0.35,0.85-nera*0.7*lsize,0.35+2*0.20,0.85);
  leg->SetNColumns(2);
  leg->SetTextSize(lsize);

  //int nb = bin["2024BC"];//bin.size();
  int nb = bin["2024I"];//bin.size();
  TH1D *h100 = new TH1D("h100",";Era;JER SF",nb,0.5,nb+0.5);
  TH1D *h250 = new TH1D("h250",";Era;JER SF",nb,0.5,nb+0.5);
  TH1D *h1000 = new TH1D("h1000",";Era;JER SF",nb,0.5,nb+0.5);
  TH1D *h2000 = new TH1D("h2000",";Era;JER SF",nb,0.5,nb+0.5);

  vector<TH1D*> vh1(nera);
  vector<TH1D*> vh1r(nera);
  // Load results
  for (int i = 0; i != nera; ++i) {

    string run = vera[i];
    const char *cr = run.c_str();
    const char *cm = "Summer23MGBPix";
  
    TH2D *h2(0), *h2r(0);
    h2 = (TH2D*)f->Get(Form("Fits/h2jersf_%s_%s",cr,cm));
    h2r = (TH2D*)f->Get(Form("Dijet/h2jersfRaw_%s_%s",cr,cm));
    assert(h2);
    assert(h2r);

    int ix1 = h2->GetXaxis()->FindBin(0.);
    int ix2 = h2->GetXaxis()->FindBin(etamax-0.05);
    TH1D *h1 = profileY(Form("h1_%s_%d",cr,i),h2,ix1,ix2,true);
    TH1D *h1r = profileY(Form("h1r_%s_%d",cr,i),h2r,ix1,ix2);

    vh1[i] = h1;
    vh1r[i] = h1r;
  } // for i (load)

  // Draw fit bands at the back
  for (int i = 0; i != nera; ++i) {

    string run = vera[i];
    const char *cr = run.c_str();
    TH1D *h1 = vh1[i];
    
    tdrDraw(h1,"E3",kNone,color[cr],kSolid,-1,1001,color[cr]-9);
    tdrDraw(new TGraph(h1),"L",kNone,color[cr],kSolid,-1);
    h1->SetFillColorAlpha(color[cr]-9,0.70);
  } // for i (back)
  
  // Draw raw data in front
  for (int i = 0; i != nera; ++i) {

    string run = vera[i];
    const char *cr = run.c_str();
    TH1D *h1 = vh1[i];
    TH1D *h1r = vh1r[i];
    
    tdrDraw(h1r,"Pz",marker[cr],color[cr],kSolid,-1,1001,color[cr]-9);
    h1r->SetMarkerSize(0.5);

    //leg->AddEntry(h1r,label[cr],"PLEF");
    leg->AddEntry(h1r,label[cr],"PLE");

    int ib = bin[cr];
    int ix100 = h1r->FindBin(100.);
    int ix250 = h1r->FindBin(250.);
    int ix1000 = h1r->FindBin(1000.);
    int ix2000 = h1r->FindBin(2000.);
    if (h100->GetBinContent(ib)==0) {
      h100->SetBinContent(ib, h1r->GetBinContent(ix100));
      h100->SetBinError(ib, h1r->GetBinError(ix100));
      h100->GetXaxis()->SetBinLabel(ib, label[cr]);

      h250->SetBinContent(ib, h1r->GetBinContent(ix250));
      h250->SetBinError(ib, h1r->GetBinError(ix250));
      h1000->SetBinContent(ib, h1r->GetBinContent(ix1000));
      h1000->SetBinError(ib, h1r->GetBinError(ix1000));
      h2000->SetBinContent(ib, h1->GetBinContent(ix2000));
      h2000->SetBinError(ib, h1->GetBinError(ix2000));

      cout << "cr " << cr << ": h1="<<h1->GetBinContent(ix2000)
	   << " (bin="<<ix2000<<")" << endl;
    }
    else {
      if (!(h100->GetXaxis()->GetBinLabel(ib)==label[cr]))
	cout << h100->GetXaxis()->GetBinLabel(ib) << " vs " << label[cr]
	     << endl << flush;
      assert(string(h100->GetXaxis()->GetBinLabel(ib))==string(label[cr]));
    }
      
  } // for i (front)
  
  gPad->RedrawAxis();

  c1->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparisonV2.pdf");
  //c1->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparisonV2.png");

  h100->GetYaxis()->SetRangeUser(0.9,1.5);//0.8,2.0);//0.90,1.5);  
  h100->GetXaxis()->SetLabelSize(0.025);
  TCanvas *c2 = tdrCanvas("c2",h100,8,11,kSquare);
  tdrDraw(h250,"Pz",kFullSquare,kGreen+2);
  tdrDraw(h100,"Pz",kOpenCircle,kBlue);
  tdrDraw(h1000,"Pz",kFullCircle,kRed);
  //tdrDraw(h1000,"Pz",kFullCircle,kOrange+2);
  //tdrDraw(h2000,"Pz",kFullDiamond,kRed); h2000->SetMarkerSize(1.5);

  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+2);
  l->DrawLine(0.5,1,nb+0.5,1);//20.5,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(0.5,1.1,nb+0.5,1.1);//20.5,1.1);

  double y1 = 0.9;//0.8
  double y2 = 1.35;//2.0
  double y3 = 1.44;//1.5;//2.0
  l->SetLineStyle(kSolid);
  l->SetLineColor(kGray+1);
  l->DrawLine(bin["2022E_19Dec2023"]+0.5,y1,bin["2022E_19Dec2023"]+0.5,y2);
  l->DrawLine(bin["2022G_19Dec2023"]+0.5,y1,bin["2022G_19Dec2023"]+0.5,y2);
  l->DrawLine(bin["2023D_19Dec2023"]+0.5,y1,bin["2023D_19Dec2023"]+0.5,y3);
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray);
  l->DrawLine(bin["2022C_22Sep2023"]-0.5,y1,bin["2022C_22Sep2023"]-0.5,y2);
  l->DrawLine(bin["2022C_19Dec2023"]-0.5,y1,bin["2022C_19Dec2023"]-0.5,y2);
  l->DrawLine(bin["2022F_19Dec2023"]-0.5,y1,bin["2022F_19Dec2023"]-0.5,y2);
  l->DrawLine(bin["2023Cv123_19Dec2023"]-0.5,y1,bin["2023Cv123_19Dec2023"]-0.5,y2);
  l->DrawLine(bin["2024F"]-0.5,y1,bin["2024F"]-0.5,y2);

  tex->SetTextSize(0.045);
  tex->DrawLatex(0.70,0.87,Form("|#eta| < %1.3f",etamax));

  TLegend *leg2 = tdrLeg(0.35,0.90-0.04*3,0.65,0.90);
  leg2->AddEntry(h1000,"1000 GeV","PLE");
  leg2->AddEntry(h100,"100 GeV","PLE");
  leg2->AddEntry(h250,"250 GeV","PLE");

  tex->SetNDC(kFALSE);
  double d = 0.35;
  tex->DrawLatex(1,1.36,"2022(EE)");
  tex->DrawLatex(10-d,1.36,"22EE");//"2022EE");
  tex->DrawLatex(14-d,1.36,"23(BPix)");//"2023(BPix)");
  tex->DrawLatex(14-d,1.36,"23(BPix)");//"2023(BPix)");
  tex->DrawLatex(20-d,1.36,"2024");//"2023(BPix)");
  
  tex->SetTextSize(0.030);
  //tex->DrawLatex(1,1.36,"2022(EE)");
  //tex->DrawLatex(1,1.18,"Prompt");
  //tex->DrawLatex(1,1.16,"C");
  //tex->DrawLatex(2,1.16,"D");
  tex->DrawLatex(1,1.34,"Prompt");
  tex->DrawLatex(1,1.32,"C");
  tex->DrawLatex(2,1.32,"D");

  tex->DrawLatex(4-d,1.18,"22Sep");
  tex->DrawLatex(4-d,1.14,"C");
  tex->DrawLatex(5-d,1.15,"D");
  tex->DrawLatex(6-d,1.16,"E");

  tex->DrawLatex(7-d,1.19,"19Dec");
  tex->DrawLatex(7-d,1.15,"C");
  tex->DrawLatex(8-d,1.16,"D");
  tex->DrawLatex(9-d,1.17,"E");

  //tex->DrawLatex(10,1.36,"2022EE");
  tex->DrawLatex(10-d,1.18,"22Sep");
  tex->DrawLatex(10-d,1.14,"F");
  tex->DrawLatex(11-d,1.15,"G");
  
  tex->DrawLatex(12-d,1.20,"19Dec");
  tex->DrawLatex(12-d,1.16,"F");
  tex->DrawLatex(13-d,1.17,"G");

  //tex->DrawLatex(14,1.36,"2023(BPix)");
  tex->DrawLatex(14-d,1.24,"22Sep");
  tex->DrawLatex(14-d,1.17,"C3");
  tex->DrawLatex(15-d,1.21,"C4");
  tex->DrawLatex(16-d,1.22,"D");
  
  tex->DrawLatex(17-d,1.25,"19Dec");
  tex->DrawLatex(17-d,1.20,"C3");
  tex->DrawLatex(18-d,1.22,"C4");
  tex->DrawLatex(19-d,1.23,"D");

  gPad->RedrawAxis();
  c2->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparisonV2_eras.pdf");
  //c2->SaveAs("pdf/drawJERSFcomparison/drawJERSFcomparisonV2_eras.png");
  
} // void drawJERSFcomparisonV2
