// Purpose: draw PF variations for HF scale
//          from files containing them in 2D (pT-eta)
#include "TFile.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <set>

#include "../tdrstyle_mod22.C"

void subtractAndScale(TH1D  *h, double off=1, double scale=1) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    if (!(h->GetBinError(i)==0 && h->GetBinContent(i)==0)) {
      h->SetBinContent(i, scale*(h->GetBinContent(i)-off));
    }
  }
} // subtractOne

// Implement basic fitShapes to check if can fit ratio
#include "../globalFitSettings.h"

// Generic fit function for JES, adapted from globalFit.C
vector<fitShape> _vshape; // obs->shapes
void loadShapes();
Double_t jesFit(Double_t *x, Double_t *p) {

  // Choose JES (var~1) or PF composition (var~0) as baseline
  double var = 1;
  double pt = x[0];

  // Load available shapes for this observable
  //vector<fitShape> &v = _mshape[_obs];
  if (_vshape.size()==0) loadShapes();
  vector<fitShape> &v = _vshape; // Adapted: Rjet only
  
  for (unsigned int i = 0; i != v.size(); ++i) {

    // Calculate variation and add it to total
    int idx = v[i].idx;   assert(idx>=0);
    TF1 *f1 = v[i].func;  assert(f1);
    double par = p[idx];
    if (v[i].ispos) par = max(par,0.);

    var += par * f1->Eval(pt) * 0.01; // fullSimShapes in %'s
  } // for i in mshape

  //return var;
  return (100.*(var-1));
} // jesFit

void draw2DvariationsHF() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_2Dvariations-HF.root","READ");
  //TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_2Dvariations-HF2.root","READ");
  //TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_2Dvariations-HF2-noRespCut.root","READ");
  //TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_2Dvariations-HF3-noRespCut.root","READ");
  TFile *f = new TFile("rootfiles/JME-Run3Summer22_1M_2Dvariations-HF4-noRespCut.root","READ");
  assert(f && !f->IsZombie());

  string year = "2022";//"2023"
  const char *cy = year.c_str();
  TFile *f2(0);
  //if (year=="2022") f2 = new TFile("rootfiles/L2Res_Summer22.root","READ");
  //if (year=="2022") f2 = new TFile("rootfiles/L2Res_Summer22_v1b.root","READ");
  //if (year=="2022") f2 = new TFile("rootfiles/L2Res_Summer22_v2_ECMixup_flatHF.root","READ"); // Bad for EC but ok for HF
  if (year=="2022") f2 = new TFile("rootfiles/L2Res_Summer22_v2_ECMixup_slopeHF.root","READ"); // Bad for EC but ok for HF
  //if (year=="2023") f2 = new TFile("rootfiles/L2Res_Summer23.root","READ");
  if (year=="2023") f2 = new TFile("rootfiles/L2Res_Summer23_v2.root","READ");
  assert(f2 && !f2->IsZombie());

  curdir->cd();

  TH2D *hr(0), *hr_1(0), *hr_2(0), *hr_3(0), *hr_4(0);
  if (year=="2023") {
    hr_1 = (TH2D*)f2->Get("h2jes_2023Cv123"); assert(hr_1);
    hr_2 = (TH2D*)f2->Get("h2jes_2023Cv4"); assert(hr_2);
    hr_3 = (TH2D*)f2->Get("h2jes_2023D"); assert(hr_3);
    hr = (TH2D*)hr_3->Clone("hr_2023");
    hr->Add(hr_1,hr_2);
    hr->Add(hr_3);
    hr->Scale(1./3.);
  }
  if (year=="2022") {
    hr_1 = (TH2D*)f2->Get("h2jes_2022CD"); assert(hr_1);
    hr_2 = (TH2D*)f2->Get("h2jes_2022E"); assert(hr_2);
    hr_3 = (TH2D*)f2->Get("h2jes_2022F"); assert(hr_3);
    hr_4 = (TH2D*)f2->Get("h2jes_2022G"); assert(hr_4);
    hr = (TH2D*)hr_4->Clone("hr_2022"); hr->Reset();
    hr->Add(hr_1,hr_2);
    hr->Add(hr_3);
    hr->Add(hr_4);
    hr->Scale(1./4.);
    //hr = (TH2D*)hr_1->Clone("hr_2022");
  }
  
  TProfile2D *pr(0);
  //pr = (TProfile2D*)f->Get("RjetPuppi_HFCalib23"); assert(pr);
  //pr = (TProfile2D*)f->Get("RjetPuppi_HFCalibComb"); assert(pr);
  if (year=="2023") pr = (TProfile2D*)f->Get("RjetPuppi_HFCalib23e");
  //if (year=="2022") pr = (TProfile2D*)f->Get("RjetPuppi_HFCalib22e");//_v1
  if (year=="2022") pr = (TProfile2D*)f->Get("RjetPuppi_HFCalib22e_v2");
  assert(pr);

  int color[] =
    {kMagenta+2, kMagenta+1, kBlue, kCyan+1, kCyan+2, kGreen+2,
     kYellow+1, kYellow+2, kOrange+1, kOrange+2, kRed, kBlack, kGray+1};
  const int ncolor = sizeof(color)/sizeof(color[0]);
  
  //TH1D *h = tdrHist("h","HF response variation (%)",-48,18);//15);//-23,7);
  TH1D *h = tdrHist("h","HF JES change (%)",-48,18);
  lumi_136TeV = "Run3, 2023";
  extraText = "Private";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(15,0,3500,0);

  TLegend *leg = tdrLeg(0.63,0.15,0.88,0.15);
  leg->SetTextSize(0.030);

  int j1 = pr->GetYaxis()->FindBin(2.964+0.05);
  int j2 = pr->GetYaxis()->FindBin(5.191-0.05);
  for (int j = j1; j != j2+1; ++j) {
    TH1D *h1 = pr->ProjectionX(Form("h1_%d",j),j,j);
    subtractAndScale(h1,1,100);
    tdrDraw(h1,"HIST][",kNone,color[(j-j1)%ncolor],kSolid,-1,kNone);
    leg->AddEntry(h1,Form("%1.3f#leq|#eta|<%1.3f",
			  pr->GetYaxis()->GetBinLowEdge(j),
			  pr->GetYaxis()->GetBinLowEdge(j+1)), "PLE");
    leg->SetY2(leg->GetY2()+0.030);
  } // for j
  int k1 = hr->GetXaxis()->FindBin(2.964+0.05);//3.664+0.05);//3.489+0.05);
  int k2 = hr->GetXaxis()->FindBin(5.191-0.05);//3.839-0.05);//4.013-0.05);
  for (int k = k1; k != k2+1; ++k) {
    TH1D *hr1 = hr->ProjectionY(Form("hr1_%d",k),k,k);//k1,k2);
    //hr1->Scale(1./(k2-k1+1));
    subtractAndScale(hr1,1,100);
    int j = pr->GetYaxis()->FindBin(hr->GetXaxis()->GetBinLowEdge(k)+0.05);
    tdrDraw(hr1,"Pz",kFullCircle,color[(j-j1)%ncolor],kDotted);//kBlack);
    hr1->SetMarkerSize(0.7);
    hr1->GetXaxis()->SetRangeUser(15.,3500.);
    if (k==k2) {
      leg->AddEntry(hr,Form("Data [%1.3f,%1.3f]",
			    hr->GetXaxis()->GetBinLowEdge(k1),
			    hr->GetXaxis()->GetBinLowEdge(k2+1)),"PLE");
      leg->SetY2(leg->GetY2()+0.030);
    }
  } // for k
  
  gPad->RedrawAxis();
  gPad->Update();
  
  c1->SaveAs(Form("pdf/draw2DvariationsHF/draw2DvariationsHF_vsPt_%s.pdf",cy));

  
  //TH1D  *h_2 = tdrHist("h_2","HF response variation (%)",-23,7,
  //		       "|#eta|",2.964,5.191);
  //TH1D *h_2 = new TH1D("h_2",";|#eta|;HF response variation (%);",
  TH1D *h_2 = new TH1D("h_2",";|#eta|;HF JES change (%);",
		       10*12,2.964,5.191);
  h_2->GetYaxis()->SetRangeUser(-48,29.999);//-40,25);//-23,7);
  TCanvas *c2 = tdrCanvas("c2",h_2,8,11,kSquare);

  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+1);
  l->DrawLine(2.964,0,5.191,0);

  TLegend *leg2 = tdrLeg(0.50,0.88,0.75,0.88);
  leg2->SetTextSize(0.030);

  int i1 = pr->GetXaxis()->FindBin(21+0.1);//15.+0.1);
  int i2 = pr->GetXaxis()->FindBin(468.-0.1);
  for (int i = i1; i != i2+1; ++i) {
    TH1D *h2 = pr->ProjectionY(Form("h2_%d",i),i,i);
    subtractAndScale(h2,1,100);
    tdrDraw(h2,"HIST][",kNone,color[(i-i1)%ncolor],kSolid,-1,kNone);
    leg2->AddEntry(h2,Form("%1.0f#leqp_{T}<%1.0f GeV",
			   pr->GetXaxis()->GetBinLowEdge(i),
			   pr->GetXaxis()->GetBinLowEdge(i+1)), "PLE");
    leg2->SetY1(leg2->GetY1()-0.030);
  } // for j
  int m1 = hr->GetYaxis()->FindBin(21+0.1);
  int m2 = hr->GetYaxis()->FindBin(468.-0.1);
  for (int m = m1; m != m2+1; ++m) {
    TH1D *hr2 = hr->ProjectionX(Form("hr2_%d",m),m,m);
    subtractAndScale(hr2,1,100);
    int i = pr->GetXaxis()->FindBin(hr->GetYaxis()->GetBinLowEdge(m)+0.05);
    tdrDraw(hr2,"Pz",kFullCircle,color[(i-i1)%ncolor],kDotted);
    hr2->SetMarkerSize(0.7);
    hr2->GetXaxis()->SetRangeUser(15.,3500.);
    if (m==m2) {
      leg2->AddEntry(hr2,Form("Data [%1.0f,%1.0f]",
			      hr->GetYaxis()->GetBinLowEdge(m1),
			      hr->GetYaxis()->GetBinLowEdge(m2+1)),"PLE");
      leg2->SetY2(leg2->GetY2()+0.030);
    }
  } // for k
  
  gPad->RedrawAxis();
  gPad->Update();
  
  c2->SaveAs(Form("pdf/draw2DvariationsHF/draw2DvariationsHF_vsEta_%s.pdf",cy));
  c2->SaveAs(Form("pdf/draw2DvariationsHF/draw2DvariationsHF_vsEta_%s.png",cy));

} // draw2Dvariations



// Adapted from globalFit.C
void loadShapes() {

  // Set whitelist for quickly selecting only subset of shapes
  set<string> whitelistshape;
  for (unsigned int i = 0; i != _gf_shapes_whitelist.size(); ++i) {
    if (_gf_shapes_whitelist[i]!="")
      whitelistshape.insert(_gf_shapes_whitelist[i]);
  }

  // Set list of positive-definite sources
  set<string> posdeflist;
  for (unsigned int i = 0; i != _gf_posdef.size(); ++i) {
      posdeflist.insert(_gf_posdef[i]);
  }
 
  // Create listing and mapping of active response/composition shapes
  set<string> shapes;
  map<string, int> mshapeidx;
  for (unsigned int i = 0; i != _gf_shapes.size(); ++i) {
    
    string name       = _gf_shapes[i][0];
    string appliesTo  = _gf_shapes[i][1];
    string funcstring = _gf_shapes[i][2];
    if (appliesTo!="Rjet") continue; // Adapted: Rjet only
    
    if (name=="") continue; // ignore extra empty (commented out) elements
    // Check if input dataset is whitelisted
    if (!whitelistshape.empty() &&
	whitelistshape.find(name)==whitelistshape.end()) continue;

    // Use running indexing for active shapes
    if (shapes.find(name)==shapes.end())
      mshapeidx[name] = shapes.size();

    //if (debug) cout << "..." << name << "->" << appliesTo << endl << flush;
    cout << (i==0? "Load shapes: " : ", ") << name;
    //<< "(" << funcstring << ")" << endl << flush;

    // Create fitShape object
    fitShape shape;
    shape.idx = mshapeidx[name];
    shape.name = name;
    shape.appliesTo = appliesTo;
    shape.ispos = (posdeflist.find(name)!=posdeflist.end());
    shape.func = new TF1(Form("f1_%s_%s",name.c_str(),appliesTo.c_str()),
			 funcstring.c_str(),15.,6500.);

    // Save fitShapes to list and map
    shapes.insert(name);
    //_mshape[appliesTo].push_back(shape);
    _vshape.push_back(shape); // Adapted: Rjet only
  } // for i in shapes
  cout << endl;
} // loadShapes
