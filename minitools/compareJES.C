// Purpose: compare JES from two different MC (or data)
// In this first case, Summer23 vs Summer22
#include "../tdrstyle_mod22.C"

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

  return var;
} // jesFit


void compareJES() {

  gROOT->ProcessLine(".! mkdir pdf/compareJES");
  
  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //TFile *f18 = new TFile(" ../gamjet/files/GamHistosFill_mc_2018P8_v20.root","READ");
  //assert(f18 && !f18->IsZombie());
  //TProfile *pjes18 = (TProfile*)f18->Get("resp_JES_MC_a100_eta00_13");
  // doesn't seem to be stored for Run2 yet 
  
  TFile *f23 = new TFile("rootfiles/Summer23_L2ResOnly/gamjet_all/GamHistosFill_mc_2023P8-BPix_w4.root","READ");
  assert(f23 && !f23->IsZombie());
  TProfile *pjes23 = (TProfile*)f23->Get("resp_JES_MC_a100_eta00_13");

  TFile *f22 = new TFile("../gamjet/rootfiles/GamHistosFill_mc_2022P8_v32.root","READ");
  assert(f22 && !f22->IsZombie());
  TProfile *pjes22 = (TProfile*)f22->Get("resp_JES_MC_a100_eta00_13");
  
  TH1D *h = tdrHist("h","JES",0.885,1.00,"p_{T} (GeV)",15,4500);
  TH1D *hd = tdrHist("hd","Ratio",1.0,1.06,"p_{T} (GeV)",15,4500);
  lumi_136TeV = "Summer23 vs Summer22";
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);
  c1->cd(1);
  gPad->SetLogx();

  tdrDraw(pjes22,"HISTE][",kNone,kBlack,kSolid,-1,kNone);
  tdrDraw(pjes23,"HISTE][",kNone,kRed,kSolid,-1,kNone);

  pjes22->SetLineWidth(2);
  pjes23->SetLineWidth(2);
  
  TLegend *leg = tdrLeg(0.60,0.90-2*0.05,0.85,0.90);
  leg->AddEntry(pjes22,"Summer22","LE");
  leg->AddEntry(pjes23,"Summer23","LE");

  TF1 *fr22 = new TF1("fr22","[0]+[1]*log(x/1000.)+[2]*pow(x/2000.,1.5)",20,4500);
  TF1 *fr23 = new TF1("fr23","[0]+[1]*log(x/1000.)+[2]*pow(x/2000.,1.5)",20,4500);

  fr22->SetParameters(0.933,0.004,0.02);
  fr23->SetParameters(0.905,0.000,0.02);
  pjes22->Fit(fr22,"QRNW");
  pjes23->Fit(fr23,"QRNW");

  fr22->SetLineColor(kBlack);
  fr22->Draw("SAME");
  fr23->SetLineColor(kRed);
  fr23->Draw("SAME");

  TLegend *legb = tdrLeg(0.17,0.70-4*0.035,0.42,0.70);
  legb->SetTextSize(0.035);
  legb->AddEntry(fr22,fr22->GetExpFormula(),"L");
  legb->AddEntry(fr22,Form("%1.3f %+12.4f %+14.4f",
			   fr22->GetParameter(0),fr22->GetParameter(1),
			   fr22->GetParameter(2)),"");
  legb->AddEntry(fr23,fr22->GetExpFormula(),"L");
  legb->AddEntry(fr23,Form("%1.3f %+12.4f %+14.4f",
			   fr23->GetParameter(0),fr23->GetParameter(1),
			   fr23->GetParameter(2)),"");
  
  c1->cd(2);
  gPad->SetLogx();

  TH1D *hr = pjes22->ProjectionX();
  hr->Divide(pjes23);
  tdrDraw(hr,"HISTE][",kNone,kBlack,kSolid,-1,kNone);
  hr->SetLineWidth(2);

  TF1 *frr = new TF1("frr","([0]+[1]*log(x/1000.)+[2]*pow(x/2000.,1.5))/"
		     "([3]+[4]*log(x/1000.)+[5]*pow(x/2000.,1.5))",20,4500);
  frr->SetParameters(fr22->GetParameter(0),fr22->GetParameter(1),
		     fr22->GetParameter(2),
		     fr23->GetParameter(0),fr23->GetParameter(1),
		     fr23->GetParameter(2));
  frr->SetLineColor(kBlack);
  frr->Draw("SAME");
  
  // Edit shapes to be loaded in globalFitSettings.h
  loadShapes();
  TF1 *f1 = new TF1("f1",jesFit,15,4500,_vshape.size());
  TF1 *f2 = new TF1("f2",jesFit,15,4500,_vshape.size());
  TF1 *f3 = new TF1("f3",jesFit,15,4500,_vshape.size());
  for (int i = 0; i != f1->GetNpar(); ++i) {
    f1->SetParameter(i, 0);
    f1->SetParLimits(i, -1, 1);
    f2->SetParameter(i, 0);
    f2->SetParLimits(i, -2, 2);
    f3->SetParameter(i, 0);
    f3->SetParLimits(i, -3, 3);
  }
  hr->Fit(f3,"QRN");
  f3->SetLineColor(kRed);
  f3->Draw("SAME");
  hr->Fit(f2,"QRN");
  f2->SetLineColor(kOrange+2);
  f2->Draw("SAME");
  hr->Fit(f1,"QRN");
  f1->SetLineColor(kGreen+2);
  f1->Draw("SAME");
  
  TLegend *legd = tdrLeg(0.17,0.88-3*0.035*2,0.42,0.88);
  legd->SetTextSize(2*0.035);
  legd->AddEntry(f1,"ParLimits #pm1#sigma (full JES shape)","L");
  legd->AddEntry(f2,"ParLimits #pm2#sigma","L");
  legd->AddEntry(f3,"ParLimits #pm3#sigma","L");

  cout << "Fitted:" << endl;
  for (int i = 0; i != f1->GetNpar(); ++i) {
    cout << Form("%8s: %6.2f +/- %4.2f | %6.2f +/- %4.2f | %6.2f +/- %4.2f\n",
		 _vshape[i].name.c_str(),
		 f1->GetParameter(i), f1->GetParError(i),
    		 f2->GetParameter(i), f2->GetParError(i),
		 f3->GetParameter(i), f3->GetParError(i));
  }
  cout << endl;
  
  c1->SaveAs("pdf/compareJES/compareJES_Summer23_vs_Summer22.pdf");
}

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
