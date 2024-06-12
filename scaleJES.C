// Purpose: scale JES in rootfiles/jecdataX.root
// Input:   get reco-to-reco ratio from dijet/compareLite.C,
//          processed into an averaged graph
// Output:  directly modify mpfchs1_Y and ptchs_Y in jecdataX.root
#include "TFile.h"
#include "TGraphErrors.h"

#include <iostream>

void scaleJES(string run, string channel) {

  cout << "scalesJES("<<run<<", "<<channel<<")..." << endl << flush;
  
  const char *crun = run.c_str();  
  TDirectory *curdir = gDirectory;

  // Load reference reco-to-reco
  TFile *fr(0);
  if (run=="Run24CS")
    fr = new TFile("rootfiles/compareLite_2024CR_HCALDI_vs_ECALRATIO.root",
		   "READ");
  else
    fr = new TFile("rootfiles/compareLite_22Sep2023_vs_19Dec2023.root",
		   "READ");
  assert(fr && !fr->IsZombie());

  TString srun = crun;
  if (run!="Run3") {
    srun.ReplaceAll("Run","20"); // e.g. Run22CD->2022CD
    srun.ReplaceAll("2023C","2023Cv"); // e.g. 2023C123->2023Cv123
    //srun.ReplaceAll("2022FG","2022F"); // temporary fix
    srun.ReplaceAll("2024CS","2024C");
  }
  const char *refname = Form("djes_%s",srun.Data());
  cout << refname << endl << flush;
  TGraphErrors *gref = (TGraphErrors*)fr->Get(refname);
  assert(gref);
  gref->SetBit(TGraph::kIsSortedX); // speed up Eval()
  
  // Load graphs to be modified
  TFile *f = new TFile(Form("rootfiles/jecdata%s.root",crun),"UPDATE");
  assert(f && !f->IsZombie());

  string vdir[] = {"data","ratio"};
  const int ndir = sizeof(vdir)/sizeof(vdir[0]);
  string vset[] = {"mpfchs1","ptchs"};
  const int nset = sizeof(vset)/sizeof(vset[0]);

  for (int idir = 0; idir != ndir; ++idir) {
    for (int iset = 0; iset != nset; ++iset) {

      const char *cdir = Form("%s/eta00-13",vdir[idir].c_str());
      
      f->cd(cdir);
      TDirectory *d = gDirectory;
      cout << cdir << endl << flush;

      const char *cs = vset[iset].c_str();
      const char *cc = channel.c_str();
      const char *gname = Form("%s_%s_a100",cs,cc);
      cout << gname << endl << flush;
      TGraphErrors *g = (TGraphErrors*)d->Get(gname);
      assert(g);

      if (channel=="multijet") {

	TGraphErrors *gc = (TGraphErrors*)d->Get("crecoil_multijet_a100");
	assert(gc);
	
	for (int i = 0; i != g->GetN(); ++i) {

	  double pt = g->GetX()[i];
	  double r = gref->Eval(pt); // linear interpolation of <B>/<A>
	  double c = (r ? 1./r : 1.); // <A>/<B> for <19Dec>/<22Sep>
	  
	  //double ptr = gc->GetX()[i] * pt;
	  double ptr = gc->Eval(pt) * pt;
	  double rr = gref->Eval(ptr);
	  double cr = (rr ? 1./rr : 1.);
	  
	  g->SetPoint(i, pt, c/cr*g->GetY()[i]);
	  // To-do: for MPF, should only scale ptchs component in data
	  //        gets a bit tricky to do this for data/MC ratio
	} // for i
	
      } // multijet
      else {
	
	for (int i = 0; i != g->GetN(); ++i) {
	  double pt = g->GetX()[i];
	  double r = gref->Eval(pt); // linear interpolation of <B>/<A>
	  double c = (r ? 1/r : 1); // <A>/<B> for <19Dec>/<22Sep>
	  if (run=="Run24CS") {
	    r = gref->Eval(max(pt,40.)); // low pT unstable
	    c = (r ? r : 1.); // <A>/<B> for <HCALDI>/<X>
	  }
	  g->SetPoint(i, pt, c*g->GetY()[i]); // <B>/Z vs Z to <A>/Z vs Z
	  // To-do: for MPF, should only scale ptchs component in data
	  //        gets a bit tricky to do this for data/MC ratio
	} // for i
	
      } // !multijet

      g->Write(g->GetName(),TObject::kOverwrite);
    } // for iset
  } // for idir

  f->Write();
  f->Close();
} // scaleJES
