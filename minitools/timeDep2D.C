// L2L3Res time dependence with inclusive jets
// Use Incjet/h2pteta_lumi
// Projection sanity check with Incjet/p2pteta vs hpt13
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"

#include <fstream>
#include <cstdio>
#include <string>

#include "../tdrstyle_mod22.C"

bool debug2 = true;//false;

// 2D cross section parameterization (new)
Double_t fx(Double_t *x, Double_t *p) {

  double eta = x[0];
  double pt = x[1];
  double N = p[0];
  double alpha = p[1];
  double beta0 = p[2];
  double beta2 = p[3];
  double beta5 = p[4];
  double pt0 = p[5];
  double etaref = 2.5;
  double etaref2 = 5.0;

  double y = fabs(eta);
  double z = etaref;
  double z2 = etaref2;
  //double beta  = beta0 + (beta2-beta0) * min(y,z)/z; // linear in |eta|
  //double beta  = beta0 + (beta2-beta0) * y/z; // linear in |eta|
  //double beta  = beta0 + (beta2-beta0) * sqrt(y/z); // sqrt in |eta|
  double beta  = beta0 + (beta2-beta0) * min(y,z)/z
    + (beta5-beta2) * (max(min(y,z2),z)-z)/(z2-z); // lin x3 in |eta|
  double erg = pt*cosh(eta);
  double sqrts = 13600.;
  if (erg>=sqrts*0.5) return 0;
  double f = p[0]*exp(-pt0/pt)*pow(pt,alpha)*pow(1-2.*erg/sqrts,beta);

  return f;
}

// 2D cross section parameterization (old)
Double_t fx_old(Double_t *x, Double_t *p) {

  double eta = x[0];
  double pt = x[1];
  double N = p[0];
  double alpha = p[1];
  double beta0 = p[2];
  double beta2 = p[3];
  double etaref = 2.5;

  double y = fabs(eta);
  double z = etaref;
  double beta  = beta0 + (beta2-beta0) * min(y,z)/z; // linear in |eta|
  double erg = pt*cosh(eta);
  double sqrts = 13600.;
  if (erg>=sqrts*0.5) return 0;
  double f = p[0]*pow(pt,alpha)*pow(1-2.*erg/sqrts,beta);

  return f;
}

// Calculate xsec ratio for given JES and [etamin,etamax,ptmin,ptmax]
TF2 *_f2(0);
Double_t fjes(Double_t *x, Double_t *p) {

  double jes = (*x);
  double etamin = p[0];
  double etamax = p[1];
  double ptmin = min(p[2],0.5*13600./cosh(etamin));
  double ptmax = min(p[3],0.5*13600./cosh(etamin));

  assert(_f2);
  double xsec = _f2->Integral(etamin,etamax,ptmin,ptmax);
  double xsec2 = _f2->Integral(etamin,etamax,ptmin/jes,ptmax/jes);

  return (xsec>0 ? xsec2/xsec : 0);
}

/*
// Simple struct to hold the region boundaries for each trigger
struct TriggerBox {
  std::string triggerName;
  double ptMin;
  double ptMax;
  double absyMin;
  double absyMax;
};
// Utility to check if a bin center is in range
inline bool inRange(double val, double minVal, double maxVal) {
  return (val >= minVal) && (val < maxVal);
}
*/
// Include DijetHistosCombine.C to define TriggerBox, incjetBoxes, inRange
#include "../minitools/DijetHistosCombine.C"
string triggerString(double pt, double eta) {
  
  // Quick printout from files:
  // for (int i = 1; i != hpt13->GetNbinsX()+1; ++i) cout << hpt13->GetBinLowEdge(i) << (hpt13->GetBinContent(i)>hpt13->GetBinContent(i-1) ? "*":"") << ", "; cout << endl;
  
  // Inclusive jet thresholds, mi[]
  /*
  // Forward triggers
  if (pt>=686. &&           fabs(eta)>2.964) return "HLT_PFJetFwd500";
  if (pt>=548. && pt<686. && fabs(eta)>2.964) return "HLT_PFJetFwd450";
  if (pt>=468. && pt<548. && fabs(eta)>2.964) return "HLT_PFJetFwd400";
  if (pt>=395. && pt<468. && fabs(eta)>2.964) return "HLT_PFJetFwd320";
  if (pt>=330. && pt<395. && fabs(eta)>2.964) return "HLT_PFJetFwd260";
  if (pt>=272. && pt<330. && fabs(eta)>2.964) return "HLT_PFJetFwd200";
  if (pt>=196. && pt<272. && fabs(eta)>3.139) return "HLT_PFJetFwd140";
  if (pt>=114. && pt<196. && fabs(eta)>3.139) return "HLT_PFJetFwd80";
  if (pt>=84.  && pt<114. && fabs(eta)>3.239) return "HLT_PFJetFwd60";
  if (pt>=49.  && pt<84.  && fabs(eta)>2.965) return "HLT_PFJetFwd40";

  // Normal inclusive triggers
  if (pt>=686.) return "HLT_PFJet500";
  if (pt>=548.) return "HLT_PFJet450";
  if (pt>=468.) return "HLT_PFJet400";
  if (pt>=395.) return "HLT_PFJet320";
  if (pt>=330.) return "HLT_PFJet260";
  if (pt>=272.) return "HLT_PFJet200";
  if (pt>=196.) return "HLT_PFJet140";
  //if (pt>196.) return "HLT_PFJet110";
  if (pt>=114.) return "HLT_PFJet80";
  if (pt>=84.)  return "HLT_PFJet60";
  if (pt>=49.)  return "HLT_PFJet40";
  return "HLT_ZeroBias";
  */
  /*
  // Inclusive jet thresholds, incjetBoxes from DijetHistosCombine.C
  // Inclusive jets pT bins in relevan trigger range
  //28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245,
  //272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846,
  // Example: define a vector of TriggerBoxes for the inclusive-jet folder
  std::vector<TriggerBox> incjetBoxes = {
    //{"HLT_ZeroBias", 0,   49,  0, 5.2},
    //{"HLT_ZeroBias", 0,   56,  0.0, 0.5},
    //{"HLT_ZeroBias", 0,   56,  0.5, 1.0},
    {"HLT_ZeroBias", 0,   64,  0.0, 0.5},
    {"HLT_ZeroBias", 0,   64,  0.5, 1.0},
    {"HLT_ZeroBias", 0,   64,  1.0, 1.5},
    {"HLT_ZeroBias", 0,   64,  1.5, 2.0},
    {"HLT_ZeroBias", 0,  114,  2.0, 2.5},
    {"HLT_ZeroBias", 0,   97,  2.5, 3.0},
    //{"HLT_PFJet40",  49,  84,  0, 3.0},
    //{"HLT_PFJet40", 56,  84,  0.0, 0.5},
    //{"HLT_PFJet40", 56,  84,  0.5, 1.0},
    {"HLT_PFJet40", 64,  84,  0.0, 0.5},
    {"HLT_PFJet40", 64,  84,  0.5, 1.0},
    {"HLT_PFJet40", 64,  84,  1.0, 1.5},
    {"HLT_PFJet40", 64,  84,  1.5, 2.0},
    {"HLT_PFJet40", 114, 153,  2.0, 2.5},
    //{"HLT_PFJet40", 114, 133,  2.0, 2.5}, // 25_V2M
    {"HLT_PFJet40", 97,  114,  2.5, 3.0},  
    //{"HLT_PFJet60",  84,  114, 0, 3.0},
    {"HLT_PFJet60", 84,  114,  0.0, 0.5},
    {"HLT_PFJet60", 84,  114,  0.5, 1.0},
    {"HLT_PFJet60", 84,  114,  1.0, 1.5},
    {"HLT_PFJet60", 84,  114,  1.5, 2.0},
    {"HLT_PFJet60", 153, 174,  2.0, 2.5},
    //{"HLT_PFJet60", 133, 153,  2.0, 2.5}, // 25_V2M
    {"HLT_PFJet60", 114, 133,  2.5, 3.0},  
    //{"HLT_PFJet80",  114, 196, 0, 3.0},
    
    // 2024_V9M
    {"HLT_PFJet80", 114, 196,  0.0, 0.5},
    {"HLT_PFJet80", 114, 196,  0.5, 1.0},
    {"HLT_PFJet80", 114, 196,  1.0, 1.5},
    {"HLT_PFJet80", 114, 196,  1.5, 2.0},
    {"HLT_PFJet80", 174, 196,  2.0, 2.5},
    {"HLT_PFJet80", 133, 196,  2.5, 3.0},  

    {"HLT_PFJet140", 196, 272, 0, 3.0},
    {"HLT_PFJet200", 272, 330, 0, 3.0},
    {"HLT_PFJet260", 330, 395, 0, 3.0},
    {"HLT_PFJet320", 395, 468, 0, 3.0},
    {"HLT_PFJet400", 468, 548, 0, 3.0},
    {"HLT_PFJet450", 548, 686, 0, 3.0},
    {"HLT_PFJet500", 686,7000, 0, 3.0},
    
    //{"HLT_ZeroBias", 0,   49,  0, 5.2}, 
    {"HLT_ZeroBias", 0,   64,  3.0, 3.5},
    {"HLT_ZeroBias", 0,   64,  3.5, 4.0},
    {"HLT_ZeroBias", 0,   64,  4.0, 4.5},
    {"HLT_ZeroBias", 0,   64,  4.5, 5.2},
    //{"HLT_PFJetFwd40",  49,  84,  3.0, 5.2},
    {"HLT_PFJetFwd40",  64,   97,  3.0, 3.5},
    {"HLT_PFJetFwd40",  64,   97,  3.5, 4.0},
    {"HLT_PFJetFwd40",  64,   97,  4.0, 4.5},
    {"HLT_PFJetFwd40",  64,   97,  4.5, 5.2},
    //{"HLT_PFJetFwd60",  84,  114, 3.0, 5.2},
    {"HLT_PFJetFwd60",  97,  220,  3.0, 3.5},
    {"HLT_PFJetFwd60",  97,  220,  3.5, 4.0},
    {"HLT_PFJetFwd60",  97,  220,  4.0, 4.5},
    {"HLT_PFJetFwd60",  97,  153,  4.5, 5.2},
    //{"HLT_PFJetFwd80",  114, 196, 3.0, 5.2},
    // Fwd80 somehow botched up, very bad turn-on at ~196 (last 133?)
    {"HLT_PFJetFwd80", 220,  220,  3.0, 3.5},
    {"HLT_PFJetFwd80", 220,  220,  3.5, 4.0},
    {"HLT_PFJetFwd80", 220,  220,  4.0, 4.5},
    {"HLT_PFJetFwd80", 153,  153,  4.5, 5.2},
    //{"HLT_PFJetFwd140", 196, 272, 3.0, 5.2},
    {"HLT_PFJetFwd140", 220,  272,  3.0, 3.5},
    {"HLT_PFJetFwd140", 220,  272,  3.5, 4.0},
    {"HLT_PFJetFwd140", 220,  272,  4.0, 4.5},
    {"HLT_PFJetFwd140", 153, 7000,  4.5, 5.2}, // E<150 for hpt50
    
    //{"HLT_PFJetFwd200", 272, 330, 3.0, 5.2},
    {"HLT_PFJetFwd200", 272, 330, 3.0, 4.0},
    {"HLT_PFJetFwd200", 272, 300, 4.0, 4.5},
    {"HLT_PFJetFwd200",7000,7000, 4.5, 5.2}, // 150-200 for hpt50
    //{"HLT_PFJetFwd260", 330, 395, 3.0, 5.2},
    {"HLT_PFJetFwd260", 330, 362, 3.0, 4.0},
    {"HLT_PFJetFwd260", 300,7000, 4.0, 4.5}, // 300 for hpt45
    {"HLT_PFJetFwd260",7000,7000, 4.5, 5.2}, // 150-200 for hpt50
    //{"HLT_PFJetFwd320", 395, 468, 3.0, 5.2},
    {"HLT_PFJetFwd320", 362, 430, 3.0, 4.0},
    {"HLT_PFJetFwd320",7000,7000, 4.0, 4.5},
    {"HLT_PFJetFwd320",7000,7000, 4.5, 5.2},
    //{"HLT_PFJetFwd400", 468, 7000, 3.0, 5.2},
    {"HLT_PFJetFwd400", 430,7000, 3.0, 4.0},
    {"HLT_PFJetFwd400",7000,7000, 4.0, 4.5},
    {"HLT_PFJetFwd400",7000,7000, 4.5, 5.2},
    //{"HLT_PFJetFwd400", 468, 548, 3.0, 5.2},
    //{"HLT_PFJetFwd450", 548, 686, 3.0, 5.2},
    //{"HLT_PFJetFwd500", 686,7000, 3.0, 5.2},

  };
*/
  string trigger("none");
  bool covered = false;
  for (auto& box : incjetBoxes) {
    // Check if (pt, eta) is in the region
    if (inRange(pt, box.ptMin, box.ptMax) &&
	inRange(fabs(eta), box.absyMin, box.absyMax)) {
      
      // Copy bin content for this trigger
      trigger = box.triggerName;
      
      // Sanity check for overlapping boxes
      if (covered) {
	std::cerr << "Warning: trigger " << trigger
		  << " bin (pt=" << pt << ", eta=" << eta << ")"
		  << " already covered by another trigger!\n";
      }
      
      // Mark as covered and break if you expect only one coverage
      covered = true;
      //break;
    }
  }

  // Sanity check for missing phase space corners
  if (!covered) {
    std::cerr << "Warning: bin (pt=" << pt
	      << ", eta=" << eta << ") not covered by any trigger!\n";
  }
  
  return trigger;
} // triggerString

map<string, map<string, double> > _trgLumCache;
double triggerLumi(string run, string trg) {
  if (_trgLumCache[run][trg]!=0) return _trgLumCache[run][trg];

  const char *cr = run.c_str();
  //string sf = Form("rootfiles/Prompt2024/lumiText/lumi_HLT_Golden385863_%s_PFJet.csv",cr);
  const char *ct = (trg=="HLT_ZeroBias" ? "ZB" : "PFJet");
  string sf = Form("rootfiles/NestorLumiJSON/lumi_HLT_Golden_%s_%s.csv",cr,ct);
  if (debug2) cout << "Opening " << sf << endl << flush;
  ifstream f(sf.c_str());
  assert(f.is_open());

  double lumi(0);
  float rec(0), del(0); // in /fb
  int ver,nfill,nrun,ncms,nruntot(0);
  string s;
  string sRegExp = Form("#%s_v%%d,%%d,%%d,%%d,%%f,%%f",trg.c_str());
  const char *csr = sRegExp.c_str();
  if (debug2) cout << "Reading regExp " << csr << endl << flush;
  while (f >> s) {
    const char *cs = s.c_str();
    if (sscanf(cs,csr,&ver,&nfill,&nrun,&ncms,&del,&rec)==6) {
      if (debug2) cout << "Parsing string " << cs << endl << flush;
      if (debug2) cout << trg << "_v" << ver<<": " << rec << endl << flush;
      lumi += rec;
      nruntot += nrun;
    }
  } // while f

  if (debug2) cout << trg << ": " << lumi << " ("<<nruntot<<" runs)" << endl << flush;
  if (lumi==0) lumi = -1;

  // Ad-hoc factors to match eras (brilcalc by hand over .csv)
  //if (run=="2024BCD") lumi *= 15.3/14.7;
  //if (run=="2024E") lumi *= 11.3/11.3;
  //if (run=="2024F") lumi *= 27.8/27.8;
  //if (run=="2024G") lumi *= 37.8/33.1;// * 0.97;

  // Ad-hoc factors to match eras
  bool doAdhocLumi = true;
  if (doAdhocLumi) {
    if (run=="2024C_nib1") lumi *= 1.06;
    if (run=="2024E_nib1") lumi *= 0.94;
    if (run=="2024I_nib1") lumi *= 0.50*0.96;

    if (run=="2024F_nib1" && trg=="HLT_ZeroBias") lumi *= 0.613;
    if (run=="2024F_nib2" && trg=="HLT_ZeroBias") lumi *= 0.206;
    if (run=="2024F_nib3" && trg=="HLT_ZeroBias") lumi *= 0.352;
    if (run=="2024G_nib1" && trg=="HLT_ZeroBias") lumi *= 0.287;
    if (run=="2024G_nib2" && trg=="HLT_ZeroBias") lumi *= 0.386;
    if (run=="2024H_nib1" && trg=="HLT_ZeroBias") lumi *= 0.945;

    if (run=="2025D" && trg=="HLT_ZeroBias") lumi *= 0.90;
    if (run=="2025D" && trg=="HLT_PFJet40")  lumi *= 0.90;
    if (run=="2025D" && trg=="HLT_PFJet60")  lumi *= 0.90;
    if (run=="2025D" && trg=="HLT_PFJet80")  lumi *= 0.90;

    /*
    double k140 = 1.03;
    if (run=="2025C" && trg=="HLT_PFJet140") lumi *= k140*0.96;
    if (run=="2025D" && trg=="HLT_PFJet140") lumi *= k140*0.95;//0.96;
    if (run=="2025E" && trg=="HLT_PFJet140") lumi *= k140*0.96;
    if (run=="2025F" && trg=="HLT_PFJet140") lumi *= k140*0.97;
    if (run=="2025G" && trg=="HLT_PFJet140") lumi *= k140*0.99;//1.01;
    if (run=="2025CDEFG" && trg=="HLT_PFJet140") lumi *= k140*0.985;
    */
    double k140 = 1.00;//1.03
    if (run=="2025C" && trg=="HLT_PFJet140") lumi *= k140*0.99;//0.96;
    if (run=="2025D" && trg=="HLT_PFJet140") lumi *= k140*0.98;//0.95;
    if (run=="2025E" && trg=="HLT_PFJet140") lumi *= k140*0.99;//0.96;
    if (run=="2025F" && trg=="HLT_PFJet140") lumi *= k140*1.00;//0.97;
    if (run=="2025G" && trg=="HLT_PFJet140") lumi *= k140*1.02;//0.99;
    if (run=="2025CDEFG" && trg=="HLT_PFJet140") lumi *= k140*1.030;//1.025;//1.045;//0.985;
    
    //if (run=="2025G" && trg=="HLT_PFJet80") lumi *= 1.01;

    double w25(109.5), wc(21.8), wd(25.5), we(14.0), wf(26.4), wg(21.8);
    if (run=="2025CDEFG" &&
	(trg=="HLT_ZeroBias" || trg=="HLT_PFJet40" || trg=="HLT_PFJet60" ||
	 trg=="HLT_PFJet80"))
      lumi *= (0.90*wd + 1.00*(wc+we+wf+wg)) / w25;
    
    // v155+v156 issues
    if (run=="2025F" && trg!="HLT_ZeroBias") lumi *= 0.853;//0.86;//0.88;//0.82;//0.84;
    if (run=="2025CDEFG" && trg!="HLT_ZeroBias") lumi *= 0.965;//0.97;//0.95; 
    //if (run=="2025CDEFG" && trg=="HLT_PFJet140") lumi *= 1.03;
    
    //if (run=="2025CDEFG" && trg=="HLT_PFJet140")
    //lumi *= (0.96*wc + 0.96*wd + 0.96*we + 0.97*wf + 1.01*wg) / w25;
    //if (run=="2025CDEFG" && trg=="HLT_PFJet140")
    //lumi *= (0.96*wc + 0.95*wd + 0.96*we + 0.97*wf + 0.99*wg) / w25;
  }

  _trgLumCache[run][trg] = lumi * 1e3; // fb-1 to pb-1

  return lumi;
} //triggerLumi

void normalizeHistByLumiAndBin(TH1D *h, string run) {
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double eta = 0.;
    double pt = h->GetBinCenter(i);
    double dpt = h->GetXaxis()->GetBinWidth(i);
    string trg = triggerString(pt, eta);
    double lumi = triggerLumi(run, trg);
    h->SetBinContent(i, h->GetBinContent(i)/lumi/dpt);
    h->SetBinError(i, h->GetBinError(i)/lumi/dpt);
  }
} // normalizeHistByLumi (TH1D)

void normalizeHistByLumiAndBin(TH2D *h2, string run) {
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      double eta = h2->GetXaxis()->GetBinCenter(i);
      double deta = h2->GetXaxis()->GetBinWidth(i);
      double pt = h2->GetYaxis()->GetBinCenter(j);
      double dpt = h2->GetYaxis()->GetBinWidth(j);
      string trg = triggerString(pt, eta);
      double lumi = triggerLumi(run, trg);
      h2->SetBinContent(i, j, h2->GetBinContent(i,j)/lumi/deta/dpt);
      h2->SetBinError(i, j, h2->GetBinError(i,j)/lumi/deta/dpt);
    }
  }
} // normalizeHistByLumi (TH2D)

TH1D* timeDep1Ds(string run="2025CDEFG", double lum = 109) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  const char *cr = run.c_str();
  TString t(cr);
  //TFile *f = new TFile(Form("rootfiles/Prompt2024/v110_2024/jmenano_data_cmb_%s_JME_v110_2024.root",cr),"READ");
  TFile *f(0);
  //if (t.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_%s_JME_v153_copy.root",cr),"READ");
  if (t.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Jet_v156/jmenano_data_cmb_%s_JME_v156_copy.root",cr),"READ");
  //if (t.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_%s_JME_v155.root",cr),"READ");
  if (t.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_%s_JME_v155_copy.root",cr),"READ");
  assert(f && !f->IsZombie());

  TH1D *hpt13 = (TH1D*)f->Get("Incjet/hpt13"); assert(hpt13);
  TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta"); assert(h2pteta);
  //TH2D *h2pteta_lum = (TH2D*)f->Get("Incjet/h2pteta_lumi"); assert(h2pteta_lum);

  hpt13 = (TH1D*)hpt13->Clone(Form("hpt13_%s",cr));
  //hpt13->Scale(0.001/lum,"width");
  normalizeHistByLumiAndBin(hpt13, run);
  
  int i1 = h2pteta->GetXaxis()->FindBin(-1.305+1e-2);
  int i2 = h2pteta->GetXaxis()->FindBin(+1.305-1e-2);
  TH1D *hpt13p = h2pteta->ProjectionY(Form("hpt13p_%s",cr),i1,i2);
  //hpt13p->Scale(0.001/lum,"width");
  normalizeHistByLumiAndBin(hpt13p, run);

  TH1D *hpt13l = h2pteta->ProjectionY(Form("hpt13l_%s",cr),i1,i2);
  //TH1D *hpt13l = h2pteta_lum->ProjectionY(Form("hpt13l_%s",cr),i1,i2);
  //double lum500 = triggerLumi(run,"HLT_PFJet500");
  //double lum500 = 109;
  //hpt13l->Scale(22.5/lum500*1e-5,"width"); // ??
  //hpt13l->Scale(0.001*lum/lum500,"width");
  //hpt13l->Scale(1./lum/lum,"width");
  //hpt13l->Scale(0.001/lum,"width");
  normalizeHistByLumiAndBin(hpt13l, run);
  //if (run=="2025CDEFG") hpt13l->Scale(1./(2+1./7.));
  //if (run=="2025D")     hpt13l->Scale(1./7.);
  
  TH1D *h = tdrHist(Form("h1_%s",cr),"Cross section (pb / GeV)",0.2e-7,2e11,
		    "p_{T,jet} (GeV)",15,4500.);
  TH1D *hd = tdrHist(Form("hd_%s",cr),"Ratio",0.95,1.05,
		     "p_{T,jet} (GeV)",15,4500.);
  lumi_136TeV = Form("%s, %1.1f fb^{-1}",cr,lum);//lum500);
  extraText = "Private";
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cr),h,hd,8,11);

  //normalizeHistByLumi(hpt13,cr);
  //normalizeHistByLumi(hpt13p,cr);
  
  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();
  tdrDraw(hpt13,"Pz",kOpenCircle);
  tdrDraw(hpt13p,"Pz",kFullDiamond);
  tdrDraw(hpt13l,"Pz",kOpenDiamond);

  gPad->RedrawAxis();
  
  c1->cd(2);
  gPad->SetLogx();
  TH1D *hpt13r = (TH1D*)hpt13p->Clone(Form("hpt13r_%s",cr));
  hpt13r->Divide(hpt13);
  tdrDraw(hpt13r,"Pz",kFullDiamond);

  TH1D *hpt13lr = (TH1D*)hpt13l->Clone(Form("hpt13lr_%s",cr));
  hpt13lr->Divide(hpt13p);
  //hpt13lr->Scale(1.04);
  tdrDraw(hpt13lr,"Pz",kOpenDiamond);
  
  gPad->RedrawAxis();

  //c1->SaveAs(Form("pdf/timeDep2D/timeDep2D_1DBB_%s.pdf",cr));
  c1->Close();

  //return hpt13p;
  return hpt13l;
} // TH1D* timeDep1Ds

void timeDep2Ds(string run="2025G", string ref="2025CDEFG", TH1D *hr=0) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  // NB: if trigger binning unknown, use minitools/DijetHistosCombine.C to re-merge a copy
  
  const char *cr = run.c_str();
  TString t(cr);
  //TFile *f = new TFile(Form("rootfiles/Prompt2024/v110_2024/jmenano_data_cmb_%s_JME_v110_2024.root",cr),"READ");
  TFile *f(0);
  //if (t.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_%s_JME_v153_copy.root",cr),"READ");
  if (t.Contains("2025")) f = new TFile(Form("rootfiles/Prompt2025/Jet_v156/jmenano_data_cmb_%s_JME_v156_copy.root",cr),"READ");
  //if (t.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_%s_JME_v155.root",cr),"READ");
  if (t.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_%s_JME_v155_copy.root",cr),"READ");
  assert(f && !f->IsZombie());

  const char *cref = ref.c_str();
  //TFile *fref = new TFile(Form("rootfiles/Prompt2024/v110_2024/jmenano_data_cmb_%s_JME_v110_2024.root",cref),"READ");
  //TFile *fref = new TFile(Form("rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_%s_JME_v153_copy.root",cref),"READ");
  TFile *fref = new TFile(Form("rootfiles/Prompt2025/Jet_v156/jmenano_data_cmb_%s_JME_v156_copy.root",cref),"READ");
  assert(fref && !fref->IsZombie());

  TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta"); assert(h2pteta);
  TH2D *h2pteta_ref = (TH2D*)fref->Get("Incjet/h2pteta"); assert(h2pteta_ref);
  //TH2D *h2pteta = (TH2D*)f->Get("Incjet/h2pteta_lumi"); assert(h2pteta);
  //TH2D *h2pteta_ref = (TH2D*)fref->Get("Incjet/h2pteta_lumi"); assert(h2pteta_ref);
  curdir->cd();
  
  double lum500 = triggerLumi(run,"HLT_PFJet500") * 1e-3;
  TH1D *h = tdrHist(Form("h2_%s",cr),"p_{T} (GeV)",15,4500,
		    "#eta (rad)",-5.2,5.2);
  h->GetYaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetNoExponent();

  lumi_136TeV = Form("%s, %1.1f fb^{-1}",cr,lum500);
  extraText = "Private";
  TCanvas *c2 = tdrCanvas(Form("c2_%s",cr),h,8,11,kSquare);
  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.20);
  h->GetYaxis()->SetTitleOffset(1.5);
  
  TH2D *h2 = (TH2D*)h2pteta->Clone(Form("h2pteta_%s",cr));
  normalizeHistByLumiAndBin(h2,cr);

  TH2D *h2ref = (TH2D*)h2pteta_ref->Clone(Form("h2pteta_%s_%s",cr,cref));
  normalizeHistByLumiAndBin(h2ref,cref);

  gPad->SetLogy();
  gPad->SetLogz();
  h2->Draw("SAME COLZ");
  h2->GetZaxis()->SetRangeUser(0.2e-6,2.5e8);
  h2->GetZaxis()->SetTitle("Cross section (pb / GeV / rad)");
  h2->GetZaxis()->SetTitleOffset(1.6);
  h2->GetZaxis()->SetTitleSize(0.05);
  h2->GetZaxis()->SetLabelSize(0.05);

  gPad->RedrawAxis();

  // Draw lines of "safe space"
  TLine *l = new TLine();
  l->SetLineStyle(kDotted);
  l->SetLineColor(kGray+2);
  l->DrawLine(-2.5,30,+2.5,30);
  l->DrawLine(-2.5,30,-2.5,1109);
  l->DrawLine(+2.5,30,+2.5,1109);
  TF1 *f1 = new TF1(Form("f1_%s",cr),"13600.*0.5/cosh(x)",-5.191,+5.191);
  f1->SetLineStyle(kDotted);
  f1->SetLineColor(kGray+2);
  f1->Draw("SAME");

  gPad->Update();
  c2->SaveAs(Form("pdf/timeDep2D/timeDep2D_xsec_%s.pdf",cr));

  
  TH1D *h_3 = tdrHist(Form("h_3_%s",cr),"p_{T} (GeV)",15,4500,
		      "#eta (rad)",-5.2,5.2);
  h_3->GetYaxis()->SetMoreLogLabels();
  h_3->GetYaxis()->SetNoExponent();
  TCanvas *c3 = tdrCanvas(Form("c3_%s",cr),h_3,8,11,kSquare);
  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.20);
  h_3->GetYaxis()->SetTitleOffset(1.5);

  TH2D *h2r = (TH2D*)h2->Clone(Form("h2r_%s",cr));
  h2r->Divide(h2ref);

  gPad->SetLogy();
  h2r->Draw("SAME COLZ");
  double eps = 1e-4;
  h2r->GetZaxis()->SetRangeUser(0.75+eps,1.25-eps);
  h2r->GetZaxis()->SetTitle(Form("%s / %s",cr,cref));
  h2r->GetZaxis()->SetTitleOffset(1.6);
  h2r->GetZaxis()->SetTitleSize(0.05);
  h2r->GetZaxis()->SetLabelSize(0.05);

  gPad->RedrawAxis();

  // Draw lines of "safe space"
  l->DrawLine(-2.5,30,+2.5,30);
  l->DrawLine(-2.5,30,-2.5,1109);
  l->DrawLine(+2.5,30,+2.5,1109);
  f1->Draw("SAME");
  
  gPad->Update();
  c3->SaveAs(Form("pdf/timeDep2D/timeDep2D_ratio_%s.pdf",cr));

  c2->Close();
  c3->Close();
  

  TH1D *h_4 = tdrHist(Form("h_4_%s",cr),"p_{T} (GeV)",15,4500,
		      "#eta (rad)",-5.2,5.2);
  h_4->GetYaxis()->SetMoreLogLabels();
  h_4->GetYaxis()->SetNoExponent();
  TCanvas *c4 = tdrCanvas(Form("c4_%s",cr),h_4,8,11,kSquare);
  gPad->SetLeftMargin(0.18);
  gPad->SetRightMargin(0.20);
  h_4->GetYaxis()->SetTitleOffset(1.5);

  //TF2 *f2 = new TF2(Form("f2_%s",cr),fx_old,-5.2,5.2,15,4500,4);
  //f2->SetParameters(2e14,-5,4.,4.); // pT^alpha
  ///f2->FixParameter(1,-5);
  //f2->FixParameter(2,8.);
  //f2->FixParameter(3,4.);
  //TF2 *f2 = new TF2(Form("f2_%s",cr),fx,-5.2,5.2,15,4500,5);
  //f2->SetParameters(3e14,-5.00,9.5,4.75,38.); // pT^alpha, lin in |eta|
  TF2 *f2 = new TF2(Form("f2_%s",cr),fx,-5.2,5.2,15,4500,6);
  f2->SetParameters(3e14,-5.0,9.5,4.75,4.,38.); // pT^alpha, lin x2 in |eta|
  //h2ref->Fit(f2);
  TH2D *h2rf = (TH2D*)h2->Clone(Form("h2rf_%s",cr)); h2rf->Reset();
  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      double eta = h2->GetXaxis()->GetBinCenter(i);
      double pt = h2->GetYaxis()->GetBinCenter(j);
      double fit = f2->Eval(eta, pt);
      double xsec = h2->GetBinContent(i,j);
      double ex = h2->GetBinError(i,j);
      if (pt*cosh(eta)<13600.*0.5) {
	h2rf->SetBinContent(i, j, xsec / fit);
	h2rf->SetBinError(i, j, ex / fit);
      }
    } // for j
  } // for i

  gPad->SetLogy();
  h2rf->Draw("SAME COLZ");
  h2rf->GetZaxis()->SetRangeUser(0.,2.);//0.75+eps,1.25-eps);
  h2rf->GetZaxis()->SetTitle(Form("%s / Fit",cr));
  h2rf->GetZaxis()->SetTitleOffset(1.6);
  h2rf->GetZaxis()->SetTitleSize(0.05);
  h2rf->GetZaxis()->SetLabelSize(0.05);

  gPad->RedrawAxis();

  // Draw lines of "safe space"
  l->DrawLine(-2.5,30,+2.5,30);
  l->DrawLine(-2.5,30,-2.5,1109);
  l->DrawLine(+2.5,30,+2.5,1109);
  f1->Draw("SAME");
  
  gPad->Update();
  c4->SaveAs(Form("pdf/timeDep2D/timeDep2D_fitratio_%s.pdf",cr));


  // Now turn the results into actual JECs
  // We need:
  // - reference 2D spectrum fitted from 2025CDEFG
  //   - f2 in this script
  // - xsec ratio of a given era to 2025CDEFG
  //   - h2r in this script
  // - JES applied to 2025CDEFG as baseline
  //   - jecdata2025CDEFG.root:ratio/eta00-13/run3/hFit_Rjet [was herr_l2l3res]
  //   - L2Res.root:h2jes1_(const_eta)_2025CDEFG
  // Store result as TProfile2D with same statistics as inclusive jet
  // cross section in era to easy slicing and dicing later
  bool extractJES(true);
  if (extractJES) {

    cout << "Entering extractJES  for run " << run << endl << flush;
    TFile *fl2res = new TFile("rootfiles/L2Res_2024_V9M_2025_V3M_noincjet.root","READ");
    TH2D *h2res = (TH2D*)fl2res->Get("h2res1_2025CDEFG");
    TH2D *h2refit = (TH2D*)fl2res->Get("h2jes1_const_eta_2025CDEFG");
    assert(h2res);

    //TFile  *fl3res = new TFile("rootfiles/jecdata2025CDEFG.root","READ");
    //TFile  *fl3res = new TFile("rootfiles/jecdata2025CDEFG_V3M_noincjet_v2.root","READ");
    TFile  *fl3res = new TFile("rootfiles/jecdata2025CDEFG_V3M_noincjet_v3.root","READ");
    TH1D *h1res = (TH1D*)fl3res->Get("ratio/eta00-13/herr_l2l3res");
    assert(h1res);
    TH1D *h1refit = (TH1D*)fl3res->Get("ratio/eta00-13/run3/hFit_Rjet");
    assert(h1refit);

    TFile *fout = new TFile("rootfiles/timeDep2D.root","UPDATE");
    TH1D *h1rel = (hr!=0 ? (TH1D*)hr->Clone(Form("h1rel_%s",cr)) : 0);
    if (h1rel) h1rel->Reset();
    TH1D *h1jes = (hr!=0 ? (TH1D*)hr->Clone(Form("h1jes_%s",cr)) : 0);
    if (h1jes) h1jes->Reset();
    TH2D *h2rel = (TH2D*)h2r->Clone(Form("h2rel_%s",cr)); h2rel->Reset();
    TH2D *h2jes = (TH2D*)h2r->Clone(Form("h2jes_%s",cr)); h2jes->Reset();
    curdir->cd();

    _f2 = f2;
    // Calculate xsec ratio given JES shift and [etamin,etamax,ptmin,ptmax]
    TF1 *f1jes = new TF1("f1jes",fjes,0.5,1.5,4);
    // Full 1D JES in barrel
    for (int ipt = 1; h1jes && ipt != h1jes->GetNbinsX()+1; ++ipt) {
      
      double eta = 0;
      double pt = h1jes->GetXaxis()->GetBinCenter(ipt);
      double etamax = 1.305;
      double ptmax = h1jes->GetXaxis()->GetBinLowEdge(ipt+1);
      double ergmax = pt*cosh(eta);//ptmax*cosh(etamax);
      double sqrts = 13600.;
      
      // Don't try to solve this for bins extending out of phase space
      if (ergmax>0.5*sqrts) continue;
      
      double r = hr->GetBinContent(ipt);
      double er = hr->GetBinError(ipt);
      
      f1jes->SetParameters(-1.305,1.305,
			   h1jes->GetXaxis()->GetBinLowEdge(ipt),
			   h1jes->GetXaxis()->GetBinLowEdge(ipt+1));
      
      // Skip bins that are clearly out of bounds (e.g. due to low stats)
      double r_max = f1jes->Eval(1.5);
      double r_min = f1jes->Eval(0.5);
      if (r-er<r_min || r+er>r_max) continue;
      
      double jes = f1jes->GetX(r,0.5,1.5);//,1e-4);
      double jes_up = f1jes->GetX(r+er,0.5,1.5,1e-4);
      double jes_dw = f1jes->GetX(r-er,0.5,1.5,1e-4);
      double ejes = 0.5*(fabs(jes_up-jes)+fabs(jes_dw-jes));

      // Residual corrections from previous round
      double res2 = 1.0;
      //double res3 = h1res->GetBinContent(h1res->GetXaxis()->FindBin(pt));
      double res3 = h1res->Interpolate(pt);

      // Residual corrections from refit
      double refit2 = 1.0;
      double refit3 = h1refit->Interpolate(pt);

      // Relative scale with respect to 2025CDEFG with latest JEC
      h1rel->SetBinContent(ipt, jes);
      h1rel->SetBinError(ipt, ejes);

      // Same, but undoing previous JEC so can do both closure and refit
      h1jes->SetBinContent(ipt, jes * (refit2 * refit3) / (res2 * res3));
      h1jes->SetBinError(ipt, ejes * (refit2 * refit3) / (res2 * res3));
    } // for ipt
    
    // Full 2D JES
    for (int ieta = 1; ieta != h2jes->GetNbinsX()+1; ++ieta) {
      for (int ipt = 1; ipt != h2jes->GetNbinsY()+1; ++ipt) {
	
	double eta = h2jes->GetXaxis()->GetBinCenter(ieta);
	double pt = h2jes->GetYaxis()->GetBinCenter(ipt);
	double etamax = h2jes->GetXaxis()->GetBinLowEdge(ieta+1);
	double ptmax = h2jes->GetYaxis()->GetBinLowEdge(ipt+1);
	double ergmax = pt*cosh(eta);//ptmax*cosh(etamax);
	double sqrts = 13600.;

	// Don't try to solve this for bins extending out of phase space
	if (ergmax>0.5*sqrts) continue;
	
	double r = h2r->GetBinContent(ieta, ipt);
	double er = h2r->GetBinError(ieta, ipt);

	f1jes->SetParameters(h2jes->GetXaxis()->GetBinLowEdge(ieta),
			     h2jes->GetXaxis()->GetBinLowEdge(ieta+1),
			     h2jes->GetYaxis()->GetBinLowEdge(ipt),
			     h2jes->GetYaxis()->GetBinLowEdge(ipt+1));

	/*
	// debugging tests
	if (ptmax>97 && ptmax<=114 && etamax>0 && etamax<=0.087) {
	  cout << "f1jes->Eval(0.9) = " << f1jes->Eval(0.9) << endl << flush;
	  cout << "f1jes->Eval(1.0) = " << f1jes->Eval(1.0) << endl << flush;
	  cout << "f1jes->Eval(1.1) = " << f1jes->Eval(1.1) << endl << flush;
	}
	else continue;
	*/
	
	// Skip bins that are clearly out of bounds (e.g. due to low stats)
	double r_max = f1jes->Eval(1.5);
	double r_min = f1jes->Eval(0.5);
	if (r-er<r_min || r+er>r_max) continue;
	
	double jes = f1jes->GetX(r,0.5,1.5);//,1e-4);
	double jes_up = f1jes->GetX(r+er,0.5,1.5,1e-4);
	double jes_dw = f1jes->GetX(r-er,0.5,1.5,1e-4);
	double ejes = 0.5*(fabs(jes_up-jes)+fabs(jes_dw-jes));

	// Residual corrections from previoius round
	//double res2 = h2res->GetBinContent(h2res->GetXaxis()->FindBin(eta),
	//				   h2res->GetYaxis()->FindBin(pt));
	double eta2 = min(max(eta, h2res->GetXaxis()->GetBinLowEdge(1)),
			  h2res->GetXaxis()->GetBinLowEdge(h2res->GetNbinsX()+1));
	double pt2 = min(max(pt, h2res->GetYaxis()->GetBinLowEdge(1)),
			  h2res->GetYaxis()->GetBinLowEdge(h2res->GetNbinsY()+1));
	double res2 = h2res->Interpolate(eta2,pt2);
	//double res3 = h1res->GetBinContent(h1res->GetYaxis()->FindBin(pt));
	double res3 = h1res->Interpolate(pt);

	// Residual corrections from refit
	double eta3 = min(max(eta, h2refit->GetXaxis()->GetBinLowEdge(1)),
			  h2refit->GetXaxis()->GetBinLowEdge(h2refit->GetNbinsX()+1));
	double pt3 = min(max(pt, h2refit->GetYaxis()->GetBinLowEdge(1)),
			 h2refit->GetYaxis()->GetBinLowEdge(h2refit->GetNbinsY()+1));
	double refit2 = h2refit->Interpolate(eta3,pt3);
	double refit3 = h1refit->Interpolate(pt);
					   
	h2rel->SetBinContent(ieta, ipt, jes);
	h2rel->SetBinError(ieta, ipt, ejes);
	h2jes->SetBinContent(ieta, ipt, jes * (refit2*refit3) / (res2*res3));
	h2jes->SetBinError(ieta, ipt, ejes * (refit2*refit3) / (res2*res3));
      } // for ipt
    } // for ieta
    fout->Write("",TObject::kOverwrite);
    fout->Close();
  }
  
} // timeDep2Ds

void timeDep2D() {

  gROOT->ProcessLine(".! mkdir pdf");
  gROOT->ProcessLine(".! mkdir pdf/timeDep2D");
  gROOT->ProcessLine(".! touch pdf");
  gROOT->ProcessLine(".! touch pdf/timeDep2D");

  ///////////////////////////////////////
  // Start with 1D spectra as a warmup //
  ///////////////////////////////////////

  string runs[] =
    // v155+v153_copy
    {"2025CDEFG",
     "2025C","2025D","2025E","2025F","2025G",
     "2024C_nib1","2024D_nib1",              "2024E_nib1",
     "2024F_nib1","2024F_nib2","2024F_nib3","2024G_nib1","2024G_nib2",
     "2024H_nib1","2024I_nib1"};
  const int nrun = sizeof(runs)/sizeof(runs[0]);

  map<string,int> marker;
  marker["2025CDEFG"] = kNone;
  marker["2025C"] = kFullDiamond;
  marker["2025D"] = kFullDiamond;
  marker["2025E"] = kFullDiamond;
  marker["2025F"] = kFullDiamond;
  marker["2025G"] = kFullDiamond;

  marker["2024C_nib1"] = kOpenDiamond;
  marker["2024D_nib1"] = kOpenDiamond;
  marker["2024E_nib1"] = kOpenDiamond;
  marker["2024F_nib1"] = kOpenDiamond;
  marker["2024F_nib2"] = kOpenDiamond;
  marker["2024F_nib3"] = kOpenDiamond;
  marker["2024G_nib1"] = kOpenDiamond;
  marker["2024G_nib2"] = kOpenDiamond;
  marker["2024H_nib1"] = kOpenDiamond;
  marker["2024I_nib1"] = kOpenDiamond;
  
  map<string,int> color;
  color["2025CDEFG"] = kBlack;
  color["2025C"] = kRed;
  color["2025D"] = kOrange+2;
  color["2025E"] = kGreen+2;
  color["2025F"] = kMagenta+2;
  color["2025G"] = kBlack;

  color["2024C_nib1"] = kRed;
  color["2024D_nib1"] = kOrange+2;
  color["2024E_nib1"] = kGreen+2;
  color["2024F_nib1"] = kBlue-9;
  color["2024F_nib2"] = kMagenta;
  color["2024F_nib3"] = kMagenta+2;
  color["2024G_nib1"] = kCyan;
  color["2024G_nib2"] = kCyan+2;
  color["2024H_nib1"] = kBlue;
  color["2024I_nib1"] = kBlack;
  
  // Normalized spectra for each eta
  TH1D *href(0);
  vector<TH1D*> hrun(nrun);
  for (int i = 0; i != nrun; ++i) {
    hrun[i] = timeDep1Ds(runs[i],triggerLumi(runs[i],"HLT_PFJet500")*1e-3);
    if (!href) href = (TH1D*)hrun[0]->Clone("href");
  }
  
  TH1D *h = tdrHist("h","Cross section (pb / GeV)",0.2e-7,2e11,
		    "p_{T,jet} (GeV)",15,4500.);
  //TH1D *hd = tdrHist("hd","Ratio to 2024",0.50,1.50,
  //TH1D *hd = tdrHist("hd","Ratio to 2025",0.50,1.50,
  TH1D *hd = tdrHist("hd","Ratio to 2025",0.90,1.15,
		     "p_{T,jet} (GeV)",15,4500.);
  extraText = "Private";
  //lumi_136TeV = Form("2024, %1.1f fb^{-1}",lumsum);
  lumi_136TeV = Form("2025, %1.1f fb^{-1}", triggerLumi("2025CDEFG","HLT_PFJet500")*1e-3);
  TCanvas *c1 = tdrDiCanvas("c1",h,hd,8,11);

  c1->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  const int nmax1(5);
  TLegend *leg1 = tdrLeg(0.60,0.90-min(nmax1,(nrun-1))*0.05,0.85,0.90);
  TLegend *leg2 = tdrLeg(0.17,0.03,0.42,0.03+max(0,(nrun-1)-nmax1)*0.045);

  TLatex *tex = new TLatex();
  tex->SetNDC(); tex->SetTextSize(0.045);
  tex->DrawLatex(0.35,0.85,"|#eta| < 1.305");

  c1->cd(2);
  gPad->SetLogx();

  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(15,1,4500,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(15,1.05,4500,1.05);
  l->DrawLine(15,0.95,4500,0.95);

  c1->cd(1);
  TH1D *href2 = (TH1D*)href->Clone("href2");
  for (int i = 1; i != href2->GetNbinsX()+1; ++i) {
    double pt = href2->GetBinCenter(i);
    if (href2->GetBinContent(i)!=0) {
      double err = href->GetBinError(i);
      double y = href->GetBinContent(i);
      double k = 0.0075;
      href2->SetBinError(i, sqrt(pow(err,2)+pow(k*y,2)));
    }
  }
  /*
  TF1 *fref = new TF1("fref","[0]*exp(-[7]/x)*pow(x,[1]+[3]*log10(0.01*x)+"
		      "[4]*pow(log10(0.01*x),2)+[5]*pow(log10(0.01*x),3))*"
		      "pow(1-2*x*cosh([6]+[8]*log10(0.001*x))/13600.,[2])",15,3500);
  fref->SetParameters(7.4e14,-5,9.8,0.,0.,0.,0.44,38.,0.);
  //fref->FixParameter(6,0);
  fref->FixParameter(5,0);
  fref->FixParameter(4,0);
  fref->FixParameter(3,0);
  */
  TF1 *fref = new TF1("fref","[0]*exp(-[1]/x)*pow(x,[2]+[3]*log(0.01*x))*"
		      "pow(1-2*x*cosh([4]+[5]*log(0.001*x))/13600.,[6])",
		      18,3500);
  //fref->SetParameters(7.4e14, 38., -5,0., 0.44,0., 9.8);
  fref->SetParameters(5.76e14, 38.5, -4.82,0., 1.01,-0.47, 10.77); // [15,3500]
  fref->FixParameter(3,0);
  //fref->FixParameter(5,0);
  href2->Fit(fref,"RN");
  fref->Draw("SAME");

  c1->cd(2);
  TH1D *hfrefr = (TH1D*)href->Clone("hfrefr");
  hfrefr->Divide(href);
  for (int i = 1; i != hfrefr->GetNbinsX()+1; ++i) {
    double pt = hfrefr->GetBinCenter(i);
    if (href->GetBinContent(i)!=0)
      hfrefr->SetBinContent(i, fref->Eval(pt)/href->GetBinContent(i));
  }
  tdrDraw(hfrefr,"HIST",kNone,kRed,kSolid,-1,kNone,0);

  TH1D* vhr[nrun];
  for (int i = 0; i != nrun; ++i) {
    string r = runs[i];
    const char *cr = r.c_str();
    double lum = triggerLumi(r, "HLT_PFJet500")*1e-3;
    TH1D *h = hrun[i];
    
    c1->cd(1);
    tdrDraw(h, i==0 ? "HIST" : "Pz", marker[r], color[r], kSolid, -1, kNone);
    TLegend *leg = (i<nmax1+1 ? leg1 : leg2);
    if (i!=0)
      leg->AddEntry(h,Form("%s, %1.1f fb^{-1}",cr,lum), i==0 ? "FL" : "PLE");

    c1->cd(2);
    TH1D *hr = (TH1D*)h->Clone(Form("hr_%s",cr));
    hr->Divide(href);
    tdrDraw(hr, i==0 ? "HIST" : "Pz", marker[r], color[r], kSolid, -1, kNone);

    vhr[i] = hr;
  }
  
  c1->cd(1);
  gPad->RedrawAxis();

  c1->cd(2);
  gPad->RedrawAxis();

  c1->SaveAs("pdf/timeDep2D/timeDep2D_1D.pdf");


  TH1D *h1_b = tdrHist("h1_b","Cross section (pb / GeV)",0.2e-7,2e11,
		       "p_{T,jet} (GeV)",15,4500.);
  TH1D *h1d_b = tdrHist("h1d_b","Ratio to fi",0.90,1.15,
			"p_{T,jet} (GeV)",15,4500.);
  extraText = "Private";
  lumi_136TeV = Form("2025, %1.1f fb^{-1}", triggerLumi("2025CDEFG","HLT_PFJet500")*1e-3);
  TCanvas *c1b = tdrDiCanvas("c1b",h1_b,h1d_b,8,11);

  c1b->cd(1);
  gPad->SetLogx();
  gPad->SetLogy();

  tex->DrawLatex(0.35,0.85,"|#eta| < 1.305");

  c1b->cd(2);
  gPad->SetLogx();
  /*
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->SetLineColor(kGray+2);
  l->DrawLine(15,1,4500,1);
  l->SetLineStyle(kDotted);
  l->DrawLine(15,1.05,4500,1.05);
  l->DrawLine(15,0.95,4500,0.95);

  for (int i = 0; i != nrun; ++i) {
    string r = runs[i];
    const char *cr = r.c_str();
    double lum = triggerLumi(r, "HLT_PFJet500")*1e-3;
    TH1D *h = hrun[i];
    
    c1->cd(1);
    tdrDraw(h, i==0 ? "HIST" : "Pz", marker[r], color[r], kSolid, -1, kNone);
    if (i!=0) leg->AddEntry(h,Form("%s, %1.1f fb^{-1}",cr,lum), i==0 ? "FL" : "PLE");

    c1->cd(2);
    TH1D *hr = (TH1D*)h->Clone(Form("hr_%s",cr));
    hr->Divide(href);
    tdrDraw(hr, i==0 ? "HIST" : "Pz", marker[r], color[r], kSolid, -1, kNone);
  }
  
  c1->cd(1);
  gPad->RedrawAxis();

  c1->cd(2);
  gPad->RedrawAxis();

  c1->SaveAs("pdf/timeDep2D/timeDep2D_1D.pdf");
  */


  

  ///////////////////
  // On to full 2D //
  ///////////////////
  /*
  timeDep2Ds("2024BCD");
  timeDep2Ds("2024E");
  timeDep2Ds("2024F");
  timeDep2Ds("2024G");
  //timeDep2Ds("2024FG");
  */

  for (int i = 0; i != nrun; ++i) {
    timeDep2Ds(runs[i],"2025CDEFG",vhr[i]);
  }
}
