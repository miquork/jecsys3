#pragma once

// Purpose: Perform Run 3 or Run 2 Legacy global fit
//
// Pre-requisites:
// - reprocess.C : produce rootfiles/jecdata[X].root input file
// - softrad3.C : produce FSR+ISR corrections for HDM (MPF+DB)
// - globalFitSyst.C : produce uncertainty shapes
// - globalFitSettings.h : input definitions and configurable settings
// Post-processing:
// - globalFitPulls.C : plot pull distributions
// - createL2L3ResTextFile.C : produce simplified (L2)L3Res text file
// [- minitools/mergerL2L3ResTextFiles.C : combine L3Res with L2Res]
//
// Author: Mikko Voutilainen, Andrea Malara
//
// Notes: enable systematic source offsetting?

#include "constants.hpp"
#include "config.hpp"
#include "Containers.hpp"

#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFitter.h"
#include "TLine.h"

// #include "tools.C"
// #include "tdrstyle_mod22.C"
// #include "globalFitSettings.h"



class GlobalFit {

public:
  GlobalFit(std::string runName_, std::string mode_, double eta_min_,double eta_max_);
  ~GlobalFit();
  void OpenFiles();
  void Run();
  void LoadInputs();
  void LoadSystematics();
  void LoadShapes();
  void SetupFitFunction();
  void DoGlobalFit();
  void StoreFitOutput();
  // Utils
  TString ReplaceDefault(TString name);

  const int number_of_fit_iterations = 1;
  double fit_min = 15.; // Define func range
  double fit_max = 2000; // Define func range
  static const bool useJESref = true;//false; // Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
  static const bool penalizeFitPars = true;
  static constexpr double globalErrMin = 0;
  static TString current_obs;
  static constexpr double ScaleFullSimShape = 0.01;

  std::string runName, mode, eta_min, eta_max;

  std::map<TString, TFile*> input_files;
  std::map<TString, TFile*> output_files;
  std::vector<TString> input_hnames;
  std::set<TString> shape_types;
  std::map<TString, TH1D*> reference_objects;

  TFitter *fitter;
  std::unique_ptr<TMatrixD> error_matrix;

  int nFitPars, nNuisancePars, nTotPars, nTotalPoints=0;

};


Double_t jesFit(Double_t *x, Double_t *p); // Generic fit function for JES and composition
std::function<double(Double_t *, Double_t *)> jesFit_wrapper(TH1D* g, std::map<TString, ShapeContainer*> shapes);

void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par, Int_t flag=1);
// std::function<void(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par, Int_t flag)> jesFitter_wrapper(std::map<TString, ShapeContainer*> shapes, std::map<TString, DataContainer*> my_data, std::map<TString, SystematicContainer*> sources);

static std::map<TString, ShapeContainer*> shapes;
static std::map<TString, DataContainer*> my_data;
static std::map<TString, SystematicContainer*> sources;

static TF1* _jesFit;

static int nFittedDataPoints(0), fit_counter(0);
