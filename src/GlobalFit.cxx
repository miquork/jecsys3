#include "GlobalFit.hpp"
#include "Containers.hpp"
#include "utils.hpp"
#include "config.hpp"

using namespace std;

// Set to a dummy value
TString GlobalFit::current_obs = "";

GlobalFit::GlobalFit(std::string runName_, std::string mode_, double eta_min_,double eta_max_) {
  // Constructor set the run, mode, eta range and open files
  runName = runName_; mode = mode_;
  eta_min = Form("%02.0f",10*eta_min_);
  eta_max = Form("%02.0f",10*eta_max_);
  PrintLine("Running on: run="+runName+", mode="+mode+", eta_min="+eta_min+", eta_max="+eta_max, green);
  Print("  --> samples:", green); PrintLine(samples);
  Print("  --> hdm_methods:", green); PrintLine(hdm_methods);
  Print("  --> types:", green); PrintLine(types);
  for (auto sample: samples) {
    for (auto method: hdm_methods) {
      for (auto type: types) {
        input_hnames.push_back(type+"_"+sample+"_"+method);
      }
    }
  }

  OpenFiles();
}

GlobalFit::~GlobalFit() {
  // Destructor takes care of closing files
  for (auto [fname, f_]: input_files) {
    if (debug) PrintLine("Closing file for "+fname, yellow);
    f_->Close();
  }
  for (auto [fname, f_]: output_files) {
    if (debug) PrintLine("Closing file for "+fname, yellow);
    f_->Close();
  }
}


void GlobalFit::OpenFiles(){
  // Open input/output files
  for (auto [name, fname]: input_fnames) {
    if (debug) PrintLoading("file", name, fname);
    input_files[name] = new TFile(ReplaceDefault(fname),"UPDATE");
    assert(input_files[name] && !input_files[name]->IsZombie());
    if (debug) PrintLine("successfully", yellow);
  }
  output_files["output"] = new TFile(ReplaceDefault(output_fname), "RECREATE");
}


TString GlobalFit::ReplaceDefault(TString name) {
  // Function to update the default strings
  TString res = name;
  res = res.ReplaceAll("RUN",runName).ReplaceAll("MODE",mode);
  res = res.ReplaceAll("ETAMIN",eta_min).ReplaceAll("ETAMAX",eta_max);
  return res;
};


void GlobalFit::LoadInputs(){
  for (auto [name, info]: input_hnames_map) {
    if (FindInVector(input_hnames,name)<0) {
      PrintLine("Skipping hist: "+name, yellow);
      continue;
    }
    TString hname = ReplaceDefault(info["hname"]);
    if (debug) PrintLoading("hist", name, hname);
    TGraphErrors* g;
    if (info["type"]=="Resp"){
      unique_ptr<TH1D> h; h.reset((TH1D*)input_files[info["fname"].Data()]->Get(hname));
      assert(h);
      g = new TGraphErrors(h.get());
    } else {
      g = (TGraphErrors*)input_files[info["fname"].Data()]->Get(hname);
    }
    RemoveZerosFromGraph(g);
    my_data[name] = new DataContainer(name, info["type"], hname, g);
    nTotalPoints += my_data[name]->GetN();
    if (debug) {PrintLine("successfully", yellow); cout << *my_data[name] << endl;}
  }


  for (auto [name, info]: reference_obj_map) {
    TString hname = ReplaceDefault(info["hname"]);
    if (debug) PrintLoading("hist", name, hname);
    reference_objects[name] =(TH1D*)input_files[info["fname"].Data()]->Get(hname)->Clone(name);
    assert(reference_objects[name]);
    reference_objects[name]->SetDirectory(0);
    if (debug) {PrintLine("successfully", yellow);}
  }
}

void GlobalFit::LoadSystematics(){
  for (auto [name, info]: sources_hnames_map) {
    TString appliesTo = info["appliesTo"];
    if (FindInVector(input_hnames,appliesTo)<0) {
      PrintLine("Skipping source: "+name+" index: "+to_string(FindInVector(input_hnames,appliesTo)), yellow);
      continue;
    }
    TString hname = ReplaceDefault(info["hname"]);
    if (debug) PrintLoading("source", name, hname);
    unique_ptr<TH1D> h; h.reset((TH1D*)input_files[info["fname"].Data()]->Get(hname));
    sources[name] = new SystematicContainer(name,appliesTo, atoi(info["index"]), hname, h.get());
    assert(sources[name]);
    if (debug) {PrintLine("successfully", yellow); cout << *sources[name] << endl;};
  }
}

void GlobalFit::LoadShapes(){
  for (auto [name, info]: shapes_map) {
    if (debug) PrintLoading("shape", name, "");
    shapes[name] = new ShapeContainer(name,info["form"],info["appliesTo"],atoi(info["index"]),atoi(info["ispositive"]));
    assert(shapes[name]);
    if (debug) {PrintLine("successfully", yellow); cout << *shapes[name] << endl;};
    shape_types.insert(info["type"]);
  }
}


void GlobalFit::SetupFitFunction(){
  // Set the minimizer tool for the global fit
  if (debug) PrintLine("Loaded "+to_string(my_data.size())+" hists", green);
  if (debug) PrintLine("Loaded "+to_string(shapes.size())+" shapes", green);
  if (debug) PrintLine("Loaded "+to_string(sources.size())+" sources", green);

  nFitPars = shape_types.size();
  nNuisancePars = sources.size();
  nTotPars = nFitPars+nNuisancePars;
  PrintLine("Global fit has "+to_string(nTotPars)+" total parameters:", blue);
  PrintLine("  --> "+to_string(nFitPars)+" fit parameters", blue);
  PrintLine("  --> "+to_string(nNuisancePars)+" nuisance parameters", blue);
  if (penalizeFitPars) PrintLine("  --> Fit parameters have Gaussian prior", blue);

  auto jesFit_wrapper_ = jesFit_wrapper( reference_objects["hjesref"],shapes);
  _jesFit = new TF1("jesFit", jesFit_wrapper_ ,fit_min,fit_max,nFitPars);
  // TMinuit *fitter = new TMinuit(nTotPars);
  fitter = new TFitter(nTotPars);
  fitter->SetFCN(jesFitter);
  for (int i = 0; i != nTotPars; ++i) {
    fitter->SetParameter(i, "", 0, (i<nFitPars ? 0.01 : 1),-100, 100);
  }

}

void GlobalFit::DoGlobalFit(){

  // Run fitter (multiple times if needed)
  for (int i = 0; i != number_of_fit_iterations; ++i) {
    fitter->ExecuteCommand("MINI", 0, 0);
  }

  // Verify that the degrees of freedom make sense. Important to check immediately
  nFittedDataPoints -= nNuisancePars;
  if (penalizeFitPars) nFittedDataPoints -= nFitPars;
  assert(nFittedDataPoints==nTotalPoints);

  // Set the error matrix
  error_matrix.reset(new TMatrixD(nTotPars, nTotPars));
  gMinuit->mnemat(error_matrix->GetMatrixArray(), nTotPars);

  // Retrieve the chi2 for the individual components
  Double_t tmp_par[nTotPars], grad[nTotPars];
  Double_t chi2_gbl(0);
  for (int i = 0; i < nTotPars; ++i) tmp_par[i] = fitter->GetParameter(i);
  jesFitter(nTotPars, grad, chi2_gbl, tmp_par, 1);

  // Retrieve the chi2 for the individual components
  Double_t chi2_src(0), chi2_par(0), chi2_data(0);
  int npar_true(0), nsrc_true(0);

  for (int i = 0; i != nTotPars; ++i) {
    double val = fitter->GetParameter(i);
    double err = fitter->GetParError(i);
    if (fabs(val)!=0 || fabs(err-1)>1e-2) {
      if (i < nFitPars) {
        ++npar_true;
        chi2_par += pow(val,2);
      } else {
        ++nsrc_true;
        chi2_src += pow(val,2);
      }
    }
  }

  assert(nFitPars==npar_true);
  assert(nNuisancePars==nsrc_true);

  for (auto [dt_name,dt]: my_data){
    TGraphErrors *graph_output   = dt->output();  assert(graph_output);
    current_obs = dt->type(); // this is needed inside _jesFit
    for (int j = 0; j != graph_output->GetN(); ++j) {
      double x = graph_output->GetX()[j];
      double y = graph_output->GetY()[j];
      double ey = graph_output->GetEY()[j];
      chi2_data += pow((y - _jesFit->Eval(x)) / ey, 2);
    }
  }
  // Verify that chi2 and d.o.f. make sense
  assert(chi2_gbl=chi2_data+chi2_src+chi2_par);
  assert(nFittedDataPoints=nTotalPoints+nFitPars+nNuisancePars);

  // Summary of the global fit
  PrintLine("Output GlobalFit", blue);
  PrintLine("  --> "+to_string(nTotalPoints)+" data points used", blue);
  PrintLine("  --> "+to_string(nFitPars)+" fit parameters", blue);
  PrintLine("  --> "+to_string(nsrc_true)+" nuisances", blue);
  PrintLine(Form("  --> fitting range [%1.0f,%1.0f]", _jesFit->GetXmin(),_jesFit->GetXmax()), blue);
  PrintLine(Form("  --> Total     chi2/NDF  = %1.1f / %d",chi2_gbl, nFittedDataPoints), blue);
  PrintLine(Form("  --> Data      chi2/NDF  = %1.1f / %d",chi2_data, nTotalPoints), blue);
  PrintLine(Form("  --> Nuisance  chi2/Nsrc = %1.1f / %d",chi2_src, nFitPars), blue);
  PrintLine(Form("  --> Parameter chi2/Npar = %1.1f / %d",chi2_par, nNuisancePars), blue);

  PrintLine("Fitted parameters (for Resp):", blue);
  for (auto [name,shape]: shapes){
    if (shape->appliesTo()!="Resp") continue;
    double val = fitter->GetParameter(shape->index());
    double err = fitter->GetParError(shape->index());
    PrintParameter(name, val, err);
  }

  PrintLine("Nuisance parameters:", blue);
  for (auto [name,src]: sources){
    double val = fitter->GetParameter(nFitPars+src->index());
    double err = fitter->GetParError(nFitPars+src->index());
    PrintParameter(name, val, err);
  }
}

void GlobalFit::StoreFitOutput(){
  // Store input and output into file
  output_files["output"]->cd();

  // Set range and granularity for plotting
  _jesFit->SetRange(func_range_min,func_range_max);
  _jesFit->SetNpx(func_range_max-func_range_min);

  // Convert global fit function into hist and graph
  double binning[] = { /*1, 5, 6, 8,*/ 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
  const double nbins = sizeof(binning)/sizeof(binning[0])-1;
  for (auto type: all_types){
    current_obs = type; // this is needed inside _jesFit
    // convert TF1 to TGraphErrors
    TGraphErrors* graph = new TGraphErrors();
    FuncToGraph(_jesFit,*error_matrix.get(),graph);
    // convert TF1 to TGraphErrors
    TH1D *hist = new TH1D("jesFit_hist_"+type,";p_{T} (GeV);"+type,nbins,binning);
    hist->SetDirectory(0);
    FuncToHist(_jesFit,*error_matrix.get(),hist);
    // Store
    graph->Write("jesFit_graph_"+type,TObject::kOverwrite);
    hist->Write(hist->GetName(),TObject::kOverwrite);
    _jesFit->Write("jesFit_func_"+type,TObject::kOverwrite);
  }
  error_matrix->Write("error_matrix",TObject::kOverwrite);


  // Store all jes response
  for (auto [name,jes]: my_data){
    TGraphErrors* input = jes->input();
    TGraphErrors* output = jes->output();
    TGraphErrors* variation = jes->variation();
    multiplyGraph(input,1./ScaleFullSimShape); // transform it in percentage
    multiplyGraph(output,1./ScaleFullSimShape); // transform it in percentage
    multiplyGraph(variation,1./ScaleFullSimShape); // transform it in percentage
    std::vector<TF1*> funcs;
    for (auto [name,shape]: shapes){
      if (shape->appliesTo()==jes->type()) funcs.push_back(shape->func());
    }
    PropagateErrorToGraph(variation, funcs, *error_matrix.get());
    input->Write(name+"_prefit",TObject::kOverwrite);
    output->Write(name+"_postfit",TObject::kOverwrite);
    variation->Write(name+"_variation",TObject::kOverwrite);
  }

  // Store all sources
  for (auto [name,src]: sources){
    src->hist()->SetYTitle(name);
    src->hist()->Write("source_"+name,TObject::kOverwrite);
  }

  // Store all shapes: pre/post fit
  for (auto [name,shape]: shapes){
    shape->func()->Write("shape_input_"+name,TObject::kOverwrite);
    TF1* prefit = new TF1("shape_prefit_"+name,"[0]*("+shape->form()+")",func_range_min,func_range_max);
    TF1* postfit = new TF1("shape_postfit_"+name,"[0]*("+shape->form()+")",func_range_min,func_range_max);
    prefit->SetParameter(0,ScaleFullSimShape);
    postfit->SetParameter(0,ScaleFullSimShape*_jesFit->GetParameter(shape->index()));
    prefit->Write(prefit->GetName(),TObject::kOverwrite);
    postfit->Write(postfit->GetName(),TObject::kOverwrite);
  }

  // Store reference objects
  for (auto [name,obj]: reference_objects){
    obj->Write(name,TObject::kOverwrite);
  }


}

void GlobalFit::Run(){
  LoadInputs();
  LoadSystematics();
  LoadShapes();
  SetupFitFunction();
  DoGlobalFit();
  StoreFitOutput();
}


function<double(Double_t *, Double_t *)> jesFit_wrapper(TH1D* hjesref, map<TString, ShapeContainer*> shapes){

  return [hjesref, shapes](Double_t *x, Double_t *p) -> double {

    double var = (GlobalFit::current_obs=="Resp" ? 1. : 0.);
    double pt = x[0];

    double jesref = 1;
    // if (debug) PrintLine("current_obs "+GlobalFit::current_obs,red);
    if (GlobalFit::useJESref && GlobalFit::current_obs=="Resp") {
      assert(hjesref);
      jesref = hjesref->Interpolate(pt);
      // if (debug) PrintLine("update jesref "+to_string(jesref),red);
    }

    // // Load available shapes for this observable
    for (auto [name,shape]: shapes){
      if (shape->appliesTo()!="Resp") continue;
      // Calculate variation and add it to total
      int index = shape->index();   assert(index>=0);
      TF1 *f1 = shape->func();  assert(f1);
      double par = p[index];
      if (shape->ispositive()) par = max(par,0.);
      //if (shape->ispositive()) par = min(1.,max(par,0.));
      var += par * f1->Eval(pt) * GlobalFit::ScaleFullSimShape; // fullSimShapes in %'s
    }

    return (var / jesref);;
  };
}

// Dummy value: not used at the moment but it can be changes if needed
Double_t jesFit(Double_t *x, Double_t *p) {return -1000;}


void jesFitter(Int_t& npar, Double_t* grad, Double_t& chi2, Double_t* par, Int_t flag) {

  // Basic checks
  assert(_jesFit);

  // Parametes for nuisances (sources)
  int _nFitPars = _jesFit->GetNpar();
  int _nNuisancePars = npar - _nFitPars;
  Double_t *fit_pars = &par[0];
  Double_t *nuisance_pars = &par[_nFitPars];

  if (flag) {
    // do the following calculation:
    // chi2 = sum_i( (x_i+sum_s(a_s y_si) -fit)^2/sigma_i^2) + sum_s(a_s^2)
    chi2 = 0;
    nFittedDataPoints = 0;

    // Loop over input data (graphs x bins)
    // - for each point, add up source eigenvectors x nuisance parameters
    // - then calculate chi2 adding up residuals + nuisance parameters
    for (auto [dt_name,dt]: my_data){
      cout << *dt << endl;
      TString dt_type         = dt->type();
      TGraphErrors *graph_input   = dt->input();  assert(graph_input);
      // TGraphErrors *graph_input2  = dt->input2;
      TGraphErrors *graph_output    = dt->output(); assert(graph_output);
      TGraphErrors *graph_variation = dt->variation(); assert(graph_variation);
      // TGraphErrors *graph_output2 = dt->output2;

      for (int bin = 0; bin < graph_input->GetN(); ++bin) {

        // Retrieve central value and uncertainty for this point
        double pt = graph_input->GetX()[bin];
        double data_point = graph_input->GetY()[bin];
        double sigma = graph_input->GetEY()[bin];

        // Calculate fit value at this point
        GlobalFit::current_obs = dt_type;
        _jesFit->SetParameters(fit_pars);
        double fit = _jesFit->EvalPar(&pt,par);

        // // For multijet balancing, multiply data by reference JES
        // if (dt_name.Contains("multijet")) {
        //   // assert(graph_input2);
        //   double ptref = graph_input->GetY()[bin] * pt;
        //   double fitRef = _jesFit->EvalPar(&ptref,par);
        //   data_point *= fitRef;
        //   sigma *= fitRef;
        // } // multijet
        //
        // Calculate total shift caused by all nuisance parameters
        double shift = 0;
        // vector<SystematicContainer> &_vsrc = _msrc[dt_name];
        for (auto [src_name_,source]: sources){
          if (dt_name!=source->appliesTo()) continue;
          if (debug) PrintLine("adding "+src_name_+" to shift for "+ dt_name+"("+source->appliesTo()+")", yellow);
          TH1D* hsrc = source->hist(); assert(hsrc);
          shift += nuisance_pars[source->index()] * hsrc->GetBinContent(hsrc->FindBin(pt));
        }

        for (auto [shape_name,shape]: shapes){
          if ("Resp"==shape->appliesTo()) continue;
          if (dt_type!=shape->appliesTo()) continue;
          if (pt<40 || pt>600) continue;
          if (debug) PrintLine("adding "+shape_name+" to shift for "+ dt_type+"("+shape->appliesTo()+")", yellow);
          TF1* func = shape->func(); assert(func);
          cout << "Adding to shift: " << nuisance_pars[shape->index()] << " " <<  func->Eval(pt) << " " << nuisance_pars[shape->index()] * func->Eval(pt)* GlobalFit::ScaleFullSimShape << " " << nuisance_pars[shape->index()] * func->Eval(pt) / GlobalFit::ScaleFullSimShape <<  endl;;
          shift += nuisance_pars[shape->index()] * func->Eval(pt) * GlobalFit::ScaleFullSimShape;
        }

        // Add chi2 from residual
        cout << "data_point: " << data_point << " shift: " << shift << " fit: " << fit << " sigma: " << sigma << " err: " << oplus(sigma,GlobalFit::globalErrMin) << endl;
        if (debug) PrintLine("chi2 before point: "+to_string_with_precision(chi2),yellow);
        double chi = (data_point + shift - fit) / oplus(sigma,GlobalFit::globalErrMin);
        chi2 += chi * chi;
        ++nFittedDataPoints;
        //
        // // Store shifted data
        assert(graph_output->GetN()==graph_input->GetN() && graph_output->GetX()[bin]==pt);
        graph_output->SetPoint(bin, pt, data_point + shift);
        graph_variation->SetPoint(bin, pt, shift);
        //
        // // For multijets, store also downward extrapolation
        // if (dt_name.Contains("multijet")) {
        //   double jes = _jesFit->EvalPar(&pt, par);
        //   double ptref = pt * graph_input->GetY()[bin];
        //   double jesref = _jesFit->EvalPar(&ptref, par);
        //   // MJB = jes / jesref
        //   // data = MJB*jesref => "jesref" = jes / MJB = jes * jesref/data
        //   // NB: double check the logic here
        //   graph_output->SetPoint(bin, ptref, jes * jesref / (data_point + shift));
        //   graph_output->SetPointError(bin, graph_output->GetEX()[bin], graph_output->GetEY()[bin]);
        // } // multijet
      } // for point in graph
    } // for graph

    // Add chi2 from nuisance parameters
    if (debug) PrintLine("chi2 before nuis: "+to_string_with_precision(chi2),yellow);
    for (int ipar = 0; ipar < _nNuisancePars; ++ipar) {
      chi2 += nuisance_pars[ipar]*nuisance_pars[ipar];
      ++nFittedDataPoints;
    }
    // if (debug) PrintLine("chi2 before pars: "+to_string_with_precision(chi2),yellow);
    // Add penalty for fit parameters (Bayesian prior, essentially)
    if (GlobalFit::penalizeFitPars) {
      for (int ipar = 0; ipar < _nFitPars; ++ipar) {
        chi2 += fit_pars[ipar] * fit_pars[ipar];
        ++nFittedDataPoints;
      }
    } // penalizeFitPars

    // if (debug) PrintLine("chi2 before next point+to_string_with_precision(chi2),yellow);

    // Give some feedback on progress in case loop gets stuck
    if ((++fit_counter)%100==0) cout << "." << flush;
  } // if flag
  else {
    if (grad) {}; // suppress warning;
    return;
  }

  cout << yellow << "chi2: " << chi2 << reset << endl;

} // jesFitter
