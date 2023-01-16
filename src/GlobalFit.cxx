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
    unique_ptr<TGraphErrors> g; g.reset((TGraphErrors*)input_files[info["fname"].Data()]->Get(hname));
    my_data[name] = new DataContainer(name, info["type"], hname, g.get());
    nTotalPoints += my_data[name]->GetN();
    if (debug) {PrintLine("successfully", yellow); cout << *my_data[name] << endl;}
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
    unique_ptr<TH1D> h;

    if (FindInString("scale", name.Data())){
      // h.reset((TH1D*)input_files[info["fname"].Data()]->Get(info["appliesTo"]));
      h.reset(new TH1D(name,name, 10000,10,5000));
      double scale = stod(hname.ReplaceAll("scale","").Data());
      hname = name;
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
        h->SetBinContent(i, scale);
        h->SetBinError(i, scale);
      }
    } else {
      h.reset((TH1D*)input_files[info["fname"].Data()]->Get(hname));
    }
    // sources[name] = new SystematicContainer(name,appliesTo, atoi(info["index"]), hname, h.get());
    sources[name] = new SystematicContainer(name,appliesTo, sources.size() , hname, h.get());
    assert(sources[name]);
    if (debug) {PrintLine("successfully", yellow); cout << *sources[name] << endl;};
  }
}

void GlobalFit::LoadShapes(){
  for (auto [name, info]: shapes_map) {
    TString appliesTo = info["appliesTo"];
    if (FindInVector(types,appliesTo)<0) {
      PrintLine("Skipping source: "+name+" index: "+to_string(FindInVector(types,appliesTo)), yellow);
      continue;
    }
    if (debug) PrintLoading("shape", name, "");
    // shapes[name] = new ShapeContainer(name,info["form"],info["appliesTo"],atoi(info["index"]),atoi(info["ispositive"]));
    int index = shapes.size();
    if (info["appliesTo"]!="Resp") index = shapes[info["type"]]->index();
    shapes[name] = new ShapeContainer(name,info["form"],info["appliesTo"],index,atoi(info["ispositive"]));
    assert(shapes[name]);
    if (debug) {PrintLine("successfully", yellow); cout << *shapes[name] << endl;};
    shape_types.insert(info["type"]);
  }
}

void GlobalFit::LoadReference(){
  for (auto [name, info]: reference_obj_map) {
    TString hname = ReplaceDefault(info["hname"]);
    if (debug) PrintLoading("hist", name, hname);
    reference_objects[name] =(TH1D*)input_files[info["fname"].Data()]->Get(hname)->Clone(name);
    assert(reference_objects[name]);
    reference_objects[name]->SetDirectory(0);
    if (debug) {PrintLine("successfully", yellow);}
  }
}

void GlobalFit::LoadFSR(){
  for (auto [name, info]: kfsr_hnames_map) {
    TString hname = ReplaceDefault(info["hname"]);
    if (debug) PrintLoading("kfsr", name, hname);
    unique_ptr<TH1D> h; h.reset((TH1D*)input_files[info["fname"].Data()]->Get(hname));
    fsrs[name] = new SystematicContainer(name,info["appliesTo"], -1, hname, h.get());
    assert(fsrs[name]);
    if (debug) {PrintLine("successfully", yellow); cout << *fsrs[name] << endl;};
  }
}

void GlobalFit::ScaleFSR(){
  for (auto [name, fsr]: fsrs) {
    auto dt = my_data[fsr->appliesTo()];
    TGraphErrors* raw = dt->raw();
    TGraphErrors* input = dt->input();
    TGraphErrors* output = dt->output();
    TGraphErrors* variation = dt->variation();
    for (int i = 0; i != raw->GetN(); ++i) {
      double pt = raw->GetX()[i];
      double r = raw->GetY()[i];
      double shift = fsr->hist()->GetBinContent(fsr->hist()->FindBin(pt));
      input->SetPoint(i, pt, r+shift);
      output->SetPoint(i, pt, r+shift);
      variation->SetPoint(i, pt, r+shift);
    }
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
    fitter->SetParameter(i, "", 1, (i<nFitPars ? 0.01 : 1),-100, 100);
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
  for (auto [name,dt]: my_data){
    TGraphErrors* raw = dt->raw();
    TGraphErrors* input = dt->input();
    TGraphErrors* output = dt->output();
    TGraphErrors* variation = dt->variation();
    double scale = (dt->type()=="Resp") ? 1.: 1./ScaleFullSimShape; // transform PF comp in percentage
    multiplyGraph(input,scale);
    multiplyGraph(output,scale);
    multiplyGraph(variation,scale);
    std::vector<TF1*> funcs;
    for (auto [name,shape]: shapes){
      if (shape->appliesTo()==dt->type()) funcs.push_back(shape->func());
    }
    PropagateErrorToGraph(variation, funcs, *error_matrix.get());
    raw->Write(name+"_raw", TObject::kOverwrite);
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
    prefit->SetParameter(0,1);
    postfit->SetParameter(0,_jesFit->GetParameter(shape->index()));
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
  LoadReference();
  LoadFSR();
  ScaleFSR();
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
        for (auto [src_name_,source]: sources){
          if (dt_name!=source->appliesTo()) continue;
          if (debug) PrintLine("adding1 "+src_name_+" to shift for "+ dt_name+"("+source->appliesTo()+")", yellow);
          TH1D* hsrc = source->hist(); assert(hsrc);
          shift += nuisance_pars[source->index()] * hsrc->GetBinContent(hsrc->FindBin(pt));
        }

        for (auto [shape_name,shape]: shapes){
          if ("Resp"==shape->appliesTo()) continue;
          if (dt_type!=shape->appliesTo()) continue;
          if (pt<40 || pt>600) continue;
          if (debug) PrintLine("adding2 "+shape_name+" to shift for "+ dt_type+"("+shape->appliesTo()+") index:"+to_string(shape->index()), yellow);

          TF1* func = shape->func(); assert(func);
          shift += fit_pars[shape->index()] * func->Eval(pt) * GlobalFit::ScaleFullSimShape;
          cout << "  --> index:" << shape->index() << " " << fit_pars[shape->index()] << " " << func->Eval(pt) << " " << shift << endl;
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
