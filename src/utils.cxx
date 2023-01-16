#include "utils.hpp"


std::string NWhiteSpaces(int n) { return std::string(std::max(0,n), ' ' );};
TString CompleteWhiteSpaces(TString text, int ntot) { return text+NWhiteSpaces(ntot-text.Length());};

bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

void PrintLoading(TString type, TString name, TString objname, std::string color) {
  Print("Loading "+type+" for "+CompleteWhiteSpaces(name)+": "+objname+"..."+reset, color);
}

void PrintParameter(TString name, double val, double err, std::string color){
  PrintLine("  --> "+CompleteWhiteSpaces(name)+Form(": %+5.2f +/- %5.2f", val, err), color);
}

void RemoveZerosFromGraph(TGraphErrors *graph) {
  for (int i = graph->GetN()-1; i != -1; --i) {// Must go backwards to keep point ordering
    if (graph->GetY()[i]==0 && graph->GetEY()[i]==0) {
      graph->RemovePoint(i);
    }
  }
}

double oplus(double a, double b) {return sqrt(a*a + b*b);};

double fitError(TF1* func, TMatrixD err_matrix, double x, double k){
  int n_pars = func->GetNpar();
  // Partial derivatives as differentials with 10% step size
  std::vector<double> grad(n_pars);
  func->GradientPar(&x, &grad[0]);
  // Perform standard error propagation
  double sumerr2(0);
  for (int i = 0; i != n_pars; ++i) {
    for (int j = 0; j != n_pars; ++j) {
      sumerr2 += err_matrix[i][j]*grad[i]*grad[j];
    }
  }
  return (func->Eval(x) + k*sqrt(sumerr2));
}


void FuncToGraph(TF1* func, TMatrixD err_matrix, TGraphErrors *graph, double k){
  assert(func); assert(graph);
  std::unique_ptr<TGraph>gr; gr.reset(new TGraph(func));
  for (int i = 0; i < gr->GetN(); ++i) {
    double x = gr->GetX()[i];
    graph->AddPoint(x, gr->GetY()[i]);
    graph->SetPointError(i, 0., fitError(func, err_matrix, x,k) - gr->GetY()[i]);
  }
}

void FuncToHist(TF1* func, TMatrixD err_matrix, TH1D* hist, double k){
  assert(func); assert(hist);
  for (int j = 1; j != hist->GetNbinsX()+1; ++j) {
    double x = hist->GetBinCenter(j);
    hist->SetBinContent(j, x, func->Eval(x));
    hist->SetBinError(j, 0., fitError(func, err_matrix, x,k) - func->Eval(x));
  }
}


void multiplyGraph(TGraphErrors *graph, double scale) {
  for (int i = 0; i != graph->GetN(); ++i) {
    graph->SetPoint(i, graph->GetX()[i], scale*graph->GetY()[i]);
    graph->SetPointError(i, graph->GetEX()[i], scale*graph->GetEY()[i]);
  }
}


void multiplyGraph(TGraphErrors *graph, TF1 *func) {
  for (int i = 0; i != graph->GetN(); ++i) {
    double val = func->Eval(graph->GetX()[i]);
    graph->SetPoint(i, graph->GetX()[i], val*graph->GetY()[i]);
    graph->SetPointError(i, graph->GetEX()[i], val*graph->GetEY()[i]);
  }
}

void PropagateErrorToGraph(TGraphErrors *graph, std::vector<TF1*> funcs, TMatrixD err_matrix) {
  int nfuncs = funcs.size();
  // int n_pars = func->GetNpar();
  for (int i = 0; i != graph->GetN(); ++i) {
    double x = graph->GetX()[i];
    double xerr = graph->GetEX()[i];

    double sumerr2(0);
    for (int ipar = 0; ipar != nfuncs; ++ipar) {
      for (int jpar = 0; jpar != nfuncs; ++jpar) {
        sumerr2 += funcs[ipar]->Eval(x)*funcs[jpar]->Eval(x)*err_matrix[ipar][jpar];
      }
    }
    graph->SetPointError(i, xerr, sqrt(sumerr2));
  }
}
