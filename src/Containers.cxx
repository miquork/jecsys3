#include "Containers.hpp"
#include "constants.hpp"

using namespace std;

DataContainer::DataContainer(TString name_, TString type_, TString hname_, TGraphErrors* graph_){
  set_name(name_);
  set_type(type_);
  set_hname(hname_);
  RemoveZerosFromGraph(graph_);
  set_raw((TGraphErrors*)graph_->Clone(name_+"_raw"));
  set_input((TGraphErrors*)graph_->Clone(name_+"_in"));
  set_output((TGraphErrors*)graph_->Clone(name_+"_out"));
  set_variation((TGraphErrors*)graph_->Clone(name_+"_variation"));
}

ostream& operator<<(ostream& os, const DataContainer& fitshape) {
  PrintLine("name: "+fitshape.name(), fitshape.color());
  PrintLine("  --> type: "+fitshape.type(), fitshape.color());
  Print("  --> hname: "+fitshape.hname()+reset, fitshape.color());
  return os;
}

ShapeContainer::ShapeContainer(TString name_, TString form_, TString appliesTo_, int index_, bool ispositive_){
  set_name(name_);
  set_form(form_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_ispositive(ispositive_);
  set_func(new TF1("f1_"+name_+"_"+appliesTo_, form_, func_range_min,func_range_max));
}

ostream& operator<<(ostream& os, const ShapeContainer& fitshape) {
  PrintLine("name: "+fitshape.name(), fitshape.color());
  PrintLine("  --> form: "+fitshape.form(), fitshape.color());
  PrintLine("  --> appliesTo: "+fitshape.appliesTo(), fitshape.color());
  PrintLine("  --> index: "+to_string(fitshape.index()), fitshape.color());
  Print("  --> ispositive: "+to_string(fitshape.ispositive())+reset, fitshape.color());
  return os;
}


SystematicContainer::SystematicContainer(TString name_, TString appliesTo_, int index_, TString hname_, TH1D* hist_){
  set_name(name_);
  set_appliesTo(appliesTo_);
  set_index(index_);
  set_hname(hname_);
  set_hist((TH1D*)(hist_)->Clone(name_));
}

ostream& operator<<(ostream& os, const SystematicContainer& fitshape) {
  PrintLine("name: "+fitshape.name(), fitshape.color());
  PrintLine("  --> appliesTo: "+fitshape.appliesTo(), fitshape.color());
  PrintLine("  --> hname: "+fitshape.hname(), fitshape.color());
  Print("  --> index: "+to_string(fitshape.index())+reset, fitshape.color());
  return os;
}
