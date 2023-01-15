#pragma once
// Purpose: Container objects for GlobalFit inputs
// Author: Andrea Malara
#include "TF1.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "utils.hpp"

class DataContainer {
public:
  DataContainer(TString name_, TString type_, TString hname_, TGraphErrors* graph_);
  void set_name(TString x) {name_ = x;}
  void set_type(TString x) {type_ = x;}
  void set_hname(TString x) {hname_ = x;}
  void set_input(TGraphErrors* x) {input_.reset(x);}
  void set_output(TGraphErrors* x) {output_.reset(x);}
  void set_variation(TGraphErrors* x) {variation_.reset(x);}
  std::string color() const {return color_;}
  TString name() const {return name_;}
  TString type() const {return type_;}
  TString hname() const {return hname_;}
  TGraphErrors* input() const {return input_.get();}
  TGraphErrors* output() const {return output_.get();}
  TGraphErrors* variation() const {return variation_.get();}
  int GetN() { assert(input()->GetN()==output()->GetN()); return input()->GetN();}

  friend std::ostream& operator<<(std::ostream& os, const DataContainer&);

private:
  std::string color_=blue;
  TString name_, type_, hname_;
  std::unique_ptr<TGraphErrors> input_, output_, variation_;
};


class ShapeContainer {

public:
  ShapeContainer(TString name_, TString form_, TString appliesTo_, int index_, bool ispositive_);
  void set_name(TString x) {name_ = x;}
  void set_form(TString x) {form_ = x;}
  void set_appliesTo(TString x) {appliesTo_ = x;}
  void set_index(int x) {index_ = x;}
  void set_ispositive(bool x) {ispositive_ = x;}
  void set_func(TF1* x) {func_.reset(x);}
  std::string color() const {return color_;}
  TString name() const {return name_;}
  TString form() const {return form_;}
  TString appliesTo() const {return appliesTo_;}
  int index() const {return index_;}
  bool ispositive() const {return ispositive_;}
  TF1* func() const {return func_.get();}

  friend std::ostream& operator<<(std::ostream& os, const ShapeContainer&);

private:
  std::string color_=magenta;
  TString name_, form_, appliesTo_;
  int index_;
  bool ispositive_;
  std::unique_ptr<TF1> func_;
};


class SystematicContainer {

public:
  SystematicContainer(TString name_, TString appliesTo_, int index_, TString hname_, TH1D* hist_);
  void set_name(TString x) {name_ = x;}
  void set_appliesTo(TString x) {appliesTo_ = x;}
  void set_hname(TString x) {hname_ = x;}
  void set_index(int x) {index_ = x;}
  void set_hist(TH1D* x) {hist_.reset(x); hist_->SetDirectory(0);}
  std::string color() const {return color_;}
  TString name() const {return name_;}
  TString appliesTo() const {return appliesTo_;}
  TString hname() const {return hname_;}
  int index() const {return index_;}
  TH1D* hist() const {return hist_.get();}

  friend std::ostream& operator<<(std::ostream& os, const SystematicContainer&);

  int GetN() { return hist()->GetNbinsX();}

private:
  std::string color_=cyan;
  TString name_, appliesTo_, hname_;
  int index_;
  std::unique_ptr<TH1D> hist_;
};
