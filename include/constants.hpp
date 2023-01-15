#pragma once
#include "utils.hpp"
#include "Containers.hpp"

static const bool debug=true;
static constexpr double func_range_min = 10.;  // Define fitting range
static constexpr double func_range_max = 6500.; // Define fitting range

static const std::vector<TString> all_types = {"Resp","chf","nef","nhf"};

static const std::map<TString, TString> input_fnames = {
  {"jes", "rootfiles/jecdataRUN.root"},
};

static const std::map<TString, std::map<TString,TString>> reference_obj_map = {
  {"hjesref",       { {"fname", "jes"}, {"type", "jes"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr"},}},

  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"hrun1",         { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hrun1"},}},
  {"hjes",          { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/hjes"},}},
  {"herr_l2l3res",  { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_l2l3res"},}},
  {"herr_ref",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_ref"},}},
  {"herr_noflv",    { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_noflv"},}},
  {"herr_spr",      { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_spr"},}},
  {"herr_pu",       { {"fname", "jes"}, {"type", "err"},  {"hname","MODE/etaETAMIN-ETAMAX/herr_pu"},}},
};

static const std::map<TString, std::map<TString,TString>> input_hnames_map = {
  {"Resp_zjet_mpf", { {"fname", "jes"}, {"type", "Resp"}, {"hname","MODE/etaETAMIN-ETAMAX/hdm_mpfchs1_zjet"},}},
  {"chf_zjet_mpf",  { {"fname", "jes"}, {"type", "chf"}, {"hname","MODE/etaETAMIN-ETAMAX/chf_zjet_a100"},}},
  {"nhf_zjet_mpf",  { {"fname", "jes"}, {"type", "nhf"}, {"hname","MODE/etaETAMIN-ETAMAX/nhf_zjet_a100"},}},
  {"nef_zjet_mpf",  { {"fname", "jes"}, {"type", "nef"}, {"hname","MODE/etaETAMIN-ETAMAX/nef_zjet_a100"},}},
};

static const std::map<TString, std::map<TString,TString>> sources_hnames_map = {
  {"Resp_uncl_zjet",    { {"index", "0"}, {"appliesTo", "Resp_zjet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfu2"},}},
  {"Resp_add_jet_zjet", { {"index", "1"}, {"appliesTo", "Resp_zjet_mpf"}, {"fname", "jes"}, {"hname","MODE/etaETAMIN-ETAMAX/fsr/hkfsr3_mpfchs1_zjet_mpfn2"},}},
};

static const std::map<TString, std::map<TString,TString>> shapes_map = {
  {"hhp3_Resp",    {{"index", "0"}, {"ispositive", "0"}, {"type","hhp3"},    {"appliesTo","Resp"}, {"form","2.817+log(x)*(-2.491+log(x)*(0.698+log(x)*(-0.04045+log(x)*(-0.008571+log(x)*(0.001366+log(x)*-5.705e-05)))))"},}},
  {"hhp3_chf",     {{"index", "0"}, {"ispositive", "1"}, {"type","hhp3"},    {"appliesTo","chf"},  {"form","1.363+log(x)*(-1.602+log(x)*(0.5867+log(x)*(-0.06816+log(x)*(-0.00732+log(x)*(0.002+log(x)*-0.0001039)))))"},}},
  {"hhp3_nef",     {{"index", "0"}, {"ispositive", "1"}, {"type","hhp3"},    {"appliesTo","nef"},  {"form","-7.619+log(x)*(7.721+log(x)*(-2.651+log(x)*(0.2725+log(x)*(0.03317+log(x)*(-0.00842+log(x)*0.000441)))))"},}},
  {"hhp3_nhf",     {{"index", "0"}, {"ispositive", "1"}, {"type","hhp3"},    {"appliesTo","nhf"},  {"form","7.178+log(x)*(-6.88+log(x)*(2.246+log(x)*(-0.1991+log(x)*(-0.03391+log(x)*(0.007537+log(x)*-0.0003856)))))"},}},
  // {"tv300pn_Resp", {{"index", "1"}, {"ispositive", "1"}, {"type","tv300pn"}, {"appliesTo","Resp"}, {"form","-0.8581+log(max(30.,x))*(0.5263+log(max(30.,x))*(-0.06701+log(max(30.,x))*(-0.008783+log(max(30.,x))*(0.0008284+log(max(30.,x))*(0.0002308+log(max(30.,x))*-2.113e-05)))))"},}},
};
