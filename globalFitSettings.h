#ifndef globalFitSettings
#define globalFitSettings
#include <array>
#include <string>
#include <map>

using std::array;
using std::string;
using std::map;

// Global minimum (stat.) uncertainty for input data
double globalErrMin = 0.003;//0.002559; // chi2/NDF=101.0/101

// Input data
struct fitData {
  string name;
  string type;
  TGraphErrors *input;
  TGraphErrors *input2;
  TGraphErrors *output;
  TGraphErrors *output2;
};

// Systematics for input data
struct fitSyst {
  int idx;
  string name;
  string appliesTo;
  TH1D *hist;
};

// Shapes possible for JES and PF composition
struct fitShape {
  int idx;
  string name;
  string appliesTo;
  bool ispos;
  TF1 *func;
};

// Gaussian prior for fit parameters
bool _gf_penalizeFitPars = true;

// Downweight PF compositin data to get broad features
double _gf_fitPFcomp_minErr = 0.001;
// Do PFcomp, Z+jet HDM fit (will still plot these, if enabled)
bool _gf_fitPFcomp = false;
bool _gf_fitZjetHDM = true;//false;
// Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
bool _gf_useJESref = false;//true;//false;
// Alternatively, undo JES completely to show pre-closure plots
bool _gf_undoJESref = false;//true;

// Listing of all available 'datasets'
// {name, type, name2}
// 'name' must match graph name in rootfiles/jecdata[X].root:ratio/eta_[Y]/
// 'type' decides which of 'shapes' below is applied, e.g. 'Rjet' or 'chf'
// 'name2' is secondary input, so far only used for 'multijet' modes
// How to use: comment out datasets that are not needed to switch them off
const unsigned int ndt = 44;
const array<array<string,3>,ndt> _gf_datasets = {{
    {"ptchs_zjet_a100", "Rjet",""},
    {"mpfchs1_zjet_a100", "Rjet",""},
    {"hdm_mpfchs1_zjet", "Rjet",""},
    {"chf_zjet_a100",  "chf",""},
    {"nhf_zjet_a100",  "nhf",""},
    {"nef_zjet_a100",  "nef",""},
    {"ptchs_gamjet_a100", "Rjet",""},
    {"mpfchs1_gamjet_a100", "Rjet",""},
    {"hdm_mpfchs1_gamjet", "Rjet",""},
    {"chf_gamjet_a100",  "chf",""},
    {"nhf_gamjet_a100",  "nhf",""},
    {"nef_gamjet_a100",  "nef",""},
    {"ptchs_multijet_a100", "Rjet","crecoil_multijet_a100"},
    {"mpfchs1_multijet_a100", "Rjet","crecoil_multijet_a100"},
    {"hdm_mpfchs1_multijet", "Rjet","crecoil_multijet_a100"},
    {"chf_multijet_a100",  "chf",""},
    {"nhf_multijet_a100",  "nhf",""},
    {"nef_multijet_a100",  "nef",""},
    {"chf_incjet_a100",  "chf",""},
    {"nhf_incjet_a100",  "nhf",""},
    {"nef_incjet_a100",  "nef",""},
  }};

// Use only these data sets (empty for all)
// {name}
// 'name' should match the dataset name in the list above
// How to use: uncomment individual datasets to use only those
const array<string,29> _gf_datasets_whitelist = {
  //"ptchs_zjet_a100",
  //"mpfchs1_zjet_a100",
  "hdm_mpfchs1_zjet",
  "chf_zjet_a100",
  "nhf_zjet_a100",
  "nef_zjet_a100",
  //"ptchs_gamjet_a100",
  //"mpfchs1_gamjet_a100",
  "hdm_mpfchs1_gamjet",
  "chf_gamjet_a100",
  "nhf_gamjet_a100",
  "nef_gamjet_a100",
  //"ptchs_multijet_a100",
  //"mpfchs1_multijet_a100",
  "hdm_mpfchs1_multijet",
  //"chf_multijet_a100",
  //"nhf_multijet_a100",
  //"nef_multijet_a100",
  "chf_incjet_a100",
  "nhf_incjet_a100",
  "nef_incjet_a100",
};

// Use only these shapes (empty for all)
// {name}
// 'name' should match the dataset name in the list above
// How to use: uncomment individual datasets to use only those
const array<string,29> _gf_shapes_whitelist = {
  "em3",
  "tv402","tv3n1","tv300pn",
  "hhp3",
  "hhpfc", // Run23
  "hhnoise", // Run22
  //"x0p5",
  "x1",
  //"x1p5",
  //"x2",
  //"x3",
  //"x4", // test TeV scale
  "hbtime"
};

// Listing one-sided (positive-definite) sources
const int npos = 19;
const array<string,npos> _gf_posdef =
  {"tv4","tv402","tv404","tv405","tv410","tv416","tv420","tv430",
   "tv3n1","tv300pn",
   //"hhp3",
   //"hhpfc",
   "hhred103","hhred100","hhred097", "hhblue103","hhblue100","hhblue097",
   /*"x2",*/ "hbtime"
  };

// Listing source limits
// (How to best implement this?)

// Listing of all available uncertainty 'sources'
// {name, appliesTo, histname}
// 'name' can repeat, e.g. 'FSRuncl' for multiple inputs, and uses same nuisance
// 'appliesTo' must match the dataset name above
// 'histname' must match histogram name in rootfiles/jecdata[X].root:[Y]/sys/
// How to use: list all 'sources'. Ones not matching active inputs are ignored
// NB: For one-sided (positive) sources, add them to _gf_posdef list
const int nsrc = 13;
const array<array<string,3>,nsrc> _gf_sources = {{
    //{"FSRunclDTMJ","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu1"},
    //{"FSRunclMCMJ","hdm_mpfchs1_multijet","hkfsr3_mpfchs1_multijet_mpfu2"},
    //{"Runcl1","hdm_mpfchs1_zjet","hkfsr3_mpfchs1_zjet_mpfu1"},
    //{"Runcl2","hdm_mpfchs1_zjet","hkfsr3_mpfchs1_zjet_mpfu2"},
  }};


// Listing of all available JEC fit 'shapes'
// {name, appliesTo, funcstring}
// 'name' can repeat e.g. for 'Rjet' and 'chf', and uses same fit parameter
// 'appliesTo' must match one of the 'type' in dataset listing
// 'funcstring' is a TF1-style string for a fixed functions. Use 'x' for pT
// How this works: each unique 'name' is assigned a fit parameter that
// multiplies 'funcstring' when applied to dataset with same 'appliesTo' type
//
// Fits are produced with minitools/fullSimShapes.C with doProduction=true
// Copy them from pdf/fullSimShapes/txt and /prod to /handpicked
// Each plot/txt pair is identified by unique MD5 checksum calculated from
// concatenated Rjet,chf,nef,nhf functions strings
const int nshp = 112;
const array<array<string,3>, nshp> _gf_shapes = {{

    {"x0p5","Rjet","10.*pow(x/3000.,0.5)"},
    {"x1","Rjet","10.*pow(x/3000.,1)"},
    {"x1p5","Rjet","10.*pow(x/3000.,1.5)"},
    {"x2","Rjet","10.*pow(x/3000.,2)"},
    {"x3","Rjet","10.*pow(x/3000.,3)"},
    {"x4","Rjet","10.*pow(x/3000.,4)"},
    // Quick tests for 10% pT^2 and pT^3 uncertainty at 3 TeV

    {"hbtime","Rjet","0.1394+log(x)*(-0.08167+log(x)*(-0.001402+log(x)*(0.005871+log(x)*(0.0003952+log(x)*(-0.0003595+log(x)*2.291e-05)))))"},
    {"hbtime","chf","1.029+log(x)*(-0.9924+log(x)*(0.2604+log(x)*(0.03122+log(x)*(-0.02631+log(x)*(0.004353+log(x)*-0.0002299)))))"},
    {"hbtime","nef","-1.396+log(x)*(1.355+log(x)*(-0.3602+log(x)*(-0.04115+log(x)*(0.03584+log(x)*(-0.005978+log(x)*0.0003229)))))"},
    {"hbtime","nhf","0.5664+log(x)*(-0.5401+log(x)*(0.1418+log(x)*(0.01537+log(x)*(-0.01331+log(x)*(0.002179+log(x)*-0.0001194)))))"},
    // "hbtime" checksum: af96d36f092413556b2d757953e90fae
    /*
    // Add 2.* by hand (5ns -> 10ns equivalent)
    {"hbtime","Rjet","2.*log(x/15.)*(-0.001278+log(x/15.)*(0.01576+log(x/15.)*-0.01456))"},
    {"hbtime","chf","2.*log(x/15.)*(6.107e-05+log(x/15.)*(-0.004846+log(x/15.)*0.006862))"},
    {"hbtime","nef","2.*log(x/15.)*(-0.007979+log(x/15.)*(0.01026+log(x/15.)*-0.003737))"},
    {"hbtime","nhf","2.*log(x/15.)*(0.005884+log(x/15.)*(-0.001467+log(x/15.)*-0.004223))"},
    // "hbtime" checksum: cc2b0bb93a82fa64200db36d81c05b1e
    */
    
    {"em3","Rjet","-0.3222+log(x)*(0.1973+log(x)*(-0.0774+log(x)*(-0.006124+log(x)*(0.001035+log(x)*(0.0002058+log(x)*-2.089e-05)))))"},
    {"em3","chf","2.59+log(x)*(-1.188+log(x)*(0.1324+log(x)*(0.01962+log(x)*(-0.001522+log(x)*(-0.0004592+log(x)*3.848e-05)))))"},
    {"em3","nef","-2.715+log(x)*(1.225+log(x)*(-0.1397+log(x)*(-0.02023+log(x)*(0.001699+log(x)*(0.0004898+log(x)*-4.382e-05)))))"},
    {"em3","nhf","2.168+log(x)*(-1.761+log(x)*(0.4903+log(x)*(-0.0334+log(x)*(-0.007356+log(x)*(0.001323+log(x)*-5.843e-05)))))"},
    // "em3" checksum: c8cb4a54ecee28c98c7fe286c0beafa9


    {"tv402","Rjet","-3.512+log(x)*(3.975+log(x)*(-1.433+log(x)*(0.1443+log(x)*(0.01588+log(x)*(-0.003781+log(x)*0.0001845)))))"},
    {"tv402","chf","-14.72+log(x)*(12.93+log(x)*(-3.671+log(x)*(0.2582+log(x)*(0.04383+log(x)*(-0.007821+log(x)*0.0003428)))))"},
    {"tv402","nef","2.082+log(x)*(-0.8289+log(x)*(-0.2106+log(x)*(0.1027+log(x)*(0.00176+log(x)*(-0.002916+log(x)*0.0002057)))))"},
    {"tv402","nhf","14.14+log(x)*(-13.34+log(x)*(4.169+log(x)*(-0.3458+log(x)*(-0.06092+log(x)*(0.01287+log(x)*-0.0006435)))))"},
    // "tv402" checksum: 6b76d4d2922fa185b8fe54ae5a5a56a5

    {"tv3n1","Rjet","17.03+log(x)*(-14.94+log(x)*(4.142+log(x)*(-0.3059+log(x)*(-0.047+log(x)*(0.008868+log(x)*-0.0003968)))))"},
    {"tv3n1","chf","17.76+log(x)*(-16.17+log(x)*(4.619+log(x)*(-0.3464+log(x)*(-0.05633+log(x)*(0.01074+log(x)*-0.0004876)))))"},
    {"tv3n1","nef","-8.682+log(x)*(7.874+log(x)*(-2.232+log(x)*(0.164+log(x)*(0.02771+log(x)*(-0.00521+log(x)*0.0002358)))))"},
    {"tv3n1","nhf","-9.61+log(x)*(8.698+log(x)*(-2.455+log(x)*(0.1624+log(x)*(0.03711+log(x)*(-0.006591+log(x)*0.0002973)))))"},
    // "tv3n1" checksum: 348f764eaddd02005d2c66cabf4fc863

    {"tv300pn","Rjet","-0.8581+log(max(30.,x))*(0.5263+log(max(30.,x))*(-0.06701+log(max(30.,x))*(-0.008783+log(max(30.,x))*(0.0008284+log(max(30.,x))*(0.0002308+log(max(30.,x))*-2.113e-05)))))"},
    {"tv300pn","chf","6.831+log(max(30.,x))*(-6.417+log(max(30.,x))*(2.039+log(max(30.,x))*(-0.181+log(max(30.,x))*(-0.02781+log(max(30.,x))*(0.005889+log(max(30.,x))*-0.0002787)))))"},
    {"tv300pn","nef","-14.02+log(max(30.,x))*(12.55+log(max(30.,x))*(-3.805+log(max(30.,x))*(0.3161+log(max(30.,x))*(0.05375+log(max(30.,x))*(-0.01122+log(max(30.,x))*0.0005451)))))"},
    {"tv300pn","nhf","6.388+log(max(30.,x))*(-5.397+log(max(30.,x))*(1.524+log(max(30.,x))*(-0.1035+log(max(30.,x))*(-0.026+log(max(30.,x))*(0.00501+log(max(30.,x))*-0.0002474)))))"},
    // "tv300pn" checksum: a25af6468140b2a27a00c0fc1bab83df

    {"hhp3","Rjet","2.817+log(x)*(-2.491+log(x)*(0.698+log(x)*(-0.04045+log(x)*(-0.008571+log(x)*(0.001366+log(x)*-5.705e-05)))))"},
    {"hhp3","chf","1.363+log(x)*(-1.602+log(x)*(0.5867+log(x)*(-0.06816+log(x)*(-0.00732+log(x)*(0.002+log(x)*-0.0001039)))))"},
    {"hhp3","nef","-7.619+log(x)*(7.721+log(x)*(-2.651+log(x)*(0.2725+log(x)*(0.03317+log(x)*(-0.00842+log(x)*0.000441)))))"},
    {"hhp3","nhf","7.178+log(x)*(-6.88+log(x)*(2.246+log(x)*(-0.1991+log(x)*(-0.03391+log(x)*(0.007537+log(x)*-0.0003856)))))"},
    // "hhp3" checksum: 57569f81de14365a30d012311898ae2a

    {"hhpfc","Rjet","15.64+log(x)*(-11.27+log(x)*(2.735+log(x)*(-0.2815+log(x)*0.01065)))"},
    {"hhpfc","chf","-2.061+log(x)*(1.798+log(x)*(-0.3602+log(x)*0.02139))"},
    {"hhpfc","nef","-2.86+log(x)*(1.904+log(x)*(-0.3753+log(x)*0.02323))"},
    {"hhpfc","nhf","13.04+log(x)*(-10.17+log(x)*(2.611+log(x)*(-0.2802+log(x)*0.01085)))"},
    // "hhpfc" checksum: 2adb48123dcde1aa50a8bba8d138f7ca

    // RunGC
    /*
    {"hhnoise","Rjet","5.28+log(x)*(-6.859+log(x)*(2.04+log(x)*(-0.09325+log(x)*(-0.04397+log(x)*(0.006839+log(x)*-0.0002948)))))"},
    {"hhnoise","chf","-2.11+log(x)*(3.659+log(x)*(-1.207+log(x)*(0.07304+log(x)*(0.02414+log(x)*(-0.004134+log(x)*0.0001881)))))"},
    {"hhnoise","nef","-7.498+log(x)*(6.744+log(x)*(-1.875+log(x)*(0.1007+log(x)*(0.03691+log(x)*(-0.006191+log(x)*0.0002811)))))"},
    {"hhnoise","nhf","8.571+log(x)*(-9.248+log(x)*(2.631+log(x)*(-0.1182+log(x)*(-0.05626+log(x)*(0.008792+log(x)*-0.0003823)))))"},
    // "hhnoise" checksum: bdbd8c2eba49736fa0f627106e80d40e

    {"hhnoise3","Rjet","(5.28+log(x)*(-6.859+log(x)*(2.04+log(x)*(-0.09325+log(x)*(-0.04397+log(x)*(0.006839+log(x)*-0.0002948)))))-1)*3+1"},
    {"hhnoise3","chf","(-2.11+log(x)*(3.659+log(x)*(-1.207+log(x)*(0.07304+log(x)*(0.02414+log(x)*(-0.004134+log(x)*0.0001881))))))*3"},
    {"hhnoise3","nef","(-7.498+log(x)*(6.744+log(x)*(-1.875+log(x)*(0.1007+log(x)*(0.03691+log(x)*(-0.006191+log(x)*0.0002811))))))*3"},
    {"hhnoise3","nhf","(8.571+log(x)*(-9.248+log(x)*(2.631+log(x)*(-0.1182+log(x)*(-0.05626+log(x)*(0.008792+log(x)*-0.0003823))))))*3"},
    // "hhnoise" checksum: bdbd8c2eba49736fa0f627106e80d40e, hand modified
    */
    // RunGfix
    {"hhnoise","Rjet","-2.822+log(x)*(-0.05219+log(x)*(0.1295+log(x)*(0.003666+log(x)*(-0.001746+log(x)*(-0.000216+log(x)*2.772e-05)))))"},
    {"hhnoise","chf","-1.805+log(x)*(2.518+log(x)*(-0.8305+log(x)*(0.08099+log(x)*(0.006106+log(x)*(-0.001582+log(x)*7.705e-05)))))"},
    {"hhnoise","nef","-3.583+log(x)*(2.995+log(x)*(-0.7649+log(x)*(0.03897+log(x)*(0.01168+log(x)*(-0.001744+log(x)*6.837e-05)))))"},
    {"hhnoise","nhf","9.528+log(x)*(-9.178+log(x)*(2.683+log(x)*(-0.2096+log(x)*(-0.03239+log(x)*(0.006442+log(x)*-0.0003013)))))"},
    // "hhnoise" checksum: a6b125bd41febeacfb42c782eaa93449

    /*
    {"hhp18","Rjet","((2.817+log(x)*(-2.491+log(x)*(0.698+log(x)*(-0.04045+log(x)*(-0.008571+log(x)*(0.001366+log(x)*-5.705e-05))))))-1)*6+1"},
    {"hhp18","chf","(1.363+log(x)*(-1.602+log(x)*(0.5867+log(x)*(-0.06816+log(x)*(-0.00732+log(x)*(0.002+log(x)*-0.0001039))))))*6"},
    {"hhp18","nef","(-7.619+log(x)*(7.721+log(x)*(-2.651+log(x)*(0.2725+log(x)*(0.03317+log(x)*(-0.00842+log(x)*0.000441))))))*6"},
    {"hhp18","nhf","(7.178+log(x)*(-6.88+log(x)*(2.246+log(x)*(-0.1991+log(x)*(-0.03391+log(x)*(0.007537+log(x)*-0.0003856))))))*6"},
    // custom variant
    */
    /*
    {"hhblue103","Rjet","-52.44+log(max(55.,x))*(36.13+log(max(55.,x))*(-7.941+log(max(55.,x))*(0.2816+log(max(55.,x))*(0.116+log(max(55.,x))*(-0.01448+log(max(55.,x))*0.0004894)))))"},
    {"hhblue103","chf","6.354+log(max(55.,x))*(-3.394+log(max(55.,x))*(0.3485+log(max(55.,x))*(0.06261+log(max(55.,x))*(-0.004791+log(max(55.,x))*(-0.001639+log(max(55.,x))*0.0001499)))))"},
    {"hhblue103","nef","-0.5177+log(max(55.,x))*(0.4554+log(max(55.,x))*(-0.07566+log(max(55.,x))*(-0.008602+log(max(55.,x))*(0.001275+log(max(55.,x))*(0.0002701+log(max(55.,x))*-3.214e-05)))))"},
    {"hhblue103","nhf","-5.249+log(max(55.,x))*(2.559+log(max(55.,x))*(-0.1941+log(max(55.,x))*(-0.05668+log(max(55.,x))*(0.002404+log(max(55.,x))*(0.001512+log(max(55.,x))*-0.0001229)))))"},
    // "hhblue103" checksum: b4332baa28b900ace3c09361eee1850e
    */
  }};

#endif
