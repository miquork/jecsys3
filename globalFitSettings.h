#ifndef globalFitSettings
#define globalFitSettings
#include <array>
#include <string>
#include <map>

using std::array;
using std::string;
using std::map;

// Global minimum (stat.) uncertainty for input data
double globalErrMin = 0.005;//0.003;//0.002559; // chi2/NDF=101.0/101

//double scaleJZ = 0.98;  // default 0.98 for 2022-23
//double scaleJZA = 0.99;  // default for 2022-23
bool scaleJZperEra = true; // Enable custom scales per era (hard-coded)
bool scaleJZAperEra = true; // Enable custom scales per era (hard-coded)
double scaleJZ = 0.985;  // default 0.985 for 2024CS
double scaleJZA = 0.9925;  // default 0.9925 for 2024CS

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
bool _gf_fitPFcomp = false;//true;//false;
bool _gf_fitZjetHDM = true;//false;
// Permit baseline shift for JES with L2L3Res inputs. Needed for PFcomp fit
bool _gf_useJESref = false;//true;//false;
// Alternatively, undo JES completely to show pre-closure plots
bool _gf_undoJESref = true;//false;//true; // Set 'false' to produce "*_closure.pdf"

// Listing of all available 'datasets'
// {name, type, name2}
// 'name' must match graph name in rootfiles/jecdata[X].root:ratio/eta_[Y]/
// 'type' decides which of 'shapes' below is applied, e.g. 'Rjet' or 'chf'
// 'name2' is secondary input, so far only used for 'multijet' modes
// How to use: comment out datasets that are not needed to switch them off
const unsigned int ndt = 44;
const array<array<string,3>,ndt> _gf_datasets = {{
    {"mpfchs1_wqq_a100", "Rjet",""},
    {"ptchs_jetz_a100", "Rjet",""},
    {"mpfchs1_jetz_a100", "Rjet",""},
    {"hdm_mpfchs1_jetz", "Rjet",""},
    {"ptchs_zjav_a100", "Rjet",""},
    {"mpfchs1_zjav_a100", "Rjet",""},
    {"hdm_mpfchs1_zjav", "Rjet",""},
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
  "mpfchs1_wqq_a100",
  //"ptchs_zjet_a100",
  //"mpfchs1_zjet_a100",
  //"hdm_mpfchs1_jetz", // V8M
  //"hdm_mpfchs1_zjav", // V8M
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
  //"ptchs_multijet_a100", // use ptchs before v34d bug fixed
  //"mpfchs1_multijet_a100",
  "hdm_mpfchs1_multijet", // mpfchs->hdm botched before v34d
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
  "em3", // EM scale -3%
  "tv402","tv3n1","tv300pn", // tracking -1% Ntrk=1, Ntrk>=1; eff=0.998^{N-1}
  "hhp3", // HCAL scale -3%
  //"hhpfc", // Run23
  //"hhnoise", // Run22
  //"x0p5",
  //"x1", // 22Sep2023 V3, not 19Dec2023
  //"g1", // NEF drop above 610 GeV
  //"p5log", // log^5
  //"p4log", // log^4
  //"p3log", // log^3
  //"quadlog", // log^2
  //"loglin", // log^1
  //"const", // log^0 constant scale factor
  //"off", // Summer23
  //"off_ho", // V9M test
  "off_nhf", // Prompt24, no hhpfc+hhnoise+hbtime+hbsipm+off with this
  "dd_nhf", // Prompt25, model HCAL problems in 2025C
  //"ecalcc", // ECAL cc timing for 2024BCD, not 2024CR/CS/F
  //"qie11"
  //"x1p5",
  //"x1v4", // Summer23 variant
  //"x1p5v4", // Summer23 high pT? just tad steeper than x1
  //"x2v4", // Summer23 high pT? bit more conservative extrapolation?
  //"x2p5", // Summer23 high pT? pretty good
  //"x3", // Summer23 hight pT? Slightly too sharp?
  //"x4", // test TeV scale
  //"hbtime", // test adding for Summer23; did nothing really on top of x2p5
  //"hbsipm", // 22Sep2023 V3, not 19Dec2023, off for 2025
  "hhblue103",
  "hbd1","hbd2","hbd3","hbd4" // Depth variations by 20% (or 17%?) for 2025
};

// Listing one-sided (positive-definite) sources
const int npos = 19;
const array<string,npos> _gf_posdef =
  {"tv4","tv402","tv404","tv405","tv410","tv416","tv420","tv430",
   "tv3n1","tv300pn",
   //"hhp3",
   //"hhpfc",
   /*"hhred103","hhred100","hhred097", "hhblue103","hhblue100","hhblue097",*/
   /*"x2",*/ /*"hbtime"*/ /*, "hbsipm"*/
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

    // 400e3(max)/1300.(per GeV per channel)*4.(channels)/0.5(HCAL fraction)
    // 1500e3(max)/1300.(per GeV per channel)/0.5(HCAL fraction)
    //{"qie11","Rjet","-0.5*max(x-1150.,0.)/(3000.-1150.)"},
    // R = (0.5*x+1150.)/x
    //{"qie11","Rjet","-0.5*max(x-1150.,0.)/(3000.-1150.)"},
    {"off","Rjet","1./x"},
    {"off_ho","Rjet","5./(x/15.)"},
    //{"off_nhf","Rjet","18.46*pow(x/15.6,1.840)/(1+0.5*pow(x/15.6,3.679))"},
    //{"off_nhf","Rjet","19.15*pow(x/15.6,1.912)/(1+0.5*pow(x/15.6,3.824))"},
  //{"off_nhf","Rjet","19.26*pow(x/15.6,1.927)/(1+0.5*pow(x/15.6,3.853))"},//V8M
    {"off_nhf","Rjet","15.49*pow(x/15.1,1.939)/(1+0.5*pow(x/15.1,3.877))"},//V9M
    //{"dd_nhf","Rjet","374.11/x-83.36*pow(x,-0.5104)"}, // Prompt2025, v1
    {"dd_nhf","Rjet","370.96/x-83.46*pow(x,-0.5103)"},
    {"g1","Rjet","-6.5*max(log(x/610.),0.)"},
    {"p5log","Rjet","1.0*pow(log(x/300.),5)"},
    {"p4log","Rjet","1.0*pow(log(x/300.),4)"},
    {"p3log","Rjet","1.0*pow(log(x/300.),3)"},
    {"quadlog","Rjet","1.0*pow(log(x/300.),2)"},
    {"loglin","Rjet","1.0*log(x/300.)"},
    {"const","Rjet","1.0"},
    {"x0p5","Rjet","10.*pow(x/3000.,0.5)"},
    {"x1","Rjet","10.*pow(x/3000.,1)"},
    {"x1v4","Rjet","10.*pow(x/4000.,1)"},
    {"x1p5","Rjet","10.*pow(x/3000.,1.5)"},
    {"x1p5v4","Rjet","10.*pow(x/4000.,1.5)"},
    {"x2","Rjet","10.*pow(x/3000.,2)"},
    {"x2v4","Rjet","10.*pow(x/4000.,2)"},
    {"x2p5","Rjet","10.*pow(x/3000.,2.5)"},
    {"x2p5cap","Rjet","10.*pow(min(x/3000.,1.),2.5)"},
    {"x3","Rjet","10.*pow(x/3000.,3)"},
    {"x4","Rjet","10.*pow(x/3000.,4)"},
    // Quick tests for 10% pT^2 and pT^3 uncertainty at 3 TeV

    // HB depth variations for depths 1,2,3,4
    {"hbd1","Rjet","-3.295+log(x)*(2.364+log(x)*(-0.6283+log(x)*(0.03461+log(x)*(0.009515+log(x)*(-0.001753+log(x)*9.142e-05)))))"},
    {"hbd1","chf","0.6206+log(x)*(-0.1393+log(x)*(-0.001886+log(x)*(0.003798+log(x)*(0.0006411+log(x)*(1.277e-05+log(x)*-1.549e-05)))))"},
    {"hbd1","nef","7.78+log(x)*(-7.391+log(x)*(2.32+log(x)*(-0.1615+log(x)*(-0.04895+log(x)*(0.009229+log(x)*-0.0004418)))))"},
    {"hbd1","nhf","-10.4+log(x)*(9.327+log(x)*(-2.848+log(x)*(0.1941+log(x)*(0.05795+log(x)*(-0.01106+log(x)*0.0005452)))))"},
    // "hbd1" checksum: 605291bc9220bea47529cf5500ceb299

    {"hbd2","Rjet","-1.473+log(x)*(0.7143+log(x)*(-0.07459+log(x)*(-0.018+log(x)*(-0.0001996+log(x)*(0.0002706+log(x)*-7.195e-06)))))"},
    {"hbd2","chf","-4.981+log(x)*(4.784+log(x)*(-1.447+log(x)*(0.09612+log(x)*(0.02997+log(x)*(-0.005195+log(x)*0.0002244)))))"},
    {"hbd2","nef","13.8+log(x)*(-12.99+log(x)*(3.973+log(x)*(-0.2471+log(x)*(-0.08869+log(x)*(0.0158+log(x)*-0.000728)))))"},
    {"hbd2","nhf","-14.01+log(x)*(12.73+log(x)*(-3.826+log(x)*(0.2391+log(x)*(0.0809+log(x)*(-0.01473+log(x)*0.0007006)))))"},
    // "hbd2" checksum: 0d1494f63e15b067458feb1bd1d23f01

    {"hbd3","Rjet","-3.36+log(x)*(2.678+log(x)*(-0.6936+log(x)*(0.02975+log(x)*(0.0134+log(x)*(-0.002282+log(x)*0.0001095)))))"},
    {"hbd3","chf","-1.347+log(x)*(1.321+log(x)*(-0.4023+log(x)*(0.02565+log(x)*(0.008512+log(x)*(-0.001317+log(x)*4.75e-05)))))"},
    {"hbd3","nef","5.869+log(x)*(-5.448+log(x)*(1.624+log(x)*(-0.08804+log(x)*(-0.03836+log(x)*(0.006432+log(x)*-0.0002824)))))"},
    {"hbd3","nhf","-8.257+log(x)*(7.358+log(x)*(-2.136+log(x)*(0.1191+log(x)*(0.0468+log(x)*(-0.008154+log(x)*0.000378)))))"},
    // "hbd3" checksum: c982fae4d80f57c2e6e2a5477bedab01

    {"hbd4","Rjet","-2.433+log(x)*(1.973+log(x)*(-0.5149+log(x)*(0.02303+log(x)*(0.01015+log(x)*(-0.001673+log(x)*7.395e-05)))))"},
    {"hbd4","chf","1.292+log(x)*(-1.099+log(x)*(0.3116+log(x)*(-0.01734+log(x)*(-0.007365+log(x)*(0.001474+log(x)*-8.072e-05)))))"},
    {"hbd4","nef","1.916+log(x)*(-1.709+log(x)*(0.4718+log(x)*(-0.01473+log(x)*(-0.0123+log(x)*(0.001667+log(x)*-5.559e-05)))))"},
    {"hbd4","nhf","-3.991+log(x)*(3.467+log(x)*(-0.9667+log(x)*(0.04508+log(x)*(0.02209+log(x)*(-0.003581+log(x)*0.0001556)))))"},
    // "hbd4" checksum: b0013414f696dc10245acf35c3385aea
    
    // ECAL cc timing, dijet/drawCompareLite_Promp24_vs_ECALRATIO_TnP_2024B.pdf
    // Maximum drop put in by hand, now -9%:
    // 1) photon xsec max drop 75%: 18.75% scale if pT^-4,
    // 2) this impacts at most all ECAL energy or 60% of jet pT => 11.25% (10%)
    //{"ecalcc","Rjet","max(-1.59-6.68*(-0.798-0.5798*pow(x/396.1,1.412)/(1+pow(x/396.1,1.412))*(1-pow(x/396.1,-1.412)))-31.2*x/3000.,-11.25)"}, // 2024B
    {"ecalcc","Rjet","max(-0.77-4.58*(-0.798-0.5798*pow(x/396.1,1.412)/(1+pow(x/396.1,1.412))*(1-pow(x/396.1,-1.412)))-22*x/3000.,-10.)"}, // 2023C
    
    // Run3 wrong SiPM non-linearity corrections for data (1M variant)
    {"hbsipm","Rjet","-17.81+log(x)*(15.9+log(x)*(-4.719+log(x)*(0.3464+log(x)*(0.08054+log(x)*(-0.01553+log(x)*0.0007183)))))"},
    {"hbsipm","chf","7.678+log(x)*(-7.629+log(x)*(2.519+log(x)*(-0.2166+log(x)*(-0.04713+log(x)*(0.0104+log(x)*-0.0005435)))))"},
    {"hbsipm","nef","-1.717+log(x)*(2.522+log(x)*(-1.163+log(x)*(0.1623+log(x)*(0.01883+log(x)*(-0.006616+log(x)*0.0004364)))))"},
    {"hbsipm","nhf","-1.346+log(x)*(0.9761+log(x)*(-0.1417+log(x)*(-0.02259+log(x)*(0.003066+log(x)*(0.0009335+log(x)*-0.0001261)))))"},
    // "hbsipm" checksum: 36dd66b29cf0d8aec0b2a4cb145728ea
    
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

    // Add back for Prompt24
    {"hhblue103","Rjet","-52.44+log(max(55.,x))*(36.13+log(max(55.,x))*(-7.941+log(max(55.,x))*(0.2816+log(max(55.,x))*(0.116+log(max(55.,x))*(-0.01448+log(max(55.,x))*0.0004894)))))"},
    {"hhblue103","chf","6.354+log(max(55.,x))*(-3.394+log(max(55.,x))*(0.3485+log(max(55.,x))*(0.06261+log(max(55.,x))*(-0.004791+log(max(55.,x))*(-0.001639+log(max(55.,x))*0.0001499)))))"},
    {"hhblue103","nef","-0.5177+log(max(55.,x))*(0.4554+log(max(55.,x))*(-0.07566+log(max(55.,x))*(-0.008602+log(max(55.,x))*(0.001275+log(max(55.,x))*(0.0002701+log(max(55.,x))*-3.214e-05)))))"},
    {"hhblue103","nhf","-5.249+log(max(55.,x))*(2.559+log(max(55.,x))*(-0.1941+log(max(55.,x))*(-0.05668+log(max(55.,x))*(0.002404+log(max(55.,x))*(0.001512+log(max(55.,x))*-0.0001229)))))"},
    // "hhblue103" checksum: b4332baa28b900ace3c09361eee1850e

  }};

#endif
