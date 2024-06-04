#ifndef __GLOBALFITSTYLES__
#define __GLOBALFITSTYLES__

#include "globalFitSettings.h"

// Purpose: Define plotting styles for globalFitRun2.C
//          Keep these in a separate file not to clutter main code

map<string, int> _gf_color;
_gf_color["ptchs_zjet_a100"]    = kRed;//kPink+2;
_gf_color["mpfchs1_zjet_a100"]    = kRed;//kPink+2;
_gf_color["ptchs_zjav_a100"]    = kRed+2;//kPink+2;
_gf_color["mpfchs1_zjav_a100"]    = kRed+2;//kPink+2;
_gf_color["ptchs_gamjet_a100"]  = kBlue;//kCyan+2;
_gf_color["mpfchs1_gamjet_a100"]  = kBlue;//kCyan+2;
_gf_color["hdm_mpfchs1_gamjet"]   = kBlue;
_gf_color["hdm_mpfchs1_gamjet2"]   = kBlack; // For PF
_gf_color["hdm_gamjet"]           = kBlue;
_gf_color["hdm_mpfchs1_jetz"]     = kRed+1;//+1;//0
_gf_color["hdm_mpfchs1_zjav"]     = kRed+2;//+1;//0
_gf_color["hdm_mpfchs1_zjet"]     = kRed;//+1;//0
_gf_color["hdm_mpfchs1_zlljet"]   = kRed;//+1;
_gf_color["mpfchs1_hadw_a30"]     = kGreen+2;
_gf_color["hdm_hadw"]             = kGreen+2;
_gf_color["hdm_mpfchs1_multijet"] = kBlack;//kGray+2;//kBlack;
_gf_color["ptchs_multijet_a100"]  = kBlack;
_gf_color["mpfchs1_multijet_a100"]  = kBlack;
_gf_color["hdm_cmb"]              = kGray+1;
_gf_color["hdm_cmb_mj"]           = kBlack;
_gf_color["hdm_incjet"]           = kOrange+2;
_gf_color["chf_gamjet_a100"]      = kRed+2;
_gf_color["nhf_gamjet_a100"]      = kGreen+3;
_gf_color["nef_gamjet_a100"]      = kBlue+2;
_gf_color["cef_gamjet_a100"]      = kCyan+2;
_gf_color["muf_gamjet_a100"]      = kMagenta+2;
_gf_color["chf_pfjet_a30"]      = kRed;
_gf_color["nhf_pfjet_a30"]      = kGreen+1;
_gf_color["nef_pfjet_a30"]      = kBlue;
_gf_color["cef_pfjet_a30"]      = kCyan;
_gf_color["muf_pfjet_a30"]      = kMagenta;
_gf_color["chf_zjet_a100"]      = kRed;//+3;//2;
_gf_color["nhf_zjet_a100"]      = kGreen+2;//4;//3;
_gf_color["nef_zjet_a100"]      = kBlue;//+3;//2;
_gf_color["cef_zjet_a100"]      = kCyan+1;//+3;//2;
_gf_color["muf_zjet_a100"]      = kMagenta+1;//+3;//2;
_gf_color["chf_zlljet_a100"]      = kRed+2;//+3;
_gf_color["nhf_zlljet_a100"]      = kGreen+3;//+4;
_gf_color["nef_zlljet_a100"]      = kBlue+2;//+3;
_gf_color["cef_zlljet_a100"]      = kCyan+2;//+3;
_gf_color["muf_zlljet_a100"]      = kMagenta+2;//3;
_gf_color["cef_zmmjet_a100"]      = kCyan+4;
_gf_color["muf_zmmjet_a100"]      = kMagenta+4;
//
_gf_color["chf_multijet_a100"]      = kRed+1;
_gf_color["nhf_multijet_a100"]      = kGreen+3;
_gf_color["nef_multijet_a100"]      = kBlue+1;
_gf_color["cef_multijet_a100"]      = kCyan+2;
_gf_color["muf_multijet_a100"]      = kMagenta+2;
//
_gf_color["chf_incjet_a100"]      = kRed+1;
_gf_color["nhf_incjet_a100"]      = kGreen+3;
_gf_color["nef_incjet_a100"]      = kBlue+1;
_gf_color["cef_incjet_a100"]      = kCyan+2;
_gf_color["muf_incjet_a100"]      = kMagenta+2;

_gf_color["chf_gamjet_ren"]      = kRed+1;
_gf_color["nhf_gamjet_ren"]      = kGreen+2;
_gf_color["nef_gamjet_ren"]      = kBlue+1;
_gf_color["chf_pfjet_ren"]      = kRed;
_gf_color["nhf_pfjet_ren"]      = kGreen+1;
_gf_color["nef_pfjet_ren"]      = kBlue;
_gf_color["chf_zjet_ren"]      = kRed+3;//2;
_gf_color["nhf_zjet_ren"]      = kGreen+4;//3;
_gf_color["nef_zjet_ren"]      = kBlue+3;//2;
_gf_color["chf_zlljet_ren"]      = kRed+2;//3;
_gf_color["nhf_zlljet_ren"]      = kGreen+3;//4;
_gf_color["nef_zlljet_ren"]      = kBlue+2;//3;

_gf_color["chf_cmb_ren"]      = kRed;
_gf_color["nhf_cmb_ren"]      = kGreen+2;
_gf_color["nef_cmb_ren"]      = kBlue;

map<string, int> _gf_marker;
_gf_marker["hdm_mpfchs1_multijet"] = kFullTriangleUp;
_gf_marker["ptchs_multijet_a100"]  = kFullTriangleUp;
_gf_marker["mpfchs1_multijet_a100"]  = kFullTriangleUp;
_gf_marker["ptchs_gamjet_a100"]   = kFullSquare;
_gf_marker["mpfchs1_gamjet_a100"] = kFullSquare;
_gf_marker["hdm_mpfchs1_gamjet"]  = kFullSquare;
_gf_marker["hdm_mpfchs1_gamjet2"] = kOpenDiamond; // For PF
_gf_marker["hdm_gamjet"]          = kFullSquare;
_gf_marker["hdm_mpfchs1_jetz"]     = kFullCircle;//kFullStar;
_gf_marker["hdm_mpfchs1_jetz2"]    = kOpenSquare; // For PF
_gf_marker["hdm_mpfchs1_zjav"]     = kFullCircle;//kFullStar;
_gf_marker["hdm_mpfchs1_zjav2"]    = kOpenSquare; // For PF
_gf_marker["hdm_mpfchs1_zjet"]     = kFullCircle;//kFullStar;
_gf_marker["hdm_mpfchs1_zjet2"]    = kOpenSquare; // For PF
_gf_marker["ptchs_zjet_a100"]     = kFullCircle;//kFullStar;
_gf_marker["mpfchs1_zjet_a100"]     = kFullCircle;//kFullStar;
_gf_marker["hdm_mpfchs1_zlljet"]   = kFullStar;
_gf_marker["hdm_hadw"]             = kFullCircle;
_gf_marker["hdm_cmb"]              = kFullDiamond;
_gf_marker["hdm_cmb_mj"]           = kFullDiamond;
_gf_marker["hdm_incjet"]           = kFullDiamond;
_gf_marker["chf_zjet_a100"]      = kFullSquare;
_gf_marker["nhf_zjet_a100"]      = kFullSquare;
_gf_marker["nef_zjet_a100"]      = kFullSquare;
_gf_marker["cef_zjet_a100"]      = kFullDiamond;
_gf_marker["muf_zjet_a100"]      = kFullDiamond;
_gf_marker["chf_gamjet_a100"]      = kFullDiamond;
_gf_marker["nhf_gamjet_a100"]      = kFullDiamond;
_gf_marker["nef_gamjet_a100"]      = kFullDiamond;
_gf_marker["chf_zlljet_a100"]      = kFullSquare;
_gf_marker["nhf_zlljet_a100"]      = kFullSquare;
_gf_marker["nef_zlljet_a100"]      = kFullSquare;
_gf_marker["cef_zlljet_a100"]      = kFullDiamond;
_gf_marker["muf_zlljet_a100"]      = kFullDiamond;
_gf_marker["cef_zmmjet_a100"]      = kFullDiamond;
_gf_marker["muf_zmmjet_a100"]      = kFullDiamond;
_gf_marker["cef_gamjet_a100"]      = kFullDiamond;
_gf_marker["muf_gamjet_a100"]      = kFullDiamond;
_gf_marker["cef_pfjet_a30"]      = kFullDiamond;
_gf_marker["muf_pfjet_a30"]      = kFullDiamond;

_gf_marker["chf_multijet_a100"]      = kFullSquare;
_gf_marker["nhf_multijet_a100"]      = kFullSquare;
_gf_marker["nef_multijet_a100"]      = kFullSquare;
_gf_marker["cef_multijet_a100"]      = kFullDiamond;
_gf_marker["muf_multijet_a100"]      = kFullDiamond;

_gf_marker["chf_incjet_a100"]      = kFullSquare;
_gf_marker["nhf_incjet_a100"]      = kFullSquare;
_gf_marker["nef_incjet_a100"]      = kFullSquare;
_gf_marker["cef_incjet_a100"]      = kFullDiamond;
_gf_marker["muf_incjet_a100"]      = kFullDiamond;

_gf_marker["chf_zjet_ren"]      = kFullSquare;
_gf_marker["nhf_zjet_ren"]      = kFullSquare;
_gf_marker["nef_zjet_ren"]      = kFullSquare;
_gf_marker["chf_zlljet_ren"]      = kFullSquare;
_gf_marker["nhf_zlljet_ren"]      = kFullSquare;
_gf_marker["nef_zlljet_ren"]      = kFullSquare;

_gf_marker["chf_cmb_ren"]      = kFullCircle;
_gf_marker["nhf_cmb_ren"]      = kFullDiamond;
_gf_marker["nef_cmb_ren"]      = kFullSquare;

map<string, double> _gf_size;
_gf_size["hdm_mpfchs1_multijet"] = 1.0;
_gf_size["ptchs_multijet_a100"] = 0.8;
_gf_size["mpfchs1_multijet_a100"] = 0.8;
_gf_size["ptchs_gamjet_a100"]     = 0.8;
_gf_size["mpfchs1_gamjet_a100"]   = 0.8;
_gf_size["hdm_mpfchs1_gamjet"]    = 0.8;
_gf_size["hdm_gamjet"]           = 0.8;
_gf_size["hdm_mpfchs1_jetz"]     = 1.0;//1.4;
_gf_size["hdm_mpfchs1_zjav"]     = 1.0;//1.4;
_gf_size["hdm_mpfchs1_zjet"]     = 1.0;//1.4;
_gf_size["ptchs_zjet_a100"]     = 1.0;//1.4;
_gf_size["mpfchs1_zjet_a100"]     = 1.0;//1.4;
_gf_size["hdm_mpfchs1_zlljet"]   = 1.4;
_gf_size["hdm_hadw"]             = 1.0;

_gf_size["chf_multijet_a100"] = 0.8;
_gf_size["nef_multijet_a100"] = 0.8;
_gf_size["nhf_multijet_a100"] = 0.8;

_gf_size["chf_incjet_a100"] = 0.6;
_gf_size["nef_incjet_a100"] = 0.6;
_gf_size["nhf_incjet_a100"] = 0.6;

map<string, const char*> _gf_label;
_gf_label["hdm_mpfchs1_multijet"] = "Multijet (HDM)";//(p_{T,HDM}^{avp})";//(p_{T}^{leading})";
_gf_label["ptchs_multijet_a100"] = "Multijet (MJB)";//(p_{T,MJB}^{avp})";
_gf_label["mpfchs1_multijet_a100"] = "Multijet (MPF)";//(p_{T,MPF}^{avp})";
_gf_label["ptchs_gamjet_a100"] = "#gamma+jet (DB)";
_gf_label["mpfchs1_gamjet_a100"] = "#gamma+jet (MPF)";
_gf_label["hdm_mpfchs1_gamjet"] = "#gamma+jet (HDM)";
_gf_label["hdm_gamjet"] = "#gamma+jet (HDM)";
_gf_label["hdm_mpfchs1_jetz"] = "X";//(_gf_fitZjetHDM ? "Z+jet (HDM)" : "Not fit: Z+jet (HDM)");// (UH)";
_gf_label["hdm_mpfchs1_zjav"] = "X";//(_gf_fitZjetHDM ? "Z+jet (HDM)" : "Not fit: Z+jet (HDM)");// (UH)";
_gf_label["hdm_mpfchs1_zjet"] = (_gf_fitZjetHDM ? "Z+jet (HDM)" : "Not fit: Z+jet (HDM)");// (UH)";
_gf_label["ptchs_zjet_a100"] = "Z+jet (DB)";
_gf_label["mpfchs1_zjet_a100"] = "Z+jet (MPF)";
_gf_label["hdm_mpfchs1_zlljet"] = "Z+jet (KIT)";
_gf_label["hdm_hadw"] = "W#rightarrowqq'";
_gf_label["hdm_cmb"] = "Combined";
_gf_label["hdm_cmb_mj"] = "Combined";//"Cmb+MJ";
_gf_label["hdm_incjet"] = "Incl. jet";//"Cmb+MJ";

_gf_label["hdm_cmb2"] = "R_{jet}"; // _pf.pdf
_gf_label["hdm_cmb_mj2"] = "R_{jet}"; // _pf.pdf
_gf_label["chf_cmb_ren"] = "CHF";
_gf_label["nhf_cmb_ren"] = "NHF";
_gf_label["nef_cmb_ren"] = "NEF";

_gf_label["hdm_mpfchs1_jetz2"] = "X";//"R_{jet}";
_gf_label["hdm_mpfchs1_zjav2"] = "X";//"R_{jet}";
_gf_label["ptchs_zjet_a1002"] = "R_{jet}";
_gf_label["mpfchs1_zjet_a1002"] = "R_{jet}";
_gf_label["hdm_mpfchs1_zjet2"] = "R_{jet}";
_gf_label["hdm_mpfchs1_gamjet2"] = "X";//R_{jet}";

_gf_label["chf_zjet_a100"] = "CHF";
_gf_label["nhf_zjet_a100"] = "NHF";
_gf_label["nef_zjet_a100"] = "NEF";
_gf_label["chf_gamjet_a100"] = "X";//"CHF";
_gf_label["nhf_gamjet_a100"] = "X";//"NHF";
_gf_label["nef_gamjet_a100"] = "X";//"NEF";
_gf_label["chf_multijet_a100"] = "X";//"CHF";
_gf_label["nhf_multijet_a100"] = "X";//"NHF";
_gf_label["nef_multijet_a100"] = "X";//"NEF";

#endif
