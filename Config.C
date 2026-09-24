// Purpose: common configurations for L2Res.C, L3Res.C/globalFit, JERSF.C
//#ifndef __CONFIG_C__
//#define __CONFIG_C__

#include <map>
#include <string>
std::map<std::string, std::string> mlum;
std::map<std::string, std::string> mfile;

// old and obsolete
mlum["2024B"] = "0.13 fb^{-1}";
mlum["2024C"] = "7.5 fb^{-1}";
mlum["2024D"] = "8.3 fb^{-1}";
mlum["2024BC"] = "7.7 fb^{-1}";
mlum["2024BCD"] = "15.3 fb^{-1}";
mlum["2024CP"] = "7.5 fb^{-1}";
mlum["2024CR"] = "7.5 fb^{-1}";
mlum["2024CS"] = "7.5 fb^{-1}";
mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024FG"] = "65.5 fb^{-1}";
mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024H"] = "5.4 fb^{-1}";
mlum["2024I"] = "11.6 fb^{-1}";

// rootfiles/brilcalc/sumfibs.py

//NIB-level:
/*
mlum["2022B_nib1"] = "0.097 fb^{-1}";
mlum["2022C_nib1"] = "5.0 fb^{-1}";
mlum["2022D_nib1"] = "2.97 fb^{-1}";
mlum["2022E_nib1"] = "5.8 fb^{-1}";
mlum["2022F_nib1"] = "17.8 fb^{-1}";
mlum["2022G_nib1"] = "3.08 fb^{-1}";
mlum["2023B_nib1"] = "0.64 fb^{-1}";
mlum["2023C_nib1"] = "7.2 fb^{-1}";
mlum["2023Cv4_nib1"] = "10.4 fb^{-1}";
mlum["2023Cv4_nib2"] = "0.407 fb^{-1}";
mlum["2023D_nib1"] = "9.7 fb^{-1}";
mlum["2024B_nib1"] = "0.132 fb^{-1}";
mlum["2024C_nib1"] = "7.3 fb^{-1}";
mlum["2024D_nib1"] = "8.0 fb^{-1}";
mlum["2024Ev1_nib1"] = "6.3 fb^{-1}";
mlum["2024Ev2_nib1"] = "5.1 fb^{-1}";
mlum["2024F_nib1"] = "0.90 fb^{-1}";
mlum["2024F_nib2"] = "14.6 fb^{-1}";
mlum["2024F_nib3"] = "12.6 fb^{-1}";
mlum["2024G_nib1"] = "17.1 fb^{-1}";
mlum["2024G_nib2"] = "21.0 fb^{-1}";
mlum["2024H_nib1"] = "5.5 fb^{-1}";
mlum["2024I_nib1"] = "11.8 fb^{-1}";
mlum["2025Cv1_nib1"] = "13.5 fb^{-1}";
mlum["2025Cv2_nib1"] = "8.3 fb^{-1}";
mlum["2025D_nib1"] = "25.3 fb^{-1}";
mlum["2025E_nib1"] = "14.0 fb^{-1}";
mlum["2025F_nib1"] = "26.4 fb^{-1}";
mlum["2025G_nib1"] = "21.8 fb^{-1}";
mlum["2026A_nib1"] = "0.64 fb^{-1}";
mlum["2026B_nib1"] = "15.3 fb^{-1}";
mlum["2026C_nib1"] = "1.18 fb^{-1}";
*/
mlum["2022B_nib1"] = "0.097 fb^{-1}";
mlum["2022C_nib1"] = "5.0 fb^{-1}";
mlum["2022D_nib1"] = "2.97 fb^{-1}";
mlum["2022E_nib1"] = "5.8 fb^{-1}";
mlum["2022F_nib1"] = "17.8 fb^{-1}";
mlum["2022G_nib1"] = "3.09 fb^{-1}";
mlum["2023B_nib1"] = "0.64 fb^{-1}";
mlum["2023C_nib1"] = "7.2 fb^{-1}";
mlum["2023Cv4_nib1"] = "10.3 fb^{-1}";
mlum["2023Cv4_nib2"] = "0.407 fb^{-1}";
mlum["2023D_nib1"] = "9.7 fb^{-1}";
mlum["2024B_nib1"] = "0.132 fb^{-1}";
mlum["2024C_nib1"] = "7.3 fb^{-1}";
mlum["2024D_nib1"] = "8.0 fb^{-1}";
mlum["2024Ev1_nib1"] = "6.3 fb^{-1}";
mlum["2024Ev2_nib1"] = "5.1 fb^{-1}";
mlum["2024F_nib1"] = "0.90 fb^{-1}";
mlum["2024F_nib2"] = "14.6 fb^{-1}";
mlum["2024F_nib3"] = "12.6 fb^{-1}";
mlum["2024G_nib1"] = "17.1 fb^{-1}";
mlum["2024G_nib2"] = "21.0 fb^{-1}";
mlum["2024H_nib1"] = "5.5 fb^{-1}";
mlum["2024I_nib1"] = "11.8 fb^{-1}";
mlum["2025Cv1_nib1"] = "13.3 fb^{-1}";
mlum["2025Cv2_nib1"] = "8.3 fb^{-1}";
mlum["2025D_nib1"] = "26.0 fb^{-1}";
mlum["2025E_nib1"] = "14.1 fb^{-1}";
mlum["2025F_nib1"] = "26.9 fb^{-1}";
mlum["2025G_nib1"] = "22.4 fb^{-1}";
mlum["2026A_nib1"] = "0.64 fb^{-1}";
mlum["2026B_nib1"] = "15.3 fb^{-1}";
mlum["2026C_nib1"] = "2.11 fb^{-1}";
mlum["2026D_nib1"] = "9.9 fb^{-1}";


// ERA-level:
/*
mlum["2022B"] = "0.097 fb^{-1}";
mlum["2022C"] = "5.0 fb^{-1}";
mlum["2022D"] = "2.97 fb^{-1}";
mlum["2022E"] = "5.8 fb^{-1}";
mlum["2022F"] = "17.8 fb^{-1}";
mlum["2022G"] = "3.08 fb^{-1}";
mlum["2023B"] = "0.64 fb^{-1}";
mlum["2023C"] = "18.1 fb^{-1}";
mlum["2023D"] = "9.7 fb^{-1}";
mlum["2024B"] = "0.132 fb^{-1}";
mlum["2024C"] = "7.3 fb^{-1}";
mlum["2024D"] = "8.0 fb^{-1}";
mlum["2024E"] = "11.4 fb^{-1}";
mlum["2024F"] = "28.1 fb^{-1}";
mlum["2024G"] = "38.1 fb^{-1}";
mlum["2024H"] = "5.5 fb^{-1}";
mlum["2024I"] = "11.8 fb^{-1}";
mlum["2025C"] = "21.8 fb^{-1}";
mlum["2025D"] = "25.3 fb^{-1}";
mlum["2025E"] = "14.0 fb^{-1}";
mlum["2025F"] = "26.4 fb^{-1}";
mlum["2025G"] = "21.8 fb^{-1}";
mlum["2025CDEFG"] = "109 fb^{-1}";
mlum["2025DEFG"] = "87.5 fb^{-1}";
*/
mlum["2022B"] = "0.097 fb^{-1}";
mlum["2022C"] = "5.0 fb^{-1}";
mlum["2022D"] = "2.97 fb^{-1}";
mlum["2022E"] = "5.8 fb^{-1}";
mlum["2022F"] = "17.8 fb^{-1}";
mlum["2022G"] = "3.09 fb^{-1}";
mlum["2023B"] = "0.64 fb^{-1}";
mlum["2023C"] = "18.0 fb^{-1}";
mlum["2023D"] = "9.7 fb^{-1}";
mlum["2024B"] = "0.132 fb^{-1}";
mlum["2024C"] = "7.3 fb^{-1}";
mlum["2024D"] = "8.0 fb^{-1}";
mlum["2024E"] = "11.4 fb^{-1}";
mlum["2024F"] = "28.0 fb^{-1}";
mlum["2024G"] = "38.1 fb^{-1}";
mlum["2024H"] = "5.5 fb^{-1}";
mlum["2024I"] = "11.8 fb^{-1}";
mlum["2025C"] = "21.6 fb^{-1}";
mlum["2025D"] = "26.0 fb^{-1}";
mlum["2025E"] = "14.1 fb^{-1}";
mlum["2025F"] = "26.9 fb^{-1}";
mlum["2025G"] = "22.4 fb^{-1}";
mlum["2026A"] = "0.64 fb^{-1}";
mlum["2026B"] = "15.3 fb^{-1}";
mlum["2026C"] = "2.11 fb^{-1}";
mlum["2026D"] = "9.9 fb^{-1}";
mlum["2025CDEFG"] = "111 fb^{-1}";
mlum["2025DEFG"] = "89.4 fb^{-1}";

// EXTRA semi-year levels:
mlum["2024CDE_nib"] = "26.7 fb^{-1}";
mlum["2024FGHI_nib"] = "83.3 fb^{-1}";
mlum["2024_nib"] = "110 fb^{-1}";

// YEAR-level:
/*
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.4 fb^{-1}";
mlum["2024"] = "110 fb^{-1}";
mlum["2025"] = "109 fb^{-1}";
mlum["2026"] = "17.1 fb^{-1}";
*/
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.3 fb^{-1}";
mlum["2024"] = "110 fb^{-1}";
mlum["2025"] = "111 fb^{-1}";
mlum["2026"] = "27.9 fb^{-1}";
mlum["2026BD"] = "25.2 fb^{-1}";

// TOTAL-level:
//mlum["Run3"] = "283 fb^{-1}";
mlum["Run3"] = "312 fb^{-1}";


////////////////////////////////////////////////////////////
// Run2 (Legacy/UL, 13 TeV) luminosities                  //
////////////////////////////////////////////////////////////
// Standard CMS UL golden-JSON values (brilcalc, normtag).
// TODO: cross-check with rootfiles/brilcalc/sumfibs.py before use in
// publication-level plots; per-era values below sum to the year totals
// (2016: 19.51+16.81=36.32, 2017: 41.47, 2018: 59.83).
// 2016 convention here: B-F = preVFP (HIPM), Fpost,G,H = postVFP (no HIPM)

// ERA-level:
mlum["2016Bv1"] = "0.0 fb^{-1}";   // TODO: ver1 has negligible certified lumi, check
mlum["2016Bv2"] = "5.83 fb^{-1}";
mlum["2016B"]   = "5.83 fb^{-1}";  // = ver2, ver1 negligible
mlum["2016C"]   = "2.60 fb^{-1}";
mlum["2016D"]   = "4.29 fb^{-1}";
mlum["2016E"]   = "4.07 fb^{-1}";
mlum["2016F"]   = "2.72 fb^{-1}";  // preVFP (HIPM) part of era F
mlum["2016Fpost"] = "0.42 fb^{-1}"; // postVFP (no HIPM) part of era F
mlum["2016G"]   = "7.65 fb^{-1}";
mlum["2016H"]   = "8.74 fb^{-1}";
//
mlum["2017B"] = "4.80 fb^{-1}";
mlum["2017C"] = "9.57 fb^{-1}";
mlum["2017D"] = "4.25 fb^{-1}";
mlum["2017E"] = "9.31 fb^{-1}";
mlum["2017F"] = "13.54 fb^{-1}";
//
mlum["2018A"] = "14.03 fb^{-1}";
mlum["2018B"] = "7.07 fb^{-1}";
mlum["2018C"] = "6.90 fb^{-1}";
mlum["2018D"] = "31.83 fb^{-1}";

// ERA-COMBINATION level:
mlum["2016BCDEF"]   = "19.5 fb^{-1}";  // preVFP (HIPM)
mlum["2016FGH"]     = "16.8 fb^{-1}";  // postVFP (no HIPM)
mlum["2016BCDEFGH"] = "36.3 fb^{-1}";
mlum["2017BCDEF"]   = "41.5 fb^{-1}";
mlum["2018ABCD"]    = "59.8 fb^{-1}";

// YEAR-level:
mlum["2016"] = "36.3 fb^{-1}";
mlum["2017"] = "41.5 fb^{-1}";
mlum["2018"] = "59.8 fb^{-1}";

// TOTAL-level:
mlum["Run2"] = "138 fb^{-1}";



mlum["24to26C"] = "236 fb^{-1}"; // obsolote
mlum["24to26"] = "249 fb^{-1}";

mlum["2025C0"] = "20.8 fb^{-1}";
mlum["2025CT"] = "TrkRadDam<<20.8 fb^{-1}";

// Jet spike test
mlum["2026BJS"] = "15.3 fb^{-1}";
mlum["2026BNS"] = "15.3 fb^{-1}";
mlum["2026CJS"] = "1.18 fb^{-1}";
mlum["2026CNS"] = "1.18 fb^{-1}";

// FLAVOR
mlum["2024FLAVOR"] = "110 fb^{-1}";
mlum["2025FLAVOR"] = "111 fb^{-1}";
mlum["2026FLAVOR"] = "25.2 fb^{-1}";
mlum["RUN3FLAVOR"] = "246 fb^{-1}";

////////////////////
// File listings  //
////////////////////

//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_v134.root";
//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MG_v150.root";
mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root"; // no JERSF
//mfile["JERC_Summer24MG_MC_NOJERSF"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root"; // no JERSF
mfile["JERC_Summer24MG_NOJERSF_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root"; // no JERSF
//
//mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MC_Flat_v150.root";
//mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MC_Flat_v155.root"; // non-JMENANO?
mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt/Jet_v159/jmenano_mc_out_Summer24MC_Flat22_JME_v159.root";
mfile["JERC_Summer24MC_Flat22_Herwig_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MC_Flat22_Herwig_v155.root";
mfile["JERC_Summer24MC_FlatBase_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_Base_v151_v3.root";
mfile["JERC_Summer24MC_FlatNoDeepCore_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_NoDeepCore_v151_v3.root";
mfile["JERC_Winter25MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Winter25MG_v150.root";
//
mfile["JERC_Summer24MC_Flat22_NoPU_MC"] = "rootfiles/Prompt2025/Jet_v155/jmenano_mc_out_Summer24MC_Flat22_NoPU_v155.root";
//
mfile["JERC_Summer24MC_Flat22_NoPFHC_MC"] = "rootfiles/Prompt/Jet_v163/jmenano_mc_cmb_Summer24MC_Flat22_NoPFH_v2_v163.root";
//
mfile["JERC_Winter25MC_Flat22_PU120_MC"] = "rootfiles/NestorEEZS25/jmenano_mc_out_Winter25MC_Flat22_v154.root"; // PU 0-120, too high as reference
mfile["JERC_Winter25MC_Flat22_MC"] = "rootfiles/NestorEEZS25/jmenano_mc_out_Winter25MC_Flat22_Realistic_v154.root";
mfile["JERC_Winter25MC_Flat_EEZS9p5_MC"] = "rootfiles/NestorEEZS25/jmenano_mc_out_Winter25MC_Flat_EEZS9p5_v154.root";

// Multijet files
// v141 (MC v128): input to Prompt25_V2M
// v143 (MC v143 JRV2M: closure test of Prompt25_V2M, input to Prompt25_V3M
// 143->146: fix Dijet folder binning and improve triggers (Jet110)
// 146/147->150(_v2): more MC's, latest data, what else?
// v151
// v152: golden JSON, rhovsmu tc.
// v153: V3M closure
// v155: new 2024 files (previous was v134)
// 2024/Jet_v155->Jet_v159 -> v165
//mfile["JET_2024_nib_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2022C_DATA_OUT"] = "rootfiles/Prompt2024/v113_2022/jmenano_data_out_2022C_nib1_JME_v113_2022.root";
mfile["JET_2022D_DATA_OUT"] = "rootfiles/Prompt2024/v113_2022/jmenano_data_out_2022D_nib1_JME_v113_2022.root";
mfile["JET_2022E_DATA_OUT"] = "rootfiles/Prompt2024/v113_2022/jmenano_data_out_2022E_nib1_JME_v113_2022.root";
mfile["JET_2022F_DATA_OUT"] = "rootfiles/Prompt2024/v113_2022/jmenano_data_out_2022F_nib1_JME_v113_2022.root";
mfile["JET_2022G_DATA_OUT"] = "rootfiles/Prompt2024/v113_2022/jmenano_data_out_2022G_nib1_JME_v113_2022.root";
//
mfile["JET_2023C1_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Cv1_nib1_v113_2023.root";
mfile["JET_2023C2_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Cv2_nib1_v113_2023.root";
mfile["JET_2023C3_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Cv3_nib1_v113_2023.root";
mfile["JET_2023C4_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Cv4_nib1_v113_2023.root";
mfile["JET_2023D1_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Dv1_nib1_v113_2023.root";
mfile["JET_2023D2_DATA_OUT"] = "rootfiles/Prompt2024/v113_2023/jmenano_data_out_2023Dv2_nib1_v113_2023.root";

//mfile["JET_2024_MC"] = "rootfiles/Prompt/Jet_v161/jmenano_mc_out_Summer24MG_JMENANO_2024CDEFGHI_JERSF_v161.root";
mfile["JET_2024_MC"] = "rootfiles/Prompt/Jet_v170/jmenano_mc_out_Summer24MG_JMENANO_JERSF2024_v170.root";
//
mfile["JET_2024_nib_MC"] = mfile["JET_2024_MC"];
mfile["JET_2024_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024CDEFGHI_JME_v170.root";
mfile["JET_2024_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024CDEFGHI_JME_v170.root";
mfile["JET_2024FGHI_nib_MC"] = mfile["JET_2024_MC"];
mfile["JET_2024FGHI_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024FGHI_JME_v170.root";
mfile["JET_2024FGHI_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024FGHI_JME_v170.root";
mfile["JET_2024CDE_nib_MC"]  = mfile["JET_2024_MC"];
mfile["JET_2024CDE_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024CDE_JME_v170.root";
mfile["JET_2024CDE_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024CDE_JME_v170.root";
//
mfile["JET_2024C_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024C_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024C_Rp_JME_v170.root";
mfile["JET_2024C_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024C_Rp_JME_v170.root";
//
mfile["JET_2024D_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024D_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024D_Rp_JME_v170.root";
mfile["JET_2024D_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024D_Rp_JME_v170.root";
//
mfile["JET_2024E_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024E_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024E_Rp_JME_v170.root";
mfile["JET_2024E_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024E_Rp_JME_v170.root";
//
mfile["JET_2024F_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024F_nib1_JME_v170.root";
mfile["JET_2024F_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024F_nib1_JME_v170.root";
mfile["JET_2024F_nib2_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib2_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024F_nib2_JME_v170.root";
mfile["JET_2024F_nib2_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024F_nib2_JME_v170.root";
mfile["JET_2024F_nib3_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib3_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024F_nib3_JME_v170.root";
mfile["JET_2024F_nib3_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024F_nib3_JME_v170.root";
//
mfile["JET_2024G_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024G_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024G_nib1_JME_v170.root";
mfile["JET_2024G_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024G_nib1_JME_v170.root";
mfile["JET_2024G_nib2_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024G_nib2_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024G_nib2_JME_v170.root";
mfile["JET_2024G_nib2_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024G_nib2_JME_v170.root";
//
mfile["JET_2024H_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024H_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024H_nib1_JME_v170.root"; // was v163_v4,v159
mfile["JET_2024H_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024H_nib1_JME_v170.root";
//
mfile["JET_2024I_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024I_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024I_nib1_JME_v170.root"; // was v163_v4
mfile["JET_2024I_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024I_nib1_JME_v170.root";

//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root"; // no JER SF
//mfile["JET_2025_MC"] = "rootfiles/Prompt2025/Jet_v146/jmenano_mc_out_Winter25MG_v146.root"; // no JER SF
//mfile["JET_2025_MC"] = "rootfiles/Prompt/Jet_v161/jmenano_mc_out_Summer24MG_JMENANO_2025CDEFG_JERSF_v161.root";
mfile["JET_2025_MC"] = "rootfiles/Prompt/Jet_v170/jmenano_mc_out_Summer24MG_JMENANO_JERSF2025_v170.root";
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v143/jmenano_mc_out_Winter25MG_JRSF2025CDE_v143.root"; // with JER SF
// 2025/Jet_v153 -> Jet_v159 -> v166
mfile["JET_2025C_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025C_JME_v170.root";
mfile["JET_2025C_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025C_JME_v170.root";
mfile["JET_2025C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025D_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025D_JME_v170.root";
mfile["JET_2025D_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025D_JME_v170.root";
mfile["JET_2025D_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025E_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025E_JME_v170.root";
mfile["JET_2025E_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025E_JME_v170.root";
mfile["JET_2025E_MC"]       = mfile["JET_2025_MC"];
// v146->v147 for more F
mfile["JET_2025F_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025F_JME_v170.root";
mfile["JET_2025F_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025F_JME_v170.root";
mfile["JET_2025F_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025G_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025G_JME_v170.root";
mfile["JET_2025G_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025G_JME_v170.root"; //placeholder
mfile["JET_2025G_MC"]       = mfile["JET_2025_MC"];
// v141 CDE, v143 CDEF, v146 CDEF as v147 CDEFG placeholder
mfile["JET_2025CDEFG_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2025CDEFG_JME_v170.root";
mfile["JET_2025CDEFG_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2025CDEFG_JME_v170.root";
mfile["JET_2025CDEFG_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025JER_DATA_OUT"] = mfile["JET_2025CDEFG_DATA_OUT"];
mfile["JET_2025JER_DATA_CMB"] = mfile["JET_2025CDEFG_DATA_CMB"];
mfile["JET_2025JER_MC"]       = "rootfiles/Prompt2026/Jet_v158/jmenano_mc_out_Summer24MC_Flat_JMENANO_JERSF_v158.root";

//mfile["JET_2025DEFG_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025DEFG_JME_v153.root";
//mfile["JET_2025DEFG_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025DEFG_JME_v153.root";
//mfile["JET_2025DEFG_MC"]       = mfile["JET_2025_MC"];
//
//mfile["JET_2025C0_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_out_2025C_JME_v145.root";
//mfile["JET_2025C0_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_cmb_2025C_JME_v145.root";
//mfile["JET_2025C0_MC"] = "rootfiles/Prompt2025/Jet_v145/jmenano_mc_out_Winter25MG_v145.root";
//mfile["JET_2025CT_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_out_2025C_Trk_JME_v145.root";
//mfile["JET_2025CT_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_cmb_2025C_Trk_JME_v145.root";
//mfile["JET_2025CT_MC"] = mfile["JET_2025C0_MC"]; 


// Wqq files
// v1m: input to Prompt25_V2M
// V2M: closure of Prompt25_V2M, input to Prompt25_V3M
// V3M: closure of Prompt25_V3M
// e2: updated jet veto map
//mfile["WQQ_2024_nib_MC"] = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root";
mfile["WQQ_2024_nib_MC"] = "rootfiles/Prompt/Wqq_e7/Summer24_TTtoLNu2Q_JMENano_V11M_JER2024nib_e7.root";
mfile["WQQ_2024_nib_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2024CDEFGHI_V11M_Golden_e7.root";
//
mfile["WQQ_2024FGHI_nib_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024FGHI_nib_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024FGHI_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024CDE_nib_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024CDE_nib_DATA"]  = "rootfiles/Prompt/Wqq_e7/Muon_Run2024CDE_ReReco_V11M_Golden_e7.root";
//
mfile["WQQ_2024C_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024C_nib1_DATA"]  = "rootfiles/Prompt/Wqq_e7/Muon_Run2024C_ReReco_V11M_Golden_e7.root";
mfile["WQQ_2024D_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024D_nib1_DATA"]  = "rootfiles/Prompt/Wqq_e7/Muon_Run2024D_ReReco_V11M_Golden_e7.root";
mfile["WQQ_2024E_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024E_nib1_DATA"]  = "rootfiles/Prompt/Wqq_e7/Muon_Run2024E_ReReco_V11M_Golden_e7.root";
//
mfile["WQQ_2024F_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib1_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024F_nib1_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024F_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib2_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024F_nib2_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024F_nib3_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib3_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024F_nib3_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024G_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib1_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024G_nib1_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024G_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib2_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024G_nib2_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024H_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024H_nib1_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024H_Prompt_V11M_Golden_e7.root";
mfile["WQQ_2024I_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024I_nib1_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024I_Prompt_V11M_Golden_e7.root";
//
//mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root"; // Summer24
mfile["WQQ_2025_MC"]     = "rootfiles/Prompt/Wqq_e7/Summer24_TTtoLNu2Q_JMENano_V5M_JER2025CDEFG_e7.root"; // Summer24
mfile["WQQ_2025C_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025C_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025C_Prompt_V5M_Golden_e7.root";
mfile["WQQ_2025D_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025D_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025D_Prompt_V5M_Golden_e7.root";
mfile["WQQ_2025E_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025E_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025E_Prompt_V5M_Golden_e7.root";
mfile["WQQ_2025F_MC"]     = mfile["WQQ_2025_MC"];
//mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt/Wqq_V3M/Muon_Run2025Fv1v2_Prompt_V3M.root"; // was 2025F
mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025F_Prompt_V5M_Golden_e7.root";
//
mfile["WQQ_2025G_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025G_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025G_Prompt_V5M_Golden_e7.root";
//
mfile["WQQ_2025CDEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025CDEFG_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2025CDEFG_Prompt_V5M_Golden_e7.root";
//
mfile["WQQ_2025JER_MC"]     = "rootfiles/Prompt2024/Wqq_e2/Summer24_TTtoLNu2Q_V9M_24V10MCSF.root";
mfile["WQQ_2025JER_DATA"]   = mfile["WQQ_2025CDEFG_DATA"];
//
//mfile["WQQ_2025DEFG_MC"]     = mfile["WQQ_2025_MC"];
//mfile["WQQ_2025DEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025DEFG_Prompt_V3M.root";


// Z+jet files
// v100: input to Prompt25_V2M
// v101: closure of Prompt25_V2M, input to Prompt25_V3M
// v102 + _nomu: add eta asymmetry Rochester, revert to Summer24
// v103 : 2025G full, Summer25
// v104: 2025 final golden, recovered runs
// V3M_v105: V3M closure
// V3M_v107: fixed Q,G tagging
// V3M_v110: moved to JMENano (since v108)
// V9M_v113: 2024 data updated to JMENano
mfile["ZMM_2022C_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2022CD_Zmm_sync_v84.root";
mfile["ZMM_2022D_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2022CD_Zmm_sync_v84.root";
mfile["ZMM_2022E_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2022E_Zmm_sync_v84.root";
mfile["ZMM_2022F_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2022FG_Zmm_sync_v78.root";
mfile["ZMM_2022G_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2022FG_Zmm_sync_v78.root";
//
mfile["ZMM_2023C1_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023C123_Zmm_sync_v84.root";
mfile["ZMM_2023C2_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023C123_Zmm_sync_v84.root";
mfile["ZMM_2023C3_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023C123_Zmm_sync_v84.root";
mfile["ZMM_2023C4_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023C4_Zmm_sync_v84.root";
mfile["ZMM_2023D1_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023D_Zmm_sync_v84.root";
mfile["ZMM_2023D2_DATA"] = "rootfiles/Prompt2024/jme_bplusZ_2023D_Zmm_sync_v84.root";

// v113 -> V10MV4MV1M_v115: closure
//mfile["ZMM_2024_DATAMC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root";
//mfile["ZMM_2024_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//mfile["ZMM_2024_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_DATA"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024_Zmm_V11M_v116.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v115/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2024Smearing.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v115/jme_Zj_2024DYTT_nib1_Zmm_V10M_v115.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116.root"; // noTT
//mfile["ZMM_2024_nib_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
mfile["ZMM_2024_nib_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT, better DY xsec
//mfile["ZMM_2024_nib_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DYTT_Zmm_V9M_v113.root"; // includes JER SF?
//
//mfile["ZMM_2024FGHI_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
//mfile["ZMM_2024FGHI_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024FGHI_nib_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024FGHI_Zmm_V11M_v116.root";
mfile["ZMM_2024FGHI_nib_MC"]   = mfile["ZMM_2024_nib_MC"];
//
//mfile["ZMM_2024CDE_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
//mfile["ZMM_2024CDE_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024CDE_nib_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024CDE_Zmm_V11M_v116.root";
mfile["ZMM_2024CDE_nib_MC"]   = mfile["ZMM_2024_nib_MC"];
//
//mfile["ZMM_2024C_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024C_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024C_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024D_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024D_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024D_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024D_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024D_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024E_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024E_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024E_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024E_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024E_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024F_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024F_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024F_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
//mfile["ZMM_2024F_nib2_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib2_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024F_nib2_Zmm_V11M_v116.root";
mfile["ZMM_2024F_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib3_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root"; // bug 2024-04-22, C instead of F_nib3!
mfile["ZMM_2024F_nib3_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024F_nib3_Zmm_V11M_v116.root";
mfile["ZMM_2024F_nib3_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024G_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root"; // bug, C instead of G_nib1
mfile["ZMM_2024G_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024G_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024G_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024G_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root"; // bug, C instead of G_nib2
mfile["ZMM_2024G_nib2_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024G_nib2_Zmm_V11M_v116.root";
mfile["ZMM_2024G_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024H_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root"; // bug, C instead of H
mfile["ZMM_2024H_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024H_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024H_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024I_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root"; // bug, C instead of I
mfile["ZMM_2024I_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024I_nib1_Zmm_V11M_v116.root";
mfile["ZMM_2024I_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

// v113 (V4MV10MV1M) -> v115 (closure)
//mfile["ZMM_Summer24_MC"]   = "rootfiles/Prompt2025/Zmm_v109/jme_Zj_2025_Zmm_V3M_v109.root";
//mfile["ZMM_Winter25_MC"]   = "rootfiles/Prompt2025/Zmm_v103/jme_Zj_2025_Zmm_v103_nomu.root";
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/jme_Zj_2024DY/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing.root"; // Summer24 MC JEC + 2025 JER SF
// v113: JMENANO+Summer24MC JEC for MC+Winter25MC JEC for data+JER SF
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing.root"; // Summer24 MC JEC + 2025 JER SF
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt/Zmm_v115/jme_Zj_2025DYTT_Zmm_V4M_v115_nomu.root"; // Also DY without TT available
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025DY_Zmm_V5M_v116_nomu.root"; // No TT
//mfile["ZMM_2025MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
mfile["ZMM_2025MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
//mfile["ZMM_2025MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
//
mfile["ZMM_2025C_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025C_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025C_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025D_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025D_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025D_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025E_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025E_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025E_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025F_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025F_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025F_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025G_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025G_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025G_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025CDEFG_DATA"]   = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025_Zmm_V5M_v116_nomu.root";
mfile["ZMM_2025CDEFG_MC"] = mfile["ZMM_2025MC"];
//
/*
mfile["ZMM_2025JER_DATAMC"] = mfile["ZMM_2025CDEFG_DATAMC"]; 
mfile["ZMM_2025JER_DATA"] = mfile["ZMM_2025CDEFG_DATAMC"]; 
mfile["ZMM_2025JER_MC"] = "rootfiles/Prompt2026/Zmm_v110/jme_Zj_2024DY_Zmm_2025V3M_v110_nomu_Smearing_TEST.root";
//
mfile["ZMM_2025DEFG_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025DEFG_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025DEFG_DATA"] = mfile["ZMM_2025DEFG_DATAMC"]; 
mfile["ZMM_2025DEFG_MC"] = mfile["ZMM_2025MC"];
//
mfile["ZMM_2025C0_DATAMC"]   = mfile["ZMM_2025C_DATAMC"];
mfile["ZMM_2025CT_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v101/jme_Zj_2025C_TrkRadDamage_Zmm_v101.root"; // for L2Res
mfile["ZMM_2025CT_DATA"]   = "rootfiles/Prompt2025/Zmm_v101/jme_Zj_2025C_TrkRadDamage_Zmm_v101.root";
mfile["ZMM_2025CT_MC"]   = mfile["ZMM_2025C_DATAMC"];
*/
// v113: JMENANO+Summer24MC KEC for MC+Winter26MC JEC for data+2025 JER SF
// (v114 2026D)
// v115: closure
//mfile["ZMM_2026MC"]     = "rootfiles/Prompt/Zmm_v115/jme_Zj_2026DYTT_Zmm_V1M_v115_nomu.root";
//mfile["ZMM_2026MC"]     = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026DY_Zmm_V2M_v116_nomu.root";
//mfile["ZMM_2026MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
mfile["ZMM_2026MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
//mfile["ZMM_2026MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
mfile["ZMM_2026B_DATA"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026B_Zmm_V2M_v116_nomu.root";
mfile["ZMM_2026B_MC"]     = mfile["ZMM_2026MC"];
//mfile["ZMM_2026C_DATA"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026C_Zmm_V2M_v116_nomu.root"; // was v113 V0M 16042026 (was this incomplete?)
//mfile["ZMM_2026C_DATA"] = "rootfiles/Prompt/Zmm_v118/jme_Zj_Run2026C_Zmm_v118_2026LowPU.root"; // old MC truth+26B res :(
mfile["ZMM_2026C_DATA"] = "rootfiles/Prompt/Zmm_v119/jme_Zj_Run2026C_Zmm_v119_2026LowPU.root"; // new MC truth+26B res
//mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026DY_JER2026C_Zmm_V2M_v116_nomu.root"; // noTT
//mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v119/jme_Zj_2026DY_Zmm_v119_2026LowPU.root"; // noTT
mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v119/jme_Zj_2026DYTT_Zmm_v119_2026LowPU.root"; // withTT
//mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v115/jme_Zj_2026DY_JER2026C_Zmm_V1M_v115_nomu.root"; // noTT
//mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026DY_Zmm_V2M_v116_nomu.root"; // noTT TODO: missing JER2026C
//mfile["ZMM_2026C_MC"]     = "rootfiles/Prompt/Zmm_v115/jme_Zj_2026DY_Zmm_V1M_v115_nomu.root"; // noTT TODO: missing JER2026C
//mfile["ZMM_2026D_DATA"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2026D_Zmm_V0M_v113_nomu_07052026.root";
mfile["ZMM_2026D_DATA"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026D_Zmm_V2M_v116_nomu.root";
mfile["ZMM_2026D_MC"]     = mfile["ZMM_2026MC"];
//
mfile["ZMM_2026BD_DATA"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026BD_Zmm_V2M_v116_nomu.root";
mfile["ZMM_2026BD_MC"]     = mfile["ZMM_2026MC"];



// Photon+jet files
// w60 (MC w54): input to Prompt25_V2M
// w62: closure of Prompt25_V2M, input to Prompt_V3M
// w64: what was this again? added G?
// w65: added Gamjet1 folder for eta asymmetry
// w66: full 2025
// w67: final 2025 golden JSON (remove high PU runs?)
// w68: V3M closure test
//mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024_V9M_w56.root";
// w65,w68,w71->w73 fix RhoVsNpv, add JES vs NHF
mfile["GAM_2022CDE_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022CDE_v32.root";
mfile["GAM_2022CD_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022CD_v32.root";
mfile["GAM_2022C_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022C_v32.root";
mfile["GAM_2022D_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022D_v32.root";
mfile["GAM_2022E_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022E_v32.root";
mfile["GAM_2022FG_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022FG_v32.root";
mfile["GAM_2022F_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022F_v32.root";
mfile["GAM_2022G_DATA"] = "../gamjet/rootfiles/GamHistosFill_data_2022G_v32.root";
//
mfile["GAM_2023C1_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023Cv123_w8.root";
mfile["GAM_2023C2_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023Cv123_w8.root";
mfile["GAM_2023C3_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023Cv123_w8.root";
mfile["GAM_2023C4_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023Cv4_w8.root";
mfile["GAM_2023D1_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023D_w8.root";
mfile["GAM_2023D2_DATA"]= "rootfiles/Summer23_L2L3Res/GamHistosFill_data_2023D_w8.root";
  
//mfile["Gam_2024_MC"] = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root";
mfile["GAM_2024_nib_MC"]= //"rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
  //"rootfiles/Prompt2024/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU   // TO BE UPDATED with proper JER SF
  //"rootfiles/Prompt/Gam_w85/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_w85.root"; // TO BE UPDATED with proper JER SF
  "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root"; // TO BE UPDATED with proper JER SF
  //"rootfiles/Prompt/Gam_w80/GamHistosFill_mc_summer2024P8-jmenano_no-pu_JERSFX_w80.root"; // wrong bin weights, breaks L2Res
//"rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
//
// w81->w83_jmenano->w83_prompt->w84->w87
// Prompt2024/w73->Prompt/Gam_w83-jmenano
//mfile["GAM_2024FGHI_nib_DATA"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024FGHI_V9M_w56.root";
//mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024CDEFGHI-rereco_w73.root";
//mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024CDEFGHI-rereco-prompt_w87.root";//-jmenano_w83.root";
mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024CDEFGHI_w87.root";//-jmenano_w83.root";
//mfile["GAM_2024FGHI_nib_DATA"] = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024FGHI-nib_w73.root";
mfile["GAM_2024FGHI_nib_DATA"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024FGHI_w87.root";//jmenano_w83.root";
mfile["GAM_2024FGHI_nib_MC"]   = mfile["GAM_2024_nib_MC"];
//mfile["GAM_2024CDE_nib_DATA"]  = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024CDE_V9M_w56.root";
//mfile["GAM_2024CDE_nib_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024CDE-rereco_w73.root";
mfile["GAM_2024CDE_nib_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024CDE-rereco_w87.root";//jmenano_w83.root";
mfile["GAM_2024CDE_nib_MC"]    = mfile["GAM_2024_nib_MC"];

//mfile["GAM_2024C_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024C-rereco_w73.root";
mfile["GAM_2024C_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024C-rereco_w87.root";//jmenano_w83.root";
mfile["GAM_2024C_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024D_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024D-rereco_w87.root";//jmenano_w83.root";
mfile["GAM_2024D_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024E_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024E-rereco_w87.root";//jmenano_w83.root";
mfile["GAM_2024E_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Fnib1_w87.root";
mfile["GAM_2024F_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib2_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Fnib2_w87.root";
mfile["GAM_2024F_nib2_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib3_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Fnib3_w87.root";
mfile["GAM_2024F_nib3_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024G_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Gnib1_w87.root";
mfile["GAM_2024G_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024G_nib2_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Gnib2_w87.root";
mfile["GAM_2024G_nib2_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024H_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Hnib1_w87.root";//-jmenano_w83.root";
mfile["GAM_2024H_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024I_nib1_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024Inib1_w87.root";//-jmenano_w83.root";
mfile["GAM_2024I_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

// Prompt2025_w73 -> Promptw83 (closure) -> w87 (flavor fix)
//mfile["GAM_2025_MC"]       = "rootfiles/Prompt2025/Gam_w54/GamHistosFill_mc_winter2025P8_no-pu_w54.root";
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC, withPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
mfile["GAM_2025_MC"]    =
  //"rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
  mfile["GAM_2024_nib_MC"]; // TO BE UPDATED with proper JER SF
  //"rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
//mfile["GAM_Summer24_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
//mfile["GAM_Winter25_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_winter2025P8_no-pu_w73.root"; // Winter24 MC, noPU
  // "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_winter2025P8_no-pu_w65.root"; // Winter24 MC, noPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024QCD_no-pu_w65.root"; // Summer24 MC QCD, noPU, TEST ONLY!!
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosMix_mc_summer2024P8_Summer2024QCD_no-pu_w65.root"; // Summer24 MC, noPU
//mfile["GAM_2025_MIX"]      = "rootfiles/Prompt2025/Gam_w54/GamHistosMix_mc_winter2025P8_Winter2025QCD_no-pu_w54.root";
mfile["GAM_2025_MIX"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_DATA"]    = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025C-jmenano_w87.root";
mfile["GAM_2025C_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025D_DATA"]    = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025D-jmenano_w87.root";
mfile["GAM_2025D_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025D_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025E_DATA"]    = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025E-jmenano_w87.root";
mfile["GAM_2025E_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025E_MIX"]     = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025F_DATA"]    = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025F-jmenano_w87.root";
mfile["GAM_2025F_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025F_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025G_DATA"]    = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025G-jmenano_w87.root";
mfile["GAM_2025G_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025G_MIX"]     = mfile["GAM_2025_MIX"];
// w60: CDE, w62: CDEF, w64: CDEFG, w65: CDEFG (more G), w67: all (B)CDEFG
//mfile["GAM_2025CDEFG_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025CDEFG-jmenano_w87.root";
mfile["GAM_2025CDEFG_DATA"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025CDEFG_w87.root";
//mfile["GAM_2025CDE_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025CDEFG_MC"]      = mfile["GAM_2025_MC"];
//mfile["GAM_2025CDEFG_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC
mfile["GAM_2025CDEFG_MIX"]   = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025JER_DATA"]   = mfile["GAM_2025CDEFG_DATA"]; 
mfile["GAM_2025JER_MC"]     = mfile["GAM_2025_MC"]; // placeholder
mfile["GAM_2025JER_MIX"]    = mfile["GAM_2025_MIX"]; // placeholder
//
//
//mfile["GAM_2025DEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w83/GamHistosFill_data_2025DEFG-jmenano_w83.root";
//mfile["GAM_2025DEFG_MC"]    = mfile["GAM_2025_MC"];
//mfile["GAM_2025DEFG_MIX"]   = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025C0_DATA"]    = "rootfiles/Prompt2025/Gam_w62/GamHistosFill_data_2025C_w62.root";
mfile["GAM_2025C0_MC"]   = mfile["GAM_2025_MC"];
mfile["GAM_2025C0_MIX"]   = mfile["GAM_2025_MIX"];
mfile["GAM_2025CT_DATA"]    = "rootfiles/Prompt2025/Gam_w63/GamHistosFill_data_2025C-TrkRadDamage_w63.root";
mfile["GAM_2025CT_MC"]   = mfile["GAM_2025_MC"]; // placeholder
mfile["GAM_2025CT_MIX"]   = mfile["GAM_2025_MIX"]; // placeholder

// 2026A initial studies
/*
mfile["JET_2026A_DATA_OUT"] = "rootfiles/Prompt2026/2026A/jmenano_data_out_2026A_0_v157.root";
mfile["JET_2026A_DATA_CMB"] = "rootfiles/Prompt2026/2026A/jmenano_data_cmb_2026A_0_v157.root";
mfile["JET_2026A_MC"]       = mfile["JET_2025_MC"];
//
mfile["ZMM_2026A_DATAMC"] = "rootfiles/Prompt2026/2026A/jme_Zj_2026A_Zmm_v110_nomu.root";
mfile["ZMM_2026A_DATA"]   = mfile["ZMM_2026A_DATAMC"]; 
mfile["ZMM_2026A_MC"]     = mfile["ZMM_2025MC"];
//
mfile["GAM_2026A_DATA"]   = "rootfiles/Prompt2026/2026A/GamHistosFill_data_2026A_w74.root";
mfile["GAM_2026A_MC"]     = mfile["GAM_2025_MC"];
mfile["GAM_2026A_MIX"]    = mfile["GAM_2025_MIX"];
//
mfile["WQQ_2026A_DATA"]   = "rootfiles/Prompt2026/2026A/Muon_Run2026A_Prompt.root";
mfile["WQQ_2026A_MC"]     = mfile["WQQ_2025_MC"];
*/

// v163_v2(B), v163_v4(C), v167(D) -> v168
// 2026B first Prompt2026 JECs with new HB+HE+HF+PFHC+Winter26 MC JEC
mfile["JET_2026_MC"] = "rootfiles/Prompt/Jet_v170/jmenano_mc_out_Summer24MG_JMENANO_JERSF2026B_v170.root";
mfile["JET_2026B_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2026B_JME_v170.root";
mfile["JET_2026B_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2026B_JME_v170.root";
mfile["JET_2026B_MC"]       = mfile["JET_2026_MC"];
//mfile["JET_2026B_MC"]       = mfile["JERC_Summer24MG_MC_NOJERSF"];
// v2: old trig mix + high PU contamination, v4: both fixed
//mfile["JET_2026C_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2026C_JME_v170.root";
//mfile["JET_2026C_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2026C_JME_v170.root";
mfile["JET_2026C_DATA_OUT"] = "rootfiles/Prompt/Jet_v175/jmenano_data_out_2026C_JME_v175.root"; // new MC truth + 26B res
mfile["JET_2026C_DATA_CMB"] = "rootfiles/Prompt/Jet_v175/jmenano_data_cmb_2026C_JME_v175.root"; // new MC truth + 26B res
//mfile["JET_2026C_MC"]       = "rootfiles/Prompt/Jet_v170/jmenano_mc_out_Summer24MG_JMENANO_JERSF2026C_v170.root";
mfile["JET_2026C_MC"]       = "rootfiles/Prompt/Jet_v175/jmenano_mc_out_Summer26MG_JMENANO_v175.root";
mfile["JET_2026D_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2026D_JME_v170.root";
mfile["JET_2026D_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2026D_JME_v170.root";
mfile["JET_2026D_MC"]       = mfile["JET_2026_MC"];
//
mfile["JET_2026BD_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2026BD_JME_v170.root";
mfile["JET_2026BD_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2026BD_JME_v170.root";
mfile["JET_2026BD_MC"]       = mfile["JET_2026_MC"];

mfile["GAM_2026_MC"] = mfile["GAM_2025_MC"]; // TO BE UPDATED
mfile["GAM_2026_MIX"] = mfile["GAM_2025_MIX"]; // TO BE UPDATED
//mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt/Gam_w84/GamHistosFill_data_2026B-jmenano_w84.root"; // with L2L3Res, for V0M closure, w79->w83
mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2026B_w87.root";
mfile["GAM_2026B_MC"]     = mfile["GAM_2026_MC"];
mfile["GAM_2026B_MIX"]    = mfile["GAM_2026_MIX"];
//
//mfile["GAM_2026C_DATA"]   = "rootfiles/Prompt/Gam_w84/GamHistosFill_data_2026C-jmenano_w84.root"; // with L2L3Res, new low PU data w79->w83
//mfile["GAM_2026C_DATA"]   = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2026C_w87.root"; // old MC+26C res
mfile["GAM_2026C_DATA"]   = "rootfiles/Prompt/Gam_w90/GamHistosFill_data_2026C_w90.root"; // new MC truth+26B res
mfile["GAM_2026C_MC"]     = mfile["GAM_2026_MC"];
mfile["GAM_2026C_MIX"]    = mfile["GAM_2026_MIX"];
//
//mfile["GAM_2026D_DATA"]   = "rootfiles/Prompt/Gam_w84/GamHistosFill_data_2026D-jmenano_w84.root";
mfile["GAM_2026D_DATA"]   = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2026D_w87.root";
mfile["GAM_2026D_MC"]     = mfile["GAM_2026_MC"];
mfile["GAM_2026D_MIX"]    = mfile["GAM_2026_MIX"];
//
mfile["GAM_2026BD_DATA"]   = "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2026BD_w87.root";
mfile["GAM_2026BD_MC"]     = mfile["GAM_2026_MC"];
mfile["GAM_2026BD_MIX"]    = mfile["GAM_2026_MIX"];

//mfile["WQQ_2026B_DATA"]   = "rootfiles/Prompt/Wqq_e4/Muon_Run2026B_Prompt_V0M_GoldenJSON_e4.root"; // v2: no L2L3Res, e2,e4: L3L2Res
mfile["WQQ_2026B_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2026B_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root";
//mfile["WQQ_2026B_MC"]     = mfile["WQQ_2025_MC"];
//mfile["WQQ_2026B_MC"]     = "rootfiles/Prompt2026/Wqq_e2/Summer24_TTtoLNu2Q_V9M_26BV0MCSF.root"; // Summer24 + JER SF
mfile["WQQ_2026B_MC"]     = "rootfiles/Prompt/Wqq_e7/Summer24_TTtoLNu2Q_JMENano_V2M_JER2026BD_e7.root"; // Summer24 + JER SF
//
//mfile["WQQ_2026C_DATA"]   = "rootfiles/Prompt/Wqq_e4/Muon_Run2026C_Prompt_V0M_CombinedJSON_e4.root"; // v2: no L2L3Res, e2, e3, e4: L3L2Res
//mfile["WQQ_2026C_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2026C_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root";
//mfile["WQQ_2026C_DATA"]   = "rootfiles/Prompt/Wqq_e8/Muon_Run2026C_Prompt_V2M_lowPU._e8.root"; // ?
mfile["WQQ_2026C_DATA"]   = "rootfiles/Prompt/Wqq_e9/Muon_Run2026C_Prompt_V2M_lowPUMC_highPUL2L3_e9.root";
//mfile["WQQ_2026C_MC"]     = "rootfiles/Prompt/Wqq_e7/Summer24_TTtoLNu2Q_JMENano_V2M_JER2026C_e7.root";
mfile["WQQ_2026C_MC"]     = "rootfiles/Prompt/Wqq_e9/Summer26_TTtoLNu2Q_lowPUJMENano_V2M_JER2026C_e8.root";
//
//mfile["WQQ_2026D_DATA"]   = "rootfiles/Prompt/Wqq_e5/Muon_Run2026D_Prompt_V0M_CombinedJSON_e5.root";
mfile["WQQ_2026D_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2026D_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root";
mfile["WQQ_2026D_MC"]     = mfile["WQQ_2026B_MC"];
//
mfile["WQQ_2026BD_DATA"]   = "rootfiles/Prompt/Wqq_e7/Muon_Run2026BD_Prompt_V2M_MLEnhancedGolden_Latest1.6._e7.root";
mfile["WQQ_2026BD_MC"]     = mfile["WQQ_2026B_MC"];


// Jet spike tests
mfile["WQQ_2026BJS_DATA"]=mfile["WQQ_2026BNS_DATA"]=mfile["WQQ_2026B_DATA"];
mfile["WQQ_2026BJS_MC"]  =mfile["WQQ_2026BNS_MC"]  =mfile["WQQ_2026B_MC"];
mfile["GAM_2026BJS_DATA"]=mfile["GAM_2026BNS_DATA"]=mfile["GAM_2026B_DATA"];
mfile["GAM_2026BJS_MC"]  =mfile["GAM_2026BNS_MC"]  =mfile["GAM_2026B_MC"];
mfile["JET_2026BJS_DATA_CMB"]=mfile["JET_2026BNS_DATA_CMB"]=mfile["JET_2026B_DATA_CMB"];
mfile["JET_2026BJS_MC"]  =mfile["JET_2026BNS_MC"]  =mfile["JET_2026B_MC"];
//
mfile["ZMM_2026BJS_DATA"]="rootfiles/Prompt/Zmm_v113_jetspike/jme_Zj_2026B_Zmm_V0M_NewJetSpikeCut_nomu.root";
mfile["ZMM_2026BNS_DATA"]="rootfiles/Prompt/Zmm_v113_jetspike/jme_Zj_2026B_Zmm_V0M_NoJetSpikeCut_nomu.root";
mfile["ZMM_2026CJS_DATA"]="rootfiles/Prompt/Zmm_v113_jetspike/jme_Zj_2026C_Zmm_V0M_NewJetSpikeCut_nomu.root";
mfile["ZMM_2026CNS_DATA"]="rootfiles/Prompt/Zmm_v113_jetspike/jme_Zj_2026C_Zmm_V0M_NoJetSpikeCut_nomu.root";
mfile["ZMM_2026BJS_MC"]=mfile["ZMM_2026BNS_MC"]=mfile["ZMM_2026B_MC"];
mfile["ZMM_2026CJS_MC"]=mfile["ZMM_2026CJS_MC"]=mfile["ZMM_2026C_MC"];


// FSR studies
mfile["ZMM_NOPS_MC"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing.root";
mfile["ZMM_PSW4_MC"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing_PSWeight4.root";
//mfile["ZMM_NOPS_DATA"] = mfile["ZMM_PSW4_DATA"] = mfile["ZMM_2024_nib_DATA"];
mfile["ZMM_NOPS_DATA"] = mfile["ZMM_PSW4_DATA"] = mfile["ZMM_2025CDEFG_DATA"];

mfile["GAM_NOPS_MC"] = "rootfiles/Prompt/Gam_w85/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_w85.root";
mfile["GAM_PSW4_MC"] = "rootfiles/Prompt/Gam_w86/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_psweightIndex4_w86.root";
//mfile["GAM_NOPS_DATA"] = mfile["GAM_PSW4_DATA"] = mfile["GAM_2024_nib_DATA"];
mfile["GAM_NOPS_DATA"] = mfile["GAM_PSW4_DATA"] = mfile["GAM_2025CDEFG_DATA"];

mfile["WQQ_NOPS_MC"] = "rootfiles/Prompt2024/Wqq_e2/Summer24_TTtoLNu2Q_V9M_24V10MCSF.root";
mfile["WQQ_PSW4_MC"] = "rootfiles/Prompt2024/Wqq_e2/Summer24_TTtoLNu2Q_V9M_FSR.root"; 
//mfile["WQQ_NOPS_DATA"] = mfile["WQQ_PSW4_DATA"] = mfile["WQQ_2024_nib_DATA"];
mfile["WQQ_NOPS_DATA"] = mfile["WQQ_PSW4_DATA"] = mfile["WQQ_2025CDEFG_DATA"];

mfile["JET_NOPS_MC"] = "rootfiles/Prompt/Jet_v166/jmenano_mc_out_Summer24MG_JMENANO_JERSF2025_v166.root";
mfile["JET_PSW4_MC"] = "rootfiles/Prompt/Jet_v166/jmenano_mc_out_Summer24MG_JMENANO_P25_JERSF2025_v166.root";
//mfile["JET_NOPS_DATA_CMB"] = mfile["JET_PSW4_DATA_CMB"] = mfile["JET_2024_nib_DATA_CMB"];
mfile["JET_NOPS_DATA_CMB"] = mfile["JET_PSW4_DATA_CMB"] = mfile["JET_2025CDEFG_DATA_CMB"];


// Flavor studies
mfile["JET_2024FLAVOR_DATA_CMB"] = mfile["JET_2024_nib_DATA_CMB"];
mfile["JET_2024FLAVOR_MC"] = //mfile["JET_2024_nib_MC"];
  "rootfiles/Prompt/Jet_v171/jmenano_mc_out_Summer24MG_JMENANO_v171.root";
mfile["WQQ_2024FLAVOR_DATA"] = mfile["WQQ_2024_nib_DATA"];
mfile["WQQ_2024FLAVOR_MC"] = mfile["WQQ_2024_nib_MC"];
mfile["ZMM_2024FLAVOR_DATA"] = mfile["ZMM_2024_nib_DATA"];
//mfile["ZMM_2024FLAVOR_MC"] = mfile["ZMM_2024_nib_MC"];
//mfile["ZMM_2024FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116.root"; // with TT
//mfile["ZMM_2024FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
//mfile["ZMM_2024FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116b.root"; // no TT
mfile["ZMM_2024FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
//mfile["ZMM_2024FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
mfile["GAM_2024FLAVOR_DATA"] = mfile["GAM_2024_nib_DATA"];
mfile["GAM_2024FLAVOR_MC"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root";

mfile["JET_2025FLAVOR_DATA_CMB"] = mfile["JET_2025CDEFG_DATA_CMB"];
mfile["JET_2025FLAVOR_MC"] = //mfile["JET_2024_nib_MC"];
  "rootfiles/Prompt/Jet_v171/jmenano_mc_out_Summer24MG_JMENANO_v171.root";
mfile["WQQ_2025FLAVOR_DATA"] = mfile["WQQ_2025CDEFG_DATA"];
mfile["WQQ_2025FLAVOR_MC"] = mfile["WQQ_2025CDEFG_MC"];
mfile["ZMM_2025FLAVOR_DATA"] = mfile["ZMM_2025CDEFG_DATA"];
//mfile["ZMM_2025FLAVOR_MC"] = mfile["ZMM_2025CDEFG_MC"];
//mfile["ZMM_2025FLAVOR_MC"]  = "rootfiles/Prompt/Zmm_v116/jme_Zj_2025DYTT_Zmm_V5M_v116_nomu.root"; // with TT
//mfile["ZMM_2025FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
mfile["ZMM_2025FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
//mfile["ZMM_2025FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
mfile["GAM_2025FLAVOR_DATA"] = mfile["GAM_2025CDEFG_DATA"];
  //"rootfiles/Prompt/Gam_w87/GamHistosFill_data_2025CDEFG_w87.root";
mfile["GAM_2025FLAVOR_MC"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root";

mfile["JET_2026FLAVOR_DATA_CMB"] = mfile["JET_2026BD_DATA_CMB"];
mfile["JET_2026FLAVOR_MC"] = //mfile["JET_2026_MC];
  "rootfiles/Prompt/Jet_v171/jmenano_mc_out_Summer24MG_JMENANO_v171.root";
mfile["WQQ_2026FLAVOR_DATA"] = mfile["WQQ_2026BD_DATA"];
mfile["WQQ_2026FLAVOR_MC"] = mfile["WQQ_2026BD_MC"];
mfile["ZMM_2026FLAVOR_DATA"] = mfile["ZMM_2026BD_DATA"];
//mfile["ZMM_2026FLAVOR_MC"] = mfile["ZMM_2026BD_MC"];
//mfile["ZMM_2026FLAVOR_MC"]  = "rootfiles/Prompt/Zmm_v116/jme_Zj_2026DYTT_Zmm_V2M_v116_nomu.root"; // with TT
//mfile["ZMM_2026FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116b.root"; // with TT => v116b works best
mfile["ZMM_2026FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
//mfile["ZMM_2026FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DY_Zmm_V11M_v116c.root"; // no TT
mfile["GAM_2026FLAVOR_DATA"] = mfile["GAM_2026BD_DATA"];
mfile["GAM_2026FLAVOR_MC"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root";


mfile["JET_RUN3FLAVOR_DATA_CMB"] = //"rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024to2025_JME_v170.root";
  "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024to2026_JME_v170.root";
mfile["JET_RUN3FLAVOR_MC"] = "rootfiles/Prompt/Jet_v171/jmenano_mc_out_Summer24MG_JMENANO_v171.root";
mfile["WQQ_RUN3FLAVOR_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_2024to2025_e7.root";
mfile["WQQ_RUN3FLAVOR_MC"] = mfile["WQQ_2024_nib_MC"];
mfile["ZMM_RUN3FLAVOR_DATA"] = //"rootfiles/Prompt/Zmm_v116/jme_Zj_2024to2025_Zmm_v116.root";
  "rootfiles/Prompt/Zmm_v116/jme_Zj_Run2024CDEFGHI_Run2025CDEFG_Run2026BD_Zmm_v116.root";
mfile["ZMM_RUN3FLAVOR_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // with TT
mfile["GAM_RUN3FLAVOR_DATA"] = //"rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024to2025_w87.root";
  "rootfiles/Prompt/Gam_w87/GamHistosFill_data_2024to2026_w87.root";
mfile["GAM_RUN3FLAVOR_MC"] = "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root";


// Test new Z+jet method on 2024I+Summer24 MC (nib1->nix)
mfile["JET_2024I_nix_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024I_nix_DATA_OUT"] = "rootfiles/Prompt/Jet_v170/jmenano_data_out_2024I_nib1_JME_v170.root"; // was v163_v4
mfile["JET_2024I_nix_DATA_CMB"] = "rootfiles/Prompt/Jet_v170/jmenano_data_cmb_2024I_nib1_JME_v170.root";
//
mfile["WQQ_2024I_nix_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024I_nix_DATA"] = "rootfiles/Prompt/Wqq_e7/Muon_Run2024I_Prompt_V11M_Golden_e7.root";
//
mfile["GAM_2024I_nix_DATA"]  = "rootfiles/Prompt/Gam_w84/GamHistosFill_data_2024Inib1_w84.root";//-jmenano_w83.root";
mfile["GAM_2024I_nix_MC"]    = mfile["GAM_2024_nib_MC"];
//
mfile["ZMM_2024I_nix_DATA"]   = "../zjet/rootfiles/zjet_JMENANO_compat.root";
//mfile["ZMM_2024I_nix_DATA"]   = "../zjet/rootfiles/zjet_JMENANO_run2024i_sync_sami_20260813_v1_legacy.root";
//mfile["ZMM_2024I_nix_DATA"]   = "../zjet/rootfiles/zjet_JMENANO_run2024i_sync_alpha_20260813_v1_legacy.root";
mfile["ZMM_2024I_nix_MC"]   = mfile["ZMM_2024I_nix_DATA"];


////////////////////////////////////////////////////////////
// Run2 (Legacy/UL) file listings                         //
////////////////////////////////////////////////////////////
// Added from the file listing received from colleagues (rootfiles/Run2/...)
// 2016 convention: eras B-F are preVFP (HIPM), Fpost/G/H are postVFP (no HIPM).
// NOTE: channel coverage is incomplete, see the TODOs below. In particular
//       there is currently NO epoch for which all four channels exist:
//       2016BCDEF   : JET, ZMM        (no GAM, no WQQ)
//       2016FGH     : JET, ZMM, WQQ   (no GAM)
//       2016BCDEFGH : JET, GAM        (no ZMM, no WQQ)
//       2017BCDEF   : JET, GAM, ZMM   (no WQQ combination file)
//       2018ABCD    : JET, GAM, ZMM   (no WQQ combination file)

// Multijet files (Jet_v172)
// First round: Run2 data vs Run3 Summer24 MC (same MC as the FLAVOR epochs
// above, without era-specific JER SF). The Run2 UL MCs from the listing are
// kept commented below as a later alternative.
mfile["JET_RUN2_MC"] = "rootfiles/Prompt/Jet_v171/jmenano_mc_out_Summer24MG_JMENANO_v171.root"; // Summer24, no JER SF
mfile["JET_2016HIPM_MC"] = mfile["JET_RUN2_MC"];
mfile["JET_2016_MC"]     = mfile["JET_RUN2_MC"];
mfile["JET_2017_MC"]     = mfile["JET_RUN2_MC"];
mfile["JET_2018_MC"]     = mfile["JET_RUN2_MC"];
// Run2 UL MC (Jet_v172), mc_out as used by JET_*_MC elsewhere. NB: the
// HIPM / non-HIPM split of the year-level keys only matters once these are
// enabled; with Summer24 all four point to the same file.
//mfile["JET_2016HIPM_MC"] = "rootfiles/Run2/Jet_v172/jmenano_mc_out_Summer20UL16MG_HIPM_v172.root";
//mfile["JET_2016_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_out_Summer20UL16MG_v172.root";
//mfile["JET_2017_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_out_Summer20UL17MG_v172.root";
//mfile["JET_2018_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_out_Summer20UL18MG_v172.root";
// ... and the same as mc_cmb:
//mfile["JET_2016HIPM_MC"] = "rootfiles/Run2/Jet_v172/jmenano_mc_cmb_Summer20UL16MG_HIPM_v172.root";
//mfile["JET_2016_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_cmb_Summer20UL16MG_v172.root";
//mfile["JET_2017_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_cmb_Summer20UL17MG_v172.root";
//mfile["JET_2018_MC"]     = "rootfiles/Run2/Jet_v172/jmenano_mc_cmb_Summer20UL18MG_v172.root";
//
mfile["JET_2016B_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016B_HIPM_JME_v172.root";
mfile["JET_2016B_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016B_HIPM_JME_v172.root";
mfile["JET_2016B_MC"]       = mfile["JET_2016HIPM_MC"];
mfile["JET_2016Bv1_DATA_OUT"] = mfile["JET_2016B_DATA_OUT"]; // no ver1/ver2 split in Jet
mfile["JET_2016Bv1_DATA_CMB"] = mfile["JET_2016B_DATA_CMB"];
mfile["JET_2016Bv1_MC"]       = mfile["JET_2016HIPM_MC"];
mfile["JET_2016Bv2_DATA_OUT"] = mfile["JET_2016B_DATA_OUT"];
mfile["JET_2016Bv2_DATA_CMB"] = mfile["JET_2016B_DATA_CMB"];
mfile["JET_2016Bv2_MC"]       = mfile["JET_2016HIPM_MC"];
//
mfile["JET_2016C_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016C_HIPM_JME_v172.root";
mfile["JET_2016C_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016C_HIPM_JME_v172.root";
mfile["JET_2016C_MC"]       = mfile["JET_2016HIPM_MC"];
//
mfile["JET_2016D_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016D_HIPM_JME_v172.root";
mfile["JET_2016D_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016D_HIPM_JME_v172.root";
mfile["JET_2016D_MC"]       = mfile["JET_2016HIPM_MC"];
//
mfile["JET_2016E_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016E_HIPM_JME_v172.root";
mfile["JET_2016E_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016E_HIPM_JME_v172.root";
mfile["JET_2016E_MC"]       = mfile["JET_2016HIPM_MC"];
//
mfile["JET_2016F_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016F_HIPM_JME_v172.root";
mfile["JET_2016F_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016F_HIPM_JME_v172.root";
mfile["JET_2016F_MC"]       = mfile["JET_2016HIPM_MC"];
mfile["JET_2016Fpost_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016F_JME_v172.root";
mfile["JET_2016Fpost_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016F_JME_v172.root";
mfile["JET_2016Fpost_MC"]       = mfile["JET_2016_MC"];
//
mfile["JET_2016G_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016G_JME_v172.root";
mfile["JET_2016G_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016G_JME_v172.root";
mfile["JET_2016G_MC"]       = mfile["JET_2016_MC"];
//
mfile["JET_2016H_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016H_JME_v172.root";
mfile["JET_2016H_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016H_JME_v172.root";
mfile["JET_2016H_MC"]       = mfile["JET_2016_MC"];
//
mfile["JET_2016BCDEF_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016BCDEF_HIPM_JME_v172.root";
mfile["JET_2016BCDEF_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016BCDEF_HIPM_JME_v172.root";
mfile["JET_2016BCDEF_MC"]       = mfile["JET_2016HIPM_MC"];
mfile["JET_2016FGH_DATA_OUT"]   = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016FGH_JME_v172.root";
mfile["JET_2016FGH_DATA_CMB"]   = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016FGH_JME_v172.root";
mfile["JET_2016FGH_MC"]         = mfile["JET_2016_MC"];
mfile["JET_2016BCDEFGH_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2016BCDEFGH_JME_v172.root";
mfile["JET_2016BCDEFGH_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2016BCDEFGH_JME_v172.root";
mfile["JET_2016BCDEFGH_MC"]       = mfile["JET_RUN2_MC"]; // Summer24, as for Run3
//
mfile["JET_2017B_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017B_JME_v172.root";
mfile["JET_2017B_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017B_JME_v172.root";
mfile["JET_2017B_MC"]       = mfile["JET_2017_MC"];
mfile["JET_2017C_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017C_JME_v172.root";
mfile["JET_2017C_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017C_JME_v172.root";
mfile["JET_2017C_MC"]       = mfile["JET_2017_MC"];
mfile["JET_2017D_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017D_JME_v172.root";
mfile["JET_2017D_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017D_JME_v172.root";
mfile["JET_2017D_MC"]       = mfile["JET_2017_MC"];
mfile["JET_2017E_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017E_JME_v172.root";
mfile["JET_2017E_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017E_JME_v172.root";
mfile["JET_2017E_MC"]       = mfile["JET_2017_MC"];
mfile["JET_2017F_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017F_JME_v172.root";
mfile["JET_2017F_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017F_JME_v172.root";
mfile["JET_2017F_MC"]       = mfile["JET_2017_MC"];
mfile["JET_2017BCDEF_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2017BCDEF_JME_v172.root";
mfile["JET_2017BCDEF_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2017BCDEF_JME_v172.root";
mfile["JET_2017BCDEF_MC"]       = mfile["JET_2017_MC"];
//
mfile["JET_2018A_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2018A_JME_v172.root";
mfile["JET_2018A_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2018A_JME_v172.root";
mfile["JET_2018A_MC"]       = mfile["JET_2018_MC"];
mfile["JET_2018B_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2018B_JME_v172.root";
mfile["JET_2018B_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2018B_JME_v172.root";
mfile["JET_2018B_MC"]       = mfile["JET_2018_MC"];
mfile["JET_2018C_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2018C_JME_v172.root";
mfile["JET_2018C_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2018C_JME_v172.root";
mfile["JET_2018C_MC"]       = mfile["JET_2018_MC"];
mfile["JET_2018D_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2018D_JME_v172.root";
mfile["JET_2018D_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2018D_JME_v172.root";
mfile["JET_2018D_MC"]       = mfile["JET_2018_MC"];
mfile["JET_2018ABCD_DATA_OUT"] = "rootfiles/Run2/Jet_v172/jmenano_data_out_2018ABCD_JME_v172.root";
mfile["JET_2018ABCD_DATA_CMB"] = "rootfiles/Run2/Jet_v172/jmenano_data_cmb_2018ABCD_JME_v172.root";
mfile["JET_2018ABCD_MC"]       = mfile["JET_2018_MC"];


// Photon+jet files (Gam_w88)
// First round: Run2 data vs Run3 Summer24 MC. reprocess.C reads the Run2
// flavour histograms from flavor/ (the Run2 files have flavor/ and flavor_new/,
// no flavor_old/), so this uses the same file as the FLAVOR epochs.
mfile["GAM_RUN2_MC"]  = "rootfiles/Prompt/Gam_w87/GamHistosFill_mc_summer2024P8-jmenano_no-pu_no-jersf_no-psweight_w87.root";
//mfile["GAM_RUN2_MC"] = mfile["GAM_2025_MC"]; // Gam_w73, read from flavor_old/
mfile["GAM_RUN2_MIX"] = mfile["GAM_RUN2_MC"];
mfile["GAM_2016HIPM_MC"] = mfile["GAM_RUN2_MC"];
mfile["GAM_2016_MC"]     = mfile["GAM_RUN2_MC"];
mfile["GAM_2017_MC"]     = mfile["GAM_RUN2_MC"];
mfile["GAM_2018_MC"]     = mfile["GAM_RUN2_MC"];
mfile["GAM_2016HIPM_MIX"] = mfile["GAM_RUN2_MIX"];
mfile["GAM_2016_MIX"]     = mfile["GAM_RUN2_MIX"];
mfile["GAM_2017_MIX"]     = mfile["GAM_RUN2_MIX"];
mfile["GAM_2018_MIX"]     = mfile["GAM_RUN2_MIX"];
// No Run2 gamma+jet MC delivered; guessed paths kept for the later round:
//mfile["GAM_2016HIPM_MC"] = "rootfiles/Run2/Gam_w88/GamHistosFill_mc_summer20UL16P8_HIPM-jmenano_no-pu_w88.root"; // PLACEHOLDER
//mfile["GAM_2016_MC"]     = "rootfiles/Run2/Gam_w88/GamHistosFill_mc_summer20UL16P8-jmenano_no-pu_w88.root"; // PLACEHOLDER
//mfile["GAM_2017_MC"]     = "rootfiles/Run2/Gam_w88/GamHistosFill_mc_summer20UL17P8-jmenano_no-pu_w88.root"; // PLACEHOLDER
//mfile["GAM_2018_MC"]     = "rootfiles/Run2/Gam_w88/GamHistosFill_mc_summer20UL18P8-jmenano_no-pu_w88.root"; // PLACEHOLDER
//
mfile["GAM_2016Bv1_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016Bv1-jmenano_w88.root";
mfile["GAM_2016Bv1_MC"]   = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016Bv1_MIX"]  = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016Bv2_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016Bv2-jmenano_w88.root";
mfile["GAM_2016Bv2_MC"]   = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016Bv2_MIX"]  = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016B_DATA"]   = mfile["GAM_2016Bv2_DATA"]; // ver1 lumi negligible
mfile["GAM_2016B_MC"]     = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016B_MIX"]    = mfile["GAM_2016HIPM_MIX"];
//
mfile["GAM_2016C_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016C-jmenano_w88.root";
mfile["GAM_2016C_MC"]     = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016C_MIX"]    = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016D_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016D-jmenano_w88.root";
mfile["GAM_2016D_MC"]     = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016D_MIX"]    = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016E_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016E-jmenano_w88.root";
mfile["GAM_2016E_MC"]     = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016E_MIX"]    = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016F_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016Fhipm-jmenano_w88.root";
mfile["GAM_2016F_MC"]     = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016F_MIX"]    = mfile["GAM_2016HIPM_MIX"];
//
mfile["GAM_2016Fpost_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016Fnohipm-jmenano_w88.root";
mfile["GAM_2016Fpost_MC"]   = mfile["GAM_2016_MC"];
mfile["GAM_2016Fpost_MIX"]  = mfile["GAM_2016_MIX"];
mfile["GAM_2016G_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016G-jmenano_w88.root";
mfile["GAM_2016G_MC"]     = mfile["GAM_2016_MC"];
mfile["GAM_2016G_MIX"]    = mfile["GAM_2016_MIX"];
mfile["GAM_2016H_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016H-jmenano_w88.root";
mfile["GAM_2016H_MC"]     = mfile["GAM_2016_MC"];
mfile["GAM_2016H_MIX"]    = mfile["GAM_2016_MIX"];
//
// hadd of B(v2)+C+D+E+Fhipm resp. Fnohipm+G+H:
mfile["GAM_2016BCDEF_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016BCDEF-jmenano_w88.root"; // PLACEHOLDER
mfile["GAM_2016BCDEF_MC"]   = mfile["GAM_2016HIPM_MC"];
mfile["GAM_2016BCDEF_MIX"]  = mfile["GAM_2016HIPM_MIX"];
mfile["GAM_2016FGH_DATA"]   = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016FGH-jmenano_w88.root"; // PLACEHOLDER
mfile["GAM_2016FGH_MC"]     = mfile["GAM_2016_MC"];
mfile["GAM_2016FGH_MIX"]    = mfile["GAM_2016_MIX"];
mfile["GAM_2016BCDEFGH_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2016BCDEFGH-jmenano_w88.root";
mfile["GAM_2016BCDEFGH_MC"]   = mfile["GAM_RUN2_MC"]; // Summer24, as for Run3
mfile["GAM_2016BCDEFGH_MIX"]  = mfile["GAM_RUN2_MIX"]; // Summer24, as for Run3
//
mfile["GAM_2017B_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017B-jmenano_w88.root";
mfile["GAM_2017B_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017B_MIX"]  = mfile["GAM_2017_MIX"];
mfile["GAM_2017C_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017C-jmenano_w88.root";
mfile["GAM_2017C_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017C_MIX"]  = mfile["GAM_2017_MIX"];
mfile["GAM_2017D_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017D-jmenano_w88.root";
mfile["GAM_2017D_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017D_MIX"]  = mfile["GAM_2017_MIX"];
mfile["GAM_2017E_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017E-jmenano_w88.root";
mfile["GAM_2017E_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017E_MIX"]  = mfile["GAM_2017_MIX"];
mfile["GAM_2017F_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017F-jmenano_w88.root";
mfile["GAM_2017F_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017F_MIX"]  = mfile["GAM_2017_MIX"];
mfile["GAM_2017BCDEF_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2017BCDEF-jmenano_w88.root";
mfile["GAM_2017BCDEF_MC"]   = mfile["GAM_2017_MC"];
mfile["GAM_2017BCDEF_MIX"]  = mfile["GAM_2017_MIX"];
//
mfile["GAM_2018A_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2018A-jmenano_w88.root";
mfile["GAM_2018A_MC"]   = mfile["GAM_2018_MC"];
mfile["GAM_2018A_MIX"]  = mfile["GAM_2018_MIX"];
mfile["GAM_2018B_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2018B-jmenano_w88.root";
mfile["GAM_2018B_MC"]   = mfile["GAM_2018_MC"];
mfile["GAM_2018B_MIX"]  = mfile["GAM_2018_MIX"];
mfile["GAM_2018C_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2018C-jmenano_w88.root";
mfile["GAM_2018C_MC"]   = mfile["GAM_2018_MC"];
mfile["GAM_2018C_MIX"]  = mfile["GAM_2018_MIX"];
mfile["GAM_2018D_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2018D-jmenano_w88.root";
mfile["GAM_2018D_MC"]   = mfile["GAM_2018_MC"];
mfile["GAM_2018D_MIX"]  = mfile["GAM_2018_MIX"];
mfile["GAM_2018ABCD_DATA"] = "rootfiles/Run2/Gam_w88/GamHistosFill_data_2018ABCD-jmenano_w88.root";
mfile["GAM_2018ABCD_MC"]   = mfile["GAM_2018_MC"];
mfile["GAM_2018ABCD_MIX"]  = mfile["GAM_2018_MIX"];


// Z+jet files (Zmm_v116)
// First round: Run2 data vs Run3 Summer24 MC (same file as ZMM_RUN3FLAVOR_MC).
// The Run2 MCs that came with the data (DYLO, DYLOTT, DYNLO, DYNLOTT per epoch)
// are kept commented below as a later alternative; DYNLOTT = NLO DY + TT would
// match the "with TT" choice used for Run3 (v116c).
mfile["ZMM_RUN2_MC"] = "rootfiles/Prompt/Zmm_v116/jme_Zj_2024DYTT_Zmm_V11M_v116c.root"; // Summer24, with TT
mfile["ZMM_2016BCDEF_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEF_Zmm_V7_v116_DATA.root";
mfile["ZMM_2016BCDEF_MC"]   = mfile["ZMM_RUN2_MC"];
//mfile["ZMM_2016BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEF_Zmm_V7_v116_DYNLOTT.root"; // NLO + TT
//mfile["ZMM_2016BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEF_Zmm_V7_v116_DYLOTT.root"; // LO + TT
//mfile["ZMM_2016BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEF_Zmm_V7_v116_DYNLO.root"; // no TT
//mfile["ZMM_2016BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEF_Zmm_V7_v116_DYLO.root"; // LO, no TT
//
mfile["ZMM_2016FGH_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016FGH_Zmm_V7_v116_DATA.root";
mfile["ZMM_2016FGH_MC"]     = mfile["ZMM_RUN2_MC"];
//mfile["ZMM_2016FGH_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016FGH_Zmm_V7_v116_DYNLOTT.root"; // NLO + TT
//mfile["ZMM_2016FGH_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016FGH_Zmm_V7_v116_DYLOTT.root"; // LO + TT
//mfile["ZMM_2016FGH_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016FGH_Zmm_V7_v116_DYNLO.root"; // no TT
//mfile["ZMM_2016FGH_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016FGH_Zmm_V7_v116_DYLO.root"; // LO, no TT
//
mfile["ZMM_2017BCDEF_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017BCDEF_Zmm_V5_v116_DATA.root";
mfile["ZMM_2017BCDEF_MC"]   = mfile["ZMM_RUN2_MC"];
//mfile["ZMM_2017BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017BCDEF_Zmm_V5_v116_DYNLOTT.root"; // NLO + TT
//mfile["ZMM_2017BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017BCDEF_Zmm_V5_v116_DYLOTT.root"; // LO + TT
//mfile["ZMM_2017BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017BCDEF_Zmm_V5_v116_DYNLO.root"; // no TT
//mfile["ZMM_2017BCDEF_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017BCDEF_Zmm_V5_v116_DYLO.root"; // LO, no TT
//
mfile["ZMM_2018ABCD_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018ABCD_Zmm_V5_v116_DATA.root";
mfile["ZMM_2018ABCD_MC"]    = mfile["ZMM_RUN2_MC"];
//mfile["ZMM_2018ABCD_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018ABCD_Zmm_V5_v116_DYNLOTT.root"; // NLO + TT
//mfile["ZMM_2018ABCD_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018ABCD_Zmm_V5_v116_DYLOTT.root"; // LO + TT
//mfile["ZMM_2018ABCD_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018ABCD_Zmm_V5_v116_DYNLO.root"; // no TT
//mfile["ZMM_2018ABCD_MC"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018ABCD_Zmm_V5_v116_DYLO.root"; // LO, no TT
//
// Per-era Z+jet files: requested from the colleagues (cannot be made by hadd
// from the merged files). Paths below are guesses, fix when they arrive.
mfile["ZMM_2016B_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016B_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016B_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016Bv1_DATA"] = mfile["ZMM_2016Bv2_DATA"] = mfile["ZMM_2016B_DATA"]; // no ver1/ver2 split in Zmm
mfile["ZMM_2016Bv1_MC"]   = mfile["ZMM_2016Bv2_MC"]   = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016C_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016C_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016C_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016D_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016D_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016D_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016E_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016E_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016E_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016F_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016Fhipm_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016F_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016Fpost_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016Fnohipm_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016Fpost_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016G_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016G_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016G_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016H_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016H_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016H_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2016BCDEFGH_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2016BCDEFGH_Zmm_V7_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2016BCDEFGH_MC"]    = mfile["ZMM_RUN2_MC"];
//
mfile["ZMM_2017B_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017B_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2017B_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2017C_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017C_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2017C_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2017D_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017D_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2017D_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2017E_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017E_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2017E_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2017F_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2017F_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2017F_MC"]    = mfile["ZMM_RUN2_MC"];
//
mfile["ZMM_2018A_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018A_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2018A_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2018B_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018B_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2018B_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2018C_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018C_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2018C_MC"]    = mfile["ZMM_RUN2_MC"];
mfile["ZMM_2018D_DATA"] = "rootfiles/Run2/Zmm_v116/jme_Zj_2018D_Zmm_V5_v116_DATA.root"; // PLACEHOLDER
mfile["ZMM_2018D_MC"]    = mfile["ZMM_RUN2_MC"];


// Wqq files (Wqq_v1)
// Assumption: "Muon_*" = data, named by MC campaign + era. NanoV15.
// First round: Run2 data vs Run3 Summer24 MC (same file as WQQ_RUN3FLAVOR_MC).
mfile["WQQ_RUN2_MC"] = mfile["WQQ_2024_nib_MC"]; // Summer24 TTtoLNu2Q, carries JER2024 SF
mfile["WQQ_2016HIPM_MC"] = mfile["WQQ_RUN2_MC"];
mfile["WQQ_2016_MC"]     = mfile["WQQ_RUN2_MC"];
mfile["WQQ_2017_MC"]     = mfile["WQQ_RUN2_MC"];
mfile["WQQ_2018_MC"]     = mfile["WQQ_RUN2_MC"];
// No Run2 TTtoLNu2Q MC delivered; guessed paths kept for the later round:
//mfile["WQQ_2016HIPM_MC"] = "rootfiles/Run2/Wqq_v1/Summer20UL16preVFP_TTtoLNu2Q_NanoV15.root"; // PLACEHOLDER
//mfile["WQQ_2016_MC"]     = "rootfiles/Run2/Wqq_v1/Summer20UL16postVFP_TTtoLNu2Q_NanoV15.root"; // PLACEHOLDER
//mfile["WQQ_2017_MC"]     = "rootfiles/Run2/Wqq_v1/Summer20UL17_TTtoLNu2Q_NanoV15.root"; // PLACEHOLDER
//mfile["WQQ_2018_MC"]     = "rootfiles/Run2/Wqq_v1/Summer20UL18_TTtoLNu2Q_NanoV15.root"; // PLACEHOLDER
//
// preVFP data missing entirely, to be requested from the colleagues:
mfile["WQQ_2016BCDEF_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16BCDEFpreVFP_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016BCDEF_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016FGH_DATA"]   = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16FGHpostVFP_NanoV15.root";
mfile["WQQ_2016FGH_MC"]     = mfile["WQQ_2016_MC"];
mfile["WQQ_2016BCDEFGH_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16BCDEFGH_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016BCDEFGH_MC"]   = mfile["WQQ_RUN2_MC"]; // Summer24, as for Run3
//
mfile["WQQ_2017B_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17B_NanoV15.root";
mfile["WQQ_2017B_MC"]   = mfile["WQQ_2017_MC"];
mfile["WQQ_2017C_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17C_NanoV15.root";
mfile["WQQ_2017C_MC"]   = mfile["WQQ_2017_MC"];
mfile["WQQ_2017D_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17D_NanoV15.root";
mfile["WQQ_2017D_MC"]   = mfile["WQQ_2017_MC"];
mfile["WQQ_2017E_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17E_NanoV15.root";
mfile["WQQ_2017E_MC"]   = mfile["WQQ_2017_MC"];
mfile["WQQ_2017F_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17F_NanoV15.root";
mfile["WQQ_2017F_MC"]   = mfile["WQQ_2017_MC"];
mfile["WQQ_2017BCDEF_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL17BCDEF_NanoV15.root"; // PLACEHOLDER (hadd B-F)
mfile["WQQ_2017BCDEF_MC"]   = mfile["WQQ_2017_MC"];
//
mfile["WQQ_2018A_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL18A_NanoV15.root";
mfile["WQQ_2018A_MC"]   = mfile["WQQ_2018_MC"];
mfile["WQQ_2018B_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL18B_NanoV15.root";
mfile["WQQ_2018B_MC"]   = mfile["WQQ_2018_MC"];
mfile["WQQ_2018C_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL18C_NanoV15.root";
mfile["WQQ_2018C_MC"]   = mfile["WQQ_2018_MC"];
mfile["WQQ_2018D_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL18D_NanoV15.root";
mfile["WQQ_2018D_MC"]   = mfile["WQQ_2018_MC"];
mfile["WQQ_2018ABCD_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL18ABCD_NanoV15.root"; // PLACEHOLDER (hadd A-D)
mfile["WQQ_2018ABCD_MC"]   = mfile["WQQ_2018_MC"];
//
// Per-era 2016 Wqq files: requested from the colleagues (cannot be made by
// hadd from Muon_Summer20UL16FGHpostVFP). Paths below are guesses.
mfile["WQQ_2016B_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16B_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016B_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016Bv1_DATA"] = mfile["WQQ_2016Bv2_DATA"] = mfile["WQQ_2016B_DATA"]; // no ver1/ver2 split in Wqq
mfile["WQQ_2016Bv1_MC"]   = mfile["WQQ_2016Bv2_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016C_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16C_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016C_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016D_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16D_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016D_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016E_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16E_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016E_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016F_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16FpreVFP_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016F_MC"]   = mfile["WQQ_2016HIPM_MC"];
mfile["WQQ_2016Fpost_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16FpostVFP_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016Fpost_MC"]   = mfile["WQQ_2016_MC"];
mfile["WQQ_2016G_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16G_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016G_MC"]   = mfile["WQQ_2016_MC"];
mfile["WQQ_2016H_DATA"] = "rootfiles/Run2/Wqq_v1/Muon_Summer20UL16H_NanoV15.root"; // PLACEHOLDER
mfile["WQQ_2016H_MC"]   = mfile["WQQ_2016_MC"];


//#endif
