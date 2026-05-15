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

// ERA-level:
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

// EXTRA semi-year levels:
mlum["2024CDE_nib"] = "26.7 fb^{-1}";
mlum["2024FGHI_nib"] = "83.3 fb^{-1}";
mlum["2024_nib"] = "110 fb^{-1}";

// YEAR-level:
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.4 fb^{-1}";
mlum["2024"] = "110 fb^{-1}";
mlum["2025"] = "109 fb^{-1}";
mlum["2026"] = "17.1 fb^{-1}";

// TOTAL-level:
mlum["Run3"] = "283 fb^{-1}";

mlum["24to26C"] = "236 fb^{-1}";

mlum["2025C0"] = "20.8 fb^{-1}";
mlum["2025CT"] = "TrkRadDam<<20.8 fb^{-1}";

mlum["2026A"] = "0.64 fb^{-1}";
mlum["2026B"] = "15.3 fb^{-1}";
mlum["2026C"] = "1.18 fb^{-1}";

// Jet spike test
mlum["2026BJS"] = "15.3 fb^{-1}";
mlum["2026BNS"] = "15.3 fb^{-1}";
mlum["2026CJS"] = "1.18 fb^{-1}";
mlum["2026CNS"] = "1.18 fb^{-1}";


////////////////////
// File listings  //
////////////////////

//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_v134.root";
//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MG_v150.root";
mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root"; // no JERSF
mfile["JERC_Summer24MG_MC_NOJERSF"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root"; // no JERSF
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
mfile["JET_2024_MC"] = "rootfiles/Prompt/Jet_v161/jmenano_mc_out_Summer24MG_JMENANO_2024CDEFGHI_JERSF_v161.root";
//
mfile["JET_2024_nib_MC"] = mfile["JET_2024_MC"];
mfile["JET_2024_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024CDEFGHI_JME_v165.root";
mfile["JET_2024_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024CDEFGHI_JME_v165.root";
mfile["JET_2024FGHI_nib_MC"] = mfile["JET_2024_MC"];
mfile["JET_2024FGHI_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024FGHI_JME_v165.root";
mfile["JET_2024FGHI_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024FGHI_JME_v165.root";
mfile["JET_2024CDE_nib_MC"]  = mfile["JET_2024_MC"];
mfile["JET_2024CDE_nib_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024CDE_JME_v165.root";
mfile["JET_2024CDE_nib_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024CDE_JME_v165.root";
//
mfile["JET_2024C_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024C_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024C_Rp_JME_v165.root";
mfile["JET_2024C_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024C_Rp_JME_v165.root";
//
mfile["JET_2024D_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024D_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024D_Rp_JME_v165.root";
mfile["JET_2024D_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024D_Rp_JME_v165.root";
//
mfile["JET_2024E_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024E_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024E_Rp_JME_v165.root";
mfile["JET_2024E_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024E_Rp_JME_v165.root";
//
mfile["JET_2024F_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024F_nib1_JME_v165.root";
mfile["JET_2024F_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024F_nib1_JME_v165.root";
mfile["JET_2024F_nib2_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib2_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024F_nib2_JME_v165.root";
mfile["JET_2024F_nib2_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024F_nib2_JME_v165.root";
mfile["JET_2024F_nib3_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024F_nib3_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024F_nib3_JME_v165.root";
mfile["JET_2024F_nib3_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024F_nib3_JME_v165.root";
//
mfile["JET_2024G_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024G_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024G_nib1_JME_v165.root";
mfile["JET_2024G_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024G_nib1_JME_v165.root";
mfile["JET_2024G_nib2_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024G_nib2_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024G_nib2_JME_v165.root";
mfile["JET_2024G_nib2_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024G_nib2_JME_v165.root";
//
mfile["JET_2024H_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024H_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024H_nib1_JME_v159.root"; // was v163_v4
mfile["JET_2024H_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024H_nib1_JME_v165.root";
//
mfile["JET_2024I_nib1_MC"]     = mfile["JET_2024_MC"];
mfile["JET_2024I_nib1_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2024I_nib1_JME_v165.root"; // was v163_v4
mfile["JET_2024I_nib1_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2024I_nib1_JME_v165.root";

//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root"; // no JER SF
//mfile["JET_2025_MC"] = "rootfiles/Prompt2025/Jet_v146/jmenano_mc_out_Winter25MG_v146.root"; // no JER SF
mfile["JET_2025_MC"] = "rootfiles/Prompt/Jet_v161/jmenano_mc_out_Summer24MG_JMENANO_2025CDEFG_JERSF_v161.root";
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v143/jmenano_mc_out_Winter25MG_JRSF2025CDE_v143.root"; // with JER SF
// 2025/Jet_v153 -> Jet_v159 -> v166
mfile["JET_2025C_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025C_JME_v166.root";
mfile["JET_2025C_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025C_JME_v166.root";
mfile["JET_2025C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025D_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025D_JME_v166.root";
mfile["JET_2025D_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025D_JME_v166.root";
mfile["JET_2025D_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025E_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025E_JME_v166.root";
mfile["JET_2025E_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025E_JME_v166.root";
mfile["JET_2025E_MC"]       = mfile["JET_2025_MC"];
// v146->v147 for more F
mfile["JET_2025F_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025F_JME_v166.root";
mfile["JET_2025F_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025F_JME_v166.root";
mfile["JET_2025F_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025G_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025G_JME_v166.root";
mfile["JET_2025G_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025G_JME_v166.root"; //placeholder
mfile["JET_2025G_MC"]       = mfile["JET_2025_MC"];
// v141 CDE, v143 CDEF, v146 CDEF as v147 CDEFG placeholder
mfile["JET_2025CDEFG_DATA_OUT"] = "rootfiles/Prompt/Jet_v166/jmenano_data_out_2025CDEFG_JME_v166.root";
mfile["JET_2025CDEFG_DATA_CMB"] = "rootfiles/Prompt/Jet_v166/jmenano_data_cmb_2025CDEFG_JME_v166.root";
mfile["JET_2025CDEFG_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025JER_DATA_OUT"] = mfile["JET_2025CDEFG_DATA_OUT"];
mfile["JET_2025JER_DATA_CMB"] = mfile["JET_2025CDEFG_DATA_CMB"];
mfile["JET_2025JER_MC"]       = "rootfiles/Prompt2026/Jet_v158/jmenano_mc_out_Summer24MC_Flat_JMENANO_JERSF_v158.root";
//
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
mfile["WQQ_2024_nib_MC"] = "rootfiles/Prompt2024/Wqq_e2/Summer24_TTtoLNu2Q_V9M_24V10MCSF.root";
mfile["WQQ_2024_nib_DATA"]   = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024CDEFGHI_ReReco_V9M.root";
//
mfile["WQQ_2024FGHI_nib_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024FGHI_nib_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024FGHI_Prompt_V9M.root";
mfile["WQQ_2024CDE_nib_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024CDE_nib_DATA"]  = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024CDE_ReReco_V9M.root";
//
mfile["WQQ_2024C_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024C_nib1_DATA"]  = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024C_ReReco_V9M.root";
mfile["WQQ_2024D_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024D_nib1_DATA"]  = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024D_ReReco_V9M.root";
mfile["WQQ_2024E_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024E_nib1_DATA"]  = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024E_ReReco_V9M.root";
//
mfile["WQQ_2024F_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib1_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024F_nib1_Prompt_V9M.root";
mfile["WQQ_2024F_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib2_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024F_nib2_Prompt_V9M.root";
mfile["WQQ_2024F_nib3_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib3_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024F_nib3_Prompt_V9M.root";
mfile["WQQ_2024G_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib1_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024G_nib1__V9M.root";
mfile["WQQ_2024G_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib2_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024G_nib2_Prompt_V9M.root";
mfile["WQQ_2024H_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024H_nib1_DATA"] = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2024H_nib1_Prompt_V9M.root";
mfile["WQQ_2024I_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024I_nib1_DATA"] = "rootfiles/Prompt2024/Wqq_e2/Muon_Run2024I_nib1_Prompt_V9M.root";
//
//mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root"; // Summer24
mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_e2/Summer24_TTtoLNu2Q_V9M_25V3MCSF.root"; // Summer24
mfile["WQQ_2025C_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025C_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025C_Prompt_V3M.root";
mfile["WQQ_2025D_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025D_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025D_Prompt_V3M.root";
mfile["WQQ_2025E_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025E_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025E_Prompt_V3M.root";
mfile["WQQ_2025F_MC"]     = mfile["WQQ_2025_MC"];
//mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025Fv1v2_Prompt_V3M.root"; // was 2025F
mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025F_Prompt_V3M.root";
//
mfile["WQQ_2025G_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025G_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025G_Prompt_V3M.root";
//
mfile["WQQ_2025CDEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025CDEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_e2/Muon_Run2025CDEFG_Prompt_V3M.root";
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
//mfile["ZMM_2024_DATAMC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root";
//mfile["ZMM_2024_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//mfile["ZMM_2024_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_DATA"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024_Zmm_V9M_v113.root";
mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2024Smearing.root";
//mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DYTT_Zmm_V9M_v113.root"; // includes JER SF?
//
//mfile["ZMM_2024FGHI_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
//mfile["ZMM_2024FGHI_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024FGHI_nib_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024FGHI_Zmm_V9M_v113.root"; // FGHI_nib broken for L3Res due to hadd, was used for L2Res
mfile["ZMM_2024FGHI_nib_MC"]   = mfile["ZMM_2024_nib_MC"];
//
//mfile["ZMM_2024CDE_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
//mfile["ZMM_2024CDE_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024CDE_nib_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024CDE_Zmm_V9M_v113.root"; // CDE_nib broken for L3Res due to hadd, was used for L2Res
mfile["ZMM_2024CDE_nib_MC"]   = mfile["ZMM_2024_nib_MC"];
//
//mfile["ZMM_2024C_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024C_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024C_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024D_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024D_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024D_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024D_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024D_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024E_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024E_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024E_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024E_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024E_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024F_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024F_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024F_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
//mfile["ZMM_2024F_nib2_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib2_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024F_nib2_Zmm_V9M_v113.root";
mfile["ZMM_2024F_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024F_nib3_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root"; // bug 2024-04-22, C instead of F_nib3!
mfile["ZMM_2024F_nib3_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024F_nib3_Zmm_V9M_v113.root";
mfile["ZMM_2024F_nib3_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024G_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root"; // bug, C instead of G_nib1
mfile["ZMM_2024G_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024G_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024G_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024G_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root"; // bug, C instead of G_nib2
mfile["ZMM_2024G_nib2_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024G_nib2_Zmm_V9M_v113.root";
mfile["ZMM_2024G_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024H_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root"; // bug, C instead of H
mfile["ZMM_2024H_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024H_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024H_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2024I_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root";
//mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root"; // bug, C instead of I
mfile["ZMM_2024I_nib1_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024I_nib1_Zmm_V9M_v113.root";
mfile["ZMM_2024I_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_Summer24_MC"]   = "rootfiles/Prompt2025/Zmm_v109/jme_Zj_2025_Zmm_V3M_v109.root";
//mfile["ZMM_Winter25_MC"]   = "rootfiles/Prompt2025/Zmm_v103/jme_Zj_2025_Zmm_v103_nomu.root";
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/jme_Zj_2024DY/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing.root"; // Summer24 MC JEC + 2025 JER SF
// v113: JMENANO+Summer24MC JEC for MC+Winter25MC JEC for data+JER SF
mfile["ZMM_2025MC"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2024DY_Zmm_2024V2_v112_nomu_2025Smearing.root"; // Summer24 MC JEC + 2025 JER SF
//
mfile["ZMM_2025C_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025C_Zmm_V3M_v113_nomu_nAOD.root";
mfile["ZMM_2025C_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025D_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025D_Zmm_V3M_v113_nomu_nAOD.root";
mfile["ZMM_2025D_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025E_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025E_Zmm_V3M_v113_nomu_nAOD.root";
mfile["ZMM_2025E_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025F_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025F_Zmm_V3M_v113_nomu_nAOD.root";
mfile["ZMM_2025F_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025G_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025G_Zmm_V3M_v113_nomu_nAOD.root";
mfile["ZMM_2025G_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2025CDEFG_DATA"]   = "rootfiles/Prompt/Zmm_v113/jme_Zj_2025_Zmm_V3M_v113_nomu_nAOD.root";
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
mfile["ZMM_2026B_DATA"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2026B_Zmm_V0M_v113_nomu.root";
mfile["ZMM_2026B_MC"]     = mfile["ZMM_2025MC"];
mfile["ZMM_2026C_DATA"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2026C_Zmm_V0M_v113_nomu_16042026.root";
mfile["ZMM_2026C_MC"]     = mfile["ZMM_2025MC"];
//mfile["ZMM_2026D_DATA"] = "rootfiles/Prompt/Zmm_v113/jme_Zj_2026D_Zmm_V0M_v113_nomu_07052026.root";
mfile["ZMM_2026D_DATA"] = "rootfiles/Prompt/Zmm_v114/jme_Zj_2026D_Zmm_V0M_v114_nomu_10052026.root";
mfile["ZMM_2026D_MC"]     = mfile["ZMM_2025MC"];



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
mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024CDEFGHI-rereco_w73.root";
//mfile["Gam_2024_MC"] = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root";
mfile["GAM_2024_nib_MC"]= //"rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
  "rootfiles/Prompt2024/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
  //"rootfiles/Prompt/Gam_w80/GamHistosFill_mc_summer2024P8-jmenano_no-pu_JERSFX_w80.root"; // wrong bin weights, breaks L2Res
//"rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
//
//mfile["GAM_2024FGHI_nib_DATA"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024FGHI_V9M_w56.root";
mfile["GAM_2024FGHI_nib_DATA"] = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024FGHI-nib_w73.root";
mfile["GAM_2024FGHI_nib_MC"]   = mfile["GAM_2024_nib_MC"];
//mfile["GAM_2024CDE_nib_DATA"]  = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024CDE_V9M_w56.root";
mfile["GAM_2024CDE_nib_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024CDE-rereco_w73.root";
mfile["GAM_2024CDE_nib_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024C_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024C-rereco_w73.root";
mfile["GAM_2024C_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024D_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024D-rereco_w73.root";
mfile["GAM_2024D_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024E_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024E-rereco_w73.root";
mfile["GAM_2024E_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Fnib1_w73.root";
mfile["GAM_2024F_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib2_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Fnib2_w73.root";
mfile["GAM_2024F_nib2_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024F_nib3_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Fnib3_w73.root";
mfile["GAM_2024F_nib3_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024G_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Gnib1_w73.root";
mfile["GAM_2024G_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024G_nib2_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Gnib2_w73.root";
mfile["GAM_2024G_nib2_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024H_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Hnib1_w73.root";
mfile["GAM_2024H_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

mfile["GAM_2024I_nib1_DATA"]  = "rootfiles/Prompt2024/Gam_w73/GamHistosFill_data_2024Inib1_w73.root";
mfile["GAM_2024I_nib1_MC"]    = mfile["GAM_2024_nib_MC"];

//mfile["GAM_2025_MC"]       = "rootfiles/Prompt2025/Gam_w54/GamHistosFill_mc_winter2025P8_no-pu_w54.root";
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC, withPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
mfile["GAM_2025_MC"]    =
  //"rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
  mfile["GAM_2024_nib_MC"];
  //"rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
//mfile["GAM_Summer24_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
//mfile["GAM_Winter25_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_winter2025P8_no-pu_w73.root"; // Winter24 MC, noPU
  // "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_winter2025P8_no-pu_w65.root"; // Winter24 MC, noPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024QCD_no-pu_w65.root"; // Summer24 MC QCD, noPU, TEST ONLY!!
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosMix_mc_summer2024P8_Summer2024QCD_no-pu_w65.root"; // Summer24 MC, noPU
//mfile["GAM_2025_MIX"]      = "rootfiles/Prompt2025/Gam_w54/GamHistosMix_mc_winter2025P8_Winter2025QCD_no-pu_w54.root";
mfile["GAM_2025_MIX"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_DATA"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025C_w73.root";
mfile["GAM_2025C_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025D_DATA"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025D_w73.root";
mfile["GAM_2025D_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025D_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025E_DATA"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025E_w73.root";
mfile["GAM_2025E_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025E_MIX"]     = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025F_DATA"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025F_w73.root";
mfile["GAM_2025F_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025F_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025G_DATA"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025G_w73.root";
mfile["GAM_2025G_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025G_MIX"]     = mfile["GAM_2025_MIX"];
// w60: CDE, w62: CDEF, w64: CDEFG, w65: CDEFG (more G), w67: all (B)CDEFG
mfile["GAM_2025CDEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025CDEFG_w73.root";
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
mfile["GAM_2025DEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025DEFG_w73.root";
mfile["GAM_2025DEFG_MC"]    = mfile["GAM_2025_MC"];
mfile["GAM_2025DEFG_MIX"]   = mfile["GAM_2025_MIX"];
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

// 2026B first Prompt2026 JECs with new HB+HE+HF+PFHC+Winter26 MC JEC
mfile["JET_2026B_DATA_OUT"] = "rootfiles/Prompt/Jet_v163_v2/jmenano_data_out_2026B_JME_v163_v2.root";
mfile["JET_2026B_DATA_CMB"] = "rootfiles/Prompt/Jet_v163_v2/jmenano_data_cmb_2026B_JME_v163_v2.root";
mfile["JET_2026B_MC"]       = mfile["JET_2025_MC"];
//mfile["JET_2026B_MC"]       = mfile["JERC_Summer24MG_MC_NOJERSF"];
// v2: old trig mix + high PU contamination, v4: both fixed
mfile["JET_2026C_DATA_OUT"] = "rootfiles/Prompt/Jet_v163_v4/jmenano_data_out_2026C_JME_v163_v4.root";
mfile["JET_2026C_DATA_CMB"] = "rootfiles/Prompt/Jet_v163_v4/jmenano_data_cmb_2026C_JME_v163_v4.root";
mfile["JET_2026C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2026D_DATA_OUT"] = "rootfiles/Prompt/Jet_v165/jmenano_data_out_2026D_JME_v165.root";
mfile["JET_2026D_DATA_CMB"] = "rootfiles/Prompt/Jet_v165/jmenano_data_cmb_2026D_JME_v165.root";
mfile["JET_2026D_MC"]       = mfile["JET_2025_MC"];

mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt/Gam_w79/GamHistosFill_data_2026B_w79.root"; // with L2L3Res, for V0M closure
mfile["GAM_2026B_MC"]     = mfile["GAM_2025_MC"];
mfile["GAM_2026B_MIX"]    = mfile["GAM_2025_MIX"];
//
mfile["GAM_2026C_DATA"]   = "rootfiles/Prompt/Gam_w79/GamHistosFill_data_2026C_w79.root"; // with L2L3Res, new low PU data
mfile["GAM_2026C_MC"]     = mfile["GAM_2025_MC"];
mfile["GAM_2026C_MIX"]    = mfile["GAM_2025_MIX"];
//
mfile["GAM_2026D_DATA"]   = "rootfiles/Prompt/Gam_w81/GamHistosFill_data_2026D_w81.root";
mfile["GAM_2026D_MC"]     = mfile["GAM_2025_MC"];
mfile["GAM_2026D_MIX"]    = mfile["GAM_2025_MIX"];

mfile["WQQ_2026B_DATA"]   = "rootfiles/Prompt/Wqq_e4/Muon_Run2026B_Prompt_V0M_GoldenJSON_e4.root"; // v2: no L2L3Res, e2,e4: L3L2Res
//mfile["WQQ_2026B_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2026B_MC"]     = "rootfiles/Prompt2026/Wqq_e2/Summer24_TTtoLNu2Q_V9M_26BV0MCSF.root"; // Summer24 + JER SF
//
mfile["WQQ_2026C_DATA"]   = "rootfiles/Prompt/Wqq_e4/Muon_Run2026C_Prompt_V0M_CombinedJSON_e4.root"; // v2: no L2L3Res, e2, e3, e4: L3L2Res
mfile["WQQ_2026C_MC"]     = mfile["WQQ_2026B_MC"];
//
mfile["WQQ_2026D_DATA"]   = "rootfiles/Prompt/Wqq_e5/Muon_Run2026D_Prompt_V0M_CombinedJSON_e5.root";
mfile["WQQ_2026D_MC"]     = mfile["WQQ_2026B_MC"];

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

//#endif
