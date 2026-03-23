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

// YEAR-level:
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.4 fb^{-1}";
mlum["2024"] = "110 fb^{-1}";
mlum["2025"] = "109 fb^{-1}";

// TOTAL-level:
mlum["Run3"] = "283 fb^{-1}";

mlum["2025C0"] = "20.8 fb^{-1}";
mlum["2025CT"] = "TrkRadDam<<20.8 fb^{-1}";

mlum["2026A"] = "0.4 fb^{-1}";
mlum["2026B"] = "1? fb^{-1}";

////////////////////
// File listings  //
////////////////////

//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_v134.root";
//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MG_v150.root";
mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MG_v155.root";
//mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MC_Flat_v150.root";
mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MC_Flat_v155.root";
mfile["JERC_Summer24MC_Flat22_Herwig_MC"] = "rootfiles/Prompt2024/Jet_v155/jmenano_mc_out_Summer24MC_Flat22_Herwig_v155.root";
mfile["JERC_Summer24MC_FlatBase_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_Base_v151_v3.root";
mfile["JERC_Summer24MC_FlatNoDeepCore_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_NoDeepCore_v151_v3.root";
mfile["JERC_Winter25MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Winter25MG_v150.root";
//
mfile["JERC_Summer24MC_Flat22_NoPU_MC"] = "rootfiles/Prompt2025/Jet_v155/jmenano_mc_out_Summer24MC_Flat22_NoPU_v155.root";
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
mfile["JET_2024_nib_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024_nib_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024CDEFGHI_JME_v155.root";
mfile["JET_2024_nib_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024CDEFGHI_JME_v155.root";
mfile["JET_2024FGHI_nib_MC"] = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024FGHI_nib_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024FGHI_JME_v155.root";
mfile["JET_2024FGHI_nib_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024FGHI_JME_v155.root";
mfile["JET_2024CDE_nib_MC"]  = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024_nib_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024CDE_JME_v155.root";
mfile["JET_2024CDE_nib_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024CDE_JME_v155.root";
//
mfile["JET_2024C_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024C_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024C_nib1_JME_v155.root";
mfile["JET_2024C_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024C_nib1_JME_v155.root";
//
mfile["JET_2024D_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024D_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024D_nib1_JME_v155.root";
mfile["JET_2024D_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024D_nib1_JME_v155.root";
//
mfile["JET_2024E_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024E_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024E_nib1_JME_v155.root";
mfile["JET_2024E_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024E_nib1_JME_v155.root";
//
mfile["JET_2024F_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024F_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib1_JME_v155.root";
mfile["JET_2024F_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib1_JME_v155.root";
mfile["JET_2024F_nib2_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024F_nib2_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib2_JME_v155.root";
mfile["JET_2024F_nib2_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib2_JME_v155.root";
mfile["JET_2024F_nib3_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024F_nib3_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib3_JME_v155.root";
mfile["JET_2024F_nib3_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024F_nib3_JME_v155.root";
//
mfile["JET_2024G_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024G_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024G_nib1_JME_v155.root";
mfile["JET_2024G_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024G_nib1_JME_v155.root";
mfile["JET_2024G_nib2_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024G_nib2_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024G_nib2_JME_v155.root";
mfile["JET_2024G_nib2_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024G_nib2_JME_v155.root";
//
mfile["JET_2024H_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024H_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024H_nib1_JME_v155.root";
mfile["JET_2024H_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024H_nib1_JME_v155.root";
//
mfile["JET_2024I_nib1_MC"]     = mfile["JERC_Summer24MG_MC"];
mfile["JET_2024I_nib1_DATA_OUT"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024I_nib1_JME_v155.root";
mfile["JET_2024I_nib1_DATA_CMB"] = "rootfiles/Prompt2024/Jet_v155/jmenano_data_cmb_2024I_nib1_JME_v155.root";

//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root"; // no JER SF
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v146/jmenano_mc_out_Winter25MG_v146.root"; // no JER SF
mfile["JET_2025_MC"]        = mfile["JERC_Summer24MG_MC"];
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v143/jmenano_mc_out_Winter25MG_JRSF2025CDE_v143.root"; // with JER SF
mfile["JET_2025C_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025C_JME_v153.root";
mfile["JET_2025C_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025C_JME_v153.root";
mfile["JET_2025C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025D_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025D_JME_v153.root";
mfile["JET_2025D_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025D_JME_v153.root";
mfile["JET_2025D_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025E_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025E_JME_v153.root";
mfile["JET_2025E_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025E_JME_v153.root";
mfile["JET_2025E_MC"]       = mfile["JET_2025_MC"];
// v146->v147 for more F
mfile["JET_2025F_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025F_JME_v153.root";
mfile["JET_2025F_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025F_JME_v153.root";
mfile["JET_2025F_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025G_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025G_JME_v153.root";
mfile["JET_2025G_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025G_JME_v153.root"; //placeholder
mfile["JET_2025G_MC"]       = mfile["JET_2025_MC"];
// v141 CDE, v143 CDEF, v146 CDEF as v147 CDEFG placeholder
mfile["JET_2025CDEFG_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025CDEFG_JME_v153.root";
mfile["JET_2025CDEFG_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025CDEFG_JME_v153.root";
mfile["JET_2025CDEFG_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025DEFG_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_out_2025DEFG_JME_v153.root";
mfile["JET_2025DEFG_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v153/jmenano_data_cmb_2025DEFG_JME_v153.root";
mfile["JET_2025DEFG_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025C0_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_out_2025C_JME_v145.root";
mfile["JET_2025C0_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_cmb_2025C_JME_v145.root";
mfile["JET_2025C0_MC"] = "rootfiles/Prompt2025/Jet_v145/jmenano_mc_out_Winter25MG_v145.root";
mfile["JET_2025CT_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_out_2025C_Trk_JME_v145.root";
mfile["JET_2025CT_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v145/jmenano_data_cmb_2025C_Trk_JME_v145.root";
mfile["JET_2025CT_MC"] = mfile["JET_2025C0_MC"]; 


// Wqq files
// v1m: input to Prompt25_V2M
// V2M: closure of Prompt25_V2M, input to Prompt25_V3M
// V3M: closure of Prompt25_V3M
mfile["WQQ_2024_nib_MC"] = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root";
mfile["WQQ_2024_nib_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDEFGHI_ReReco_V9M.root";
//
mfile["WQQ_2024FGHI_nib_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024FGHI_nib_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root";
mfile["WQQ_2024CDE_nib_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024CDE_nib_DATA"]  = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDE_ReReco_V9M.root";
//
mfile["WQQ_2024C_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024C_nib1_DATA"]  = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDE_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024D_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024D_nib1_DATA"]  = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDE_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024E_nib1_MC"]    = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024E_nib1_DATA"]  = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDE_ReReco_V9M.root"; // placeholder
//
mfile["WQQ_2024F_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib1_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024F_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib2_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024F_nib3_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024F_nib3_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024G_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib1_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024G_nib2_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024G_nib2_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024H_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024H_nib1_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
mfile["WQQ_2024I_nib1_MC"]   = mfile["WQQ_2024_nib_MC"];
mfile["WQQ_2024I_nib1_DATA"] = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024FGHI_ReReco_V9M.root"; // placeholder
//
mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root"; // Summer24
mfile["WQQ_2025C_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025C_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025C_Prompt_V3M.root";
mfile["WQQ_2025D_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025D_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025D_Prompt_V3M.root";
mfile["WQQ_2025E_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025E_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025E_Prompt_V3M.root";
mfile["WQQ_2025F_MC"]     = mfile["WQQ_2025_MC"];
//mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025Fv1v2_Prompt_V3M.root"; // was 2025F
mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025F_Prompt_V3M.root"; // was 2025F
//
mfile["WQQ_2025G_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025G_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025G_Prompt_V3M.root"; // partial G
//
mfile["WQQ_2025CDEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025CDEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025CDEFG_Prompt_V3M.root";
//
mfile["WQQ_2025DEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025DEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_V3M/Muon_Run2025DEFG_Prompt_V3M.root";


// Z+jet files
// v100: input to Prompt25_V2M
// v101: closure of Prompt25_V2M, input to Prompt25_V3M
// v102 + _nomu: add eta asymmetry Rochester, revert to Summer24
// v103 : 2025G full, Summer25
// v104: 2025 final golden, recovered runs
// V3M_v105: V3M closure
// V3M_v107: fixed Q,G tagging
// V3M_v110: moved to JMENano (since v108)
//mfile["ZMM_2024_DATAMC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root";
mfile["ZMM_2024_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//
mfile["ZMM_2024FGHI_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024FGHI_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024FGHI_nib_MC"]   = mfile["ZMM_2024_nib_MC"];//"rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024FGHI_nib_Zmm_V9M_v103.root";
//
mfile["ZMM_2024CDE_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024CDE_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
mfile["ZMM_2024CDE_nib_MC"]   = mfile["ZMM_2024_nib_MC"];//"rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024CDE_nib_Zmm_V9M_v103.root";
//
mfile["ZMM_2024C_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024D_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024D_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024D_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024E_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024E_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024E_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024F_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib1_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024F_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib2_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024F_nib3_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024F_nib3_Zmm_V9M_v103.root";
mfile["ZMM_2024F_nib3_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024G_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024G_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024G_nib2_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024G_nib2_Zmm_V9M_v103.root";
mfile["ZMM_2024G_nib2_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024H_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024H_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024H_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

mfile["ZMM_2024I_nib1_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024C_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024I_nib1_Zmm_V9M_v103.root";
mfile["ZMM_2024I_nib1_MC"]   = mfile["ZMM_2024_nib_MC"];

//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024_Zmm_V9M_v100.root"; // Summer24 to replace Winter25, gluontag not broken
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root"; // Summer24 to replace Winter25, gluontag broken
mfile["ZMM_Summer24_MC"]   = "rootfiles/Prompt2025/Zmm_v109/jme_Zj_2025_Zmm_V3M_v109.root"; // Summer24 MC now
//"rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_Winter25_MC"]   = "rootfiles/Prompt2025/Zmm_v103/jme_Zj_2025_Zmm_v103_nomu.root";
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root"; // Summer24 to replace Winter25, gluontag broken (also in data)
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2025/Zmm_v107/jme_Zj_2025_Zmm_V3M_v107_nomu.root"; // Summer24 or Winter25 MC?
mfile["ZMM_2025MC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025_Zmm_V3M_v110_nomu.root"; // Summer24 JMENano (TBC=
//
mfile["ZMM_2025C_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025C_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025C_DATA"] = mfile["ZMM_2025C_DATAMC"]; 
mfile["ZMM_2025C_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025D_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025D_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025D_DATA"] = mfile["ZMM_2025D_DATAMC"]; 
mfile["ZMM_2025D_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025E_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025E_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025E_DATA"] = mfile["ZMM_2025E_DATAMC"]; 
mfile["ZMM_2025E_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025F_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025F_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025F_DATA"] = mfile["ZMM_2025F_DATAMC"]; 
mfile["ZMM_2025F_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025G_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025G_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025G_DATA"] = mfile["ZMM_2025G_DATAMC"]; 
mfile["ZMM_2025G_MC"] = mfile["ZMM_2025MC"];
//mfile["ZMM_2025G_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_JMENano/jme_Zj_2025G_PromptJMENano_v1_Zmm_v103_nomu.root"; // TEST of lower Z pT cut
//mfile["ZMM_2025G_MC"] = "rootfiles/Prompt2025/Zmm_v103_JMENano/jme_Zj_2024I_JMENANOv15_v1_Zmm_v103_nomu.root"; // TEST of lower Z pT cut
// v100: CDE, v102: CDEF (=CDE+F), v102_nomu: CDEFG (=2025)
mfile["ZMM_2025CDEFG_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025CDEFG_DATA"] = mfile["ZMM_2025CDEFG_DATAMC"]; 
mfile["ZMM_2025CDEFG_MC"] = mfile["ZMM_2025MC"];
//
mfile["ZMM_2025DEFG_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v110/jme_Zj_2025DEFG_Zmm_V3M_v110_nomu.root";
mfile["ZMM_2025DEFG_DATA"] = mfile["ZMM_2025DEFG_DATAMC"]; 
mfile["ZMM_2025DEFG_MC"] = mfile["ZMM_2025MC"];
//
mfile["ZMM_2025C0_DATAMC"]   = mfile["ZMM_2025C_DATAMC"];
mfile["ZMM_2025CT_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v101/jme_Zj_2025C_TrkRadDamage_Zmm_v101.root"; // for L2Res
mfile["ZMM_2025CT_DATA"]   = "rootfiles/Prompt2025/Zmm_v101/jme_Zj_2025C_TrkRadDamage_Zmm_v101.root";
mfile["ZMM_2025CT_MC"]   = mfile["ZMM_2025C_DATAMC"];


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
  "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
  //"rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
mfile["GAM_Summer24_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_summer2024P8_no-pu_w73.root"; // Summer24 MC, noPU
mfile["GAM_Winter25_MC"]    = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_mc_winter2025P8_no-pu_w73.root"; // Winter24 MC, noPU
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
mfile["GAM_2025DEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w73/GamHistosFill_data_2025DEFG_w73.root";
mfile["GAM_2025DEFG_MC"]    = mfile["GAM_2025_MC"];
mfile["GAM_2025DEFG_MIX"]   = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025C0_DATA"]    = "rootfiles/Prompt2025/Gam_w62/GamHistosFill_data_2025C_w62.root";
mfile["GAM_2025C0_MC"]   = mfile["GAM_2025_MC"];
mfile["GAM_2025C0_MIX"]   = mfile["GAM_2025_MIX"];
mfile["GAM_2025CT_DATA"]    = "rootfiles/Prompt2025/Gam_w63/GamHistosFill_data_2025C-TrkRadDamage_w63.root";
mfile["GAM_2025CT_MC"]   = mfile["GAM_2025_MC"];
mfile["GAM_2025CT_MIX"]   = mfile["GAM_2025_MIX"];



// 2026A initial studies
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

// 2026B first Prompt2026 JECs with new HB+HE+HF+PFHC+Winter26 MC JEC
//mfile["JET_2026Bnib1_DATA_OUT"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_out_2026Bnib1_JME_v157_v3.root"; // v2 no L2L3Res
//mfile["JET_2026Bnib1_DATA_CMB"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_cmb_2026Bnib1_JME_v157_v3.root";
//mfile["JET_2026Bnib1_MC"]       = mfile["JET_2025_MC"];
//
//mfile["JET_2026Bnib2_DATA_OUT"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_out_2026Bnib2_JME_v157_v3.root"; // v2 no L2L3Res
//mfile["JET_2026Bnib2_DATA_CMB"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_cmb_2026Bnib2_JME_v157_v3.root";
//mfile["JET_2026Bnib2_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2026B_DATA_OUT"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_out_2026B_JME_v157_v4.root"; // v2 no L2L3Res, v3 L2L3Res, v4 more lum
mfile["JET_2026B_DATA_CMB"] = "rootfiles/Prompt2026/Jet_v157/jmenano_data_cmb_2026B_JME_v157_v4.root"; // v2 no L2L3Res, v3 L2L3Res, v4 more lum
mfile["JET_2026B_MC"]       = mfile["JET_2025_MC"];

//mfile["ZMM_2026Bnib1_DATAMC"] = "rootfiles/Prompt2026/Zmm_v110/jme_Zj_2026B_Zmm_L2L3Res_v110_nomu.root";
//mfile["ZMM_2026Bnib1_DATA"]   = mfile["ZMM_2026Bnib1_DATAMC"]; 
//mfile["ZMM_2026Bnib1_MC"]     = mfile["ZMM_2025MC"];
//
//mfile["ZMM_2026Bnib2_DATAMC"] = "rootfiles/Prompt2026/Zmm_v110/jme_Zj_2026B_Zmm_L2L3Res_v110_nomu.root";
//mfile["ZMM_2026Bnib2_DATA"]   = mfile["ZMM_2026Bnib2_DATAMC"]; 
//mfile["ZMM_2026Bnib2_MC"]     = mfile["ZMM_2025MC"];
//
mfile["ZMM_2026B_DATAMC"] = "rootfiles/Prompt2026//Zmm_v110/jme_Zj_2026B_Zmm_v110_nomu_20032026.root"; // Zmm_L2L3Res as option
mfile["ZMM_2026B_DATA"]   = mfile["ZMM_2026B_DATAMC"]; 
mfile["ZMM_2026B_MC"]     = mfile["ZMM_2025MC"];

//mfile["GAM_2026Bnib1_DATA"]   = "rootfiles/Prompt2026/Gam_w75/GamHistosFill_data_2026B_w75_17Mar2026_Bnib1_L2Rel2026MC_L2L3Res2025G.root";
//mfile["GAM_2026Bnib1_MC"]     = mfile["GAM_2025_MC"];
//mfile["GAM_2026Bnib1_MIX"]    = mfile["GAM_2025_MIX"];
//
//mfile["GAM_2026Bnib2_DATA"]   = "rootfiles/Prompt2026/Gam_w75/GamHistosFill_data_2026B_w75_17Mar2026_Bnib2_L2Rel2026MC_L2L3Res2025G.root";
//mfile["GAM_2026Bnib2_MC"]     = mfile["GAM_2025_MC"];
//mfile["GAM_2026Bnib2_MIX"]    = mfile["GAM_2025_MIX"];
//
//mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt2026/Gam_w75/GamHistosFill_data_2026B_w75_17Mar2026_L2Rel2026MC_L2L3Res2025G.root";
//mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt2026/Gam_w75/GamHistosFill_data_2026B_w75_17Mar2026_onlyL2Rel2026MC.root";
mfile["GAM_2026B_DATA"]   = "rootfiles/Prompt2026/Gam_w76/GamHistosFill_data_2026B_w75_onlyL2Rel2026MC.root";
mfile["GAM_2026B_MC"]     = mfile["GAM_2025_MC"];
mfile["GAM_2026B_MIX"]    = mfile["GAM_2025_MIX"];

//mfile["WQQ_2026Bnib1_DATA"]   = "rootfiles/Prompt2026/Wqq_v1/Muon_Run2026Bnib1_Prompt_L2L32025G.root";
//mfile["WQQ_2026Bnib1_MC"]     = mfile["WQQ_2025_MC"];
//
//mfile["WQQ_2026Bnib2_DATA"]   = "rootfiles/Prompt2026/Wqq_v1/Muon_Run2026Bnib2_Prompt_L2L32025G.root";
//mfile["WQQ_2026Bnib2_MC"]     = mfile["WQQ_2025_MC"];
//
//mfile["WQQ_2026B_DATA"]   = "rootfiles/Prompt2026/Wqq_v1/Muon_Run2026B_Prompt.root";//_L2L32025G.root";
mfile["WQQ_2026B_DATA"]   = "rootfiles/Prompt2026/Wqq_v2/Muon_Run2026B_Prompt_CombinedJSON.root";
mfile["WQQ_2026B_MC"]     = mfile["WQQ_2025_MC"];

//#endif
