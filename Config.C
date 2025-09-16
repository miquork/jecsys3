// Purpose: common configurations for L2Res.C, L3Res.C/globalFit, JERSF.C
#ifndef __CONFIG_C__
#define __CONFIG_C__

#include <map>
#include <string>
std::map<std::string, std::string> mlum;
std::map<std::string, std::string> mfile;

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

/*
mlum["2024B_nib1"] = "0.13 fb^{-1}";
mlum["2024C_nib1"] = "7.5 fb^{-1}";
mlum["2024D_nib1"] = "8.3 fb^{-1}";
//mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024Ev1_nib1"] = "X.X fb^{-1}";
mlum["2024Ev2_nib1"] = "X.X fb^{-1}";
//mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024F_nib1"] = "0.88 fb^{-1}";
mlum["2024F_nib2"] = "X.X fb^{-1}";
mlum["2024F_nib3"] = "X.X fb^{-1}";
//mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024G_nib1"] = "X.X fb^{-1}";
mlum["2024G_nib2"] = "X.X fb^{-1}";
mlum["2024H_nib1"] = "5.4 fb^{-1}";
mlum["2024I_nib1"] = "11.6 fb^{-1}";
*/
// rootfiles/brilcalc/sumfibs.py

// NIB-level:
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
mlum["2024B_nib1"] = "0.130 fb^{-1}";
mlum["2024C_nib1"] = "7.2 fb^{-1}";
mlum["2024D_nib1"] = "8.0 fb^{-1}";
mlum["2024Ev1_nib1"] = "6.3 fb^{-1}";
mlum["2024Ev2_nib1"] = "5.0 fb^{-1}";
mlum["2024E_nib1"] = "11.3 fb^{-1}";
mlum["2024F_nib1"] = "0.88 fb^{-1}";
mlum["2024F_nib2"] = "14.4 fb^{-1}";
mlum["2024F_nib3"] = "12.5 fb^{-1}";
mlum["2024G_nib1"] = "16.9 fb^{-1}";
mlum["2024G_nib2"] = "20.9 fb^{-1}";
mlum["2024H_nib1"] = "5.4 fb^{-1}";
mlum["2024I_nib1"] = "11.5 fb^{-1}";

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
mlum["2024B"] = "0.130 fb^{-1}";
mlum["2024C"] = "7.2 fb^{-1}";
mlum["2024D"] = "8.0 fb^{-1}";
mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024H"] = "5.4 fb^{-1}";
mlum["2024I"] = "11.5 fb^{-1}";

mlum["2025C"] = "X.X fb^{-1}";

// YEAR-level:
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.4 fb^{-1}";
mlum["2024"] = "109 fb^{-1}";

mlum["2025"] = "X.X fb^{-1}";
//Cv1: 12.463961608
//Cv2: 8.320040199
mlum["2025C"] = "20.8 fb^{-1}";
//D: 24.040148153
mlum["2025D"] = "24.0 fb^{-1}";
//E: 13.896356116 (DIALS JSON)
mlum["2025E"] = "13.9 fb^{-1}"; // DIALS JSON
// CDE: 58.720506 (hybrid)
mlum["2025CDE"] = "58.7 fb^{-1}"; // HYBRID JSON

// File listings
// Multijet files
mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root";
mfile["JET_2025C_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025C_JME_v141.root";
mfile["JET_2025C_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_cmb_2025C_JME_v141.root";
mfile["JET_2025C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025D_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025D_JME_v141.root";
mfile["JET_2025D_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_cmb_2025D_JME_v141.root";
mfile["JET_2025D_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025E_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025E_JME_v141.root";
mfile["JET_2025E_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_cmb_2025E_JME_v141.root";
mfile["JET_2025E_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025CDE_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_out_2025CDE_JME_v141.root";
mfile["JET_2025CDE_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v141/jmenano_data_cmb_2025CDE_JME_v141.root";
mfile["JET_2025CDE_MC"]       = mfile["JET_2025_MC"];

// Wqq files
mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_v1m/Summer24_TTtoLNu2Q.root";
mfile["WQQ_2025C_DATA"]   = "rootfiles/Prompt2025/Wqq_v1m/Muon_Run2025C_Prompt_V1M.root";
mfile["WQQ_2025D_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025D_DATA"]   = "rootfiles/Prompt2025/Wqq_v1m/Muon_Run2025D_Prompt_V1M.root";
mfile["WQQ_2025E_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025E_DATA"]   = "rootfiles/Prompt2025/Wqq_v1m/Muon_Run2025E_Prompt_V1M.root";

// Z+jet files
mfile["ZMM_2025C_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v100/jme_Zj_2025C_Zmm_v100.root"; // Missing
mfile["ZMM_2025D_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v100/jme_Zj_2025D_Zmm_v100.root";
mfile["ZMM_2025E_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v100/jme_Zj_2025E_Zmm_v100.root";
mfile["ZMM_2025CDE_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v100/jme_Zj_2025D_Zmm_v100.root"; // Placeholder

// Photon+jet files
mfile["GAM_2025_MC"]       = "rootfiles/Prompt2025/Gam_w54/GamHistosFill_mc_winter2025P8_no-pu_w54.root";
mfile["GAM_2025_MIX"]      = "rootfiles/Prompt2025/Gam_w54/GamHistosMix_mc_winter2025P8_Winter2025QCD_no-pu_w54.root";
mfile["GAM_2025C_DATA"]    = "rootfiles/Prompt2025/Gam_w60/GamHistosFill_data_2025C_w60.root";
mfile["GAM_2025C_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025D_DATA"]    = "rootfiles/Prompt2025/Gam_w60/GamHistosFill_data_2025D_w60.root";
mfile["GAM_2025D_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025D_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025E_DATA"]    = "rootfiles/Prompt2025/Gam_w60/GamHistosFill_data_2025E_w60.root";
mfile["GAM_2025E_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025E_MIX"]     = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025CDE_DATA"]  = "rootfiles/Prompt2025/Gam_w60/GamHistosFill_data_2025CDEF_w60.root"; // NB: also F
//mfile["GAM_2025CDE_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025CDE_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC
mfile["GAM_2025CDE_MIX"]   = mfile["GAM_2025_MIX"];



#endif
