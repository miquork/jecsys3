// Purpose: common configurations for L2Res.C, L3Res.C/globalFit, JERSF.C
#ifndef __CONFIG_C__
#define __CONFIG_C__

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

////////////////////
// File listings  //
////////////////////

//mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_v134.root";
mfile["JERC_Summer24MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MG_v150.root";
mfile["JERC_Summer24MC_Flat_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Summer24MC_Flat_v150.root";
mfile["JERC_Summer24MC_FlatBase_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_Base_v151_v3.root";
mfile["JERC_Summer24MC_FlatNoDeepCore_MC"] = "rootfiles/Prompt2025/Jet_v151/jmenano_mc_out_Summer24MC_Flat22_NoDeepCore_v151_v3.root";
mfile["JERC_Winter25MG_MC"] = "rootfiles/Prompt2025/Jet_v150/jmenano_mc_out_Winter25MG_v150.root";

// Multijet files
// v141 (MC v128): input to Prompt25_V2M
// v143 (MC v143 JRV2M: closure test of Prompt25_V2M, input to Prompt25_V3M
// 143->146: fix Dijet folder binning and improve triggers (Jet110)
// 146/147->150(_v2): more MC's, latest data, what else?
// v151
// v152: golden JSON, rhovsmu tc.
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v128/jmenano_mc_out_Winter25MG_v128.root"; // no JER SF
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v146/jmenano_mc_out_Winter25MG_v146.root"; // no JER SF
mfile["JET_2025_MC"]        = mfile["JERC_Summer24MG_MC"];
//mfile["JET_2025_MC"]        = "rootfiles/Prompt2025/Jet_v143/jmenano_mc_out_Winter25MG_JRSF2025CDE_v143.root"; // with JER SF
mfile["JET_2025C_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025C_JME_v152.root";
mfile["JET_2025C_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025C_JME_v152.root";
mfile["JET_2025C_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025D_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025D_JME_v152.root";
mfile["JET_2025D_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025D_JME_v152.root";
mfile["JET_2025D_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025E_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025E_JME_v152.root";
mfile["JET_2025E_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025E_JME_v152.root";
mfile["JET_2025E_MC"]       = mfile["JET_2025_MC"];
// v146->v147 for more F
mfile["JET_2025F_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v147/jmenano_data_out_2025F_JME_v147.root";
mfile["JET_2025F_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v147/jmenano_data_cmb_2025F_JME_v147.root";
mfile["JET_2025F_MC"]       = mfile["JET_2025_MC"];
mfile["JET_2025G_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025G_JME_v152.root";
mfile["JET_2025G_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025G_JME_v152.root"; //placeholder
mfile["JET_2025G_MC"]       = mfile["JET_2025_MC"];
// v141 CDE, v143 CDEF, v146 CDEF as v147 CDEFG placeholder
mfile["JET_2025CDEFG_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025CDEFG_JME_v152.root";
mfile["JET_2025CDEFG_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025CDEFG_JME_v152.root";
mfile["JET_2025CDEFG_MC"]       = mfile["JET_2025_MC"];
//
mfile["JET_2025DEFG_DATA_OUT"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_out_2025DEFG_JME_v152.root";
mfile["JET_2025DEFG_DATA_CMB"] = "rootfiles/Prompt2025/Jet_v152/jmenano_data_cmb_2025DEFG_JME_v152.root";
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
mfile["WQQ_2024_nib_MC"]     = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root";
mfile["WQQ_2024_nib_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2024CDEFGHI_ReReco_V9M.root";
//
mfile["WQQ_2025_MC"]     = "rootfiles/Prompt2025/Wqq_V2M/Summer24_TTtoLNu2Q.root"; // Summer24
mfile["WQQ_2025C_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025C_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025C_Prompt_V2M.root";
mfile["WQQ_2025D_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025D_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025D_Prompt_V2M.root";
mfile["WQQ_2025E_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025E_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025E_Prompt_V2M.root";
mfile["WQQ_2025F_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025F_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025Fv1v2_Prompt_V2M.root"; // was 2025F
//
mfile["WQQ_2025G_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025G_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025G_Prompt_V2M.root"; // partial G
//
mfile["WQQ_2025CDEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025CDEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025CDEFG_Prompt_V2M.root"; // was 2025CDEF, but still partial G
//
mfile["WQQ_2025DEFG_MC"]     = mfile["WQQ_2025_MC"];
mfile["WQQ_2025DEFG_DATA"]   = "rootfiles/Prompt2025/Wqq_V2M/Muon_Run2025DEFG_Prompt_V2M.root"; // was 2025CDEF, but still partial G


// Z+jet files
// v100: input to Prompt25_V2M
// v101: closure of Prompt25_V2M, input to Prompt25_V3M
// v102 + _nomu: add eta asymmetry Rochester, revert to Summer24
// v103 : 2025G full, Summer25
// v104: 2025 final golden, recovered runs
//mfile["ZMM_2024_DATAMC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root";
mfile["ZMM_2024_nib_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_DATA"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_2024_nib_MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
//
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024_Zmm_V9M_v100.root"; // Summer24 to replace Winter25, gluontag not broken
//mfile["ZMM_2025MC"]   = "rootfiles/Prompt2024/Zmm_v102/jme_Zj_2024_Zmm_V9M_v102.root"; // Summer24 to replace Winter25, gluontag broken
mfile["ZMM_Summer24_MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root";
mfile["ZMM_Winter25_MC"]   = "rootfiles/Prompt2025/Zmm_v103/jme_Zj_2025_Zmm_v103_nomu.root";
mfile["ZMM_2025MC"]   = "rootfiles/Prompt2025/Zmm_v103_2024_V9M/jme_Zj_2024_Zmm_V9M_v103.root"; // Summer24 to replace Winter25, gluontag broken (also in data)
mfile["ZMM_2025C_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025C_Zmm_v104_nomu.root";
mfile["ZMM_2025C_DATA"] = mfile["ZMM_2025C_DATAMC"]; 
mfile["ZMM_2025C_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025D_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025D_Zmm_v104_nomu.root";
mfile["ZMM_2025D_DATA"] = mfile["ZMM_2025D_DATAMC"]; 
mfile["ZMM_2025D_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025E_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025E_Zmm_v104_nomu.root";
mfile["ZMM_2025E_DATA"] = mfile["ZMM_2025E_DATAMC"]; 
mfile["ZMM_2025E_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025F_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025F_Zmm_v104_nomu.root"; // hadd F1+F2
mfile["ZMM_2025F_DATA"] = mfile["ZMM_2025F_DATAMC"]; 
mfile["ZMM_2025F_MC"] = mfile["ZMM_2025MC"];
mfile["ZMM_2025G_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025G_Zmm_v104_nomu.root";
mfile["ZMM_2025G_DATA"] = mfile["ZMM_2025G_DATAMC"]; 
mfile["ZMM_2025E_MC"] = mfile["ZMM_2025MC"];
// v100: CDE, v102: CDEF (=CDE+F), v102_nomu: CDEFG (=2025)
mfile["ZMM_2025CDEFG_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025_Zmm_v104_nomu.root";
mfile["ZMM_2025CDEFG_DATA"] = mfile["ZMM_2025CDEFG_DATAMC"]; 
mfile["ZMM_2025CDEFG_MC"] = mfile["ZMM_2025MC"];
//
mfile["ZMM_2025DEFG_DATAMC"]   = "rootfiles/Prompt2025/Zmm_v104/jme_Zj_2025DEFG_Zmm_v104_nomu.root";
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
mfile["GAM_2024_nib_DATA"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024_V9M_w56.root";
//mfile["Gam_2024_MC"] = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root";
mfile["GAM_2024_nib_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
//
//mfile["GAM_2025_MC"]       = "rootfiles/Prompt2025/Gam_w54/GamHistosFill_mc_winter2025P8_no-pu_w54.root";
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC, withPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_no-pu_w48.root"; // Summer24 MC, noPU
mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
mfile["GAM_Summer24_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024P8_no-pu_w65.root"; // Summer24 MC, noPU
mfile["GAM_Winter25_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_winter2025P8_no-pu_w65.root"; // Winter24 MC, noPU
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosFill_mc_summer2024QCD_no-pu_w65.root"; // Summer24 MC QCD, noPU, TEST ONLY!!
//mfile["GAM_2025_MC"]    = "rootfiles/Prompt2025/Gam_w65/GamHistosMix_mc_summer2024P8_Summer2024QCD_no-pu_w65.root"; // Summer24 MC, noPU
//mfile["GAM_2025_MIX"]      = "rootfiles/Prompt2025/Gam_w54/GamHistosMix_mc_winter2025P8_Winter2025QCD_no-pu_w54.root";
mfile["GAM_2025_MIX"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_DATA"]    = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025C_w67.root";
mfile["GAM_2025C_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025C_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025D_DATA"]    = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025D_w67.root";
mfile["GAM_2025D_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025D_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025E_DATA"]    = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025E_w67.root";
mfile["GAM_2025E_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025E_MIX"]     = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025F_DATA"]    = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025F_w67.root";
mfile["GAM_2025F_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025F_MIX"]     = mfile["GAM_2025_MIX"];
mfile["GAM_2025G_DATA"]    = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025G_w67.root";
mfile["GAM_2025G_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025G_MIX"]     = mfile["GAM_2025_MIX"];
// w60: CDE, w62: CDEF, w64: CDEFG, w65: CDEFG (more G), w67: all (B)CDEFG
mfile["GAM_2025CDEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025CDEFG_w67.root";
//mfile["GAM_2025CDE_MC"]      = mfile["GAM_2025_MC"];
mfile["GAM_2025CDEFG_MC"]      = mfile["GAM_2025_MC"];
//mfile["GAM_2025CDEFG_MC"]    = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root"; // Summer24 MC
mfile["GAM_2025CDEFG_MIX"]   = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025DEFG_DATA"]  = "rootfiles/Prompt2025/Gam_w67/GamHistosFill_data_2025DEFG_w67.root";
mfile["GAM_2025DEFG_MC"]    = mfile["GAM_2025_MC"];
mfile["GAM_2025DEFG_MIX"]   = mfile["GAM_2025_MIX"];
//
mfile["GAM_2025C0_DATA"]    = "rootfiles/Prompt2025/Gam_w62/GamHistosFill_data_2025C_w62.root";
mfile["GAM_2025C0_MC"]   = mfile["GAM_2025_MC"];
mfile["GAM_2025C0_MIX"]   = mfile["GAM_2025_MIX"];
mfile["GAM_2025CT_DATA"]    = "rootfiles/Prompt2025/Gam_w63/GamHistosFill_data_2025C-TrkRadDamage_w63.root";
mfile["GAM_2025CT_MC"]   = mfile["GAM_2025_MC"];
mfile["GAM_2025CT_MIX"]   = mfile["GAM_2025_MIX"];


#endif
