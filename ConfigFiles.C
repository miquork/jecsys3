// Purpose: common file listings for L2Res.C, L3Res.C/globalFit, JERSF.C
#ifndef __CONFIGFILES_C__
#define __CONFIGFILES_C__

#include <map>
#include <string>
std::map<std::string, map<std::string, map<std::string, std::string> > > mf;

// 2024 partial re-reco + Summer24
string t = "24V9M";

// Options for Zmm: v100. v100pu. v100noPU_smear. v100pu_smear.root
// Use find-and-replace to switch
// v100*smear.root are corrupt, bad JRV9M file
mf[t]["Zmm"]["2024C_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDEReprocessing_v1_2024C_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024D_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDEReprocessing_v1_2024D_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024E_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDEReprocessing_v1_2024E_nib1_Zmm_V9M_v100pu.root";
//mf[t]["Zmm"]["2024Ev1_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDEReprocessing_v1_2024Ev1_nib1_Zmm_V9M_v100pu.root";
//mf[t]["Zmm"]["2024Ev2_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDEReprocessing_v1_2024Ev2_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024F_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024F_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024F_nib2"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024F_nib2_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024F_nib3"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024F_nib3_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024G_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024G_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024G_nib2"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024G_nib2_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024H_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024H_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024I_nib1"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024I_nib1_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024CDE_nib"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024CDE_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024FGHI_nib"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024FGHI_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["2024_nib"] = "rootfiles/Prompt2024/Zmm_v100/jme_Zj_2024_Zmm_V9M_v100pu.root";
mf[t]["Zmm"]["MC"] = ""; // 'mc' directory in input files

mf[t]["Gam"]["2024C_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024C-rereco_w56.root";
mf[t]["Gam"]["2024D_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024D-rereco_w56.root";
mf[t]["Gam"]["2024E_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024E-rereco_w56.root";
mf[t]["Gam"]["2024F_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Fnib1_w56.root";
mf[t]["Gam"]["2024F_nib2"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Fnib2_w56.root";
mf[t]["Gam"]["2024F_nib3"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Fnib3_w56.root";
mf[t]["Gam"]["2024G_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Gnib1_w56.root";
mf[t]["Gam"]["2024G_nib2"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Gnib2_w56.root";
mf[t]["Gam"]["2024H_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Hnib1_w56.root";
mf[t]["Gam"]["2024I_nib1"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024Inib1_w56.root";
mf[t]["Gam"]["2024CDE_nib"]  = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024CDE_V9M_w56.root";
mf[t]["Gam"]["2024FGHI_nib"] = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024FGHI_V9M_w56.root";
mf[t]["Gam"]["2024_nib"]   = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_data_2024_V9M_w56.root";
mf[t]["Gam"]["MC"]         = "rootfiles/Prompt2024/w48_Gam/GamHistosFill_mc_summer2024P8_pu-2024CDEFGHI_w48.root";
//mf[t]["Gam"]["MC"]         = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_mc_summer2024P8_no-pu_w56.root"; // no QCD, no PU RW
//mf[t]["Gam"]["MC"]         = "rootfiles/Prompt2024/Gam_w56/GamHistosFill_mc_summer2024QCD_no-pu_w56.root"; // QCD, no PU RW; unstable stats, revert

mf[t]["Jet"]["2024C_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024C_Rp_JME_v134.root";
mf[t]["Jet"]["2024D_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024D_Rp_JME_v134.root";
mf[t]["Jet"]["2024E_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024E_Rp_JME_v134.root";
mf[t]["Jet"]["2024F_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024F_nib1_JME_v134.root";
mf[t]["Jet"]["2024F_nib2"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024F_nib2_JME_v134.root";
mf[t]["Jet"]["2024F_nib3"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024F_nib3_JME_v134.root";
mf[t]["Jet"]["2024G_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024G_nib1_JME_v134.root";
mf[t]["Jet"]["2024G_nib2"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024G_nib2_JME_v134.root";
mf[t]["Jet"]["2024H_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024H_nib1_JME_v134.root";
mf[t]["Jet"]["2024I_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024I_nib1_JME_v134.root";
//mf[t]["Jet"]["2024Iv1_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024Iv1_nib1_JME_v134.root";
//mf[t]["Jet"]["2024Iv2_nib1"] = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024Iv2_nib1_JME_v134.root";
mf[t]["Jet"]["2024CDE_nib"]   = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024CDE_JME_v134.root";
mf[t]["Jet"]["2024FGHI_nib"]   = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024FGHI_JME_v134.root";
mf[t]["Jet"]["2024_nib"]   = "rootfiles/Prompt2024/Jet_v134/jmenano_data_cmb_2024CDEFGHI_JME_v134.root";
//mf[t]["Jet"]["MC"] = "rootfiles/Prompt2024/Jet_v132/jmenano_mc_out_Summer24MG_v132.root"; // no JER SF? no PU RW?
//mf[t]["Jet"]["MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_v134.root"; // no JER SF, no PU RW
mf[t]["Jet"]["MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_PU69_v134.root"; // no JER SF, with PU RW
//mf[t]["Jet"]["MC"] = "rootfiles/Prompt2024/Jet_v134/jmenano_mc_out_Summer24MG_JERSF_v134.root"; // with JER SF, no PU RW




#endif
