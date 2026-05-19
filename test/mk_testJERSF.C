// Purpose: test JER SF text file produced by JERSF.C
{
  // For JEC
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L test/testJERSF.C+g");

  /*
  testJERSF("pdf/JERSF/Summer23_2023Cv123_JRV1_MC_SF_AK4PFPuppi.txt");
  testJERSF("pdf/JERSF/Summer23_2023Cv4_JRV1_MC_SF_AK4PFPuppi.txt");
  testJERSF("pdf/JERSF/Summer23_2023D_JRV1_MC_SF_AK4PFPuppi.txt");
  */
  
  //testJERSF("textfiles/ReReco24_V9M/ReReco24_2024CDEFGHI_nib_JRV9M_MC_SF_AK4PFPuppi.txt"); // Buggy file!
  //testJERSF("textfiles/ReReco24/ReReco24_2024_nib_JRV9M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/ReReco24_V9M_goodJER/ReReco24_2024_nib_JRV9M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/ReReco24_V10M_JERSF/ReReco24_2024_nib_JRV10M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/Prompt25/Prompt25_2025C_JRV1M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/Prompt25/Prompt25_2025CDEFG_JRV2M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/Prompt25_JRV3M/Prompt25_2025CDEFG_JRV3M_MC_SF_AK4PFPuppi.txt");
  //testJERSF("textfiles/Prompt26_V0M/Prompt26_2026B_JRV0M_MC_SF_AK4PFPuppi.txt");

  testJERSF("textfiles/Prompt/Prompt24_2024_nib_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024CDE_nib_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024FGHI_nib_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024C_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024D_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024E_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024F_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024F_nib2_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024F_nib3_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024G_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024G_nib2_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024H_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt24_2024I_nib1_JRV10M_MC_SF_AK4PFPuppi.txt");
  
  testJERSF("textfiles/Prompt/Prompt25_2025CDEFG_JRV4M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt25_2025C_JRV4M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt25_2025D_JRV4M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt25_2025E_JRV4M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt25_2025F_JRV4M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt25_2025G_JRV4M_MC_SF_AK4PFPuppi.txt");
    
  testJERSF("textfiles/Prompt/Prompt26_2026B_JRV1M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt26_2026C_JRV1M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt/Prompt26_2026D_JRV1M_MC_SF_AK4PFPuppi.txt");
}
