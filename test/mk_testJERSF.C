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
  testJERSF("textfiles/ReReco24/ReReco24_2024_nib_JRV9M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/ReReco24_V9M_goodJER/ReReco24_2024_nib_JRV9M_MC_SF_AK4PFPuppi.txt");
  testJERSF("textfiles/Prompt25/Prompt25_2025C_JRV1M_MC_SF_AK4PFPuppi.txt");
}
