{
  // For JEC
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L test/testJERSF.C+g");
  
  testJERSF("pdf/JERSF/Summer23_2023Cv123_JRV1_MC_SF_AK4PFPuppi.txt");
  testJERSF("pdf/JERSF/Summer23_2023Cv4_JRV1_MC_SF_AK4PFPuppi.txt");
  testJERSF("pdf/JERSF/Summer23_2023D_JRV1_MC_SF_AK4PFPuppi.txt");
}
