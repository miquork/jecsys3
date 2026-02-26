// Purpose: test JER SF text file produced by JERSF.C
{
  // For JEC
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L test/testMCJER.C+g");

  testMCJER("textfiles/Prompt25_JRV3M/Summer24MG_JRV3M_MC_PtResolution_AK4PFPuppi.txt", "Summer24 AK4PFPuppi");
	    //"Summer24 QCD MadGraph");
}
