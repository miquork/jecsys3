{
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L test/drawCMSresponse.C+g");

  gROOT->ProcessLine(".! touch pdf/test");

  for (int i = 0; i != 2; ++i) {
  Bool_t l2res = (i!=0);
  //drawCMSresponse("Prompt24_Run2024B_nib1_V8M","Prompt24B_nib1_V8M",l2res);
  /*
  drawCMSresponse("Prompt24_Run2024C_nib1_V9M","Prompt24C_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024D_nib1_V9M","Prompt24D_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024Ev1_nib1_V9M","Prompt24Ev1_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024Ev2_nib1_V9M","Prompt24Ev2_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024F_nib1_V9M","Prompt24F_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024F_nib2_V9M","Prompt24F_nib2_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024F_nib3_V9M","Prompt24F_nib3_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024G_nib1_V9M","Prompt24G_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024G_nib2_V9M","Prompt24G_nib2_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024H_nib1_V9M","Prompt24H_nib1_V9M",l2res);
  drawCMSresponse("Prompt24_Run2024I_nib1_V9M","Prompt24I_nib1_V9M",l2res);
  */
  //drawCMSresponse("Prompt25_Run2025C_V1M","Prompt25C_V1M",l2res);
  /*
  drawCMSresponse("ReReco24_Run2024C_nib1_V9M","ReReco24C_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024D_nib1_V9M","ReReco24D_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024E_nib1_V9M","ReReco24E_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024F_nib1_V9M","ReReco24F_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024F_nib2_V9M","ReReco24F_nib2_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024F_nib3_V9M","ReReco24F_nib3_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024G_nib1_V9M","ReReco24G_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024G_nib2_V9M","ReReco24G_nib2_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024H_nib1_V9M","ReReco24H_nib1_V9M",l2res);
  drawCMSresponse("ReReco24_Run2024I_nib1_V9M","ReReco24I_nib1_V9M",l2res);
  */

  /*
  drawCMSresponse("Prompt25_Run2025CDE_V2M","Prompt25CDE_V2M",l2res);
  drawCMSresponse("Prompt25_Run2025C_V2M","Prompt25C_V2M",l2res);
  drawCMSresponse("Prompt25_Run2025D_V2M","Prompt25D_V2M",l2res);
  drawCMSresponse("Prompt25_Run2025E_V2M","Prompt25E_V2M",l2res);
  */

  drawCMSresponse("Prompt25_Run2025CDEFG_V3M","Prompt25CDEFG_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025DEFG_V3M","Prompt25DEFG_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025C_V3M","Prompt25C_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025D_V3M","Prompt25D_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025E_V3M","Prompt25E_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025F_V3M","Prompt25F_V3M",l2res);
  drawCMSresponse("Prompt25_Run2025G_V3M","Prompt25G_V3M",l2res);
  }
  
  //drawCMSresponse("Winter24Run3_V1_MC","Winter24Run3_V1_MC", false);
  //drawCMSresponse("Winter24Run3_V1_MC","Winter24Run3_V1_MC", true);
  
  
}
