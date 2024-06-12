{
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");

  gROOT->ProcessLine(".L test/drawCMSresponse.C+");

  //drawCMSresponse();
  /*
  drawCMSresponse("Summer22-22Sep2023_Run2022CD_V3","22CD-22Sep2023_V3",false);
  drawCMSresponse("Summer22-22Sep2023_Run2022CD_V3","22CD-22Sep2023_V3",true);
  drawCMSresponse("Summer22EE-22Sep2023_Run2022E_V3","22E-22Sep2023_V3",false);
  drawCMSresponse("Summer22EE-22Sep2023_Run2022E_V3","22E-22Sep2023_V3",true);
  drawCMSresponse("Summer22EEPrompt22_Run2022F_V3","22F-Prompt_V3",false);
  drawCMSresponse("Summer22EEPrompt22_Run2022F_V3","22F-Prompt_V3",true);
  drawCMSresponse("Summer22EEPrompt22_Run2022G_V3","22G-Prompt_V3",false);
  drawCMSresponse("Summer22EEPrompt22_Run2022G_V3","22G-Prompt_V3",true);
  drawCMSresponse("Summer22Prompt23_Run2023Cv123_V3","23C123-Prompt_V3",false);
  drawCMSresponse("Summer22Prompt23_Run2023Cv123_V3","23C123-Prompt_V3",true);
  drawCMSresponse("Summer22Prompt23_Run2023Cv4_V3","23C4-Prompt_V3",false);
  drawCMSresponse("Summer22Prompt23_Run2023Cv4_V3","23C4-Prompt_V3",true);
  drawCMSresponse("Summer22Prompt23_Run2023D_V3","23D-Prompt_V3",false);
  drawCMSresponse("Summer22Prompt23_Run2023D_V3","23D-Prompt_V3",true);
  */
  /*
  drawCMSresponse("Run22CD-22Sep2023_DATA","22CD-22Sep2023");
  drawCMSresponse("Run22E-22Sep2023_DATA","22E-22Sep2023");
  drawCMSresponse("Run22F-Prompt_DATA","22F-Prompt");
  drawCMSresponse("Run22G-Prompt_DATA","22G-Prompt");
  drawCMSresponse("Run23C123-Prompt_DATA","23C123-Prompt");
  drawCMSresponse("Run23C4-Prompt_DATA","23C4-Prompt");
  drawCMSresponse("Run23D-Prompt_DATA","23D-Prompt");
  drawCMSresponse("Run23C4D-Prompt_DATA","23C4D-Prompt");
  */

  /*
  drawCMSresponse("Summer23Prompt23_Run2023Cv123_V1","23C123-Summer23_V2",true);
  drawCMSresponse("Summer23Prompt23_Run2023Cv4_V1","23C4-Summer23_V2",true);
  drawCMSresponse("Summer23Prompt23_Run2023D_V1","23D-Summer23_V2",true);

  drawCMSresponse("Summer23Prompt23_Run2023Cv123_V2","23C123-Summer23_V2",false);
  drawCMSresponse("Summer23Prompt23_Run2023Cv4_V2","23C4-Summer23_V2",false);
  drawCMSresponse("Summer23Prompt23_Run2023D_V2","23D-Summer23_V2",false);
  */

  //drawCMSresponse("Prompt24_Run2024BC_V1M","Prompt24BC_V1M",true);
  //drawCMSresponse("Prompt24_Run2024BC_V1M","Prompt24BC_V1M",false);
  //drawCMSresponse("Prompt24_Run2024BC_V2M","Prompt24BC_V2M",true);
  //drawCMSresponse("Prompt24_Run2024BC_V2M","Prompt24BC_V2M",false);

  //drawCMSresponse("Summer23Prompt23_Run2023D_V2xHFscale","23D_V2xHF",false);

  /*
  drawCMSresponse("Prompt24_Run2024BCD_V3M","Prompt24BCD_V3M",true);
  drawCMSresponse("Prompt24_Run2024BCD_V3M","Prompt24BCD_V3M",false);
  drawCMSresponse("Prompt24_Run2024CR_V3M","Prompt24CR_V3M",true);
  drawCMSresponse("Prompt24_Run2024CR_V3M","Prompt24CR_V3M",false);
  */
  drawCMSresponse("Prompt24_Run2024BCD_V4M","Prompt24BCD_V4M",true);
  drawCMSresponse("Prompt24_Run2024BCD_V4M","Prompt24BCD_V4M",false);
  drawCMSresponse("Prompt24_Run2024E_V4M","Prompt24E_V4M",true);
  drawCMSresponse("Prompt24_Run2024E_V4M","Prompt24E_V4M",false);
  drawCMSresponse("Prompt24_Run2024CR_V4M","Prompt24CR_V4M",true);
  drawCMSresponse("Prompt24_Run2024CR_V4M","Prompt24CR_V4M",false);
  drawCMSresponse("Prompt24_Run2024CS_V4M","Prompt24CS_V4M",true);
  drawCMSresponse("Prompt24_Run2024CS_V4M","Prompt24CS_V4M",false);
}
