{
  // make sure asserts are run
  #undef NDEBUG
  // does not seem to be enough, also need to combile with +g

  // For JEC residual (and pile-up)
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/Utilities.cc+");
  gROOT->ProcessLine(".L CondFormats//JetMETObjects/src/JetCorrectorParameters.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrector.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc+");
  //
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+");
  gROOT->ProcessLine(".L CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc+");

  // Compile with +g to make sure asserts are run
  gROOT->ProcessLine(".L tools.C+g");
  gROOT->ProcessLine(".L Flavor.C+g");
  gROOT->ProcessLine(".L reprocess.C+g");
  gROOT->ProcessLine(".L scaleJES.C+g");
  gROOT->ProcessLine(".L softrad3.C+g");
  //gROOT->ProcessLine(".L globalFitSyst.C+g");
  //gROOT->ProcessLine(".L globalFitRenormPF.C+g");
  gROOT->ProcessLine(".L globalFit.C+g");

  // Merge inputs from separate groups
  // NB: this does not need to be run, if the merged inputs
  //     are already available in 'rootfiles/jecdata.root'
  string epoch = "Run22C";//"Run23B";//"Run3";//"RunCD";
  #ifdef epochname
  std::cout << epoch.c_str()<< std::endl;
  std::cout << inputepoch.c_str()<< std::endl;
  epoch = inputepoch;
  #endif

  // Read in files from different groups and merge them in jecdata[epoch].root
  reprocess(epoch); // Comment out if using archived jecdata[epoch].root

  // Code for 19Dec2023 use only!
  // Scale Z+jet, gamma+jet (and multijet) from 22Sep2023 to 19Dec2023
  //if (epoch!="Run3" && epoch!="Run22FG") { // Run3,22FG inputs already scaled
  //scaleJES(epoch, "zjet");
  //}
  //if (epoch!="Run22FG") { // Run22FG inputs already scaled (19Dec)
  //scaleJES(epoch, "gamjet");
  //}
  //scaleJES(epoch, "multijet"); // 19Dec no scaling needed
  if (epoch=="Run24CS") {
    scaleJES(epoch, "zjet");
    scaleJES(epoch, "zjav");
    scaleJES(epoch, "jetz");
  }
  
  // HDM method: use HT decomposition (lead, soft jets, unclustered) for FSR
  softrad3(0.0,1.3,epoch); // 3-point FSR

  // Produce central systematic uncertainties for globalFitL3Res
  //globalFitSyst(epoch);     // also for globalFitRun2.C
  //globalFitRenormPF(epoch); // for globalFitRun2.C

  // Run global fit
  /////////////////
  //globalFitEtaBin(0.0, 1.3, epoch, "22Sep2023reV3");
  //globalFitEtaBin(0.0, 1.3, epoch, "19Dec2023");
  //globalFitEtaBin(0.0, 1.3, epoch, "Summer23");
  //globalFitEtaBin(0.0, 1.3, epoch, "Summer23_V1test");
  //globalFitEtaBin(0.0, 1.3, epoch, "Summer23_V2test2");
  //globalFitEtaBin(0.0, 1.3, epoch, "Summer23_V2closure");
  //globalFitEtaBin(0.0, 1.3, epoch, "Summer23_V2closure2");
  //globalFitEtaBin(0.0, 1.3, epoch, "Prompt24_3fb");
  //globalFitEtaBin(0.0, 1.3, epoch, "Prompt24_3fb_V2Mclosure");
  //globalFitEtaBin(0.0, 1.3, epoch, "Prompt24_12fb_V3M");
  globalFitEtaBin(0.0, 1.3, epoch, "Prompt24_V3Mclosure");
  
  exit(0); // Avoid page full of THastList::Delete errors
}
