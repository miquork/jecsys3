void mk_reprocess_RunEpoch(std::string inepoch="", int doClosure=-1,
			   int doRecompile=-1){
  #define epochname
  gROOT->ProcessLine(Form("string inputepoch(\"%s\");",inepoch.c_str()));
  #define closurebool
  gROOT->ProcessLine(Form("int inputclosure(%d);",doClosure));
  #define recompilebool
  gROOT->ProcessLine(Form("int recompile(%d);",doRecompile));
  gROOT->ProcessLine(".x mk_reprocess.C");
}
