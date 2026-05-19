# jecsys3
CMS Run 3 jet energy scale

Start rebuilding jecsys for Run3 based on Run2 jecsys UL2017 branch.
Simplify reprocess.C to only work on alpha<1.0.
FSR corrections with HDM method only.
Base global fit on globalFitRun2.C.

Quick starter guide:

`emacs -nw Config.C`
[edit file to add new inputs files; NB: use ++g to recompile codes after!]

`root -l -b -q JetVeto.C++g`
[veto maps stored rootfiles/jetveto20*.root]

`root -l -b -q 'L2Res.C++g(0,"2025CDEFG","",-1)'` [single run test]
`python3 minitools/runAllL2Res.py` [batch process]
[L2Res corrections stored in textfiles/Prompt/*L2Residual[VsPtRef[Asymm]]*.txt]

`root -l -b -q 'mk_reprocess_RunEpoch.C("2025CDEFG",0,1)` [single run test]
`python3 minitools/runAllIOVs.py` [batch process]
[L3Res results stored in rootfiles/jecdata20*.root]

`root -l -b -q createL2L3ResTextFileV2.C++g` [batch process]
[L2L3Res corrections stored in textfiles/Prompt/*L2L3Residual[VsPtRefAsymmm]*.txt]

`root -l -b -q JERSF.C++g(0,"2025CDEFG","Summer24MG_NOJERSF")` [single run test]
`python3 minitools/runAllJERSF.py` [batch process]
[JER SF corrections stored in textfiles/Prompt/*SF*.txt]

Core plots:

`root -l -b -q minitools/drawL2ResVsTime.C`
`root -l -b -q minitools/drawL3ResVsTime.C`
`root -l -b -q minitools/drawJERSFvsTime.C`

Text file validation and plots:

`root -l -b -q test/mk_drawCMSresponse.C` [L2L3Res, L2Res]
`root -l -b -q test/mk_testJERSF.C` [JER SF]