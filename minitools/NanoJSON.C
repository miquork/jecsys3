// Purpose: run over a set of NanoAOD files and produce output JSON for brilcalc
#include "TChain.h"
#include "TTree.h"

#include <iostream>
#include <map>
#include <set>

void NanoJSON() {

  // Header file for the classes stored in the TTree if any.
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  // Declaration of leaf types
  UInt_t          run;
  UInt_t          luminosityBlock;

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_luminosityBlock;   //!

  // Open NanoAODs in a new chain
  TChain *c = new TChain("Events");
  c->AddFile("/Users/voutila/Downloads/Data_2024I_DiJet_ZeroBiasv1_2024I-nib1-fib3-386592-386594.root");
  
  // Set branch addresses and branch pointers
  fChain = c;
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);

  // List of all Run,LS pairs encountered
  map<int, set<int> > runls;

  // Loop over all events to find active lumisections
  Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry = 0; jentry != nentries; jentry++) {
    //Long64_t ientry = LoadTree(jentry);
    Long64_t ientry = fChain->LoadTree(jentry);
    fChain->GetEntry(jentry);

    if (jentry%10000==0) cout << "." << flush;
    runls[run].insert(luminosityBlock);
    
  }
  cout << endl;

  // List all run,LS pairs in JSON format
  typedef map<int, set<int> >::const_iterator IT;
  typedef set<int>::const_iterator JT;
  for (IT it = runls.begin(); it != runls.end(); ++it) {
    int run = it->first;
    set<int> const& sls = it->second;
    int ls_begin = (*(sls.begin()));
    
    for (JT jt = sls.begin(); jt != sls.end(); ++jt) {

      // Get the current LS
      int ls = (*jt);

      // Peek ahead at the next element
      JT nxt = next(jt);
      
      // End of run: printout
      if (nxt == sls.end()) {
	cout << Form("%d: [%d, %d]\n",run,ls_begin,ls);
      }
      // Break in contiguous LS numbers: printout and start new block
      else if ((*nxt)-ls>2) {
	cout << Form("%d: [%d, %d]\n",run,ls_begin,ls);
	ls_begin = (*nxt);
      }
    } // for jt
  } // for it
  
}
