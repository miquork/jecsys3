// Purpose: dump Rho, NPV etc. from root files into a text file
//          used by minitools/drawSigmaMBV2.C
#include "TFile.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>

void dumpRhoDijet();

void dumpRho() {
  dumpRhoDijet();
}

void dumpRhoDijet() {

  string vera[] = {
    "2022C","2022D","2022E","2022F","2022G",
    "2023Bv1","2023Cv1","2023Cv2","2023Cv3","2023Cv4","2023Dv1","2024Dv2"
  };
  const int nera = sizeof(vera)/sizeof(vera[0]);

  string vobs[] = {"RhoAll","Rho_C","Rho_CCPU","NPV","NPVGood"};//,"PUProfile"};
  const int nobs = sizeof(vobs)/sizeof(vobs[0]);
  
  string vmc[] = {"data","mc","mcpuoff"};
  const int nmc = sizeof(vmc)/sizeof(vmc[0]);

  ofstream fout("rootfiles/histogram_means_2022_2023.csv");
  
  for (int iera = 0; iera != nera; ++iera) {
    
    string sera = vera[iera];
    const char *cera = sera.c_str();
    //TString tera(cera);

    for (int iobs = 0; iobs != nobs; ++iobs) {

      const char *co = vobs[iobs].c_str();
    
      for (int imc = 0; imc != nmc; ++imc) {
	
	string sd = vmc[imc];
	const char *cd = sd.c_str();;
	
	TFile *f(0);
	if ((sera=="2022C"||sera=="2022D") && sd=="mc") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_mc_out_Summer22MG_v113_%s.root",cera),"READ");
	if ((sera=="2022C"||sera=="2022D") && sd=="mcpuoff") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_mc_out_Summer22MG_v113_2022.root"),"READ");
	if ((sera=="2022C"||sera=="2022D") && sd=="data") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_data_out_%s_nib1_JME_v113_2022.root",cera),"READ");
	//
	if ((sera=="2022E"||sera=="2022F"||sera=="2022G") && sd=="mc") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_mc_out_Summer22EEMG_v113_%s.root",cera),"READ");
	if ((sera=="2022E"||sera=="2022F"||sera=="2022G") && sd=="mcpuoff") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_mc_out_Summer22EEMG_v113_2022.root"),"READ");
	if ((sera=="2022E"||sera=="2022F"||sera=="2022G") && sd=="data") 
	  f = new TFile(Form("rootfiles/Prompt2024/v113_2022/jmenano_data_out_%s_nib1_JME_v113_2022.root",cera),"READ");
	//
	assert(f && !f->IsZombie());

	TH1D *h(0);
	if (sd=="data")
	  h = (TH1D*)f->Get(Form("HLT_PFJet500/Pileup/h_%s",co));
	else
	  h = (TH1D*)f->Get(Form("HLT_MC/Pileup/h_%s",co));
	assert(h);
	
	//Category,Source,Variable,Mean,MeanError
	//2024F,data,RhoAll,39.39716880,0.00273497
	//2024F,mc,RhoAll,42.35479369,0.00143170
	//2024F,mcpuoff,RhoAll,37.80844056,0.00163237
	
	cout << Form("%s,%s,%s,%1.4f,%1.4f",cera,cd,co,h->GetMean(),h->GetMeanError()) << endl;
	fout << Form("%s,%s,%s,%1.4f,%1.4f",cera,cd,co,h->GetMean(),h->GetMeanError()) << endl;
      } // for imc
    } // for iobs
  } // for iera
} // dumpRhoDijet
