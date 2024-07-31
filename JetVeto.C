// Purpose: Produce jet veto maps from JMENANO analyzer on dijets
//          Updated version for 2024. Sum absolute pulls to avoid cancellation
// Author:  Mikko Voutilainen (at) cern (dot) ch
// Date:    2024-06-14
// NB: Update of minitools/doJetVetoV2.C to match L2Res.C, JERSF.C, L3Res.C
#include "TFile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TMath.h"

#include <string>
#include <vector>

#include "tdrstyle_mod22.C"

// Threshold for veto maps
double pullThreshold = 70; // default: 70

// Separate threshold for outer HF
double pullThresholdHF45 = 70; // default: 70

// Minimum of non-empty towers to consider eta strip
int nMinTowers = 70; // default: 70 (out of 72)

// Use pulls instead of relative changes
bool doPull = true; // default: true

// plot results from each step
bool plotJetVeto = true; // default: true
bool plotJetVetoOnJES = true;
bool plotJetVetoOnJESnorm = false;

// If this is non-zero (not ""), will use just this one trigger
string oneTrig = ""; // default: ""
// Possible options to comment out
//string oneTrig = "HLT_PFJet500";
//string oneTrig = "HLT_PFJet450";
//string oneTrig = "HLT_PFJet320";
//string oneTrig = "HLT_PFJet40";
//string oneTrig = "HLT_PFJet60";
//string oneTrig = "HLT_PFJet80";
//string oneTrig = "HLT_PFJet140";
//string oneTrig = "HLT_PFJet260";

// If this is non-zero (not ""), will use just this one histogram
string oneHist = ""; // default: ""
// Possible options to comment out
//string oneHist = "p2asymm";
//string oneHist = "h2pt";
//string oneHist = "p2chf";
//string oneHist = "p2nhf";


void JetVetos(string run, string version);
void JetVeto(string run = "", string version = "vx") {

  if (run!="") { JetVeto(run,version); exit(0); }
  /*
  JetVetos("2022CD",version);
  JetVetos("2022EFG",version);
  JetVetos("2023BC",version);
  JetVetos("2023D",version);
  */
  
  JetVetos("2024BCD",version);
  JetVetos("2024E",version);
}

void JetVetos(string run, string version) {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *fout = new TFile(Form("rootfiles/jetveto%s.root",run.c_str()),
			  "RECREATE");
  fout->mkdir("trigs");
  
  TFile *f(0);
  if (run=="2022CD") { // Summer22
    lumi_136TeV = "Run2022CD, 8.1 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230816/jmenano_data_out_2022CD_v1.root","READ");
  }
  if (run=="2022EFG") { // Summer22BPix
    lumi_136TeV = "Run2022EFG, 27.0 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230816/jmenano_data_out_2022EF_v1.root","READ");
  }  
  if (run=="2023BC") { // Summer23
    lumi_136TeV = "Run2023BC, 13.6 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230814/nano_data_out_2023BC_v1.root","READ");
  }
  if (run=="2023D") { // Summer23BPix
    lumi_136TeV = "Run2023D, 9.5 fb^{-1}";
    f = new TFile("rootfiles/Iita_20230814/nano_data_out_2023D_v1.root","READ");
  }
  if (run=="2024BCD") {
    lumi_136TeV = "Run2024BCD, 12.3 fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v50_2024/jmenano_data_cmb_2024BCD_JME_v50_2024.root","READ"); // May 16 golden, 12.3/fb => some problem with cmb?
    f = new TFile("rootfiles/Prompt2024/v50_2024/jmenano_data_out_2024BCD_JME_v50_2024.root","READ"); // May 16 golden, 12.3/fb
    //pullThreshold = 250;//200;//150;//100;//85;//100;//70;//50;
    //pullThresholdHF45 = 300;
    nMinTowers = 50; // for BPix hole
  }
  if (run=="2024E") {
    lumi_136TeV = "Run2024E, X.X fb^{-1}";
    //f = new TFile("rootfiles/Prompt2024/v76_2024/jmenano_data_cmb_2024E_JME_v76_2024.root","READ"); // June 6 hybrid, 27.0/fb => some problem with cmb?
    f = new TFile("rootfiles/Prompt2024/v76_2024/jmenano_data_out_2024E_JME_v76_2024.root","READ"); // June 6 hybrid, 27.0/fb
    //pullThreshold = 250;//200;//150;//100;//85;//100;//70;//50;
    //pullThresholdHF45 = 300;
    nMinTowers = 50; // for BPix hole
  }

  // Don't use *cmb* files, trigger folder uncertainties messed up!!
  assert(!TString(f->GetName()).Contains("_cmb_"));

  assert(f && !f->IsZombie());

  vector<string> vtrg;
  //vtrg.push_back("HLT_ZeroBias"); // analyze correct PD

  vtrg.push_back("HLT_PFJet40");
  vtrg.push_back("HLT_PFJet60");
  vtrg.push_back("HLT_PFJet80");
  vtrg.push_back("HLT_PFJet140");
  vtrg.push_back("HLT_PFJet200");
  vtrg.push_back("HLT_PFJet260");
  vtrg.push_back("HLT_PFJet320");
  vtrg.push_back("HLT_PFJet450");
  vtrg.push_back("HLT_PFJet500");
  //vtrg.push_back("HLT_PFJet550");

  //if (run!="2024BC") {
  vtrg.push_back("HLT_PFJetFwd40");
  vtrg.push_back("HLT_PFJetFwd60");
  vtrg.push_back("HLT_PFJetFwd80");
  // NB: these had wrong |eta| threshold in v12+v13. To be fixed in v14
  vtrg.push_back("HLT_PFJetFwd140");
  vtrg.push_back("HLT_PFJetFwd200");
  vtrg.push_back("HLT_PFJetFwd260");
  vtrg.push_back("HLT_PFJetFwd320");
  vtrg.push_back("HLT_PFJetFwd400");
  vtrg.push_back("HLT_PFJetFwd450");
  vtrg.push_back("HLT_PFJetFwd500");
  //} 
  //else { // 2024BCD
  vtrg.push_back("HLT_ZeroBias");
  //}
  
  // Add dijet average triggers for better h2jes
  vtrg.push_back("HLT_DiPFJetAve40");
  vtrg.push_back("HLT_DiPFJetAve60");
  vtrg.push_back("HLT_DiPFJetAve80");
  vtrg.push_back("HLT_DiPFJetAve140");
  vtrg.push_back("HLT_DiPFJetAve200");
  vtrg.push_back("HLT_DiPFJetAve260");
  vtrg.push_back("HLT_DiPFJetAve320");
  vtrg.push_back("HLT_DiPFJetAve400");
  vtrg.push_back("HLT_DiPFJetAve500");
  
  vtrg.push_back("HLT_DiPFJetAve60_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve80_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve100_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve160_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve220_HFJEC");
  vtrg.push_back("HLT_DiPFJetAve300_HFJEC");
  
  if (oneTrig!="") {
    vtrg.clear();
    vtrg.push_back(oneTrig);
  }
  int ntrg = vtrg.size();

  vector<string> vh;
  vh.push_back("p2asymm");
  vh.push_back("h2phieta");//h2pt");
  vh.push_back("p2nef");
  vh.push_back("p2chf");
  vh.push_back("p2nhf");
  if (oneHist!="") {
    vh.clear();
    vh.push_back(oneHist);
  }
  int nh = vh.size();

  TH2D *h2nomsums(0), *h2nomrefsum(0), *h2jes(0);
  TH2D *h2abssums(0), *h2absrefsum(0);
  for (int ih = 0; ih != nh; ++ih) {

    string hname = vh[ih];
    TH2D *h2nomsum(0), *h2abssum(0);
  
    for (int itrg = 0; itrg != ntrg; ++itrg) {

      string trg = vtrg[itrg];
      //string trg = "HLT_PFJet500";
      //string trg = "HLT_ZeroBias";
      //string hname = "p2nef";
      //isProfile2D = (hname[0]='p' && hname[1]='2');

      string objname = Form("%s/Jetveto/%s",trg.c_str(),hname.c_str());
      TObject *obj = f->Get(objname.c_str());
      if (!obj) {
	cout << "Missing " << objname << endl << flush;
      }
      assert(obj);
      assert(obj->InheritsFrom("TH2D"));
      bool isProf2D = obj->InheritsFrom("TProfile2D");
      
      TH2D *h2 = (isProf2D ? ((TProfile2D*)obj)->ProjectionXY() : (TH2D*)obj);
      h2->UseCurrentStyle();
      TH2D *h2nom = (TH2D*)h2->Clone(Form("h2nom_%s",hname.c_str()));
      TH2D *h2abs = (TH2D*)h2->Clone(Form("h2abs_%s",hname.c_str()));
      
      // Calculate average JES shift
      if (hname=="p2asymm") {
	if (!h2jes) {
	  h2jes = (TH2D*)h2->Clone("h2jes");
	}
	else {
	  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
	    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
	      if (h2->GetBinError(i, j)!=0) {
		if (h2jes->GetBinError(i, j)!=0) {
		  double val1 = h2jes->GetBinContent(i,j);
		  double err1 = h2jes->GetBinError(i,j);
		  double n1 = 1./pow(err1,2);
		  double val2 = h2->GetBinContent(i,j);
		  double err2 = h2->GetBinError(i,j);
		  double n2 = 1./pow(err2,2);
		  double val = (n1*val1 + n2*val2) / (n1+n2);
		  double err = (n1*err1 + n2*err2) / (n1+n2);
		  h2jes->SetBinContent(i,j,val);
		  h2jes->SetBinError(i,j,err);
		}
		else if (h2->GetBinContent(i,j)!=0) {
		  h2jes->SetBinContent(i,j,h2->GetBinContent(i,j));
		  h2jes->SetBinError(i,j,h2->GetBinError(i,j));
		}
	      }
	    } // for j
	  } // for i
	}
      } // p2asymm
      
      // Normalize eta strips vs phi
      for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
	
	// Recreate histogram with automatic binning for each eta bin
	TH1D *htmp = new TH1D("htmp","Distribution of values",100,-1,-1);
	htmp->SetBuffer(72.);
    
	int ny = h2->GetNbinsY();
	int n(0);
	for (int j = 1; j != ny+1; ++j) {
	  if (h2->GetBinError(i,j)!=0) {
	    htmp->Fill(h2->GetBinContent(i, j));	
	    ++n;
	  }
	} // for j
	//htmp->Draw();
	
	// If not enough entries, clear column and continue
	//if (htmp->GetEntries()<70. || htmp->Integral()!=0) {
	if (n<nMinTowers) {
	  for (int j = 1; j != ny+1; ++j) {
	    h2->SetBinContent(i, j, 0.);
	    h2->SetBinError(i, j, 0.);
	  } // for j
	}
	// 
	else {
	  // Determine median and 68% quantile ([0.16,0.84] range divided by 2)
	  const Int_t nprobSum = 3;
	  Double_t probSum[nprobSum] = {0.16, 0.50, 0.84};
	  Double_t q[nprobSum];
	  htmp->GetQuantiles(nprobSum, q, probSum);
	  double median = q[1];
	  double q68 = 0.5*(q[2]-q[0]);
	  
	  // Determine pull width
	  for (int j = 1; j != ny+1; ++j) {
	    if (h2->GetBinError(i,j)!=0 && q68!=0) {
	      if (doPull) {
		h2abs->SetBinContent(i, j, fabs(h2->GetBinContent(i, j)-median) / q68 - 1); // -1 to model the normal cancellation of +1 and -1 fluctuations
		h2abs->SetBinError(i, j, h2->GetBinError(i, j) / q68);
		h2nom->SetBinContent(i, j, (h2->GetBinContent(i, j)-median) / q68);
		h2nom->SetBinError(i, j, h2->GetBinError(i, j) / q68);
		
	      }
	      else {
		// TBD: code not checked for h2abs
		if (median<0.5) {
		  h2abs->SetBinContent(i,j,fabs(h2->GetBinContent(i,j)-median)/
				       (1+fabs(median)));
		  h2abs->SetBinError(i, j, h2->GetBinError(i, j) /
				     (1+fabs(median)));
		  h2nom->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/
				       (1+fabs(median)));
		  h2nom->SetBinError(i, j, h2->GetBinError(i, j) /
				     (1+fabs(median)));

		}
		else {
		  h2abs->SetBinContent(i,j,fabs(h2->GetBinContent(i,j)-median)/
				       median);
		  h2abs->SetBinError(i, j, h2->GetBinError(i, j) / median);
		  h2nom->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/median);
		  h2nom->SetBinError(i, j, h2->GetBinError(i, j) / median);
		}
	      }
	    }
	  } // for j
	}
	
	delete htmp;
      } // for i
      
      //h2->Draw("COLZ");
      if (!h2nomsum && !h2abssum) {
	h2nomsum = (TH2D*)h2nom->Clone(Form("h2nomsum_%s",hname.c_str()));
	h2abssum = (TH2D*)h2abs->Clone(Form("h2abssum_%s",hname.c_str()));
      }
      else {
	h2nomsum->Add(h2nom);
	h2abssum->Add(h2abs);
      }
      fout->cd("trigs");
      h2nom->Write(Form("jetpullmap_nom_%s_%s",hname.c_str(),trg.c_str()));
      curdir->cd();
      delete h2;
					
    } // for itrg

    double ntrg2 = (oneTrig!="" ? 1 : max(1, ntrg/2));

    const char *c1nomname = Form("c1nom_%s",hname.c_str());
    TH1D *h1nom = tdrHist(Form("h1nom_%s_%s",oneTrig.c_str(),hname.c_str()),
			  "#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
    if (doPull) h1nom->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h1nom->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c1nom = tdrCanvas(c1nomname,h1nom,8,11,kRectangular);
    c1nom->SetRightMargin(0.15);
    h2nomsum->Draw("COLZ SAME");
    if (doPull) h2nomsum->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h2nomsum->GetZaxis()->SetRangeUser(-0.5,0.5);

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.035);
    tex->DrawLatex(0.15,0.87,hname.c_str());
    tex->DrawLatex(0.15,0.83,oneTrig.c_str());

    gPad->RedrawAxis();
    gPad->Update();

    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_Run%s.pdf",
		       hname.c_str(), doPull ? "nompull" : "nomrel", run.c_str()));


    const char *c1absname = Form("c1abs_%s",hname.c_str());
    TH1D *h1abs = tdrHist(Form("h1abs_%s_%s",oneTrig.c_str(),hname.c_str()),
			  "#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
    if (doPull) h1nom->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h1nom->GetZaxis()->SetRangeUser(-0.50,0.50);
    //if (doPull) h1abs->GetZaxis()->SetRangeUser(0,5*ntrg2);
    //else        h1abs->GetZaxis()->SetRangeUser(0.,0.50);

    TCanvas *c1abs = tdrCanvas(c1absname,h1abs,8,11,kRectangular);
    c1abs->SetRightMargin(0.15);
    h2abssum->Draw("COLZ SAME");
    if (doPull) h2abssum->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h2abssum->GetZaxis()->SetRangeUser(-0.5,0.5);
    //if (doPull) h2abssum->GetZaxis()->SetRangeUser(0,5*ntrg2);
    //else        h2abssum->GetZaxis()->SetRangeUser(0,0.5);

    tex->DrawLatex(0.15,0.87,hname.c_str());
    tex->DrawLatex(0.15,0.83,oneTrig.c_str());

    gPad->RedrawAxis();
    gPad->Update();

    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_Run%s.pdf",
		       hname.c_str(), doPull ? "abspull" : "absrel",
		       run.c_str()));


    
    if (!h2nomrefsum && !h2absrefsum && hname=="p2asymm") {
      h2nomrefsum = (TH2D*)h2nomsum->Clone("h2nomrefsum");
      h2absrefsum = (TH2D*)h2abssum->Clone("h2absrefsum");
    }

    if (!h2nomsums && !h2abssums) {
      h2nomsums = (TH2D*)h2nomsum->Clone("h2nomsums");
      h2abssums = (TH2D*)h2abssum->Clone("h2abssums");
    }
    else {
      h2nomsums->Add(h2nomsum);
      h2abssums->Add(h2abssum);
    }
  } // for ih

  TH1D *h1nom = tdrHist("h1nom","#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
  TH1D *h1abs = tdrHist("h1abs","#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
  if (doPull) {
    h1nom->GetZaxis()->SetRangeUser(-100,100);
    h1abs->GetZaxis()->SetRangeUser(0,100);
    //if (run=="2024BCD") {
    //h1nom->GetZaxis()->SetRangeUser(-150,150);
    //h1abs->GetZaxis()->SetRangeUser(0,150);
    //}
  }
  else {
    h1nom->GetZaxis()->SetRangeUser(-0.50,0.50);
    h1abs->GetZaxis()->SetRangeUser(0.,0.50);
  }

  TCanvas *c1nom = tdrCanvas("c1nom",h1nom,8,11,kRectangular);
  gPad->SetRightMargin(0.15);
  h2nomsums->Draw("COLZ SAME");

  TCanvas *c1abs = tdrCanvas("c1abs",h1abs,8,11,kRectangular);
  gPad->SetRightMargin(0.15);
  h2abssums->Draw("COLZ SAME");

  if (doPull) {
    h2nomsums->GetZaxis()->SetRangeUser(-100,100);
    h2abssums->GetZaxis()->SetRangeUser(-100,100);
    //h2abssums->GetZaxis()->SetRangeUser(0,100);
    //if (run=="2024BCD") {
    //h2nomsums->GetZaxis()->SetRangeUser(-150,150);
    //h2abssums->GetZaxis()->SetRangeUser(-150,150);
    //}
  }
  else {
    h2nomsums->GetZaxis()->SetRangeUser(-0.50,0.50);
    h2abssums->GetZaxis()->SetRangeUser(-0.50,0.50);
    //h2abssums->GetZaxis()->SetRangeUser(0,0.50);
  }

  
  TH2D *h2veto = (TH2D*)h2nomsums->Clone("jetvetomap"); // default map
  TH2D *h2hot = (TH2D*)h2nomsums->Clone("jetvetomap_hot");
  TH2D *h2cold = (TH2D*)h2nomsums->Clone("jetvetomap_cold");
  TH2D *h2old = (TH2D*)h2nomsums->Clone("jetvetomap_old");
  TH2D *h2hotandcold = (TH2D*)h2nomsums->Clone("jetvetomap_hotandcold");
  TH2D *h2eep = (TH2D*)h2nomsums->Clone("jetvetomap_eep");
  TH2D *h2bpix = (TH2D*)h2nomsums->Clone("jetvetomap_bpix");
  TH2D *h2all = (TH2D*)h2nomsums->Clone("jetvetomap_all");
  h2veto->Reset();
  h2hot->Reset();
  h2cold->Reset();
  h2old->Reset();
  h2hotandcold->Reset();
  h2eep->Reset();
  h2bpix->Reset();
  h2all->Reset();

  h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold");
  h2hot ->SetTitle("Hot zones. Use for steep jet pT spectra and MET tails");
  h2cold->SetTitle("Cold zones. Mostly ECAL holes. Use for MET tails");
  h2old->SetTitle("Old zones. Already in previous map.");
  h2hotandcold->SetTitle("Union of hot and cold zone maps");
  h2eep->SetTitle("EE+ water leak region. Complete loss of ECAL energy in 2022EFG");
  h2bpix->SetTitle("Barrel pixel failure. Reduction of tracking efficiency in 2023D and later");
  h2all->SetTitle("Union of hot and cold maps");
  if (run=="2022E" || run=="2022F" || run=="2022G" ||
      run=="2022EF" || run=="2022EFG") {
    h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold+EEP");
    h2all->SetTitle("Union of hot, cold and EEP maps");
  }
  if (run=="2023D") {
    h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold+BPIX");
    h2all->SetTitle("Union of hot, cold and BPIX maps");
  }
  if (run=="2024BCD") {
    h2veto->SetTitle("JME recommended map. Hot+Cold(+not BPIX)");
    h2all->SetTitle("Union of all maps, used for JEC. Hot+cold+BPIX");
  }
  h2nomsums->SetTitle("Raw nominal pull map. NomPull=(x_i-mu)/sigma. Summed over triggers");
  h2abssums->SetTitle("Raw absolute pull map. AbsPull=|(x_i-mu)/sigma|-1. Summed over triggers");

  for (int i = 1; i != h2nomsums->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2nomsums->GetNbinsY()+1; ++j) {
      double eta = h2nomsums->GetXaxis()->GetBinCenter(i);
      double phi = h2nomsums->GetYaxis()->GetBinCenter(j);

      // Match positive asymmetry with hot region
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (!h2nomrefsum || h2nomrefsum->GetBinContent(i,j)>0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2hot->SetBinContent(i,j,100);
      }

      // Match negative asymmetry with cold region
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (h2nomrefsum && h2nomrefsum->GetBinContent(i,j)<0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2cold->SetBinContent(i,j,-100);
      }

      // Match zero asymmetry with old region (probably had p2asymm masked)
      if ((h2abssums->GetBinContent(i,j)>pullThreshold && fabs(eta)<4.5 ||
	   h2abssums->GetBinContent(i,j)>pullThresholdHF45 && fabs(eta)>4.5) &&
	  (h2nomrefsum && h2nomrefsum->GetBinContent(i,j)==0)) {
	h2veto->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2old->SetBinContent(i,j,-100);
      }

      // Extra manual margin for EE+ water leak region (2022E,F,G)
      if (eta>1.5 && phi>1.85 &&
	  eta<2.2 && phi<2.7 &&
	  phi<2.7-(2.7-2.1)/(2.1-1.6)*(eta-1.6)) {
	if (run=="2022E" || run=="2022F" || run=="2022G" ||
	    run=="2022EF" || run=="2022EFG") {
	  h2veto->SetBinContent(i,j,100);
	  h2eep->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
      } // EEP

      // BPIX reference region (2023D and later)
      if (eta>-1.5 && phi>-1.22 &&
	  eta<0.1 && phi<-0.87) {
	if (run=="2023D") { // Remove BPix
	  h2veto->SetBinContent(i,j,100);
	  h2bpix->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
	if (run=="2024BCD") { // Keep BPix for recommended, drop for JEC
	  h2veto->SetBinContent(i,j,0);
	  h2bpix->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	}
      } // BPIX

    } // for j
  } // for j

  if (plotJetVeto) {
    h2veto->SetLineColor(kRed);
    h2hotandcold->SetLineColor(kOrange+1);
    h2hot->SetLineColor(kRed);

    c1nom->cd();
    h2veto->Draw("SAME BOX");
    h2hotandcold->Draw("SAME BOX");
    h2hot->Draw("SAME BOX");

    c1abs->cd();
    h2veto->Draw("SAME BOX");
    h2hotandcold->Draw("SAME BOX");
    h2hot->Draw("SAME BOX");
  }
  gPad->RedrawAxis();
  gPad->Update();
  
  if (oneTrig!="" && oneHist!="") {
    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_%s_Run%s.pdf",
		       oneTrig.c_str(), oneHist.c_str(),
		       doPull ? "nompull" : "absrel",run.c_str()));
    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_%s_%s_%s_Run%s.pdf",
		       oneTrig.c_str(), oneHist.c_str(),
		       doPull ? "abspull" : "absrel",run.c_str()));
  }
  else {
    c1nom->SaveAs(Form("pdf/JetVeto/JetVeto_h2nommap_Run%s.pdf",run.c_str()));
    c1abs->SaveAs(Form("pdf/JetVeto/JetVeto_h2absmap_Run%s.pdf",run.c_str()));
  }

  TH2D *h2jesnorm(0);
  if (h2jes) {
    TH1D *h2 = tdrHist("h2","#phi",-TMath::Pi(),+TMath::Pi(),"#eta",-5.2,5.2);
    h2->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c2 = tdrCanvas("c2",h2,8,11,kRectangular);
    gPad->SetRightMargin(0.15);
    h2jes->Draw("SAME COLZ");
    h2jes->GetZaxis()->SetRangeUser(-0.15,0.15);
    gPad->RedrawAxis();
    gPad->Update();

    if (plotJetVetoOnJES) {
      h2veto->SetLineColor(kRed);
      h2veto->Draw("SAME BOX");
    }
      
    c2->SaveAs(Form("pdf/JetVeto/JetVeto_x2jes_%s.pdf",run.c_str()));

    h2jes->SetTitle("Dijet asymmetry map, (pTprobe-pTtag)/pTave for |eta,tag|<1.3");


    // Normalize eta strips
    h2jesnorm = (TH2D*)h2jes->Clone("h2jesnorm");
    for (int i = 1; i != h2jes->GetNbinsX()+1; ++i) {
      double norm = h2jes->Integral(i,i,1,72) / 72.;
      for (int j = 1; j != h2jes->GetNbinsY()+1; ++j) {
	h2jesnorm->SetBinContent(i,j,(1+h2jes->GetBinContent(i,j))/(1+norm)-1);
	h2jesnorm->SetBinError(i,j,h2jes->GetBinError(i,j)/(1+norm));
      }
    }
    
    TH1D *h2norm = tdrHist("h2norm","#phi",-TMath::Pi(),+TMath::Pi(),
			   "#eta",-5.2,5.2);
    h2norm->GetZaxis()->SetRangeUser(-0.50,0.50);

    TCanvas *c2norm = tdrCanvas("c2norm",h2norm,8,11,kRectangular);
    gPad->SetRightMargin(0.15);
    h2jesnorm->Draw("SAME COLZ");
    h2jesnorm->GetZaxis()->SetRangeUser(-0.15,0.15);
    gPad->RedrawAxis();
    gPad->Update();

    if (plotJetVetoOnJESnorm) {
      h2veto->SetLineColor(kRed);
      h2veto->Draw("SAME BOX");
    }
      
    c2norm->SaveAs(Form("pdf/JetVeto/JetVeto_x2jesnorm_%s.pdf",run.c_str()));

    h2jesnorm->SetTitle("Dijet asymmetry map, (pTprobe-pTtag)/pTave for |eta,tag|<1.3 normalized vs #phi_{probe}");
  } // if (h2jes)

  fout->cd();
  h2veto->Write("jetvetomap");
  h2hot->Write("jetvetomap_hot");
  h2cold->Write("jetvetomap_cold");
  h2old->Write("jetvetomap_old");
  h2hotandcold->Write("jetvetomap_hotandcold");
  if (h2eep->Integral()!=0) h2eep->Write("jetvetomap_eep");
  if (h2bpix->Integral()!=0) h2bpix->Write("jetvetomap_bpix");
  h2all->Write("jetvetomap_all");
  if (h2jes) h2jes->Write("jetasymmetrymap");
  if (h2jesnorm) h2jesnorm->Write("jetasymmetrymap_norm");
  h2nomsums->Write("jetpullsummap_nom");
  h2abssums->Write("jetpullsummap_abs");
  if (h2nomrefsum) h2nomrefsum->Write("jetpullsummap_ref");
  fout->Write();
  fout->Close();
} // JetVeto
