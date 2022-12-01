// Purpose: Produce jet veto maps from JMENANO analyzer on dijets
#include "TFile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TMath.h"

#include <string>
#include <vector>

#include "../tdrstyle_mod22.C"

// use pulls instead of relative changes
bool doPull = true;
bool plotJetVeto = true;
// If these are non-zero (not ""), will use just this
string oneTrig = "";
//string oneTrig = "HLT_PFJet500";
//string oneTrig = "HLT_PFJet450";
//string oneTrig = "HLT_PFJet320";
//string oneTrig = "HLT_PFJet40";
//string oneTrig = "HLT_PFJet60";
//string oneTrig = "HLT_PFJet80";
//string oneTrig = "HLT_PFJet140";
//string oneTrig = "HLT_PFJet260";
string oneHist = "";
//string oneHist = "p2asymm";
//string oneHist = "p2nhf";
//string oneHist = "p2nhf";

void doJetVeto(string run = "CD") {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f(0);
  if (run=="CD") {
    lumi_136TeV = "RunCD, 7.6 fb^{-1}";
    f = new TFile("../dijet/rootfiles/jmenano_data_out_v13cd.root","READ");
  }
  if (run=="E") {
    lumi_136TeV = "RunE, 5.5 fb^{-1}";
    f = new TFile("../dijet/rootfiles/jmenano_data_out_v13e.root","READ"); // new v6
  }
 ///TFile *f = new TFile("../dijet/rootfiles/jmenano_data_out_v12e.root","READ");
  //lumi_136TeV = "RunF, 1.4 fb^{-1}";
  //TFile *f = new TFile("../dijet/rootfiles/jmenano_data_out_v12f.root","READ");

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
  vtrg.push_back("HLT_PFJet550");

  vtrg.push_back("HLT_PFJetFwd40");
  vtrg.push_back("HLT_PFJetFwd60");
  vtrg.push_back("HLT_PFJetFwd80");
  // NB: these had wrong |eta| threshold in v12+v13. To be fixed in v14
  /*
  vtrg.push_back("HLT_PFJetFwd140");
  vtrg.push_back("HLT_PFJetFwd200");
  vtrg.push_back("HLT_PFJetFwd260");
  vtrg.push_back("HLT_PFJetFwd320");
  vtrg.push_back("HLT_PFJetFwd400");
  vtrg.push_back("HLT_PFJetFwd450");
  vtrg.push_back("HLT_PFJetFwd500");
  */

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
  vh.push_back("h2pt");
  vh.push_back("p2nef");
  //vh.push_back("p2chf");
  //vh.push_back("p2nhf");
  if (oneHist!="") {
    vh.clear();
    vh.push_back(oneHist);
  }
  int nh = vh.size();

  TH2D *h2sums(0), *h2refsum(0), *h2jes(0);
  for (int ih = 0; ih != nh; ++ih) {

    string hname = vh[ih];
    TH2D *h2sum(0);
  
    for (int itrg = 0; itrg != ntrg; ++itrg) {

      string trg = vtrg[itrg];
      //string trg = "HLT_PFJet500";
      //string trg = "HLT_ZeroBias";
      //string hname = "p2nef";
      //isProfile2D = (hname[0]='p' && hname[1]='2');
      
      TObject *obj = f->Get(Form("%s/Jetveto/%s",trg.c_str(),hname.c_str()));
      assert(obj);
      assert(obj->InheritsFrom("TH2D"));
      bool isProf2D = obj->InheritsFrom("TProfile2D");
      
      TH2D *h2 = (isProf2D ? ((TProfile2D*)obj)->ProjectionXY() : (TH2D*)obj);

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
	//int i = int(h2->GetNbinsX()/2.); {
	
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
	if (n<70) {
	  for (int j = 1; j != ny+1; ++j) {
	    h2->SetBinContent(i, j, 0.);
	    h2->SetBinError(i, j, 0.);
	  } // for j
	}
	// 
	else {
	  // Determine median and 68% quantile ([0.16,0.84] range divided by 2)
	  const Int_t nprobSum = 3;//5;
	  //Double_t probSum[nprobSum] = {0.025, 0.16, 0.50, 0.84, 0.975};
	  Double_t probSum[nprobSum] = {0.16, 0.50, 0.84};
	  Double_t q[nprobSum];
	  htmp->GetQuantiles(nprobSum, q, probSum);
	  double median = q[1];//q[2];
	  double q68 = 0.5*(q[2]-q[0]);
	  //double q68 = 0.5*(q[3]-q[1]);
	  //double q95 = 0.25*(q[4]-q[0]);
	  //cout << "Bin " << i << ", eta=" << h2->GetXaxis()->GetBinCenter(i)<<endl;
	  //cout << "Mean = " << htmp->GetMean()
	  //   << ", RMS = " << htmp->GetRMS() << endl;
	  //cout << "Median = " << median << ", q68% = " << q68 << endl;
	  //<< ", q95% = " << q95 << endl;
	  
	  // Determine pull width
	  for (int j = 1; j != ny+1; ++j) {
	    if (h2->GetBinError(i,j)!=0 && q68!=0) {
	      //if (h2->GetBinError(i, j) < q68) {
	      if (doPull) {
		h2->SetBinContent(i, j, (h2->GetBinContent(i, j)-median) / q68);
		h2->SetBinError(i, j, h2->GetBinError(i, j) / q68);
	      }
	      else {
		//if (fabs(median)<0.25) {
		if (median<0.5) {
		  h2->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/
				    (1+fabs(median)));
		  h2->SetBinError(i, j, h2->GetBinError(i, j) /
				  (1+fabs(median)));
		}
		else {
		  h2->SetBinContent(i,j,(h2->GetBinContent(i,j)-median)/median);
		  h2->SetBinError(i, j, h2->GetBinError(i, j) / median);
		}
	      }
	      //}
	      //else {
	      //h2->SetBinContent(i, j, 0.);
	      //h2->SetBinError(i, j, 0.);
	      //}
	    }
	  } // for j
	}
	
	delete htmp;
      } // for i
      
      //h2->Draw("COLZ");
      if (!h2sum) {
	h2sum = (TH2D*)h2->Clone(Form("h2sum_%s",hname.c_str()));
      }
      else {
	h2sum->Add(h2);
      }
      delete h2;
		
    } // for itrg

    double ntrg2 = (oneTrig!="" ? 1 : max(1, ntrg/2));
    const char *c1name = Form("c1_%s",hname.c_str());
    //TCanvas *c1 = new TCanvas(c1name,c1name,800,600);
    TH1D *h1 = tdrHist(Form("h1_%s_%s",oneTrig.c_str(),hname.c_str()),
		       "#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
    if (doPull) h1->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h1->GetZaxis()->SetRangeUser(-0.50,0.50);
    //lumi_136TeV = "RunE, 5.5 fb^{-1}";
    TCanvas *c1 = tdrCanvas(c1name,h1,8,11,kRectangular);
    c1->SetRightMargin(0.15);
    h2sum->Draw("COLZ SAME");
    if (doPull) h2sum->GetZaxis()->SetRangeUser(-5*ntrg2,5*ntrg2);
    else        h2sum->GetZaxis()->SetRangeUser(-0.5,0.5);
    //h2sum->Draw("COLZ");

    TLatex *tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.035);
    tex->DrawLatex(0.15,0.87,hname.c_str());//oneHist.c_str());
    tex->DrawLatex(0.15,0.83,oneTrig.c_str());

    gPad->RedrawAxis();
    gPad->Update();

    c1->SaveAs(Form("pdf/doJetVeto/doJetVeto_%s_%s_Run%s.pdf",
		    hname.c_str(), doPull ? "pull" : "rel", run.c_str()));
    
    if (!h2refsum && hname=="p2asymm") {
      h2refsum = (TH2D*)h2sum->Clone("h2refsum");
    }

    if (!h2sums) {
      h2sums = (TH2D*)h2sum->Clone("h2sums");
    }
    else {
      h2sums->Add(h2sum);
    }
  } // for ih

  TH1D *h1 = tdrHist("h1","#phi",-TMath::Pi(),TMath::Pi(),"#eta",-5.2,5.2);
  if (doPull)
    h1->GetZaxis()->SetRangeUser(-100,100);
  else
    h1->GetZaxis()->SetRangeUser(-0.50,0.50);
  //lumi_136TeV = "RunE, 5.5 fb^{-1}";
  //TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TCanvas *c1 = tdrCanvas("c1",h1,8,11,kRectangular);
  gPad->SetRightMargin(0.15);
  h2sums->Draw("COLZ SAME");
  if (doPull)
    h2sums->GetZaxis()->SetRangeUser(-100,100);
  else
    h2sums->GetZaxis()->SetRangeUser(-0.50,0.50);
  //h2sums->GetZaxis()->SetRangeUser(-50,5);
  //h2sums->Draw("COLZ");

  TH2D *h2veto = (TH2D*)h2sums->Clone("jetvetomap"); // default map
  TH2D *h2hot = (TH2D*)h2sums->Clone("jetvetomap_hot");
  TH2D *h2cold = (TH2D*)h2sums->Clone("jetvetomap_cold");
  //TH2D *h2both = (TH2D*)h2sums->Clone("jetvetomap_both");
  TH2D *h2hotandcold = (TH2D*)h2sums->Clone("jetvetomap_hotandcold");
  TH2D *h2eep = (TH2D*)h2sums->Clone("jetvetomap_eep");
  TH2D *h2all = (TH2D*)h2sums->Clone("jetvetomap_all");
  h2veto->Reset();
  h2hot->Reset();
  h2cold->Reset();
  h2hot->Reset();
  //h2both->Reset();
  h2hotandcold->Reset();
  h2eep->Reset();
  h2all->Reset();

  h2veto->SetTitle("JME recommended map, used for JEC. Hot+Cold+EEP");
  h2hot ->SetTitle("Hot zones. Use for steep jet pT spectra and MET tails");
  h2cold->SetTitle("Cold zones. Mostly ECAL holes. Use for MET tails");
  h2hotandcold->SetTitle("Union of hot and cold zone maps");
  h2eep->SetTitle("EE+ water leak region. Complete loss of ECAL energy");
  h2all->SetTitle("Union of hot, cold and EEP maps");
  h2sums->SetTitle("Raw pull map. Pull=(x_i-mu)/sigma. Summed over triggers");

  for (int i = 1; i != h2sums->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2sums->GetNbinsY()+1; ++j) {
      double eta = h2sums->GetXaxis()->GetBinCenter(i);
      double phi = h2sums->GetYaxis()->GetBinCenter(j);
      //if ((h2sums->GetBinContent(i,j)>80 && fabs(eta)<2.5) || // v3-4
	  //h2sums->GetBinContent(i,j)>100) { // v3
	  //h2sums->GetBinContent(i,j)>50) { // looser for hot regions, v4
      if (h2sums->GetBinContent(i,j)>70) { // v5+
	h2veto->SetBinContent(i,j,100);
	//h2both->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2hot->SetBinContent(i,j,100);
      }
      //if ((h2sums->GetBinContent(i,j)<-50 && fabs(eta)<2.5) || // v3-4
      //  h2sums->GetBinContent(i,j)<-100) { // v3-4
      if (h2sums->GetBinContent(i,j)<-80) { // v5+
	h2veto->SetBinContent(i,j,100);
	//h2both->SetBinContent(i,j,100);
	h2hotandcold->SetBinContent(i,j,100);
	h2all->SetBinContent(i,j,100);
	h2cold->SetBinContent(i,j,-100);
      }
      // extra manual margin for EE+ water leak region
      if (eta>1.5 && phi>1.85 &&
	  eta<2.2 && phi<2.7 &&
	  phi<2.7-(2.7-2.1)/(2.1-1.6)*(eta-1.6)) {
	if (run=="E" || run=="F") {
	  h2veto->SetBinContent(i,j,100);
	  h2all->SetBinContent(i,j,100);
	  h2eep->SetBinContent(i,j,100);
	}
      }
    } // for j
  } // for j
  if (plotJetVeto) {
    h2veto->SetLineColor(kRed);
    h2veto->Draw("SAME BOX");
    //h2both->SetLineColor(kOrange+1);
    //h2both->Draw("SAME BOX");
    h2hotandcold->SetLineColor(kOrange+1);
    h2hotandcold->Draw("SAME BOX");
    h2hot->SetLineColor(kRed);
    h2hot->Draw("SAME BOX");
  }
  gPad->RedrawAxis();
  gPad->Update();
  if (oneTrig!="" && oneHist!="") {
    c1->SaveAs(Form("pdf/doJetVeto/doJetVeto_%s_%s_%s_Run%s.pdf",
		    oneTrig.c_str(), oneHist.c_str(),
		    doPull ? "pull" : "rel",run.c_str()));
  }
  else
    c1->SaveAs(Form("pdf/doJetVeto/doJetVeto_h2map_Run%s.pdf",run.c_str()));

  //h2eep->SetLineColor(kOrange+1);
  //h2eep->Draw("SAME BOX");

  if (h2jes) {
    TH1D *h2 = tdrHist("h2","#phi",-TMath::Pi(),+TMath::Pi(),"#eta",-5.2,5.2);
    h2->GetZaxis()->SetRangeUser(-0.50,0.50);
    //TCanvas *c2 = new TCanvas("c2","c2",800,600);
    TCanvas *c2 = tdrCanvas("c2",h2,8,11,kRectangular);
    gPad->SetRightMargin(0.15);
    h2jes->Draw("SAME COLZ");
    h2jes->GetZaxis()->SetRangeUser(-0.50,0.50);
    //h2jes->Draw("COLZ");
    gPad->RedrawAxis();
    gPad->Update();
    c2->SaveAs(Form("pdf/doJetVeto/doJetVeto_h2jes_%s.pdf",run.c_str()));

    h2jes->SetTitle("Dijet asymmetry map, (pTprobe-pTtag)/pTave for |eta,tag|<1.3");
  }

  //TFile *fout = new TFile("rootfiles/jetveto2022E.root","RECREATE");
  TFile *fout = new TFile(Form("rootfiles/jetveto2022%s.root",run.c_str()),
			  "RECREATE");
  fout->cd();
  h2veto->Write("jetvetomap");
  h2hot->Write("jetvetomap_hot");
  h2cold->Write("jetvetomap_cold");
  //h2both->Write("jetvetomap_both");
  h2hotandcold->Write("jetvetomap_hotandcold");
  h2eep->Write("jetvetomap_eep");
  h2all->Write("jetvetomap_all");
  h2jes->Write("jetasymmetrymap");
  h2sums->Write("jetpullsummap");
  fout->Write();
  fout->Close();
} // doJetveto
