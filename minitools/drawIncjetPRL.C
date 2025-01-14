// Purpose: Draw two PRL plots for inclusive jet xsec at 13.6 TeV (late 2024)
#include "TFile.h"
#include "TH1D.h"

#include "../tdrstyle_mod22.C"

// Normalize histogram for pT and lumi
void normalizeHist(TH1D *h, double absy);
void cleanHist(TH1D *h, double xmin, double xmax);

void drawIncjetPRL() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  TFile *f = new TFile("rootfiles/Prompt2024/v113_2024/jmenano_data_cmb_2024FGHI_JME_v113.root","READ");
  assert(f && !f->IsZombie());
  curdir->cd();

  //TH1D *h = tdrHist("h","Cross section d#sigma/dp_{T}dy (pb/GeV)",1e-8,1e11,
  //		    "Jet p_{T} (GeV)",10-1e-3,5000);
  TH1D *h = tdrHist("h","Cross section d#sigma/dp_{T}dy (pb/GeV)",
		    1.001e-8,0.999e8,
		    "Jet p_{T} (GeV)",50,5000);
  extraText = "Private";
  lumi_136TeV = "Late 2024, 82.4 fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,8,11,kSquare);
  h->GetXaxis()->SetMoreLogLabels(kFALSE);
  c1->SetLogx();
  c1->SetLogy();

  const int nybins(10);
  int marker[nybins] =
    {kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullDiamond,
     kFullSquare, kOpenSquare, kFullCircle, kOpenCircle, kFullDiamond};
  int color[nybins] =
    {kBlack, kBlack, kBlack, kBlack, kBlack,
     kRed, kBlue, kBlue, kBlue, kBlue};
  double size[nybins] =
    {0.5, 0.5, 0.5, 0.5, 0.8,
     0.5, 0.5, 0.5, 0.5, 0.8};

  TLegend *leg = tdrLeg(0.85-0.25,0.90-0.035*nybins,0.85,0.90);
  leg->SetTextSize(0.040);
  
  for (int ybin = 0; ybin != nybins; ++ybin) {

    int ieta = 5*(ybin+1);
    double y1 = 0+ybin*0.5;
    double y2 = 0.5+ybin*0.5;
    double y = 0.5*(y1+y2);
    
    TH1D *hpt = (TH1D*)f->Get(Form("Incjet/hpt%02d",ieta)); assert(hpt);
    normalizeHist(hpt, y);
    cleanHist(hpt, 50., 13000.);

    tdrDraw(hpt,"HP][",marker[ybin],color[ybin],kSolid,-1,kNone,0,size[ybin]);

    if (y1==0) leg->AddEntry(hpt,Form("|y| < %1.1f",y2));
    else leg->AddEntry(hpt,Form("%1.1f #leq |y| < %1.1f",y1,y2));
  } // for ybin

  gPad->RedrawAxis();
  c1->SaveAs("pdf/drawIncjetPRL/drawIncjetPRL_spectra.pdf");
}


void normalizeHist(TH1D *h, double absy) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {

    double x = h->GetBinCenter(i);
    double ex = h->GetBinWidth(i);
    double y = h->GetBinContent(i);
    double ey = h->GetBinError(i);

    /*
    mi["HLT_PFJet40"]  = range{49,  84,  0, fwdeta0};
    mi["HLT_PFJet60"]  = range{84,  114, 0, fwdeta};
    mi["HLT_PFJet80"]  = range{114, 196, 0, fwdeta};
    mi["HLT_PFJet140"] = range{196, 272, 0, fwdeta};
    mi["HLT_PFJet200"] = range{272, 330, 0, fwdeta0};
    mi["HLT_PFJet260"] = range{330, 395, 0, fwdeta0};
    mi["HLT_PFJet320"] = range{395, 468, 0, fwdeta0};
    mi["HLT_PFJet400"] = range{468, 548, 0, fwdeta0};
    mi["HLT_PFJet450"] = range{548, 686, 0, fwdeta0};
    mi["HLT_PFJet500"] = range{686,6500, 0, fwdeta0};

    double fwdeta = 3.139;  // was 2.853. 80% (100%) on negative (positive) side
    double fwdeta0 = 2.964; // 2.853; // 40 and 260 up
    mi["HLT_PFJetFwd40"]  = range{49,  84,  fwdeta0, 5.2};
    mi["HLT_PFJetFwd60"]  = range{84,  114, fwdeta, 5.2};
    mi["HLT_PFJetFwd80"]  = range{114, 196, fwdeta, 5.2};
    mi["HLT_PFJetFwd140"] = range{196, 272, fwdeta, 5.2};
    mi["HLT_PFJetFwd200"] = range{272, 330, fwdeta0, 5.2};
    mi["HLT_PFJetFwd260"] = range{330, 395, fwdeta0, 5.2};
    mi["HLT_PFJetFwd320"] = range{395, 468, fwdeta0, 5.2};
    mi["HLT_PFJetFwd400"] = range{468, 548, fwdeta0, 5.2};
    mi["HLT_PFJetFwd450"] = range{548, 686, fwdeta0, 5.2};
    mi["HLT_PFJetFwd500"] = range{686,6500, fwdeta0, 5.2};
    */
    
    double lum(82438.188401497); // 2024FGHI_all.txt online
    // Prescale trigger luminosities calculated with brilcalc, e.g.:
    // brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json --datatag online -u /pb --begin 382229 -i Cert_Collisions2024_378981_386951_Golden.json --hltpath "HLT_PFJet500_v*" > 2024FGHI_500.txt
    if (absy<3.0) {
      if (x>686) { // HLT_PFJet500_v*
	lum = 82419.755946457;
      } else if (x>548) { // HLT_PFJet450_v*
	lum = 10302.469493307;
      } else if (x>468) { // HLT_PFJet400_v*
	lum = 5151.234746654;
      } else if (x>395) { // HLT_PFJet320_v*
	lum = 2575.617373327;
      } else if (x>330) { // HLT_PFJet260_v*
	lum = 643.904343332;
      } else if (x>272) { // HLT_PFJet200_v*
	lum = 235.485016990;
      } else if (x>196) { // HLT_PFJet140_v*
	lum = 54.946503964;
      } else if (x>114) {  // HLT_PFJet80_v*
	lum = 4.488592109;
      } else if (x>84) { // HLT_PFJet60_v*
	lum = 1.180534644;
      } else if (x>49) { // HLT_PFJet40_v*
	lum = 0.157404619;
      } else if (x>0) { // HLT_ZeroBias_v*
	lum = 0.157432157;
      }
    }
    else {
      if (x>686) { // HLT_PFJetFwe500_v*
	lum = 82419.755946457;//82419.755946457;
      } else if (x>548) { // HLT_PFJetFwd450_v*
	lum = 82419.755946457;//10302.469493307;
      } else if (x>468) { // HLT_PFJetFwd400_v*
	lum = 82419.755946457;//5151.234746654;
      } else if (x>395) { // HLT_PFJetFwd320_v*
	lum = 20604.938986614;//2575.617373327;
      } else if (x>330) { // HLT_PFJetFwd260_v*
	lum = 10302.469493307;//643.904343332;
      } else if (x>272) { // HLT_PFJetFwd200_v*
	lum = 1831.550132143;//235.485016990;
      } else if (x>196) { // HLT_PFJetFwd140_v*
	lum = 312.196045252;//54.946503964;
      } else if (x>114) {  // HLT_PFJetFwd80_v*
	lum = 15.095193397;//4.488592109;
      } else if (x>84) { // HLT_PFJetFwd60_v*
	lum = 2.361069288;//1.180534644;
      } else if (x>49) { // HLT_PFJetFwd40_v*
	lum = 0.190857498;//0.157404619;
      } else if (x>0) { // HLT_ZeroBias_v*
	lum = 0.157432157;
      }
    }
      
    h->SetBinContent(i, y / ex / lum);
    h->SetBinError(i, ey / ex / lum);
  } // for i
  
} // normalizeHist


void cleanHist(TH1D *h, double xmin, double xmax) {

  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    double x = h->GetBinCenter(i);
    if (x<xmin || x>xmax) {
      h->SetBinContent(i, 0.);
      h->SetBinError(i, 0.);
    }
  } // for i
} // cleanHist
