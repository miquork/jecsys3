// Purpose: draw jetvetomaps for Run3
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TBox.h"

#include "../tdrstyle_mod22.C"

// Set negative values to zero in TH2D for better "BOX" drawing
void rezero(TH2D* h2, double thr=0, double max=10) {
  
  assert(h2);

  for (int i = 1; i != h2->GetNbinsX()+1; ++i) {
    for (int j = 1; j != h2->GetNbinsY()+1; ++j) {
      if (h2->GetBinContent(i,j)<thr) h2->SetBinContent(i,j,0);
      if (h2->GetBinContent(i,j)>max) h2->SetBinContent(i,j,max);
    }
  }
} // void rezero


void hotjetsRun3() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  //string eras[] = {"22CD", "22EFG", "23BC", "23D", "24BCD", "24E"};
  //string eras[] = {"22CD", "22EFG", "23BC", "23D", "24BCDE", "24F"};
  //string eras[] = {"22CD", "22EFG", "23BC", "23D", "24BCDE", "24FG"};
  //string eras[] = {"22CD", "22EFG", "23BC", "23D", "24BCDE", "24BCDEFGHI"};
  string eras[] = {"22CD", "22EFG", "23BC", "23D", "24BCDEFGHI", "25C"};
  const int nera = sizeof(eras)/sizeof(eras[0]);

  map<string, const char*> mera;
  mera["22CD"] = "22";
  mera["22EFG"] = "22EE";
  mera["23BC"] = "23";
  mera["23D"] = "23BPix";
  //mera["24BCD"] = "24";
  //mera["24BCDE"] = "24BCDE";
  mera["24BCDE"] = "24B-E";
  //mera["24E"] = "24E";
  //mera["24F"] = "24F";
  mera["24FG"] = "24FG";
  mera["24BCDEFGHI"] = "24";
  mera["25C"] = "25C";

  map<string, int> mcolor;
  mcolor["22CD"] = kGreen+2;
  mcolor["22EFG"] = kBlue;
  mcolor["23BC"] = kMagenta+2;//kOrange;
  mcolor["23D"] = kCyan+2;
  //mcolor["24BCD"] = kRed;
  mcolor["24BCDE"] = kOrange +1;
  //mcolor["24E"] = kYellow+2;
  //mcolor["24F"] = kRed;
  mcolor["24FG"] = kRed;
  //mcolor["24BCDEFGHI"] = kRed;
  mcolor["24BCDEFGHI"] = kSpring;
  mcolor["25C"] = kRed;

  TH1D *h = new TH1D("h",";#eta_{jet};#phi_{jet}",100,-4.7,4.7);
  h->SetMaximum(+TMath::Pi());
  h->SetMinimum(-TMath::Pi());

  //lumi_13TeV = "Run2024BC, 0.74 fb^{-1}";
  //lumi_13TeV = "Run2024BC, 3.3 fb^{-1}";
  //lumi_13TeV = "Run2024BCD, 12.3 fb^{-1}";
  //lumi_13TeV = "Run2024BCDE, 27.0 fb^{-1}";
  //lumi_13TeV = "Run2024BCDEF, 46.0 fb^{-1}";
  //lumi_13TeV = "Run2024BCDEFG, 92.2 fb^{-1}";
  //lumi_13TeV = "Run2024B-I, 109.2 fb^{-1}";
  lumi_13TeV = "Run2025C, X.X fb^{-1}";
  TCanvas *c1 = tdrCanvas("c1",h,4,0,kRectangular);

  TLine *l = new TLine();

  l->SetLineStyle(kSolid);
  double etahf = 2.964;
  l->DrawLine(-etahf,-TMath::Pi(),-etahf,+TMath::Pi());
  l->DrawLine(+etahf,-TMath::Pi(),+etahf,+TMath::Pi());
  l->SetLineStyle(kDashed);
  double etatr = 2.5;
  l->DrawLine(-etatr,-TMath::Pi(),-etatr,+TMath::Pi());
  l->DrawLine(+etatr,-TMath::Pi(),+etatr,+TMath::Pi());
  l->SetLineStyle(kDotted);
  double etaec = 1.305;
  l->DrawLine(-etaec,-TMath::Pi(),-etaec,+TMath::Pi());
  l->DrawLine(+etaec,-TMath::Pi(),+etaec,+TMath::Pi());

  //TLegend *leg = tdrLeg(0.15,0.90-nera*0.05,0.35,0.90);
  TLegend *leg = tdrLeg(0.13,0.90-(nera+1)*0.05,0.33,0.90);

  TH2D *h2sum(0);
  for (int i = 0; i != nera; ++i) {
    
    string era = eras[i];
    const char *ce = era.c_str();
    
    TFile *f = new TFile(Form("rootfiles/jetveto20%s.root",ce),"READ");
    assert(f && !f->IsZombie());

    curdir->cd();

    TH2D *h2 = (TH2D*)f->Get("jetvetomap"); assert(h2);

    rezero(h2);
    h2->GetZaxis()->SetRangeUser(-10,10);
    h2->SetLineColor(mcolor[era]);
    h2->SetFillStyle(1001);
    h2->SetFillColor(kNone);
    h2->SetFillColorAlpha(mcolor[era]-9, 0.35); // 65% transparent
    h2->Draw("SAMEBOX");

    //if (era=="24BCD")
    //h2->SetFillColorAlpha(mcolor[era], 0.60); // 40% transparent
    if (era=="24BCDE" || era=="24F")
      h2->SetFillColorAlpha(mcolor[era], 0.50); // 50% transparent
    if (era=="25C")
      h2->SetFillColorAlpha(mcolor[era], 0.35); // 60% transparent
    
    leg->AddEntry(h2,mera[era],"F");
    
    // combination of regions
    if (!h2sum) h2sum = (TH2D*)h2->Clone("h2sum");
    else h2sum->Add(h2);
  }
  
  //rezero(h2sum,4*10,10); // overlap min. 4
  rezero(h2sum,3*10,10); // overlap min. 3
  //rezero(h2sum,2*10,10); // overlap min. 2

  h2sum->SetLineColor(kBlack);
  h2sum->SetLineStyle(kNone);
  h2sum->SetFillStyle(1001);
  h2sum->SetFillColor(kNone);
  h2sum->DrawClone("SAMEBOX");

  //leg->AddEntry(h2sum,"Min. 4","F");
  leg->AddEntry(h2sum,"Min. 3","F");

  gPad->Paint();

  //c1->SaveAs("pdf/hotjetsRun3_V4M.pdf");
  //c1->SaveAs("pdf/hotjetsRun3_V5M.pdf");
  //c1->SaveAs("pdf/hotjetsRun3_V6M.pdf");
  //c1->SaveAs("pdf/hotjetsRun3_V7M.pdf");
  c1->SaveAs("pdf/hotjetsRun3_Prompt25_V1M.pdf");

  //TFile *fout = new TFile("rootfiles/hotjets-Run3.root","RECREATE");
  //h2sum->Write();
  //fout->Close();
}
