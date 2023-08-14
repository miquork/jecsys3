// Purpose: Plot jet veto maps from minitools/doJetVetoV2.C
//          Used for Run3 plots (2022CD,EF(G), 2023BC,D)
//          based on jecsys/hotjets-Run2.C
// Author:  Mikko Voutilainen (at) cern (dot) ch
// Date:    2023-08-14
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TBox.h"

#include "../tdrstyle_mod22.C"

bool excludeHEP17 = false;

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


void drawJetVetoV2() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *f22cd = new TFile("rootfiles/jetveto2022CD.root","READ");
  assert(f22cd && !f22cd->IsZombie());

  TFile *f22ef = new TFile("rootfiles/jetveto2022EF.root","READ");
  assert(f22ef && !f17->IsZombie());

  TFile *f23bc = new TFile("rootfiles/jetveto2023BC.root","READ");
  assert(f23bc && !f23bc->IsZombie());

  TFile *f23d = new TFile("rootfiles/jetveto2023D.root","READ");
  assert(f23d && !f23d->IsZombie());


  curdir->cd();

  TH2D *h2_22cd = (TH2D*)f22cd->Get("jetvetomap"); assert(h2_22cd);
  TH2D *h2_22ef = (TH2D*)f22ef->Get("jetvetomap"); assert(h2_22ef);
  TH2D *h2_23bc = (TH2D*)f23bc->Get("jetvetomap"); assert(h2_23bc);
  TH2D *h2_23d  = (TH2D*)f23d->Get("jetvetomap"); assert(h2_23d);

  TH1D *h = new TH1D("h",";#eta_{jet};#phi_{jet}",100,-4.7,4.7);
  h->SetMaximum(+TMath::Pi());
  h->SetMinimum(-TMath::Pi());

  lumi_13TeV = "Run3, 45.6 fb^{-1} (+23D, not 22G)";
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
  l->SetLineStyle(kDashDotted);
  double etaec = 1.305;
  l->DrawLine(-etaec,-TMath::Pi(),-etaec,+TMath::Pi());
  l->DrawLine(+etaec,-TMath::Pi(),+etaec,+TMath::Pi());
  l->SetLineStyle(kDotted);
  double etazr = 0.0;
  l->DrawLine(-etazr,-TMath::Pi(),-etazr,+TMath::Pi());
  l->DrawLine(+etazr,-TMath::Pi(),+etazr,+TMath::Pi());

  rezero(h2_22cd);
  h2_22cd->GetZaxis()->SetRangeUser(-10,10);
  h2_22cd->SetFillStyle(1001);
  h2_22cd->SetLineColor(kYellow+2);
  h2_22cd->SetFillColorAlpha(kYellow+2, 0.35); // 35% transparent
  h2_22cd->Draw("SAMEBOX");

  rezero(h2_22ef);
  h2_22ef->GetZaxis()->SetRangeUser(-10,10);
  h2_22ef->SetFillStyle(1001);
  h2_22ef->SetLineColor(kBlue);
  h2_22ef->SetFillColorAlpha(kBlue, 0.35); // 35% transparent
  h2_22ef->Draw("SAMEBOX");

  rezero(h2_23bc);
  h2_23bc->GetZaxis()->SetRangeUser(-10,10);
  h2_23bc->SetFillStyle(1001);
  h2_23bc->SetLineColor(kGreen+1);
  h2_23bc->SetFillColorAlpha(kGreen+1, 0.35); // 35% transparent
  h2_23bc->Draw("SAMEBOX");

  rezero(h2_23d);
  h2_23d->GetZaxis()->SetRangeUser(-10,10);
  h2_23d->SetFillStyle(1001);
  h2_23d->SetLineColor(kRed);
  h2_23d->SetFillColorAlpha(kRed, 0.35); // 35% transparent
  h2_23d->Draw("SAMEBOX");
  
  TH2D *h2sum = (TH2D*)h2_22cd->Clone("h2_run3");
  h2sum->Add(h2_22ef);
  h2sum->Add(h2_23bc);
  h2sum->Add(h2_23d);
  TH2D *h2all = (TH2D*)h2sum->Clone("h2all");
  rezero(h2sum,30,10); // overlap min. 3
  //rezero(h2sum); // no overlap needed
  rezero(h2all,40,10); // overlap all 4

  h2sum->SetLineColor(kMagenta+2);//kGray+2);
  h2sum->SetLineStyle(kNone);
  h2sum->SetFillStyle(1001);
  h2sum->SetFillColor(kNone);
  h2sum->DrawClone("SAMEBOX");

  h2all->SetLineColor(kBlack);
  h2all->SetLineStyle(kNone);
  h2all->SetFillStyle(1001);
  h2all->SetFillColor(kNone);
  h2all->DrawClone("SAMEBOX");
  
  //h2em->DrawClone("SAMEBOX");
  //h2em->SetFillColorAlpha(kAzure-9, 0.35); // 35% transparent
  //h2em->Draw("SAMEBOX");

  TLegend *leg = tdrLeg(0.14,0.90-6*0.05,0.34,0.90);
  //leg->AddEntry(h2em,"Cold (GH)","F");
  leg->AddEntry(h2_22cd,"22CD","F");
  leg->AddEntry(h2_22ef,"22EF","F");
  leg->AddEntry(h2_23bc,"23BC","F");
  leg->AddEntry(h2_23d, "23D","F");
  //leg->AddEntry(h2em,"EM mask","F");
  leg->AddEntry(h2sum,"(min. 3)","F");
  leg->AddEntry(h2all,"(all 4)","F");

  gPad->Paint();

  c1->SaveAs("pdf/doJetVetoV2/drawJetVetoV2_Run3.pdf");

  //TFile *fout = new TFile("rootfiles/hotjets-UL16.root","RECREATE");
  //h2sum->Write();
  //fout->Close();
}
