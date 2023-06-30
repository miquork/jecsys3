// Purpose: reformat DATA_ReReco-Run2022C-JetHT_PFCutVariation.root and
//          DATA_PromptReco-Run2022G-JetMET_PFCutVariation.root
//          to use with fullSimShapes.C
//          New: merge as patch
TH1D *extendHist(TH1D* h) {
  vector<double> x(h->GetNbinsX()+1+2);
  for (int i = 1; i != h->GetNbinsX()+2; ++i) {
    x[i-1] = h->GetBinLowEdge(i);
  }
  x[h->GetNbinsX()+1+0] = 3000;
  x[h->GetNbinsX()+1+1] = 3500;
  string name = h->GetName(); h->SetName(Form("%sX",name.c_str()));
  TH1D *hnew = new TH1D(name.c_str(),h->GetTitle(),x.size()-1,&x[0]);
  for (int i = 1; i != h->GetNbinsX()+1; ++i) {
    hnew->SetBinContent(i, h->GetBinContent(i));
    hnew->SetBinError(i, h->GetBinError(i));
  }
  return hnew;
}

void reformatHcalNoise() {

  TFile *f = new TFile("rootfiles/DATA_PromptReco-Run2022G-JetMET_PFCutVariation.root","READ");
  TFile *f2 = new TFile("rootfiles/DATA_ReReco-Run2022C-JetHT_PFCutVariation.root","READ");
  assert(f && !f->IsZombie());
  assert(f2 && !f2->IsZombie());
  TFile *fout = new TFile("rootfiles/DATA_Run2022GC_HBnoise.root","RECREATE");
  //DATA_Run2022C_HBnoise.root","RECREATE");
  assert(fout && !fout->IsZombie());

  TProfile *pr = (TProfile*)f->Get("Rjet_newBinning");
  assert(pr);
  TProfile *pr2 = (TProfile*)f2->Get("Rjet_newBinning");
  assert(pr2);

  TProfile *pchf = (TProfile*)f2->Get("chf");
  assert(pchf);
  TH1D *hchf = pchf->ProjectionX("chf_HBnoise");
  hchf = extendHist(hchf);

  TProfile *pnhf = (TProfile*)f2->Get("nhf");
  assert(pnhf);
  TH1D *hnhf = pnhf->ProjectionX("nhf_HBnoise");
  hnhf = extendHist(hnhf);

  TProfile *pnef = (TProfile*)f2->Get("gammaf");
  assert(pnef);
  TH1D *hnef = pnef->ProjectionX("gammaf_HBnoise");
  hnef = extendHist(hnef);

  TH1D *hr = pr->ProjectionX("Rjet_HBnoise");
  TH1D *hr2 = pr2->ProjectionX("Rjet_HBnoise2");
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {
    double a = hr->GetBinContent(i);
    double ae = hr->GetBinError(i);
    // a=(c-b)/(c+b) => ac+ab=c-b => (1+a)b=(1-a)c => c/b=(1+a)/(1-a)~1+2a
    double r = (1+a)/(1-a);
    hr->SetBinContent(i, r);
    hr->SetBinError(i, 2*ae);

    double a2 = hr2->GetBinContent(i);
    double ae2 = hr2->GetBinError(i);
    double r2 = (1+a2)/(1-a2);
    hr2->SetBinContent(i, r2);
    hr2->SetBinError(i, 2*ae2);
  }
  hr = extendHist(hr);
  hr2 = extendHist(hr2);

  // Fix CHF in 2022G
  /*
  for (int i = 1; i != hchf->GetNbinsX()+1; ++i) {
    double nef = hnef->GetBinContent(i);
    double enef = hnef->GetBinError(i);
    double nhf = hnhf->GetBinContent(i);
    double enhf = hnhf->GetBinError(i);
    double chf = 1 - nef - nhf;
    double echf = sqrt(pow(enef,2)+pow(enhf,2));
    hchf->SetBinContent(i, chf);
    hchf->SetBinError(i, echf);
  }
  */
  // Scale RunG Rjet by x1/3 and match CHF,NHF,NEF from RunC
  for (int i = 1; i != hr->GetNbinsX()+1; ++i) {

    if (hr->GetBinContent(i)==0 && hr->GetBinError(i)==0) continue;

    double dr = (1-hr->GetBinContent(i)) / 3.;
    double r = 1 - dr;
    double er = hr->GetBinError(i) / 3.;
    double dr2 = (1-hr2->GetBinContent(i));

    double k = (dr2!=0 ? dr / dr2 : 1);
    double chf = hchf->GetBinContent(i) * k;
    double echf = hchf->GetBinError(i) * k;
    double nhf = hnhf->GetBinContent(i) * k;
    double enhf = hnhf->GetBinError(i) * k;
    double nef = hnef->GetBinContent(i) * k;
    double enef = hnef->GetBinError(i) * k;

    hr->SetBinContent(i, r);
    hr->SetBinError(i, er);
    hchf->SetBinContent(i, chf);
    hchf->SetBinError(i, echf);
    hnhf->SetBinContent(i, nhf);
    hnhf->SetBinError(i, enhf);
    hnef->SetBinContent(i, nef);
    hnef->SetBinError(i, enef);
  }

  fout->cd();
  hr->Write("Rjet_HBnoise",TObject::kOverwrite);
  hchf->Write("chf_HBnoise",TObject::kOverwrite);
  hnhf->Write("nhf_HBnoise",TObject::kOverwrite);
  hnef->Write("gammaf_HBnoise",TObject::kOverwrite);
  fout->Close();
}
