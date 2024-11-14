// Purpose: Draw Z+jet, gamma+jet and inclusive jet xsec and MPF vs time
//          Use breakup made by rootfiles/brilcalc/clusterRuns.C
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"

#include <iostream>

#include "../tdrstyle_mod22.C"

void rebinProfileCustom(TProfile* p, TProfile* p_new);

void drawTimeStability() {

  setTDRStyle();
  TDirectory *curdir = gDirectory;

  TFile *flum = new TFile("rootfiles/brilcalc/clusterRuns.root","READ");
  assert(flum && !flum->IsZombie());

  TH1D *hlum = (TH1D*)flum->Get("hlum"); assert(hlum); // per run
  TH1D *hlum2 = (TH1D*)flum->Get("hlum2"); assert(hlum2); // per range
  TH1D *hcumlum = (TH1D*)flum->Get("hcumlum"); assert(hcumlum); // per run
  TH1D *hcumlum2 = (TH1D*)flum->Get("hcumlum2"); assert(hcumlum2); // per range
  TH1D *hbins = (TH1D*)flum->Get("hcumlumbins2"); assert(hbins); // per range

  string vf[] =
    {"2022CDE_v32","2022FG_v32",
     "2023Cv123_w8","2023Cv4_w8","2023D_w8",
     "2024BCD_w39","2024E_w39","2024F_w39","2024G_w39","2024H_w40","2024I_w40"};
  const int nf = sizeof(vf)/sizeof(vf[0]);
    
  string vh[] =
    {//"pr50n","pr110n","pr230n",
      //"pr50m","pr110m"};//,"pr230m"};
      "pr110m"};
  const int nh = sizeof(vh)/sizeof(vh[0]);

  map<string, int> mcolor;
  mcolor["pr50n"] = kRed;
  mcolor["pr110n"] = kBlue;
  mcolor["pr230n"] = kGreen+2;
  mcolor["pr50m"] = kRed;
  mcolor["pr110m"] = kBlue;
  mcolor["pr230m"] = kGreen+2;
  
  // Loop over files to retrieve stuff
  for (int ih = 0; ih != nh; ++ih) {

    const char *ch = vh[ih].c_str();
    cout << "Analyuzing " << ch << endl << flush;

    TH1D *hsum = (TH1D*)hlum2->Clone(Form("hsum_%s",ch)); hsum->Reset();
    //double vx[hsum->GetNbinsX()+1];
    //for (int i = 1; i != hsum->GetNbinsX()+2; ++i) vx[i-1] = hsum->GetBinLowEdge(i);
    //TProfile *psum = new TProfile(Form("psum_%s",ch),"",hsum->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());
    double vx[hbins->GetNbinsX()+1];
    for (int i = 1; i != hbins->GetNbinsX()+2; ++i) vx[i-1] = hbins->GetBinLowEdge(i);
    TProfile *psum = new TProfile(Form("psum_%s",ch),"",hbins->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());
    TProfile *psumjes = new TProfile(Form("psumjes_%s",ch),"",hbins->GetNbinsX(),&vx[0]);//hsum->GetXaxis()->GetXbins()->GetArray());

    // Keep track of JES (1=L2L3Res) for gamma+jet MPF
    //TH1D *hjes = new TH1D(Form("hjes_%s",ch),";JES;Cum.Lum(/fb)",hbins->GetNbinsX(),&vx[0]);
    
    for (int iff = 0; iff != nf; ++iff) {

      const char *cf = vf[iff].c_str();
      TString tf(cf);
      TFile *f(0);
      if (tf.Contains("2024")) f = new TFile(Form("rootfiles/Prompt2024/GamHistosFill_data_%s.root",cf),"READ");
      if (tf.Contains("2023")) f = new TFile(Form("rootfiles/Summer23_L2L3Res/GamHistosFill_data_%s.root",cf),"READ");
      if (tf.Contains("2022")) f = new TFile(Form("../gamjet/rootfiles/GamHistosFill_data_%s.root",cf),"READ");
      if (!f || f->IsZombie()) cout << "Missing " << cf << endl << flush;
      assert(f && !f->IsZombie());
      curdir->cd();

      TH1D *h(0);
      TProfile *p(0);
      TObject *o = f->Get(Form("runs/%s",ch));
      // Detect object type automatically
      if (o) {
	if (o->InheritsFrom("TProfile"))  p = (TProfile*)o;
	else if (o->InheritsFrom("TH1D")) h = (TH1D*)o;
      }

      // Determine the average JES that was applied to each run
      TProfile2D *p2res = (TProfile2D*)f->Get("Gamjet2/p2res");
      if (!p2res) cout << "Missing Gamjet2/p2res in " << cf << endl << flush;
      double jes(1.);
      if (p2res) {
	// Slice out [110,230]x[-1.3,1.3] region to get mean JES here
	int ipt1 = p2res->GetYaxis()->FindBin(110.);
	int ipt2 = p2res->GetYaxis()->FindBin(230.);
	TProfile *px = p2res->ProfileX("px",ipt1,ipt2);
	int ieta1 = p2res->GetXaxis()->FindBin(-1.3);
	int ieta2 = p2res->GetXaxis()->FindBin(+1.3);
	px->GetXaxis()->SetRangeUser(ieta1,ieta2);
	jes = px->GetMean(2);
	delete px;
      }
      
      if (!h && !p) cout << "Missing " << ch << " in " << cf << endl << flush;
      if (h) {
	for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	  int j = hsum->FindBin(h->GetBinCenter(i));
	  hsum->SetBinContent(j, hsum->GetBinContent(j) + h->GetBinContent(i));
	  hsum->SetBinError(j, sqrt(pow(hsum->GetBinError(j),2) + pow(h->GetBinError(i),2)));
	} // for i
      } // if h
	       
      if (p) {

	cout << "Fill p for file " << cf << endl << flush;
	//rebinProfileCustom(psum, p);

	// Scale out L2L3Res
	TProfile *pjes = (TProfile*)p->Clone(Form("pjes_%s",ch));
	pjes->Scale(jes);
	
	for (int i = 1; i != p->GetNbinsX()+1; ++i) {
	  int run = p->GetBinCenter(i);
	  double cumlum = hcumlum2->GetBinContent(hcumlum2->FindBin(run));
	  //double y = p->GetBinContent(i);
	  //double w = p->GetBinEntries(i);
	  //psum->Fill(cumlum, y, w);

	  // Improved calculation to keep track of RMS
	  int j = psum->FindBin(cumlum);

	  // Add content and entries from the original bin to the new one
	  (*psum)[j] = (*p)[i] + (*psum)[j]; // sumwy
	  (*psum->GetSumw2())[j] = (*p->GetSumw2())[i] + (*psum->GetSumw2())[j]; // sumwy2
	  psum->SetBinEntries(j, p->GetBinEntries(i) + psum->GetBinEntries(j)); // sumw

	  (*psumjes)[j] = (*pjes)[i] + (*psumjes)[j]; // sumwy
	  (*psumjes->GetSumw2())[j] = (*pjes->GetSumw2())[i] + (*psumjes->GetSumw2())[j]; // sumwy2
	  psumjes->SetBinEntries(j, pjes->GetBinEntries(i) + psumjes->GetBinEntries(j)); // sumw
      
	  // Copy (if needed) bin sum of weight square
	  if (p->GetBinSumw2()->fN > i) {
	    psum->Sumw2();
	    (*psum->GetBinSumw2())[j] = (*p->GetBinSumw2())[i] + (*psum->GetBinSumw2())[j]; // sum2
	  }
	  
	  if (pjes->GetBinSumw2()->fN > i) {
	    psumjes->Sumw2();
	    (*psumjes->GetBinSumw2())[j] = (*pjes->GetBinSumw2())[i] + (*psumjes->GetBinSumw2())[j]; // sum2
	  }
      
	  // Accumulate overall profile entries
	  psum->SetEntries(psum->GetEntries() + p->GetEntries());
	  psumjes->SetEntries(psumjes->GetEntries() + pjes->GetEntries());
	} // for i
      } // if p
    } // for iff
    hsum->Divide(hlum2);

    double xmin = hbins->GetXaxis()->GetXmin();
    double xmax = hbins->GetXaxis()->GetXmax();
    double ymin = -2.5;
    double ymax = +5.5;//+2.5;
    //TH1D *h1 = tdrHist(Form("h1_%s",ch),"JES",0.975,1.025,"Cumulative luminosity (fb^{-1})",0,hbins->GetXaxis()->GetXmax());
    TH1D *h1 = tdrHist(Form("h1_%s",ch),"JES-1 (%)",ymin,ymax,"Cumulative luminosity (fb^{-1})",xmin,xmax);
    lumi_136TeV = "Run3, 2022-24";
    extraText = "Private";
    TCanvas *c1 = tdrCanvas(Form("c1_%s",ch),h1,8,11,kRectangular);

    TF1 *f1 = new TF1("f1","[0]",0,1e6);
    if (hsum->Integral()!=0) {
      hsum->Fit(f1,"QRN");
      hsum->Scale(1./f1->GetParameter(0));
      hsum->SetLineColor(mcolor[ch]);
      //hsum->Draw(ih==0 ? "HE" : "HESAME");
      tdrDraw(hsum,"PE",kFullCircle,kBlack,kSolid,-1,kNone,0,0.6);
    }
    if (psum->Integral()!=0) {

      // Normalize average to unity to focus on time dependence
      psum->Fit(f1,"QRN");
      TH1D *h = psum->ProjectionX(Form("h_%s",ch));
      h->Scale(1./f1->GetParameter(0));

      psumjes->Fit(f1,"QRN");
      TH1D *hjes = psumjes->ProjectionX(Form("hjes_%s",ch));
      hjes->Scale(1./f1->GetParameter(0));

      // Turn JES into JES-1 (%)
      for (int i = 1; i != h->GetNbinsX()+1; ++i) {
	h->SetBinContent(i, (h->GetBinContent(i)-1)*100.);
	h->SetBinError(i, h->GetBinError(i)*100.);

	hjes->SetBinContent(i, (hjes->GetBinContent(i)-1)*100.);
	hjes->SetBinError(i, hjes->GetBinError(i)*100.);
      }

      // Horizontal lines at 0 and +/-1% for reference
      TLine *l = new TLine();
      l->SetLineStyle(kDashed);
      l->SetLineColor(kGray+1);
      l->DrawLine(xmin, 0, xmax, 0);
      l->SetLineStyle(kDotted);
      l->DrawLine(xmin, +1, xmax, +1);
      l->DrawLine(xmin, -1, xmax, -1);

      // Lines at era breaks
      // Define era boundaries as pairs of (start_run, end_run)
      std::vector<std::pair<int, int>> era_boundaries = {
        // 2022 (CDE) (FG)
        {355100, 355793}, /*{355794, 357486}, {357487, 357733},*/ {357734, 358219}, /*{358220, 359021}, {359022, 360331},*/
        //{360332, 362180}, {362181, 362349}, {362350, 362760},
        // 2023
        {366442, 367079}, //{367080, 367515}, {367516, 367620}, {367621, 367763}, {367765, 369802}, {369803, 370602}, {370603, 370790},
        // 2024
        {378971, 379411}, /*{379412, 380252}, {380253, 380947},*/ {380948, 381383}, {381384, 381943}, {381944, 383779},
        {383780, 385813}, {385814, 386408}, {386409, 386951}
      };
      // List of eras as (start_run, era name)
      std::vector<std::pair<int, string>> eras = {
        // 2022
        //{355100, "22"},/* {355794, "22C"}, {357487, "22D"}, {357734, "22Dv2"}, {358220, "22Dv3"},*/ {359022, "E"},
        //{360332, "F"}, /*{362181, "22HI"}, {362350, "G"},*/
	{355374, "22"}, {359569, "E"}, {360390, "F"}, /*{362181, "22HI"},*/ {362437, "G"}, // Updated actual first run
        // 2023
        //{366442, "23"}, /*{367080, "23C"}, {367516, "23Cv2"}, {367621, "23Cv3"},*/ {367765, "Cv4"}, {369803, "D"}, /*{370603, "23Dv2"},*/
	{366727, "23"}, {367770, "Cv4"}, {369927, "D"}, // Updated actual first run
        // 2024
        //{378971, "24"}, /*{379412, "24C"}, {380253, "24D"},*/ {380948, "E"}, /*{381384, "24Ev2"},*/ {381944, "F"},
	{378985, "24"}, {380963, "E"}, {382229, "F"}, // Updated actual first run
        //{383780, "G"}, {385814, "H"}, {386409, "I"}
	{383811, "G"}, {385836, "H"}, {386478, "I"} // Updated actual first run
      };

      // Additional HCAL breaks
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalRespCorrsTagsRun3
      // https://twiki.cern.ch/twiki/bin/view/CMS/HcalRespCorrsTags2011
      //eras.push_back(pair<int,string>(383195,"24_v2.0")); // HLT
      //eras.push_back(pair<int,string>(383219,"24_v2.1")); // HLT
      //eras.push_back(pair<int,string>(386401,"24_v3.0")); // HLT, eraH

      // => https://cms-talk.web.cern.ch/t/fast-track-validation-hlt-prompt-hcal-respcorrs-condition-update-from-hb-time-adjustment/25302/5
      //eras.push_back(pair<int,string>(368775,"V1.0")); // 368765, 368822 =>
      eras.push_back(pair<int,string>(368822,"V1.0")); // (actual run, 23D-)
      
      //eras.push_back(pair<int,string>(367765,"V1.0")); // 23_v1.0 => 23Cv4
      //eras.push_back(pair<int,string>(380637,"V1.0")); // Something in mid-24C?
      //eras.push_back(pair<int,string>(380852,"V1.0")); // 24E
      //eras.push_back(pair<int,string>(382287,"V2.0")); // 24_v2.0 =>
      eras.push_back(pair<int,string>(382298,"V2.0")); // 24_v2.0 (actual first run)
      //eras.push_back(pair<int,string>(383219,"V2.1")); // 24_v2.1 =>
      eras.push_back(pair<int,string>(383247,"V2.1")); // 24_v2.1 (actual first run)
      //eras.push_back(pair<int,string>(386401,"V3.0")); // 24_v3.0 => 24I =>
      //eras.push_back(pair<int,string>(386478,"V3.0")); // 24_v3.0 / 24I (actual first run)
      //eras.push_back(pair<int,string>(383811,"VG")); // 24G tes
      
      // Additional ECAL intercalibration updates
      // https://cms-talk.web.cern.ch/c/ppd/alca/108
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-368782-w24/25407
      //era_boundaries.push_back(pair<int,int>(368919,368919)); // 23D-, 368823, 369927 =>
      //eras.push_back(pair<int,string>(369927,"EP")); // =? 23D (actual run)
      // => https://cms-talk.web.cern.ch/t/gt-online-hlt-express-prompt-update-of-ecal-pedestals-conditions-run-384719-w34/46354/10
      //eras.push_back(pair<int,string>(384719,"IC")); // mid-24G! =>
      eras.push_back(pair<int,string>(384933,"IC")); // mid-24G! (actual first run) 384644, 384933
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386025/57306
      //eras.push_back(pair<int,string>(386025,"IC")); // also actual run? or deployed later?
      // => https://cms-talk.web.cern.ch/t/l1-pre-announcement-of-ecal-intercalibration-update-at-l1-run-386945/61017
      // eras.push_back(pair<int,string>(386945,"IC")); // post 24I

      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379660/39889/5
      //eras.push_back(pair<int,string>(380115,"IC")); // mid-24C (too late)
      // => https://cms-talk.web.cern.ch/t/full-track-validation-hlt-prompt-ecalintercalibconstants-conditions-from-runs-378981-379616/39588/4
      //eras.push_back(pair<int,string>(379956,"IC")); // mid-24C 379866, 379984 =>
      eras.push_back(pair<int,string>(379984,"IC")); // (actual run)
      
      //eras.push_back(pair<int,string>(368824,"XX")); // 23D- mystery run
      //eras.push_back(pair<int,string>(368765,"XX")); // 23D- mystery run
      //eras.push_back(pair<int,string>(380030,"XX")); // mid-24C mystery runx
	
      TLatex *tex = new TLatex();
      tex->SetTextSize(0.045);
      int ks(0), kh(0);
      for (int i = 0; i != eras.size(); ++i) {

	int run = eras[i].first;
	string s = eras[i].second;
	const char *cn = s.c_str();
	TString t = TString(cn);

	int j = hcumlum2->FindBin(run);
	//if (hcumlum2->GetBinLowEdge(j)<run) ++j;
	if (hcumlum2->GetBinLowEdge(j)!=run) cout << "run="<<run<<", bin edge="<<hcumlum2->GetBinLowEdge(j)<<endl<<flush;
	double cumlum = hcumlum2->GetBinContent(j);
	//int k = hcumlum->FindBin(run)+1;
	//if (hcumlum->GetBinLowEdge(k)<run) ++k;
	//double cumlum1 = hcumlum->GetBinContent(k);

	l->SetLineColor(kGray);
	if (t.Contains("V")) l->SetLineColor(kRed);
	if (t.Contains("IC")) l->SetLineColor(kBlue);
	if (t.Contains("EP")) l->SetLineColor(kBlue);
	if (t.Contains("XX")) l->SetLineColor(kMagenta+1);
	l->SetLineStyle(kSolid);
	//if (t.Contains("V")) l->DrawLine(cumlum1,ymin,cumlum1,ymax);
	//else
	l->DrawLine(cumlum,ymin,cumlum,ymax);
	if (s=="22"||s=="23"||s=="24") ks = 0;
	//tex->DrawLatex(cumlum+1,-1.2-0.2*(ks++),cn);
	//if (t.Contains("V")) tex->DrawLatex(cumlum1+1,-1.6-0.4*(kh++),cn);
	if (t.Contains("V") || t.Contains("IC") || t.Contains("XX") || t.Contains("EP")) tex->DrawLatex(cumlum+1,-1.6-0.4*(kh++%3),cn);
	else tex->DrawLatex(cumlum+1,+3.5-0.4*(ks++),cn);
      }
            
      tdrDraw(hjes,"PE",kOpenSquare,kRed,kSolid,-1,kNone,0,0.6);
      tdrDraw(h,"PE",kFullCircle,kGreen+2,kSolid,-1,kNone,0,0.6);

      TLegend *leg = tdrLeg(0.65,0.85-0.05*3,0.90,0.85);
      leg->SetHeader("#gamma+jet MPF 110EB");
      leg->AddEntry(hjes,"Before JEC");
      leg->AddEntry(h,"After JEC");
    }

    c1->RedrawAxis();
    c1->SaveAs(Form("pdf/drawTimeStability/drawTimeStability_%s.pdf",ch));
  } // for ih
  
  
} // drawTimeStability



void rebinProfileCustom(TProfile* p, TProfile* p_new) {
  // Iterate through the bins of the original TProfile
  for (int ix = 1; ix <= p->GetNbinsX(); ix++) {
                
    // Original bin index
    int bin_orig = ix;
      
    // Get bin centers from original profile
    double x_center = p->GetXaxis()->GetBinCenter(ix);

    // Find corresponding bin in the new profile
    int new_ix = p_new->GetXaxis()->FindBin(x_center);
      
    // Get the new bin index
    int bin_new = new_ix;
      
    // profile keeps track of sumy, sumwy, sumwy2, sumw2
    // sumw=fArray, sumwy=fBinEntries.fArray, 
    // sumwy2 = fBinSumw2.fArray, sumw2 = fSum2.fArray
    // GetBinContent = sumwy/sumw
    // https://root-forum.cern.ch/t/copy-entries-of-tprofile/11828
      
    // Add content and entries from the original bin to the new one
    (*p_new)[bin_new] = (*p)[bin_orig] + (*p_new)[bin_new]; // sumwy
    (*p_new->GetSumw2())[bin_new] = (*p->GetSumw2())[bin_orig] +
      (*p_new->GetSumw2())[bin_new]; // sumwy2
    p_new->SetBinEntries(bin_new, p->GetBinEntries(bin_orig) + 
                          p_new->GetBinEntries(bin_new)); // sumw
      
    // Copy (if needed) bin sum of weight square
    if (p->GetBinSumw2()->fN > bin_orig) {
      p_new->Sumw2();
      (*p_new->GetBinSumw2())[bin_new] = (*p->GetBinSumw2())[bin_orig] +
        (*p_new->GetBinSumw2())[bin_new]; // sum2
    }
      
    // Accumulate overall profile entries
    p_new->SetEntries(p_new->GetEntries() + p->GetEntries());
  } // for ix
} // rebinProfileCustom

