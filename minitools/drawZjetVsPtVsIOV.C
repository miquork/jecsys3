// Purpose: draw stability of Z+jet results vs pT over time
//          compare data and MC for response etc 
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TLine.h"

#include "../tools.C"
#include "../tdrstyle_mod22.C"

//bool addMPFu2n = true;
// v61 (no L2L3Res -> v63/v64 L2L3Res->v66(L2L3Res V3))
string id = "v66";
//bool drawFullIOVList = false;//true;

// Forward declaration of call
void drawZjetVsPtVsIOVs(string so, string var, string name,
			double y1, double y2, double z1, double z2);

// Multiple calls to draw function
void drawZjetVsPtVsIOV() {

  drawZjetVsPtVsIOVs("%s/eta_00_13/rmpf_zmmjet_a100","MPF","MPF",0.90,1.15,0.89,1.09);
  drawZjetVsPtVsIOVs("%s/eta_00_13/rmpfjet1_zmmjet_a100","DB","DB",0.75,1.20,0.89,1.09);

} // drawZjetVsPtVsIOV

void drawZjetVsPtVsIOVs(string so, string var, string name,
			double y1, double y2, double z1, double z2) {
  setTDRStyle();
  TDirectory *curdir = gDirectory;
  
  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD","Run2"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8","Run2P8"};
  //string iovs[] = {"2016BCDEF","2016FGH","2017BCDEF","2018ABCD"};
  //string mcs[] = {"2016APVP8","2016P8","2017P8","2018P8"};
  //string iovs_long[] = {
  //"2022C","2022D","2022E","2022F","2022G",
  //"2023Cv123","2023Cv4","2023D"
  //};
  //string iovs_short[] = {
  //"2018ABCD","Run3",
  //"2022CD","2022E","2022FG", //hadd
  //"2023Cv123","2023Cv4D", //hadd
  //};

  string iovs[] = {
    "Run3",
    "2022CD","2022E","2022FG",
    //"2022CD","2022E","2022F","2022G",
    //"2023C123","2023C4D"
    "2023C123","2023C4","2023D",
  };
  const int niov = sizeof(iovs)/sizeof(iovs[0]);
  
  //string mcs_long[] = {
  //"2022P8","2022P8","2022EEP8","2022EEP8","2022EEP8", //orig
  //"2022P8","2022P8","2022P8"
  //};
  //string mcs_short[] = {
  //"2018P8","Run3P8",
  //"2022P8","2022EEP8","2022EEP8", //hadd
  //"2022P8","2022P8"
  //};
  //const int niov_long = sizeof(iovs_long)/sizeof(iovs_long[0]);
  //const int nmc_long = sizeof(mcs_long)/sizeof(mcs_long[0]);
  //const int niov_short = sizeof(iovs_short)/sizeof(iovs_short[0]);
  //const int nmc_short = sizeof(mcs_short)/sizeof(mcs_short[0]);

  //string *iovs = (drawFullIOVList ? &iovs_long[0] : &iovs_short[0]);
  //string *mcs = (drawFullIOVList ? &mcs_long[0] : &mcs_short[0]);
  //const int niov = (drawFullIOVList ? niov_long : niov_short);
  //const int nmc = (drawFullIOVList ? nmc_long : nmc_short);

  //assert(niov==nmc);
  
  map<string,int> mcolor;
  mcolor["2016BCDEF"] = kBlue;
  mcolor["2016FGH"] = kCyan+2;
  mcolor["2017BCDEF"] = kGreen+2;
  mcolor["2018ABCD"] = kGray+2;//kRed;
  mcolor["Run2"] = kBlack;
  //
  mcolor["2022C"] = kGreen+2;//kBlue;
  mcolor["2022D"] = kGreen+3;//kCyan+1;
  mcolor["2022E"] = kCyan+2;
  mcolor["2022CD"] = kGreen+2;//kBlue;
  mcolor["2022CDE"] = kGreen+2;//kBlue;
  mcolor["2022F"] = kRed;
  mcolor["2022G"] = kRed+1;//kOrange+2;
  mcolor["2022FG"] = kRed;
  mcolor["2023C123"] = kOrange+1;//kYellow+2;
  mcolor["2023C4"] = kBlue;//kGreen+2;
  mcolor["2023D"] = kMagenta+2;
  mcolor["2023C4D"] = kBlue;//kGreen+2;
  mcolor["Run3"] = kBlack;

  map<string,int> mmarker;
  mmarker["2022C"] = kFullSquare;
  mmarker["2022D"] = kOpenSquare;
  mmarker["2022E"] = kOpenSquare;
  mmarker["2022CD"] = kFullSquare;
  mmarker["2022CDE"] = kFullSquare;
  mmarker["2022F"] = kFullTriangleUp;
  mmarker["2022G"] = kOpenTriangleUp;
  mmarker["2022FG"] = kFullTriangleUp;
  mmarker["2023C123"] = kFullDiamond;//kFullCircle;
  mmarker["2023BC123"] = kFullDiamond;//kFullCircle;
  mmarker["2023C4"] = kFullTriangleDown;//kFullDiamond;
  mmarker["2023D"] = kOpenTriangleDown;//kOpenDiamond;
  mmarker["2023C4D"] = kFullTriangleDown;//kFullDiamond;
  mmarker["Run3"] = kFullSquare;
  
  const char *cvar = var.c_str();
  const char *cname = name.c_str();

  TH1D *h = tdrHist("h",cvar,y1,y2);
  TH1D *h2 = tdrHist("h2","Data/MC",z1,z2);
  //lumi_13TeV = "Run2 v16";
  //lumi_13TeV = "Run3 v21";
  lumi_13TeV = "Run2"; // 4=13 TeV
  if (id!="") lumi_136TeV = Form("Run3 %s",id.c_str()); // 8=13.6 TeV
  TCanvas *c1 = tdrDiCanvas(Form("c1_%s",cname),h,h2,8,11);


  c1->cd(1);
  TLine *l = new TLine();
  l->SetLineStyle(kDashed);
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
  gPad->SetLogx();

  //TLegend *leg = tdrLeg(0.60,0.90-niov*0.06,0.80,0.90);
  TLegend *leg = tdrLeg(0.60,0.90-niov*0.045,0.80,0.90);

  c1->cd(2);
  gPad->SetLogx();
  l->DrawLine(h->GetXaxis()->GetXmin(),1,h->GetXaxis()->GetXmax(),1);
  l->DrawLine(h->GetXaxis()->GetXmin(),0,h->GetXaxis()->GetXmax(),0);
  
  for (int i = 0; i != niov; ++i) {
    
    string iov = iovs[i];
    const char *ciov = iov.c_str();
    //const char *cmc = mcs[i].c_str();
    const char *cid = id.c_str();
    
    TFile *f(0);
    if (iovs[i]=="2018ABCD") {
      assert(false);
    }
    else {
      f = new TFile(Form("rootfiles/jme_bplusZ_%s_Zmm_sync_%s.root",ciov,cid));
    }
    assert(f && !f->IsZombie());
    
    curdir->cd();
    
    TObject *od = f->Get(Form(so.c_str(),"data")); assert(od);
    TObject *om = f->Get(Form(so.c_str(),"mc")); assert(om);

    /*
    TH1D *hd(0), *hm(0);
    if (od->InheritsFrom("TProfile")) {
      hd = ((TProfile*)od)->ProjectionX(Form("hd_%s_%s",cvar,ciov));
      hm = ((TProfile*)om)->ProjectionX(Form("hm_%s_%s",cvar,ciov));

      if (name=="MPFn" && addMPFu2n) {
	const char *co = "resp_MpfRuchs_%s_a100_eta00_13";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,1.2);
	hm->Add(pm,1.2);
      }
      if ((name=="MPF" || name=="MPF1") && addG1toMPF) {
	const char *co = "control/pgain1vspt";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,0.02);
      }
      if ((name=="MPF" || name=="MPF1") && addG12toMPF) {
	const char *co = "control/pgain12vspt";
	TProfile *pd = (TProfile*)fd->Get(Form(co,"DATA")); assert(pd);
	TProfile *pm = (TProfile*)fm->Get(Form(co,"MC")); assert(pm);
	hd->Add(pd,0.0);
      }
    }
    else {
      hd = (TH1D*)od;
      hm = (TH1D*)om;
    }
    assert(hd);
    assert(hm);
    
    TH1D *hr = (TH1D*)hd->Clone(Form("hr_%s_%s",cvar,ciov));
    if (name=="MPFn" || name=="MPFu" || name=="Leakage" || name=="Leakage0") {
      hr->Add(hm,-1);
      h2->SetYTitle("Data-MC");
    }
    else
      hr->Divide(hm);
    */

    TGraphErrors *gd(0), *gm(0), *gr(0);
    if (od->InheritsFrom("TGraphErrors")) {
      gd = (TGraphErrors*)od;
      gm = (TGraphErrors*)om;
      gr = tools::ratioGraphs(gd, gm);
    }
    else {
      assert(false);
    }
    
    c1->cd(1);
    gPad->SetLogx();
    //tdrDraw(hm,"H",kNone,(mcolor[iov] ? mcolor[iov] : kBlack),kSolid,-1,kNone);
    //tdrDraw(hd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
    //	    (mcolor[iov] ? mcolor[iov] : kBlack));
    
    //leg->AddEntry(hd,ciov,"PLE");
    
    tdrDraw(gm,"L0z",kNone,(mcolor[iov] ? mcolor[iov] : kBlack),kSolid,-1,kNone);
    tdrDraw(gd,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
	    (mcolor[iov] ? mcolor[iov] : kBlack));
    
    leg->AddEntry(gd,ciov,"PLE");

    c1->cd(2);
    gPad->SetLogx();
    //tdrDraw(hr,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
    //	    (mcolor[iov] ? mcolor[iov] : kBlack));
    tdrDraw(gr,"Pz",(mmarker[iov] ? mmarker[iov] : kFullCircle),
    	    (mcolor[iov] ? mcolor[iov] : kBlack));
  } // for iov

  if (id!="")
      c1->SaveAs(Form("pdf/drawZjetVsPtVsIOVs/drawZjetVsPtVsIOVs_%s_%s.pdf",
		      name.c_str(),id.c_str()));
  else
    c1->SaveAs(Form("pdf/drawZjetVsPtVsIOVs/drawZjetVsPtVsIOVs_%s.pdf",name.c_str()));
} // void drawPhotonJetVsPtVsIOVs
