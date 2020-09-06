#ifndef __FITX_FIT_H_
#define __FITX_FIT_H_

#include "xjjrootuti.h"
#include "xjjcuti.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <RooGlobalFunc.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include <RooAddPdf.h>
#include <RooConstVar.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooStats/SPlot.h>
#include <RooMsgService.h>
#include <RooStats/ModelConfig.h>
#include <string>
#include <iostream>
#include <iomanip>

#include "fitX.h"

/* ----------------------------------------
// ==> Usage <==
#include "fit.h"

TH1F* h;                         // data binned histogram
TH1F *hmcp_a, *hmcp_b;           // mc mass template histogram
RooDataSet* dsh;                 // data unbinned RooDataSet
RooDataSet *dshmcp_a, *dshmcp_b; // mc mass template RooDataSet

// ====>
std::map<std::string, fitX::fitXresult*> result = fitX::fit(h, 0, hmcp_a, hmcp_b, dsh, dshmcp_a, dshmcp_b, "directory/", false, (bool)saveornot, "name");
float N_psi2s = result["unbinned"]->ysig_a();
float Nerr_psi2s = result["unbinned"]->ysigerr_a();
float N_psi2s = result["unbinned"]->ysig_a();
float Nerr_psi2s = result["unbinned"]->ysigerr_a();
// <====
---------------------------------------- */

namespace fitX
{
  class fitXresult
  {
  public:
    fitXresult(float ysig_a, float ysig_b, float ysigerr_a, float ysigerr_b, 
               float msig_a, float msig_b, float msigerr_a, float msigerr_b) : 
      fysig_a(ysig_a), fysig_b(ysig_b), fysigerr_a(ysigerr_a), fysigerr_b(ysigerr_b),
      fmsig_a(msig_a), fmsig_b(msig_b), fmsigerr_a(msigerr_a), fmsigerr_b(msigerr_b) { ; }
    float ysig_a() { return fysig_a; }
    float ysig_b() { return fysig_b; }
    float ysigerr_a() { return fysigerr_a; }
    float ysigerr_b() { return fysigerr_b; }
    float msig_a() { return fmsig_a; }
    float msig_b() { return fmsig_b; }
    float msigerr_a() { return fmsigerr_a; }
    float msigerr_b() { return fmsigerr_b; }
    TF1* f() { return ff; }
    TF1* fsig_a() { return ffsig_a; }
    TF1* fsig_b() { return ffsig_b; }
    TF1* fbkg() { return ffbkg; }
    float minNll() { return fminNll; }
    TH2F* hcorr() { return fhcorr; }
    RooWorkspace* ww() { return fww; }
    void setf(TF1* f, TF1* fsig_a, TF1* fsig_b, TF1* fbkg) { ff = f; ffsig_a = fsig_a; ffsig_b = fsig_b; ffbkg = fbkg; }
    void setw(RooWorkspace* w, float minNll, TH2F* hcorr) { fww = w; fminNll = minNll; fhcorr = hcorr; }

  private:
    float fysig_a, fysig_b, fysigerr_a, fysigerr_b;
    float fmsig_a, fmsig_b, fmsigerr_a, fmsigerr_b;
    float fminNll;
    TF1 *ff, *ffsig_a, *ffsig_b, *ffbkg;
    RooWorkspace* fww;
    TH2F* fhcorr;
  };
  // ===>
  std::map<std::string, fitX::fitXresult*> fit(TH1F* hh, TH1F* hh_ss, TH1F* hhmc_a, TH1F* hhmc_b, 
                                               RooDataSet* dshh, RooDataSet* dshhmc_a, RooDataSet* dshhmc_b,
                                               std::string outputdir, float fixmean, bool saveplot, std::string name="", std::string title="", std::string option="default", bool silence=true,
                                               std::string cmsleft="#scale[1.25]{#bf{CMS}}", std::string cmsright="1.7 nb^{-1} (2018 PbPb 5.02 TeV)");
  // <===
  void setmasshist(TH1* h, float xoffset=0, float yoffset=0, Color_t pcolor=kBlack);
  void setmasshist(RooPlot* h, float xoffset=0, float yoffset=0, Color_t pcolor=kBlack);

  int NBIN = 38, NBIN_L = 50, NBIN_H = 50;
  float BIN_MIN = 3.62, BIN_MAX = 4.0, BIN_MIN_L = 3.64, BIN_MAX_L = 3.74, BIN_MIN_H = 3.82, BIN_MAX_H = 3.92;
  // int NBIN = 30, NBIN_L = 50, NBIN_H = 50;
  // float BIN_MIN = 3.62, BIN_MAX = 3.92, BIN_MIN_L = 3.64, BIN_MAX_L = 3.74, BIN_MIN_H = 3.82, BIN_MAX_H = 3.92;
  float BIN_WIDTH = (BIN_MAX-BIN_MIN)/NBIN*1.0, BIN_WIDTH_L = (BIN_MAX_L-BIN_MIN_L)/NBIN_L*1.0, BIN_WIDTH_H = (BIN_MAX_H-BIN_MIN_H)/NBIN_H*1.0;
  const int NMVA_BIN = 30;
  const float MVA_BIN_MIN = 0.06, MVA_BIN_MAX = 0.36;
  float MVA_BIN_WIDTH = (MVA_BIN_MAX-MVA_BIN_MIN)/NMVA_BIN*1.0;

  const float PDG_MASS_X = 3.87169, PDG_MASS_X_ERR = 0.00017, FIT_MASS_X = 3.867, FIT_MASS_X_WIN = 0.01;
  const float PDG_MASS_PSI2S = 3.686097, PDG_MASS_PSI2S_ERR = 0.000010, FIT_MASS_PSI2S = 3.686097, FIT_MASS_PSI2S_WIN = 0.01;

  void drawpull(TH1* hmc, TF1* f, Color_t color);

  double getparmin(TF1* f, int ipar) { double parmin, parmax; f->GetParLimits(ipar, parmin, parmax); return parmin; } 
  double getparmax(TF1* f, int ipar) { double parmin, parmax; f->GetParLimits(ipar, parmin, parmax); return parmax; } 

  void labelsmc(std::string label, double mean, double sigma1, double sigma2);
  void labelsdata(std::string label, std::string cmsleft, std::string cmsright);
  void labelsdata_params(double mean_a, double mean_a_err, double yield_a, double yield_a_err,
                         double mean_b, double mean_b_err, double yield_b, double yield_b_err,
                         double chi2prob);
  void zeroparameters(TF1* f, std::vector<int> ipars);
  std::map<std::string, TF1*> resolvef(TF1* f, std::string name);
  void plotonmc(RooPlot* fremptymc, RooDataSet* dshmc, RooAbsPdf* sig_mc, RooAbsPdf* sig_mc1, RooAbsPdf* sig_mc2, RooAbsPdf* sig_mc3, Color_t color, std::string name, std::string option);
  TF1* astf(RooAbsPdf* pdf, TF1* f, std::string name, float width=fitX::BIN_WIDTH);
  TF1* astfsig(RooAbsPdf* pdf, TF1* f, std::string name, float norm, std::vector<int> pars);

  void printfit(TF1* f, RooAbsPdf* pdf);
}


// --->
std::map<std::string, fitX::fitXresult*> fitX::fit(TH1F* hh, TH1F* hh_ss, TH1F* hhmc_a, TH1F* hhmc_b, 
                                                   RooDataSet* dshh, RooDataSet* dshhmc_a, RooDataSet* dshhmc_b, 
                                                   std::string outputdir, float fixmean, bool saveplot, std::string name, std::string title, std::string option, bool silence,
                                                   std::string cmsleft, std::string cmsright)
{
  if(option=="range")
    { NBIN = 30; BIN_MAX = 3.92; }

  std::string uniqstr(xjjc::currenttime());
  if(saveplot) gSystem->Exec(Form("mkdir -p %s", outputdir.c_str()));
  if(silence) RooMsgService::instance().setSilentMode(true);

  std::map<std::string, fitX::fitXresult*> fitresult;
  RooWorkspace* ww = new RooWorkspace(Form("ww%s", name.c_str()));

  /***************************************************
    Preparations
  ***************************************************/

  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooRealVar* massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  RooRealVar* massmc_b = new RooRealVar("Bmass", "massmc_b", fitX::BIN_MIN_H, fitX::BIN_MAX_H);

  RooPlot* frempty = mass->frame(RooFit::Title(""));
  frempty->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  RooPlot* fremptymc_a = massmc_a->frame(RooFit::Title(""));
  fremptymc_a->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH_L*1.e+3));
  RooPlot* fremptymc_b = massmc_b->frame(RooFit::Title(""));
  fremptymc_b->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH_H*1.e+3));

  TH1F* h = new TH1F(*hh); h->SetName(std::string("h_"+uniqstr).c_str());
  TH1F* h_ss = 0; if(hh_ss) { h_ss = new TH1F(*hh_ss); h_ss->SetName(std::string("h_ss_"+uniqstr).c_str()); }
  TH1F* hmc_a = new TH1F(*hhmc_a); hmc_a->SetName(std::string("hmc_a_"+uniqstr).c_str());
  TH1F* hmc_b = new TH1F(*hhmc_b); hmc_b->SetName(std::string("hmc_b_"+uniqstr).c_str());
  RooDataSet* dsh = new RooDataSet(*dshh, std::string("dsh_"+uniqstr).c_str()); 
  RooDataSet* dshmc_a = new RooDataSet(*dshhmc_a, std::string("dshmc_a_"+uniqstr).c_str());
  RooDataSet* dshmc_b = new RooDataSet(*dshhmc_b, std::string("dshmc_b_"+uniqstr).c_str());

  dsh->Print();
  dshmc_a->Print();
  dshmc_b->Print();

  std::vector<TF1*> funs;
  xjjroot::setgstyle(2);
  gStyle->SetPadLeftMargin(gStyle->GetPadLeftMargin()*0.7);

  fitX::setmasshist(h, 0, -0.2);
  h->SetMinimum(0);
  fitX::setmasshist(frempty, 0, -0.2);
  frempty->SetMinimum();
  if(h_ss) fitX::setmasshist(h_ss, 0, -0.2, color_ss);
  fitX::setmasshist(hmc_a, 0, 0.05);
  hmc_a->SetMaximum(hmc_a->GetMaximum()*1.2);
  hmc_a->SetMinimum(0 - hmc_a->GetMaximum()*0.1);
  hmc_a->SetNdivisions(-505);
  fitX::setmasshist(hmc_b, 0, 0.05);
  hmc_b->SetMaximum(hmc_b->GetMaximum()*1.2);
  hmc_b->SetMinimum(0 - hmc_b->GetMaximum()*0.1);
  hmc_b->SetNdivisions(-505);
  fitX::setmasshist(fremptymc_a, 0, 0.05);
  fremptymc_a->SetMaximum(hmc_a->GetMaximum());
  fremptymc_a->SetMinimum(0 - hmc_a->GetMaximum()*0.1);
  fremptymc_a->SetNdivisions(-505);  
  fitX::setmasshist(fremptymc_b, 0, 0.05);
  fremptymc_b->SetMaximum(hmc_b->GetMaximum());
  fremptymc_b->SetMinimum(0 - hmc_b->GetMaximum()*0.1);
  fremptymc_b->SetNdivisions(-505);

  TCanvas* c = new TCanvas("c", "", 700, 600);
  TCanvas* cmc = new TCanvas("cmc", "", 1200, 600);
  cmc->Divide(2, 1);
  TCanvas* cr = new TCanvas("cr", "", 700, 600);
  TCanvas* crnofit = new TCanvas("crnofit", "", 700, 600);
  TCanvas* crmc = new TCanvas("crmc", "", 1200, 600);
  crmc->Divide(2, 1);
  
  TString str_bkg = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x";
  if(xjjc::str_contains(option,"cheb")) { str_bkg = option; }
  // if(xjjc::str_contains(option,"cheb")) { str_bkg = Form("[19]*(%s)", option.c_str()); }
  TString str_sig_a = "[5]*( [9]*TMath::Gaus(x, [6], [7])/(TMath::Sqrt(2*3.14159)* [7]) + (1- [9])*( [16]*TMath::Gaus(x,[6], [8])/(TMath::Sqrt(2*3.14159)* [8]) + (1- [16])*TMath::Gaus(x,[6], [15])/(TMath::Sqrt(2*3.14159)*[15])))";
  TString str_sig_b = "[10]*([14]*TMath::Gaus(x,[11],[12])/(TMath::Sqrt(2*3.14159)*[12]) + (1-[14])*([18]*TMath::Gaus(x,[11],[13])/(TMath::Sqrt(2*3.14159)*[13]) + (1-[18])*TMath::Gaus(x,[11],[17])/(TMath::Sqrt(2*3.14159)*[17])))";

  TF1* f_sig_a = new TF1("f_sig_a", "[6]*([4]*TMath::Gaus(x, [0], [1])/(TMath::Sqrt(2*3.14159)*[1]) + (1-[4])*([5]*TMath::Gaus(x,[0], [2])/(TMath::Sqrt(2*3.14159)*[2]) + (1-[5])*TMath::Gaus(x,[0], [3])/(TMath::Sqrt(2*3.14159)*[3])))", BIN_MIN_L, BIN_MAX_L);
  TF1* f_sig_b = new TF1("f_sig_b", "[6]*([4]*TMath::Gaus(x, [0], [1])/(TMath::Sqrt(2*3.14159)*[1]) + (1-[4])*([5]*TMath::Gaus(x,[0], [2])/(TMath::Sqrt(2*3.14159)*[2]) + (1-[5])*TMath::Gaus(x,[0], [3])/(TMath::Sqrt(2*3.14159)*[3])))", BIN_MIN_H, BIN_MAX_H);
  f_sig_a->SetNpx(1000);
  f_sig_b->SetNpx(1000);
  TF1* f = new TF1("f", str_bkg+"+"+str_sig_a+"+"+str_sig_b, BIN_MIN, BIN_MAX);
  const double parinit[] = {0., 0.     , 0.     , 0.  , 0. , // 0-4
                            1., 3.686  , 0.00357, 0.005, 0.5, // 5-9
                            1., 3.87169, 0.001  , 0.008, 0.5, // 10-14
                            0.005, 1, 0.005, 1, // 15-18 
  };
  f->SetNpx(1000);
  f->SetParameters(parinit);
  f->SetParLimits(9, 0, 1);
  f->SetParLimits(16, 0, 1); if(option != "3gaus") { f->FixParameter(16, 1); }
  f->SetParLimits(14, 0, 1);
  f->SetParLimits(18, 0, 1); if(option != "3gaus") { f->FixParameter(18, 1); }
  f->SetParLimits(7, 1.e-3, 0.01);
  f->SetParLimits(8, 1.e-3, 0.01);
  f->SetParLimits(15, 1.e-3, 0.01);
  f->SetParLimits(12, 1.e-3, 0.01);
  f->SetParLimits(13, 1.e-3, 0.01);
  f->SetParLimits(17, 1.e-3, 0.01);
  f->FixParameter(0, 0);
  f->FixParameter(1, 0);
  f->FixParameter(2, 0);
  f->FixParameter(3, 0);
  f->FixParameter(4, 0);
  // if(xjjc::str_contains(option,"cheb")) { f->FixParameter(19, 0); f->FixParameter(0, 1);}

  std::map<int, RooRealVar*> mcpars;
  // mcpars[5] = new RooRealVar("mcpar5", "", parinit[5], 0, 1.e+5);
  mcpars[6] = new RooRealVar("mcpar6", "", parinit[6], fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
  mcpars[7] = new RooRealVar("mcpar7", "", parinit[7], getparmin(f, 7), getparmax(f, 7));
  mcpars[8] = new RooRealVar("mcpar8", "", parinit[8], getparmin(f, 8), getparmax(f, 8));
  mcpars[15] = new RooRealVar("mcpar15", "", parinit[15], getparmin(f, 15), getparmax(f, 15));
  mcpars[9] = new RooRealVar("mcpar9", "", parinit[9], getparmin(f, 9), getparmax(f, 9));
  mcpars[16] = new RooRealVar("mcpar16", "", parinit[16], getparmin(f, 16), getparmax(f, 16));
  if(option!="3gaus") { mcpars[16]->setVal(1); mcpars[16]->setConstant(); }
  RooGaussian* sig_mc_a1 = new RooGaussian("sig_mc_a1", "", *massmc_a, *(mcpars[6]), *(mcpars[7]));
  RooGaussian* sig_mc_a2 = new RooGaussian("sig_mc_a2", "", *massmc_a, *(mcpars[6]), *(mcpars[8]));
  RooGaussian* sig_mc_a3 = new RooGaussian("sig_mc_a3", "", *massmc_a, *(mcpars[6]), *(mcpars[15]));
  RooAddPdf* sig_mc_a = new RooAddPdf("sig_mc_a", "", RooArgList(*sig_mc_a1, *sig_mc_a2, *sig_mc_a3), RooArgList(*(mcpars[9]), *(mcpars[16])), true);
  // mcpars[10] = new RooRealVar("mcpar10", "", parinit[10], 0, 1.e+5);
  mcpars[11] = new RooRealVar("mcpar11", "", parinit[11], fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
  mcpars[12] = new RooRealVar("mcpar12", "", parinit[12], getparmin(f, 12), getparmax(f, 12));
  mcpars[13] = new RooRealVar("mcpar13", "", parinit[13], getparmin(f, 13), getparmax(f, 13));
  mcpars[17] = new RooRealVar("mcpar17", "", parinit[17], getparmin(f, 17), getparmax(f, 17));
  mcpars[14] = new RooRealVar("mcpar14", "", parinit[14], getparmin(f, 14), getparmax(f, 14));
  mcpars[18] = new RooRealVar("mcpar18", "", parinit[18], getparmin(f, 18), getparmax(f, 18));
  if(option!="3gaus") { mcpars[18]->setVal(1); mcpars[18]->setConstant(); }
  RooGaussian* sig_mc_b1 = new RooGaussian("sig_mc_b1", "", *massmc_b, *(mcpars[11]), *(mcpars[12]));
  RooGaussian* sig_mc_b2 = new RooGaussian("sig_mc_b2", "", *massmc_b, *(mcpars[11]), *(mcpars[13]));
  RooGaussian* sig_mc_b3 = new RooGaussian("sig_mc_b3", "", *massmc_b, *(mcpars[11]), *(mcpars[17]));
  RooAddPdf* sig_mc_b = new RooAddPdf("sig_mc_b", "", RooArgList(*sig_mc_b1, *sig_mc_b2, *sig_mc_b3), RooArgList(*(mcpars[14]), *(mcpars[18])), true);

  /***************************************************
    Start Fitting -- Constraint by MC
  ***************************************************/

  // unbinned fit
  crmc->cd(1);
  std::cout<<std::endl<<"\e[34;2m";
  RooFitResult* fitr_mc_a;
  if(silence) fitr_mc_a = sig_mc_a->fitTo(*dshmc_a, RooFit::Save(), RooFit::SumW2Error(kTRUE), RooFit::PrintEvalErrors(-1));
  else        fitr_mc_a = sig_mc_a->fitTo(*dshmc_a, RooFit::Save(), RooFit::SumW2Error(kTRUE));
  fitX::plotonmc(fremptymc_a, dshmc_a, sig_mc_a, sig_mc_a1, sig_mc_a2, sig_mc_a3, color_a, "__a", option);
  TF1* fsig_mc_a = astfsig(sig_mc_a, f_sig_a, "fsig_mc_a", hmc_a->Integral()*fitX::BIN_WIDTH_L, std::vector<int>({6, 7, 8, 15, 9, 16}));
  fitX::drawpull(hmc_a, fsig_mc_a, color_a);
  fitX::labelsmc("Gen-matched #psi(2S)", mcpars[6]->getVal(), mcpars[7]->getVal(), mcpars[8]->getVal());
  crmc->cd(2);
  RooFitResult* fitr_mc_b; 
  if(silence) fitr_mc_b = sig_mc_b->fitTo(*dshmc_b, RooFit::Save(), RooFit::SumW2Error(kTRUE), RooFit::PrintEvalErrors(-1));
  else        fitr_mc_b = sig_mc_b->fitTo(*dshmc_b, RooFit::Save(), RooFit::SumW2Error(kTRUE));
  fitX::plotonmc(fremptymc_b, dshmc_b, sig_mc_b, sig_mc_b1, sig_mc_b2, sig_mc_b3, color_b, "__b", option);
  TF1* fsig_mc_b = astfsig(sig_mc_b, f_sig_b, "fsig_mc_b", hmc_b->Integral()*fitX::BIN_WIDTH_H, std::vector<int>({11, 12, 13, 17, 14, 18}));
  fitX::drawpull(hmc_b, fsig_mc_b, color_b);
  fitX::labelsmc("Gen-matched X(3872)", mcpars[11]->getVal(), mcpars[12]->getVal(), mcpars[13]->getVal());
  if(saveplot) crmc->SaveAs(Form("%s/chmassmcr%s.pdf", outputdir.c_str(), name.c_str()));
  std::cout<<"\e[0m";
  std::cout<<"\e[34;1m";
  fitr_mc_a->Print("v");
  fitr_mc_b->Print("v");
  std::cout<<"\e[0m";

  // binned fit
  cmc->cd(1);
  hmc_a->Draw("pe");
  xjjroot::settfstyle(f, color_a, 1, 2);
  f->FixParameter(10, 0);
  f->FixParameter(6, 3.68610);
  hmc_a->Fit("f", "Nq");
  hmc_a->Fit("f", "NLLq");
  f->ReleaseParameter(6);
  hmc_a->Fit("f", "NLLq");
  float mean_a = f->GetParameter(6), meanerr_a = f->GetParError(6);
  f->FixParameter(6, f->GetParameter(6));
  f->SetParLimits(7, f->GetParameter(7)-f->GetParError(7), f->GetParameter(7)+f->GetParError(7));
  f->SetParLimits(8, f->GetParameter(8)-f->GetParError(8), f->GetParameter(8)+f->GetParError(8));
  f->SetParLimits(15, f->GetParameter(15)-f->GetParError(15), f->GetParameter(15)+f->GetParError(15));
  if(option != "floatwidth")
    {
      f->FixParameter(7, f->GetParameter(7));
      f->FixParameter(8, f->GetParameter(8));
      f->FixParameter(15, f->GetParameter(15));
    }
  f->FixParameter(9, f->GetParameter(9));
  f->FixParameter(16, f->GetParameter(16));
  f->ReleaseParameter(10);
  fitX::drawpull(hmc_a, f, color_a);
  xjjroot::copyobject(f, "fmc_a")->Draw("same");
  hmc_a->Draw("pesame");
  fitX::labelsmc("Gen-matched #psi(2S)", f->GetParameter(6), f->GetParameter(7), f->GetParameter(8));
  // 
  cmc->cd(2);
  hmc_b->Draw("pe");
  xjjroot::settfstyle(f, color_b, 1, 2);
  f->SetLineColor(color_b);
  f->FixParameter(5, 0);
  f->FixParameter(11, 3.87169);
  hmc_b->Fit("f", "Nq");
  hmc_b->Fit("f", "NLLq");
  f->ReleaseParameter(11);
  hmc_b->Fit("f", "NLLq");
  float mean_b = f->GetParameter(11), meanerr_b = f->GetParError(11);
  f->FixParameter(11, f->GetParameter(11));
  f->SetParLimits(12, f->GetParameter(12)-f->GetParError(12), f->GetParameter(12)+f->GetParError(12));
  f->SetParLimits(13, f->GetParameter(13)-f->GetParError(13), f->GetParameter(13)+f->GetParError(13));
  f->SetParLimits(17, f->GetParameter(17)-f->GetParError(17), f->GetParameter(17)+f->GetParError(17));
  if(option != "floatwidth")
    {
      f->FixParameter(12, f->GetParameter(12));
      f->FixParameter(13, f->GetParameter(13));
      f->FixParameter(17, f->GetParameter(17));
    }
  f->FixParameter(14, f->GetParameter(14));
  f->FixParameter(18, f->GetParameter(18));
  f->ReleaseParameter(5);
  fitX::drawpull(hmc_b, f, color_b);
  xjjroot::copyobject(f, "fmc_b")->Draw("same");
  hmc_b->Draw("pesame");
  fitX::labelsmc("Gen-matched X(3872)", f->GetParameter(11), f->GetParameter(12), f->GetParameter(13));
  cmc->RedrawAxis();
  if(saveplot) cmc->SaveAs(Form("%s/chmassmc%s.pdf", outputdir.c_str(), name.c_str()));

  /***************************************************
    Start Fitting -- Data
  ***************************************************/
  // ---> binned fit
  xjjroot::settfstyle(f, color_data, 1, 3);
  f->ReleaseParameter(0);
  f->ReleaseParameter(1);
  f->ReleaseParameter(2);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  // if(xjjc::str_contains(option,"cheb")) { f->ReleaseParameter(19); f->FixParameter(0, 1); }
  if(option.back()=='3') { f->FixParameter(4, 0); }
  if(option.back()=='2') { f->FixParameter(4, 0); f->FixParameter(3, 0); }
  if(option.back()=='1') { f->FixParameter(4, 0); f->FixParameter(3, 0); f->FixParameter(2, 0); }
  f->SetParLimits(5, 0, 1.e+5);
  f->SetParLimits(10, 0, 1.e+5);
  // if(xjjc::str_contains(option,"cheb"))
  //   {
  //     f->SetParLimits(0, 0, 1.e+10);
  //     f->SetParameter(0, 0);
  //     // !!!
  //   }
  c->cd();
  h->Draw("pe");
  h->Fit("f","Nq");
  h->Fit("f","NLLq");
  h->Fit("f","NLLq");

  cr->cd();
  // ---> unbinned fit
  std::map<int, RooRealVar*> pars;
  pars[5] = new RooRealVar("par5", "", f->GetParameter(5)/fitX::BIN_WIDTH, -1.e+4, 1.e+4); // !!
  pars[6] = new RooRealVar("par6", "", mcpars[6]->getVal(), mcpars[6]->getVal(), mcpars[6]->getVal());
  pars[7] = new RooRealVar("par7", "", mcpars[7]->getVal(), mcpars[7]->getVal()-mcpars[7]->getError(), mcpars[7]->getVal()+mcpars[7]->getError()); if(option!="floatwidth") { pars[7]->setConstant(); }
  pars[8] = new RooRealVar("par8", "", mcpars[8]->getVal(), mcpars[8]->getVal()-mcpars[8]->getError(), mcpars[8]->getVal()+mcpars[8]->getError()); if(option!="floatwidth") { pars[8]->setConstant(); }
  pars[15] = new RooRealVar("par15", "", mcpars[15]->getVal(), mcpars[15]->getVal()-mcpars[15]->getError(), mcpars[15]->getVal()+mcpars[15]->getError()); pars[15]->setConstant();
  pars[9] = new RooRealVar("par9", "", mcpars[9]->getVal(), mcpars[9]->getVal(), mcpars[9]->getVal()); pars[9]->setConstant();
  pars[16] = new RooRealVar("par16", "", mcpars[16]->getVal(), mcpars[16]->getVal(), mcpars[16]->getVal()); pars[16]->setConstant();
  RooGaussian sig_a1("sig_a1", "", *mass, *(pars[6]), *(pars[7]));
  RooGaussian sig_a2("sig_a2", "", *mass, *(pars[6]), *(pars[8]));
  RooGaussian sig_a3("sig_a3", "", *mass, *(pars[6]), *(pars[15]));
  RooAddPdf* sig_a = new RooAddPdf("sig_a", "", RooArgList(sig_a1, sig_a2, sig_a3), RooArgList(*(pars[9]), *(pars[16])), true);
  pars[10] = new RooRealVar("par10", "", f->GetParameter(10)/fitX::BIN_WIDTH, -1.e+4, 1.e+4); // !
  RooRealVar* fsyst = new RooRealVar("fsyst", "", 1, 0, 2); // !
  RooFormulaVar* nsig = new RooFormulaVar("nsig", "@0*@1", RooArgList(*(pars[10]), *fsyst)); // <---------- include systematics
  pars[11] = new RooRealVar("par11", "", fitX::FIT_MASS_X, fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
  pars[12] = new RooRealVar("par12", "", mcpars[12]->getVal(), mcpars[12]->getVal()-mcpars[12]->getError(), mcpars[12]->getVal()+mcpars[12]->getError()); if(option!="floatwidth") { pars[12]->setConstant(); }
  pars[13] = new RooRealVar("par13", "", mcpars[13]->getVal(), mcpars[13]->getVal()-mcpars[13]->getError(), mcpars[13]->getVal()+mcpars[13]->getError()); if(option!="floatwidth") { pars[13]->setConstant(); }
  pars[17] = new RooRealVar("par17", "", mcpars[17]->getVal(), mcpars[17]->getVal()-mcpars[17]->getError(), mcpars[17]->getVal()+mcpars[17]->getError()); pars[17]->setConstant();
  pars[14] = new RooRealVar("par14", "", mcpars[14]->getVal(), mcpars[14]->getVal(), mcpars[14]->getVal()); pars[14]->setConstant();
  pars[18] = new RooRealVar("par18", "", mcpars[18]->getVal(), mcpars[18]->getVal(), mcpars[18]->getVal()); pars[18]->setConstant();
  RooGaussian sig_b1("sig_b1", "", *mass, *(pars[11]), *(pars[12]));
  RooGaussian sig_b2("sig_b2", "", *mass, *(pars[11]), *(pars[13]));
  RooGaussian sig_b3("sig_b3", "", *mass, *(pars[11]), *(pars[17]));
  RooAddPdf* sig_b = new RooAddPdf("sig_b", "", RooArgList(sig_b1, sig_b2, sig_b3), RooArgList(*(pars[14]), *(pars[18])), true);
  if(xjjc::str_contains(option,"cheb"))
    {
      pars[0] = new RooRealVar("par0", "", f->GetParameter(0) / f->GetParameter(0), -1.e+6, 1.e+6);
      pars[1] = new RooRealVar("par1", "", f->GetParameter(1) / f->GetParameter(0), -1.e+6, 1.e+6);
      pars[2] = new RooRealVar("par2", "", f->GetParameter(2) / f->GetParameter(0), -1.e+6, 1.e+6);
      pars[3] = new RooRealVar("par3", "", f->GetParameter(3) / f->GetParameter(0), -1.e+6, 1.e+6);
      pars[4] = new RooRealVar("par4", "", f->GetParameter(4) / f->GetParameter(0), -1.e+6, 1.e+6);
    }
  else
    {
      pars[0] = new RooRealVar("par0", "", f->GetParameter(0), (f->GetParameter(0)<-1.e+5?f->GetParameter(0)*2:-1.e+5), (f->GetParameter(0)>1.e+5?f->GetParameter(0)*2:1.e+5));
      pars[1] = new RooRealVar("par1", "", f->GetParameter(1), (f->GetParameter(1)<-1.e+5?f->GetParameter(1)*2:-1.e+5), (f->GetParameter(1)>1.e+5?f->GetParameter(1)*2:1.e+5));
      pars[2] = new RooRealVar("par2", "", f->GetParameter(2), (f->GetParameter(2)<-1.e+5?f->GetParameter(2)*2:-1.e+5), (f->GetParameter(2)>1.e+5?f->GetParameter(2)*2:1.e+5));
      pars[3] = new RooRealVar("par3", "", f->GetParameter(3), (f->GetParameter(3)<-1.e+5?f->GetParameter(3)*2:-1.e+5), (f->GetParameter(3)>1.e+5?f->GetParameter(3)*2:1.e+5));
      pars[4] = new RooRealVar("par4", "", f->GetParameter(4), (f->GetParameter(4)<-1.e+5?f->GetParameter(4)*2:-1.e+5), (f->GetParameter(4)>1.e+5?f->GetParameter(4)*2:1.e+5));
    }
  if(option.back()=='3') { pars[4]->setVal(0); pars[4]->setConstant(); }
  if(option.back()=='2') { pars[4]->setVal(0); pars[4]->setConstant(); pars[3]->setVal(0); pars[3]->setConstant(); }
  if(option.back()=='1') { pars[4]->setVal(0); pars[4]->setConstant(); pars[3]->setVal(0); pars[3]->setConstant(); pars[2]->setVal(0); pars[2]->setConstant(); }
  RooAbsPdf* bkg;
  if(xjjc::str_contains(option,"cheb")) { bkg = new RooChebychev("bkg", "", *mass, RooArgSet(*(pars[1]), *(pars[2]), *(pars[3]), *(pars[4]))); }
  else                                  { bkg = new RooPolynomial("bkg", "", *mass, RooArgSet(*(pars[0]), *(pars[1]), *(pars[2]), *(pars[3]), *(pars[4])), 0); }
  RooRealVar* nbkg = new RooRealVar("nbkg", "", 0, -1.e+6, 1.e+6);
  RooAddPdf* pdf = new RooAddPdf("pdf", "",
                                 RooArgList(*sig_a    , *sig_b     , *bkg),
                                 RooArgList(*(pars[5]), *(pars[10]), *nbkg)); // <---------- switch w/ and w/o syst
                                 // RooArgList(*(pars[5]), *nsig, *nbkg));  // <---------- switch w/ and w/o syst
  // >>> fix mean
  if(!fixmean)
    {
      f->ReleaseParameter(6);
      f->ReleaseParameter(11);
      f->SetParameter(6, fitX::FIT_MASS_PSI2S);
      f->SetParameter(11, fitX::FIT_MASS_X);
      f->SetParLimits(6, fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      f->SetParLimits(11, fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
      pars[6]->setVal(fitX::FIT_MASS_PSI2S);
      pars[6]->setRange(fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      pars[11]->setVal(fitX::FIT_MASS_X);
      pars[11]->setRange(fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
    }
  else
    { // !! tricky fixmean !! todo !!
      f->ReleaseParameter(6);
      f->SetParameter(6, fitX::FIT_MASS_PSI2S);
      f->SetParLimits(6, fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      pars[6]->setVal(fitX::FIT_MASS_PSI2S);
      pars[6]->setRange(fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      f->FixParameter(11, fixmean);
      pars[11]->setVal(fixmean);
      pars[11]->setConstant();
    }
  // <<<

  c->cd();
  // TFitResultPtr r = histo->Fit(func, "S");
  std::cout<<std::endl<<"\e[35;1m";
  h->Fit("f","NLL", "");
  std::cout<<"\e[0m";
  if(h_ss) h_ss->Draw("pe same");

  cr->cd();
  std::cout<<std::endl<<"\e[35;2m";
  RooFitResult* fitr;
  // RooGaussian fconssyst("fconssyst", "fconssyst", *fsyst, RooFit::RooConst(1), RooFit::RooConst(0.043)); // yield extraction syst
  // fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::Constrain(fconssyst), RooFit::Strategy(0)); // <---------- switch
  // fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::Constrain(fconssyst), RooFit::Strategy(1)); // <---------- switch
  // fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::ExternalConstraints(fconssyst), RooFit::Strategy(0)); // <---------- switch
  // fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::ExternalConstraints(fconssyst), RooFit::Strategy(1)); // <---------- switch
  fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::Strategy(0)); // <---------- switch
  fitr = pdf->fitTo(*dsh, RooFit::Save(), RooFit::PrintEvalErrors(-1), RooFit::Strategy(1)); // <---------- switch
  std::cout<<"\e[0m";
  std::cout<<"\e[35;1m";
  fitr->Print("v");
  std::cout<<"\e[0m";
  TF1* frf = fitX::astf(pdf, f, "fr");

  std::cout<<std::endl<<"\e[7m"<<"==> Finished fitting"<<"\e[0m"<<std::endl<<std::endl;

  /***************************************************
    Plotting
  ***************************************************/
  std::map<std::string, TF1*> fs = fitX::resolvef(f, Form("%s-binned", name.c_str()));
  std::map<std::string, TF1*> fr = fitX::resolvef(frf, Form("%s-unbinned", name.c_str()));

  float maxy = std::max(h->GetMaximum(), frempty->GetMaximum());
  float maxscale = xjjc::str_contains(hh->GetName(), "Benr")?1.5:1.25;
  std::cout<<maxscale<<" "<<hh->GetName()<<std::endl;
  h->SetMaximum(maxy*maxscale);
  frempty->SetMaximum(maxy*maxscale);
  TH1F* hclone = xjjroot::copyobject(h, "hclone");

  // ---> binned fit
  c->cd();
  hclone->Draw("pe");
  fitX::drawpull(h, f, color_data);
  fs["fsig_a"]->Draw("same");
  fs["fsig_b"]->Draw("same");
  fs["fbkg"]->Draw("same");
  fs["f"]->Draw("same");
  hclone->Draw("pe same");
  c->RedrawAxis();
  float ysig_a = fs["fsig_a"]->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  float ysig_b = fs["fsig_b"]->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  if(!fixmean)
    {
      xjjroot::drawbox(fitX::PDG_MASS_PSI2S - fitX::PDG_MASS_PSI2S_ERR, h->GetMinimum(), fitX::PDG_MASS_PSI2S + fitX::PDG_MASS_PSI2S_ERR, h->GetMaximum(), color_a, 0.8, 1001, 0, 0, 0);
      xjjroot::drawbox(fitX::PDG_MASS_X - fitX::PDG_MASS_X_ERR, h->GetMinimum(), fitX::PDG_MASS_X + fitX::PDG_MASS_X_ERR, h->GetMaximum(), color_b, 0.8, 1001, 0, 0, 0);
    }
  fitX::labelsdata(title, cmsleft, cmsright);
  if(saveplot) c->SaveAs(Form("%s/chmass%s.pdf", outputdir.c_str(), name.c_str()));
  fitX::labelsdata_params(f->GetParameter(6), f->GetParError(6), ysig_a, f->GetParError(5)*ysig_a/f->GetParameter(5),
                          f->GetParameter(11), f->GetParError(11), ysig_b, f->GetParError(10)*ysig_b/f->GetParameter(10),
                          TMath::Prob(f->GetChisquare(), f->GetNDF()));
  if(saveplot) c->SaveAs(Form("%s/chmass%s_params.pdf", outputdir.c_str(), name.c_str()));

  // ---> no fit histogram
  crnofit->cd();
  dsh->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  frempty->Draw();
  crnofit->RedrawAxis();
  fitX::labelsdata(title, cmsleft, cmsright);
  if(saveplot) crnofit->SaveAs(Form("%s/chmassrnofit%s.pdf", outputdir.c_str(), name.c_str()));

  // ---> unbinned fit
  cr->cd();
  dsh->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  pdf->plotOn(frempty, RooFit::Name("sig_a"), RooFit::Components(*sig_a), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_a), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("sig_b"), RooFit::Components(*sig_b), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_b), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("bkg"), RooFit::Components(*bkg), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_data), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_data), RooFit::LineWidth(3));
  dsh->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  frempty->Draw();
  fitX::drawpull(h, frf, color_data);
  frf->SetLineColor(color_data);
  frf->Draw("same");
  hclone->Draw("same pe");
  int32_t ndof = fitr->floatParsFinal().getSize();
  RooRealVar rndof("ndof", "", ndof);
  float minNll = fitr->minNll();
  std::cout<<"\e[32;7m"<<minNll<<"\e[0m"<<std::endl;
  cr->RedrawAxis();
  fitX::labelsdata(title, cmsleft, cmsright);
  if(saveplot) cr->SaveAs(Form("%s/chmassr%s.pdf", outputdir.c_str(), name.c_str()));
  fitX::labelsdata_params(pars[6]->getVal(), pars[6]->getError(), pars[5]->getVal(), pars[5]->getError(),
                          pars[11]->getVal(), pars[11]->getError(), pars[10]->getVal(), pars[10]->getError(), // <---------- switch
                          // pars[11]->getVal(), pars[11]->getError(), pars[10]->getVal()*fsyst->getVal(), pars[10]->getError()*fsyst->getVal(), // <---------- switch
                          TMath::Prob(frempty->chiSquare("pdf", "dshist", ndof)*ndof, ndof));
  if(saveplot) cr->SaveAs(Form("%s/chmassr%s_params.pdf", outputdir.c_str(), name.c_str()));

  printfit(f, pdf);

  // ---> for sPlot
  pars[6]->setConstant();
  pars[11]->setConstant();
  pars[0]->setConstant();
  pars[1]->setConstant();
  pars[2]->setConstant();
  pars[3]->setConstant();
  pars[4]->setConstant();

  RooStats::SPlot* sData = new RooStats::SPlot("sData", "An SPlot", *dsh, pdf, RooArgList(*(pars[5]), *(pars[10]), *nbkg));
  // RooStats::SPlot* sData = new RooStats::SPlot("sData", "An SPlot", *dsh, pdf, RooArgList(*(pars[5]), *(nsig), *nbkg));
  ww->import(*dsh, RooFit::Rename(Form("%s_ws", dshh->GetName())), RooFit::Silence());
  ww->import(*pdf, RooFit::Rename(Form("pdf_%s", name.c_str())), RooFit::Silence());
  ww->import(rndof);

  delete cmc;
  delete c;
  delete crmc;
  delete cr;

  fitresult["binned"] = new fitX::fitXresult(ysig_a, ysig_b, f->GetParError(5)*ysig_a/f->GetParameter(5), f->GetParError(10)*ysig_a/f->GetParameter(10), 
                                             f->GetParameter(6), f->GetParameter(11), f->GetParError(6), f->GetParError(11));
  fitresult["binned"]->setf(fs["f"], fs["fsig_a"], fs["fsig_b"], fs["fbkg"]);
  fitresult["binned"]->setw(0, 0, new TH2F(Form("hcorr-binned-%s", name.c_str()), ";;", 1, 0, 1, 1, 0, 1));
  fitresult["unbinned"] = new fitX::fitXresult(pars[5]->getVal(), pars[10]->getVal()*fsyst->getVal(), pars[5]->getError(), pars[10]->getError()*fsyst->getVal(),
                                               pars[6]->getVal(), pars[11]->getVal(), pars[6]->getError(), pars[11]->getError());
  fitresult["unbinned"]->setf(fr["f"], fr["fsig_a"], fr["fsig_b"], fr["fbkg"]);
  fitresult["unbinned"]->setw(ww, minNll, (TH2F*)fitr->correlationHist(Form("hcorr-unbinned-%s", name.c_str())));

  return fitresult;
}

void fitX::setmasshist(TH1* h, float xoffset/*=0*/, float yoffset/*=0*/, Color_t pcolor/*=kBlack*/)
{
  xjjroot::sethempty(h, xoffset, yoffset);
  h->Sumw2();
  h->SetLineColor(pcolor);
  h->SetMarkerColor(pcolor);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.0);
}

void fitX::setmasshist(RooPlot* h, float xoffset/*=0*/, float yoffset/*=0*/, Color_t pcolor/*=kBlack*/)
{
  xjjroot::sethempty(h, xoffset, yoffset);
  h->SetTitle("");
  h->SetXTitle("m_{#mu#mu#pi#pi} (GeV/c^{2})");
}

void fitX::drawpull(TH1* hmc, TF1* f, Color_t color)
{
  int nbin = hmc->GetXaxis()->GetNbins();
  float binmin = hmc->GetBinCenter(1) - hmc->GetBinWidth(1)/2.;
  float binmax = hmc->GetBinCenter(nbin) + hmc->GetBinWidth(nbin)/2.;
  for(int bb=0; bb<nbin; bb++)
    { 
      float realval = hmc->GetBinError(bb+1)==0?0:(hmc->GetBinContent(bb+1)-f->Eval(hmc->GetBinCenter(bb+1)))/hmc->GetBinError(bb+1);
      float fillval = ((realval+4)/(4*2))*hmc->GetMaximum();
      xjjroot::drawbox(hmc->GetBinCenter(bb+1)-hmc->GetBinWidth(bb+1)/2., hmc->GetMaximum()/2., hmc->GetBinCenter(bb+1)+hmc->GetBinWidth(bb+1)/2., fillval, color, 0.1, 1001);
    }
  xjjroot::drawline(binmin, hmc->GetMaximum()/2., binmax, hmc->GetMaximum()/2., kGray, 2, 2, 0.5);
  xjjroot::drawaxis(binmax, 0, binmax, hmc->GetMaximum(), -4, 4, color, 1, 2, "+L");
  xjjroot::drawtex(0.93, 0.51, "Pull", 0.04, 33, 62, color);
  xjjroot::drawtex(1.00, 0.53, "#sigma", 0.04, 33, 62, color);
}

void fitX::labelsmc(std::string label, double mean, double sigma1, double sigma2)
{
  xjjroot::drawtex(0.18, 0.85, label.c_str(), 0.04, 12, 62);
  xjjroot::drawtex(0.18, 0.85-0.05, pttag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*2, ytag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*3, centtag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", mean), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma_{1} = %.4f GeV", std::min(sigma1, sigma2)), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05*2, Form("#sigma_{2} = %.4f GeV", std::max(sigma1, sigma2)), 0.038);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
}

void fitX::labelsdata(std::string title, std::string cmsleft, std::string cmsright)
{
  xjjroot::drawtex(0.92, 0.84, Form("%.0f < p_{T} < %.0f GeV/c", fitX::ptmincut, fitX::ptmaxcut), 0.038, 32, 62);
  xjjroot::drawtex(0.92, 0.79, Form("%s|y| < %.1f", (fitX::ymincut?Form("%.1f < ", fitX::ymincut):""), fitX::ymaxcut), 0.038, 32, 62);
  xjjroot::drawtex(0.92, 0.74, Form("Cent. %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%"), 0.038, 32, 62);
  // xjjroot::drawtex(0.17, 0.84-0.04, title.c_str(), 0.038, 13, 62);
  xjjroot::drawtex(0.17, 0.86, title.c_str(), 0.045, 13, 52);

  xjjroot::drawCMSleft(); //preliminary
  // xjjroot::drawCMSleft(cmsleft.c_str());
  xjjroot::drawCMSright(cmsright.c_str());
}

void fitX::labelsdata_params(double mean_a, double mean_a_err, double yield_a, double yield_a_err,
                             double mean_b, double mean_b_err, double yield_b, double yield_b_err,
                             double chi2prob)
{
  xjjroot::drawtex(0.32, 0.36, Form("#bar{m}_{#psi(2S)} = %.4f #pm %.4f GeV", mean_a, mean_a_err), 0.03, 12, 62, color_a);
  xjjroot::drawtex(0.32, 0.31, Form("N_{#psi(2S)} = %.0f #pm %.0f", yield_a, yield_a_err), 0.03, 12, 62, color_a);
  xjjroot::drawtex(0.32, 0.26, Form("#bar{m}_{X(3872)} = %.4f #pm %.4f GeV", mean_b, mean_b_err), 0.03, 12, 62, color_b);
  xjjroot::drawtex(0.32, 0.21, Form("N_{X(3872)} = %.0f #pm %.0f", yield_b, yield_b_err), 0.03, 12, 62, color_b);

  xjjroot::drawtex(0.17, 0.84-0.04*2, Form("#chi^{2} Prob = %.1f%s", chi2prob*100., "%"), 0.042, 12);
}

void fitX::printfit(TF1* f, RooAbsPdf* pdf)
{
  int stdwid = 15;
  std::map<std::string, TF1*> frr = resolvef(f, "tmp");
  float ysig_a = frr["fsig_a"]->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  float ysig_b = frr["fsig_b"]->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  float ysigerr_a = f->GetParError(5)*ysig_a/f->GetParameter(5);
  float ysigerr_b = f->GetParError(10)*ysig_b/f->GetParameter(10);
  RooArgSet* variables = pdf->getVariables();

  std::cout << std::left         << "\e[33;1m"   << std::endl         << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "Val"        << std::setw(stdwid) << "Binned Fit"        << std::setw(stdwid) << "Unbinned Fit"                                             << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "m(psi')"    << std::setw(stdwid) << f->GetParameter(6)  << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 6)))->getVal()    << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "errm(psi')" << std::setw(stdwid) << f->GetParError(6)   << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 6)))->getError()  << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "m(X)"       << std::setw(stdwid) << f->GetParameter(11) << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 11)))->getVal()   << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "errm(X)"    << std::setw(stdwid) << f->GetParError(11)  << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 11)))->getError() << std::endl << std::string(stdwid*3, '-') << std::endl
            << "\e[33;7m"
            << std::setw(stdwid) << "N(psi')"    << std::setw(stdwid) << ysig_a              << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 5)))->getVal()    << "\e[0m"<<"\e[33;1m" << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "errN(psi')" << std::setw(stdwid) << ysigerr_a           << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 5)))->getError()  << std::endl << std::string(stdwid*3, '-') << std::endl
            << "\e[33;7m"
            << std::setw(stdwid) << "N(X)"       << std::setw(stdwid) << ysig_b              << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 10)))->getVal()   << "\e[0m"<<"\e[33;1m" << std::endl << std::string(stdwid*3, '-') << std::endl
            << std::setw(stdwid) << "errN(X)"    << std::setw(stdwid) << ysigerr_b           << std::setw(stdwid) << ((RooRealVar*)variables->find(Form("par%d", 10)))->getError() << std::endl << std::string(stdwid*3, '-') << std::endl
            << "\e[0m"           << std::endl;
}

void fitX::zeroparameters(TF1* f, std::vector<int> ipars)
{
  for(auto& i : ipars)
    f->SetParameter(i, 0);
}

std::map<std::string, TF1*> fitX::resolvef(TF1* f, std::string name)
{
  std::map<std::string, TF1*> fs;

  fs["f"] = xjjroot::copyobject(f, Form("f-%s", name.c_str()));

  fs["fsig_a"] = xjjroot::copyobject(f, Form("fsig_a-%s", name.c_str()));
  fs["fsig_a"]->SetRange(BIN_MIN_L, BIN_MAX_L);
  xjjroot::settfstyle(fs["fsig_a"], color_a, 7, 3, 0, 0, 1001);
  fitX::zeroparameters(fs["fsig_a"], std::vector<int>({10, 0, 1, 2, 3, 4}));

  fs["fsig_b"] = xjjroot::copyobject(f, Form("fsig_b-%s", name.c_str()));
  fs["fsig_b"]->SetRange(BIN_MIN_H, BIN_MAX_H);
  xjjroot::settfstyle(fs["fsig_b"], color_b, 7, 3, 0, 0, 1001);
  fitX::zeroparameters(fs["fsig_b"], std::vector<int>({5, 0, 1, 2, 3, 4}));

  fs["fbkg"] = xjjroot::copyobject(f, Form("fbkg-%s", name.c_str()));
  fs["fbkg"]->SetRange(BIN_MIN, BIN_MAX);
  xjjroot::settfstyle(fs["fbkg"], color_bkg, 7, 3, 0, 0, 1001);
  fitX::zeroparameters(fs["fbkg"], std::vector<int>({5, 10}));

  return fs;
}

void fitX::plotonmc(RooPlot* fremptymc, RooDataSet* dshmc, RooAbsPdf* sig_mc, RooAbsPdf* sig_mc1, RooAbsPdf* sig_mc2, RooAbsPdf* sig_mc3, Color_t color, std::string name, std::string option)
{
  dshmc->plotOn(fremptymc, RooFit::Name(Form("dshmc%s", name.c_str())), RooFit::Binning(fitX::NBIN_L), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  sig_mc->plotOn(fremptymc, RooFit::Name(Form("sig_mc1%s", name.c_str())), RooFit::Components(*sig_mc1), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color), RooFit::LineWidth(2));
  sig_mc->plotOn(fremptymc, RooFit::Name(Form("sig_mc2%s", name.c_str())), RooFit::Components(*sig_mc2), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color), RooFit::LineWidth(2));
  if(option=="3gaus") sig_mc->plotOn(fremptymc, RooFit::Name(Form("sig_mc3%s", name.c_str())), RooFit::Components(*sig_mc3), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color), RooFit::LineWidth(2));
  sig_mc->plotOn(fremptymc, RooFit::Name(Form("sig_mc%s", name.c_str())), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color), RooFit::LineWidth(2));
  dshmc->plotOn(fremptymc, RooFit::Name(Form("dshmc%s", name.c_str())), RooFit::Binning(fitX::NBIN_L), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  fremptymc->Draw();    
}

// TF1* frf = (TF1*)(pdf->asTF(RooArgList(*mass))->Clone("frf"));
TF1* fitX::astf(RooAbsPdf* pdf, TF1* f, std::string name, float width)
{
  TF1* frf = (TF1*)f->Clone(name.c_str());
  RooArgSet* variables = pdf->getVariables();
  RooArgSet* components = pdf->getComponents();
  // components->Print();
  RooArgSet* params = pdf->getParameters((*variables)["Bmass"]);
  std::cout<<std::endl;
  // params->Print("v");
  double polynorm = ((RooAbsPdf*)(components->find("bkg")))->getNormIntegral((*variables)["Bmass"])->getVal();
  for(int i=0;i<f->GetNpar();i++)
    {
      // float ipar = i?params->getRealValue(Form("par%d",i)):1; // roopoly par0 = 1
      float ipar = params->getRealValue(Form("par%d",i));
      if(i<=4) ipar *= (params->getRealValue("nbkg")/polynorm);
      if(i<=4 || i==5 || i==10) ipar *= width;
      frf->SetParameter(i, ipar);
    }
  return frf;
}

TF1* fitX::astfsig(RooAbsPdf* pdf, TF1* f, std::string name, float norm, std::vector<int> pars)
{
  TF1* frf = (TF1*)f->Clone(name.c_str());
  RooArgSet* variables = pdf->getVariables();
  RooArgSet* components = pdf->getComponents();
  // components->Print();
  RooArgSet* params = pdf->getParameters((*variables)["Bmass"]);
  // std::cout<<"\e[31;1m"<<std::endl;
  // variables->Print("v");
  // components->Print("v");
  // params->Print("v");
  // double norm = pdf->getNormIntegral((*variables)["Bmass"])->getVal();
  // double norm = pdf->getNorm((*variables)["Bmass"]);
  // double norm = pdf->getVal();
  frf->SetParameter(f->GetNpar()-1, norm);
  for(int i=0;i<f->GetNpar()-1;i++)
    {
      float ipar = params->getRealValue(Form("mcpar%d",pars[i]));
      frf->SetParameter(i, ipar);
    }
  // std::cout<<frf->GetExpFormula("P CLING")<<std::endl;
  // std::cout<<"\e[0m";
  return frf;
}

#endif


  // RooStats::ModelConfig* mcpdf = new RooStats::ModelConfig("ModelConfig", ww);
  // mcpdf->SetPdf(*pdf);
  // mcpdf->SetParametersOfInterest(*(ww->var("par10")));
  // mcpdf->SetObservables(*(ww->var("Bmass")));
  // ww->import(*mcpdf);
