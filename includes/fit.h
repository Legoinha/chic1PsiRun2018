#ifndef __FITX_FIT_H_
#define __FITX_FIT_H_

#include "xjjrootuti.h"
#include "xjjcuti.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
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
#include <RooFitResult.h>
#include <RooPlot.h>

#include <string>

#include "fitX.h"

/* ----------------------------------------
// ==> Usage <==
#include "fitX.h"

// hdata: data mass; hmc_psi, hmc_x: mc template;
std::vector<TF1*> funs = fitX::fit(hdata, 0, hmc_a, hmc_b, "plots/", true, true, "name"); 
// psi' yield:
float ysig1 = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
float ysig1err = funs[0]->GetParError(5)*ysig1/funs[0]->GetParameter(5);
// X(3872) yield:
float ysig2 = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
float ysig2err = funs[0]->GetParError(10)*ysig2/funs[0]->GetParameter(10);

---------------------------------------- */

namespace fitX
{
  std::vector<TF1*> fit(TH1F* hh, TH1F* hh_ss, TH1F* hhmc_a, TH1F* hhmc_b, 
                        RooDataSet* dshh, RooDataSet* dshhmc_a, RooDataSet* dshhmc_b,
                        std::string outputdir, bool fixmean, bool saveplot, std::string name="", std::string option="default");
  void setmasshist(TH1* h, float xoffset=0, float yoffset=0, Color_t pcolor=kBlack);
  void setmasshist(RooPlot* h, float xoffset=0, float yoffset=0, Color_t pcolor=kBlack);

  const int NBIN = 38, NBIN_L = 50, NBIN_H = 50;
  // const int NBIN = 57, NBIN_L = 50, NBIN_H = 50;
  const float BIN_MIN = 3.62, BIN_MAX = 4.0, BIN_MIN_L = 3.64, BIN_MAX_L = 3.74, BIN_MIN_H = 3.82, BIN_MAX_H = 3.92;
  float BIN_WIDTH = (BIN_MAX-BIN_MIN)/NBIN*1.0, BIN_WIDTH_L = (BIN_MAX_L-BIN_MIN_L)/NBIN_L*1.0, BIN_WIDTH_H = (BIN_MAX_H-BIN_MIN_H)/NBIN_H*1.0;
  void drawpull(TH1* hmc, TF1* f, Color_t color);

  Color_t color_data = kRed-3, color_a = kAzure+4, color_b = kGreen-1, color_ss = kGray+1, color_bkg = color_data;
  int ibin_a = 2, ibin_b = 4, nbin = 5;
  std::string title_a = "#psi(2S)", title_b = "X(3872)";

  double getparmin(TF1* f, int ipar) { double parmin, parmax; f->GetParLimits(ipar, parmin, parmax); return parmin; } 
  double getparmax(TF1* f, int ipar) { double parmin, parmax; f->GetParLimits(ipar, parmin, parmax); return parmax; } 
  
}

// --->
std::vector<TF1*> fitX::fit(TH1F* hh, TH1F* hh_ss, TH1F* hhmc_a, TH1F* hhmc_b, 
                            RooDataSet* dshh, RooDataSet* dshhmc_a, RooDataSet* dshhmc_b, 
                            std::string outputdir, bool fixmean, bool saveplot, std::string name, std::string option)
{
  std::string uniqstr(xjjc::currenttime());
  if(saveplot) gSystem->Exec(Form("mkdir -p %s", outputdir.c_str()));

  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooRealVar* massmc_a = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_L, fitX::BIN_MAX_L);
  RooRealVar* massmc_b = new RooRealVar("Bmass", "massmc_a", fitX::BIN_MIN_H, fitX::BIN_MAX_H);
  RooPlot* frempty = mass->frame(RooFit::Title(""));
  frempty->SetYTitle(Form("Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3));
  RooPlot* fremptymc_a = massmc_a->frame(RooFit::Title(""));
  fremptymc_a->SetYTitle(Form("Entries / %.0f MeV", fitX::BIN_WIDTH_L*1.e+3));
  RooPlot* fremptymc_b = massmc_b->frame(RooFit::Title(""));
  fremptymc_b->SetYTitle(Form("Entries / %.0f MeV", fitX::BIN_WIDTH_H*1.e+3));

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
  TCanvas* crmc = new TCanvas("crmc", "", 1200, 600);
  crmc->Divide(2, 1);
  
  TString str_bkg = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x";
  if(xjjc::str_contains(option,"cheb")) { str_bkg = option; }
  TString str_sig_a = "[5]*([9]*TMath::Gaus(x,[6],[7])/(TMath::Sqrt(2*3.14159)*[7]) + (1-[9])*TMath::Gaus(x,[6],[8])/(TMath::Sqrt(2*3.14159)*[8]))";
  TString str_sig_b = "[10]*([14]*TMath::Gaus(x,[11],[12])/(TMath::Sqrt(2*3.14159)*[7]) + (1-[14])*TMath::Gaus(x,[11],[13])/(TMath::Sqrt(2*3.14159)*[13]))";
  TF1* f = new TF1("f", str_bkg+"+"+str_sig_a+"+"+str_sig_b, BIN_MIN, BIN_MAX);
  const double parinit[15] = {0., 0.     , 0.     , 0.  , 0. , // 0-4
                              1., 3.686  , 0.00357, 0.005, 0.5, // 5-9
                              1., 3.87169, 0.001  , 0.008, 0.5, // 10-14
  };

  f->SetNpx(1000);
  f->SetParameters(parinit);
  f->SetParLimits(9, 0, 1);
  f->SetParLimits(14, 0, 1);
  f->SetParLimits(7, 1.e-6, 0.01);
  f->SetParLimits(8, 1.e-6, 0.01);
  f->SetParLimits(12, 1.e-6, 0.01);
  f->SetParLimits(13, 1.e-6, 0.01);
  f->FixParameter(0, 0);
  f->FixParameter(1, 0);
  f->FixParameter(2, 0);
  f->FixParameter(3, 0);
  f->FixParameter(4, 0);

  std::map<int, RooRealVar*> mcpars;
  mcpars[5] = new RooRealVar("mcpar5", "", parinit[5], 0, 1.e+5);
  mcpars[6] = new RooRealVar("mcpar6", "", parinit[6], fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
  mcpars[7] = new RooRealVar("mcpar7", "", parinit[7], 1.e-6, 0.01);
  mcpars[8] = new RooRealVar("mcpar8", "", parinit[8], 1.e-6, 0.01);
  mcpars[9] = new RooRealVar("mcpar9", "", parinit[9], 0, 1);
  RooGaussian sig_mc_a1("sig_mc_a1", "", *massmc_a, *(mcpars[6]), *(mcpars[7]));
  RooGaussian sig_mc_a2("sig_mc_a2", "", *massmc_a, *(mcpars[6]), *(mcpars[8]));
  RooAddPdf* sig_mc_a = new RooAddPdf("sig_mc_a", "", RooArgList(sig_mc_a1, sig_mc_a2), *(mcpars[9]), true);
  // RooAddPdf* pdf_mc_a = new RooAddPdf("pdf_mc_a", "", RooArgList(*sig_mc_a), RooArgList(*(mcpars[5])));
  std::cout<<"test0 "<<sig_mc_a<<std::endl;
  RooFitResult* fitr_mc_a = sig_mc_a->fitTo(*dshmc_a, RooFit::Save());
  std::cout<<"test"<<std::endl;
  mcpars[10] = new RooRealVar("mcpar10", "", parinit[10], 0, 1.e+5);
  mcpars[11] = new RooRealVar("mcpar11", "", parinit[11], fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
  mcpars[12] = new RooRealVar("mcpar12", "", parinit[12], getparmin(f, 12), getparmax(f, 12));
  mcpars[13] = new RooRealVar("mcpar13", "", parinit[13], getparmin(f, 13), getparmax(f, 13));
  mcpars[14] = new RooRealVar("mcpar14", "", parinit[14], getparmin(f, 14), getparmax(f, 14));
  RooGaussian sig_mc_b1("sig_mc_b1", "", *massmc_b, *(mcpars[11]), *(mcpars[12]));
  RooGaussian sig_mc_b2("sig_mc_b2", "", *massmc_b, *(mcpars[11]), *(mcpars[13]));
  RooAddPdf* sig_mc_b = new RooAddPdf("sig_mc_b", "", RooArgList(sig_mc_b1, sig_mc_b2), *(mcpars[14]), true);
  // RooAddPdf* pdf_mc_b = new RooAddPdf("pdf_mc_b", "", RooArgList(*sig_mc_b), RooArgList(*(mcpars[10])));
  RooFitResult* fitr_mc_b = sig_mc_b->fitTo(*dshmc_b, RooFit::Save());

  crmc->cd(1);
  dshmc_a->plotOn(fremptymc_a, RooFit::Name("dshmc_a"), RooFit::Binning(fitX::NBIN_L), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  sig_mc_a->plotOn(fremptymc_a, RooFit::Name("sig_mc_a1"), RooFit::Components(sig_mc_a1), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_a), RooFit::LineWidth(2));
  sig_mc_a->plotOn(fremptymc_a, RooFit::Name("sig_mc_a2"), RooFit::Components(sig_mc_a2), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_a), RooFit::LineWidth(2));
  sig_mc_a->plotOn(fremptymc_a, RooFit::Name("sig_mc_a"), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_a), RooFit::LineWidth(2));
  dshmc_a->plotOn(fremptymc_a, RooFit::Name("dshmc_a"), RooFit::Binning(fitX::NBIN_L), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  fremptymc_a->Draw();
  crmc->cd(2);
  dshmc_b->plotOn(fremptymc_b, RooFit::Name("dshmc_b"), RooFit::Binning(fitX::NBIN_H), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  sig_mc_b->plotOn(fremptymc_b, RooFit::Name("sig_mc_b1"), RooFit::Components(sig_mc_b1), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_b), RooFit::LineWidth(2));
  sig_mc_b->plotOn(fremptymc_b, RooFit::Name("sig_mc_b2"), RooFit::Components(sig_mc_b2), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_b), RooFit::LineWidth(2));
  sig_mc_b->plotOn(fremptymc_b, RooFit::Name("sig_mc_b"), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_b), RooFit::LineWidth(2));
  dshmc_b->plotOn(fremptymc_b, RooFit::Name("dshmc_b"), RooFit::Binning(fitX::NBIN_H), RooFit::MarkerSize(1.), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(2));
  fremptymc_b->Draw();

  //
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
  f->FixParameter(7, f->GetParameter(7));
  f->FixParameter(8, f->GetParameter(8));
  f->FixParameter(9, f->GetParameter(9));
  f->ReleaseParameter(10);
  fitX::drawpull(hmc_a, f, color_a);
  xjjroot::copyobject(f, "fmc_a")->Draw("same");
  hmc_a->Draw("pesame");
  xjjroot::drawtex(0.18, 0.85, "Gen-matched #psi(2S)", 0.04, 12, 62);
  xjjroot::drawtex(0.18, 0.85-0.05, pttag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*2, ytag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*3, centtag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", f->GetParameter(6)), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma_{1} = %.4f GeV", std::min(f->GetParameter(7),f->GetParameter(8))), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05*2, Form("#sigma_{2} = %.4f GeV", std::max(f->GetParameter(7),f->GetParameter(8))), 0.038);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
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
  f->FixParameter(12, f->GetParameter(12));
  f->FixParameter(13, f->GetParameter(13));
  f->FixParameter(14, f->GetParameter(14));
  f->ReleaseParameter(5);
  fitX::drawpull(hmc_b, f, color_b);
  xjjroot::copyobject(f, "fmc_b")->Draw("same");
  hmc_b->Draw("pesame");
  xjjroot::drawtex(0.18, 0.85, "Gen-matched X(3872)", 0.042, 12, 62);
  xjjroot::drawtex(0.18, 0.85-0.05, pttag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*2, ytag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.18, 0.85-0.05*3, centtag().c_str(), 0.04, 12, 42);
  xjjroot::drawtex(0.65, 0.86, Form("#bar{m} = %.4f GeV", f->GetParameter(11)), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05, Form("#sigma_{1} = %.4f GeV", std::min(f->GetParameter(12),f->GetParameter(13))), 0.038);
  xjjroot::drawtex(0.65, 0.86-0.05*2, Form("#sigma_{2} = %.4f GeV", std::max(f->GetParameter(12),f->GetParameter(13))), 0.038);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  cmc->RedrawAxis();
  if(saveplot) cmc->SaveAs(Form("%s/chmassmc%s.pdf", outputdir.c_str(), name.c_str()));
  if(saveplot) crmc->SaveAs(Form("%s/chmassmcr%s.pdf", outputdir.c_str(), name.c_str()));

  //
  xjjroot::settfstyle(f, color_data, 1, 3);
  f->ReleaseParameter(0);
  f->ReleaseParameter(1);
  f->ReleaseParameter(2);
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  if(option.back()=='3') { f->FixParameter(4, 0); }
  if(option.back()=='2') { f->FixParameter(4, 0); f->FixParameter(3, 0); }

  // >>>
  if(!fixmean)
    {
      f->ReleaseParameter(6);
      f->ReleaseParameter(11);
      f->SetParameter(6, fitX::FIT_MASS_PSI2S);
      f->SetParameter(11, fitX::FIT_MASS_X);
      f->SetParLimits(6, fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      f->SetParLimits(11, fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
    }

  // <<<
  f->SetParLimits(5, 0, 1.e+5);
  f->SetParLimits(10, 0, 1.e+5);

  std::map<int, RooRealVar*> pars;
  pars[5] = new RooRealVar("par5", "", mcpars[5]->getVal(), 0, 1.e+5);
  pars[6] = new RooRealVar("par6", "", mcpars[6]->getVal(), mcpars[6]->getVal(), mcpars[6]->getVal());
  pars[7] = new RooRealVar("par7", "", mcpars[7]->getVal(), mcpars[7]->getVal(), mcpars[7]->getVal()); pars[7]->setConstant();
  pars[8] = new RooRealVar("par8", "", mcpars[8]->getVal(), mcpars[8]->getVal(), mcpars[8]->getVal()); pars[8]->setConstant();
  pars[9] = new RooRealVar("par9", "", mcpars[9]->getVal(), mcpars[9]->getVal(), mcpars[9]->getVal()); pars[9]->setConstant();
  RooGaussian sig_a1("sig_a1", "", *mass, *(pars[6]), *(pars[7]));
  RooGaussian sig_a2("sig_a2", "", *mass, *(pars[6]), *(pars[8]));
  RooAddPdf* sig_a = new RooAddPdf("sig_a", "", RooArgList(sig_a1, sig_a2), *(pars[9]), true);
  pars[10] = new RooRealVar("par10", "", mcpars[10]->getVal(), 0, 1.e+5);
  pars[11] = new RooRealVar("par11", "", mcpars[11]->getVal(), mcpars[11]->getVal(), mcpars[11]->getVal());
  pars[12] = new RooRealVar("par12", "", mcpars[12]->getVal(), mcpars[12]->getVal(), mcpars[12]->getVal()); pars[12]->setConstant();
  pars[13] = new RooRealVar("par13", "", mcpars[13]->getVal(), mcpars[13]->getVal(), mcpars[13]->getVal()); pars[13]->setConstant();
  pars[14] = new RooRealVar("par14", "", mcpars[14]->getVal(), mcpars[14]->getVal(), mcpars[14]->getVal()); pars[14]->setConstant();
  RooGaussian sig_b1("sig_b1", "", *mass, *(pars[11]), *(pars[12]));
  RooGaussian sig_b2("sig_b2", "", *mass, *(pars[11]), *(pars[13]));
  RooAddPdf* sig_b = new RooAddPdf("sig_b", "", RooArgList(sig_b1, sig_b2), *(pars[14]), true);
  pars[0] = new RooRealVar("par0", "", 0, -1.e+5, 1.e+5); // parlimits !!! todo
  pars[1] = new RooRealVar("par1", "", 0, -1.e+5, 1.e+5);
  pars[2] = new RooRealVar("par2", "", 0, -1.e+5, 1.e+5);
  pars[3] = new RooRealVar("par3", "", 0, -1.e+5, 1.e+5);
  pars[4] = new RooRealVar("par4", "", 0, -1.e+5, 1.e+5);
  RooPolynomial bkg_poly("bkg_poly", "", *mass, RooArgSet(*(pars[0]), *(pars[1]), *(pars[2]), *(pars[3]), *(pars[4])));
  RooChebychev bkg_cheb("bkg_cheb", "", *mass, RooArgSet(*(pars[0]), *(pars[1]), *(pars[2]), *(pars[3]), *(pars[4])));
  if(option.back()=='3') { pars[4]->setConstant(); }
  if(option.back()=='2') { pars[4]->setConstant(); pars[3]->setConstant(); }
  RooRealVar* nbkg = new RooRealVar("nbkg", "", 1, 0, 1.e+6);
  RooAddPdf* pdf;
  if(xjjc::str_contains(option,"cheb")) 
    { pdf = new RooAddPdf("pdf", "", 
                          // RooArgList(*sig_a    , *sig_b     , bkg_cheb)); }
                          RooArgList(*sig_a    , *sig_b     , bkg_cheb),
                          RooArgList(*(pars[5]), *(pars[10]), *nbkg)); }
  else 
    { pdf = new RooAddPdf("pdf", "",
                          // RooArgList(*sig_a    , *sig_b     , bkg_poly)); }
                          RooArgList(*sig_a    , *sig_b     , bkg_poly), 
                          RooArgList(*(pars[5]), *(pars[10]), *nbkg)); }
  if(!fixmean)
    {
      pars[6]->setVal(fitX::FIT_MASS_PSI2S);
      pars[6]->setRange(fitX::FIT_MASS_PSI2S - fitX::FIT_MASS_PSI2S_WIN, fitX::FIT_MASS_PSI2S + fitX::FIT_MASS_PSI2S_WIN);
      pars[11]->setVal(fitX::FIT_MASS_X);
      pars[11]->setRange(fitX::FIT_MASS_X - fitX::FIT_MASS_X_WIN, fitX::FIT_MASS_X + fitX::FIT_MASS_X_WIN);
    }
  else
    {
      pars[6]->setConstant();
      pars[11]->setConstant();
    }
  RooFitResult* fitr = pdf->fitTo(*dsh, RooFit::Save());
  TF1* roof = xjjroot::copyobject(f, "roof");
  for(int i=0; i<15; i++) 
    { roof->SetParameter(i, pars[i]->getVal()); 
      roof->SetParError(i, pars[i]->getError()); }
  roof->SetLineColor(color_data);
  roof->SetLineStyle(9);
  
  c->cd();
  f->SetLineColor(color_data);
  h->SetMaximum(h->GetMaximum()*1.25);
  frempty->SetMaximum(h->GetMaximum());
  h->Draw("pe");
  if(!fixmean)
    {
      xjjroot::drawbox(fitX::PDG_MASS_PSI2S - fitX::PDG_MASS_PSI2S_ERR, h->GetMinimum(), fitX::PDG_MASS_PSI2S + fitX::PDG_MASS_PSI2S_ERR, h->GetMaximum(), color_a, 0.8, 1001, 0, 0, 0);
      xjjroot::drawbox(fitX::PDG_MASS_X - fitX::PDG_MASS_X_ERR, h->GetMinimum(), fitX::PDG_MASS_X + fitX::PDG_MASS_X_ERR, h->GetMaximum(), color_b, 0.8, 1001, 0, 0, 0);
    }
  h->Fit("f","Nq");
  h->Fit("f","NLLq");
  // TFitResultPtr r = histo->Fit(func, "S");
  h->Fit("f","NLLq");
  roof->Draw("same");

  fitX::drawpull(h, f, color_data);
  if(h_ss) h_ss->Draw("pe same");

  TF1* fsig_a = xjjroot::copyobject(f, "fsig_a");
  fsig_a->SetRange(BIN_MIN_L, BIN_MAX_L);
  fsig_a->FixParameter(10, 0);
  fsig_a->FixParameter(0, 0);
  fsig_a->FixParameter(1, 0);
  fsig_a->FixParameter(2, 0);
  fsig_a->FixParameter(3, 0);
  fsig_a->FixParameter(4, 0);
  fsig_a->SetParError(5, f->GetParError(5));
  // xjjroot::settfstyle(fsig_a, color_a, 1, 3, color_a, 0.3, 1001);
  xjjroot::settfstyle(fsig_a, color_a, 7, 3, 0, 0, 1001);
  fsig_a->Draw("same");
  
  TF1* fsig_b = xjjroot::copyobject(f, "fsig_b");
  fsig_b->SetRange(BIN_MIN_H, BIN_MAX_H);
  fsig_b->FixParameter(5, 0);
  fsig_b->FixParameter(0, 0);
  fsig_b->FixParameter(1, 0);
  fsig_b->FixParameter(2, 0);
  fsig_b->FixParameter(3, 0);
  fsig_b->FixParameter(4, 0);
  fsig_b->SetParError(10, f->GetParError(10));
  // xjjroot::settfstyle(fsig_b, color_b, 1, 3, color_b, 0.3, 1001);
  xjjroot::settfstyle(fsig_b, color_b, 7, 3, 0, 0, 1001);
  fsig_b->Draw("same");

  TF1* fbkg = xjjroot::copyobject(f, "fbkg");
  fbkg->SetRange(BIN_MIN, BIN_MAX);
  fbkg->FixParameter(5, 0);
  fbkg->FixParameter(10, 0);
  xjjroot::settfstyle(fbkg, color_bkg, 7, 3);
  fbkg->Draw("same");

  f->Draw("same");
  xjjroot::copyobject(h, "hclone")->Draw("pesame");  
  c->RedrawAxis();
  float ysig_a = fsig_a->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;
  float ysig_b = fsig_b->Integral(BIN_MIN, BIN_MAX)/BIN_WIDTH;

  std::cout<<"\e[33;1m"<<ysig_a<<"  "<<ysig_b<<"\e[0m"<<std::endl;
  std::cout<<"\e[33;1m"<<pars[5]->getVal()<<"  "<<pars[10]->getVal()<<"\e[0m"<<std::endl;
  fitr->Print("v");

  xjjroot::drawtex(0.92, 0.84, Form("%.0f < p_{T} < %.0f GeV/c", fitX::ptmincut, fitX::ptmaxcut), 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.79, Form("%s|y| < %.1f", (fitX::ymincut?Form("%.1f < ", fitX::ymincut):""), fitX::ymaxcut), 0.042, 32, 62);
  xjjroot::drawtex(0.92, 0.74, Form("Cent. %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%"), 0.042, 32, 62);
  xjjroot::drawtex(0.17, 0.84, Form("#chi^{2} Prob = %.1f%s", TMath::Prob(f->GetChisquare(), f->GetNDF())*100., "%"), 0.042, 12, 62);

  xjjroot::drawtex(0.32, 0.36, Form("#bar{m}_{#psi(2S)} = %.4f GeV", f->GetParameter(6)), 0.034, 12, 62, color_a);
  xjjroot::drawtex(0.32, 0.31, Form("N_{#psi(2S)} = %.0f #pm %.0f", ysig_a, f->GetParError(5)*ysig_a/f->GetParameter(5)), 0.034, 12, 62, color_a);
  xjjroot::drawtex(0.68, 0.36, Form("#bar{m}_{X(3872)} = %.4f GeV", f->GetParameter(11)), 0.034, 12, 62, color_b);
  xjjroot::drawtex(0.68, 0.31, Form("N_{X(3872)} = %.0f #pm %.0f", ysig_b, f->GetParError(10)*ysig_b/f->GetParameter(10)), 0.034, 12, 62, color_b);
  xjjroot::drawCMS();

  if(saveplot) c->SaveAs(Form("%s/chmass%s.pdf", outputdir.c_str(), name.c_str()));

  cr->cd();
  dsh->plotOn(frempty, RooFit::Name("dsh"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  pdf->plotOn(frempty, RooFit::Name("sig_a"), RooFit::Components(*sig_a), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_a), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("sig_b"), RooFit::Components(*sig_b), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_b), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("bkg_poly"), RooFit::Components(bkg_poly), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(color_data), RooFit::LineWidth(3));
  pdf->plotOn(frempty, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(color_data), RooFit::LineWidth(3));
  dsh->plotOn(frempty, RooFit::Name("dsh"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(0.9), RooFit::MarkerStyle(20), RooFit::LineColor(1), RooFit::LineWidth(1));
  frempty->Draw();
  if(!fixmean)
    {
      xjjroot::drawbox(fitX::PDG_MASS_PSI2S - fitX::PDG_MASS_PSI2S_ERR, h->GetMinimum(), fitX::PDG_MASS_PSI2S + fitX::PDG_MASS_PSI2S_ERR, h->GetMaximum(), color_a, 0.8, 1001, 0, 0, 0);
      xjjroot::drawbox(fitX::PDG_MASS_X - fitX::PDG_MASS_X_ERR, h->GetMinimum(), fitX::PDG_MASS_X + fitX::PDG_MASS_X_ERR, h->GetMaximum(), color_b, 0.8, 1001, 0, 0, 0);
    }
  if(saveplot) cr->SaveAs(Form("%s/chmassr%s.pdf", outputdir.c_str(), name.c_str()));

  delete cmc;
  delete c;
  delete crmc;
  delete cr;

  funs.push_back(f);
  funs.push_back(fsig_a);
  funs.push_back(fsig_b);
  funs.push_back(fbkg);
  return funs;
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
  // h->GetXaxis()->SetNdivisions(510);
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
  xjjroot::drawline(binmin, hmc->GetMaximum()/2., binmax, hmc->GetMaximum()/2., kGray, 2, 2);
  xjjroot::drawaxis(binmax, 0, binmax, hmc->GetMaximum(), -4, 4, color, 1, 2, "+L");
  xjjroot::drawtex(0.93, 0.55, "Pull", 0.04, 33, 62, color);
}

#endif


// RooGaussian sig_a1("sig_a1", "", *mass, *(mcpars[6]), *(mcpars[7]));
// RooGaussian sig_a2("sig_a2", "", *mass, *(mcpars[6]), *(mcpars[8]));
// RooGaussian sig_mc_a1("sig_mc_a1", "", *mass, *(mcpars[6]), *(mcpars[7]));
// RooGaussian sig_mc_a2("sig_mc_a2", "", *mass, *(mcpars[6]), *(mcpars[8]));
// RooGaussian sig_b1("sig_b1", "", *mass, *(mcpars[11]), *(mcpars[12]));
// RooGaussian sig_b2("sig_b2", "", *mass, *(mcpars[11]), *(mcpars[13]));
// RooGaussian sig_mc_b1("sig_mc_b1", "", *mass, *(mcpars[11]), *(mcpars[12]));
// RooGaussian sig_mc_b2("sig_mc_b2", "", *mass, *(mcpars[11]), *(mcpars[13]));
