#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TPad.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fit.h"
#include "MCefficiency.h"
#include "systematics.h"

void drawkinematic();
void drawlabels();
void drawval(TEfficiency* greff_a, TEfficiency* greff_b, float ymax);
void fitX_fithist(std::string input, std::string output, std::string inputtnp_a, std::string inputtnp_b, std::string fitopt="")
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
  if(fitopt!="") { output += ("/"+fitopt); }
  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  RooDataSet* dsh = (RooDataSet*)ww->data("dsh");
  RooDataSet* dshBenr = (RooDataSet*)ww->data("dshBenr");
  RooDataSet* dshmcp_a = (RooDataSet*)ww->data("dshmcp_a");
  RooDataSet* dshmcp_b = (RooDataSet*)ww->data("dshmcp_b");
  TH1F* h = (TH1F*)inf->Get("h");
  TH1F* hBenr = (TH1F*)inf->Get("hBenr");
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  // hmcp_a->Scale(hmcp_a->GetEntries()/hmcp_a->Integral());
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  // hmcp_b->Scale(hmcp_b->GetEntries()/hmcp_b->Integral());
  TH1F* hlxymcnp_a = (TH1F*)inf->Get("hlxymcnp_a");
  TH1F* hlxymcnp_b = (TH1F*)inf->Get("hlxymcnp_b");
  TH1F* hlxymcp_a = (TH1F*)inf->Get("hlxymcp_a");
  TH1F* hlxymcp_b = (TH1F*)inf->Get("hlxymcp_b");
  MCeff::MCefficiency* mceff_a = new MCeff::MCefficiency(inf, "_a");
  MCeff::MCefficiency* mceff_b = new MCeff::MCefficiency(inf, "_b");

  // fit + yield
  TH1F* hyield_a = new TH1F("hyield_a", "", 5, 0, 5); //hyield_a->Sumw2();
  TH1F* hyield_b = new TH1F("hyield_b", "", 5, 0, 5); //hyield_b->Sumw2();
  TH1F* hBenryield_a = new TH1F("hBenryield_a", "", 5, 0, 5); //hBenryield_a->Sumw2();
  TH1F* hBenryield_b = new TH1F("hBenryield_b", "", 5, 0, 5); //hBenryield_b->Sumw2();
  std::vector<TH1F*> vh = {h, hBenr};
  std::vector<RooDataSet*> vdsh = {dsh, dshBenr};
  std::vector<TH1F*> vhyield_a = {hyield_a, hBenryield_a};
  std::vector<TH1F*> vhyield_b = {hyield_b, hBenryield_b};
  std::vector<std::string> vname = {"th", "thBenr"};
  std::vector<std::string> vtitle = {Form("#splitline{Inclusive}{%s}", fitopt.c_str()), Form("#splitline{B-enriched (l#scale[0.60]{xy} > 0.1 mm)}{%s}", fitopt.c_str())};
  std::vector<Style_t> mstyle = {20, 24};

  float mm = 0;
  xjjroot::setgstyle();
  std::vector<RooAbsPdf*> pdf(2, 0), bkg(2, 0);
  std::vector<TF1*> frf(2, 0);
  std::vector<float> maxy(2, 0);
  TCanvas* cy = new TCanvas("cy", "", 1200, 600);
  cy->Divide(2, 1);
  for(int l=0; l<2; l++)
    {
      // ====>
      std::map<std::string, fitX::fitXresult*> result = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, 
                                                                  vdsh[l], dshmcp_a, dshmcp_b,
                                                                  Form("plots/%s", output.c_str()), mm, true, "_"+vname[l], vtitle[l], fitopt);
      cy->cd(l+1);
      xjjroot::setgstyle();
      // <====
      float ysig_a = result["unbinned"]->ysig_a();
      float ysigerr_a = result["unbinned"]->ysigerr_a();
      float ysig_b = result["unbinned"]->ysig_b();
      float ysigerr_b = result["unbinned"]->ysigerr_b();
      mm = result["unbinned"]->msig_b();
      RooWorkspace* wwr = result["unbinned"]->ww();
      pdf[l] = wwr->pdf("pdf");
      bkg[l] = wwr->pdf("bkg");
      TString str_bkg = "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x";
      TString str_sig_a = "[5]*( [9]*TMath::Gaus(x, [6], [7])/(TMath::Sqrt(2*3.14159)* [7]) + (1- [9])*( [16]*TMath::Gaus(x,[6], [8])/(TMath::Sqrt(2*3.14159)* [8]) + (1- [16])*TMath::Gaus(x,[6], [15])/(TMath::Sqrt(2*3.14159)*[15])))";
      TString str_sig_b = "[10]*([14]*TMath::Gaus(x,[11],[12])/(TMath::Sqrt(2*3.14159)*[12]) + (1-[14])*([18]*TMath::Gaus(x,[11],[13])/(TMath::Sqrt(2*3.14159)*[13]) + (1-[18])*TMath::Gaus(x,[11],[17])/(TMath::Sqrt(2*3.14159)*[17])))";
      TF1* f = new TF1("f", str_bkg+"+"+str_sig_a+"+"+str_sig_b, fitX::BIN_MIN, fitX::BIN_MAX);
      frf[l] = fitX::astf(pdf[l], f, "fr");
      maxy[l] = result["unbinned"]->get("maxy");

      // yield
      xjjroot::setthgrstyle(vhyield_a[l], fitX::color_a, mstyle[l], 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
      vhyield_a[l]->SetBinContent(fitX::ibin_a, ysig_a);
      vhyield_a[l]->SetBinError(fitX::ibin_a, ysigerr_a);
      xjjroot::setthgrstyle(vhyield_b[l], fitX::color_b, mstyle[l], 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
      vhyield_b[l]->SetBinContent(fitX::ibin_b, ysig_b);
      vhyield_b[l]->SetBinError(fitX::ibin_b, ysigerr_b);
      // float ymax = 300.; // !
      float ymax = vhyield_a[l]->GetMaximum()*2.5; // !
      TH2F* hempty = new TH2F(Form("hempty%d",l), ";;Raw Yield", 5, 0, 5, 10, 0, ymax);
      xjjroot::sethempty(hempty, 0, 0.3);
      hempty->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
      hempty->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
      hempty->GetXaxis()->SetLabelSize(hempty->GetXaxis()->GetLabelSize()*1.5);
      hempty->Draw();
      vhyield_a[l]->Draw("ple same");
      vhyield_b[l]->Draw("ple same");
      xjjroot::drawtex(0.42, ysig_a/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", ysig_a, ysigerr_a), 0.042, 22, 62, fitX::color_a);
      xjjroot::drawtex(0.72, ysig_b/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", ysig_b, ysigerr_b), 0.042, 22, 62, fitX::color_b);
      drawkinematic();
      xjjroot::drawtex(0.22, 0.84, vtitle[l].c_str(), 0.042, 12, 62);
      xjjroot::drawCMS();
    }
  xjjroot::mkdir(Form("plots/%s/cyield.pdf", output.c_str()));
  cy->SaveAs(Form("plots/%s/cyield.pdf", output.c_str()));

  // fprompt
  hlxymcnp_a->GetYaxis()->SetTitle("Events");
  hlxymcnp_b->GetYaxis()->SetTitle("Events");
  hlxymcp_a->GetYaxis()->SetTitle("Events");
  hlxymcp_b->GetYaxis()->SetTitle("Events");
  hlxymcnp_a->Scale(1./hlxymcnp_a->Integral(), "width");
  hlxymcnp_b->Scale(1./hlxymcnp_b->Integral(), "width");
  hlxymcp_a->Scale(1./hlxymcp_a->Integral(), "width");
  hlxymcp_b->Scale(1./hlxymcp_b->Integral(), "width");
  xjjroot::sethempty(hlxymcnp_a, 0, 0);
  xjjroot::sethempty(hlxymcp_a, 0, 0);
  xjjroot::setthgrstyle(hlxymcp_a, fitX::color_a, 21, 1, fitX::color_a, 1, 2);
  xjjroot::setthgrstyle(hlxymcp_b, fitX::color_b, 21, 1, fitX::color_b, 1, 2);
  xjjroot::setthgrstyle(hlxymcnp_a, fitX::color_a, 21, 1, fitX::color_a, 1, 2);
  xjjroot::setthgrstyle(hlxymcnp_b, fitX::color_b, 21, 1, fitX::color_b, 1, 2);
  hlxymcp_a->SetMaximum(hlxymcp_a->GetMaximum()*10.);
  hlxymcnp_a->SetMaximum(hlxymcnp_a->GetMaximum()*10.);
  TCanvas* clxy = new TCanvas("clxy", "", 1200, 600);
  clxy->Divide(2, 1);
  clxy->cd(1);
  gPad->SetLogy();
  hlxymcp_a->Draw("histe");
  hlxymcp_b->Draw("histe same");
  // xjjroot::drawline(0.1, 0, 0.1, hlxymcp_a->GetMaximum(), kGray+1, 2, 2);
  xjjroot::drawline(0.1, 0, 0.1, 2., kGray+1, 2, 2);
  drawkinematic();
  xjjroot::drawtex(0.23, 0.85, fitX::title_a.c_str(), 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, fitX::title_b.c_str(), 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.85, "Prompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  clxy->cd(2);
  gPad->SetLogy();
  hlxymcnp_a->Draw("histe");
  hlxymcnp_b->Draw("histe same");
  // xjjroot::drawline(0.1, 0, 0.1, hlxymcnp_a->GetMaximum(), kGray+1, 2, 2);
  xjjroot::drawline(0.1, 0, 0.1, 2., kGray+1, 2, 2);
  drawkinematic();
  xjjroot::drawtex(0.23, 0.85, fitX::title_a.c_str(), 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, fitX::title_b.c_str(), 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.85, "Nonprompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  xjjroot::mkdir(Form("plots/%s/clxy.pdf", output.c_str()));
  clxy->SaveAs(Form("plots/%s/clxy.pdf", output.c_str()));
  
  std::vector<double> vlxyfrac = lxydis::nplxyfrac(hlxymcnp_a, hlxymcnp_b);
  TH1F* hlxyfrac_a = new TH1F("hlxyfrac_a", "", 5, 0, 5); //hlxyfrac_a->Sumw2();
  hlxyfrac_a->SetBinContent(fitX::ibin_a, vlxyfrac[0]);
  hlxyfrac_a->SetBinError(fitX::ibin_a, vlxyfrac[1]);
  TH1F* hlxyfrac_b = new TH1F("hlxyfrac_b", "", 5, 0, 5); //hlxyfrac_b->Sumw2();
  hlxyfrac_b->SetBinContent(fitX::ibin_b, vlxyfrac[2]);
  hlxyfrac_b->SetBinError(fitX::ibin_b, vlxyfrac[3]);
  TH1F *hyieldprompt_a, *hyieldprompt_b;
  TEfficiency* grfprompt_a = lxydis::calclxyfprompt(hyield_a, hBenryield_a, hlxyfrac_a, "grfprompt_a", &hyieldprompt_a);
  TEfficiency* grfprompt_b = lxydis::calclxyfprompt(hyield_b, hBenryield_b, hlxyfrac_b, "grfprompt_b", &hyieldprompt_b);

  xjjroot::setthgrstyle(grfprompt_a, fitX::color_a, 25, 1.5, fitX::color_a, 2, 3);
  xjjroot::setthgrstyle(grfprompt_b, fitX::color_b, 25, 1.5, fitX::color_b, 2, 3);
  xjjroot::setgstyle(2);
  TCanvas* cfprompt = new TCanvas("cfprompt", "", 600, 600);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", 5, 0, 5, 10, 0, 1.4);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyfprompt->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.5);
  hemptyfprompt->Draw();
  xjjroot::drawline(0, 1, 5, 1, kGray+1, 9, 2);
  grfprompt_a->Draw("same");
  grfprompt_b->Draw("same");
  drawkinematic();
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::drawCMS();
  xjjroot::mkdir(Form("plots/%s/cfprompt.pdf", output.c_str()));
  cfprompt->SaveAs(Form("plots/%s/cfprompt.pdf", output.c_str()));
  std::cout<<"f_prompt(psi') = "<<grfprompt_a->GetEfficiency(fitX::ibin_a)<<" + "<<grfprompt_a->GetEfficiencyErrorUp(fitX::ibin_a)<<" - "<<grfprompt_a->GetEfficiencyErrorLow(fitX::ibin_a)<<std::endl;
  std::cout<<"f_prompt(X) = "<<grfprompt_b->GetEfficiency(fitX::ibin_b)<<" + "<<grfprompt_b->GetEfficiencyErrorUp(fitX::ibin_b)<<" - "<<grfprompt_b->GetEfficiencyErrorLow(fitX::ibin_b)<<std::endl;

  // efficiency
  mceff_a->calceff();
  mceff_a->calcacc();
  mceff_a->setstyle(fitX::color_a);
  mceff_b->calceff();
  mceff_b->calcacc();
  mceff_b->setstyle(fitX::color_b);
  float ymaxeff = 0.2, ymaxacc = 1.0, ymaxeffpre = 1.0, ymaxeffcut = 1.0, ymaxeffbdt = 1.3, ymaxeffqvl = 1.3;
  TH2F* hemptyeff = MCeff::createhempty("hemptyeff", "#alpha #times #epsilon_{reco} #times #epsilon_{sel}", ymaxeff);
  TH2F* hemptyeff_incl = MCeff::createhempty_incl("hemptyeff_incl", "#alpha #times #epsilon_{reco} #times #epsilon_{sel}", ymaxeff);
  TH2F* hemptyacc = MCeff::createhempty("hemptyacc", "#alpha", ymaxacc);
  TH2F* hemptyacc_incl = MCeff::createhempty_incl("hemptyacc_incl", "#alpha", ymaxacc);
  TH2F* hemptyeffpre = MCeff::createhempty("hemptyeffpre", "#epsilon_{reco}", ymaxeffpre);
  TH2F* hemptyeffpre_incl = MCeff::createhempty_incl("hemptyeffpre_incl", "#epsilon_{reco}", ymaxeffpre);
  TH2F* hemptyeffcut = MCeff::createhempty("hemptyeffcut", "#epsilon_{sel}", ymaxeffcut);
  TH2F* hemptyeffcut_incl = MCeff::createhempty_incl("hemptyeffcut_incl", "#epsilon_{sel}", ymaxeffcut);
  TH2F* hemptyeffbdt = MCeff::createhempty("hemptyeffbdt", "#epsilon_{sel} (BDT)", ymaxeffbdt);
  TH2F* hemptyeffbdt_incl = MCeff::createhempty_incl("hemptyeffbdt_incl", "#epsilon_{sel} (BDT)", ymaxeffbdt);
  TH2F* hemptyeffqvl = MCeff::createhempty("hemptyeffqvl", "#epsilon_{sel} (Q value)", ymaxeffqvl);
  TH2F* hemptyeffqvl_incl = MCeff::createhempty_incl("hemptyeffqvl_incl", "#epsilon_{sel} (Q value)", ymaxeffqvl);
  TLegend* legeff = new TLegend(0.70, 0.20, 1.10, 0.32);
  xjjroot::setleg(legeff, 0.042);
  legeff->AddEntry(mceff_a->greff(), fitX::title_a.c_str(), "fl");
  legeff->AddEntry(mceff_b->greff(), fitX::title_b.c_str(), "fl");
  TLegend* legacc = new TLegend(0.70, 0.40, 1.10, 0.52);
  xjjroot::setleg(legacc, 0.042);
  legacc->AddEntry(mceff_a->gracc(), fitX::title_a.c_str(), "fl");
  legacc->AddEntry(mceff_b->gracc(), fitX::title_b.c_str(), "fl");

  xjjroot::setgstyle();

  TCanvas* ceff = new TCanvas("ceff", "", 1200, 600);
  ceff->Divide(2, 1);
  ceff->cd(1);
  hemptyeff->Draw();
  mceff_a->greff()->Draw("same3");
  mceff_a->greff()->Draw("samelX");
  mceff_b->greff()->Draw("same3");
  mceff_b->greff()->Draw("samelX");
  legeff->Draw();
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  xjjroot::drawcomment(output.c_str(), "r");
  ceff->cd(2);
  hemptyeff_incl->Draw();
  mceff_a->greff_incl()->Draw("same ple");
  mceff_b->greff_incl()->Draw("same ple");
  float effval_a = mceff_a->greff_incl()->GetEfficiency(fitX::ibin_a);
  xjjroot::drawtex(0.42, effval_a/ymaxeff*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, 
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", mceff_a->greff_incl()->GetEfficiency(fitX::ibin_a)*1000, mceff_a->greff_incl()->GetEfficiencyErrorUp(fitX::ibin_a)*1000, mceff_a->greff_incl()->GetEfficiencyErrorLow(fitX::ibin_a)*1000), 0.042, 22, 62, fitX::color_a);
  float effval_b = mceff_b->greff_incl()->GetEfficiency(fitX::ibin_b);
  xjjroot::drawtex(0.72, effval_b/ymaxeff*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1,
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", mceff_b->greff_incl()->GetEfficiency(fitX::ibin_b)*1000, mceff_b->greff_incl()->GetEfficiencyErrorUp(fitX::ibin_b)*1000, mceff_b->greff_incl()->GetEfficiencyErrorLow(fitX::ibin_b)*1000), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  drawkinematic();
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  xjjroot::mkdir(Form("plots/%s/ceff.pdf", output.c_str()));
  ceff->SaveAs(Form("plots/%s/ceff.pdf", output.c_str()));

  TCanvas* cacc = new TCanvas("cacc", "", 1800, 1200);
  cacc->Divide(3, 2);
  cacc->cd(1);
  hemptyacc->Draw();
  mceff_a->gracc()->Draw("same3");
  mceff_a->gracc()->Draw("samelX");
  mceff_b->gracc()->Draw("same3");
  mceff_b->gracc()->Draw("samelX");
  legacc->Draw();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  drawlabels();
  cacc->cd(2);
  hemptyeffpre->Draw();
  mceff_a->greffpre()->Draw("same3");
  mceff_a->greffpre()->Draw("samelX");
  mceff_b->greffpre()->Draw("same3");
  mceff_b->greffpre()->Draw("samelX");
  legacc->Draw();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  drawlabels();
  cacc->cd(3);
  hemptyeffcut->Draw();
  mceff_a->greffcut()->Draw("same3");
  mceff_a->greffcut()->Draw("samelX");
  mceff_b->greffcut()->Draw("same3");
  mceff_b->greffcut()->Draw("samelX");
  legacc->Draw();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  drawlabels();
  cacc->cd(4);
  hemptyacc_incl->Draw();
  mceff_a->gracc_incl()->Draw("same ple");
  mceff_b->gracc_incl()->Draw("same ple");
  drawval(mceff_a->gracc_incl(), mceff_b->gracc_incl(), hemptyacc_incl->GetYaxis()->GetXmax());
  drawkinematic();
  drawlabels();
  cacc->cd(5);
  hemptyeffpre_incl->Draw();
  mceff_a->greffpre_incl()->Draw("same ple");
  mceff_b->greffpre_incl()->Draw("same ple");
  drawval(mceff_a->greffpre_incl(), mceff_b->greffpre_incl(), hemptyeffpre_incl->GetYaxis()->GetXmax());
  drawkinematic();
  drawlabels();
  cacc->cd(6);
  hemptyeffcut_incl->Draw();
  mceff_a->greffcut_incl()->Draw("same ple");
  mceff_b->greffcut_incl()->Draw("same ple");
  drawval(mceff_a->greffcut_incl(), mceff_b->greffcut_incl(), hemptyeffcut_incl->GetYaxis()->GetXmax());
  drawkinematic();
  xjjroot::drawcomment(output.c_str(), "r");
  drawlabels();
  xjjroot::mkdir(Form("plots/%s/cacc.pdf", output.c_str()));
  cacc->SaveAs(Form("plots/%s/cacc.pdf", output.c_str()));

  TCanvas* ceffsel = new TCanvas("ceffsel", "", 1200, 1200);
  ceffsel->Divide(2, 2);
  ceffsel->cd(1);
  hemptyeffbdt->Draw();
  mceff_a->greffbdt()->Draw("same3");
  mceff_a->greffbdt()->Draw("samelX");
  mceff_b->greffbdt()->Draw("same3");
  mceff_b->greffbdt()->Draw("samelX");
  legacc->Draw();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  drawlabels();
  ceffsel->cd(2);
  hemptyeffqvl->Draw();
  mceff_a->greffqvl()->Draw("same3");
  mceff_a->greffqvl()->Draw("samelX");
  mceff_b->greffqvl()->Draw("same3");
  mceff_b->greffqvl()->Draw("samelX");
  legacc->Draw();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  drawlabels();
  ceffsel->cd(3);
  hemptyeffbdt_incl->Draw();
  mceff_a->greffbdt_incl()->Draw("same ple");
  mceff_b->greffbdt_incl()->Draw("same ple");
  drawval(mceff_a->greffbdt_incl(), mceff_b->greffbdt_incl(), hemptyeffbdt_incl->GetYaxis()->GetXmax());
  drawkinematic();
  drawlabels();
  ceffsel->cd(4);
  hemptyeffqvl_incl->Draw();
  mceff_a->greffqvl_incl()->Draw("same ple");
  mceff_b->greffqvl_incl()->Draw("same ple");
  drawval(mceff_a->greffqvl_incl(), mceff_b->greffqvl_incl(), hemptyeffqvl_incl->GetYaxis()->GetXmax());
  drawkinematic();
  drawlabels();
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::mkdir(Form("plots/%s/ceffsel.pdf", output.c_str()));
  ceffsel->SaveAs(Form("plots/%s/ceffsel.pdf", output.c_str()));

  // tnp
  TFile* inftnp_a = TFile::Open(inputtnp_a.c_str());
  TH1D* hscale_htnp_total_nominal_a = (TH1D*)inftnp_a->Get("scale_htnp_total_nominal_1");
  hscale_htnp_total_nominal_a->SetName("hscale_htnp_total_nominal_a");
  hyieldprompt_a->Scale(hscale_htnp_total_nominal_a->GetBinContent(1));
  TFile* inftnp_b = TFile::Open(inputtnp_b.c_str());
  TH1D* hscale_htnp_total_nominal_b = (TH1D*)inftnp_b->Get("scale_htnp_total_nominal_1");
  hscale_htnp_total_nominal_b->SetName("hscale_htnp_total_nominal_b");
  hyieldprompt_b->Scale(hscale_htnp_total_nominal_b->GetBinContent(1));

  // correct eff
  xjjroot::setgstyle(2);
  TH1F* hyieldpromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldpromptCorr_a");
  hyieldpromptCorr_a->Divide(mceff_a->heff_incl());
  TH1F* hyieldpromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldpromptCorr_b");
  hyieldpromptCorr_b->Divide(mceff_b->heff_incl());
  std::vector<float> xx_a, yy_a, xel_a, xeh_a, yel_a, yeh_a;
  xx_a.push_back(hyieldpromptCorr_a->GetBinCenter(fitX::ibin_a));
  xel_a.push_back(0.1);
  xeh_a.push_back(0.1);
  yy_a.push_back(hyieldpromptCorr_a->GetBinContent(fitX::ibin_a));
  yel_a.push_back(syst::getsyst(0, "d")*hyieldpromptCorr_a->GetBinContent(fitX::ibin_a));
  yeh_a.push_back(syst::getsyst(0, "u")*hyieldpromptCorr_a->GetBinContent(fitX::ibin_a));
  TGraphAsymmErrors* gsyst_a = new TGraphAsymmErrors(1, xx_a.data(), yy_a.data(), xel_a.data(), xeh_a.data(), yel_a.data(), yeh_a.data());
  gsyst_a->SetName("gr_a_syst");
  xjjroot::setthgrstyle(gsyst_a, fitX::color_a, 20, 1.2, 0, 0, 0, fitX::color_a, 0.3, 1001);
  std::vector<float> xx_b, yy_b, xel_b, xeh_b, yel_b, yeh_b;
  xx_b.push_back(hyieldpromptCorr_b->GetBinCenter(fitX::ibin_b));
  xel_b.push_back(0.1);
  xeh_b.push_back(0.1);
  yy_b.push_back(hyieldpromptCorr_b->GetBinContent(fitX::ibin_b));
  yel_b.push_back(syst::getsyst(1, "d")*hyieldpromptCorr_b->GetBinContent(fitX::ibin_b));
  yeh_b.push_back(syst::getsyst(1, "u")*hyieldpromptCorr_b->GetBinContent(fitX::ibin_b));
  TGraphAsymmErrors* gsyst_b = new TGraphAsymmErrors(1, xx_b.data(), yy_b.data(), xel_b.data(), xeh_b.data(), yel_b.data(), yeh_b.data());
  gsyst_b->SetName("gr_b_syst");
  xjjroot::setthgrstyle(gsyst_b, fitX::color_b, 20, 1.2, 0, 0, 0, fitX::color_b, 0.3, 1001);
  TLegend* leg = new TLegend(0.24, 0.85-0.05*2, 0.65, 0.85);
  xjjroot::setleg(leg, 0.042);
  leg->SetNColumns(2);
  leg->AddEntry(hyieldpromptCorr_a, " ", "pe");
  leg->AddEntry(hyieldpromptCorr_b, "Stat.", "pe");
  leg->AddEntry(gsyst_a, " ", "pf");
  leg->AddEntry(gsyst_b, "Syst.", "pf");
  float ymaxyieldpromptCorr = std::max(hyieldpromptCorr_a->GetMaximum(), hyieldpromptCorr_b->GetMaximum())*2.5; // !
  TH2F* hemptyyieldpromptCorr = new TH2F("hemptyyieldpromptCorr", ";;N_{signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", 5, 0, 5, 10, 0, ymaxyieldpromptCorr);
  xjjroot::sethempty(hemptyyieldpromptCorr, 0, 0.3);
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(fitX::ibin_a, fitX::title_a.c_str());
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(fitX::ibin_b, fitX::title_b.c_str());
  hemptyyieldpromptCorr->GetXaxis()->SetLabelSize(hemptyyieldpromptCorr->GetXaxis()->GetLabelSize()*1.5);
  TCanvas* cyieldpromptCorr = new TCanvas("cyieldpromptCorr", "", 600, 600);
  hemptyyieldpromptCorr->Draw();
  gsyst_a->Draw("2 same");
  gsyst_b->Draw("2 same");
  hyieldpromptCorr_a->Draw("ple same");
  hyieldpromptCorr_b->Draw("ple same");
  leg->Draw();
  float yyieldpromptCorr_a = hyieldpromptCorr_a->GetBinContent(fitX::ibin_a);
  float yerryieldpromptCorr_a = hyieldpromptCorr_a->GetBinError(fitX::ibin_a);
  xjjroot::drawtex(0.42, yyieldpromptCorr_a/ymaxyieldpromptCorr*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.2, Form("%.0f #pm %.0f", yyieldpromptCorr_a, yerryieldpromptCorr_a), 0.042, 22, 62, fitX::color_a);
  float yyieldpromptCorr_b = hyieldpromptCorr_b->GetBinContent(fitX::ibin_b);
  float yerryieldpromptCorr_b = hyieldpromptCorr_b->GetBinError(fitX::ibin_b);
  xjjroot::drawtex(0.72, yyieldpromptCorr_b/ymaxyieldpromptCorr*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.2, Form("%.0f #pm %.0f", yyieldpromptCorr_b, yerryieldpromptCorr_b), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawCMS();
  drawkinematic();
  xjjroot::drawcomment(output.c_str(), "r");
  xjjroot::mkdir(Form("plots/%s/cyieldpromptCorr.pdf", output.c_str()));
  cyieldpromptCorr->SaveAs(Form("plots/%s/cyieldpromptCorr.pdf", output.c_str()));

  std::cout<<"-->"<<std::endl;
  std::cout<<"_a: "<<xx_a[0]<<" "<<yy_a[0]<<" "<<hyieldpromptCorr_a->GetBinContent(fitX::ibin_a)<<" "<<hyieldpromptCorr_a->GetBinError(fitX::ibin_a)<<" "<<yel_a[0]<<" "<<yeh_a[0]<<std::endl;
  std::cout<<"_b: "<<xx_b[0]<<" "<<yy_b[0]<<" "<<hyieldpromptCorr_b->GetBinContent(fitX::ibin_b)<<" "<<hyieldpromptCorr_b->GetBinError(fitX::ibin_b)<<" "<<yel_b[0]<<" "<<yeh_b[0]<<std::endl;
  std::cout<<"<--"<<std::endl;

  // merge a + b
  TH1F* hyield = (TH1F*)hyield_a->Clone("hyield");
  hyield->Add(hyield_b);
  TH1F* hBenryield = (TH1F*)hBenryield_a->Clone("hBenryield");
  hBenryield->Add(hBenryield_b);
  TH1F* heff_incl = (TH1F*)mceff_a->heff_incl()->Clone("heff_incl");
  heff_incl->Add(mceff_b->heff_incl());
  TH1F* hyieldprompt = (TH1F*)hyieldprompt_a->Clone("hyieldprompt");
  hyieldprompt->Add(hyieldprompt_b);
  TH1F* hyieldpromptCorr = (TH1F*)hyieldpromptCorr_a->Clone("hyieldpromptCorr");
  hyieldpromptCorr->Add(hyieldpromptCorr_b);

  TH1F* hratio_a = new TH1F("hratio_a", ";;", 1, fitX::ptmincut, fitX::ptmaxcut);
  hratio_a->SetBinContent(1, hyieldpromptCorr_a->GetBinContent(fitX::ibin_a));
  hratio_a->SetBinError(1, hyieldpromptCorr_a->GetBinError(fitX::ibin_a));
  TH1F* hratio_b = new TH1F("hratio_b", ";;", 1, fitX::ptmincut, fitX::ptmaxcut);
  hratio_b->SetBinContent(1, hyieldpromptCorr_b->GetBinContent(fitX::ibin_b));
  hratio_b->SetBinError(1, hyieldpromptCorr_b->GetBinError(fitX::ibin_b));
  TH1F* hratio = (TH1F*)hratio_b->Clone("hratio");
  hratio->Divide(hratio_a);

  // merge inclusive + B-enr
  TF1 *f_dump_incl = new TF1(), *bkg_dump_incl = new TF1(), *f_dump_Benr = new TF1(), *bkg_dump_Benr = new TF1();
  TH1F* h_dump = new TH1F("h_dump", "", 1, 0, 1);
  xjjroot::setthgrstyle(h_dump, kBlack, 20, 0.9, kBlack, 1, 1);
  xjjroot::settfstyle(f_dump_incl, fitX::color_data, 1, 2);
  xjjroot::settfstyle(f_dump_Benr, fitX::color_data2, 1, 2);
  xjjroot::settfstyle(bkg_dump_incl, fitX::color_data, 7, 2);
  xjjroot::settfstyle(bkg_dump_Benr, fitX::color_data2, 7, 2);
  h->SetBinErrorOption(TH1::kPoisson);
  hBenr->SetBinErrorOption(TH1::kPoisson);
  // TH1F* h_create = (TH1F*)vdsh[0]->createHistogram("h_create", );
  xjjroot::setgstyle(2);
  gStyle->SetPadLeftMargin(gStyle->GetPadLeftMargin()*0.7);
  gStyle->SetPadRightMargin(gStyle->GetPadRightMargin()*1.4);  
  TCanvas* cr = new TCanvas("cr", "", 700, 600);
  RooRealVar* mass = new RooRealVar("Bmass", "Bmass", fitX::BIN_MIN, fitX::BIN_MAX);
  RooPlot* frempty = mass->frame(RooFit::Title(""));
  frempty->SetXTitle("m_{#mu#mu#pi#pi} (GeV/c^{2})");
  frempty->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  fitX::setmasshist(frempty, 0, -0.2);
  // frempty->SetTitle("");
  vdsh[0]->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()));
  pdf[0]->plotOn(frempty, RooFit::Name("bkg"), RooFit::Components(*(bkg[0])), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(fitX::color_data), RooFit::LineWidth(3));
  pdf[0]->plotOn(frempty, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(fitX::color_data), RooFit::LineWidth(3));
  vdsh[0]->plotOn(frempty, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()));
  vdsh[1]->plotOn(frempty, RooFit::Name("dshist_Benr"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()));
  pdf[1]->plotOn(frempty, RooFit::Name("bkg_Benr"), RooFit::Components(*(bkg[1])), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(7), RooFit::LineColor(fitX::color_data2), RooFit::LineWidth(3));
  pdf[1]->plotOn(frempty, RooFit::Name("pdf_Benr"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(1), RooFit::LineColor(fitX::color_data2), RooFit::LineWidth(3));
  vdsh[1]->plotOn(frempty, RooFit::Name("dshist_Benr"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()));
  frempty->SetMinimum(0); h->SetMinimum(0);
  frempty->SetMaximum(maxy[0]); h->SetMaximum(maxy[0]);
  frempty->Draw();
  // fitX::drawpull(h, frf[0], fitX::color_data);
  xjjroot::drawtex(0.24, 0.41, fitX::title_a.c_str(), 0.038);
  xjjroot::drawtex(0.61, 0.52, fitX::title_b.c_str(), 0.038);
  xjjroot::drawtex(0.93, 0.46, "Inclusive", 0.035, 32, 42, fitX::color_data);
  xjjroot::drawtex(0.93, 0.23, "b-enriched", 0.035, 32, 42, fitX::color_data2);
  fitX::labelsdata("", "#scale[1.25]{#bf{CMS}}", "1.7 nb^{-1} (PbPb 5.02 TeV)");
  cr->SaveAs(Form("plots/%s/chmassr_both.pdf", output.c_str()));
  delete cr;

  xjjroot::setgstyle(2);
  float yincl = 0.5, ybenr = 1-yincl;
  RooPlot* frempty_incl = mass->frame(RooFit::Title(""));
  frempty_incl->SetXTitle("m_{#mu#mu#pi#pi} (GeV/c^{2})");
  frempty_incl->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  fitX::setmasshist(frempty_incl, 0, -0.6);
  frempty_incl->GetXaxis()->SetTitleSize(frempty_incl->GetXaxis()->GetTitleSize()*1.2);
  frempty_incl->GetYaxis()->SetTitleSize(frempty_incl->GetYaxis()->GetTitleSize()*1.2);
  frempty_incl->GetXaxis()->SetLabelSize(frempty_incl->GetXaxis()->GetLabelSize()*1.2);
  frempty_incl->GetYaxis()->SetLabelSize(frempty_incl->GetYaxis()->GetLabelSize()*1.2); h->GetYaxis()->SetLabelSize(frempty_incl->GetYaxis()->GetLabelSize());
  RooPlot* frempty_benr = mass->frame(RooFit::Title(""));
  frempty_benr->SetXTitle("m_{#mu#mu#pi#pi} (GeV/c^{2})");
  frempty_benr->SetYTitle(Form("Entries / (%.0f MeV/c^{2})", fitX::BIN_WIDTH*1.e+3));
  fitX::setmasshist(frempty_benr, 0, -0.6);
  frempty_benr->GetXaxis()->SetTitleSize(frempty_benr->GetXaxis()->GetTitleSize()*1.2);
  frempty_benr->GetYaxis()->SetTitleSize(frempty_benr->GetYaxis()->GetTitleSize()*1.2);
  frempty_benr->GetXaxis()->SetLabelSize(frempty_benr->GetXaxis()->GetLabelSize()*1.2);
  frempty_benr->GetYaxis()->SetLabelSize(frempty_benr->GetYaxis()->GetLabelSize()*1.2); hBenr->GetYaxis()->SetLabelSize(frempty_benr->GetYaxis()->GetLabelSize());
  TLegend* leg_both = new TLegend(0.70, 0.95-0.06*3, 0.95, 0.95);
  xjjroot::setleg(leg_both, 0.055);
  leg_both->SetNColumns(2);
  leg_both->AddEntry((TObject*)0, "", NULL);
  leg_both->AddEntry(h_dump, "data", "pe");
  leg_both->AddEntry(f_dump_incl, "", "l");
  leg_both->AddEntry(f_dump_Benr, "total fit", "l");
  leg_both->AddEntry(bkg_dump_incl, "", "l");
  leg_both->AddEntry(bkg_dump_Benr, "background", "l");
  cr = new TCanvas("cr", "", 700, 700);
  TPad* pincl = new TPad("pincl", "", 0, 1-yincl, 1, 1);
  pincl->SetMargin(xjjroot::margin_pad_left*0.5, xjjroot::margin_pad_right*1.5, 0, xjjroot::margin_pad_top*1.3);
  pincl->Draw();
  pincl->cd();
  vdsh[0]->plotOn(frempty_incl, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  pdf[0]->plotOn(frempty_incl, RooFit::Name("bkg"), RooFit::Components(*(bkg[0])), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(bkg_dump_incl->GetLineStyle()), RooFit::LineColor(bkg_dump_incl->GetLineColor()), RooFit::LineWidth(bkg_dump_incl->GetLineWidth()));
  pdf[0]->plotOn(frempty_incl, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(f_dump_incl->GetLineStyle()), RooFit::LineColor(f_dump_incl->GetLineColor()), RooFit::LineWidth(f_dump_incl->GetLineWidth()));
  vdsh[0]->plotOn(frempty_incl, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  frempty_incl->SetMinimum(130);
  frempty_incl->SetMaximum(440);
  frempty_incl->Draw();
  TH1F* hforpull_incl = fitX::createhistforpull(frempty_incl, "dshist", frf[0]);
  fitX::drawpull(hforpull_incl, frf[0], fitX::color_data);
  xjjroot::drawbox(4.02, h->GetMinimum(), 5, h->GetMaximum(), kWhite, 1);
  xjjroot::drawtex(0.99, 0.39, "Pull", frempty_incl->GetYaxis()->GetTitleSize(), 33, 42, fitX::color_data, 270);
  xjjroot::drawtex(0.89, 0.79, Form("%.0f < p_{T} < %.0f GeV/c", fitX::ptmincut, fitX::ptmaxcut), 0.05, 32, 62);
  xjjroot::drawtex(0.89, 0.79-0.06, Form("%s|y| < %.1f", (fitX::ymincut?Form("%.1f < ", fitX::ymincut):""), fitX::ymaxcut), 0.05, 32, 62);
  xjjroot::drawtex(0.89, 0.79-0.06*2, Form("Cent. %.0f-%.0f%s", fitX::centmincut, fitX::centmaxcut, "%"), 0.05, 32, 62);
  xjjroot::drawtex(0.21, 0.12, fitX::title_a.c_str(), 0.055);
  xjjroot::drawtex(0.60, 0.31, fitX::title_b.c_str(), 0.055);
  xjjroot::drawtex(0.14, 0.79, "#scale[1.25]{#bf{CMS}}", 0.06, 12);
  xjjroot::drawtex(0.93, 0.92, "1.7 nb^{-1} (PbPb 5.02 TeV)", 0.055, 32);
  xjjroot::drawtex(0.14, 0.70, "Inclusive", 0.06, 12, 52, kBlack);
  cr->cd();
  TPad* pbenr = new TPad("pbenr", "", 0, 0, 1, ybenr);
  pbenr->SetMargin(xjjroot::margin_pad_left*0.5, xjjroot::margin_pad_right*1.5, xjjroot::margin_pad_bottom, 0);
  pbenr->Draw();
  pbenr->cd();
  vdsh[1]->plotOn(frempty_benr, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  pdf[1]->plotOn(frempty_benr, RooFit::Name("bkg"), RooFit::Components(*(bkg[1])), RooFit::Precision(1e-6), RooFit::DrawOption("L"), RooFit::LineStyle(bkg_dump_Benr->GetLineStyle()), RooFit::LineColor(bkg_dump_Benr->GetLineColor()), RooFit::LineWidth(bkg_dump_Benr->GetLineWidth()));
  pdf[1]->plotOn(frempty_benr, RooFit::Name("pdf"), RooFit::Precision(1e-6), RooFit::Normalization(1.0, RooAbsReal::RelativeExpected), RooFit::DrawOption("L"), RooFit::LineStyle(f_dump_Benr->GetLineStyle()), RooFit::LineColor(f_dump_Benr->GetLineColor()), RooFit::LineWidth(f_dump_Benr->GetLineWidth()));
  vdsh[1]->plotOn(frempty_benr, RooFit::Name("dshist"), RooFit::Binning(fitX::NBIN), RooFit::MarkerSize(h_dump->GetMarkerSize()), RooFit::MarkerStyle(h_dump->GetMarkerStyle()), RooFit::LineColor(h_dump->GetLineColor()), RooFit::LineWidth(h_dump->GetLineWidth()), RooFit::XErrorSize(0));
  frempty_benr->SetMinimum(0); hBenr->SetMinimum(0);
  frempty_benr->SetMaximum(75); hBenr->SetMaximum(75); // maxy[1]
  frempty_benr->Draw();
  leg_both->Draw();
  TH1F* hforpull_benr = fitX::createhistforpull(frempty_benr, "dshist", frf[1]);
  fitX::drawpull(hforpull_benr, frf[1], fitX::color_data2);
  // fitX::drawpull(hBenr, frf[1], fitX::color_data2);
  xjjroot::drawbox(4.02, hBenr->GetMinimum(), 5, hBenr->GetMaximum(), kWhite, 1);
  xjjroot::drawtex(0.99, 0.53, "Pull", frempty_incl->GetYaxis()->GetTitleSize(), 33, 42, fitX::color_data2, 270);
  xjjroot::drawtex(0.14, 0.90, "B-enriched (l#scale[0.8]{xy} > 0.1 mm)", 0.06, 12, 52, kBlack);
  cr->SaveAs(Form("plots/%s/chmassr_both_2p.pdf", output.c_str()));

  // write
  std::string outputname(Form("rootfiles/%s/fitX_fithist.root", output.c_str()));
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  h->Write();
  hBenr->Write();
  hmcp_a->Write();
  hmcp_b->Write();
  hyield_a->Write();
  hyield_b->Write();
  hyield->Write();
  hBenryield_a->Write();
  hBenryield_b->Write();
  hBenryield->Write();
  hyieldprompt_a->Write();
  hyieldprompt_b->Write();
  hyieldprompt->Write();
  hyieldpromptCorr_a->Write();
  hyieldpromptCorr_b->Write();
  hyieldpromptCorr->Write();
  mceff_a->heff_incl()->Write();
  mceff_b->heff_incl()->Write();
  heff_incl->Write();
  mceff_a->greff_incl()->Write();
  mceff_b->greff_incl()->Write();
  grfprompt_a->Write();
  grfprompt_b->Write();
  hscale_htnp_total_nominal_a->Write();
  hscale_htnp_total_nominal_b->Write();
  hratio->Write();
  fitX::write();
  outf->Close();
  std::cout<<"output: "<<outputname<<std::endl;
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  fitX::init(TFile::Open(argv[1]));
  std::string inputnametnp_a = "rootfiles/"+std::string(argv[2])+fitX::tagname()+"/drawtnp_a.root";
  std::string inputnametnp_b = "rootfiles/"+std::string(argv[2])+fitX::tagname()+"/drawtnp_b.root";
  std::string outputname = std::string(argv[2])+fitX::tagname();
  if(argc==4) { fitX_fithist(argv[1], outputname, inputnametnp_a, inputnametnp_b, argv[3]); return 0; }
  if(argc==3) { fitX_fithist(argv[1], outputname, inputnametnp_a, inputnametnp_b); return 0; }
  return 1;
}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::pttag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.74, fitX::centtag().c_str(), 0.04, 32, 42);
}

void drawlabels()
{
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
}

void drawval(TEfficiency* greff_a, TEfficiency* greff_b, float ymax)
{
  float effval_a = greff_a->GetEfficiency(fitX::ibin_a);
  xjjroot::drawtex(0.42, effval_a/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1,
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-2}", greff_a->GetEfficiency(fitX::ibin_a)*1.e+2, greff_a->GetEfficiencyErrorUp(fitX::ibin_a)*1.e+2, greff_a->GetEfficiencyErrorLow(fitX::ibin_a)*1.e+2), 0.042, 22, 62, fitX::color_a);
  float effval_b = greff_b->GetEfficiency(fitX::ibin_b);
  xjjroot::drawtex(0.72, effval_b/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1,
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-2}", greff_b->GetEfficiency(fitX::ibin_b)*1.e+2, greff_b->GetEfficiencyErrorUp(fitX::ibin_b)*1.e+2, greff_b->GetEfficiencyErrorLow(fitX::ibin_b)*1.e+2), 0.042, 22, 62, fitX::color_b);
}
