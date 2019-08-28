#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fit.h"
#include "MCefficiency.h"
#include "systematics.h"

void drawkinematic();

void fitX_fithist(std::string input, std::string output, std::string inputtnp_a, std::string inputtnp_b)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;
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
  hmcp_a->Scale(hmcp_a->GetEntries()/hmcp_a->Integral());
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
  hmcp_b->Scale(hmcp_b->GetEntries()/hmcp_b->Integral());
  TH1F* hlxymcnp_a = (TH1F*)inf->Get("hlxymcnp_a");
  TH1F* hlxymcnp_b = (TH1F*)inf->Get("hlxymcnp_b");
  TH1F* hlxymcp_a = (TH1F*)inf->Get("hlxymcp_a");
  TH1F* hlxymcp_b = (TH1F*)inf->Get("hlxymcp_b");
  MCeff::MCefficiency mceff_a(inf, "_a");
  MCeff::MCefficiency mceff_b(inf, "_b");

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
  std::vector<std::string> vtitle = {"Inclusive", "B-enriched (l_{xy} > 0.1 mm)"};
  std::vector<Style_t> mstyle = {20, 24};

  xjjroot::setgstyle();
  TCanvas* cy = new TCanvas("cy", "", 1200, 600);
  cy->Divide(2, 1);
  for(int l=0; l<2; l++)
    {
      // ====>
      std::map<std::string, fitX::fitXresult*> result = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, 
                                                                  vdsh[l], dshmcp_a, dshmcp_b,
                                                                  Form("plots/%s", output.c_str(), vname[l].c_str()), l==1, true, "_"+vname[l], vtitle[l]);
      cy->cd(l+1);
      xjjroot::setgstyle();
      // <====
      float ysig_a = result["unbinned"]->ysig_a();
      float ysigerr_a = result["unbinned"]->ysigerr_a();
      float ysig_b = result["unbinned"]->ysig_b();
      float ysigerr_b = result["unbinned"]->ysigerr_b();

      // yield
      xjjroot::setthgrstyle(vhyield_a[l], fitX::color_a, mstyle[l], 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
      vhyield_a[l]->SetBinContent(2, ysig_a);
      vhyield_a[l]->SetBinError(2, ysigerr_a);
      xjjroot::setthgrstyle(vhyield_b[l], fitX::color_b, mstyle[l], 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
      vhyield_b[l]->SetBinContent(4, ysig_b);
      vhyield_b[l]->SetBinError(4, ysigerr_b);
      // float ymax = 300.; // !
      float ymax = vhyield_a[l]->GetMaximum()*2.5; // !
      TH2F* hempty = new TH2F(Form("hempty%d",l), ";;Raw Yield", 5, 0, 5, 10, 0, ymax);
      xjjroot::sethempty(hempty, 0, 0.3);
      hempty->GetXaxis()->SetBinLabel(2, "#psi(2S)");
      hempty->GetXaxis()->SetBinLabel(4, "X(3872)");
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
  xjjroot::drawline(0.1, 0, 0.1, hlxymcp_a->GetMaximum(), kGray+1, 2, 2);
  drawkinematic();
  xjjroot::drawtex(0.23, 0.85, "#psi(2S)", 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, "X(3872)", 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.85, "Prompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  clxy->cd(2);
  gPad->SetLogy();
  hlxymcnp_a->Draw("histe");
  hlxymcnp_b->Draw("histe same");
  xjjroot::drawline(0.1, 0, 0.1, hlxymcnp_a->GetMaximum(), kGray+1, 2, 2);
  drawkinematic();
  xjjroot::drawtex(0.23, 0.85, "#psi(2S)", 0.042, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.23, 0.78, "X(3872)", 0.042, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.57, 0.85, "Nonprompt", 0.042, 22, 62);
  xjjroot::drawCMS("Simulation");
  xjjroot::mkdir(Form("plots/%s/clxy.pdf", output.c_str()));
  clxy->SaveAs(Form("plots/%s/clxy.pdf", output.c_str()));
  
  std::vector<double> vlxyfrac = lxydis::nplxyfrac(hlxymcnp_a, hlxymcnp_b);
  TH1F* hlxyfrac_a = new TH1F("hlxyfrac_a", "", 5, 0, 5); //hlxyfrac_a->Sumw2();
  hlxyfrac_a->SetBinContent(2, vlxyfrac[0]);
  hlxyfrac_a->SetBinError(2, vlxyfrac[1]);
  TH1F* hlxyfrac_b = new TH1F("hlxyfrac_b", "", 5, 0, 5); //hlxyfrac_b->Sumw2();
  hlxyfrac_b->SetBinContent(4, vlxyfrac[2]);
  hlxyfrac_b->SetBinError(4, vlxyfrac[3]);
  TH1F *hyieldprompt_a, *hyieldprompt_b;
  TEfficiency* grfprompt_a = lxydis::calclxyfprompt(hyield_a, hBenryield_a, hlxyfrac_a, "grfprompt_a", &hyieldprompt_a);
  TEfficiency* grfprompt_b = lxydis::calclxyfprompt(hyield_b, hBenryield_b, hlxyfrac_b, "grfprompt_b", &hyieldprompt_b);

  xjjroot::setthgrstyle(grfprompt_a, fitX::color_a, 25, 1.5, fitX::color_a, 2, 3);
  xjjroot::setthgrstyle(grfprompt_b, fitX::color_b, 25, 1.5, fitX::color_b, 2, 3);
  xjjroot::setgstyle(2);
  TCanvas* cfprompt = new TCanvas("cfprompt", "", 600, 600);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", 5, 0, 5, 10, 0, 1.4);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  hemptyfprompt->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptyfprompt->GetXaxis()->SetBinLabel(4, "X(3872)");
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.5);
  hemptyfprompt->Draw();
  xjjroot::drawline(0, 1, 5, 1, kGray+1, 9, 2);
  grfprompt_a->Draw("same");
  grfprompt_b->Draw("same");
  drawkinematic();
  fitX::drawcomment(output.c_str());
  xjjroot::drawCMS();
  xjjroot::mkdir(Form("plots/%s/cfprompt.pdf", output.c_str()));
  cfprompt->SaveAs(Form("plots/%s/cfprompt.pdf", output.c_str()));

  // efficiency
  mceff_a.calceff();
  mceff_a.setstyle(fitX::color_a);
  mceff_b.calceff();
  mceff_b.setstyle(fitX::color_b);
  float ymaxeff = 0.2;
  TH2F* hemptyeff = new TH2F("hemptyeff", ";p_{T} (GeV/c);#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 10, MCeff::ptBins[0], MCeff::ptBins[MCeff::nPtBins], 10, 0, ymaxeff);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  TH2F* hemptyeff_incl = new TH2F("hemptyeff_incl", ";;#alpha #times #epsilon_{reco} #times #epsilon_{sel}", 5, 0, 5, 10, 0, ymaxeff);
  xjjroot::sethempty(hemptyeff_incl, 0, 0.3);
  hemptyeff_incl->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptyeff_incl->GetXaxis()->SetBinLabel(4, "X(3872)");
  hemptyeff_incl->GetXaxis()->SetLabelSize(hemptyeff_incl->GetXaxis()->GetLabelSize()*1.5);
  TLegend* legeff = new TLegend(0.70, 0.20, 1.10, 0.32);
  xjjroot::setleg(legeff, 0.042);
  legeff->AddEntry(mceff_a.greff, "#psi(2S)", "fl");
  legeff->AddEntry(mceff_b.greff, "X(3872)", "fl");
  xjjroot::setgstyle();
  TCanvas* ceff = new TCanvas("ceff", "", 1200, 600);
  ceff->Divide(2, 1);
  ceff->cd(1);
  hemptyeff->Draw();
  mceff_a.greff->Draw("same3");
  mceff_a.greff->Draw("samelX");
  mceff_b.greff->Draw("same3");
  mceff_b.greff->Draw("samelX");
  legeff->Draw();
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);

  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::centtag().c_str(), 0.04, 32, 42);
  fitX::drawcomment(output.c_str());
  ceff->cd(2);
  hemptyeff_incl->Draw();
  mceff_a.greff_incl->Draw("same ple");
  mceff_b.greff_incl->Draw("same ple");
  float effval_a = mceff_a.greff_incl->GetEfficiency(fitX::ibin_a);
  xjjroot::drawtex(0.42, effval_a/ymaxeff*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, 
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", mceff_a.greff_incl->GetEfficiency(fitX::ibin_a)*1000, mceff_a.greff_incl->GetEfficiencyErrorUp(fitX::ibin_a)*1000, mceff_a.greff_incl->GetEfficiencyErrorLow(fitX::ibin_a)*1000), 0.042, 22, 62, fitX::color_a);
  float effval_b = mceff_b.greff_incl->GetEfficiency(fitX::ibin_b);
  xjjroot::drawtex(0.72, effval_b/ymaxeff*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1,
                   Form("%.1f {}^{+ %.1f}_{-  %.1f} #times 10^{-3}", mceff_b.greff_incl->GetEfficiency(fitX::ibin_b)*1000, mceff_b.greff_incl->GetEfficiencyErrorUp(fitX::ibin_b)*1000, mceff_b.greff_incl->GetEfficiencyErrorLow(fitX::ibin_b)*1000), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawtex(0.23, 0.84, "PYTHIA8 + HYDJET", 0.042, 12, 62);
  xjjroot::drawtex(0.23, 0.77, "Prompt", 0.042, 12, 62);
  drawkinematic();
  xjjroot::drawCMSleft("#scale[1.25]{#bf{CMS}} #it{Simulation}");
  xjjroot::drawCMSright();
  xjjroot::mkdir(Form("plots/%s/ceff.pdf", output.c_str()));
  ceff->SaveAs(Form("plots/%s/ceff.pdf", output.c_str()));

  // tnp
  TFile* inftnp_a = TFile::Open(inputtnp_a.c_str());
  TH1D* hscale_htnp_total_nominal_a = (TH1D*)inftnp_a->Get("scale_htnp_total_nominal");
  hyieldprompt_a->Scale(hscale_htnp_total_nominal_a->GetBinContent(1));
  TFile* inftnp_b = TFile::Open(inputtnp_b.c_str());
  TH1D* hscale_htnp_total_nominal_b = (TH1D*)inftnp_b->Get("scale_htnp_total_nominal");
  hyieldprompt_b->Scale(hscale_htnp_total_nominal_b->GetBinContent(1));
  // correct eff
  xjjroot::setgstyle(2);
  TH1F* hyieldpromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldpromptCorr_a");
  hyieldpromptCorr_a->Divide(mceff_a.heff_incl);
  TH1F* hyieldpromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldpromptCorr_b");
  hyieldpromptCorr_b->Divide(mceff_b.heff_incl);
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
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(4, "X(3872)");
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
  fitX::drawcomment(output.c_str());
  xjjroot::mkdir(Form("plots/%s/cyieldpromptCorr.pdf", output.c_str()));
  cyieldpromptCorr->SaveAs(Form("plots/%s/cyieldpromptCorr.pdf", output.c_str()));

  // merge a + b
  TH1F* hyield = (TH1F*)hyield_a->Clone("hyield");
  hyield->Add(hyield_b);
  TH1F* hBenryield = (TH1F*)hBenryield_a->Clone("hBenryield");
  hBenryield->Add(hBenryield_b);
  TH1F* heff_incl = (TH1F*)mceff_a.heff_incl->Clone("heff_incl");
  heff_incl->Add(mceff_b.heff_incl);
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

  // write
  TFile* outf = new TFile(Form("rootfiles/%s/fitX_fithist.root", output.c_str()), "recreate");
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
  mceff_a.heff_incl->Write();
  mceff_b.heff_incl->Write();
  heff_incl->Write();
  mceff_a.greff_incl->Write();
  mceff_b.greff_incl->Write();
  grfprompt_a->Write();
  grfprompt_b->Write();
  hratio->Write();
  fitX::write();
  outf->Close();
  std::cout<<"output: "<<Form("rootfiles/%s/fitX_fithist.root", output.c_str())<<std::endl;
  std::cout<<std::endl;
}

int main(int argc, char* argv[])
{
  fitX::init(TFile::Open(argv[1]));
  std::string inputnametnp_a = "rootfiles/"+std::string(argv[2])+fitX::tagname()+"/drawtnp_a.root";
  std::string inputnametnp_b = "rootfiles/"+std::string(argv[2])+fitX::tagname()+"/drawtnp_b.root";
  std::string outputname = std::string(argv[2])+fitX::tagname();
  if(argc==3) { fitX_fithist(argv[1], outputname, inputnametnp_a, inputnametnp_b); return 0; }
  return 1;
}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, fitX::ytag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.79, fitX::pttag().c_str(), 0.04, 32, 42);
  xjjroot::drawtex(0.92, 0.74, fitX::centtag().c_str(), 0.04, 32, 42);
}
