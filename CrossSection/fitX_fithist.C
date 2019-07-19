#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"
#include "lxydis.h"
#include "fitX.h"
#include "MCefficiency.h"

float weight_ss = 1./0.681;

void drawkinematic();

void fitX_fithist(std::string input, std::string output)
{
  TFile* inf = new TFile(Form("%s.root", input.c_str()));
  TH1F* h = (TH1F*)inf->Get("h");
  TH1F* hBenr = (TH1F*)inf->Get("hBenr");
  TH1F* hmcp_a = (TH1F*)inf->Get("hmcp_a");
  TH1F* hmcp_b = (TH1F*)inf->Get("hmcp_b");
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
      // std::vector<TF1*> funsft = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, "plots/fltm", false, false); // fix mean = false
      std::vector<TF1*> funs   = fitX::fit(vh[l], 0, hmcp_a, hmcp_b, Form("plots/%s", output.c_str(), vname[l].c_str()), true, true, "_"+vname[l]); // fix mean = true
      cy->cd(l+1);
      xjjroot::setgstyle();
      // <====
      float ysig1 = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
      float ysig1err = funs[0]->GetParError(5)*ysig1/funs[0]->GetParameter(5);
      float ysig2 = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH;
      float ysig2err = funs[0]->GetParError(10)*ysig2/funs[0]->GetParameter(10);
      // float ybkg1 = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH;
      // float ybkg2 = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH;

      // yield
      xjjroot::setthgrstyle(vhyield_a[l], fitX::color_a, mstyle[l], 1.2, fitX::color_a, 2, 3, fitX::color_a, 0.1, 1001);
      vhyield_a[l]->SetBinContent(2, ysig1);
      vhyield_a[l]->SetBinError(2, ysig1err);
      xjjroot::setthgrstyle(vhyield_b[l], fitX::color_b, mstyle[l], 1.2, fitX::color_b, 2, 3, fitX::color_b, 0.1, 1001);
      vhyield_b[l]->SetBinContent(4, ysig2);
      vhyield_b[l]->SetBinError(4, ysig2err);
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
      xjjroot::drawtex(0.42, ysig1/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", ysig1, ysig1err), 0.042, 22, 62, fitX::color_a);
      xjjroot::drawtex(0.72, ysig2/ymax*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", ysig2, ysig2err), 0.042, 22, 62, fitX::color_b);
      drawkinematic();
      xjjroot::drawtex(0.22, 0.84, vtitle[l].c_str(), 0.042, 12, 62);
      xjjroot::drawCMS();
    }
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
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";;f_{prompt}", 5, 0, 5, 10, 0, 1.3);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  hemptyfprompt->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptyfprompt->GetXaxis()->SetBinLabel(4, "X(3872)");
  hemptyfprompt->GetXaxis()->SetLabelSize(hemptyfprompt->GetXaxis()->GetLabelSize()*1.5);
  hemptyfprompt->Draw();
  xjjroot::drawline(0, 1, 5, 1, kGray+1, 9, 2);
  grfprompt_a->Draw("same");
  grfprompt_b->Draw("same");
  drawkinematic();
  xjjroot::drawCMS();
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
  xjjroot::drawtex(0.92, 0.84, Form("|y| < %s", xjjc::number_remove_zero(fitX::ycut).c_str()), 0.042, 32, 42);
  ceff->SaveAs(Form("plots/%s/ceff.pdf", output.c_str()));

  // correct eff
  xjjroot::setgstyle(2);
  TH1F* hyieldpromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldpromptCorr_a");
  hyieldpromptCorr_a->Divide(mceff_a.heff_incl);
  TH1F* hyieldpromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldpromptCorr_b");
  hyieldpromptCorr_b->Divide(mceff_b.heff_incl);
  float ymaxyieldpromptCorr = hyieldpromptCorr_b->GetMaximum()*2.5; // !
  TH2F* hemptyyieldpromptCorr = new TH2F("hemptyyieldpromptCorr", ";;N_{signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", 5, 0, 5, 10, 0, ymaxyieldpromptCorr);
  xjjroot::sethempty(hemptyyieldpromptCorr, 0, 0.3);
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(2, "#psi(2S)");
  hemptyyieldpromptCorr->GetXaxis()->SetBinLabel(4, "X(3872)");
  hemptyyieldpromptCorr->GetXaxis()->SetLabelSize(hemptyyieldpromptCorr->GetXaxis()->GetLabelSize()*1.5);
  TCanvas* cyieldpromptCorr = new TCanvas("cyieldpromptCorr", "", 600, 600);
  hemptyyieldpromptCorr->Draw();
  hyieldpromptCorr_a->Draw("ple same");
  hyieldpromptCorr_b->Draw("ple same");
  float yyieldpromptCorr_a = hyieldpromptCorr_a->GetBinContent(fitX::ibin_a);
  float yerryieldpromptCorr_a = hyieldpromptCorr_a->GetBinError(fitX::ibin_a);
  xjjroot::drawtex(0.42, yyieldpromptCorr_a/ymaxyieldpromptCorr*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", yyieldpromptCorr_a, yerryieldpromptCorr_a), 0.042, 22, 62, fitX::color_a);
  float yyieldpromptCorr_b = hyieldpromptCorr_b->GetBinContent(fitX::ibin_b);
  float yerryieldpromptCorr_b = hyieldpromptCorr_b->GetBinError(fitX::ibin_b);
  xjjroot::drawtex(0.72, yyieldpromptCorr_b/ymaxyieldpromptCorr*(1-gStyle->GetPadBottomMargin()-gStyle->GetPadTopMargin()) + gStyle->GetPadBottomMargin() + 0.1, Form("%.0f #pm %.0f", yyieldpromptCorr_b, yerryieldpromptCorr_b), 0.042, 22, 62, fitX::color_b);
  xjjroot::drawCMS();
  drawkinematic();
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

  TH1F* hratio_a = new TH1F("hratio_a", ";;", 1, 15, 50);
  hratio_a->SetBinContent(1, hyieldpromptCorr_a->GetBinContent(fitX::ibin_a));
  hratio_a->SetBinError(1, hyieldpromptCorr_a->GetBinError(fitX::ibin_a));
  TH1F* hratio_b = new TH1F("hratio_b", ";;", 1, 15, 50);
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
  outf->Close();

  std::cout<<std::endl<<"output: "<<Form("rootfiles/%s/fitX_fithist.root", output.c_str())<<std::endl;
}

int main(int argc, char* argv[])
{
  if(argc==3) { fitX_fithist(argv[1], argv[2]); return 0; }
  return 1;
}

void drawkinematic()
{
  xjjroot::drawtex(0.92, 0.84, Form("|y| < %s", xjjc::number_remove_zero(fitX::ycut).c_str()), 0.042, 32, 42);
  xjjroot::drawtex(0.92, 0.77, Form("p_{T} > %s GeV/c", xjjc::number_remove_zero(fitX::ptcut).c_str()), 0.042, 32, 42);
}
