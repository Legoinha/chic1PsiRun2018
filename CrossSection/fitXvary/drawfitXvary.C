#include "fitX.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitXvary.h"
#include "lxydis.h"

bool saveornot = true; // save separate fit or not
void drawfitXvary(std::string output)
{
  xjjroot::setgstyle(3);

  fitX::varycut vc(bdtg, "BDTG");

  TFile* inf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()));
  // read
  for(int l=0, lp=0; l<nbdtg; l++)
    {
      vc.hdata[l] = (TH1F*)inf->Get(Form("hdata_%d", l));
      fitX::setmasshist(vc.hdata[l], 0, -0.1);
      vc.hdata[l]->SetMinimum(0);
      xjjroot::setthgrstyle(vc.hdata[l], xjjroot::colorlist_dark[lp], 20, 0.9, xjjroot::colorlist_dark[lp], 1, 1);
      vc.hdataBenr[l] = (TH1F*)inf->Get(Form("hdataBenr_%d", l));
      fitX::setmasshist(vc.hdataBenr[l], 0, -0.1);
      vc.hdataBenr[l]->SetMinimum(0);
      xjjroot::setthgrstyle(vc.hdataBenr[l], xjjroot::colorlist_dark[lp], 20, 0.9, xjjroot::colorlist_dark[lp], 1, 1);
      vc.hmc_a[l] = (TH1F*)inf->Get(Form("hmc_a_%d", l));
      vc.hmc_b[l] = (TH1F*)inf->Get(Form("hmc_b_%d", l));

      vc.hlxymcnp_a[l] = (TH1F*)inf->Get(Form("hlxymcnp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_a[l], xjjroot::colorlist_dark[lp], 20, 0.2, xjjroot::colorlist_dark[lp], 1, 1);
      vc.hlxymcnp_b[l] = (TH1F*)inf->Get(Form("hlxymcnp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_b[l], xjjroot::colorlist_dark[lp], 20, 0.2, xjjroot::colorlist_dark[lp], 1, 1);
      vc.hlxymcp_a[l] = (TH1F*)inf->Get(Form("hlxymcp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_a[l], xjjroot::colorlist_dark[lp], 20, 0.2, xjjroot::colorlist_dark[lp], 1, 1);
      vc.hlxymcp_b[l] = (TH1F*)inf->Get(Form("hlxymcp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_b[l], xjjroot::colorlist_dark[lp], 20, 0.2, xjjroot::colorlist_dark[lp], 1, 1);

      if(pbdtg[l]) lp++;
    }
  vc.heff_a = (TH1F*)inf->Get("heff_a");
  vc.heff_b = (TH1F*)inf->Get("heff_b");
  vc.greff_a = (TEfficiency*)inf->Get("greff_a");
  vc.greff_b = (TEfficiency*)inf->Get("greff_b");
  vc.hsideband_a = (TH1F*)inf->Get("hsideband_a");
  vc.hsideband_b = (TH1F*)inf->Get("hsideband_b");

  // ==> yield
  TH1F* hyieldbdtg_a = new TH1F("hyieldbdtg_a"        , ";BDTG;N_{Signal}"                 , bdtg.size()-1, bdtg.data());
  TH1F* hyieldbdtg_b = new TH1F("hyieldbdtg_b"        , ";BDTG;N_{Signal}"                 , bdtg.size()-1, bdtg.data());
  TH1F* hbkgbdtg_a = new TH1F("hbkgbdtg_a"            , ";BDTG;N_{Background}"             , bdtg.size()-1, bdtg.data());
  TH1F* hbkgbdtg_b = new TH1F("hbkgbdtg_b"            , ";BDTG;N_{Background}"             , bdtg.size()-1, bdtg.data());
  TH1F* hsigfitbdtg_a = new TH1F("hsigfitbdtg_a"      , ";BDTG;S / #sqrt{S+B}"             , bdtg.size()-1, bdtg.data());
  TH1F* hsigfitbdtg_b = new TH1F("hsigfitbdtg_b"      , ";BDTG;S / #sqrt{S+B}"             , bdtg.size()-1, bdtg.data());
  TH1F* hyieldbdtgBenr_a = new TH1F("hyieldbdtgBenr_a", ";BDTG;N_{Signal} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
  TH1F* hyieldbdtgBenr_b = new TH1F("hyieldbdtgBenr_b", ";BDTG;N_{Signal} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
  TH1F* hlxyfracbdtg_a = new TH1F("hlxyfracbdtg_a"    , ";BDTG;"                           , bdtg.size()-1, bdtg.data());
  TH1F* hlxyfracbdtg_b = new TH1F("hlxyfracbdtg_b"    , ";BDTG;"                           , bdtg.size()-1, bdtg.data());

  TCanvas* c = new TCanvas("cvary", "", 800, 1200);
  c->Divide(1, 2);
  std::vector<std::vector<TH1F*>*> vhdata        = {&vc.hdata,       &vc.hdataBenr};
  std::vector<TH1F*>               vhyieldbdtg_a = {hyieldbdtg_a, hyieldbdtgBenr_a};
  std::vector<TH1F*>               vhyieldbdtg_b = {hyieldbdtg_b, hyieldbdtgBenr_b};
  std::vector<std::string>         vtitle        = {"Inclusive",  "l_{xy} > 100 #mum"};
  for(int cc=0; cc<vhdata.size(); cc++)
    {
      c->cd(cc+1);
      for(int l=0; l<nbdtg; l++) { if(pbdtg[l]) { vhdata[cc]->at(l)->Draw("pe"); break; } }
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          c->cd(cc+1);
          if(pbdtg[l]) vhdata[cc]->at(l)->Draw("pe same");
          std::vector<TF1*> funs = fitX::fit(xjjroot::copyobject(vhdata[cc]->at(l), Form("%s_%d", vhdata[cc]->at(l)->GetName(), l)), 0, vc.hmc_a[0], vc.hmc_b[0], 
                                             Form("idx/%s/%s", output.c_str(), vhdata[cc]->at(l)->GetName()), true, saveornot); // fix mean
          c->cd(cc+1);
          xjjroot::settfstyle(funs[0], xjjroot::colorlist_dark[lp], 2, 2);
          if(pbdtg[l])
            {
              xjjroot::copyobject(funs[0], Form("ff_%s_%d", vhdata[cc]->at(l)->GetName(), l))->Draw("l same");
              xjjroot::drawtex(0.70+0.08*(lp%3), 0.75-0.04*(lp/3), Form("> %.2f%s", bdtg[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[lp]);
              lp++;
            }
          xjjroot::drawtex(0.90, 0.85, vtitle[cc].c_str(), 0.04, 33, 42, kBlack);
          if(l==nbdtg-1) continue;
          float ysig_a = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysig_a < 0) ysig_a = 0;
          float ysig_aerr = funs[0]->GetParError(5)*ysig_a/funs[0]->GetParameter(5);
          float ysig_b = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysig_b < 0) ysig_b = 0;
          float ysig_berr = funs[0]->GetParError(10)*ysig_b/funs[0]->GetParameter(10);
          float ybkg_a = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH; if(ybkg_a < 0) ybkg_a = 0;
          float ybkg_b = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH; if(ybkg_b < 0) ybkg_b = 0;
          float ysigf_a = (ysig_a+ybkg_a>0)?(ysig_a/TMath::Sqrt(ysig_a+ybkg_a)):0;
          float ysigf_b = (ysig_b+ybkg_b>0)?(ysig_b/TMath::Sqrt(ysig_b+ybkg_b)):0;
          vhyieldbdtg_a[cc]->SetBinContent(l+1, ysig_a);
          vhyieldbdtg_a[cc]->SetBinError(l+1, ysig_aerr);
          vhyieldbdtg_b[cc]->SetBinContent(l+1, ysig_b);
          vhyieldbdtg_b[cc]->SetBinError(l+1, ysig_berr);

          // only for inclusive
          if(cc) continue; 
          hbkgbdtg_a->SetBinContent(l+1, ybkg_a);
          hbkgbdtg_b->SetBinContent(l+1, ybkg_b);
          hsigfitbdtg_a->SetBinContent(l+1, ysigf_a);
          hsigfitbdtg_b->SetBinContent(l+1, ysigf_b);
        }
      xjjroot::setgstyle(1);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
    }
  c->SaveAs(Form("plots/%s/cmass_varybdtg.pdf", output.c_str()));

  // ==> significance
  const int iref = 16;
  float yieldnocut_a = hyieldbdtg_a->GetBinContent(iref+1) / vc.heff_a->GetBinContent(iref+1);
  float yieldnocut_b = hyieldbdtg_b->GetBinContent(iref+1) / vc.heff_b->GetBinContent(iref+1);
  TH1F* hsigbdtg_a = new TH1F("hsigbdtg_a", ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());  
  TH1F* hsigbdtg_b = new TH1F("hsigbdtg_b", ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());
  for(int l=0; l<nbdtg-1; l++)
    {
      float nsig_a = yieldnocut_a * vc.heff_a->GetBinContent(l+1);
      float nsig_b = yieldnocut_b * vc.heff_b->GetBinContent(l+1);
      float nbkg_a = vc.hsideband_a->GetBinContent(l+1) / (masswinH-masswinL) * 0.0070*3;
      float nbkg_b = vc.hsideband_b->GetBinContent(l+1) / (masswinH-masswinL) * 0.0093*3;
      float significance_a = nsig_a / TMath::Sqrt(nsig_a + nbkg_a);
      float significance_b = nsig_b / TMath::Sqrt(nsig_b + nbkg_b);
      hsigbdtg_a->SetBinContent(l+1, significance_a);
      hsigbdtg_b->SetBinContent(l+1, significance_b);
    }

  // ==> treatment for low-stat X
  TH1F* hyieldbdtgerase_b = new TH1F("hyieldbdtgerase_b", ";BDTG;N_{Signal}", nbdtg-1, bdtg.data());
  TH1F* hyieldbdtgeraseBenr_b = new TH1F("hyieldbdtgeraseBenr_b", ";BDTG;N_{Signal}", nbdtg-1, bdtg.data());
  for(int l=0; l<hyieldbdtg_b->GetNbinsX(); l++)
    {
      if(bdtg[l] < 0.5) continue;
      hyieldbdtgerase_b->SetBinContent(l+1, hyieldbdtg_b->GetBinContent(l+1));
      hyieldbdtgerase_b->SetBinError(l+1, hyieldbdtg_b->GetBinError(l+1));
      hyieldbdtgeraseBenr_b->SetBinContent(l+1, hyieldbdtgBenr_b->GetBinContent(l+1));
      hyieldbdtgeraseBenr_b->SetBinError(l+1, hyieldbdtgBenr_b->GetBinError(l+1));
    }

  // ==> fprompt
  TCanvas* cp = new TCanvas("cp", "", 1600, 1200);
  cp->Divide(2, 2);
  for(int l=0; l<nbdtg-1; l++)
    {
      vc.hlxymcp_a[l]->Scale(1./vc.hlxymcp_a[l]->Integral(), "width");
      vc.hlxymcp_b[l]->Scale(1./vc.hlxymcp_b[l]->Integral(), "width");
      vc.hlxymcnp_a[l]->Scale(1./vc.hlxymcnp_a[l]->Integral(), "width");
      vc.hlxymcnp_b[l]->Scale(1./vc.hlxymcnp_b[l]->Integral(), "width");
      std::vector<double> vlxyfrac = lxydis::nplxyfrac(vc.hlxymcnp_a[l], vc.hlxymcnp_b[l]);
      hlxyfracbdtg_a->SetBinContent(l+1, vlxyfrac[0]);
      hlxyfracbdtg_a->SetBinError(l+1, vlxyfrac[1]);
      hlxyfracbdtg_b->SetBinContent(l+1, vlxyfrac[2]);
      hlxyfracbdtg_b->SetBinError(l+1, vlxyfrac[3]);

      if(pbdtg[l])
        {
          cp->cd(1);
          vc.hlxymcp_a[l]->Draw((l?"histe same":"histe"));
          cp->cd(2);
          vc.hlxymcp_b[l]->Draw((l?"histe same":"histe"));
          cp->cd(3);
          vc.hlxymcnp_a[l]->Draw((l?"histe same":"histe"));
          cp->cd(4);
          vc.hlxymcnp_b[l]->Draw((l?"histe same":"histe"));
        }
    }
  cp->SaveAs(Form("plots/%s/clxy_varybdtg.pdf", output.c_str()));
  TH1F *hyieldprompt_a, *hyieldprompt_b;
  TEfficiency* grfprompt_a = lxydis::calclxyfprompt(hyieldbdtg_a, hyieldbdtgBenr_a, hlxyfracbdtg_a, "grfprompt_a", &hyieldprompt_a);
  TEfficiency* grfprompt_b = lxydis::calclxyfprompt(hyieldbdtgerase_b, hyieldbdtgeraseBenr_b, hlxyfracbdtg_b, "grfprompt_b", &hyieldprompt_b);

  // ==> corrected yield
  TH1F* hyieldCorr_a = (TH1F*)hyieldbdtg_a->Clone("hyieldCorr_a");
  hyieldCorr_a->Divide(vc.heff_a);
  TH1F* hyieldCorr_b = (TH1F*)hyieldbdtgerase_b->Clone("hyieldCorr_b");
  hyieldCorr_b->Divide(vc.heff_b);
  TH1F* hyieldPromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldCorr_a");
  hyieldPromptCorr_a->Divide(vc.heff_a);
  TH1F* hyieldPromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldCorr_b");
  hyieldPromptCorr_b->Divide(vc.heff_b);

  /*================== Draw ====================*/

  // draw
  TH2F* hemptyyd = new TH2F("hemptyyd", ";BDTG;Raw N_{Signal} by Fit", 10, -1.1, 1, 10, 0, hyieldbdtg_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyyd, 0, 0.4);
  TGraphErrors* gr_hyieldbdtg_a = xjjroot::shifthistcenter(hyieldbdtg_a, "gr_hyieldbdtg_a");
  TGraphErrors* gr_hyieldbdtg_b = xjjroot::shifthistcenter(hyieldbdtg_b, "gr_hyieldbdtg_b");
  xjjroot::setthgrstyle(gr_hyieldbdtg_a, fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldbdtg_b, fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyeff = new TH2F("hemptyeff", ";BDTG;(#alpha #times #epsilon )_{prompt}", 10, -1.1, 1, 10, 0, 0.1);
  xjjroot::sethempty(hemptyeff, 0, 0.4);
  TGraphAsymmErrors* gr_greff_a = xjjroot::shifthistcenter(vc.greff_a, "gr_greff_a");
  TGraphAsymmErrors* gr_greff_b = xjjroot::shifthistcenter(vc.greff_b, "gr_greff_b");
  xjjroot::setthgrstyle(gr_greff_a, fitX::color_a, 34, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_greff_b, fitX::color_b, 34, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  TH2F* hemptysigf = new TH2F("hemptysigf", ";BDTG;S / #sqrt{S+B}", 10, -1.1, 1, 10, 0, hsigbdtg_b->GetMaximum()*1.4);
  xjjroot::sethempty(hemptysigf, 0, 0.4);
  TGraphErrors* gr_hsigfitbdtg_a = xjjroot::shifthistcenter(hsigfitbdtg_a, "gr_hsigfitbdtg_a");
  TGraphErrors* gr_hsigfitbdtg_b = xjjroot::shifthistcenter(hsigfitbdtg_b, "gr_hsigfitbdtg_b");
  TGraphErrors* gr_hsigbdtg_a = xjjroot::shifthistcenter(hsigbdtg_a, "gr_hsigbdtg_a");
  TGraphErrors* gr_hsigbdtg_b = xjjroot::shifthistcenter(hsigbdtg_b, "gr_hsigbdtg_b");
  xjjroot::setthgrstyle(gr_hsigfitbdtg_a, fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hsigfitbdtg_b, fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  xjjroot::setthgrstyle(gr_hsigbdtg_a, fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hsigbdtg_b, fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyBenr = new TH2F("hemptyBenr", ";BDTG;N_{Signal} (l_{xy} > 100#mum) by Fit", 10, -1.1, 1, 10, 0, hyieldbdtgBenr_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyBenr, 0, 0.4);
  TGraphErrors* gr_hyieldbdtgBenr_a = xjjroot::shifthistcenter(hyieldbdtgBenr_a, "gr_hyieldbdtgBenr_a");
  TGraphErrors* gr_hyieldbdtgBenr_b = xjjroot::shifthistcenter(hyieldbdtgBenr_b, "gr_hyieldbdtgBenr_b");
  xjjroot::setthgrstyle(gr_hyieldbdtgBenr_a, fitX::color_a, 24, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldbdtgBenr_b, fitX::color_b, 24, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyydprompt = new TH2F("hemptyydprompt", ";BDTG;N_{Signal} #times f_{prompt}", 10, -1.1, 1, 10, 0, hyieldprompt_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyydprompt, 0, 0.4);
  TGraphErrors* gr_hyieldprompt_a = xjjroot::shifthistcenter(hyieldprompt_a, "gr_hyieldprompt_a");
  TGraphErrors* gr_hyieldprompt_b = xjjroot::shifthistcenter(hyieldprompt_b, "gr_hyieldprompt_b");
  xjjroot::setthgrstyle(gr_hyieldprompt_a, fitX::color_a, 21, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldprompt_b, fitX::color_b, 21, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";BDTG;f_{prompt} After Cuts", 10, -1.1, 1, 10, 0, 1.2);
  xjjroot::sethempty(hemptyfprompt, 0, 0.4);
  TGraphAsymmErrors* gr_grfprompt_a = xjjroot::shifthistcenter(grfprompt_a, "gr_grfprompt_a");
  TGraphAsymmErrors* gr_grfprompt_b = xjjroot::shifthistcenter(grfprompt_b, "gr_grfprompt_b");
  xjjroot::setthgrstyle(gr_grfprompt_a, fitX::color_a, 45, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_grfprompt_b, fitX::color_b, 45, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  TH2F* hemptyCorr = new TH2F("hemptyCorr", ";BDTG;N_{Signal} / (#alpha #times #epsilon )_{prompt}", 10, -1.1, 1, 10, 0, hyieldCorr_b->GetMaximum()*1.5);
  xjjroot::sethempty(hemptyCorr, 0, 0.4);
  TGraphErrors* gr_hyieldCorr_a = xjjroot::shifthistcenter(hyieldCorr_a, "gr_hyieldCorr_a");
  TGraphErrors* gr_hyieldCorr_b = xjjroot::shifthistcenter(hyieldCorr_b, "gr_hyieldCorr_b");
  xjjroot::setthgrstyle(gr_hyieldCorr_a, fitX::color_a, 3, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldCorr_b, fitX::color_b, 3, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyPromptCorr = new TH2F("hemptyPromptCorr", ";BDTG;N_{Signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", 10, -1.1, 1, 10, 0, hyieldPromptCorr_b->GetMaximum()*1.5);
  xjjroot::sethempty(hemptyPromptCorr, 0, 0.4);
  TGraphErrors* gr_hyieldPromptCorr_a = xjjroot::shifthistcenter(hyieldPromptCorr_a, "gr_hyieldPromptCorr_a");
  TGraphErrors* gr_hyieldPromptCorr_b = xjjroot::shifthistcenter(hyieldPromptCorr_b, "gr_hyieldPromptCorr_b");
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_a, fitX::color_a, 5, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_b, fitX::color_b, 5, 1.2, fitX::color_b, 1, 1);

  // do draw
  xjjroot::setgstyle(1);
  TCanvas* c3 = new TCanvas("c3", "", 2400, 1200);
  c3->Divide(4, 2);
  c3->cd(1);
  hemptyyd->Draw();
  gr_hyieldbdtg_a->Draw("pe same");
  gr_hyieldbdtg_b->Draw("pe same");
  drawalltext();
  c3->cd(2);
  hemptyeff->Draw();
  gr_greff_a->Draw("p3e same");
  gr_greff_b->Draw("p3e same");
  drawalltext();
  c3->cd(3);
  hemptysigf->Draw();
  // gr_hsigfitbdtg_b->Draw("plX0 same");
  // gr_hsigfitbdtg_a->Draw("plX0 same");
  gr_hsigbdtg_b->Draw("plX0 same");
  gr_hsigbdtg_a->Draw("plX0 same");
  drawalltext();
  c3->cd(4);
  hemptyCorr->Draw();
  gr_hyieldCorr_a->Draw("pe same");
  gr_hyieldCorr_b->Draw("pe same");
  drawalltext();
  c3->cd(5);
  hemptyBenr->Draw();
  gr_hyieldbdtgBenr_a->Draw("pe same");
  gr_hyieldbdtgBenr_b->Draw("pe same");
  drawalltext();
  c3->cd(6);
  hemptyydprompt->Draw();
  gr_hyieldprompt_a->Draw("pe same");
  gr_hyieldprompt_b->Draw("pe same");
  drawalltext();
  c3->cd(7);
  hemptyfprompt->Draw();
  gr_grfprompt_a->Draw("pe same");
  gr_grfprompt_b->Draw("pe same");
  drawalltext();
  c3->cd(8);
  hemptyPromptCorr->Draw();
  gr_hyieldPromptCorr_a->Draw("pe same");
  gr_hyieldPromptCorr_b->Draw("pe same");
  drawalltext();
  c3->SaveAs(Form("plots/%s/csig_varybdtg.pdf", output.c_str()));

  TFile* outf = new TFile(Form("rootfiles/root_fitXvary_yield_%s.root", output.c_str()), "recreate");
  outf->cd();
  hyieldbdtg_a->Write();
  hyieldbdtg_b->Write();
  hbkgbdtg_a->Write();
  hbkgbdtg_b->Write();
  hsigfitbdtg_a->Write();
  hsigfitbdtg_b->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==2) { drawfitXvary(argv[1]); return 0; }
  return 1;
}
