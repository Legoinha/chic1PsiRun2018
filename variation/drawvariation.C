#ifndef __VARIABLE__
#define __VARIABLE__ __INPUT__FROM__OUTSIDE__
#endif

#include "fitX.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitXvary.h"
#include "lxydis.h"

#include "variation.h"

bool saveornot = false; // save separate fit or not
void drawvariaiton(std::string output, std::string vartitle, int iscutordis)
{
  xjjroot::setgstyle(3);

  std::vector<float> vvector(variation::varvectors[(iscutordis?vartitle+"cut":vartitle)]);
  fitX::varycut vc(vvector, vartitle, iscutordis);

  TFile* inf = new TFile(Form("rootfiles/root_variation_%s.root", output.c_str()));
  // read
  for(int l=0; l<vc.getnv(); l++)
    {
      vc.hdata[l] = (TH1F*)inf->Get(Form("hdata_%d", l));
      fitX::setmasshist(vc.hdata[l], 0, -0.1);
      vc.hdata[l]->SetMinimum(0);
      xjjroot::setthgrstyle(vc.hdata[l], xjjroot::colorlist_dark[l], 20, 0.9, xjjroot::colorlist_dark[l], 1, 1);
      vc.hdataBenr[l] = (TH1F*)inf->Get(Form("hdataBenr_%d", l));
      fitX::setmasshist(vc.hdataBenr[l], 0, -0.1);
      vc.hdataBenr[l]->SetMinimum(0);
      xjjroot::setthgrstyle(vc.hdataBenr[l], xjjroot::colorlist_dark[l], 20, 0.9, xjjroot::colorlist_dark[l], 1, 1);
      vc.hmc_a[l] = (TH1F*)inf->Get(Form("hmc_a_%d", l));
      vc.hmc_b[l] = (TH1F*)inf->Get(Form("hmc_b_%d", l));
      vc.hlxymcnp_a[l] = (TH1F*)inf->Get(Form("hlxymcnp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_a[l], xjjroot::colorlist_dark[l], 20, 0.2, xjjroot::colorlist_dark[l], 1, 1);
      vc.hlxymcnp_b[l] = (TH1F*)inf->Get(Form("hlxymcnp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_b[l], xjjroot::colorlist_dark[l], 20, 0.2, xjjroot::colorlist_dark[l], 1, 1);
      vc.hlxymcp_a[l] = (TH1F*)inf->Get(Form("hlxymcp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_a[l], xjjroot::colorlist_dark[l], 20, 0.2, xjjroot::colorlist_dark[l], 1, 1);
      vc.hlxymcp_b[l] = (TH1F*)inf->Get(Form("hlxymcp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_b[l], xjjroot::colorlist_dark[l], 20, 0.2, xjjroot::colorlist_dark[l], 1, 1);
    }
  vc.heff_a = (TH1F*)inf->Get("heff_a");
  vc.heff_b = (TH1F*)inf->Get("heff_b");
  vc.greff_a = (TEfficiency*)inf->Get("greff_a");
  vc.greff_b = (TEfficiency*)inf->Get("greff_b");
  vc.hsideband_a = (TH1F*)inf->Get("hsideband_a");
  vc.hsideband_b = (TH1F*)inf->Get("hsideband_b");
  vc.hmcdisp_a = (TH1F*)inf->Get("hmcdisp_a");
  vc.hmcdisp_b = (TH1F*)inf->Get("hmcdisp_b");
  vc.hmcdisnp_a = (TH1F*)inf->Get("hmcdisnp_a");
  vc.hmcdisnp_b = (TH1F*)inf->Get("hmcdisnp_b");

  // ==> yield
  TH1F* hyield_a = new TH1F("hyield_a"        , Form(";%s;N_{Signal};", vartitle.c_str())                 , vvector.size()-1, vvector.data());
  TH1F* hyield_b = new TH1F("hyield_b"        , Form(";%s;N_{Signal};", vartitle.c_str())                 , vvector.size()-1, vvector.data());
  TH1F* hbkg_a = new TH1F("hbkg_a"            , Form(";%s;N_{Background};", vartitle.c_str())             , vvector.size()-1, vvector.data());
  TH1F* hbkg_b = new TH1F("hbkg_b"            , Form(";%s;N_{Background};", vartitle.c_str())             , vvector.size()-1, vvector.data());
  TH1F* hsigfit_a = new TH1F("hsigfit_a"      , Form(";%s;S / #sqrt{S+B};", vartitle.c_str())             , vvector.size()-1, vvector.data());
  TH1F* hsigfit_b = new TH1F("hsigfit_b"      , Form(";%s;S / #sqrt{S+B};", vartitle.c_str())             , vvector.size()-1, vvector.data());
  TH1F* hyieldBenr_a = new TH1F("hyieldBenr_a", Form(";%s;N_{Signal} (l_{xy} > 0.1mm);", vartitle.c_str()), vvector.size()-1, vvector.data());
  TH1F* hyieldBenr_b = new TH1F("hyieldBenr_b", Form(";%s;N_{Signal} (l_{xy} > 0.1mm);", vartitle.c_str()), vvector.size()-1, vvector.data());
  TH1F* hlxyfrac_a = new TH1F("hlxyfrac_a"    , Form(";%s;;", vartitle.c_str())                           , vvector.size()-1, vvector.data());
  TH1F* hlxyfrac_b = new TH1F("hlxyfrac_b"    , Form(";%s;;", vartitle.c_str())                           , vvector.size()-1, vvector.data());
  std::vector<TF1*> ffun, ffunBenr;

  TCanvas* c = new TCanvas("cvary", "", 800, 1200);
  c->Divide(1, 2);
  std::vector<std::vector<TH1F*>*> vhdata    = {&vc.hdata   , &vc.hdataBenr};
  std::vector<TH1F*>               vhyield_a = {hyield_a    , hyieldBenr_a};
  std::vector<TH1F*>               vhyield_b = {hyield_b    , hyieldBenr_b};
  std::vector<std::string>         vtitle    = {"Inclusive" , "l_{xy} > 100 #mum"};
  std::vector<std::vector<TF1*>*>  vffun     = {&ffun       , &ffunBenr};
  for(int cc=0; cc<vhdata.size(); cc++)
    {
      c->cd(cc+1);
      float maxydata = 0;
      for(int l=0; l<vc.getnv(); l++) { if(vhdata[cc]->at(l)->GetMaximum() > maxydata) { maxydata = vhdata[cc]->at(l)->GetMaximum(); } }
      for(int l=0; l<vc.getnv(); l++)
        {
          c->cd(cc+1);
          if(!l) { vhdata[cc]->at(l)->SetMaximum(maxydata*1.2); vhdata[cc]->at(l)->Draw("pe"); }
          else vhdata[cc]->at(l)->Draw("pe same");
          std::vector<TF1*> funs = fitX::fit(xjjroot::copyobject(vhdata[cc]->at(l), Form("%s_%d", vhdata[cc]->at(l)->GetName(), l)), 0, vc.hmc_a[0], vc.hmc_b[0], 
                                             Form("plots/%s/idx", output.c_str()), true, saveornot, vhdata[cc]->at(l)->GetName()); // fix mean
          c->cd(cc+1);
          xjjroot::settfstyle(funs[0], xjjroot::colorlist_dark[l], 2, 2);
          vffun[cc]->push_back(xjjroot::copyobject(funs[0], Form("ff_%s_%d", vhdata[cc]->at(l)->GetName(), l)));
          vffun[cc]->back()->Draw("l same");

          // if(l==vc.getnv()-1) continue;
          if(vhdata[cc]->at(l)->GetEntries() < 1) continue;

          float ysig_a = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysig_a < 0) ysig_a = 0;
          float ysig_aerr = funs[0]->GetParError(5)*ysig_a/funs[0]->GetParameter(5);
          float ysig_b = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysig_b < 0) ysig_b = 0;
          float ysig_berr = funs[0]->GetParError(10)*ysig_b/funs[0]->GetParameter(10);
          float ybkg_a = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH; if(ybkg_a < 0) ybkg_a = 0;
          float ybkg_b = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH; if(ybkg_b < 0) ybkg_b = 0;
          float ysigf_a = (ysig_a+ybkg_a>0)?(ysig_a/TMath::Sqrt(ysig_a+ybkg_a)):0;
          float ysigf_b = (ysig_b+ybkg_b>0)?(ysig_b/TMath::Sqrt(ysig_b+ybkg_b)):0;
          vhyield_a[cc]->SetBinContent(l+1, ysig_a);
          vhyield_a[cc]->SetBinError(l+1, ysig_aerr);
          vhyield_b[cc]->SetBinContent(l+1, ysig_b);
          vhyield_b[cc]->SetBinError(l+1, ysig_berr);

          // only for inclusive
          if(cc) continue; 
          hbkg_a->SetBinContent(l+1, ybkg_a);
          hbkg_b->SetBinContent(l+1, ybkg_b);
          hsigfit_a->SetBinContent(l+1, ysigf_a);
          hsigfit_b->SetBinContent(l+1, ysigf_b);
        }
      xjjroot::drawtex(0.23, 0.85, vtitle[cc].c_str(), 0.035, 12, 62, kBlack);
      variation::drawtexlist(iscutordis, vvector, vc.getnv(), vartitle);
      xjjroot::setgstyle(1);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
    }
  c->SaveAs(Form("plots/%s/cmass_vary.pdf", output.c_str()));

  // ==> significance
  const int iref = xjjc::findibin(vvector, (float)0.75);
  float yieldnocut_a = hyield_a->GetBinContent(iref+1) / vc.heff_a->GetBinContent(iref+1);
  float yieldnocut_b = hyield_b->GetBinContent(iref+1) / vc.heff_b->GetBinContent(iref+1);
  TH1F* hsig_a = new TH1F("hsig_a", Form(";%s;S / #sqrt{S+B};", vartitle.c_str()), vvector.size()-1, vvector.data());  
  TH1F* hsig_b = new TH1F("hsig_b", Form(";%s;S / #sqrt{S+B};", vartitle.c_str()), vvector.size()-1, vvector.data());
  for(int l=0; l<vc.getnv()-1; l++)
    {
      float nsig_a = yieldnocut_a * vc.heff_a->GetBinContent(l+1);
      float nsig_b = yieldnocut_b * vc.heff_b->GetBinContent(l+1);
      float nbkg_a = vc.hsideband_a->GetBinContent(l+1) / (vc.masswinH-vc.masswinL) * 0.0070*3;
      float nbkg_b = vc.hsideband_b->GetBinContent(l+1) / (vc.masswinH-vc.masswinL) * 0.0093*3;
      float significance_a = nsig_a / TMath::Sqrt(nsig_a + nbkg_a);
      float significance_b = nsig_b / TMath::Sqrt(nsig_b + nbkg_b);
      hsig_a->SetBinContent(l+1, significance_a);
      hsig_b->SetBinContent(l+1, significance_b);
    }

  // ==> treatment for low-stat X
  // TH1F* hyielderase_b = new TH1F("hyielderase_b", Form(";%s;N_{Signal};", vartitle.c_str()), vc.getnv()-1, vvector.data());
  // TH1F* hyieldbdtgeraseBenr_b = new TH1F("hyieldbdtgeraseBenr_b", Form(";%s;N_{Signal};", vartitle.c_str()), vc.getnv()-1, vvector.data());
  // for(int l=0; l<hyieldbdtg_b->GetNbinsX(); l++)
  //   {
  //     if(bdtg[l] < 0.5) continue;
  //     hyieldbdtgerase_b->SetBinContent(l+1, hyieldbdtg_b->GetBinContent(l+1));
  //     hyieldbdtgerase_b->SetBinError(l+1, hyieldbdtg_b->GetBinError(l+1));
  //     hyieldbdtgeraseBenr_b->SetBinContent(l+1, hyieldbdtgBenr_b->GetBinContent(l+1));
  //     hyieldbdtgeraseBenr_b->SetBinError(l+1, hyieldbdtgBenr_b->GetBinError(l+1));
  //   }

  // ==> fprompt
  TCanvas* cp = new TCanvas("cp", "", 1600, 1200);
  cp->Divide(2, 2);
  for(int l=0; l<vc.getnv()-1; l++)
    {
      vc.hlxymcp_a[l]->Scale(1./vc.hlxymcp_a[l]->Integral(), "width");
      vc.hlxymcp_b[l]->Scale(1./vc.hlxymcp_b[l]->Integral(), "width");
      vc.hlxymcnp_a[l]->Scale(1./vc.hlxymcnp_a[l]->Integral(), "width");
      vc.hlxymcnp_b[l]->Scale(1./vc.hlxymcnp_b[l]->Integral(), "width");
      std::vector<double> vlxyfrac = lxydis::nplxyfrac(vc.hlxymcnp_a[l], vc.hlxymcnp_b[l]);
      hlxyfrac_a->SetBinContent(l+1, vlxyfrac[0]);
      hlxyfrac_a->SetBinError(l+1, vlxyfrac[1]);
      hlxyfrac_b->SetBinContent(l+1, vlxyfrac[2]);
      hlxyfrac_b->SetBinError(l+1, vlxyfrac[3]);

      cp->cd(1);
      vc.hlxymcp_a[l]->Draw((l?"histe same":"histe"));
      cp->cd(2);
      vc.hlxymcp_b[l]->Draw((l?"histe same":"histe"));
      cp->cd(3);
      vc.hlxymcnp_a[l]->Draw((l?"histe same":"histe"));
      cp->cd(4);
      vc.hlxymcnp_b[l]->Draw((l?"histe same":"histe"));
    }
  cp->SaveAs(Form("plots/%s/clxy_vary.pdf", output.c_str()));
  TH1F *hyieldprompt_a, *hyieldprompt_b;
  TEfficiency* grfprompt_a = lxydis::calclxyfprompt(hyield_a, hyieldBenr_a, hlxyfrac_a, "grfprompt_a", &hyieldprompt_a);
  TEfficiency* grfprompt_b = lxydis::calclxyfprompt(hyield_b, hyieldBenr_b, hlxyfrac_b, "grfprompt_b", &hyieldprompt_b);

  // ==> corrected yield
  TH1F* hyieldCorr_a = (TH1F*)hyield_a->Clone("hyieldCorr_a");
  hyieldCorr_a->Divide(vc.heff_a);
  TH1F* hyieldCorr_b = (TH1F*)hyield_b->Clone("hyieldCorr_b");
  hyieldCorr_b->Divide(vc.heff_b);
  TH1F* hyieldPromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldCorr_a");
  hyieldPromptCorr_a->Divide(vc.heff_a);
  TH1F* hyieldPromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldCorr_b");
  hyieldPromptCorr_b->Divide(vc.heff_b);

  // ==> data-MC comparison
  TH1F* hyieldpromptscale_a = (TH1F*)hyieldprompt_a->Clone("hyieldpromptscale_a");
  vc.hmcdisp_a->Scale(1./vc.hmcdisp_a->Integral(), "width");
  hyieldpromptscale_a->Scale(1./hyieldpromptscale_a->Integral(), "width");
  xjjroot::setthgrstyle(vc.hmcdisp_a, kRed+1, 20, 1.2, kRed+1, 2, 2);
  xjjroot::setthgrstyle(hyieldpromptscale_a, kBlack, 25, 1.2, kBlack, 1, 1);
  TH1F* hyieldpromptscale_b = (TH1F*)hyieldprompt_b->Clone("hyieldpromptscale_b");
  vc.hmcdisp_b->Scale(1./vc.hmcdisp_b->Integral(), "width");
  hyieldpromptscale_b->Scale(1./hyieldpromptscale_b->Integral(), "width");
  xjjroot::setthgrstyle(vc.hmcdisp_b, kRed+1, 20, 1.2, kRed+1, 2, 2);
  xjjroot::setthgrstyle(hyieldpromptscale_b, kBlack, 25, 1.2, kBlack, 1, 1);

  TLegend* leg = new TLegend(0.20, 0.70, 0.50, 0.70+0.04*2);
  xjjroot::setleg(leg, 0.038);
  leg->AddEntry(vc.hmcdisp_a, "Prompt MC", "p");
  leg->AddEntry(hyieldpromptscale_a, "Data N_{signal} #times f_{prompt}", "p");

  /*================== Draw ====================*/

  // draw
  int shiftoption = 0 - iscutordis;
  TH2F* hemptyyd = new TH2F("hemptyyd", Form(";%s;Raw N_{Signal} by Fit;", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hyield_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyyd, 0, 0.4);
  TGraphErrors* gr_hyield_a = xjjroot::shifthistcenter(hyield_a, "gr_hyield_a", shiftoption);
  TGraphErrors* gr_hyield_b = xjjroot::shifthistcenter(hyield_b, "gr_hyield_b", shiftoption);
  xjjroot::setthgrstyle(gr_hyield_a, fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyield_b, fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyeff = new TH2F("hemptyeff", Form(";%s;(#alpha #times #epsilon )_{prompt};", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, 0.1);
  xjjroot::sethempty(hemptyeff, 0, 0.4);
  TGraphAsymmErrors* gr_greff_a = xjjroot::shifthistcenter(vc.greff_a, "gr_greff_a", shiftoption);
  TGraphAsymmErrors* gr_greff_b = xjjroot::shifthistcenter(vc.greff_b, "gr_greff_b", shiftoption);
  xjjroot::setthgrstyle(gr_greff_a, fitX::color_a, 34, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_greff_b, fitX::color_b, 34, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  TH2F* hemptysigf = new TH2F("hemptysigf", Form(";%s;S / #sqrt{S+B};", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hsig_b->GetMaximum()*1.4);
  xjjroot::sethempty(hemptysigf, 0, 0.4);
  TGraphErrors* gr_hsigfit_a = xjjroot::shifthistcenter(hsigfit_a, "gr_hsigfit_a", shiftoption);
  TGraphErrors* gr_hsigfit_b = xjjroot::shifthistcenter(hsigfit_b, "gr_hsigfit_b", shiftoption);
  TGraphErrors* gr_hsig_a = xjjroot::shifthistcenter(hsig_a, "gr_hsig_a", shiftoption);
  TGraphErrors* gr_hsig_b = xjjroot::shifthistcenter(hsig_b, "gr_hsig_b", shiftoption);
  xjjroot::setthgrstyle(gr_hsigfit_a, fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hsigfit_b, fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  xjjroot::setthgrstyle(gr_hsig_a, fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hsig_b, fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyBenr = new TH2F("hemptyBenr", Form(";%s;N_{Signal} (l_{xy} > 100#mum) by Fit;", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hyieldBenr_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyBenr, 0, 0.4);
  TGraphErrors* gr_hyieldBenr_a = xjjroot::shifthistcenter(hyieldBenr_a, "gr_hyieldBenr_a", shiftoption);
  TGraphErrors* gr_hyieldBenr_b = xjjroot::shifthistcenter(hyieldBenr_b, "gr_hyieldBenr_b", shiftoption);
  xjjroot::setthgrstyle(gr_hyieldBenr_a, fitX::color_a, 24, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldBenr_b, fitX::color_b, 24, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyydprompt = new TH2F("hemptyydprompt", Form(";%s;N_{Signal} #times f_{prompt};", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hyieldprompt_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyydprompt, 0, 0.4);
  TGraphErrors* gr_hyieldprompt_a = xjjroot::shifthistcenter(hyieldprompt_a, "gr_hyieldprompt_a", shiftoption);
  TGraphErrors* gr_hyieldprompt_b = xjjroot::shifthistcenter(hyieldprompt_b, "gr_hyieldprompt_b", shiftoption);
  xjjroot::setthgrstyle(gr_hyieldprompt_a, fitX::color_a, 21, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldprompt_b, fitX::color_b, 21, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", Form(";%s;f_{prompt} After Cuts;", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, 1.2);
  xjjroot::sethempty(hemptyfprompt, 0, 0.4);
  TGraphAsymmErrors* gr_grfprompt_a = xjjroot::shifthistcenter(grfprompt_a, "gr_grfprompt_a", shiftoption);
  TGraphAsymmErrors* gr_grfprompt_b = xjjroot::shifthistcenter(grfprompt_b, "gr_grfprompt_b", shiftoption);
  xjjroot::setthgrstyle(gr_grfprompt_a, fitX::color_a, 45, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_grfprompt_b, fitX::color_b, 45, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  TH2F* hemptyCorr = new TH2F("hemptyCorr", Form(";%s;N_{Signal} / (#alpha #times #epsilon )_{prompt};", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hyieldCorr_b->GetMaximum()*1.5);
  xjjroot::sethempty(hemptyCorr, 0, 0.4);
  TGraphErrors* gr_hyieldCorr_a = xjjroot::shifthistcenter(hyieldCorr_a, "gr_hyieldCorr_a", shiftoption);
  TGraphErrors* gr_hyieldCorr_b = xjjroot::shifthistcenter(hyieldCorr_b, "gr_hyieldCorr_b", shiftoption);
  xjjroot::setthgrstyle(gr_hyieldCorr_a, fitX::color_a, 3, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldCorr_b, fitX::color_b, 3, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyPromptCorr = new TH2F("hemptyPromptCorr", Form(";%s;N_{Signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt};", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, hyieldPromptCorr_b->GetMaximum()*1.5);
  xjjroot::sethempty(hemptyPromptCorr, 0, 0.4);
  TGraphErrors* gr_hyieldPromptCorr_a = xjjroot::shifthistcenter(hyieldPromptCorr_a, "gr_hyieldPromptCorr_a", shiftoption);
  TGraphErrors* gr_hyieldPromptCorr_b = xjjroot::shifthistcenter(hyieldPromptCorr_b, "gr_hyieldPromptCorr_b", shiftoption);
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_a, fitX::color_a, 5, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_b, fitX::color_b, 5, 1.2, fitX::color_b, 1, 1);
  
  TH2F* hemptyfit_a = new TH2F("hemptyfit_a", ";m_{#mu#mu#pi#pi} (GeV/c);", variation::nmass_a, variation::massmin_a, variation::massmax_a, 10, 0, vc.hdata[0]->GetMaximum());
  xjjroot::sethempty(hemptyfit_a);
  hemptyfit_a->GetXaxis()->SetNdivisions(505);
  TH2F* hemptyfit_b = new TH2F("hemptyfit_b", ";m_{#mu#mu#pi#pi} (GeV/c);", variation::nmass_b, variation::massmin_b, variation::massmax_b, 10, 0, vc.hdata[0]->GetMaximum());
  xjjroot::sethempty(hemptyfit_b);
  hemptyfit_b->GetXaxis()->SetNdivisions(505);
  TH2F* hemptyydpromptscale_a = new TH2F("hemptyydpromptscale_a", Form(";%s;Probability;", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, vc.hmcdisp_a->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyydpromptscale_a, 0, 0.4);
  TH2F* hemptyydpromptscale_b = new TH2F("hemptyydpromptscale_b", Form(";%s;Probability;", vartitle.c_str()), 10, vvector.front()-0.1*(vvector.back()-vvector.front()), vvector.back(), 10, 0, vc.hmcdisp_b->GetMaximum()*1.4);
  xjjroot::sethempty(hemptyydpromptscale_b, 0, 0.4);

  // do draw
  xjjroot::setgstyle(1);

  TCanvas* c2_a = new TCanvas("c2_a", "", 1200, 600);
  c2_a->Divide(2, 1);
  c2_a->cd(1);
  hemptyfit_a->Draw();
  for(auto& hh : vc.hdata) hh->Draw("pe same");
  for(auto& ff : ffun) ff->Draw("l same");
  variation::drawtext("#psi(2S)");
  variation::drawtexlist(iscutordis, vvector, vc.getnv(), vartitle, -0.1, -0.5);
  c2_a->cd(2);
  if(iscutordis)
    {
      hemptyPromptCorr->Draw();
      gr_hyieldPromptCorr_a->Draw("pe same");
    }
  else
    {
      hemptyydpromptscale_a->Draw();
      vc.hmcdisp_a->Draw("histe same");
      hyieldpromptscale_a->Draw("pe same");
      leg->Draw();
    }
  variation::drawtext("#psi(2S)");
  c2_a->SaveAs(Form("plots/%s/cyieldcorr_vary_a.pdf", output.c_str()));

  TCanvas* c2_b = new TCanvas("c2_b", "", 1200, 600);
  c2_b->Divide(2, 1);
  c2_b->cd(1);
  hemptyfit_b->Draw();
  for(auto& hh : vc.hdata) hh->Draw("pe same");
  for(auto& ff : ffun) ff->Draw("l same");
  variation::drawtext("X(3872)");
  variation::drawtexlist(iscutordis, vvector, vc.getnv(), vartitle, -0.1, -0.5);
  c2_b->cd(2);
  if(iscutordis)
    {
      hemptyPromptCorr->Draw();
      gr_hyieldPromptCorr_b->Draw("pe same");
    }
  else
    {
      hemptyydpromptscale_b->Draw();
      vc.hmcdisp_b->Draw("histe same");
      hyieldpromptscale_b->Draw("pe same");
      leg->Draw();
    }
  variation::drawtext("X(3872)");
  c2_b->SaveAs(Form("plots/%s/cyieldcorr_vary_b.pdf", output.c_str()));

  TCanvas* csig = new TCanvas("csig", "", 600, 600);
  hemptysigf->Draw();
  // gr_hsigfit_b->Draw("plX0 same");
  // gr_hsigfit_a->Draw("plX0 same");
  gr_hsig_b->Draw("plX0 same");
  gr_hsig_a->Draw("plX0 same");
  variation::drawalltext();
  csig->SaveAs(Form("plots/%s/csig_vary.pdf", output.c_str()));

  //
  TCanvas* c3 = new TCanvas("c3", "", 2400, 1200);
  c3->Divide(4, 2);
  c3->cd(1);
  hemptyyd->Draw();
  gr_hyield_a->Draw("pe same");
  gr_hyield_b->Draw("pe same");
  variation::drawalltext();
  c3->cd(2);
  hemptyeff->Draw();
  gr_greff_a->Draw("p3e same");
  gr_greff_b->Draw("p3e same");
  variation::drawalltext_simulation();
  c3->cd(3);
  hemptysigf->Draw();
  // gr_hsigfit_b->Draw("plX0 same");
  // gr_hsigfit_a->Draw("plX0 same");
  gr_hsig_b->Draw("plX0 same");
  gr_hsig_a->Draw("plX0 same");
  variation::drawalltext();
  c3->cd(4);
  hemptyCorr->Draw();
  gr_hyieldCorr_a->Draw("pe same");
  gr_hyieldCorr_b->Draw("pe same");
  variation::drawalltext();
  c3->cd(5);
  hemptyBenr->Draw();
  gr_hyieldBenr_a->Draw("pe same");
  gr_hyieldBenr_b->Draw("pe same");
  variation::drawalltext();
  c3->cd(6);
  hemptyydprompt->Draw();
  gr_hyieldprompt_a->Draw("pe same");
  gr_hyieldprompt_b->Draw("pe same");
  variation::drawalltext();
  c3->cd(7);
  hemptyfprompt->Draw();
  gr_grfprompt_a->Draw("pe same");
  gr_grfprompt_b->Draw("pe same");
  variation::drawalltext();
  c3->cd(8);
  hemptyPromptCorr->Draw();
  gr_hyieldPromptCorr_a->Draw("pe same");
  gr_hyieldPromptCorr_b->Draw("pe same");
  variation::drawalltext();
  c3->SaveAs(Form("plots/%s/csig_vary_details.pdf", output.c_str()));

  TFile* outf = new TFile(Form("rootfiles/root_variation_yield_%s.root", output.c_str()), "recreate");
  outf->cd();
  hyield_a->Write();
  hyield_b->Write();
  hbkg_a->Write();
  hbkg_b->Write();
  hsigfit_a->Write();
  hsigfit_b->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==4) { drawvariaiton(argv[1], argv[2], atoi(argv[3])); return 0; }
  return 1;
}
