#include "fitX.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitXvaryfit.h"
#include "fitXvary.h"
#include "lxydis.h"

bool saveornot = true; // save separate fit or not
float setcorrerr(TH1F* hCorr, int ibinsig);
void drawfitXvary(std::string input, std::string output, std::string cutvar)
{
  std::cout<<"\e[32;1m -- "<<__FUNCTION__<<"\e[0m"<<std::endl;

  xjjroot::setgstyle(3);
  fitX::varymva* mvas = fitX::initvarycut(cutvar);
  if(!mvas) { std::cout<<__FUNCTION__<<"error: "<<cutvar<<std::endl; return; }
  fitX::varycut vc(mvas->mva(), cutvar.c_str());
  xjjroot::mkdir(Form("plots/%s/%s/", output.c_str(), cutvar.c_str()));

  TFile* inf = new TFile(input.c_str());
  fitX::init(inf);
  RooWorkspace* ww = (RooWorkspace*)inf->Get("ww");
  // read
  for(int l=0, lp=0; l<mvas->n(); l++)
    {
      vc.hdata[l] = (TH1F*)inf->Get(Form("hdata_%d", l));
      fitX::setmasshist(vc.hdata[l], 0, -0.1);
      vc.hdata[l]->SetMinimum(0);
      vc.hdata[l]->SetMaximum(vc.hdata[l]->GetMaximum()*1.2);
      xjjroot::setthgrstyle(vc.hdata[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.9, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 2);
      vc.hdataBenr[l] = (TH1F*)inf->Get(Form("hdataBenr_%d", l));
      fitX::setmasshist(vc.hdataBenr[l], 0, -0.1);
      vc.hdataBenr[l]->SetMinimum(0);
      vc.hdataBenr[l]->SetMaximum(vc.hdataBenr[l]->GetMaximum()*1.2);
      xjjroot::setthgrstyle(vc.hdataBenr[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.9, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 2);
      vc.hmc_a[l] = (TH1F*)inf->Get(Form("hmc_a_%d", l));
      vc.hmc_b[l] = (TH1F*)inf->Get(Form("hmc_b_%d", l));

      vc.dshdata[l] = (RooDataSet*)ww->data(Form("dshdata_%d", l));
      vc.dshdataBenr[l] = (RooDataSet*)ww->data(Form("dshdataBenr_%d", l));
      vc.dshmc_a[l] = (RooDataSet*)ww->data(Form("dshmc_a_%d", l));
      vc.dshmc_b[l] = (RooDataSet*)ww->data(Form("dshmc_b_%d", l));

      vc.hlxymcnp_a[l] = (TH1F*)inf->Get(Form("hlxymcnp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_a[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.2, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 1);
      vc.hlxymcnp_b[l] = (TH1F*)inf->Get(Form("hlxymcnp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcnp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcnp_b[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.2, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 1);
      vc.hlxymcp_a[l] = (TH1F*)inf->Get(Form("hlxymcp_a_%d", l));
      xjjroot::sethempty(vc.hlxymcp_a[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_a[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.2, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 1);
      vc.hlxymcp_b[l] = (TH1F*)inf->Get(Form("hlxymcp_b_%d", l));
      xjjroot::sethempty(vc.hlxymcp_b[l], 0, -0.1);
      xjjroot::setthgrstyle(vc.hlxymcp_b[l], xjjroot::mycolor_middle[xjjroot::cc[lp]], 20, 0.2, xjjroot::mycolor_middle[xjjroot::cc[lp]], 1, 1);

      if(mvas->ifdraw()[l]) lp++;
    }
  vc.heff_a = (TH1F*)inf->Get("heff_a");
  vc.heff_b = (TH1F*)inf->Get("heff_b");
  vc.greff_a = (TEfficiency*)inf->Get("greff_a");
  vc.greff_b = (TEfficiency*)inf->Get("greff_b");
  vc.hsideband_a = (TH1F*)inf->Get("hsideband_a");
  vc.hsideband_b = (TH1F*)inf->Get("hsideband_b");

  // ==> yield
  TH1F* hyieldmva_a = new TH1F("hyieldmva_a"        , Form(";%s;N_{Signal}", mvas->type().c_str())                 , mvas->n()-1, mvas->mva().data());
  TH1F* hyieldmva_b = new TH1F("hyieldmva_b"        , Form(";%s;N_{Signal}", mvas->type().c_str())                 , mvas->n()-1, mvas->mva().data());
  TH1F* hyieldmvaBenr_a = new TH1F("hyieldmvaBenr_a", Form(";%s;N_{Signal} (l_{xy} > 0.1mm)", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());
  TH1F* hyieldmvaBenr_b = new TH1F("hyieldmvaBenr_b", Form(";%s;N_{Signal} (l_{xy} > 0.1mm)", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());

  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("cvary", "", 700, 600);
  c->cd();
  fitXvary::fitXvaryfit(mvas, "Inclusive", 
                        vc.hdata, vc.hmc_a[0], vc.hmc_b[0],
                        vc.dshdata, vc.dshmc_a[0], vc.dshmc_b[0],
                        hyieldmva_a, hyieldmva_b, Form("%s/%s/%s", output.c_str(), cutvar.c_str(), "inclusive"), c, false); // fixmean = false
  c->SaveAs(Form("plots/%s/%s/cmass_varymva.pdf", output.c_str(), cutvar.c_str()));

  TCanvas* cBenr = new TCanvas("cvaryBenr", "", 700, 600);
  cBenr->cd();
  fitXvary::fitXvaryfit(mvas, "B-enr: l_{xy} > 0.1mm",
                        vc.hdataBenr, vc.hmc_a[0], vc.hmc_b[0],
                        vc.dshdataBenr, vc.dshmc_a[0], vc.dshmc_b[0],
                        hyieldmvaBenr_a, hyieldmvaBenr_b, Form("%s/%s/%s", output.c_str(), cutvar.c_str(), "Benr"), cBenr, true); // fixmean = true
  cBenr->SaveAs(Form("plots/%s/%s/cmass_varymva_Benr.pdf", output.c_str(), cutvar.c_str()));

  // ==> significance
  const int iref = xjjc::findibin(mvas->mva(), mvas->thatval()); // !!
  float yieldnocut_a = hyieldmva_a->GetBinContent(iref+1) / vc.heff_a->GetBinContent(iref+1);
  float yieldnocut_b = hyieldmva_b->GetBinContent(iref+1) / vc.heff_b->GetBinContent(iref+1);
  TH1F* hsigmva_a = new TH1F("hsigmva_a", Form(";%s;S / #sqrt{S+B}", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());  
  TH1F* hsigmva_b = new TH1F("hsigmva_b", Form(";%s;S / #sqrt{S+B}", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());
  float maxsig = 0, xmaxsig = 0; int ibinsig = 0;
  for(int l=0; l<mvas->n()-1; l++)
    {
      float nsig_a = yieldnocut_a * vc.heff_a->GetBinContent(l+1);
      float nsig_b = yieldnocut_b * vc.heff_b->GetBinContent(l+1);
      float nbkg_a = vc.hsideband_a->GetBinContent(l+1) / (vc.masswinH-vc.masswinL) * 0.0036*2;
      float nbkg_b = vc.hsideband_b->GetBinContent(l+1) / (vc.masswinH-vc.masswinL) * 0.0047*2;
      float significance_a = nsig_a / TMath::Sqrt(nsig_a + nbkg_a);
      float significance_b = nsig_b / TMath::Sqrt(nsig_b + nbkg_b);
      hsigmva_a->SetBinContent(l+1, significance_a);
      hsigmva_b->SetBinContent(l+1, significance_b);
      // std::cout<<l<<" "<<nsig_a<<" "<<nbkg_a<<" "<<nsig_b<<" "<<nbkg_b<<std::endl;
      // !!!!!
      if(significance_b > maxsig) { maxsig = significance_b; ibinsig = l+1; xmaxsig = hsigmva_b->GetBinCenter(l+1)-hsigmva_b->GetBinWidth(l+1)/2.; }
      // if(TMath::Abs((hsigmva_b->GetBinCenter(l+1)-hsigmva_b->GetBinWidth(l+1)/2.)-mvas->thatval())<0.01) 
      //   { maxsig = significance_b; ibinsig = l+1; xmaxsig = hsigmva_b->GetBinCenter(l+1)-hsigmva_b->GetBinWidth(l+1)/2.; }
      // !!!!!
    }
  std::cout<<"\e[33;1m"<<xmaxsig<<"\e[0m"<<std::endl;

  // ==> treatment for low-stat X
  TH1F* hyieldmvaerase_b = new TH1F("hyieldmvaerase_b", Form(";%s;N_{Signal}", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());
  TH1F* hyieldmvaeraseBenr_b = new TH1F("hyieldmvaeraseBenr_b", Form(";%s;N_{Signal}", mvas->type().c_str()), mvas->n()-1, mvas->mva().data());
  for(int l=0; l<hyieldmva_b->GetNbinsX(); l++)
    {
      if(mvas->mva()[l] < mvas->minval()) continue;
      hyieldmvaerase_b->SetBinContent(l+1, hyieldmva_b->GetBinContent(l+1));
      hyieldmvaerase_b->SetBinError(l+1, hyieldmva_b->GetBinError(l+1));
      hyieldmvaeraseBenr_b->SetBinContent(l+1, hyieldmvaBenr_b->GetBinContent(l+1));
      hyieldmvaeraseBenr_b->SetBinError(l+1, hyieldmvaBenr_b->GetBinError(l+1));
    }

  // ==> fprompt
  TH1F* hlxyfracmva_a = new TH1F("hlxyfracmva_a" , Form(";%s;", mvas->type().c_str()) , mvas->n()-1, mvas->mva().data());
  TH1F* hlxyfracmva_b = new TH1F("hlxyfracmva_b" , Form(";%s;", mvas->type().c_str()) , mvas->n()-1, mvas->mva().data());
  TCanvas* cp = new TCanvas("cp", "", 1700, 1200);
  cp->Divide(2, 2);
  for(int l=0; l<mvas->n()-1; l++)
    {
      vc.hlxymcp_a[l]->Scale(1./vc.hlxymcp_a[l]->Integral(), "width");
      vc.hlxymcp_b[l]->Scale(1./vc.hlxymcp_b[l]->Integral(), "width");
      vc.hlxymcnp_a[l]->Scale(1./vc.hlxymcnp_a[l]->Integral(), "width");
      vc.hlxymcnp_b[l]->Scale(1./vc.hlxymcnp_b[l]->Integral(), "width");
      std::vector<double> vlxyfrac = lxydis::nplxyfrac(vc.hlxymcnp_a[l], vc.hlxymcnp_b[l]);
      hlxyfracmva_a->SetBinContent(l+1, vlxyfrac[0]);
      hlxyfracmva_a->SetBinError(l+1, vlxyfrac[1]);
      hlxyfracmva_b->SetBinContent(l+1, vlxyfrac[2]);
      hlxyfracmva_b->SetBinError(l+1, vlxyfrac[3]);

      if(mvas->ifdraw()[l])
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
  cp->SaveAs(Form("plots/%s/%s/clxy_varymva.pdf", output.c_str(), cutvar.c_str()));
  TH1F *hyieldprompt_a, *hyieldprompt_b;
  TEfficiency* grfprompt_a = lxydis::calclxyfprompt(hyieldmva_a, hyieldmvaBenr_a, hlxyfracmva_a, "grfprompt_a", &hyieldprompt_a);
  TEfficiency* grfprompt_b = lxydis::calclxyfprompt(hyieldmvaerase_b, hyieldmvaeraseBenr_b, hlxyfracmva_b, "grfprompt_b", &hyieldprompt_b);

  // ==> corrected yield
  TH1F* hyieldCorr_a = (TH1F*)hyieldmva_a->Clone("hyieldCorr_a");
  hyieldCorr_a->Divide(vc.heff_a);
  TH1F* hyieldCorr_b = (TH1F*)hyieldmvaerase_b->Clone("hyieldCorr_b");
  hyieldCorr_b->Divide(vc.heff_b);
  TH1F* hyieldPromptCorr_a = (TH1F*)hyieldprompt_a->Clone("hyieldPromptCorr_a");
  hyieldPromptCorr_a->Divide(vc.heff_a);
  TH1F* hyieldPromptCorr_b = (TH1F*)hyieldprompt_b->Clone("hyieldPromptCorr_b");
  hyieldPromptCorr_b->Divide(vc.heff_b);

  float errhyieldCorr_a = setcorrerr(hyieldCorr_a, ibinsig);
  float errhyieldCorr_b = setcorrerr(hyieldCorr_b, ibinsig);
  float errhyieldPromptCorr_a = setcorrerr(hyieldPromptCorr_a, ibinsig);
  float errhyieldPromptCorr_b = setcorrerr(hyieldPromptCorr_b, ibinsig);

  auto bPromptCorr_a = xjjroot::drawbox(mvas->binmin(),  hyieldPromptCorr_a->GetBinContent(ibinsig)-errhyieldPromptCorr_a, mvas->binmax(), hyieldPromptCorr_a->GetBinContent(ibinsig)+errhyieldPromptCorr_a, fitX::color_a, 0.1, 1001);
  auto lPromptCorr_a = xjjroot::drawline(mvas->binmin(), hyieldPromptCorr_a->GetBinContent(ibinsig),                       mvas->binmax(), hyieldPromptCorr_a->GetBinContent(ibinsig),                       fitX::color_a, 2,   1, 0.8);
  auto bPromptCorr_b = xjjroot::drawbox(mvas->binmin(),  hyieldPromptCorr_b->GetBinContent(ibinsig)-errhyieldPromptCorr_b, mvas->binmax(), hyieldPromptCorr_b->GetBinContent(ibinsig)+errhyieldPromptCorr_b, fitX::color_b, 0.1, 1001);
  auto lPromptCorr_b = xjjroot::drawline(mvas->binmin(), hyieldPromptCorr_b->GetBinContent(ibinsig),                       mvas->binmax(), hyieldPromptCorr_b->GetBinContent(ibinsig),                       fitX::color_b, 2,   1, 0.8);
  auto bCorr_a       = xjjroot::drawbox(mvas->binmin(),  hyieldCorr_a->GetBinContent(ibinsig)-errhyieldCorr_a,             mvas->binmax(), hyieldCorr_a->GetBinContent(ibinsig)+errhyieldCorr_a,             fitX::color_a, 0.1, 1001);
  auto lCorr_a       = xjjroot::drawline(mvas->binmin(), hyieldCorr_a->GetBinContent(ibinsig),                             mvas->binmax(), hyieldCorr_a->GetBinContent(ibinsig),                             fitX::color_a, 2,   1, 0.8);
  auto bCorr_b       = xjjroot::drawbox(mvas->binmin(),  hyieldCorr_b->GetBinContent(ibinsig)-errhyieldCorr_b,             mvas->binmax(), hyieldCorr_b->GetBinContent(ibinsig)+errhyieldCorr_b,             fitX::color_b, 0.1, 1001);
  auto lCorr_b       = xjjroot::drawline(mvas->binmin(), hyieldCorr_b->GetBinContent(ibinsig),                             mvas->binmax(), hyieldCorr_b->GetBinContent(ibinsig),                             fitX::color_b, 2,   1, 0.8);

  /*================== Draw ====================*/

  // draw
  // --> 1
  TH2F* hemptyyd = new TH2F("hemptyyd", Form(";%s >;Raw N_{Signal} by Fit", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hyieldmva_a->GetMaximum(), hyieldmva_b->GetMaximum())*1.4);
  xjjroot::sethempty(hemptyyd, 0, 0.4); hemptyyd->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hyieldmva_a = xjjroot::shifthistcenter(hyieldmva_a, "gr_hyieldmva_a");
  TGraphErrors* gr_hyieldmva_b = xjjroot::shifthistcenter(hyieldmva_b, "gr_hyieldmva_b");
  xjjroot::setthgrstyle(gr_hyieldmva_a, fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldmva_b, fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);
  // --> 2
  TH2F* hemptyeff = new TH2F("hemptyeff", Form(";%s >;(#alpha #times #epsilon )_{prompt}", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, 0.1);
  xjjroot::sethempty(hemptyeff, 0, 0.4); hemptyeff->GetXaxis()->SetNdivisions(505);
  TGraphAsymmErrors* gr_greff_a = xjjroot::shifthistcenter(vc.greff_a, "gr_greff_a");
  TGraphAsymmErrors* gr_greff_b = xjjroot::shifthistcenter(vc.greff_b, "gr_greff_b");
  xjjroot::setthgrstyle(gr_greff_a, fitX::color_a, 34, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_greff_b, fitX::color_b, 34, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  // --> 3
  TH2F* hemptysigf = new TH2F("hemptysigf", Form(";%s >;S / #sqrt{S+B}", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hsigmva_a->GetMaximum(), hsigmva_b->GetMaximum())*1.3);
  xjjroot::sethempty(hemptysigf, 0, 0.4); hemptysigf->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hsigmva_a = xjjroot::shifthistcenter(hsigmva_a, "gr_hsigmva_a");
  TGraphErrors* gr_hsigmva_b = xjjroot::shifthistcenter(hsigmva_b, "gr_hsigmva_b");
  xjjroot::setthgrstyle(gr_hsigmva_a, fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hsigmva_b, fitX::color_b, 46, 1.2, fitX::color_b, 1, 1); 
  // --> 5
  TH2F* hemptyBenr = new TH2F("hemptyBenr", Form(";%s >;N_{Signal} (l_{xy} > 100#mum) by Fit", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hyieldmvaBenr_a->GetMaximum(), hyieldmvaBenr_b->GetMaximum())*1.4);
  xjjroot::sethempty(hemptyBenr, 0, 0.4); hemptyBenr->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hyieldmvaBenr_a = xjjroot::shifthistcenter(hyieldmvaBenr_a, "gr_hyieldmvaBenr_a");
  TGraphErrors* gr_hyieldmvaBenr_b = xjjroot::shifthistcenter(hyieldmvaBenr_b, "gr_hyieldmvaBenr_b");
  xjjroot::setthgrstyle(gr_hyieldmvaBenr_a, fitX::color_a, 24, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldmvaBenr_b, fitX::color_b, 24, 1.2, fitX::color_b, 1, 1);
  // --> 6
  TH2F* hemptyydprompt = new TH2F("hemptyydprompt", Form(";%s >;N_{Signal} #times f_{prompt}", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hyieldprompt_a->GetMaximum(), hyieldprompt_b->GetMaximum())*1.4);
  xjjroot::sethempty(hemptyydprompt, 0, 0.4); hemptyydprompt->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hyieldprompt_a = xjjroot::shifthistcenter(hyieldprompt_a, "gr_hyieldprompt_a");
  TGraphErrors* gr_hyieldprompt_b = xjjroot::shifthistcenter(hyieldprompt_b, "gr_hyieldprompt_b");
  xjjroot::setthgrstyle(gr_hyieldprompt_a, fitX::color_a, 21, 1.2, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldprompt_b, fitX::color_b, 21, 1.2, fitX::color_b, 1, 1);
  // --> 7
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", Form(";%s >;f_{prompt} After Cuts", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, 1.2);
  xjjroot::sethempty(hemptyfprompt, 0, 0.4); hemptyfprompt->GetXaxis()->SetNdivisions(505);
  TGraphAsymmErrors* gr_grfprompt_a = xjjroot::shifthistcenter(grfprompt_a, "gr_grfprompt_a");
  TGraphAsymmErrors* gr_grfprompt_b = xjjroot::shifthistcenter(grfprompt_b, "gr_grfprompt_b");
  xjjroot::setthgrstyle(gr_grfprompt_a, fitX::color_a, 45, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(gr_grfprompt_b, fitX::color_b, 45, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);
  // --> 4
  TH2F* hemptyCorr = new TH2F("hemptyCorr", Form(";%s >;N_{Signal} / (#alpha #times #epsilon )_{prompt}", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hyieldCorr_a->GetMaximum(), hyieldCorr_b->GetMaximum())*1.5);
  xjjroot::sethempty(hemptyCorr, 0, 0.4); hemptyCorr->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hyieldCorr_a = xjjroot::shifthistcenter(hyieldCorr_a, "gr_hyieldCorr_a");
  TGraphErrors* gr_hyieldCorr_b = xjjroot::shifthistcenter(hyieldCorr_b, "gr_hyieldCorr_b");
  xjjroot::setthgrstyle(gr_hyieldCorr_a, fitX::color_a, 47, 1.3, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldCorr_b, fitX::color_b, 47, 1.3, fitX::color_b, 1, 1);
  // --> 8
  TH2F* hemptyPromptCorr = new TH2F("hemptyPromptCorr", Form(";%s >;N_{Signal} #times f_{prompt} / (#alpha #times #epsilon )_{prompt}", mvas->type().c_str()), 10, mvas->binmin(), mvas->binmax(), 10, 0, std::max(hyieldPromptCorr_a->GetMaximum(), hyieldPromptCorr_b->GetMaximum())*1.5);
  xjjroot::sethempty(hemptyPromptCorr, 0, 0.4); hemptyPromptCorr->GetXaxis()->SetNdivisions(505);
  TGraphErrors* gr_hyieldPromptCorr_a = xjjroot::shifthistcenter(hyieldPromptCorr_a, "gr_hyieldPromptCorr_a");
  TGraphErrors* gr_hyieldPromptCorr_b = xjjroot::shifthistcenter(hyieldPromptCorr_b, "gr_hyieldPromptCorr_b");
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_a, fitX::color_a, 5, 1.3, fitX::color_a, 1, 1);
  xjjroot::setthgrstyle(gr_hyieldPromptCorr_b, fitX::color_b, 5, 1.3, fitX::color_b, 1, 1);

  // do draw
  xjjroot::setgstyle(1);
  TCanvas* c8 = new TCanvas("c8", "", 2400, 1200);
  c8->Divide(4, 2);
  c8->cd(1);
  hemptyyd->Draw();
  gr_hyieldmva_a->Draw("pe same");
  gr_hyieldmva_b->Draw("pe same");
  drawalltext();
  c8->cd(2);
  hemptyeff->Draw();
  gr_greff_a->Draw("p3e same");
  gr_greff_b->Draw("p3e same");
  drawalltext();
  c8->cd(3);
  hemptysigf->Draw();
  gr_hsigmva_b->Draw("plX0 same");
  gr_hsigmva_a->Draw("plX0 same");
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hsigmva_a->GetMaximum(), hsigmva_b->GetMaximum())*1.3, kGray+1, 2, 1);
  drawalltext();
  c8->cd(4);
  hemptyCorr->Draw();
  bCorr_a->Draw(); lCorr_a->Draw(); bCorr_b->Draw(); lCorr_b->Draw();
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hyieldCorr_a->GetMaximum(), hyieldCorr_b->GetMaximum())*1.5, kGray+1, 2, 1);
  gr_hyieldCorr_a->Draw("pe same");
  gr_hyieldCorr_b->Draw("pe same");
  drawalltext();
  c8->cd(5);
  hemptyBenr->Draw();
  gr_hyieldmvaBenr_a->Draw("pe same");
  gr_hyieldmvaBenr_b->Draw("pe same");
  drawalltext();
  c8->cd(6);
  hemptyydprompt->Draw();
  gr_hyieldprompt_a->Draw("pe same");
  gr_hyieldprompt_b->Draw("pe same");
  drawalltext();
  c8->cd(7);
  hemptyfprompt->Draw();
  gr_grfprompt_a->Draw("pe same");
  gr_grfprompt_b->Draw("pe same");
  drawalltext();
  c8->cd(8);
  hemptyPromptCorr->Draw();
  bPromptCorr_a->Draw(); lPromptCorr_a->Draw(); bPromptCorr_b->Draw(); lPromptCorr_b->Draw();
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hyieldPromptCorr_a->GetMaximum(), hyieldPromptCorr_b->GetMaximum())*1.5, kGray+1, 2, 1);
  gr_hyieldPromptCorr_a->Draw("pe same");
  gr_hyieldPromptCorr_b->Draw("pe same");
  drawalltext();
  c8->SaveAs(Form("plots/%s/%s/c8detail_varymva.pdf", output.c_str(), cutvar.c_str()));

  xjjroot::setgstyle(3);
  TCanvas* c83 = new TCanvas("c83", "", 600, 600);
  hemptysigf->Draw();
  gr_hsigmva_b->SetLineWidth(2);
  gr_hsigmva_a->SetLineWidth(2);
  gr_hsigmva_b->Draw("plX0 same");
  gr_hsigmva_a->Draw("plX0 same");
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hsigmva_a->GetMaximum(), hsigmva_b->GetMaximum())*1.3, kGray+1, 2, 2);
  drawalltext();
  c83->SaveAs(Form("plots/%s/%s/c83_varymva.pdf", output.c_str(), cutvar.c_str()));
  TCanvas* c84 = new TCanvas("c84", "", 600, 600);
  hemptyCorr->Draw();
  bCorr_a->Draw(); lCorr_a->Draw(); bCorr_b->Draw(); lCorr_b->Draw();
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hyieldCorr_a->GetMaximum(), hyieldCorr_b->GetMaximum())*1.5, kGray+1, 2, 1);
  gr_hyieldCorr_a->Draw("pe same");
  gr_hyieldCorr_b->Draw("pe same");
  drawalltext();
  c84->SaveAs(Form("plots/%s/%s/c84_varymva.pdf", output.c_str(), cutvar.c_str()));
  TCanvas* c88 = new TCanvas("c88", "", 600, 600);
  hemptyPromptCorr->Draw();
  bPromptCorr_a->Draw(); lPromptCorr_a->Draw(); bPromptCorr_b->Draw(); lPromptCorr_b->Draw();
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hyieldPromptCorr_a->GetMaximum(), hyieldPromptCorr_b->GetMaximum())*1.5, kGray+1, 2, 1);
  gr_hyieldPromptCorr_a->Draw("pe same");
  gr_hyieldPromptCorr_b->Draw("pe same");
  drawalltext();
  c88->SaveAs(Form("plots/%s/%s/c88_varymva.pdf", output.c_str(), cutvar.c_str()));

  xjjroot::setgstyle(1);
  TCanvas* c8124 = new TCanvas("c8124", "", 1800, 600);
  c8124->Divide(3, 1);
  c8124->cd(1);
  hemptyyd->Draw();
  gr_hyieldmva_a->Draw("pe same");
  gr_hyieldmva_b->Draw("pe same");
  drawalltext();
  c8124->cd(2);
  hemptyeff->Draw();
  gr_greff_a->Draw("p3e same");
  gr_greff_b->Draw("p3e same");
  drawalltext();
  c8124->cd(3);
  hemptyCorr->Draw();
  bCorr_a->Draw(); lCorr_a->Draw(); bCorr_b->Draw(); lCorr_b->Draw();
  xjjroot::drawline(xmaxsig, 0, xmaxsig, std::max(hyieldCorr_a->GetMaximum(), hyieldCorr_b->GetMaximum())*1.5, kGray+1, 2, 1);
  gr_hyieldCorr_a->Draw("pe same");
  gr_hyieldCorr_b->Draw("pe same");
  drawalltext();
  c8124->SaveAs(Form("plots/%s/%s/c8124_varymva.pdf", output.c_str(), cutvar.c_str()));

  std::string outputname = Form("rootfiles/%s%s/%s/root_fitXvary_yields.root", output.c_str(), fitX::tagname().c_str(), cutvar.c_str());
  xjjroot::mkdir(outputname);
  TFile* outf = new TFile(outputname.c_str(), "recreate");
  outf->cd();
  hyieldmva_a->Write();
  hyieldmva_b->Write();
  outf->Close();

  // print
  float nominal_a = hyieldCorr_a->GetBinContent(ibinsig);
  float per_a = 0;
  float nominal_b = hyieldCorr_b->GetBinContent(ibinsig);
  float per_b = 0;
  for(int i=0; i<hyieldCorr_a->GetNbinsX(); i++)
    {
      float dev;
      dev = TMath::Abs(hyieldCorr_a->GetBinContent(i+1)-nominal_a)/nominal_a;
      if(dev > per_a && hyieldCorr_a->GetBinContent(i+1)>0) per_a = dev;
      dev = TMath::Abs(hyieldCorr_b->GetBinContent(i+1)-nominal_b)/nominal_b;
      if(dev > per_b && hyieldCorr_b->GetBinContent(i+1)>0) per_b = dev;
    }
  float per_ab = TMath::Sqrt(per_a*per_a + per_b*per_b);
  std::cout<<"Efficiency & "<<Form("%.1f", per_a*1.e+2)<<"\\% & "<<Form("%.1f", per_b*1.e+2)<<"\\% & "<<Form("%.1f", per_ab*1.e+2)<<"\\% \\\\"<<std::endl;

}

int main(int argc, char* argv[])
{
  if(argc==4) { drawfitXvary(argv[1], argv[2], argv[3]); return 0; }
  return 1;
}

float setcorrerr(TH1F* hCorr, int ibinsig)
{
  float nominalerr = hCorr->GetBinError(ibinsig);
  for(int i=0; i<hCorr->GetXaxis()->GetNbins(); i++)
    {
      hCorr->SetBinError(i+1, TMath::Sqrt(fabs(nominalerr*nominalerr -
                                               hCorr->GetBinError(i+1)*hCorr->GetBinError(i+1)))
                         );
    }
  return nominalerr;
}

