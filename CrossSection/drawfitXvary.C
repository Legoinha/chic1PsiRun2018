#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <string>

#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitXvary.h"
#include "lxydis.h"

void drawfitXvary(std::string output)
{
  xjjroot::setgstyle(3);
  TFile* inf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()));
  std::vector<TH1F*> hdata(nbdtg*ndls), hdatagt(nbdtg*ndls), hmc_a(nbdtg*ndls), hmc_b(nbdtg*ndls);
  std::vector<TF1*> ff(nbdtg*ndls, 0), ffgt(nbdtg*ndls, 0);
  std::vector<TH1F*> heff_a(ndls), heff_b(ndls), hlxymcnp_a(nbdtg*ndls), hlxymcnp_b(nbdtg*ndls);
  std::vector<TEfficiency*> greff_a(ndls), greff_b(ndls);

  for(int k=0; k<ndls; k++)
    {
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx] = (TH1F*)inf->Get(Form("hdata_%d_%d", k, l));
          fitX::setmasshist(hdata[idx], 0, -0.1);
          hdata[idx]->SetMinimum(0);
          xjjroot::setthgrstyle(hdata[idx], xjjroot::colorlist_dark[lp], 20, 0.9, xjjroot::colorlist_dark[lp], 1, 1);
          hdatagt[idx] = (TH1F*)inf->Get(Form("hdatagt_%d_%d", k, l));
          fitX::setmasshist(hdatagt[idx], 0, -0.1);
          hdatagt[idx]->SetMinimum(0);
          xjjroot::setthgrstyle(hdatagt[idx], xjjroot::colorlist_dark[lp], 20, 0.9, xjjroot::colorlist_dark[lp], 1, 1);
          hmc_a[idx] = (TH1F*)inf->Get(Form("hmc_a_%d_%d", k, l));
          hmc_b[idx] = (TH1F*)inf->Get(Form("hmc_b_%d_%d", k, l));
          if(pbdtg[l]) lp++;
          hlxymcnp_a[idx] = (TH1F*)inf->Get(Form("hlxymcnp_a_%d_%d", k, l));
          hlxymcnp_b[idx] = (TH1F*)inf->Get(Form("hlxymcnp_b_%d_%d", k, l));
        }
      heff_a[k] = (TH1F*)inf->Get(Form("heff_a_%d", k));
      heff_b[k] = (TH1F*)inf->Get(Form("heff_b_%d", k));
      greff_a[k] = (TEfficiency*)inf->Get(Form("greff_a_%d", k));
      greff_b[k] = (TEfficiency*)inf->Get(Form("greff_b_%d", k));
    }

  std::vector<TH1F*> hyieldbdtgL(ndls), hyieldbdtgH(ndls), hbkgbdtgL(ndls), hbkgbdtgH(ndls), hsigfbdtgL(ndls), hsigfbdtgH(ndls);
  std::vector<TH1F*> hyieldbdtgBenrL(ndls), hyieldbdtgBenrH(ndls), hbkgbdtgBenrL(ndls), hbkgbdtgBenrH(ndls), hsigfbdtgBenrL(ndls), hsigfbdtgBenrH(ndls);
  std::vector<TH1F*> hlxyfracL(ndls), hlxyfracH(ndls);
  for(int k=0; k<ndls; k++)
    {
      hyieldbdtgL[k] = new TH1F(Form("hyieldbdtgL_%d", k), ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
      hyieldbdtgH[k] = new TH1F(Form("hyieldbdtgH_%d", k), ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
      hbkgbdtgL[k] = new TH1F(Form("hbkgbdtgL_%d", k), ";BDTG;N_{Background}", bdtg.size()-1, bdtg.data());
      hbkgbdtgH[k] = new TH1F(Form("hbkgbdtgH_%d", k), ";BDTG;N_{Background}", bdtg.size()-1, bdtg.data());
      hsigfbdtgL[k] = new TH1F(Form("hsigfbdtgL_%d", k), ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());
      hsigfbdtgH[k] = new TH1F(Form("hsigfbdtgH_%d", k), ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());
      hyieldbdtgBenrL[k] = new TH1F(Form("hyieldbdtgBenrL_%d", k), ";BDTG;N_{Signal} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hyieldbdtgBenrH[k] = new TH1F(Form("hyieldbdtgBenrH_%d", k), ";BDTG;N_{Signal} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hbkgbdtgBenrL[k] = new TH1F(Form("hbkgbdtgBenrL_%d", k), ";BDTG;N_{Background} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hbkgbdtgBenrH[k] = new TH1F(Form("hbkgbdtgBenrH_%d", k), ";BDTG;N_{Background} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hsigfbdtgBenrL[k] = new TH1F(Form("hsigfbdtgBenrL_%d", k), ";BDTG;S / #sqrt{S+B} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hsigfbdtgBenrH[k] = new TH1F(Form("hsigfbdtgBenrH_%d", k), ";BDTG;S / #sqrt{S+B} (l_{xy} > 0.1mm)", bdtg.size()-1, bdtg.data());
      hlxyfracL[k] = new TH1F(Form("hlxyfracL_%d", k), ";BDTG;", bdtg.size()-1, bdtg.data());
      hlxyfracH[k] = new TH1F(Form("hlxyfracH_%d", k), ";BDTG;", bdtg.size()-1, bdtg.data());

      TCanvas* c = new TCanvas("cvary", "", 800, 1200);
      c->Divide(1, 2);
      c->cd(1);
      hdata[k*nbdtg]->Draw("pe");
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          c->cd(1);
          if(pbdtg[l]) hdata[idx]->Draw("pe same");
          std::vector<TF1*> funs;
          bool saveornot = k==0;
          funs = fitX::fit(xjjroot::copyobject(hdata[idx], Form("hdata_%d_%d", k, l)), 0, hmc_a[0], hmc_b[0], Form("idx/hdata_%d_%d", k, l), true, saveornot); // fix mean
          ff[idx] = funs[0];
          ff[idx]->SetName(Form("ff_%d_%d", k, l));
          c->cd(1);
          xjjroot::settfstyle(ff[idx], xjjroot::colorlist_dark[lp], 2, 2);
          if(pbdtg[l]) 
            {
              ff[idx]->Draw("l same");
              xjjroot::drawtex(0.70+0.08*(lp%3), 0.75-0.04*(lp/3), Form("> %.2f%s", bdtg[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[lp]);
              lp++;
            }
          if(l<nbdtg-1)
            {
              float ysigL = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigL < 0) ysigL = 0;
              float ysigLerr = funs[0]->GetParError(5)*ysigL/funs[0]->GetParameter(5);
              float ysigH = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigH < 0) ysigH = 0;
              float ysigHerr = funs[0]->GetParError(10)*ysigH/funs[0]->GetParameter(10);
              float ybkgL = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH; if(ybkgL < 0) ybkgL = 0;
              float ybkgH = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH; if(ybkgH < 0) ybkgH = 0;
              float ysigfL = (ysigL+ybkgL>0)?(ysigL/TMath::Sqrt(ysigL+ybkgL)):0;
              float ysigfH = (ysigH+ybkgH>0)?(ysigH/TMath::Sqrt(ysigH+ybkgH)):0;
              hyieldbdtgL[k]->SetBinContent(l+1, ysigL);
              hyieldbdtgL[k]->SetBinError(l+1, ysigLerr);
              hyieldbdtgH[k]->SetBinContent(l+1, ysigH);
              hyieldbdtgH[k]->SetBinError(l+1, ysigHerr);
              hbkgbdtgL[k]->SetBinContent(l+1, ybkgL);
              hbkgbdtgH[k]->SetBinContent(l+1, ybkgH);
              hsigfbdtgL[k]->SetBinContent(l+1, ysigfL);
              hsigfbdtgH[k]->SetBinContent(l+1, ysigfH);
            }
        }
      xjjroot::setgstyle(1);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
      c->cd(2);
      hdatagt[k*nbdtg]->Draw("pe");
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          c->cd(2);
          if(pbdtg[l]) hdatagt[idx]->Draw("pe same");
          std::vector<TF1*> funs;
          bool saveornot = k==0;
          funs = fitX::fit(xjjroot::copyobject(hdatagt[idx], Form("hdatagt_%d_%d", k, l)), 0, hmc_a[0], hmc_b[0], Form("idx/hdatagt_%d_%d", k, l), true, saveornot); // fix mean
          ffgt[idx] = funs[0];
          ffgt[idx]->SetName(Form("ffgt_%d_%d", k, l));
          c->cd(2);
          xjjroot::settfstyle(ffgt[idx], xjjroot::colorlist_dark[lp], 2, 2);
          if(pbdtg[l]) 
            {
              ffgt[idx]->Draw("l same");
              xjjroot::drawtex(0.70+0.08*(lp%3), 0.75-0.04*(lp/3), Form("> %.2f%s", bdtg[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[lp]);
              lp++;
            }
          if(l<nbdtg-1)
            {
              float ysigL = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigL < 0) ysigL = 0;
              float ysigLerr = funs[0]->GetParError(5)*ysigL/funs[0]->GetParameter(5);
              float ysigH = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigH < 0) ysigH = 0;
              float ysigHerr = funs[0]->GetParError(10)*ysigH/funs[0]->GetParameter(10);
              float ybkgL = funs[3]->Integral(MASS_PSI2S - mytmva::sigwindowL, MASS_PSI2S + mytmva::sigwindowL) / fitX::BIN_WIDTH; if(ybkgL < 0) ybkgL = 0;
              float ybkgH = funs[3]->Integral(MASS_X     - mytmva::sigwindowH, MASS_X     + mytmva::sigwindowH) / fitX::BIN_WIDTH; if(ybkgH < 0) ybkgH = 0;
              float ysigfL = (ysigL+ybkgL>0)?(ysigL/TMath::Sqrt(ysigL+ybkgL)):0;
              float ysigfH = (ysigH+ybkgH>0)?(ysigH/TMath::Sqrt(ysigH+ybkgH)):0;
              hyieldbdtgBenrL[k]->SetBinContent(l+1, ysigL);
              hyieldbdtgBenrL[k]->SetBinError(l+1, ysigLerr);
              hyieldbdtgBenrH[k]->SetBinContent(l+1, ysigH);
              hyieldbdtgBenrH[k]->SetBinError(l+1, ysigHerr);
              hbkgbdtgBenrL[k]->SetBinContent(l+1, ybkgL);
              hbkgbdtgBenrH[k]->SetBinContent(l+1, ybkgH);
              hsigfbdtgBenrL[k]->SetBinContent(l+1, ysigfL);
              hsigfbdtgBenrH[k]->SetBinContent(l+1, ysigfH);
            }
        }
      c->cd(2);
      xjjroot::setgstyle(1);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();

      c->SaveAs(Form("plots/cmass_varybdtg_%s_%d.pdf", output.c_str(), k));

      TCanvas* c2 = new TCanvas("cvary2", "", 800, 600);
      c2->cd();
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          std::string opt = lp==0?"pe":"pe same";
          if(!pbdtg[l]) continue;
          if(lp >= 3) 
            {
              int idx = k*nbdtg + l;
              hdata[idx]->Draw(opt.c_str());
              ffgt[idx]->Draw("l same");
              xjjroot::drawtex(0.70+0.08*(lp%3), 0.75-0.04*(lp/3-1), Form("> %.2f%s", bdtg[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[lp]);
            }
          lp++;
        }
      xjjroot::setgstyle(3);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
      c2->SaveAs(Form("plots/cmass_varybdtg_%s_%d_zoomin.pdf", output.c_str(), k));

      delete c;
      delete c2;

      for(int l=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hlxymcnp_a[idx]->Scale(1./hlxymcnp_a[idx]->Integral(), "width");
          hlxymcnp_b[idx]->Scale(1./hlxymcnp_b[idx]->Integral(), "width");
          std::vector<double> vlxyfrac = lxydis::nplxyfrac(hlxymcnp_a[idx], hlxymcnp_b[idx]);
          std::cout<<vlxyfrac[0]<<" "<<vlxyfrac[2]<<std::endl;
          hlxyfracL[k]->SetBinContent(l+1, vlxyfrac[0]);
          hlxyfracL[k]->SetBinError(l+1, vlxyfrac[1]);
          hlxyfracH[k]->SetBinContent(l+1, vlxyfrac[2]);
          hlxyfracH[k]->SetBinError(l+1, vlxyfrac[3]);
        }

    }

  const int idls = 0;
  xjjroot::sethempty(hyieldbdtgL[idls], 0, 0.3);
  xjjroot::setthgrstyle(hyieldbdtgL[idls], fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::sethempty(hyieldbdtgH[idls], 0, 0.3);
  xjjroot::setthgrstyle(hyieldbdtgH[idls], fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);
  xjjroot::sethempty(hsigfbdtgL[idls], 0, 0.3);
  xjjroot::setthgrstyle(hsigfbdtgL[idls], fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::sethempty(hsigfbdtgH[idls], 0, 0.3);
  xjjroot::setthgrstyle(hsigfbdtgH[idls], fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  TH2F* hemptyeff = new TH2F("hemptyeff", ";BDTG;#alpha #times #epsilon", 10, -1, 1, 10, 0, 0.1);
  xjjroot::sethempty(hemptyeff, 0, 0.3);
  xjjroot::setthgrstyle(greff_a[idls], fitX::color_a, 34, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(greff_b[idls], fitX::color_b, 34, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);

  xjjroot::sethempty(hyieldbdtgBenrL[idls], 0, 0.3);
  xjjroot::setthgrstyle(hyieldbdtgBenrL[idls], fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::sethempty(hyieldbdtgBenrH[idls], 0, 0.3);
  xjjroot::setthgrstyle(hyieldbdtgBenrH[idls], fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);

  TH1F* hyieldbdtgeraseH = new TH1F("hyieldbdtgeraseH", ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
  TH1F* hyieldbdtgeraseBenrH = new TH1F("hyieldbdtgeraseBenrH", ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
  for(int l=0; l<hyieldbdtgH[idls]->GetNbinsX(); l++)
    {
      if(bdtg[l] < 0.5) continue;
      hyieldbdtgeraseH->SetBinContent(l+1, hyieldbdtgH[idls]->GetBinContent(l+1));
      hyieldbdtgeraseH->SetBinError(l+1, hyieldbdtgH[idls]->GetBinError(l+1));
      hyieldbdtgeraseBenrH->SetBinContent(l+1, hyieldbdtgBenrH[idls]->GetBinContent(l+1));
      hyieldbdtgeraseBenrH->SetBinError(l+1, hyieldbdtgBenrH[idls]->GetBinError(l+1));
    }

  TEfficiency* grfpromptL = lxydis::calclxyfprompt(hyieldbdtgL[idls], hyieldbdtgBenrL[idls], hlxyfracL[idls], "grfpromptL");
  TEfficiency* grfpromptH = lxydis::calclxyfprompt(hyieldbdtgeraseH, hyieldbdtgeraseBenrH, hlxyfracH[idls], "grfpromptH");
  TH2F* hemptyfprompt = new TH2F("hemptyfprompt", ";BDTG;f_{prompt} After Cuts", 10, -1, 1, 10, 0, 1.2);
  xjjroot::sethempty(hemptyfprompt, 0, 0.3);
  xjjroot::setthgrstyle(grfpromptL, fitX::color_a, 45, 1.2, fitX::color_a, 1, 1, fitX::color_a, 0.2, 1001);
  xjjroot::setthgrstyle(grfpromptH, fitX::color_b, 45, 1.2, fitX::color_b, 1, 1, fitX::color_b, 0.2, 1001);

  xjjroot::setgstyle(1);
  TCanvas* c3 = new TCanvas("c3", "", 1800, 1200);
  c3->Divide(3, 2);
  c3->cd(1);
  hyieldbdtgL[idls]->SetMaximum(hyieldbdtgL[idls]->GetMaximum()*1.5);
  hyieldbdtgL[idls]->SetMinimum(0);
  hyieldbdtgL[idls]->Draw("peX0");
  hyieldbdtgH[idls]->Draw("peX0 same");
  drawalltext();
  c3->cd(2);
  gStyle->SetErrorX(0);
  hemptyeff->Draw();
  greff_a[idls]->Draw("p3 same");
  greff_b[idls]->Draw("p3 same");
  drawalltext();
  c3->cd(3);
  hsigfbdtgH[idls]->SetMaximum(hsigfbdtgH[idls]->GetMaximum()*1.5);
  hsigfbdtgH[idls]->Draw("pl");
  hsigfbdtgL[idls]->Draw("pl same");
  drawalltext();
  c3->cd(4);
  hyieldbdtgBenrL[idls]->SetMaximum(hyieldbdtgBenrL[idls]->GetMaximum()*1.5);
  hyieldbdtgBenrL[idls]->SetMinimum(0);
  hyieldbdtgBenrL[idls]->Draw("peX0");
  hyieldbdtgBenrH[idls]->Draw("peX0 same");
  drawalltext();
  c3->cd(5);
  gStyle->SetErrorX(0);
  hemptyfprompt->Draw();
  grfpromptL->Draw("pe same");
  grfpromptH->Draw("pe same");
  drawalltext();
  c3->SaveAs(Form("plots/csig_varybdtg_%s.pdf", output.c_str()));

  TFile* outf = new TFile(Form("rootfiles/root_fitXvary_yield_%s.root", output.c_str()), "recreate");
  outf->cd();
  for(auto& hh : hyieldbdtgL) { hh->Write(); }
  for(auto& hh : hyieldbdtgH) { hh->Write(); }
  for(auto& hh : hbkgbdtgL) { hh->Write(); }
  for(auto& hh : hbkgbdtgH) { hh->Write(); }
  for(auto& hh : hsigfbdtgL) { hh->Write(); }
  for(auto& hh : hsigfbdtgH) { hh->Write(); }
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==2) { drawfitXvary(argv[1]); return 0; }
  return 1;
}
