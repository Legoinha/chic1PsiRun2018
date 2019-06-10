#include "fitX.h"

#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <TH1F.h>
#include <string>

#include "packtree.h"
#include "ntuple.h"
#include "xjjcuti.h"
#include "xjjrootuti.h"

#include "fitXvary.h"

void drawfitXvary(std::string output)
{
  xjjroot::setgstyle(3);
  TFile* inf = new TFile(Form("rootfiles/root_fitXvary_%s.root", output.c_str()));
  std::vector<TH1F*> hdata(nbdtg*ndls), hmc_a(nbdtg*ndls), hmc_b(nbdtg*ndls);
  std::vector<TF1*> ff(nbdtg*ndls, 0), ffL(nbdtg*ndls, 0), ffH(nbdtg*ndls, 0), ffbkg(nbdtg*ndls, 0);
  for(int k=0; k<ndls; k++)
    {
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx] = (TH1F*)inf->Get(Form("hdata_%d_%d", k, l));
          fitX::setmasshist(hdata[idx], 0, -0.1);
          hdata[idx]->SetMinimum(0);
          xjjroot::setthgrstyle(hdata[idx], xjjroot::colorlist_dark[lp], 20, 0.9, xjjroot::colorlist_dark[lp], 1, 1);
          hmc_a[idx] = (TH1F*)inf->Get(Form("hmc_a_%d_%d", k, l));
          hmc_b[idx] = (TH1F*)inf->Get(Form("hmc_b_%d_%d", k, l));
          if(pbdtg[l]) lp++;
        }
    }

  std::vector<TH1F*> hyieldbdtgL(ndls), hyieldbdtgH(ndls), hbkgbdtgL(ndls), hbkgbdtgH(ndls), hsigfbdtgL(ndls), hsigfbdtgH(ndls);
  for(int k=0; k<ndls; k++)
    {
      hyieldbdtgL[k] = new TH1F(Form("hyieldbdtgL_%d", k), ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
      hyieldbdtgH[k] = new TH1F(Form("hyieldbdtgH_%d", k), ";BDTG;N_{Signal}", bdtg.size()-1, bdtg.data());
      hbkgbdtgL[k] = new TH1F(Form("hbkgbdtgL_%d", k), ";BDTG;N_{Background}", bdtg.size()-1, bdtg.data());
      hbkgbdtgH[k] = new TH1F(Form("hbkgbdtgH_%d", k), ";BDTG;N_{Background}", bdtg.size()-1, bdtg.data());
      hsigfbdtgL[k] = new TH1F(Form("hsigfbdtgL_%d", k), ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());
      hsigfbdtgH[k] = new TH1F(Form("hsigfbdtgH_%d", k), ";BDTG;S / #sqrt{S+B}", bdtg.size()-1, bdtg.data());
      TCanvas* c = new TCanvas("cvary", "", 800, 600);
      hdata[k*nbdtg]->Draw("pe");
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          if(pbdtg[l]) hdata[idx]->Draw("pe same");
          std::vector<TF1*> funs;
          bool saveornot = k==0;
          funs = fitX::fit(xjjroot::copyobject(hdata[idx], Form("hdata_%d_%d", k, l)), 0, hmc_a[0], hmc_b[0], Form("idx/hdata_%d_%d", k, l), true, saveornot); // fix mean
          ff[idx] = funs[0];
          ff[idx]->SetName(Form("ff_%d_%d", k, l));
          ffL[idx] = funs[1];
          ffL[idx]->SetName(Form("ffL_%d_%d", k, l));
          ffH[idx] = funs[2];
          ffH[idx]->SetName(Form("ffH_%d_%d", k, l));
          ffbkg[idx] = funs[3];
          ffbkg[idx]->SetName(Form("ffbkg_%d_%d", k, l));
          c->cd();
          xjjroot::settfstyle(ff[idx], xjjroot::colorlist_dark[lp], 2, 2);
          if(pbdtg[l]) 
            {
              ff[idx]->Draw("l same");
              // hdata[idx]->Draw("pe same");
              xjjroot::drawtex(0.70+0.08*(lp%3), 0.75-0.04*(lp/3), Form("> %.2f%s", bdtg[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[lp]);
              lp++;
            }
          if(l<nbdtg-1)
            {
              float ysigL = funs[1]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigL < 0) ysigL = 0;
              float ysigLerr = funs[0]->GetParError(5)*ysigL/funs[0]->GetParameter(5);
              float ysigH = funs[2]->Integral(fitX::BIN_MIN, fitX::BIN_MAX) / fitX::BIN_WIDTH; if(ysigH < 0) ysigH = 0;
              float ysigHerr = funs[0]->GetParError(10)*ysigH/funs[0]->GetParameter(10);
              std::cout<<"\e[33;1m"<<ysigL<<" "<<ysigH<<"\e[0m"<<std::endl;
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
      xjjroot::setgstyle(3);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
      c->SaveAs(Form("plots/cmass_varybdtg_%s_%d.pdf", output.c_str(), k));

      TCanvas* c2 = new TCanvas("cvary2", "", 800, 600);
      c2->cd();
      // hdata[k*nbdtg+3]->Draw("pe");
      for(int l=0, lp=0; l<nbdtg; l++)
        {
          std::string opt = lp==0?"pe":"pe same";
          if(!pbdtg[l]) continue;
          if(lp >= 3) 
            {
              int idx = k*nbdtg + l;
              hdata[idx]->Draw(opt.c_str());
              ff[idx]->Draw("l same");
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
    }

  const int idls = 0;
  xjjroot::sethempty(hyieldbdtgL[idls]);
  xjjroot::setthgrstyle(hyieldbdtgL[idls], fitX::color_a, 20, 1.2, fitX::color_a, 1, 1);
  xjjroot::sethempty(hyieldbdtgH[idls]);
  xjjroot::setthgrstyle(hyieldbdtgH[idls], fitX::color_b, 20, 1.2, fitX::color_b, 1, 1);
  xjjroot::sethempty(hsigfbdtgL[idls]);
  xjjroot::setthgrstyle(hsigfbdtgL[idls], fitX::color_a, 46, 1.2, fitX::color_a, 1, 1);
  xjjroot::sethempty(hsigfbdtgH[idls]);
  xjjroot::setthgrstyle(hsigfbdtgH[idls], fitX::color_b, 46, 1.2, fitX::color_b, 1, 1);
  xjjroot::setgstyle(1 );
  TCanvas* c3 = new TCanvas("c3", "", 1800, 600);
  c3->Divide(3, 1);
  c3->cd(1);
  hyieldbdtgL[idls]->SetMaximum(hyieldbdtgL[idls]->GetMaximum()*1.5);
  hyieldbdtgL[idls]->SetMinimum(0);
  hyieldbdtgL[idls]->Draw("peX0");
  hyieldbdtgH[idls]->Draw("peX0 same");
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.90, 0.84, "p_{T} < 15 GeV/c", 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.042, "|y| < 1.5", 0.038, 32, 62);
  c3->cd(2);
  hsigfbdtgH[idls]->SetMaximum(hsigfbdtgH[idls]->GetMaximum()*1.5);
  hsigfbdtgH[idls]->Draw("pl");
  hsigfbdtgL[idls]->Draw("pl same");
  xjjroot::drawCMSleft();
  xjjroot::drawCMSright();
  xjjroot::drawtex(0.24, 0.84, "#psi(2S)", 0.038, 12, 62, fitX::color_a);
  xjjroot::drawtex(0.24, 0.84-0.042, "X(3872)", 0.038, 12, 62, fitX::color_b);
  xjjroot::drawtex(0.90, 0.84, "p_{T} < 15 GeV/c", 0.038, 32, 62);
  xjjroot::drawtex(0.90, 0.84-0.042, "|y| < 1.5", 0.038, 32, 62);
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
