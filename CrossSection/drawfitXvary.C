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
  std::vector<TF1*> ff(nbdtg*ndls, 0);
  for(int k=0; k<ndls; k++)
    {
      for(int l=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx] = (TH1F*)inf->Get(Form("hdata_%d_%d", k, l));
          fitX::setmasshist(hdata[idx], 0, -0.1);
          hdata[idx]->SetMinimum(0);
          xjjroot::setthgrstyle(hdata[idx], xjjroot::colorlist_dark[l], 20, 0.9, xjjroot::colorlist_dark[l], 1, 1);
          hmc_a[idx] = (TH1F*)inf->Get(Form("hmc_a_%d_%d", k, l));
          hmc_b[idx] = (TH1F*)inf->Get(Form("hmc_b_%d_%d", k, l));
        }
    }

  for(int k=0; k<ndls; k++)
    {
      TCanvas* c = new TCanvas("cvary", "", 800, 600);
      hdata[k*nbdtg]->Draw("pe");
      for(int l=0; l<nbdtg; l++)
        {
          int idx = k*nbdtg + l;
          hdata[idx]->Draw("pe same");
          ff[idx] = fitX::fit(xjjroot::copyobject(hdata[idx], Form("hdata_%d_%d", 0, l)), 0, hmc_a[0], hmc_b[0], Form("idx_%d_%d", 0, l), false);
          ff[idx]->SetName(Form("ff_%d_%d", 0, l));
          c->cd();
          xjjroot::settfstyle(ff[idx], xjjroot::colorlist_middle[l], 2, 2);
          ff[idx]->Draw("l same");
          // hdata[idx]->Draw("pe same");
          xjjroot::drawtex(0.70+0.08*(l%3), 0.75-0.04*(l/3), Form("> %.2f%s", bdtg[l], (l%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[l]);
        }
      xjjroot::setgstyle(3);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
      c->SaveAs(Form("plots/cmass_varybdtg_%s_%d.pdf", output.c_str(), k));

      TCanvas* c2 = new TCanvas("cvary2", "", 800, 600);
      c2->cd();
      hdata[k*nbdtg+3]->Draw("pe");
      for(int l=0; l<nbdtg; l++)
        {
          if(l < 3) continue;
          int idx = k*nbdtg + l;
          hdata[idx]->Draw("pe same");
          ff[idx]->Draw("l same");
          xjjroot::drawtex(0.70+0.08*(l%3), 0.75-0.04*(l/3-1), Form("> %.2f%s", bdtg[l], (l%3==2?"":", ")), 0.03, 12, 62, xjjroot::colorlist_dark[l]);
        }
      xjjroot::setgstyle(3);
      xjjroot::drawtex(0.70-0.08, 0.75, "BDTG", 0.03, 12, 62, kBlack);
      xjjroot::drawCMSleft();
      xjjroot::drawCMSright();
      c2->SaveAs(Form("plots/cmass_varybdtg_%s_%d_zoomin.pdf", output.c_str(), k));

      delete c;
      delete c2;
    }
}

int main(int argc, char* argv[])
{
  if(argc==2) { drawfitXvary(argv[1]); return 0; }
  return 1;
}
