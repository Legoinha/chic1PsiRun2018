#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>

#include <string>

#include "fitXvary.h"
#include "fit.h"

namespace fitXvary
{
  void fitXvaryfit(fitX::varymva* mvas, std::string title,
                   std::vector<TH1F*> hdata, TH1F* hmcp_a, TH1F* hmcp_b,
                   std::vector<RooDataSet*> dshdata, RooDataSet* dshmcp_a, RooDataSet* dshmcp_b,
                   TH1F* hyieldmva_a, TH1F* hyieldmva_b,
                   std::string output, TPad* pp, float fixmean)
  {
    pp->cd();
    xjjroot::setgstyle(1);
    // RooPlot* frempty = mass->frame(RooFit::Title(""));
    // frempty->SetYTitle(Form("Entries / %.0f MeV", fitX::BIN_WIDTH*1.e+3));
    for(int l=0, lp=0; l<mvas->n()-1; l++)
      {
        if(mvas->ifdraw()[l]) hdata[l]->Draw(Form("pe %s", (lp?"same":"")));
        std::map<std::string, fitX::fitXresult*> result = fitX::fit(hdata[l], 0, hmcp_a, hmcp_b,
                                                                    dshdata[l], dshmcp_a, dshmcp_b,
                                                                    Form("plots/%s/idx", output.c_str()), fixmean, true, Form("_%d", l), Form("%s > %.2f", mvas->type().c_str(), mvas->mva()[l])); // fix mean = false
        pp->cd();
        xjjroot::setgstyle(1);
        xjjroot::settfstyle(result["unbinned"]->f(), xjjroot::mycolor_middle[xjjroot::cc[lp]], 2, 4); //
        if(mvas->ifdraw()[l])
          {
            // xjjroot::copyobject(result["unbinned"]->f(), Form("ff_%s", hdata[l]->GetName()))->Draw("l same");
            result["unbinned"]->f()->Draw("l same");
            xjjroot::drawtex(0.70+0.08*(lp%3), 0.78-0.04*(lp/3), Form("> %.2f%s", mvas->mva()[l], (lp%3==2?"":", ")), 0.03, 12, 62, xjjroot::mycolor_middle[xjjroot::cc[lp]]);
            lp++;
          }
        hyieldmva_a->SetBinContent(l+1, result["unbinned"]->ysig_a());
        hyieldmva_a->SetBinError(l+1, result["unbinned"]->ysigerr_a());
        hyieldmva_b->SetBinContent(l+1, result["unbinned"]->ysig_b());
        hyieldmva_b->SetBinError(l+1, result["unbinned"]->ysigerr_b());
      }
    xjjroot::setgstyle(1);
    xjjroot::drawtex(0.90, 0.85, title.c_str(), 0.04, 33, 62, kBlack);
    xjjroot::drawtex(0.70-0.08, 0.78, mvas->type().c_str(), 0.03, 12, 62, kBlack);
    xjjroot::drawCMS();
  }

}
