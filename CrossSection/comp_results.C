#include <TFile.h>
#include <TCanvas.h>

#include <string>
#include <map>

#include "xjjrootuti.h"
#include "fitX.h"
#include "results.h"

void comp_results()
{
  std::map<std::string, std::string> options = {
    std::pair<std::string, std::string>("_wmissingfile", "PAS"),
    std::pair<std::string, std::string>("", "miss MC"),
    std::pair<std::string, std::string>("_PVz15", "PVz"),
    std::pair<std::string, std::string>("_PVz15_L2L3", "TnP L2L3"),
    std::pair<std::string, std::string>("_PVz15_newL2L3", "New TnP"),
  };
  std::vector<std::string> options_sort = {"_wmissingfile", "", "_PVz15", "_PVz15_L2L3", "_PVz15_newL2L3"};
  std::map<std::string, TFile*> inf;
  std::map<std::string, fitX::results*> rs;
  for(auto& i : options)
    {
      inf[i.first] = TFile::Open(std::string("rootfiles/trainX_20190808ptdep_sideband_tktk0p2_15p0_50p0_0-10-1-2-9_BDTQvalue" + i.first + "_pt15-50_cent090_y0p0-1p6/fitX_fithist.root").c_str());
      rs[i.first] = new fitX::results(inf[i.first], i.first);
    }
  fitX::results rs0;
  std::map<std::string, TH1F*> hdev;
  for(auto& r : rs0.vars())
    {
      hdev[r+"_a"] = new TH1F(Form("h%s_a", r.c_str()), ";;Change w.r.t. PAS", options.size(), 0, options.size());
      hdev[r+"_b"] = new TH1F(Form("h%s_b", r.c_str()), ";;Change w.r.t. PAS", options.size(), 0, options.size());
    }
  hdev["ratio"] = new TH1F("hratio", ";;Change w.r.t PAS (%)", options.size(), 0, options.size());
  for(int j=0; j<options_sort.size(); j++)
    {
      for(auto& r : rs0.vars())
        {
          hdev[r+"_a"]->SetBinContent(j+1, rs[options_sort[j]]->val(r+"_a"));
          hdev[r+"_a"]->SetBinError(j+1, rs[options_sort[j]]->err(r+"_a"));
          hdev[r+"_b"]->SetBinContent(j+1, rs[options_sort[j]]->val(r+"_b"));
          hdev[r+"_b"]->SetBinError(j+1, rs[options_sort[j]]->err(r+"_b"));
        }
      hdev["ratio"]->SetBinContent(j+1, rs[options_sort[j]]->val("ratio"));
      hdev["ratio"]->SetBinError(j+1, rs[options_sort[j]]->err("ratio"));
    }

  for(auto& h : hdev)
    {
      h.second->Scale(1./h.second->GetBinContent(1));
      xjjroot::sethempty(h.second, 0, 0.3);
      h.second->SetMinimum(0.90);
      h.second->SetMaximum(1.10);
      for(int j=0; j<options_sort.size(); j++)
        {
          h.second->GetXaxis()->SetBinLabel(j+1, options[options_sort[j]].c_str());
          h.second->GetXaxis()->SetLabelSize(h.second->GetXaxis()->GetLabelSize()*1.02);
        }
      h.second->GetYaxis()->SetNdivisions(-404);
    }
  
  xjjroot::setthgrstyle(hdev["yieldpromptCorr_a"], xjjroot::mycolor_satmiddle["azure"], 20, 1.1, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::setthgrstyle(hdev["yieldpromptCorr_b"], xjjroot::mycolor_satmiddle["green"], 20, 1.1, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::setthgrstyle(hdev["ratio"], xjjroot::mycolor_satmiddle["red"], 21, 1.1, xjjroot::mycolor_satmiddle["red"], 1, 1);
  xjjroot::setthgrstyle(hdev["yield_a"], xjjroot::mycolor_satmiddle["azure"], 20, 1.1, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::setthgrstyle(hdev["yield_b"], xjjroot::mycolor_satmiddle["green"], 20, 1.1, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::setthgrstyle(hdev["Benryield_a"], xjjroot::mycolor_satmiddle["azure"], 21, 1.1, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::setthgrstyle(hdev["Benryield_b"], xjjroot::mycolor_satmiddle["green"], 21, 1.1, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::setthgrstyle(hdev["eff_a"], xjjroot::mycolor_satmiddle["azure"], 47, 1.1, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::setthgrstyle(hdev["eff_b"], xjjroot::mycolor_satmiddle["green"], 47, 1.1, xjjroot::mycolor_satmiddle["green"], 1, 1);
  xjjroot::setthgrstyle(hdev["tnp_a"], xjjroot::mycolor_satmiddle["azure"], 29, 1.1, xjjroot::mycolor_satmiddle["azure"], 1, 1);
  xjjroot::setthgrstyle(hdev["tnp_b"], xjjroot::mycolor_satmiddle["green"], 29, 1.1, xjjroot::mycolor_satmiddle["green"], 1, 1);
  TLegend* leg = new TLegend(0.20, 0.85-0.04, 0.50, 0.85);
  xjjroot::setleg(leg, 0.035);
  leg->AddEntry(hdev["ratio"], Form("%s", rs0.leg("ratio").c_str(), fitX::title_a.c_str()), "p");
  TLegend* leg1 = new TLegend(0.20, 0.85-0.04*3, 0.50, 0.85);
  xjjroot::setleg(leg1, 0.035);
  leg1->AddEntry(hdev["yieldpromptCorr_a"], Form("%s %s", rs0.leg("yieldpromptCorr").c_str(), fitX::title_a.c_str()), "p");
  leg1->AddEntry(hdev["yieldpromptCorr_b"], Form("%s %s", rs0.leg("yieldpromptCorr").c_str(), fitX::title_b.c_str()), "p");
  leg1->AddEntry(hdev["ratio"], Form("%s", rs0.leg("ratio").c_str(), fitX::title_a.c_str()), "p");
  TLegend* leg2 = new TLegend(0.20, 0.85-0.03*8, 0.50, 0.85);
  xjjroot::setleg(leg2, 0.03);
  leg2->AddEntry(hdev["yield_a"], Form("%s %s", rs0.leg("yield").c_str(), fitX::title_a.c_str()), "p");
  leg2->AddEntry(hdev["yield_b"], Form("%s %s", rs0.leg("yield").c_str(), fitX::title_b.c_str()), "p");
  leg2->AddEntry(hdev["Benryield_a"], Form("%s %s", rs0.leg("Benryield").c_str(), fitX::title_a.c_str()), "p");
  leg2->AddEntry(hdev["Benryield_b"], Form("%s %s", rs0.leg("Benryield").c_str(), fitX::title_b.c_str()), "p");
  leg2->AddEntry(hdev["eff_a"], Form("%s %s", rs0.leg("eff").c_str(), fitX::title_a.c_str()), "p");
  leg2->AddEntry(hdev["eff_b"], Form("%s %s", rs0.leg("eff").c_str(), fitX::title_b.c_str()), "p");
  leg2->AddEntry(hdev["tnp_a"], Form("%s %s", rs0.leg("tnp").c_str(), fitX::title_a.c_str()), "p");
  leg2->AddEntry(hdev["tnp_b"], Form("%s %s", rs0.leg("tnp").c_str(), fitX::title_b.c_str()), "p");

  xjjroot::setgstyle(1);
  gStyle->SetPaintTextFormat("1.3f");
  TCanvas* c = new TCanvas("c", "", 600, 600);
  hdev["ratio"]->Draw("hist pl text0 same");
  leg->Draw();
  xjjroot::drawCMS();
  xjjroot::mkdir("plots/comp_results/x");
  c->SaveAs("plots/comp_results/cratio.pdf");

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  hdev["yieldpromptCorr_a"]->Draw("hist pl text0");
  hdev["yieldpromptCorr_b"]->Draw("hist pl text0 same");
  hdev["ratio"]->Draw("hist pl text0 same");
  leg1->Draw();
  xjjroot::drawCMS();
  xjjroot::mkdir("plots/comp_results/x");
  c1->SaveAs("plots/comp_results/cyieldpromptCorr.pdf");

  TCanvas* c2 = new TCanvas("c2", "", 600, 600);
  hdev["yield_a"]->Draw("hist pl text0");
  hdev["yield_b"]->Draw("hist pl text0 same");
  hdev["Benryield_a"]->Draw("hist pl text0 same");
  hdev["Benryield_b"]->Draw("hist pl text0 same");
  hdev["eff_a"]->Draw("hist pl text0 same");
  hdev["eff_b"]->Draw("hist pl text0 same");
  hdev["tnp_a"]->Draw("hist pl text0 same");
  hdev["tnp_b"]->Draw("hist pl text0 same");
  leg2->Draw();
  xjjroot::drawCMS();
  xjjroot::mkdir("plots/comp_results/x");
  c2->SaveAs("plots/comp_results/cyield.pdf");
}

int main()
{
  comp_results();
  return 0;
}
