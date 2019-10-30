#include <map>
#include <vector>
#include <string>

#include "fitX.h"

namespace tnpcc
{
  const int filterId = 1; // * filterId = 1: Jpsi L3 filter
  std::map<int, std::string> idxname = {
    std::pair<int, std::string>(-1, "syst_u"),
    std::pair<int, std::string>(-2, "syst_d"),
    std::pair<int, std::string>(+1, "stat_u"),
    std::pair<int, std::string>(+2, "stat_d"),
    std::pair<int, std::string>(0, "nominal"),
  };
  const std::vector<std::string> types = {"nominal", "trg", "trk", "muid", "total"};
  std::map<std::string, Color_t> typecolor = {
    std::pair<std::string, Color_t>("nominal", kBlack),
    std::pair<std::string, Color_t>("trg", xjjroot::mycolor_middle["azure"]),
    std::pair<std::string, Color_t>("trk", xjjroot::mycolor_middle["red"]),
    std::pair<std::string, Color_t>("muid", xjjroot::mycolor_middle["green"]),
    std::pair<std::string, Color_t>("total", xjjroot::mycolor_middle["yellow"]),
  };
  std::map<std::string, Style_t> idxstyle = {
    std::pair<std::string, Style_t>("syst_u", 2),
    std::pair<std::string, Style_t>("syst_d", 2),
    std::pair<std::string, Style_t>("stat_u", 4),
    std::pair<std::string, Style_t>("stat_d", 4),
    std::pair<std::string, Style_t>("nominal", 1),
  };
  std::vector<std::string> err({"stat", "syst"});

  __PTBIN_INPUT__
  int nptbins = sizeof(ptbins)/sizeof(ptbins[0]) - 1;
  float scalemin = 0.9, scalemax = 1.4;

  std::vector<float> muetabins = {0, 1.2, 2.1, 2.4}; int nmueta = muetabins.size()-1;

  class drawtnp
  {
  public:
    drawtnp(std::string inputname, std::string name, int color=0);
    std::map<std::string, std::map<std::string, TH1D*>> hh;
    std::map<std::string, std::map<std::string, TH1D*>> hhp;
    std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> gg;
    std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> ggp;
    TFile* inf() { return finf; }
    void draw();
    void drawscale();
    void drawleg() { leg->Draw(); }
    void drawkinematics(float ty);
    void setleg(float ty);
  private:
    TFile* finf;
    std::vector<float> yy, yye_d, yye_u, xx, xxe;
    std::vector<float> yyp, yype_d, yype_u;
    TLegend* leg;
  };

}

tnpcc::drawtnp::drawtnp(std::string inputname, std::string name, int color)
{
  finf = TFile::Open(inputname.c_str());
  
  for(auto& tt : tnpcc::types)
    {
      for(auto& idxk : tnpcc::idxname)
        {
          hh[tt][idxk.second] = (TH1D*)finf->Get(Form("htnp_%s_%s", tt.c_str(), idxk.second.c_str()));
          hh[tt][idxk.second]->SetName(Form("%s%s", hh[tt][idxk.second]->GetName(), name.c_str()));
          xjjroot::sethempty(hh[tt][idxk.second], 0, 0.3);
          xjjroot::setthgrstyle(hh[tt][idxk.second], (color?color:tnpcc::typecolor[tt]), 21, 0.8, (color?color:tnpcc::typecolor[tt]), tnpcc::idxstyle[idxk.second], 1, 0, 0, 1001);
        }
      hhp[tt]["nominal"] = (TH1D*)hh[tt]["nominal"]->Clone(Form("scale_%s", hh[tt]["nominal"]->GetName()));
      hhp[tt]["nominal"]->Divide(hh["nominal"]["nominal"]);
      for(int i=0; i<tnpcc::nptbins; i++) { hhp[tt]["nominal"]->SetBinError(i+1, 0); }
      hhp[tt]["nominal"]->GetYaxis()->SetTitle("#beta^{TnP}");
      for(auto& ee : tnpcc::err)
        {
          xx.clear(); xxe.clear(); yy.clear(); yye_d.clear(); yye_u.clear();
          yyp.clear(); yype_d.clear(); yype_u.clear();
          for(int i=0; i<tnpcc::nptbins; i++)
            {
              float x = (tnpcc::ptbins[i+1] + tnpcc::ptbins[i])/2.;
              float xe = (tnpcc::ptbins[i+1] - tnpcc::ptbins[i])/2.;
              float y = hh[tt]["nominal"]->GetBinContent(i+1);
              float ye_d = hh[tt]["nominal"]->GetBinContent(i+1) - hh[tt][ee+"_d"]->GetBinContent(i+1);
              float ye_u = 0-hh[tt]["nominal"]->GetBinContent(i+1) + hh[tt][ee+"_u"]->GetBinContent(i+1);
              float y0 = hh["nominal"]["nominal"]->GetBinContent(i+1);
              float yp = y/y0;
              float ype_d = ye_d/y0;
              float ype_u = ye_u/y0;
              xx.push_back(x);
              xxe.push_back(xe);
              yy.push_back(y);
              yye_d.push_back(ye_d);
              yye_u.push_back(ye_u);
              yyp.push_back(yp);
              yype_d.push_back(ype_d);
              yype_u.push_back(ype_u);
            }
          gg[tt][ee] = new TGraphAsymmErrors(tnpcc::nptbins, xx.data(), yy.data(), xxe.data(), xxe.data(), yye_d.data(), yye_u.data());
          gg[tt][ee]->SetName(Form("ggtnp_%s_%s%s", tt.c_str(), ee.c_str(), name.c_str()));
          if(ee=="stat") xjjroot::setthgrstyle(gg[tt][ee], (color?color:tnpcc::typecolor[tt]), 21, 0.8, 0, 0, 0, (color?color:tnpcc::typecolor[tt]), 0.3, 1001);
          if(ee=="syst") xjjroot::setthgrstyle(gg[tt][ee], (color?color:tnpcc::typecolor[tt]), 21, 0.8, (color?color:tnpcc::typecolor[tt]), 1, 1, 0, 0, 1001);
          ggp[tt][ee] = new TGraphAsymmErrors(tnpcc::nptbins, xx.data(), yyp.data(), xxe.data(), xxe.data(), yype_d.data(), yype_u.data());
          ggp[tt][ee]->SetName(Form("ggptnp_%s_%s%s", tt.c_str(), ee.c_str(), name.c_str()));
          if(ee=="stat") xjjroot::setthgrstyle(ggp[tt][ee], (color?color:tnpcc::typecolor[tt]), 21, 0.8, 0, 0, 0, (color?color:tnpcc::typecolor[tt]), 0.3, 1001);
          if(ee=="syst") xjjroot::setthgrstyle(ggp[tt][ee], (color?color:tnpcc::typecolor[tt]), 21, 0.8, (color?color:tnpcc::typecolor[tt]), 1, 1, 0, 0, 1001);
        }
    }

  if(tnpcc::nptbins < 3) hh["nominal"]["nominal"]->SetMaximum(hh["nominal"]["nominal"]->GetMaximum()*3.5);
  else hh["nominal"]["nominal"]->SetMaximum(hh["nominal"]["nominal"]->GetMaximum()*1.5);
  hh["nominal"]["nominal"]->SetMinimum(0.1);
  hhp["nominal"]["nominal"]->SetMaximum(1.4);
  hhp["nominal"]["nominal"]->SetMinimum(0.9);

}

void tnpcc::drawtnp::setleg(float ty)
{
  leg = new TLegend(0.60, ty-0.04-0.045*6, 0.89, ty-0.04, "Stat.          Syst.");
  xjjroot::setleg(leg, 0.04);
  leg->SetNColumns(2);
  for(auto& tt : tnpcc::types)
    {
      if(tt == "nominal") continue;
      for(auto& ee : tnpcc::err)
        {
          leg->AddEntry(ggp[tt][ee], tt.c_str(), "f");
        }
    }
  leg->AddEntry(hh["nominal"]["nominal"], "Nominal", "l");
}

void tnpcc::drawtnp::draw()
{
  hh["nominal"]["nominal"]->Draw("hist");
  for(auto& ee : tnpcc::err)
    {
      std::string drawopt = (ee=="stat"?"p2 same":"p5 same");
      gg["total"][ee]->Draw(drawopt.c_str());
      gg["trg"][ee]->Draw(drawopt.c_str());
      gg["muid"][ee]->Draw(drawopt.c_str());
      gg["trk"][ee]->Draw(drawopt.c_str());
    }
  hh["total"]["nominal"]->Draw("hist same");
  hh["nominal"]["nominal"]->Draw("hist same");
}

void tnpcc::drawtnp::drawscale()
{
  hhp["nominal"]["nominal"]->Draw("0");
  for(auto& ee : tnpcc::err)
    {
      std::string drawopt = (ee=="stat"?"2p same":"p5 same");
      ggp["total"][ee]->Draw(drawopt.c_str());
      ggp["trg"][ee]->Draw(drawopt.c_str());
      ggp["muid"][ee]->Draw(drawopt.c_str());
      ggp["trk"][ee]->Draw(drawopt.c_str());
    }
  hhp["total"]["nominal"]->Draw("hist same");
  hhp["nominal"]["nominal"]->Draw("hist same");
}

void tnpcc::drawtnp::drawkinematics(float ty)
{
  xjjroot::drawtex(0.89, ty, "PYTHIA + HYDJET", 0.04, 33, 62);
  xjjroot::drawtex(0.22, ty-0.043, "TnP Correction", 0.04, 13, 62);
  xjjroot::drawtex(0.22, ty-0.043*2, fitX::ytag().c_str(), 0.04, 13, 42);
  xjjroot::drawtex(0.22, ty-0.043*3, fitX::centtag().c_str(), 0.04, 13, 42);
}

