#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TLegend.h>
#include "xjjrootuti.h"
#include "xjjcuti.h"
#include <string>
#include "tnpcc.h"

void draw_tnp(std::string inputname)
{
  std::string tparticle("untitled");
  if(xjjc::str_contains(inputname, "romptX")) tparticle = "X(3872)";
  if(xjjc::str_contains(inputname, "romptPsi")) tparticle = "#psi(2S)";
  TFile* inf = TFile::Open(inputname.c_str());
  std::map<std::string, std::map<std::string, TH1D*>> hh;
  std::map<std::string, std::map<std::string, TH1D*>> hhp;
  std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> gg;
  std::map<std::string, std::map<std::string, TGraphAsymmErrors*>> ggp;
  std::vector<float> yy, yye_d, yye_u, xx, xxe;
  std::vector<float> yyp, yype_d, yype_u;
  for(auto& tt : tnpcc::types)
    {
      for(auto& idxk : tnpcc::idxname)
        {
          hh[tt][idxk.second] = (TH1D*)inf->Get(Form("htnp_%s_%s", tt.c_str(), idxk.second.c_str()));
          TH1D* h = hh[tt][idxk.second];
          xjjroot::sethempty(h, 0, 0.3);
          xjjroot::setthgrstyle(h, tnpcc::typecolor[tt], 21, 0.8, tnpcc::typecolor[tt], tnpcc::idxstyle[idxk.second], 1, 0, 0, 1001);
        }
      hhp[tt]["nominal"] = (TH1D*)hh[tt]["nominal"]->Clone(Form("ratio_%s", hh[tt]["nominal"]->GetName()));
      hhp[tt]["nominal"]->Divide(hh["nominal"]["nominal"]);
      for(int i=0; i<tnpcc::nptbins; i++) { hhp[tt]["nominal"]->SetBinError(i+1, 0); }
      hhp[tt]["nominal"]->GetYaxis()->SetTitle("Scale");
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
          gg[tt][ee]->SetName(Form("ggtnp_%s_%s", tt.c_str(), ee.c_str()));
          if(ee=="stat") xjjroot::setthgrstyle(gg[tt][ee], tnpcc::typecolor[tt], 21, 0.8, 0, 0, 0, tnpcc::typecolor[tt], 0.3, 1001);
          if(ee=="syst") xjjroot::setthgrstyle(gg[tt][ee], tnpcc::typecolor[tt], 21, 0.8, tnpcc::typecolor[tt], 1, 1, 0, 0, 1001);
          ggp[tt][ee] = new TGraphAsymmErrors(tnpcc::nptbins, xx.data(), yyp.data(), xxe.data(), xxe.data(), yype_d.data(), yype_u.data());
          ggp[tt][ee]->SetName(Form("ggptnp_%s_%s", tt.c_str(), ee.c_str()));
          if(ee=="stat") xjjroot::setthgrstyle(ggp[tt][ee], tnpcc::typecolor[tt], 21, 0.8, 0, 0, 0, tnpcc::typecolor[tt], 0.3, 1001);
          if(ee=="syst") xjjroot::setthgrstyle(ggp[tt][ee], tnpcc::typecolor[tt], 21, 0.8, tnpcc::typecolor[tt], 1, 1, 0, 0, 1001);
        }
    }

  float ty = 0.85, tx = 0.89;
  xjjroot::setgstyle(1);
  TCanvas* c = new TCanvas("c", "", 1200, 600);
  c->Divide(2, 1);
  c->cd(1);
  // gPad->SetLogy();
  hh["nominal"]["nominal"]->SetMaximum(hh["nominal"]["nominal"]->GetMaximum()*1.5);
  hh["nominal"]["nominal"]->SetMinimum(0.1);
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
  xjjroot::drawCMS("Simulation");
  gPad->RedrawAxis();
  c->cd(2);
  hhp["nominal"]["nominal"]->SetMaximum(1.4);
  hhp["nominal"]["nominal"]->SetMinimum(0.9);
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
  xjjroot::drawCMS("Simulation");
  gPad->RedrawAxis();

  TLegend* leg = new TLegend(0.60, ty-0.04-0.045*6, 0.89, ty-0.04, "Stat.          Syst.");
  xjjroot::setleg(leg, 0.04);
  leg->SetNColumns(2);
  for(auto& tt : tnpcc::types)
    {
      if(tt == "nominal") continue;
      for(auto& ee : tnpcc::err)
        {
          leg->AddEntry(gg[tt][ee], tt.c_str(), "f");
        }
    }
  leg->AddEntry(hh["nominal"]["nominal"], "Nominal", "l");
  c->cd(1);
  xjjroot::drawtex(0.89, ty, "PYTHIA + HYDJET", 0.04, 33, 62);
  xjjroot::drawtex(0.22, ty, tparticle.c_str(), 0.04, 13, 62);
  xjjroot::drawtex(0.22, ty-0.043, "TnP Correction", 0.04, 13, 62);
  leg->Draw();
  c->cd(2);
  xjjroot::drawtex(0.89, ty, "PYTHIA + HYDJET", 0.04, 33, 62);
  xjjroot::drawtex(0.22, ty, tparticle.c_str(), 0.04, 13, 62);
  xjjroot::drawtex(0.22, ty-0.043, "TnP Correction", 0.04, 13, 62);
  leg->Draw();
  c->cd();
  std::string outputname = xjjc::str_replaceall(xjjc::str_replaceall(inputname, "rootfiles/", "plots/c"), ".root", ".pdf");
  gSystem->Exec("mkdir -p plots");
  c->SaveAs(outputname.c_str());
  
}

int main(int argc, char* argv[])
{
  if(argc==2) { draw_tnp(argv[1]); return 0; }
  return 1;
}
